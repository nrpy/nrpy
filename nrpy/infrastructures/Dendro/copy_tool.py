#!/usr/bin/env python3
"""
Safe installer for the generated FCCZ4_GR module (whitepaper section 12.6).

This module is the single source of the generated project's
``tools/copy_into_dendro.py``: the exporter renders that file from this
source with a generated-file banner.

Default behavior is a dry run.  Supported options:
  --dendro-root PATH          Dendro-GR checkout root
  --check                     Verify the staged module and exit (no changes)
  --execute                   Apply the installation transactionally
  --apply-root-cmake-patch    Insert the root CMake marker block (idempotent)
  --backup-dir PATH           Back up an existing installed module
  --force-generated-replace   Replace a modified existing FCCZ4_GR (requires
                              --backup-dir, so the local state is preserved)
  --remove                    Remove exactly what the receipt identifies
  --verify-only               Verify an installed copy against its receipt

Transaction (on --execute):
  1. verify the staged module hashes against manifest/file_hashes.json;
  2. stage the module in a temporary sibling directory;
  3. back up an existing installed module when explicitly allowed;
  4. atomically rename the staged module into place as <root>/FCCZ4_GR;
  5. write the installation receipt (FCCZ4_GR/.nrpy-generated.json);
  6. update the root CMakeLists.txt through a temporary file;
  7. run post-install verification;
  8. print the exact changed paths.

Removal (--remove) only deletes files identified by the receipt and only
removes the exact root CMake marker block.  The complete installed file set is
validated *before* any unlink, so a refusal never leaves a half-removed
module.  All file paths in the receipt are module-relative, so the installed
copy is self-verifying independent of the generated project.
"""

import argparse
import hashlib
import json
import os
import shutil
import sys
import tempfile
from typing import Dict, List, NoReturn, Optional, cast

# Installed location under the Dendro-GR root (matches the root CMake patch).
MODULE_REL = "FCCZ4_GR"
# Location of the staged module inside the generated project.
STAGED_PREFIX = "Dendro-GR/FCCZ4_GR/"
MARKER_BEGIN = "# BEGIN NRPY FCCZ4 MODULE"
MARKER_END = "# END NRPY FCCZ4 MODULE"


def fail(message: str) -> NoReturn:
    """
    Print a failure message and exit nonzero.

    :param message: The failure description.
    """
    print(f"FAIL: {message}")
    sys.exit(1)


def sha256_of(path: str) -> str:
    """
    Return the SHA-256 hex digest of a file.

    :param path: File to digest.
    :return: Hex digest.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read(path: str) -> str:
    """
    Read a UTF-8 text file.

    :param path: File to read.
    :return: The file contents.
    """
    with open(path, "r", encoding="utf-8") as handle:
        return handle.read()


def project_root() -> str:
    """
    Return the generated project root (two levels above this tool).

    :return: The project root path.
    """
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def staged_module(project: str) -> str:
    """
    Return the staged module directory inside the generated project.

    :param project: Generated project root.
    :return: The staged module path.
    """
    return os.path.join(project, STAGED_PREFIX.rstrip("/"))


def installed_module(dendro_root: str) -> str:
    """
    Return the installed module directory under the Dendro-GR root.

    :param dendro_root: Dendro-GR checkout root.
    :return: The installed module path.
    """
    return os.path.join(dendro_root, MODULE_REL)


def preconditions(dendro_root: str) -> str:
    """
    Check that the target root is plausible and safe.

    :param dendro_root: Dendro-GR checkout root.
    :return: The resolved real path.
    """
    real = os.path.realpath(dendro_root)
    if real in (os.path.realpath("/"), os.path.realpath(os.path.expanduser("~"))):
        fail(f"refusing to install into {dendro_root}")
    if not os.path.isfile(os.path.join(real, "CMakeLists.txt")):
        fail(f"{dendro_root} does not look like a Dendro-GR root (no CMakeLists.txt)")
    return real


def read_receipt(dendro_root: str) -> Optional[Dict[str, object]]:
    """
    Return the installation receipt, or None if absent.

    :param dendro_root: Dendro-GR checkout root.
    :return: The parsed receipt, or None.
    """
    path = os.path.join(installed_module(dendro_root), ".nrpy-generated.json")
    if not os.path.isfile(path):
        return None
    return cast(Dict[str, object], json.loads(read(path)))


def receipt_files(receipt: Dict[str, object]) -> Dict[str, str]:
    """
    Return the module-relative path-to-digest map recorded in a receipt.

    :param receipt: The parsed installation receipt.
    :return: Mapping of module-relative path to SHA-256 hex digest.
    """
    return cast(Dict[str, str], receipt["files"])


def staged_file_hashes(project: str) -> Dict[str, str]:
    """
    Return the module-relative hash map, verifying the staged module.

    The manifest file_hashes.json is project-rooted; the module files are
    those under ``Dendro-GR/FCCZ4_GR/``.  They are re-expressed as
    module-relative keys (the keys the receipt uses) so the installed copy
    is self-verifying.

    :param project: The generated project root.
    :return: Mapping of module-relative path to SHA-256 hex digest.
    """
    file_hashes = json.loads(
        read(os.path.join(project, "manifest", "file_hashes.json"))
    )["files"]
    module = staged_module(project)
    out: Dict[str, str] = {}
    for rel, expected in sorted(file_hashes.items()):
        if not rel.startswith(STAGED_PREFIX):
            continue
        mod_rel = rel[len(STAGED_PREFIX) :]
        path = os.path.join(module, mod_rel)
        if not os.path.isfile(path):
            fail(f"staged file missing: {rel}")
        actual = sha256_of(path)
        if actual != expected:
            fail(f"staged hash mismatch for {rel}")
        out[mod_rel] = expected
    if not out:
        fail("no staged module files found under Dendro-GR/FCCZ4_GR")
    return out


def marker_lines() -> List[str]:
    """
    Return the exact root CMake marker block lines.

    :return: The marker block, one line per entry.
    """
    return [
        MARKER_BEGIN,
        'option(DENDRO_BUILD_FCCZ4 "Build NRPy-generated fCCZ4 solver" OFF)',
        "if(DENDRO_BUILD_FCCZ4)",
        "  add_subdirectory(FCCZ4_GR)",
        "endif()",
        MARKER_END,
    ]


def apply_root_cmake(dendro_root: str, changed: List[str]) -> None:
    """
    Insert (idempotently) the root CMake marker block via a temp file.

    :param dendro_root: Dendro-GR checkout root.
    :param changed: Accumulator of changed paths, appended to in place.
    :raises Exception: Re-raised after the temporary file is cleaned up.
    """
    root_cmake = os.path.join(dendro_root, "CMakeLists.txt")
    lines = read(root_cmake).splitlines()
    if MARKER_BEGIN in lines:
        start = lines.index(MARKER_BEGIN)
        if MARKER_END not in lines[start:] or lines.index(MARKER_END) < start:
            fail("existing root CMake marker block is malformed")
        end = lines.index(MARKER_END)
        if lines[start : end + 1] != marker_lines():
            fail("a different root CMake marker block already exists")
        print("note: root CMake marker block already present (idempotent no-op)")
        return
    fd, tmp_path = tempfile.mkstemp(dir=dendro_root, prefix=".fccz4_cmake_")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write("\n".join(lines + [""] + marker_lines()) + "\n")
        os.replace(tmp_path, root_cmake)
        changed.append(os.path.join(dendro_root, "CMakeLists.txt"))
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise


def _verify_installed(dendro_root: str, file_hashes: Dict[str, str]) -> None:
    """
    Verify every installed module file against the module-relative map.

    :param dendro_root: Dendro-GR checkout root.
    :param file_hashes: Module-relative path to SHA-256 digest.
    """
    module = installed_module(dendro_root)
    for rel, expected in sorted(file_hashes.items()):
        path = os.path.join(module, rel)
        if not os.path.isfile(path) or sha256_of(path) != expected:
            fail(f"post-install verification failed for {rel}")


def execute(project: str, dendro_root: str, args: argparse.Namespace) -> None:
    """
    Apply the transactional installation (section 12.6).

    :param project: Generated project root.
    :param dendro_root: Dendro-GR checkout root.
    :param args: Parsed command line.
    :raises Exception: Re-raised after the staged directory is cleaned up.
    """
    file_hashes = staged_file_hashes(project)
    module_dst = installed_module(dendro_root)
    if args.force_generated_replace and not args.backup_dir:
        fail(
            "--force-generated-replace overwrites a locally modified module; "
            "pass --backup-dir PATH so the existing tree is preserved"
        )
    changed: List[str] = []
    existing = read_receipt(dendro_root)
    if os.path.isdir(module_dst):
        if existing is None:
            fail("existing FCCZ4_GR without an NRPy receipt; refusing")
        if not args.force_generated_replace:
            for rel, expected in receipt_files(existing).items():
                path = os.path.join(module_dst, rel)
                if not os.path.isfile(path) or sha256_of(path) != expected:
                    fail(
                        "existing FCCZ4_GR differs from its receipt; refusing "
                        "without --force-generated-replace"
                    )
            print("note: existing FCCZ4_GR matches its receipt (idempotent)")
    if args.backup_dir and os.path.isdir(module_dst):
        backup = os.path.join(args.backup_dir, MODULE_REL)
        os.makedirs(args.backup_dir, exist_ok=True)
        shutil.rmtree(backup, ignore_errors=True)
        shutil.copytree(module_dst, backup)
        print(f"backup: {backup}")
    # Stage in a temporary sibling, then atomically rename.
    staged_dir = tempfile.mkdtemp(dir=dendro_root, prefix=".fccz4_staged_")
    try:
        shutil.rmtree(staged_dir)
        shutil.copytree(staged_module(project), staged_dir)
        receipt = {
            "module_abi_hash": json.loads(
                read(os.path.join(project, "manifest", "module.json"))
            )["module_abi_hash"],
            "files": file_hashes,
            "root_cmake_marker": marker_lines(),
        }
        with open(
            os.path.join(staged_dir, ".nrpy-generated.json"), "w", encoding="utf-8"
        ) as handle:
            json.dump(receipt, handle, sort_keys=True, indent=2)
            handle.write("\n")
        if os.path.isdir(module_dst):
            shutil.rmtree(module_dst)
        os.replace(staged_dir, module_dst)
        changed.append(os.path.join(dendro_root, MODULE_REL))
    except Exception:
        if os.path.isdir(staged_dir):
            shutil.rmtree(staged_dir, ignore_errors=True)
        raise
    if args.apply_root_cmake_patch:
        apply_root_cmake(dendro_root, changed)
    _verify_installed(dendro_root, file_hashes)
    print(f"installed FCCZ4_GR into {dendro_root}")
    for path in changed:
        print(f"changed: {path}")


def remove(dendro_root: str, _args: argparse.Namespace) -> None:
    """
    Remove exactly the receipt-identified files and marker block.

    :param dendro_root: Dendro-GR checkout root.
    :param _args: Parsed command line (unused).
    :raises Exception: Re-raised after the temporary file is cleaned up.
    """
    receipt = read_receipt(dendro_root)
    if receipt is None:
        fail("no installation receipt found; nothing to remove")
    module_dst = installed_module(dendro_root)
    # Pre-validate the COMPLETE installed file set before deleting anything:
    # every file must be receipt-identified and hash-match, and no unexpected
    # file may exist.  A modified install is refused before any unlink, so a
    # refusal never leaves a half-removed (receiptless) module (section 12.6).
    expected_rel = set(receipt_files(receipt)) | {".nrpy-generated.json"}
    on_disk_rel = set()
    for dirpath, _dirnames, filenames in os.walk(module_dst):
        for filename in filenames:
            on_disk_rel.add(
                os.path.relpath(os.path.join(dirpath, filename), module_dst).replace(
                    os.sep, "/"
                )
            )
    for rel in sorted(on_disk_rel - expected_rel):
        fail(f"unexpected file in installed module: {rel}; refusing to remove")
    for rel, expected in receipt_files(receipt).items():
        path = os.path.join(module_dst, rel)
        if not os.path.isfile(path) or sha256_of(path) != expected:
            fail(f"installed {rel} differs from the receipt; refusing to remove")
    if not os.path.isfile(os.path.join(module_dst, ".nrpy-generated.json")):
        fail("installation receipt missing; refusing to remove")
    for rel in sorted(receipt_files(receipt)):
        os.unlink(os.path.join(module_dst, rel))
    os.unlink(os.path.join(module_dst, ".nrpy-generated.json"))
    shutil.rmtree(module_dst)
    print("removed FCCZ4_GR (receipt-identified files only)")
    root_cmake = os.path.join(dendro_root, "CMakeLists.txt")
    lines = read(root_cmake).splitlines()
    if MARKER_BEGIN in lines and MARKER_END in lines:
        start = lines.index(MARKER_BEGIN)
        end = lines.index(MARKER_END)
        new_lines = lines[:start] + lines[end + 1 :]
        fd, tmp_path = tempfile.mkstemp(dir=dendro_root, prefix=".fccz4_cmake_")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write("\n".join(new_lines).rstrip() + "\n")
            os.replace(tmp_path, root_cmake)
            print(f"removed root CMake marker block from {root_cmake}")
        except Exception:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
            raise


def main() -> None:
    """Run the installer CLI (dry run by default)."""
    parser = argparse.ArgumentParser(
        description="Safe installer for the generated FCCZ4_GR module."
    )
    parser.add_argument("--dendro-root", required=True)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--apply-root-cmake-patch", action="store_true")
    parser.add_argument("--backup-dir")
    parser.add_argument("--force-generated-replace", action="store_true")
    parser.add_argument("--remove", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()

    project = project_root()
    if not os.path.isfile(os.path.join(project, "manifest", "file_hashes.json")):
        fail("project manifest not found (run from inside a generated project)")
    dendro_root = preconditions(args.dendro_root)
    if args.verify_only:
        receipt = read_receipt(dendro_root)
        if receipt is None:
            fail("no installation receipt found")
        _verify_installed(dendro_root, receipt_files(receipt))
        print(f"OK: installed FCCZ4_GR matches its receipt in {dendro_root}")
        return
    if args.remove:
        remove(dendro_root, args)
        return
    if args.check:
        staged_file_hashes(project)
        print(f"OK: staged module verifies against the manifest for {dendro_root}")
        return
    if args.execute:
        execute(project, dendro_root, args)
        return
    # Default: dry run.
    staged_file_hashes(project)
    receipt = read_receipt(dendro_root)
    print(f"DRY RUN: would install FCCZ4_GR into {dendro_root}")
    if receipt is not None:
        print("existing receipt found; re-running --execute is idempotent")
    print("use --execute (and --apply-root-cmake-patch) to apply")


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
