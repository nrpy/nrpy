#!/usr/bin/env python3
"""
Verify a generated dendro_fccz4 project (whitepaper sections 12.5/16.1).

This module is the single source of the generated project's
``tools/verify_generated_project.py``: the exporter renders that file from
this source with a generated-file banner.

Checks (all must pass, exit 0):
  * every file listed in ``manifest/file_hashes.json`` and ``manifest/SHA256SUMS``
    exists and its SHA-256 matches, and every file on disk is listed (the three
    hash artifacts cannot list themselves);
  * every generated C++/CMake artifact carries the ownership banner and the
    module ABI (section 12.7);
  * every manifest parses as JSON;
  * the CMake generated source list equals, as a set, the source path derived
    from every registered CFunction record (one source per registered
    CFunction, sections 11.5/12.1);
  * no concrete EVOL name -- bare or role-prefixed (``in_``/``rhs_``/``out_``/
    ``diag_``/``aux_``) -- and no NRPy-emitted numerical loop appears in a
    scaffold file (section 9.16).  Section 9.16's mandated scan is the
    EVOL-name one; the loop check is an additional heuristic keyed on the NRPy
    loop-helper index names.  Three of the section's prohibitions -- a physics
    parameter default, an FD coefficient, and a loop written with other index
    names -- are not mechanically enforced here and remain a review
    obligation.

Usage: ``python verify_generated_project.py [project_root]`` (the root is
auto-detected when this file sits in a generated project's ``tools/``).
"""

import hashlib
import json
import os
import re
import sys
from typing import Dict, List, NoReturn, Set

# Role prefixes that decorate an exact NRPy gridfunction name (section 5.2).
# The scan must see `in_alpha` as an occurrence of `alpha`, otherwise it is
# blind to every spelling the Dendro backend actually emits.
_ROLE_PREFIXES = ("in_", "rhs_", "out_", "diag_", "aux_")

# NRPy's loop helper emits its index variables as `i0`/`i1`/`i2` (point loops)
# and `blk_id` (block loops); a numerical loop in a fixed template is a section
# 9.16 violation whatever type it declares.  Host lifecycle and reduction loops
# owned by Dendro (section 14.6) use their own index names and are permitted.
_NRPY_LOOP_RE = re.compile(
    r"for\s*\(\s*(?:const\s+)?(?:int|unsigned|long|std::size_t|std::ptrdiff_t|auto)?"
    r"\s*(?:i[012]|blk_id)\s*="
)


def _nrpy_loop_present(text: str) -> bool:
    """
    Return True if text contains an NRPy-emitted numerical loop.

    Section 9.16 forbids an NRPy point or block loop in a fixed template.  The
    NRPy loop helper spells its index variables ``i0``/``i1``/``i2`` and
    ``blk_id``; host lifecycle and reduction loops own their index names and
    are permitted.

    :param text: Text to search.
    :return: True if an NRPy numerical loop appears.

    Doctests:
    >>> _nrpy_loop_present("for (int i0 = 0; i0 < n; i0++) {")
    True
    >>> _nrpy_loop_present("for (int blk_id = 0; blk_id < nb; blk_id++) {")
    True

    A host-owned reduction loop with its own index name is not a violation.

    >>> _nrpy_loop_present("for (unsigned bx = pad; bx < nx - pad; ++bx) {")
    False
    >>> _nrpy_loop_present("for (const VariableRef::Group group : GROUPS) {")
    False
    """
    return _NRPY_LOOP_RE.search(text) is not None


# Files that cannot list their own digest (section 12.5).
_HASH_ARTIFACTS = (
    "GENERATION_RECEIPT.json",
    "manifest/file_hashes.json",
    "manifest/SHA256SUMS",
)


def fail(message: str) -> NoReturn:
    """
    Print a failure message and exit nonzero.

    :param message: The failure description.
    """
    print(f"FAIL: {message}")
    sys.exit(1)


def _word_present(text: str, name: str) -> bool:
    """
    Return True if a name occurs as a standalone or role-prefixed identifier.

    :param text: Text to search.
    :param name: Exact registered gridfunction name.
    :return: True if the bare or role-prefixed name occurs as a whole token.

    Doctests:
    >>> _word_present("const T alpha = 1;", "alpha")
    True

    Whitepaper section 5.2: the scan must see every role-prefixed spelling the
    backend emits, or it is blind to the form scaffolds actually contain.

    >>> [_word_present(f"{p}alpha[pp] = 1;", "alpha")
    ...  for p in ("in_", "rhs_", "out_", "diag_", "aux_")]
    [True, True, True, True, True]

    It must not fire on a longer identifier that merely contains the name.

    >>> _word_present("alphabet_size", "alpha")
    False
    >>> _word_present("my_alpha_helper", "alpha")
    False
    """
    prefixes = "|".join(re.escape(prefix) for prefix in _ROLE_PREFIXES)
    pattern = (
        r"(?<![A-Za-z0-9_])(?:"
        + prefixes
        + r")?"
        + re.escape(name)
        + r"(?![A-Za-z0-9_])"
    )
    return re.search(pattern, text) is not None


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


def _relative_files(root: str) -> Set[str]:
    """
    Return every file under a directory, as project-relative slash paths.

    :param root: Directory to walk.
    :return: Set of relative paths.
    """
    found: Set[str] = set()
    for dirpath, dirnames, filenames in os.walk(root):
        # Interpreter bytecode is not a generated artifact: running the shipped
        # tools creates tools/__pycache__/, and the inventory check must not
        # then reject the very tree those tools were run from.
        dirnames[:] = [name for name in dirnames if name != "__pycache__"]
        for filename in filenames:
            if filename.endswith((".pyc", ".pyo")):
                continue
            path = os.path.join(dirpath, filename)
            found.add(os.path.relpath(path, root).replace(os.sep, "/"))
    return found


def _check_hashes(root: str) -> Dict[str, str]:
    """
    Verify the hash manifests and the on-disk inventory.

    :param root: Generated project root.
    :return: The manifest's path-to-digest mapping.
    """
    manifest_path = os.path.join(root, "manifest", "file_hashes.json")
    if not os.path.isfile(manifest_path):
        fail(
            f"{manifest_path} not found: pass the generated project root "
            "(the directory containing manifest/ and Dendro-GR/)"
        )
    file_hashes: Dict[str, str] = json.loads(read(manifest_path))["files"]
    if not file_hashes:
        fail("file_hashes.json lists no files")
    for rel, expected in sorted(file_hashes.items()):
        path = os.path.join(root, rel)
        if not os.path.isfile(path):
            fail(f"listed file missing: {rel}")
        actual = sha256_of(path)
        if actual != expected:
            fail(f"hash mismatch for {rel}: {actual} != {expected}")
    # Inventory completeness: an unlisted file would otherwise be installed by
    # the copy tool and never verified (sections 11.5/12.5).
    for rel in sorted(_relative_files(root) - set(file_hashes) - set(_HASH_ARTIFACTS)):
        fail(f"file on disk is not listed in file_hashes.json: {rel}")
    sums = read(os.path.join(root, "manifest", "SHA256SUMS"))
    listed_in_sums: Set[str] = set()
    for line in sums.splitlines():
        if not line.strip():
            continue
        digest, rel = line.split("  ", 1)
        listed_in_sums.add(rel)
        if file_hashes.get(rel) != digest:
            fail(f"SHA256SUMS disagrees with file_hashes.json for {rel}")
    if listed_in_sums != set(file_hashes):
        fail("SHA256SUMS and file_hashes.json cover different file sets")
    return file_hashes


def _check_banners(root: str, abi: str) -> None:
    """
    Require the section 12.7 ownership banner on every emitted module file.

    Exempt: JSON (no comment syntax), the markdown READMEs, and the NRPy-owned
    mock-vehicle host header (a test double, not a generated artifact).

    :param root: Generated project root.
    :param abi: Module ABI hash from module.json.
    """
    module_dir = os.path.join(root, "Dendro-GR", "FCCZ4_GR")
    for dirpath, _dirnames, filenames in os.walk(module_dir):
        for filename in filenames:
            path = os.path.join(dirpath, filename)
            if filename.endswith((".json", ".md")):
                continue
            relative = os.path.relpath(path, module_dir).replace(os.sep, "/")
            if relative == "host_mock/dendro_mock.hpp":
                continue
            text = read(path)
            is_hash_comment = (
                filename.endswith(".cmake")
                or filename == "CMakeLists.txt"
                or filename.endswith(".toml")
            )
            expected_banner = (
                "# GENERATED FILE - DO NOT EDIT"
                if is_hash_comment
                else "// GENERATED FILE - DO NOT EDIT"
            )
            if not text.startswith(expected_banner):
                fail(f"missing ownership banner: {os.path.relpath(path, root)}")
            if "Module ABI: " + abi not in text:
                fail(f"missing module ABI in banner: {os.path.relpath(path, root)}")


def _check_cmake_source_list(root: str) -> None:
    """
    Require set equality between the CMake sources and the frozen CFunctions.

    :param root: Generated project root.
    """
    module_dir = os.path.join(root, "Dendro-GR", "FCCZ4_GR")
    cmake_text = read(
        os.path.join(module_dir, "generated", "cmake", "nrpy_generated_sources.cmake")
    )
    listed = {
        line.strip().replace("${FCCZ4_MODULE_ROOT}/", "")
        for line in cmake_text.splitlines()
        if line.strip().startswith("$")
    }
    if not listed:
        fail("nrpy_generated_sources.cmake lists no sources")
    cfunctions = json.loads(read(os.path.join(root, "manifest", "cfunctions.json")))
    expected = {
        f"{record['subdirectory']}/{record['name']}.cpp".replace("./", "")
        for record in cfunctions["records"]
    }
    if listed != expected:
        fail(
            "CMake source list is not the frozen CFunction source set: "
            f"only in CMake={sorted(listed - expected)} "
            f"only in manifest={sorted(expected - listed)}"
        )
    for rel in sorted(listed):
        if not os.path.isfile(os.path.join(module_dir, rel)):
            fail(f"CMake source list references missing file: {rel}")


def _check_scaffolds(root: str) -> None:
    """
    Require scaffold files to carry no field name and no NRPy numerical loop.

    The scan runs over the delivered project's own non-generated files, so it
    works in an installed tree where the NRPy sources are absent (section 9.16).

    :param root: Generated project root.
    """
    state = json.loads(read(os.path.join(root, "manifest", "state.json")))
    evol_names: List[str] = [
        record["name"] for record in state["records"] if record["group"] == "EVOL"
    ]
    if not evol_names:
        fail("state.json has no EVOL records")
    module = os.path.join(root, "Dendro-GR", "FCCZ4_GR")
    scaffolds = [
        os.path.join(root, "README.md"),
        os.path.join(module, "README.md"),
        os.path.join(module, "CMakeLists.txt"),
        os.path.join(module, "tests", "CMakeLists.txt"),
        os.path.join(module, "src", "fccz4_main.cpp"),
        os.path.join(module, "include", "fccz4Ctx.h"),
        os.path.join(module, "src", "fccz4Ctx.cpp"),
        os.path.join(module, "generated", "include", "generated_project_preamble.hpp"),
    ]
    for dirpath, _dirnames, filenames in os.walk(os.path.join(module, "pars")):
        for filename in filenames:
            scaffolds.append(os.path.join(dirpath, filename))
    for path in scaffolds:
        if not os.path.isfile(path):
            continue
        text = read(path)
        relative = os.path.relpath(path, root)
        for name in evol_names:
            if _word_present(text, name):
                fail(
                    f"concrete EVOL name {name!r} in scaffold {relative} (section 9.16)"
                )
        if _nrpy_loop_present(text):
            match = _NRPY_LOOP_RE.search(text)
            assert match is not None
            fail(f"NRPy numerical loop in scaffold {relative}: {match.group(0)!r}")


def verify(project_root: str) -> None:
    """
    Run every verifier check on one generated project.

    :param project_root: Path of the generated project root.
    """
    root = os.path.abspath(project_root)
    if not os.path.isdir(root):
        fail(f"project root {root} is not a directory")
    file_hashes = _check_hashes(root)
    module_json = json.loads(read(os.path.join(root, "manifest", "module.json")))
    abi = module_json["module_abi_hash"]
    if len(abi) != 64:
        fail("module_abi_hash is not a SHA-256 digest")
    _check_banners(root, abi)
    for dirpath, _dirnames, filenames in os.walk(os.path.join(root, "manifest")):
        for filename in filenames:
            if not filename.endswith(".json"):
                continue
            try:
                json.loads(read(os.path.join(dirpath, filename)))
            except ValueError as exc:
                fail(f"manifest {filename} is not valid JSON: {exc}")
    _check_cmake_source_list(root)
    _check_scaffolds(root)
    print(f"OK: verified {len(file_hashes)} files, ABI {abi[:12]}...")


def _auto_detect_root() -> str:
    """
    Return the project root containing this tool.

    :return: The generated project root.
    """
    candidate = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if os.path.isfile(os.path.join(candidate, "manifest", "file_hashes.json")):
        return candidate
    fail("could not auto-detect the project root; pass it explicitly")


def main() -> None:
    """Verify the project root given on the command line, or auto-detect it."""
    if len(sys.argv) == 2:
        verify(sys.argv[1])
    elif len(sys.argv) == 1:
        verify(_auto_detect_root())
    else:
        fail("usage: verify_generated_project.py [project_root]")


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
