# nrpy/infrastructures/Dendro/project.py
"""
Pure project exporter for the generated Dendro fCCZ4 project.

Turns the frozen NRPy snapshot into the complete generated project
(whitepaper sections 9.15/12): state/parameter/CFunction artifacts, CMake,
manifests, templates, tools, and the receipt.  The exporter accepts no
formulation configuration: every scientific choice has already been
registered and frozen.  It is read-only over the registries and verifies
them before and after the staged commit (section 4.3).

Staged write: artifacts are written to a temporary sibling directory,
verified against the snapshot hashes, and then atomically replaced over the
target project directory (a non-generated existing target is refused unless
it was produced by a previous run of this exporter - detected by its
GENERATION_RECEIPT.json).

Template policy (whitepaper section 9.16).  The fixed templates this exporter
renders from ``templates/`` may contain only structures that NRPy core
registries do not model naturally: CRTP context scaffolding, Dendro vector
lifecycle calls, CMake target declarations, executable ``main()`` structure and
host adapter functions.  They must not contain a concrete state-field name, a
physics parameter default, an fCCZ4 expression, an FD coefficient, a numerical
point or block loop, or a source inventory; ``$``-style placeholders are the
only profile-dependent content.  Section 9.16 grants one exception: the
generated-project self-test template is a designated test artifact, so it may
name evolved fields and carry its own loops over a test block.  Every other
template is scanned for both by
:func:`nrpy.infrastructures.Dendro.validation._check_scaffolds`.
"""

import hashlib
import importlib
import json
import os
import shutil
import tempfile
from typing import Any, Dict, Tuple, Union, cast

from nrpy.infrastructures.Dendro import (
    CFunction_output,
    CodeParameters_output,
    cmake,
    gridfunction_output,
    manifests,
    validation,
)
from nrpy.infrastructures.Dendro.freeze import (
    FrozenNRPyDendroSnapshot,
    assert_mutable_registries_match,
)

# Module-relative paths of the host-adapter (template) files.  These are the
# only non-pure-rendered files; their content comes from fixed templates
# (section 9.16) plus the host-vehicle preamble.
_HOST_FILES = (
    "include/fccz4Ctx.h",
    "src/fccz4Ctx.cpp",
    "src/fccz4_main.cpp",
    "pars/fccz4_minkowski.toml",
    "tests/test_generated.cpp",
    "tests/CMakeLists.txt",
    "README.md",
    "host_mock/dendro_mock.hpp",
)


def _template_text(template_name: str) -> str:
    """
    Read a fixed template from the Dendro templates directory.

    :param template_name: Template file name (templates directory).
    :return: The template text.
    """
    path = os.path.join(os.path.dirname(__file__), "templates", template_name)
    with open(path, "r", encoding="utf-8") as handle:
        return handle.read()


def _render_template(template_name: str, mapping: Dict[str, str]) -> str:
    """
    Render a fixed template with simple $KEY placeholders.

    Only exact $KEY tokens for keys present in the mapping are substituted;
    all other text (including CMake ${VARIABLE} references, which use
    braces) is left untouched, so CMake templates render safely.

    :param template_name: Template file name (templates directory).
    :param mapping: Placeholder mapping.
    :return: The rendered text.
    """
    text = _template_text(template_name)
    for key, value in mapping.items():
        text = text.replace("$" + key, value)
    return text


def _module_abi(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Return the module ABI hash (banner placeholder).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: The module ABI hash string.
    """
    return snapshot.hashes.module_abi_hash


def _host_mapping(snapshot: FrozenNRPyDendroSnapshot) -> Dict[str, str]:
    """
    Build the template placeholder mapping for the host-adapter files.

    The RHS call tails are derived from the frozen CFunction records
    (section 6.2 used-parameter closure), so no parameter name is hardcoded
    in a template.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: The placeholder mapping.
    """
    gen = dict(snapshot.generation_parameters)
    used = [cp.name for cp in snapshot.codeparameters]
    # The RHS/ID CFunctions take the used CodeParameters by value.  The
    # caller supplies the values from the generated parameter table (the
    # Ctx member ``fccz4_params``, or a local ``params`` in the self-test),
    # so no parameter default is hardcoded in a template.
    rhs_tail = (
        (", " + ", ".join(f"fccz4_params.{name}" for name in used)) if used else ""
    )
    block_tail = (", " + ", ".join(f"params.{name}" for name in used)) if used else ""
    padding = list(snapshot.required_padding)
    return {
        "MODULE_ABI_HASH": _module_abi(snapshot),
        "TARGET_SCALAR": str(gen.get("fp_type", "double")),
        "NRPY_ARITHMETIC": str(gen.get("fp_type", "double")),
        "RHS_CPARAM_TAIL": rhs_tail,
        "RHS_BLOCK_CPARAM_TAIL": block_tail,
        "FD_ORDER": str(int(str(gen.get("fd_order", 0)))),
        "PADDING": str(list(padding)),
        "KO_ENABLED": "true" if bool(gen.get("Dendro_enable_KO", False)) else "false",
        "FCCZ4_PARAMS_TOML": _inert_parameter_toml(snapshot),
    }


def _inert_parameter_toml(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the sample runtime-parameter table, commented out.

    This generated profile has no TOML binding (the registered
    ``fccz4_params_parse_toml`` refuses a supplied file), so shipping a live
    table would invite edits that silently do nothing.  The table is emitted
    commented, with the reason, so the values remain documented.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: The commented TOML block.
    """
    sample = CodeParameters_output.render_parameter_toml_sample(snapshot)
    commented = "\n".join(
        f"# {line}" if line.strip() else "#"
        for line in sample.rstrip("\n").splitlines()
    )
    return (
        "# The runtime parameter table below is generated from the frozen\n"
        "# CodeParameter records and is shown for reference only: this profile\n"
        "# has no TOML binding yet, so fccz4Solver refuses a -t argument rather\n"
        "# than appear to apply these values (section 6.4).  The effective\n"
        "# values are the generated defaults, printed at startup.\n" + commented
    )


def _banner_py(source_object: str, abi: str) -> str:
    """
    Render the Python generated-file banner (section 12.7).

    :param source_object: Frozen source the file was derived from.
    :param abi: Module ABI hash.
    :return: Banner lines (comment-prefixed).
    """
    return (
        "# GENERATED FILE - DO NOT EDIT\n"
        "# Generator: nrpy.infrastructures.Dendro\n"
        "# Regenerate with: python -m nrpy.examples.dendro_fccz4\n"
        f"# Source object: {source_object}\n"
        f"# Module ABI: {abi}\n"
    )


def _render_tool_py(source_module: str, source_object: str, abi: str) -> str:
    """
    Render a generated project tool from its NRPy single-source module.

    The rendered file is the NRPy source with the generated-file banner as the
    first comment block (section 12.7: every emitted file *begins* with the
    origin banner), so the project tool and the NRPy module cannot drift
    (sections 9.15/12.6).  The NRPy module ends with the repository-standard
    ``__main__`` doctest block, which the CI static-analysis step executes;
    the rendered tool is an executable instead, so that block is replaced by a
    call to ``main()``.

    :param source_module: Dotted module name of the NRPy source module.
    :param source_object: Source object line for the banner.
    :param abi: Module ABI hash.
    :return: The rendered tool text.
    :raises ValueError: If the source module has no ``__main__`` block.
    """
    module = importlib.import_module(source_module)
    path = str(module.__file__)
    with open(path, "r", encoding="utf-8") as handle:
        source = handle.read()
    marker = '\nif __name__ == "__main__":\n'
    if marker not in source:
        raise ValueError(
            f"{source_module} has no '__main__' block to replace; the rendered "
            "project tool would not be executable."
        )
    body = source[: source.index(marker)].rstrip("\n")
    lines = body.splitlines(keepends=True)
    insert_at = 1 if lines and lines[0].startswith("#!") else 0
    banner = _banner_py(source_object, abi).splitlines(keepends=True)
    rendered = lines[:insert_at] + banner + lines[insert_at:]
    entry_point = '\n\n\nif __name__ == "__main__":\n    main()\n'
    return "".join(rendered) + entry_point


def _render_tool_wrapper(abi: str, *, docstring: str, forced_flag: str) -> str:
    """
    Render a single-action wrapper around the installer.

    Each wrapper forces exactly the action its file name promises.  Re-exposing
    the installer's own parser would make ``remove_generated_module.py
    --execute`` *install* the module -- a mutating action behind a read-only
    name (section 12.6).

    :param abi: Module ABI hash for the banner.
    :param docstring: One-line module docstring for the wrapper.
    :param forced_flag: The installer flag this wrapper always applies.
    :return: The rendered wrapper text.
    """
    return (
        _banner_py(
            "fixed tool rendered from nrpy.infrastructures.Dendro.copy_tool",
            abi,
        )
        + f'"""{docstring}"""\n'
        + "import sys\n\n"
        + "# Importing the installer must not leave a __pycache__ directory in\n"
        + "# tools/: the generated project's own verifier checks that every file\n"
        + "# on disk is a listed artifact.\n"
        + "sys.dont_write_bytecode = True\n\n"
        + "from copy_into_dendro import main  # noqa: E402\n\n"
        + f"FORCED_FLAG = {forced_flag!r}\n\n\n"
        + 'if __name__ == "__main__":\n'
        + "    if FORCED_FLAG not in sys.argv:\n"
        + "        sys.argv.insert(1, FORCED_FLAG)\n"
        + "    main()\n"
    )


def _capabilities_text(_snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Return the canonical dendrolib capabilities record text.

    The record is copied verbatim from the NRPy Dendro infrastructure
    (honest UNPINNED/UNPROVEN state while Dendrolib is not vendored).

    :param _snapshot: The frozen NRPy Dendro snapshot (unused; the record is
        the repository-pinned capability file).
    :return: Canonical JSON text.
    """
    path = os.path.join(os.path.dirname(__file__), "dendrolib_capabilities.json")
    with open(path, "r", encoding="utf-8") as handle:
        return handle.read()


def _render_artifacts(snapshot: FrozenNRPyDendroSnapshot) -> Dict[str, str]:
    """
    Render every generated project artifact keyed by project-relative path.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Mapping of project-relative path to file text.
    """
    abi = _module_abi(snapshot)
    artifacts: Dict[str, str] = {}

    # Module: generated include sources (pure renders).
    module_prefix = "Dendro-GR/FCCZ4_GR/"
    gen = dict(snapshot.generation_parameters)
    artifacts[module_prefix + "generated/include/fccz4_types.hpp"] = _render_template(
        "fccz4_types_hpp.in",
        {
            "MODULE_ABI_HASH": abi,
            "TARGET_SCALAR": str(gen.get("fp_type", "double")),
            "NRPY_ARITHMETIC": str(gen.get("fp_type", "double")),
        },
    )
    artifacts[module_prefix + "generated/include/fccz4_state.hpp"] = (
        gridfunction_output.render_fccz4_state_hpp(snapshot)
    )
    artifacts[module_prefix + "generated/include/fccz4_parameters.hpp"] = (
        CodeParameters_output.render_fccz4_parameters_hpp(snapshot)
    )
    artifacts[module_prefix + "generated/include/fccz4_provenance.hpp"] = (
        gridfunction_output.render_fccz4_provenance_hpp(snapshot)
    )
    artifacts |= CFunction_output.CFunction_artifacts(snapshot)
    artifacts[module_prefix + "generated/cmake/nrpy_generated_sources.cmake"] = (
        cmake.render_nrpy_generated_sources_cmake(snapshot)
    )
    # Generated project self-test preamble (host-vehicle include).
    artifacts[module_prefix + "generated/include/generated_project_preamble.hpp"] = (
        _render_template("generated_project_preamble_cpp.in", {"MODULE_ABI_HASH": abi})
    )

    # Module: host adapter files (fixed templates + host mock header).
    host_mapping = _host_mapping(snapshot)
    artifacts[module_prefix + "include/fccz4Ctx.h"] = _render_template(
        "fccz4Ctx_h.in", host_mapping
    )
    artifacts[module_prefix + "src/fccz4Ctx.cpp"] = _render_template(
        "fccz4Ctx_cpp.in", host_mapping
    )
    artifacts[module_prefix + "src/fccz4_main.cpp"] = _render_template(
        "fccz4_main_cpp.in", host_mapping
    )
    artifacts[module_prefix + "pars/fccz4_minkowski.toml"] = _render_template(
        "fccz4_minkowski_toml.in", host_mapping
    )
    artifacts[module_prefix + "tests/test_generated.cpp"] = _render_template(
        "generated_project_tests_cpp.in", host_mapping
    )
    artifacts[module_prefix + "tests/CMakeLists.txt"] = _render_template(
        "FCCZ4_GR_tests_CMakeLists.in", host_mapping
    )
    artifacts[module_prefix + "README.md"] = _render_template(
        "FCCZ4_GR_README_md.in", host_mapping
    )
    # Mock host header: copied verbatim from the pinned host-vehicle API that
    # this infrastructure owns (the real Dendrolib header arrives with the
    # I0-1 gates).  It is a source asset of the Dendro package, not a test
    # fixture, so it lives beside the exporter that ships it.
    mock_path = os.path.join(
        os.path.dirname(__file__),
        "host_mock",
        "dendro_mock.hpp",
    )
    with open(os.path.abspath(mock_path), "r", encoding="utf-8") as handle:
        artifacts[module_prefix + "host_mock/dendro_mock.hpp"] = handle.read()

    # Module CMake (fixed template, source list from the frozen CFunctions).
    artifacts[module_prefix + "CMakeLists.txt"] = cmake.render_module_cmake(snapshot)

    # Manifests.
    artifacts |= manifests.all_manifests(snapshot, _capabilities_text(snapshot))

    # Project README.
    artifacts["README.md"] = _render_template("project_README_md.in", host_mapping)

    # Tools: rendered from the NRPy single-source modules.
    artifacts["tools/verify_generated_project.py"] = _render_tool_py(
        "nrpy.infrastructures.Dendro.validation",
        "fixed tool rendered from nrpy.infrastructures.Dendro.validation",
        abi,
    )
    artifacts["tools/copy_into_dendro.py"] = _render_tool_py(
        "nrpy.infrastructures.Dendro.copy_tool",
        "fixed tool rendered from nrpy.infrastructures.Dendro.copy_tool",
        abi,
    )
    # verify_copy.py and remove_generated_module.py are thin wrappers over
    # the installer (single source, no duplicated logic).
    artifacts["tools/verify_copy.py"] = _render_tool_wrapper(
        abi,
        docstring="Verify an installed FCCZ4_GR copy against its receipt.",
        forced_flag="--verify-only",
    )
    artifacts["tools/remove_generated_module.py"] = _render_tool_wrapper(
        abi,
        docstring="Remove the installed FCCZ4_GR module (receipt-identified).",
        forced_flag="--remove",
    )

    return artifacts


def _write_artifacts(stage: str, artifacts: Dict[str, str]) -> None:
    """
    Write every artifact to the staged directory (LF, UTF-8).

    :param stage: Staged project directory.
    :param artifacts: Mapping of project-relative path to file text.
    """
    for rel, text in sorted(artifacts.items()):
        path = os.path.join(stage, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)


def _sha256_of_file(path: str) -> str:
    """
    Return the SHA-256 hex digest of one file.

    :param path: File to digest.
    :return: Hex digest.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_hashes(stage: str, exclude: Tuple[str, ...] = ()) -> Dict[str, str]:
    """
    Compute the SHA-256 of every staged file (project-relative paths).

    :param stage: Staged project directory.
    :param exclude: Project-relative paths to omit (the hash artifacts
        themselves, which cannot contain their own digests).
    :return: Mapping of project-relative path to SHA-256 hex digest.
    """
    digests: Dict[str, str] = {}
    for dirpath, _dirnames, filenames in os.walk(stage):
        for filename in sorted(filenames):
            path = os.path.join(dirpath, filename)
            rel = os.path.relpath(path, stage).replace(os.sep, "/")
            if rel in exclude:
                continue
            digests[rel] = _sha256_of_file(path)
    return digests


def verify_generated_project(
    project_dir: str, snapshot: FrozenNRPyDendroSnapshot
) -> None:
    """
    Run the verifier checks on a generated project directory.

    :param project_dir: The generated project root.
    :param snapshot: The frozen NRPy Dendro snapshot (for the receipt ABI).
    """
    validation.verify(project_dir)
    receipt_path = os.path.join(project_dir, "GENERATION_RECEIPT.json")
    with open(receipt_path, "r", encoding="utf-8") as handle:
        receipt = json.loads(handle.read())
    if receipt["module_abi_hash"] != _module_abi(snapshot):
        validation.fail("receipt module ABI differs from the snapshot")


def _doctest_frozen_profile(freeze: Any) -> FrozenNRPyDendroSnapshot:
    """
    Register and freeze the qualified profile, for this module's doctest.

    The exporter's contract is only meaningful against a real frozen snapshot,
    so the doctest builds the same profile the shipped example registers.

    :param freeze: The freeze entry point (injected so the doctest keeps its
        imports local).
    :return: The frozen snapshot.
    """
    # pylint: disable=import-outside-toplevel
    # These imports are deliberately deferred.  This module is the *pure*
    # exporter: it renders a frozen snapshot and must not depend on the
    # formulation builders at import time, or every export would drag the
    # fCCZ4 equation modules in and the layering would invert.  Only this
    # doctest fixture needs them, so they are local to it.
    import nrpy.params as par

    from nrpy.infrastructures.Dendro import registration as reg
    from nrpy.infrastructures.Dendro.general_relativity import (
        diagnostics,
        initial_data,
        projection,
        rhs_eval,
    )
    from nrpy.infrastructures.Dendro.runtime import parameters

    par.set_parval_from_str("Infrastructure", "Dendro")
    par.set_parval_from_str("fp_type", "double")
    par.set_parval_from_str("parallelization", "none")
    par.set_parval_from_str("fd_order", 4)
    par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    par.set_parval_from_str("Dendro_enable_KO", False)
    build = rhs_eval.build_fccz4_rhs(fd_order=4, enable_KO=False)
    for role, name, params, body, calls in (
        (
            "rhs_block",
            rhs_eval.BLOCK_CFUNCTION,
            build.block_params,
            build.block_body,
            (),
        ),
        (
            "rhs",
            rhs_eval.GLOBAL_CFUNCTION,
            build.global_params,
            build.global_body,
            (rhs_eval.BLOCK_CFUNCTION,),
        ),
        (
            "rhs_flat_block",
            rhs_eval.FLAT_BLOCK_CFUNCTION,
            build.flat_block_params,
            build.flat_block_body,
            (rhs_eval.BLOCK_CFUNCTION,),
        ),
    ):
        reg.register_dendro_CFunction(
            role=role,
            entry_point=True,
            calls=calls,
            name=name,
            desc=f"Doctest profile: {role}.",
            subdirectory="generated/src/rhs",
            params=params,
            body=body,
        )
    par.glb_extras_dict.setdefault("Dendro", {})["builder_extras"] = {
        "operator_manifest": build.operator_manifest,
        "rhs_canonical": dict(build.canonical_expression_digests),
        "provenance": {
            "rhs_symbols": list(build.rhs_symbols),
            "upwind_control_fields": list(build.upwind_control_fields),
            "used_codeparameters": list(build.used_codeparameters),
            "padding": list(build.padding),
        },
    }
    initial_data.register_minkowski_CFunctions()
    initial_data.register_perturbation_CFunctions()
    initial_data.register_initial_data_conversion_CFunctions()
    projection.register_projection_CFunctions()
    diagnostics.register_diagnostics_CFunctions()
    parameters.register_parameter_CFunctions_last()
    return cast(FrozenNRPyDendroSnapshot, freeze(profile_name="doctest_profile"))


def output_project(
    *,
    snapshot: FrozenNRPyDendroSnapshot,
    project_dir: Union[str, "os.PathLike[str]"],
) -> str:
    """
    Export the complete generated project (section 9.15).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :param project_dir: Target project directory (Path or str).
    :return: The project directory path (str).
    :raises ValueError: On drift, unsafe targets, or a non-generated
        existing target.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> import json
    >>> import tempfile
    >>> import nrpy.params as par
    >>> from nrpy.infrastructures.Dendro import validation as _val
    >>> from nrpy.infrastructures.Dendro.freeze import (
    ...     freeze_nrpy_dendro_environment,
    ... )
    >>> from nrpy.infrastructures.Dendro import project as _proj
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _snapshot = _doctest_frozen_profile(freeze_nrpy_dendro_environment)
    >>> _root = tempfile.mkdtemp(prefix="dendro_export_doctest_")
    >>> _out = os.path.join(_root, "project")
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _ = _proj.output_project(snapshot=_snapshot, project_dir=_out)

    Section 12.5: the shipped verifier accepts the tree the exporter just
    wrote -- hashes, banners, manifests and the derived CMake source list.

    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _val.verify(_out)

    Sections 11.5/12.5, inventory equality: every file on disk is hash-listed,
    apart from the three artifacts that cannot list their own digest.

    >>> with open(os.path.join(_out, "manifest", "file_hashes.json")) as _fh:
    ...     _listed = set(json.load(_fh)["files"])
    >>> _on_disk = _val._relative_files(_out)
    >>> sorted(_on_disk - _listed - set(_val._HASH_ARTIFACTS))
    []
    >>> sorted(_listed - _on_disk)
    []

    Section 9.16: no scaffold names a registered field or carries an
    NRPy-emitted numerical loop.

    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _val._check_scaffolds(_out)
    >>> shutil.rmtree(_root)
    """
    assert_mutable_registries_match(snapshot)
    project_dir = os.path.abspath(str(project_dir))
    here = os.path.abspath(__file__)
    nrpy_root = os.path.dirname(os.path.dirname(os.path.dirname(here)))
    dendro_root = os.path.dirname(os.path.dirname(here))
    if project_dir in (
        os.path.abspath("/"),
        os.path.abspath(os.path.expanduser("~")),
        nrpy_root,
        dendro_root,
    ):
        raise ValueError(f"Refusing to generate into unsafe path {project_dir}")
    if os.path.isdir(project_dir) and not os.path.isfile(
        os.path.join(project_dir, "GENERATION_RECEIPT.json")
    ):
        raise ValueError(
            f"Refusing to generate into existing non-generated tree " f"{project_dir}"
        )

    # The three hash artifacts cover every other project file; they cannot
    # contain their own digests (section 12.5: no self-referential hashes).
    hash_artifacts = (
        "GENERATION_RECEIPT.json",
        "manifest/file_hashes.json",
        "manifest/SHA256SUMS",
    )
    parent = os.path.dirname(project_dir)
    os.makedirs(parent, exist_ok=True)
    stage = tempfile.mkdtemp(prefix=".dendro_project_", dir=parent)
    try:
        artifacts = _render_artifacts(snapshot)
        _write_artifacts(stage, artifacts)
        file_hashes = _file_hashes(stage, exclude=hash_artifacts)
        receipt = manifests.render_receipt_json(snapshot, file_hashes)
        file_hashes_text = manifests.render_file_hashes_json(file_hashes)
        sums_text = manifests.render_SHA256SUMS(file_hashes)
        for rel, text in (
            ("GENERATION_RECEIPT.json", receipt),
            ("manifest/file_hashes.json", file_hashes_text),
            ("manifest/SHA256SUMS", sums_text),
        ):
            path = os.path.join(stage, rel)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w", encoding="utf-8", newline="\n") as handle:
                handle.write(text)
        # Verify the staged project against the snapshot before commit.
        verify_generated_project(stage, snapshot)
        # Re-render in memory and require a byte-identical second render
        # (section 12.5: re-render in memory and verify hashes before
        # committing the staged directory).
        second = _render_artifacts(snapshot)
        for rel, text in sorted(second.items()):
            staged_text = open(os.path.join(stage, rel), "r", encoding="utf-8").read()
            if staged_text != text:
                raise ValueError(
                    f"Staged file {rel} differs from a fresh in-memory render"
                )
        # Atomic commit: replace the existing project directory (if any).
        if os.path.isdir(project_dir):
            shutil.rmtree(project_dir)
        os.replace(stage, project_dir)
        stage = ""
    finally:
        if stage and os.path.isdir(stage):
            shutil.rmtree(stage, ignore_errors=True)
    assert_mutable_registries_match(snapshot)
    return project_dir


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
