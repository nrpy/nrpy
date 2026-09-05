# nrpy/infrastructures/Dendro/manifests.py
"""
Canonical manifest renderers for the generated Dendro project.

Every manifest is pure JSON rendered from the frozen snapshot (whitepaper
sections 12.5/15): sorted keys, stable arrays, exact rational coefficient
strings, UTF-8/LF text, no absolute source paths, and no timestamps in
hashed content.  The operator records in ``stencils.json`` come from the
builder sidecar (``snapshot.extras``), whose coefficients were computed once
by the builder from NRPy's ``compute_fdcoeffs_fdstencl`` (section 8.6: one
coefficient source, not a second stencil generator).
"""

import hashlib
import json
import re
import sys
from typing import Any, Dict, List, Tuple, cast

import sympy

import nrpy
from nrpy.infrastructures.Dendro import CodeParameters_output, gridfunction_output
from nrpy.infrastructures.Dendro.freeze import FrozenNRPyDendroSnapshot


def _gen_lookup(snapshot: FrozenNRPyDendroSnapshot) -> Dict[str, object]:
    """
    Return the frozen generation parameters as a dict.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Mapping of generation parameter name to frozen value.
    """
    return dict(snapshot.generation_parameters)


def _padding_constants(snapshot: FrozenNRPyDendroSnapshot) -> Dict[str, Any]:
    """
    Derive the generated padding constants (section 8.6).

    The required padding per axis comes from the snapshot reduction of the
    captured accesses; the minimum element order under the current Dendro
    rule (padding = elementOrder/2) follows as 2*max(padding).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: JSON-ready mapping of padding constants.
    """
    padding = list(snapshot.required_padding)
    return {
        "required_padding": padding,
        "minimum_element_order_under_current_rule": 2 * max(padding),
    }


def render_module_json(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the module.json manifest (section 15.1).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Canonical JSON text.
    """
    gen = _gen_lookup(snapshot)
    evol_count = sum(1 for fg in snapshot.gridfunctions if fg.group == "EVOL")
    payload = {
        "module_name": str(gen.get("Dendro_module_name", "FCCZ4_GR")),
        "infrastructure": "Dendro",
        "equation_system": "fCCZ4",
        "coordinate_system": str(gen.get("Dendro_fccz4_CoordSystem")),
        "conformal_factor": str(gen.get("EvolvedConformalFactor_cf", "chi")),
        "lapse": str(gen.get("Dendro_fccz4_LapseEvolutionOption")),
        "shift": str(gen.get("Dendro_fccz4_ShiftEvolutionOption")),
        "derivative_backend": str(gen.get("Dendro_derivative_backend", "full_stencil")),
        "finite_difference_order": int(str(gen.get("fd_order", 0))),
        "ko_owner": "nrpy_full_stencil",
        "ko_enabled": bool(gen.get("Dendro_enable_KO", False)),
        "evolved_variable_count": evol_count,
        "state_schema_hash": snapshot.hashes.state_schema_hash,
        "parameter_schema_hash": snapshot.hashes.parameter_schema_hash,
        "equation_hash": snapshot.hashes.equation_hash,
        "stencil_hash": snapshot.hashes.stencil_hash,
        "cfunction_api_hash": snapshot.hashes.cfunction_api_hash,
        "cfunction_source_hash": snapshot.hashes.cfunction_source_hash,
        "module_abi_hash": snapshot.hashes.module_abi_hash,
    }
    payload.update(_padding_constants(snapshot))
    return json.dumps(payload, sort_keys=True, indent=2) + "\n"


def render_generation_parameters_json(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the generation_parameters.json manifest.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Canonical JSON text.
    """
    return (
        json.dumps(
            {
                "profile_name": snapshot.profile_name,
                "parameters": dict(snapshot.generation_parameters),
            },
            sort_keys=True,
            indent=2,
        )
        + "\n"
    )


def render_stencils_json(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the stencils.json manifest from the builder sidecar.

    The operator records (name, FD order, exact rational coefficient
    strings, signed offsets, per-axis and total reach, consuming
    CFunctions, emitted fields, access count) are the single-derivation
    records produced by the PR5 operator manifest builder (section 15.4).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Canonical JSON text.
    :raises ValueError: If no operator manifest was recorded by the builder.
    """
    extras = dict(snapshot.extras)
    operator_manifest = extras.get("operator_manifest")
    if operator_manifest is None:
        raise ValueError(
            "No operator manifest recorded in the snapshot extras; the "
            "stencil records are single-derivation and cannot be re-derived "
            "by the pure exporter."
        )
    payload: Dict[str, Any] = dict(cast(Dict[str, Any], operator_manifest))
    payload["stencil_hash"] = snapshot.hashes.stencil_hash
    payload["access_padding"] = list(snapshot.required_padding)
    # The builder sidecar names only the CFunction whose kernel derived the
    # operator records.  Every frozen CFunction that captured a neighbour
    # access consumes the same operators at the same order, so the consumer
    # list is derived from the snapshot rather than inherited (section 15.4).
    payload["consuming_cfunctions"] = sorted(
        fc.name
        for fc in snapshot.cfunctions
        if any(access[1:] != (0, 0, 0) for access in fc.accesses)
    )
    return json.dumps(payload, sort_keys=True, indent=2) + "\n"


def render_cfunctions_json(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the cfunctions.json manifest (section 15.3).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Canonical JSON text.
    """
    records = []
    for fc in snapshot.cfunctions:
        accesses = sorted(set(fc.accesses))
        padding = [0, 0, 0]
        for _name, i0, i1, i2 in accesses:
            padding[0] = max(padding[0], abs(i0))
            padding[1] = max(padding[1], abs(i1))
            padding[2] = max(padding[2], abs(i2))
        records.append(
            {
                "name": fc.name,
                "subdirectory": fc.subdirectory,
                "role": fc.role,
                "entry_point": fc.entry_point,
                "lifecycle_hook": fc.lifecycle_hook,
                "calls": list(fc.calls),
                "used_codeparameters": list(fc.used_codeparameters),
                "required_padding": padding,
                "prototype_hash": hashlib.sha256(
                    fc.prototype.encode("utf-8")
                ).hexdigest(),
                "source_hash": hashlib.sha256(
                    fc.full_function.encode("utf-8")
                ).hexdigest(),
            }
        )
    payload = {
        "cfunction_api_hash": snapshot.hashes.cfunction_api_hash,
        "cfunction_source_hash": snapshot.hashes.cfunction_source_hash,
        "records": records,
    }
    return json.dumps(payload, sort_keys=True, indent=2) + "\n"


def render_provenance_json(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the provenance.json manifest (section 15.5).

    No timestamps appear in hashed content; the NRPy release string, the
    profile, and the builder provenance sidecar (when present) are recorded.
    Absolute paths are excluded.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Canonical JSON text.
    """
    extras = dict(snapshot.extras)
    payload = {
        "profile_name": snapshot.profile_name,
        "nrpy_release": getattr(nrpy, "__version__", "unknown"),
        "python": sys.version.split()[0],
        "sympy": sympy.__version__,
        "dendrolib": {
            "pinned_commit": None,
            "status": "UNPINNED - real Dendrolib gates remain open (section 13.2)",
        },
        "hashes": {
            "state_schema_hash": snapshot.hashes.state_schema_hash,
            "parameter_schema_hash": snapshot.hashes.parameter_schema_hash,
            "equation_hash": snapshot.hashes.equation_hash,
            "stencil_hash": snapshot.hashes.stencil_hash,
            "cfunction_api_hash": snapshot.hashes.cfunction_api_hash,
            "cfunction_source_hash": snapshot.hashes.cfunction_source_hash,
            "module_abi_hash": snapshot.hashes.module_abi_hash,
        },
    }
    if "provenance" in extras:
        payload["builder_provenance"] = extras["provenance"]
    return json.dumps(payload, sort_keys=True, indent=2) + "\n"


def _names_in_group(snapshot: FrozenNRPyDendroSnapshot, group: str) -> List[str]:
    """
    Return the frozen gridfunction names of one group, in registry order.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :param group: Group name (``EVOL``, ``AUXEVOL``, ``DIAG``, ``AUX``).
    :return: Ordered exact names.
    """
    return [fg.name for fg in snapshot.gridfunctions if fg.group == group]


def _written_fields(body: str, names: List[str], prefix: str) -> List[str]:
    """
    Return the names a frozen CFunction body assigns through a role pointer.

    The scan reads the frozen body, which is the authoritative registered
    text, so the manifest cannot claim a write the emitted source does not
    perform.

    :param body: The frozen ``full_function`` text.
    :param names: Candidate exact names, in registry order.
    :param prefix: Role prefix of the write target (``out_`` or ``diag_``).
    :return: The assigned names, in registry order.

    Doctests:
    >>> _written_fields("out_hDD00[pp] = x;", ["hDD00", "hDD01"], "out_")
    ['hDD00']

    A comparison is not an assignment, a longer name that merely starts with a
    candidate is not that candidate, and a pointer binding is not a write.

    >>> _written_fields("out_hDD00[pp] == x;", ["hDD00"], "out_")
    []
    >>> _written_fields("out_hDD001[pp] = x;", ["hDD00"], "out_")
    []
    >>> _written_fields("const T* const out_hDD00 = p;", ["hDD00"], "out_")
    []
    """
    # `re` rather than a plain-string method: the match must be anchored on an
    # identifier boundary (so `out_hDD00` is not found inside `out_hDD001`) and
    # must skip an arbitrary subscript before the `=`, neither of which
    # `.replace()` or `str.find` can express robustly.
    return [
        name
        for name in names
        if re.search(
            r"(?<![A-Za-z0-9_])" + re.escape(prefix + name) + r"\[[^\]]*\]\s*=(?!=)",
            body,
        )
        is not None
    ]


def _cfunction_records(
    snapshot: FrozenNRPyDendroSnapshot, roles: Tuple[str, ...]
) -> List[Dict[str, Any]]:
    """
    Return compact records for every frozen CFunction with one of the roles.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :param roles: Dendro roles to select.
    :return: JSON-ready records, in frozen (sorted name) order.
    """
    records = []
    for fc in snapshot.cfunctions:
        if fc.role not in roles:
            continue
        padding = [0, 0, 0]
        for _name, i0, i1, i2 in fc.accesses:
            padding[0] = max(padding[0], abs(i0))
            padding[1] = max(padding[1], abs(i1))
            padding[2] = max(padding[2], abs(i2))
        records.append(
            {
                "name": fc.name,
                "role": fc.role,
                "entry_point": fc.entry_point,
                "lifecycle_hook": fc.lifecycle_hook,
                "source": f"{fc.subdirectory}/{fc.name}.cpp",
                "required_padding": padding,
                "reads": sorted({access[0] for access in fc.accesses}),
            }
        )
    return records


def render_projection_json(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the projection.json manifest (whitepaper section 14.5).

    The projected and unchanged evolved fields are derived from the frozen
    CFunction bodies, not restated here.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Canonical JSON text.
    """
    records = _cfunction_records(snapshot, ("projection_block", "projection"))
    if not records:
        return (
            json.dumps(
                {
                    "status": "deferred",
                    "note": "deferred: no projection CFunction is registered in this profile",
                },
                sort_keys=True,
                indent=2,
            )
            + "\n"
        )
    evol = _names_in_group(snapshot, "EVOL")
    projected: List[str] = []
    for fc in snapshot.cfunctions:
        if fc.role == "projection_block":
            projected = _written_fields(fc.full_function, evol, "out_")
    payload = {
        "status": "registered",
        "constraints_enforced": [
            "det(gammabar)/det(gammahat) == 1",
            "gammabar^ij Atilde_ij == 0",
        ],
        "projected_fields": projected,
        "unchanged_fields": [name for name in evol if name not in set(projected)],
        "cfunctions": records,
        "structured_status": True,
        "calls_exit_on_failure": False,
        "floors_registered": [],
        "post_stage_projection": False,
        "scheduled": ["after_initial_data", "after_accepted_timestep"],
        "deferred_scheduling": ["after_remesh", "after_checkpoint_restore"],
        "note": (
            "post_stage projection is off by default (section 14.5): the "
            "supplied active host hook does no projection, and stage "
            "projection is a separate profile requiring stability evidence. "
            "Post-remesh and post-restore projection arrive with the AMR and "
            "checkpoint lifecycle in PR 11"
        ),
    }
    return json.dumps(payload, sort_keys=True, indent=2) + "\n"


def render_initial_data_json(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the initial_data.json manifest (sections 14.1/14.2/14.3).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Canonical JSON text.
    """
    roles = (
        "minkowski_initial_data_block",
        "minkowski_initial_data",
        "adm_to_evolved_block",
        "initialize_lambda_block",
    )
    records = _cfunction_records(snapshot, roles)
    evol = _names_in_group(snapshot, "EVOL")
    written_by_role = {
        fc.role: _written_fields(fc.full_function, evol, "out_")
        for fc in snapshot.cfunctions
        if fc.role in roles
    }
    payload = {
        "status": "registered" if records else "deferred",
        "cfunctions": records,
        "written_evolved_fields": {
            role: written_by_role[role] for role in sorted(written_by_role)
        },
        "adm_source_fields": _names_in_group(snapshot, "AUXEVOL"),
        "connection_pass": {
            "distinct_input_and_output_arrays": True,
            "written_fields": written_by_role.get("initialize_lambda_block", []),
        },
    }
    return json.dumps(payload, sort_keys=True, indent=2) + "\n"


def render_diagnostics_json(snapshot: FrozenNRPyDendroSnapshot) -> str:
    """
    Render the diagnostics.json manifest (section 14.6).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :return: Canonical JSON text.
    """
    records = _cfunction_records(snapshot, ("diagnostics_block", "diagnostics"))
    diag = _names_in_group(snapshot, "DIAG")
    if not records:
        return (
            json.dumps(
                {
                    "status": "deferred",
                    "note": "deferred: no diagnostic CFunction is registered in this profile",
                },
                sort_keys=True,
                indent=2,
            )
            + "\n"
        )
    written: List[str] = []
    for fc in snapshot.cfunctions:
        if fc.role == "diagnostics_block":
            written = _written_fields(fc.full_function, diag, "diag_")
    gen = _gen_lookup(snapshot)
    payload = {
        "status": "registered",
        "diagnostic_gridfunctions": diag,
        "written_fields": written,
        "cfunctions": records,
        "finite_difference_order": int(str(gen.get("fd_order", 0))),
        "recomputed_not_checkpoint_state": True,
        "exact_name_selection": True,
        "deferred": [
            "VTU output: the real writer arrives with the pinned Dendrolib "
            "host (section 13.2); this profile resolves output and refinement "
            "selections by exact NRPy name only",
            "L2 reductions: only the L-infinity reduction is exercised in the "
            "mock vehicle (section 14.6 lists further reductions as optional)",
        ],
    }
    return json.dumps(payload, sort_keys=True, indent=2) + "\n"


def render_file_hashes_json(file_hashes: Dict[str, str]) -> str:
    """
    Render the file_hashes.json manifest.

    :param file_hashes: Mapping of project-relative path to SHA-256 hex.
    :return: Canonical JSON text.
    """
    return (
        json.dumps(
            {"files": {path: file_hashes[path] for path in sorted(file_hashes)}},
            sort_keys=True,
            indent=2,
        )
        + "\n"
    )


def render_SHA256SUMS(file_hashes: Dict[str, str]) -> str:
    """
    Render the SHA256SUMS artifact (one ``<hash>  <path>`` line per file).

    :param file_hashes: Mapping of project-relative path to SHA-256 hex.
    :return: Text with one line per file, sorted by path, trailing newline.
    """
    return "".join(
        f"{digest}  {path}\n" for path, digest in sorted(file_hashes.items())
    )


def render_receipt_json(
    snapshot: FrozenNRPyDendroSnapshot, file_hashes: Dict[str, str]
) -> str:
    """
    Render the GENERATION_RECEIPT.json (and installation receipt body).

    :param snapshot: The frozen NRPy Dendro snapshot.
    :param file_hashes: Mapping of project-relative path to SHA-256 hex.
    :return: Canonical JSON text.
    """
    payload = {
        "generator": "nrpy.infrastructures.Dendro",
        "profile_name": snapshot.profile_name,
        "module_abi_hash": snapshot.hashes.module_abi_hash,
        "file_count": len(file_hashes),
        "file_hashes": {path: file_hashes[path] for path in sorted(file_hashes)},
    }
    return json.dumps(payload, sort_keys=True, indent=2) + "\n"


def all_manifests(
    snapshot: FrozenNRPyDendroSnapshot,
    dendrolib_capabilities_text: str,
) -> Dict[str, str]:
    """
    Render every fixed-name manifest keyed by project-relative path.

    :param snapshot: The frozen NRPy Dendro snapshot.
    :param dendrolib_capabilities_text: Canonical JSON text of the pinned
        capability record (copied verbatim into the project).
    :return: Mapping of ``manifest/<name>`` to file text.
    """
    manifests: Dict[str, str] = {}
    manifests["manifest/module.json"] = render_module_json(snapshot)
    manifests["manifest/generation_parameters.json"] = (
        render_generation_parameters_json(snapshot)
    )
    manifests["manifest/state.json"] = gridfunction_output.render_state_json(snapshot)
    manifests["manifest/parameters.json"] = (
        CodeParameters_output.render_parameters_json(snapshot)
    )
    manifests["manifest/cfunctions.json"] = render_cfunctions_json(snapshot)
    manifests["manifest/stencils.json"] = render_stencils_json(snapshot)
    manifests["manifest/dendrolib_capabilities.json"] = dendrolib_capabilities_text
    manifests["manifest/provenance.json"] = render_provenance_json(snapshot)
    manifests["manifest/projection.json"] = render_projection_json(snapshot)
    manifests["manifest/initial_data.json"] = render_initial_data_json(snapshot)
    manifests["manifest/diagnostics.json"] = render_diagnostics_json(snapshot)
    # Physical boundary policies are the only records this profile cannot
    # produce: they qualify in PR 10 (section 14.8).  An explicit deferral
    # record is emitted rather than invented content.
    manifests["manifest/boundaries.json"] = (
        json.dumps(
            {
                "status": "deferred",
                "note": "deferred: physical boundary policies qualify in PR 10",
            },
            sort_keys=True,
            indent=2,
        )
        + "\n"
    )
    return manifests


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
