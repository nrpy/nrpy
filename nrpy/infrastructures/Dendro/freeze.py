# nrpy/infrastructures/Dendro/freeze.py
"""
Freeze the mutable NRPy registries into an immutable Dendro snapshot.

The snapshot is a byte-stable copy of the authoritative registered objects
(gridfunctions, CodeParameters, CFunctions, Dendro roles and access
captures).  It is the only input to the pure target translators.  Freeze
validates the invariants of the frozen environment, computes independent
semantic hashes, and seals Dendro registration so that late mutation is
impossible.
"""

import hashlib
import json
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.Dendro import access_capture as cap
from nrpy.infrastructures.Dendro import naming
from nrpy.infrastructures.Dendro import registration as dendro_reg

# Host APIs that a Dendro CFunction may call without being a registered
# NRPy CFunction (whitepaper section 4.5: "a registered CFunction or an
# explicitly allowed host API").  The `ot__` prefix covers the pinned
# Dendrolib block API (`ot::Block::getOffset` and siblings); the adapter
# signature stays frozen until the I0-1 capability proof lands.
_HOST_API_ALLOWLIST = frozenset({"ot__Block__getOffset"})


def _is_allowed_host_api(callee: str) -> bool:
    """
    Return True if a callee name is an explicitly allowed host API.

    :param callee: Callee name listed in a Dendro role's `calls` tuple.
    :return: True for allowlisted names or the `ot__` Dendrolib prefix.
    """
    return callee in _HOST_API_ALLOWLIST or callee.startswith("ot__")


def body_uses_symbol(body: str, name: str) -> bool:
    """
    Return True if a C body references a symbol as a whole token.

    Word boundaries are `[A-Za-z0-9_]` on both sides, so `eta` does not
    match inside `Theta_fCCZ4` (whitepaper section 6.2 used-parameter
    closure; risk-table rule "emit only consumed CodeParameters").

    :param body: Registered CFunction body text.
    :param name: Registered CodeParameter name to look for.
    :return: True if the name occurs as a standalone token.
    """
    return (
        re.search(r"(?<![A-Za-z0-9_])" + re.escape(name) + r"(?![A-Za-z0-9_])", body)
        is not None
    )


# Types the generated translators can render end to end (struct member,
# default, validation, print, and TOML sample).  `#define` and the array
# families of whitepaper section 6.5 are deliberately absent: they have no
# implemented renderer, and accepting them here would emit either an
# ill-formed `auto` struct member or a silently zero-initialised array that is
# never validated or printed.  Failing at freeze is the honest behaviour
# (section 6.5: "unsupported or ambiguous types fail generation").
_SUPPORTED_CPARAM_TYPES = frozenset(
    {"REAL", "float", "double", "int", "bool", "DendroScalar"}
)


def _is_supported_cparam_type(cparam_type: str) -> bool:
    """
    Return True if a CodeParameter type has a specified TOML/C++ mapping.

    Whitepaper section 6.5 maps `REAL`, `float`, `double`, `int`, `bool`,
    `REAL[n]`, `int[n]`, and `#define`.  Only the scalar subset has an
    implemented renderer, so only that subset is accepted here; the frozen
    `DendroScalar` alias spelling is additionally accepted (it is the section
    4.1/6.6 scalar contract, verified against the frozen alias at freeze).
    Anything else fails generation rather than emitting code that cannot
    compile or that silently defaults to zero.

    :param cparam_type: Registered `cparam_type` string.
    :return: True if the type is supported.
    """
    return cparam_type in _SUPPORTED_CPARAM_TYPES


def _check_safe_subdirectory(subdirectory: str, name: str) -> None:
    """
    Reject unsafe registered CFunction output paths (whitepaper section 7.5).

    :param subdirectory: Registered subdirectory for one CFunction.
    :param name: CFunction name (for error messages).
    :raises ValueError: On absolute paths, `..`, or empty paths.
    """
    if not subdirectory:
        raise ValueError(f"CFunction {name} has unsafe subdirectory {subdirectory!r}.")
    normalized = os.path.normpath(subdirectory.replace("\\", "/"))
    if os.path.isabs(normalized) or normalized.startswith(".."):
        raise ValueError(f"CFunction {name} has unsafe subdirectory {subdirectory!r}.")
    # `os.path.normpath` collapses "." (the NRPy default subdirectory, which
    # means the module root) and any interior "./" segments, so only a real
    # escape reaches here.
    if any(part == ".." for part in normalized.split("/")):
        raise ValueError(f"CFunction {name} has unsafe subdirectory {subdirectory!r}.")


def _max_proven_padding() -> int:
    """
    Return the maximum block padding proven by the Dendrolib capability record.

    The gate is read from ``dendrolib_capabilities.json`` (whitepaper section
    13.2) rather than hardcoded, so extending the proof extends the gate and a
    narrowed proof narrows it.

    :return: ``padding.max_proven_padding`` from the capability record.
    :raises ValueError: If the record is missing or lacks a positive integer
        ``padding.max_proven_padding``.
    """
    path = os.path.join(os.path.dirname(__file__), "dendrolib_capabilities.json")
    try:
        with open(path, "r", encoding="utf-8") as handle:
            record = json.load(handle)
        value = record["padding"]["max_proven_padding"]
    except (OSError, KeyError, ValueError) as exc:
        raise ValueError(
            f"Cannot read padding.max_proven_padding from {path}: {exc}"
        ) from exc
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise ValueError(
            f"padding.max_proven_padding must be a positive integer, got {value!r}."
        )
    return value


def _sha256_of_lines(lines: Tuple[str, ...]) -> str:
    """
    Return the hex SHA-256 digest of canonicalized lines.

    :param lines: Canonicalized, deterministic text lines.
    :return: Hex digest.
    """
    payload = "\n".join(lines).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _canonical_json(obj: Any) -> str:
    """
    Return canonical JSON (sorted keys, compact separators) for hashing.

    :param obj: JSON-serializable object.
    :return: Canonical JSON string.
    """
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)


@dataclass(frozen=True)
class FrozenGridFunction:
    """
    Immutable copy of one registered gridfunction record.

    :param name: Exact registered NRPy gridfunction name.
    :param group: Group (EVOL, AUXEVOL, DIAG, AUX).
    :param description: Gridfunction description.
    :param rank: Tensor rank.
    :param dimension: Spatial dimension.
    :param gf_type: Generated scalar type alias.
    :param f_infinity: Value in the limit r -> infinity.
    :param wavespeed: Characteristic wave speed.
    :param is_basename: Whether this is the base name of a family.
    :param index: Group-local index assigned by freeze, in NRPy list order.
    """

    name: str
    group: str
    description: str
    rank: int
    dimension: int
    gf_type: str
    f_infinity: str
    wavespeed: str
    is_basename: bool
    index: int


@dataclass(frozen=True)
class FrozenCodeParameter:
    """
    Immutable copy of one used registered CodeParameter record.

    :param name: Exact registered CodeParameter name.
    :param module: Registering module.
    :param cparam_type: C type string.
    :param default_value: Registered default value, type-preserved.
    :param description: Registered description.
    :param commondata: Whether the parameter is common across grids.
    :param add_to_parfile: Whether the parameter appears in parfiles.
    :param consumers: Sorted names of frozen CFunctions using this parameter.
    """

    name: str
    module: str
    cparam_type: str
    default_value: object
    description: str
    commondata: bool
    add_to_parfile: bool
    consumers: Tuple[str, ...]


@dataclass(frozen=True)
class FrozenCFunction:
    """
    Immutable copy of one registered CFunction and its Dendro role.

    :param name: Exact registered CFunction name.
    :param subdirectory: Registered subdirectory (derived output path source).
    :param prototype: Registered function prototype.
    :param full_function: Registered, normalized full function text.
    :param role: Dendro scheduling role ("" if unregistered).
    :param entry_point: Whether Dendro invokes this CFunction directly.
    :param lifecycle_hook: Dendro lifecycle point, or empty.
    :param calls: Registered callee names.
    :param used_codeparameters: CodeParameter names used by this CFunction.
    :param accesses: Captured (gridfunction, i0, i1, i2) access tuples.
    """

    name: str
    subdirectory: str
    prototype: str
    full_function: str
    role: str
    entry_point: bool
    lifecycle_hook: str
    calls: Tuple[str, ...]
    used_codeparameters: Tuple[str, ...]
    accesses: Tuple[Tuple[str, int, int, int], ...]


@dataclass(frozen=True)
class FrozenHashes:
    """
    Independent semantic hashes of the frozen environment.

    :param state_schema_hash: Ordered gridfunction records and scalar contract.
    :param parameter_schema_hash: Used CodeParameter records.
    :param equation_hash: Canonical RHS expression digests and formulation
        options.
    :param stencil_hash: FD operators, exact coefficients, and access offsets.
    :param cfunction_api_hash: Names, prototypes, roles, call graph.
    :param cfunction_source_hash: Normalized registered full_function bytes.
    :param module_abi_hash: Combination of all preceding semantic hashes.
    """

    state_schema_hash: str
    parameter_schema_hash: str
    equation_hash: str
    stencil_hash: str
    cfunction_api_hash: str
    cfunction_source_hash: str
    module_abi_hash: str


@dataclass(frozen=True)
class FrozenNRPyDendroSnapshot:
    """
    Immutable snapshot of the NRPy environment for Dendro export.

    :param profile_name: Human-readable generation profile name.
    :param generation_parameters: Frozen (name, value) pairs, sorted.
    :param gridfunctions: All registered gridfunctions, in NRPy group-list
        order with group-local indices.
    :param codeparameters: Used CodeParameters, sorted by name.
    :param cfunctions: Registered CFunctions, sorted by name.
    :param required_padding: (padding_x, padding_y, padding_z) reduced from
        the captured accesses of every CFunction.
    :param hashes: Independent semantic hashes.
    :param extras: Deterministic builder sidecar (name, payload) pairs,
        sorted by name.  Payloads must be JSON-serializable.  This is the
        single-derivation channel for builder-computed audit records that
        the pure translators must render (e.g. the PR5 operator manifest,
        whose exact rational coefficients come from
        ``compute_fdcoeffs_fdstencl``); the exporter never re-derives them.
    """

    profile_name: str
    generation_parameters: Tuple[Tuple[str, object], ...]
    gridfunctions: Tuple[FrozenGridFunction, ...]
    codeparameters: Tuple[FrozenCodeParameter, ...]
    cfunctions: Tuple[FrozenCFunction, ...]
    required_padding: Tuple[int, int, int]
    hashes: FrozenHashes
    extras: Tuple[Tuple[str, object], ...] = ()


def _gridfunction_order_lines(fgs: Tuple[FrozenGridFunction, ...]) -> Tuple[str, ...]:
    """
    Canonical lines for the state schema hash.

    :param fgs: Frozen gridfunctions.
    :return: One line per record, in registry order.
    """
    return tuple(
        f"{fg.index}|{fg.group}|{fg.name}|{fg.rank}|{fg.dimension}|{fg.gf_type}"
        f"|{fg.f_infinity}|{fg.wavespeed}|{fg.is_basename}"
        for fg in fgs
    )


def canonical_expression_digest(expr: Any) -> str:
    """
    Return a canonical, DAG-aware digest of one SymPy expression.

    The digest is built bottom-up over the expression graph, memoized by node
    identity, so shared subexpressions are hashed once: the cost is
    proportional to the DAG size (a few thousand nodes for the fCCZ4 RHSs)
    rather than to the exponentially larger expanded tree that ``srepr``
    would materialize.  Any change to an operator, an operand, or a rational
    coefficient changes the digest, so this is a physics fingerprint rather
    than a fingerprint of the emitted C text (whitepaper section 9.10).

    :param expr: A SymPy expression.
    :return: Hex SHA-256 digest of the canonical expression graph.
    """
    memo: Dict[int, str] = {}
    stack: List[Tuple[Any, bool]] = [(expr, False)]
    while stack:
        node, expanded = stack.pop()
        key = id(node)
        if key in memo:
            continue
        if not expanded:
            stack.append((node, True))
            for arg in node.args:
                if id(arg) not in memo:
                    stack.append((arg, False))
            continue
        digest = hashlib.sha256()
        digest.update(type(node).__name__.encode("utf-8"))
        if not node.args:
            digest.update(repr(node).encode("utf-8"))
        for arg in node.args:
            digest.update(memo[id(arg)].encode("utf-8"))
        memo[key] = digest.hexdigest()
    return memo[id(expr)]


def _canonical_equation_lines(
    extras: Tuple[Tuple[str, object], ...],
) -> Tuple[str, ...]:
    """
    Canonical equation lines for the equation hash.

    The lines are the builder-supplied canonical expression digests (see
    :func:`canonical_expression_digest`), which fingerprint the symbolic
    right-hand sides themselves.  When a profile supplies none, an explicit
    marker is hashed instead of silently falling back to a weaker input, so
    an unfingerprinted profile is visibly distinct from a fingerprinted one.

    :param extras: The frozen builder sidecar entries.
    :return: Deterministic lines describing the canonical expressions.
    """
    canonical = dict(extras).get("rhs_canonical")
    if not isinstance(canonical, dict) or not canonical:
        return ("no_canonical_expressions",)
    return tuple(
        f"{name}|{canonical[name]}"
        for name in sorted(canonical, key=lambda n: (n.lower(), n))
    )


def _stencil_lines(
    fcs: Tuple[FrozenCFunction, ...], accesses: Tuple[Tuple[str, int, int, int], ...]
) -> Tuple[str, ...]:
    """
    Canonical finite-difference stencil lines, in deterministic order.

    The exact rational coefficients and derivative operators (lines referencing
    an ``invdxx`` stride) from every frozen CFunction body, plus the sorted
    unique signed access offsets.

    :param fcs: Frozen CFunctions.
    :param accesses: All captured (gridfunction, i0, i1, i2) access tuples.
    :return: One line per stencil fact, in deterministic order.
    """
    lines: List[str] = []
    for fc in fcs:
        for line in fc.full_function.splitlines():
            if "invdxx" in line or "_Rational_" in line:
                lines.append(line.strip())
    for name, i0, i1, i2 in sorted(set(accesses)):
        lines.append(f"offset|{name}|{i0},{i1},{i2}")
    return tuple(lines)


def _validate_builder_invariants(
    *,
    extras: Dict[str, Any],
    evol_names: Tuple[str, ...],
    padding: Tuple[int, int, int],
    gen_lookup: Dict[str, object],
    has_rhs_role: bool,
    allow_missing_builder_records: bool,
) -> None:
    """
    Enforce the section 9.9 freeze invariants that depend on builder records.

    These are steps 4 (RHS/EVOL bijection), 9 (operator-derived reach versus
    captured reach) and 11 (KO ownership and derivative backend).  They run
    here, at the one point every Dendro generation passes, rather than inside
    a single builder: a later builder (projection, diagnostics, boundaries)
    that registers evolved-field writers is then covered automatically.

    A builder that records any sidecar at all must record both of these, and
    both are then enforced: a physics profile cannot opt out of the checks by
    omitting half its records.  A profile that registers an ``rhs`` role and
    records no sidecar at all is refused unless it explicitly declares itself a
    probe (``allow_missing_builder_records``), so forgetting the records is an
    error rather than a silent skip.

    :param extras: The builder sidecar (``glb_extras_dict["Dendro"]["builder_extras"]``).
    :param evol_names: Registered EVOL names, in NRPy order.
    :param padding: Access-capture-derived (px, py, pz).
    :param gen_lookup: Frozen generation parameters.
    :param has_rhs_role: Whether any frozen CFunction carries an ``rhs`` role.
    :param allow_missing_builder_records: Permit an ``rhs``-role profile with
        no sidecar at all (toy and unit-test environments only).
    :raises ValueError: If a required record is missing or an invariant fails.
    """
    if not extras:
        if has_rhs_role and not allow_missing_builder_records:
            raise ValueError(
                "A CFunction carries an 'rhs' role but no builder records were "
                "supplied, so the section 9.9 step 4/9/11 checks cannot run.  A "
                "physics profile must record 'operator_manifest' and "
                "provenance 'rhs_symbols'; a probe kernel may pass "
                "allow_missing_builder_records=True."
            )
        return
    provenance = extras.get("provenance")
    provenance = provenance if isinstance(provenance, dict) else {}
    manifest = extras.get("operator_manifest")
    manifest = manifest if isinstance(manifest, dict) else None
    rhs_symbols = provenance.get("rhs_symbols")

    if has_rhs_role and manifest is None:
        raise ValueError(
            "A CFunction carries an 'rhs' role and builder records were "
            "supplied, but no operator manifest is among them: the section "
            "9.9 step-9 stencil cross-check cannot be performed."
        )
    if has_rhs_role and not rhs_symbols:
        raise ValueError(
            "A CFunction carries an 'rhs' role and builder records were "
            "supplied, but no rhs_symbols provenance is among them: the "
            "section 9.9 step-4 RHS/EVOL bijection cannot be verified."
        )

    # Step 4: the RHS symbols map bijectively onto the registered EVOL fields.
    if rhs_symbols:
        mapped = [
            naming.rhs_symbol_to_gridfunction_name(str(symbol))
            for symbol in rhs_symbols
        ]
        if len(set(mapped)) != len(mapped):
            raise ValueError(
                f"RHS symbols map onto duplicate gridfunctions: {sorted(mapped)}."
            )
        if set(mapped) != set(evol_names):
            raise ValueError(
                "RHS targets are not a bijection onto the registered EVOL "
                f"fields: missing={sorted(set(evol_names) - set(mapped))} "
                f"extra={sorted(set(mapped) - set(evol_names))}."
            )

    if manifest is None:
        return

    # Step 9: the operator-derived reach must equal the captured reach, per
    # axis.  One coefficient source (compute_fdcoeffs_fdstencl) produced the
    # operator records; the accesses came from the emitted memory reads.
    operator_axes = [0, 0, 0]
    for record in manifest.get("operators", ()):
        per_axis = record.get("max_offset_per_axis", (0, 0, 0))
        for axis in range(3):
            operator_axes[axis] = max(operator_axes[axis], int(per_axis[axis]))
    if tuple(operator_axes) != tuple(padding):
        raise ValueError(
            f"Operator-derived reach {tuple(operator_axes)} does not equal the "
            f"access-capture-derived padding {tuple(padding)} (section 8.6)."
        )

    # Step 11: KO ownership and derivative backend match the frozen profile.
    if bool(manifest.get("ko_enabled")) != bool(gen_lookup.get("Dendro_enable_KO")):
        raise ValueError(
            f"Operator manifest ko_enabled={manifest.get('ko_enabled')!r} does "
            f"not match Dendro_enable_KO={gen_lookup.get('Dendro_enable_KO')!r} "
            "(section 8.5: NRPy owns KO exactly once)."
        )
    if manifest.get("derivative_backend") != gen_lookup.get(
        "Dendro_derivative_backend"
    ):
        raise ValueError(
            f"Operator manifest derivative_backend="
            f"{manifest.get('derivative_backend')!r} does not match the frozen "
            f"{gen_lookup.get('Dendro_derivative_backend')!r}."
        )
    if manifest.get("ko_owner") != "nrpy_full_stencil":
        raise ValueError(
            f"Unexpected ko_owner {manifest.get('ko_owner')!r}: the canonical "
            "backend requires NRPy-owned Kreiss-Oliger dissipation."
        )


def freeze_nrpy_dendro_environment(
    profile_name: str,
    tx: Any = None,
    allow_missing_builder_records: bool = False,
) -> FrozenNRPyDendroSnapshot:
    """
    Validate, freeze, and hash the current NRPy environment for Dendro.

    :param profile_name: Human-readable generation profile name.
    :param tx: Optional generation transaction; when given it is advanced
        to ``FROZEN`` so the section 4.4 phase machine is integrated
        (invalid transitions raise).
    :param allow_missing_builder_records: Permit a profile that registers an
        ``rhs`` role to supply no builder sidecar, which skips the section 9.9
        step 4/9/11 checks.  Only a toy or unit-test environment registering a
        probe kernel may set this; a physics profile must supply the records.
    :return: The immutable snapshot.
    :raises ValueError: If Infrastructure is not Dendro or a freeze
        invariant is violated.
    """
    if tx is not None:
        advance = getattr(tx, "advance_to", None)
        if not callable(advance):
            raise ValueError(
                "freeze tx must be a GenerationTransaction with advance_to()."
            )
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            f"Infrastructure must be 'Dendro' to freeze, got "
            f"{par.parval_from_str('Infrastructure')!r}."
        )

    # Step 1: Read all four NRPy gridfunction lists and assign group-local
    # indices in the exact NRPy (deterministic, case-insensitive) order.
    # Section 9.9: gridfunction names must be exactly unique across groups.
    groups = ("EVOL", "AUXEVOL", "DIAG", "AUX")
    lists = gri.GridFunction.gridfunction_lists()
    seen_names = [name for names in lists for name in names]
    if len(set(seen_names)) != len(seen_names):
        raise ValueError("Gridfunction names are not unique across groups.")
    if set(seen_names) != set(gri.glb_gridfcs_dict):
        raise ValueError("Gridfunction lists do not cover the gridfunction registry.")
    fgs: List[FrozenGridFunction] = []
    for group, names in zip(groups, lists):
        for index, name in enumerate(names):
            gf = gri.glb_gridfcs_dict[name]
            fgs.append(
                FrozenGridFunction(
                    name=gf.name,
                    group=group,
                    description=gf.desc,
                    rank=gf.rank,
                    dimension=gf.dimension,
                    gf_type=gf.gf_type,
                    f_infinity=str(gf.f_infinity),
                    wavespeed=str(gf.wavespeed),
                    is_basename=gf.is_basename,
                    index=index,
                )
            )

    # Step 2: Per-function used CodeParameters via token-aware closure
    # (section 6.2; section 9.7 provenance).  The global closure is the union
    # of per-function sets, derived and never hand-typed.
    bodies = {
        name: cfc.CFunction_dict[name].full_function
        for name in sorted(cfc.CFunction_dict)
    }
    per_function_uses: Dict[str, Tuple[str, ...]] = {}
    for cf_name, body in bodies.items():
        per_function_uses[cf_name] = tuple(
            sorted(
                cp_name
                for cp_name in par.glb_code_params_dict
                if body_uses_symbol(body, cp_name)
            )
        )
    used_cparams = sorted({cp for uses in per_function_uses.values() for cp in uses})
    # Section 6.5: unsupported or ambiguous CodeParameter types fail here.
    for cp_name in used_cparams:
        cparam_type = par.glb_code_params_dict[cp_name].cparam_type
        if not _is_supported_cparam_type(cparam_type):
            raise ValueError(
                f"CodeParameter {cp_name!r} has unsupported type "
                f"{cparam_type!r} (section 6.5)."
            )
    consumers_by_param: Dict[str, Tuple[str, ...]] = {
        cp_name: tuple(
            sorted(
                cf_name
                for cf_name, uses in per_function_uses.items()
                if cp_name in uses
            )
        )
        for cp_name in used_cparams
    }
    frozen_cparams = tuple(
        FrozenCodeParameter(
            name=name,
            module=par.glb_code_params_dict[name].module,
            cparam_type=par.glb_code_params_dict[name].cparam_type,
            default_value=par.glb_code_params_dict[name].defaultvalue,
            description=par.glb_code_params_dict[name].description,
            commondata=par.glb_code_params_dict[name].commondata,
            add_to_parfile=par.glb_code_params_dict[name].add_to_parfile,
            consumers=consumers_by_param[name],
        )
        for name in used_cparams
    )

    # Step 3: CFunctions (sorted by name) with role/access sidecars.
    # A Dendro role marks a numerical CFunction.  Helpers are role-less
    # CFunctions without captures; they are allowed and contribute no
    # padding.  Every leaf numerical CFunction (a role with no `calls`)
    # must have a capture: a missing capture fails here instead of
    # silently yielding padding (0,0,0) (sections 8.6/18.1
    # insufficient-padding fatal).  Driver entry points whose body only
    # calls other CFunctions carry no direct accesses and need no capture.
    roles = dendro_reg.get_CFunction_roles()
    captured_names = set(cap.captured_cfunction_names())
    frozen_cfs = []
    for name in sorted(cfc.CFunction_dict):
        cf = cfc.CFunction_dict[name]
        role = roles.get(name, {})
        _check_safe_subdirectory(cf.subdirectory, name)
        if name in roles and name not in captured_names and not role.get("calls"):
            raise ValueError(
                f"Numerical CFunction {name!r} has a Dendro role but no "
                "access capture; every leaf numerical CFunction must be "
                "built inside capture_gridfunction_accesses."
            )
        frozen_cfs.append(
            FrozenCFunction(
                name=cf.name,
                subdirectory=cf.subdirectory,
                prototype=cf.function_prototype,
                full_function=cf.full_function,
                role=role.get("role", ""),
                entry_point=bool(role.get("entry_point", False)),
                lifecycle_hook=role.get("lifecycle_hook", ""),
                calls=tuple(role.get("calls", ())),
                used_codeparameters=per_function_uses[name],
                accesses=(
                    cap.get_captured_accesses(name) if name in captured_names else ()
                ),
            )
        )
    frozen_cfunctions = tuple(frozen_cfs)

    # Step 4: Role invariants: every role key is a registered CFunction;
    # every CFunction with captured (numerical) accesses has exactly one
    # role.  Role-less CFunctions without captures are non-numerical
    # helpers (section 4.5).
    for role_name in roles:
        if role_name not in cfc.CFunction_dict:
            raise ValueError(f"Dendro role {role_name!r} is not a CFunction.")
    for name in cfc.CFunction_dict:
        if name in captured_names and name not in roles:
            raise ValueError(
                f"Numerical CFunction {name!r} has captured accesses but no "
                "Dendro role."
            )

    # Step 5: Call-graph closure: every listed callee is a registered
    # CFunction or an explicitly allowed host API (section 4.5).
    for fc in frozen_cfunctions:
        for callee in fc.calls:
            if callee not in cfc.CFunction_dict and not _is_allowed_host_api(callee):
                raise ValueError(
                    f"CFunction {fc.name} lists callee {callee!r}, which is "
                    "not a registered CFunction or allowed host API."
                )

    # Step 6: Reduce captured accesses to per-axis padding; verify each
    # captured name exists in the registry.
    all_accesses = []
    for fc in frozen_cfunctions:
        for access in fc.accesses:
            all_accesses.append(access)
            if access[0] not in gri.glb_gridfcs_dict:
                raise ValueError(
                    f"Captured access names unknown gridfunction {access[0]!r}."
                )
    padding = [0, 0, 0]
    for _name, i0, i1, i2 in all_accesses:
        padding[0] = max(padding[0], abs(i0))
        padding[1] = max(padding[1], abs(i1))
        padding[2] = max(padding[2], abs(i2))

    # Sections 2.5/8.6/18.1: insufficient padding is fatal, and any padding
    # beyond what the Dendrolib capability record proves is capability-gated.
    # There is no fallback to a lower-order derivative.
    max_proven = _max_proven_padding()
    if max(padding) > max_proven:
        raise ValueError(
            f"Required padding {tuple(padding)} exceeds the maximum padding "
            f"{max_proven} recorded in dendrolib_capabilities.json: pin "
            "Dendrolib and extend that capability record before generating "
            "this profile."
        )

    # Step 7: Generation parameters (frozen copy, sorted).
    from nrpy.infrastructures.Dendro.generation_parameters import (
        get_generation_parameter_view,
        validate_generation_parameters,
    )

    validate_generation_parameters()
    gen_params = tuple(
        (name, value) for name, value in sorted(get_generation_parameter_view().items())
    )

    # Sections 4.1/6.6: both scalar spellings come from frozen parameters.
    # Every frozen gridfunction type must equal the frozen scalar alias;
    # loop emitters read the live alias at build time, so a mismatch here
    # proves build/freeze skew and fails generation.
    gen_lookup = dict(gen_params)
    frozen_scalar = gen_lookup.get("Dendro_scalar_type")
    for fg in fgs:
        if fg.gf_type != frozen_scalar:
            raise ValueError(
                f"Gridfunction {fg.name!r} has scalar type {fg.gf_type!r} "
                f"but frozen Dendro_scalar_type is {frozen_scalar!r}."
            )

    # Step 7a: Builder sidecar extras (e.g. the operator manifest with its
    # exact rational coefficients, computed once by the builder from
    # compute_fdcoeffs_fdstencl).  JSON-serializability is validated here so
    # the manifests can render them without a second derivation pass.  These
    # are hashed into the stencil hash, so they must be built before step 8.
    raw_extras = par.glb_extras_dict.get("Dendro", {}).get("builder_extras", {})
    if not isinstance(raw_extras, dict):
        raise ValueError("Dendro builder_extras must be a dict.")
    extras: List[Tuple[str, object]] = []
    for name in sorted(raw_extras):
        try:
            _canonical_json(raw_extras[name])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Dendro builder extra {name!r} is not JSON-serializable."
            ) from exc
        extras.append((name, raw_extras[name]))

    # Steps 4/9/11 (section 9.9): RHS/EVOL bijection, operator-versus-access
    # reach cross-check, and KO/derivative-backend ownership.
    _validate_builder_invariants(
        extras=raw_extras if isinstance(raw_extras, dict) else {},
        evol_names=tuple(lists[0]),
        padding=(padding[0], padding[1], padding[2]),
        gen_lookup=gen_lookup,
        has_rhs_role=any(fc.role.startswith("rhs") for fc in frozen_cfunctions),
        allow_missing_builder_records=allow_missing_builder_records,
    )

    # Step 8: Independent semantic hashes.
    state_schema_hash = _sha256_of_lines(_gridfunction_order_lines(tuple(fgs)))
    parameter_schema_hash = _sha256_of_lines(
        tuple(_canonical_json(cp.__dict__) for cp in frozen_cparams)
    )
    equation_hash = _sha256_of_lines(
        tuple(f"{name}={value}" for name, value in gen_params)
        + _canonical_equation_lines(tuple(extras))
    )
    stencil_hash = _sha256_of_lines(
        _stencil_lines(frozen_cfunctions, tuple(all_accesses))
        + tuple(f"extra|{name}|{_canonical_json(payload)}" for name, payload in extras)
    )
    cfunction_api_hash = _sha256_of_lines(
        tuple(
            f"{fc.name}|{fc.subdirectory}|{fc.prototype}|{fc.role}"
            f"|{fc.entry_point}|{fc.lifecycle_hook}|{' '.join(fc.calls)}"
            for fc in frozen_cfunctions
        )
    )
    cfunction_source_hash = _sha256_of_lines(
        tuple(f"{fc.name}\n{fc.full_function}" for fc in frozen_cfunctions)
    )
    module_abi_hash = _sha256_of_lines(
        (
            state_schema_hash,
            parameter_schema_hash,
            equation_hash,
            stencil_hash,
            cfunction_api_hash,
            cfunction_source_hash,
        )
    )
    hashes = FrozenHashes(
        state_schema_hash=state_schema_hash,
        parameter_schema_hash=parameter_schema_hash,
        equation_hash=equation_hash,
        stencil_hash=stencil_hash,
        cfunction_api_hash=cfunction_api_hash,
        cfunction_source_hash=cfunction_source_hash,
        module_abi_hash=module_abi_hash,
    )

    snapshot = FrozenNRPyDendroSnapshot(
        profile_name=profile_name,
        generation_parameters=gen_params,
        gridfunctions=tuple(fgs),
        codeparameters=frozen_cparams,
        cfunctions=frozen_cfunctions,
        required_padding=tuple(padding),  # type: ignore[arg-type]
        hashes=hashes,
        extras=tuple(extras),
    )

    # Step 9: Seal Dendro registration; further role registration raises.
    dendro_reg.set_registration_open(False)
    if tx is not None:
        tx.advance_to("FROZEN")
    return snapshot


def assert_mutable_registries_match(snapshot: FrozenNRPyDendroSnapshot) -> None:
    """
    Verify the mutable registries still match the frozen snapshot.

    Bodies, prototypes, subdirectories, roles, calls, CodeParameter
    defaults, and the token-aware use closure are compared, not just name
    sets (whitepaper section 4.3: recalculate digest before/after export
    and fail if it differs).

    :param snapshot: The snapshot to compare against.
    :raises ValueError: If any compared record drifted since freeze.
    """
    if sorted(gri.glb_gridfcs_dict) != sorted(fg.name for fg in snapshot.gridfunctions):
        raise ValueError("Gridfunction registry drifted since freeze.")
    if sorted(cfc.CFunction_dict) != sorted(fc.name for fc in snapshot.cfunctions):
        raise ValueError("CFunction registry drifted since freeze.")
    roles = dendro_reg.get_CFunction_roles()
    for fc in snapshot.cfunctions:
        live = cfc.CFunction_dict[fc.name]
        if (
            live.function_prototype != fc.prototype
            or live.full_function != fc.full_function
            or live.subdirectory != fc.subdirectory
        ):
            raise ValueError(
                f"CFunction {fc.name!r} body/prototype/path drifted since freeze."
            )
        live_role = roles.get(fc.name, {})
        if (
            live_role.get("role", "") != fc.role
            or bool(live_role.get("entry_point", False)) != fc.entry_point
            or live_role.get("lifecycle_hook", "") != fc.lifecycle_hook
            or tuple(live_role.get("calls", ())) != fc.calls
        ):
            raise ValueError(
                f"CFunction {fc.name!r} role metadata drifted since freeze."
            )
        live_uses = tuple(
            sorted(
                cp_name
                for cp_name in par.glb_code_params_dict
                if body_uses_symbol(live.full_function, cp_name)
            )
        )
        if live_uses != fc.used_codeparameters:
            raise ValueError(
                f"CFunction {fc.name!r} CodeParameter use drifted since freeze."
            )
    live_used = {
        cp_name
        for cf in cfc.CFunction_dict.values()
        for cp_name in par.glb_code_params_dict
        if body_uses_symbol(cf.full_function, cp_name)
    }
    if live_used != {cp.name for cp in snapshot.codeparameters}:
        raise ValueError("CodeParameter use closure drifted since freeze.")
    for cp in snapshot.codeparameters:
        if cp.name not in par.glb_code_params_dict:
            raise ValueError(f"CodeParameter {cp.name!r} missing since freeze.")
        if par.glb_code_params_dict[cp.name].defaultvalue != cp.default_value:
            raise ValueError(f"CodeParameter {cp.name!r} default drifted since freeze.")
    live_extras_raw = par.glb_extras_dict.get("Dendro", {}).get("builder_extras", {})
    live_extras = tuple(
        (name, live_extras_raw[name]) for name in sorted(live_extras_raw)
    )
    if live_extras != snapshot.extras:
        raise ValueError("Dendro builder extras drifted since freeze.")


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
