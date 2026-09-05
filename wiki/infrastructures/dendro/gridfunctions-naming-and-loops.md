# Gridfunctions, Naming, And Loops

> Explain the Dendro gridfunction class, the exact-name role decorations, the CFunction role sidecar, the generation parameters, and the two loop helpers. · Status: provisional · Last reconciled: 09-05-2026
> Up: [Dendro](index.md)

## Summary

The Dendro infrastructure keeps every registered NRPy name byte for byte and
derives everything machine-readable from it. Names reach generated code only
through reversible syntactic decorations, scheduling roles are a one-word
sidecar beside the CFunction registry, the non-scientific lowering choices are
ordinary NRPy parameters, and both loop helpers are thin layers over the
generic NRPy loop helper.

## Detail

### The gridfunction class

`DendroGridFunction` is a plain string formatter, like its `BHaHGridFunction`,
`ETLegacyGridFunction`, and `CarpetXGridFunction` peers. Its one-point read
helper returns a role-prefixed input pointer indexed by the base interior
index `pp` plus signed x-fastest offsets, and it rejects a SIMD read because
the qualified CPU profile is not SIMD-qualified. It imports nothing from
`nrpy.infrastructures`.

### Exact names, no aliases

The naming module permits only syntactic decorations that map one-to-one onto
the exact NRPy name. Semantic aliases are prohibited: `cf` stays `cf` and never
becomes `chi`, `lambdaU` never becomes `Gt`, and `Theta_fCCZ4` never becomes
`Theta`. The decorations are the role pointers `in_`, `rhs_`, `out_`, and
`diag_`, plus the enum member form. The four pointers are prefixes applied
to the exact registered name and the enum member is the name itself, so every
decoration is reversible by construction and a generated identifier always
traces back to one registered gridfunction. `validate_cpp_identifier` enforces
that a name is a legal C or C++ identifier.

`out_` exists because the other three roles all mean something else. Reusing
`rhs_` for the initial-data writer made its generated signature claim it fills
the right-hand-side vector, which is how a caller silently zeroes the state.

A separate helper, `rhs_symbol_to_gridfunction_name`, maps NRPy's RHS *symbol*
convention onto the registered gridfunction name — `h_rhsDD00` becomes `hDD00`
— which is the bijection the RHS builder asserts against the registry. It does
not invert the `rhs_` pointer decoration.

The recorded decision behind this is that the fCCZ4 state keeps its exact NRPy
names and the `GridFunction.gridfunction_lists()` order, stores native `hDD`
and native evolved `lambdaU` rather than substituting BSSN quantities, and
treats runtime physics parameters as the registered `CodeParameter` objects
only.

### The role sidecar

A Dendro kernel is an ordinary NRPy CFunction plus one word: its scheduling
role. `register_Dendro_CFunction` registers the CFunction through
`cfc.register_CFunction` and records the role in
`par.glb_extras_dict["Dendro"]["CFunction_roles"]`;
`CFunction_name_for_role` reads it back, which is how the host-adapter emitters
ask for "the all-block RHS entry point" without taking a dozen name arguments.
The sidecar holds no body, signature, parameter default, field declaration, or
source path, because the CFunction registry is the only body and signature
store. Duplicate names are rejected by `cfc.register_CFunction` itself.

The registered EVOL, AUXEVOL, and DIAG orders are queried from the gridfunction
registry rather than restated. The right-hand-side builder records two derived
values there as well: the ghost points its operators reach, through
`set_required_padding`, and the EVOL fields its kernel upwinds on, through
`set_upwind_control_fields`. Both readers raise when nothing is recorded, so a
lost writer fails generation rather than emitting a zero padding or an empty
upwind table. `registered_evol_order` raises on an empty EVOL registry for the
same reason: a builder that reached it with nothing registered would emit a
well-formed kernel with an empty body, which is the one failure a generated
solver cannot report.

Claim evidence:
- Claim: the Dendro role sidecar stores exactly one scheduling role string per registered CFunction name and no other metadata, and duplicate registration is rejected by `nrpy.c_function.register_CFunction` rather than by the sidecar.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/registration.py`, `register_Dendro_CFunction` and its doctest
- Corroboration: `nrpy/c_function.py`, `register_CFunction` duplicate-name `ValueError`

### Loops

The interior point loop is x-fastest with `i0` innermost over the padded block
interior. It emits the interior base index `pp`, which the gridfunction class's
one-point read helper uses for every read, along with the interior coordinates
`xx0`, `xx1`, and `xx2`. The block loop iterates the local block list supplied
by Dendro and invokes the registered per-block CFunction for each block. Both
are emitted by NRPy through the generic loop helper; the Dendro runtime supplies
only the block list and its count. `require_serial_parallelization` keeps the
emitted loop free of a parallelization directive NRPy is not entitled to choose
for the host: the point kernel runs inside Dendro's own block traversal, so an
inner OpenMP pragma would nest parallelism.

### Generation parameters

The non-scientific lowering choices are ordinary NRPy parameters registered
through `par.register_param`: the scalar alias and the Kreiss-Oliger switch.
There is no parallel configuration object; the values
live in the parameter registry and every builder reads them from there.
`validate_generation_parameters` rejects an unsupported combination before
anything is lowered, so an unqualified profile fails generation instead of
producing silently wrong output. The runtime parameter default, validation, and
print CFunctions are registered last, after the scientific CFunctions have put
every CodeParameter they use into the registry.

## Sources

- [grid.py](../../../nrpy/grid.py) - `DendroGridFunction`, `read_gf_from_memory_Ccode_onept`
- [naming.py](../../../nrpy/infrastructures/Dendro/naming.py) - `input_pointer`, `rhs_pointer`, `out_pointer`, `enum_member`, `rhs_symbol_to_gridfunction_name`, `validate_cpp_identifier`
- [registration.py](../../../nrpy/infrastructures/Dendro/registration.py) - `register_Dendro_CFunction`, `CFunction_name_for_role`, `registered_evol_order`, `set_required_padding`, `set_upwind_control_fields`
- [simple_loop.py](../../../nrpy/infrastructures/Dendro/simple_loop.py) - `simple_loop`, `require_serial_parallelization`
- [block_loop.py](../../../nrpy/infrastructures/Dendro/block_loop.py) - `block_loop`
- [generation_parameters.py](../../../nrpy/infrastructures/Dendro/generation_parameters.py) - `validate_generation_parameters`
- [parameters.py](../../../nrpy/infrastructures/Dendro/runtime/parameters.py) - `register_CFunctions_parameters`
- [ADR_dendro_names.md](../../../nrpy/infrastructures/Dendro/ADR_dendro_names.md) - `Decision`, `Consequences`

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
- See also: [Project Assembly And Emitters](project-assembly-and-emitters.md)
- Example: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
- See also: [Finite Difference](../../core/finite-difference.md)
- See also: [Infrastructure Code Style](../infrastructure-code-style.md)
