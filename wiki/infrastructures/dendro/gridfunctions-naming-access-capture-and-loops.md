# Gridfunctions, Naming, Access Capture, And Loops

> Explain the Dendro gridfunction class, exact-name role decorations, CFunction role metadata, access-derived padding, and the interior point and block loops. · Status: provisional · Last reconciled: 09-04-2026
> Up: [Dendro](index.md)

## Summary

The Dendro backend keeps every registered NRPy name byte-for-byte and derives
everything machine-readable from it. `DendroGridFunction` records each one-point
memory access as it is emitted, so the padding a kernel needs follows from the
accesses that kernel actually made rather than from a second stencil model.
Names reach generated code only through reversible syntactic decorations,
scheduling metadata is a non-authoritative sidecar beside the CFunction
registry, and both loop helpers are thin layers over the generic NRPy loop
helper.

## Detail

### Exact names, no aliases

The naming module permits only syntactic decorations that map one-to-one onto
the exact NRPy name. Semantic aliases are prohibited: `cf` stays `cf` and never
becomes `chi`, `lambdaU` never becomes `Gt`, and `Theta_fCCZ4` never becomes
`Theta`. The decorations are the role pointers `in_`, `rhs_`, `out_`, `diag_`,
and `aux_`, plus the enum member form. The five pointers are prefixes applied to
the exact registered name and the enum member is the name itself, so every
decoration is reversible by construction and a generated identity can always be
traced to one registered gridfunction. `validate_cpp_identifier`
enforces that a name is a legal C or C++ identifier. A separate helper,
`rhs_symbol_to_gridfunction_name`, maps NRPy's RHS *symbol* convention onto the
registered gridfunction name — `h_rhsDD00` becomes `hDD00` — which is the
bijection the RHS builder asserts against the registry. It does not invert the
`rhs_` pointer decoration.

The recorded decision behind this is that the fCCZ4 state keeps its exact NRPy
names and the `GridFunction.gridfunction_lists()` order, stores native `hDD` and
native evolved `lambdaU` rather than substituting BSSN quantities, and treats
runtime physics parameters as the used registered `CodeParameter` objects only.
The checkpoint consequence of that decision is recorded on
[Frozen Snapshot And Pure Translators](frozen-snapshot-and-pure-translators.md).

### Access capture and padding

During code generation, `DendroGridFunction` records every one-point access as
a gridfunction name plus signed offsets into the active capture context.
`capture_gridfunction_accesses` scopes that context to the CFunction being
built, so `required_padding` derives the block padding from the exact offsets
that function emitted. The recorded data lives in the Dendro area of
`par.glb_extras_dict` and holds only names and signed integers: no function
bodies, no field declarations, and no parallel state.

Deriving padding from captured accesses rather than from a declared stencil
width removes a class of silent mismatch, because the padding and the emitted
reads cannot disagree by construction. It also means capture must be active while a
kernel is generated: asking for the padding of a CFunction that recorded
nothing raises rather than returning zero, so a missing capture scope is an
error instead of a silent zero-padding result. A writer that reads nothing at all — the Minkowski
fill and the perturbation writer — satisfies that gate by recording an explicit
empty capture through `record_empty_capture`. A kernel that reads only its own
point reaches zero padding a different way: the algebraic projection opens an
ordinary capture whose recorded accesses all sit at offset `(0,0,0)`, and the
padding derived from them is zero.

### Role metadata beside the registry

A Dendro kernel is an ordinary NRPy CFunction plus scheduling metadata keyed by
the registered function name. That sidecar carries no function body, signature,
parameter default, field declaration, or source path, because the CFunction
registry is the only body and signature store. `register_dendro_CFunction`
writes both halves, `get_CFunction_roles` reads the metadata back, and the
registered EVOL, AUXEVOL, and DIAG orders are queried from the gridfunction
registry rather than restated. `set_registration_open` and
`assert_dendro_registration_open` bound the window in which registration is
legal, which is what the freeze boundary later seals.

### Loops

The interior point loop is x-fastest with `i0` innermost over the padded block
interior. It emits the interior base index `pp`, which the gridfunction class's
one-point read helper uses for every read, along with the interior coordinates
`xx0`, `xx1`, and `xx2`. The block loop iterates the local block list supplied
by Dendro and invokes the registered per-block CFunction for each block. Both
are emitted by NRPy through the generic loop helper; the Dendro runtime supplies
only the block list and its count. `require_serial_parallelization` keeps the
emitted loop free of a parallelization directive NRPy is not entitled to choose
for the host.

### Generation parameters

The non-scientific lowering choices are ordinary NRPy parameters registered
through `par.register_param`, with no parallel configuration object: the values
live in the parameter registry and are read back from it through
`get_generation_parameter_view`, and `validate_generation_parameters` rejects an
unsupported combination. Runtime parameter default, validation, and print
CFunctions are registered last, after the scientific CFunctions have established
the complete use closure, and are registered without a Dendro role because the
host context calls them rather than scheduling them as numerical kernels.

## Sources

- [grid.py](../../../nrpy/grid.py) - `DendroGridFunction`, `read_gf_from_memory_Ccode_onept`
- [naming.py](../../../nrpy/infrastructures/Dendro/naming.py) - `input_pointer`, `rhs_pointer`, `out_pointer`, `diag_pointer`, `aux_pointer`, `rhs_symbol_to_gridfunction_name`, `validate_cpp_identifier`
- [access_capture.py](../../../nrpy/infrastructures/Dendro/access_capture.py) - `capture_gridfunction_accesses`, `record_dendro_access`, `required_padding`, `accessed_gridfunction_names`
- [registration.py](../../../nrpy/infrastructures/Dendro/registration.py) - `register_dendro_CFunction`, `get_CFunction_roles`, `registered_evol_order`, `assert_dendro_registration_open`
- [simple_loop.py](../../../nrpy/infrastructures/Dendro/simple_loop.py) - `interior_loop`, `require_serial_parallelization`
- [block_loop.py](../../../nrpy/infrastructures/Dendro/block_loop.py) - `block_loop`
- [generation_parameters.py](../../../nrpy/infrastructures/Dendro/generation_parameters.py) - `validate_generation_parameters`, `get_generation_parameter_view`
- [parameters.py](../../../nrpy/infrastructures/Dendro/runtime/parameters.py) - `register_parameter_CFunctions_last`
- [ADR_dendro_names.md](../../../nrpy/infrastructures/Dendro/ADR_dendro_names.md) - `Decision`, `Consequences`

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
- See also: [Frozen Snapshot And Pure Translators](frozen-snapshot-and-pure-translators.md)
- Example: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
- See also: [Finite Difference](../../core/finite-difference.md)
- See also: [Infrastructure Code Style](../infrastructure-code-style.md)
