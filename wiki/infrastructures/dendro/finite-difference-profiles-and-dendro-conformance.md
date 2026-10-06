# Finite-Difference Profiles And Dendro Conformance

> Define supported centered finite-difference and Kreiss-Oliger profiles. · Status: confirmed
> Up: [Dendro](index.md)

## Summary

Dendro generation emits centered regular orders 4, 6, and 8 with KO orders 2,
4, and 6. Their stencil radii are 2, 3, and 4 points. No upwind derivative is
generated.

## Detail

| Regular order | KO finite-difference order | KO derivative order | Required padding |
| --- | --- | --- | --- |
| 4 | 2 | 4 | 2 |
| 6 | 4 | 6 | 3 |
| 8 | 6 | 8 | 4 |

One generated application contains every profile. Runtime dispatch accepts
Dendro element order 4, 6, or 8, selects the matching Ricci, RHS, constraint,
and derivative-dependent initial-data kernels, and requires padding at least
`order/2`. A mismatched KO order or unsupported element order is a configuration
error.

Kreiss-Oliger dissipation is enabled in the production profile. The generated
default uses order 6 with KO order 4. Validation must exercise each profile in
the complete application; source presence alone does not establish correct
dispatch or block indexing.

The example option `--fd-order` (4, 6, or 8; default 6) selects only the defaults
of the generated application: the constants `FD_ORDER`, `KO_FD_ORDER`, and
`REQUIRED_PADDING`, the `BSSN_ELE_ORDER` line of the generated `pars/<stem>.toml`,
and the fallback of `BSSN_ELE_ORDER` in the executable. All three orders are
emitted whatever the option. On a fresh start the order used at run time is the
`BSSN_ELE_ORDER` value of the parameter file, so a file that sets the key
overrides the option: the packaged `pars/q1.par.lowres.toml` sets
`BSSN_ELE_ORDER = 6`, and the printed run commands run FD6 whatever `--fd-order`
was. A restore requires `BSSN_ELE_ORDER` to equal the order stored in the
checkpoint and is rejected otherwise.

Claim evidence:
- Claim: `--fd-order` sets only the default order (the constants `FD_ORDER`, `KO_FD_ORDER`, `REQUIRED_PADDING`, the `BSSN_ELE_ORDER` line of the generated parameter file, and the fallback of the `BSSN_ELE_ORDER` read); all orders 4, 6, and 8 are emitted, and the order of a run is the `BSSN_ELE_ORDER` value in its parameter file, which the packaged q1 file sets to 6; a restore is rejected unless that value equals the element order stored in the checkpoint.
- Role: descriptive behavior
- Deciding authority: `nrpy/examples/dendro_bssn.py` and `nrpy/examples/dendro_fccz4.py`, `parse_args` and `main`; `nrpy/infrastructures/Dendro/param_toml.py`, `generate_default_parfile`; `nrpy/infrastructures/Dendro/constants_h.py`, `output_constants_h`; `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp` (the `BSSN_ELE_ORDER` read); `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::restore_checkpt` within `output_solver_context_cpp`, and `nrpy/infrastructures/Dendro/checkpoint.py`, `output_checkpoint_cpp` (the stored element order and the comparison with the running order).
- Corroboration: `nrpy/examples/q1.par.lowres.toml`, the `BSSN_ELE_ORDER` line; `nrpy/examples/tests/dendro_application_check.py`, `Leg.generate_and_build` and `Leg.run_variant`, which generate at FD6 and select other orders through `BSSN_ELE_ORDER`.

Every generated right-hand side is centered. `register_CFunction_rhs_eval`
replaces each directional derivative symbol of the equation modules (the `dupD`
and `ddnD` derivatives that carry shift advection) by the centered `dD`
derivative of the selected order, so advection terms use centered stencils and
no upwind or downwind stencil appears in the emitted kernels.

Claim evidence:
- Claim: The generated Dendro right-hand-side kernels replace every `dupD` and `ddnD` directional derivative of the equation modules by the centered `dD` derivative of the selected finite-difference order, so no upwind or downwind derivative is emitted.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`, `register_CFunction_rhs_eval`, the centered advection substitution.
- Corroboration: none available: no test compares the advection stencils; `nrpy/examples/tests/dendro_application_check.py`, `Leg.check_reference`, regresses only the evolved result.

## Sources

- [finite_difference.py](../../../nrpy/finite_difference.py) - finite-difference and KO stencil construction.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - order-specific RHS registration and the centered advection substitution.
- [constants_h.py](../../../nrpy/infrastructures/Dendro/constants_h.py) - emitted order constants.
- [param_toml.py](../../../nrpy/infrastructures/Dendro/param_toml.py) - the order line of the generated parameter file.
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - the run-time order read.
- [q1.par.lowres.toml](../../../nrpy/examples/q1.par.lowres.toml) - packaged parameter file with its own order.
- [Ricci_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py) - order-specific Ricci registration.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - runtime order dispatch and padding checks.
- [Dendrolib block.h](https://github.com/paralab/Dendro-5.01/blob/master/include/block.h) - regular block geometry.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Finite Difference](../../core/finite-difference.md)
- See also: [Runtime Parameter Keys](runtime-parameters.md)
- Validated by: [Production Validation And Deferred Checks](validation-standalone-host-and-deferral-gates.md)
