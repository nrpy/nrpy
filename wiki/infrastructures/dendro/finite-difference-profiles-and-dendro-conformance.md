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

One generated application contains every profile. At run time the executable
accepts only the `BSSN_ELE_ORDER` values 4, 6, and 8, and the solver context
dispatches on the element order to the matching kernel for the `ADM_to_BSSN`
conversion, the connection initialization, Ricci, the RHS, the constraints,
Psi4, and the ADM surface data; an unsupported order throws. The Ricci, RHS,
constraint, connection-initialization, Psi4, and ADM surface data kernels require
padding of at least `order/2` and throw otherwise; the pointwise `ADM_to_BSSN`
conversion checks no padding. The relations `KO_FD_ORDER = FD_ORDER - 2` and
`REQUIRED_PADDING = FD_ORDER/2` are checked when the project is generated, where
a violation raises `ValueError`. The solver reads neither `KO_FD_ORDER` nor
`REQUIRED_PADDING` at run time, reads `FD_ORDER` only as the fallback of
`BSSN_ELE_ORDER`, and prints `KO=` as `element_order - 2`.

Evidence (descriptive behavior): `Dendro/main_cpp.py::output_main_cpp`; `Dendro/solver_context.py::output_solver_context_cpp`; corroboration: `tests/dendro_application_check.py::Leg.run_negatives` (`BSSN_ELE_ORDER = 5` rejection).

Kreiss-Oliger dissipation is enabled in the production profile. The generated
default uses order 6 with KO order 4. Validation must exercise each profile in
the complete application; source presence alone does not establish correct
dispatch or block indexing.

The example option `--fd-order` (4, 6, or 8; default 6) selects only the defaults
of the generated application: the constants `FD_ORDER`, `KO_FD_ORDER`, and
`REQUIRED_PADDING`, the `BSSN_ELE_ORDER` line of the generated `pars/<stem>.toml`,
the fallback of `BSSN_ELE_ORDER` in the executable, and the `NGHOSTS` macro of the
packaged `BHaH_defines.h`, which no Dendro source uses. All three orders are
emitted whatever the option. On a fresh start the order used at run time is the
`BSSN_ELE_ORDER` value of the parameter file, so a file that sets the key
overrides the option: the packaged `pars/q1.par.lowres.toml` sets
`BSSN_ELE_ORDER = 6`, and the printed run commands run FD6 whatever `--fd-order`
was. A restore requires `BSSN_ELE_ORDER` to equal the order stored in the
checkpoint and is rejected otherwise.

Evidence (descriptive behavior): `dendro_bssn.py` and `dendro_fccz4.py::parse_args` and `main`; `Dendro/param_toml.py::generate_default_parfile`; `Dendro/constants_h.py::output_constants_h`; corroboration: `q1.par.lowres.toml`, the `BSSN_ELE_ORDER` line.

Every generated right-hand side is centered. `register_CFunction_rhs_eval`
replaces each directional derivative symbol of the equation modules (the `dupD`
and `ddnD` derivatives that carry shift advection) by the centered `dD`
derivative of the selected order, so advection terms use centered stencils and
no upwind or downwind stencil appears in the emitted kernels.

Evidence (descriptive behavior): `Dendro/general_relativity/rhs_eval.py::register_CFunction_rhs_eval`, the centered advection substitution.

## Sources

- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - order-specific RHS registration and the centered advection substitution
- [Ricci_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py) - order-specific Ricci registration and padding validation
- [constants_h.py](../../../nrpy/infrastructures/Dendro/constants_h.py) - emitted order constants
- [param_toml.py](../../../nrpy/infrastructures/Dendro/param_toml.py) - the order line of the generated parameter file
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - the run-time order read
- [q1.par.lowres.toml](../../../nrpy/examples/q1.par.lowres.toml) - packaged parameter file with its own order
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - runtime order dispatch and padding checks
- [Dendrolib block.h](https://github.com/paralab/Dendro-5.01/blob/master/include/block.h) - regular block geometry
- [checkpoint.py](../../../nrpy/infrastructures/Dendro/checkpoint.py) - the stored element order and its comparison at restore
- [BHaH_defines_h.py](../../../nrpy/infrastructures/BHaH/BHaH_defines_h.py) - the `NGHOSTS` macro of the packaged `BHaH_defines.h`
- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - `--fd-order` handling
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - `--fd-order` handling
- [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py) - the profiles that generate at FD6 and select other orders at run time

## See Also

- Depends on: [Finite Difference](../../core/finite-difference.md)
- See also: [Runtime Parameter Keys](runtime-parameters.md)
- Validated by: [Production Validation And Deferred Checks](production-validation-and-deferred-checks.md)
