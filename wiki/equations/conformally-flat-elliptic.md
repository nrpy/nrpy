# Conformally Flat Elliptic

> Map the NRPyElliptic conformally flat hyperbolic-relaxation equations and source terms. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Equations](index.md)

## Summary

The conformally flat elliptic modules build the symbolic pieces for
NRPyElliptic Bowen-York binary-puncture initial data. One module owns the hyperbolic
relaxation RHSs for the nonsingular correction `uu`; the other constructs the
background conformal factor and Bowen-York extrinsic-curvature source terms.
This is not the BHaH integration of the separate TwoPunctures library.

## Detail

`HyperbolicRelaxationCurvilinearRHSs` is the RHS owner. It selects a reference
metric by coordinate-system name, with an optional `_rfm_precompute` variant,
and builds the same reference-metric Laplacian structure used by the
curvilinear wave equation. It registers `eta_damping`, the `EVOL`
gridfunctions `uu` and `vv`, the `AUXEVOL` gridfunctions
`variable_wavespeed`, `psi_background`, and `ADD_times_AUU`, and the `AUX`
gridfunction `residual_H`.

The relaxation state is stored as `uu_rhs`, `vv_rhs`, and `residual`.
`uu_rhs` damps `uu` while advancing by `vv`; `vv_rhs` combines the
reference-metric Laplacian with the conformally flat Hamiltonian-constraint
source term using `ADD_times_AUU` and `psi_background + uu`. The class records
`residual` before multiplying `vv_rhs` by `variable_wavespeed**2`, so generated
code can separately evaluate the elliptic residual and the hyperbolized
evolution RHS. Stable generated-code names are exposed through
`NRPyElliptic_RHSs_varname_to_expr_dict`.

The expression contains `(psi_background + uu)**(-7)` and the source
construction contains puncture-distance denominators. These builders do not
guard singular evaluation points or select a numerical domain, boundary
conditions, relaxation-wave profile, or stopping criterion; those are caller
and infrastructure responsibilities.

`compute_psi_background_and_ADD_times_AUU` owns the source-term construction
for two punctures on a conformally flat background. It registers bare masses,
the puncture separation parameter `zPunc`, both punctures' linear momenta, and
both punctures' spins. It builds the singular background conformal factor with
`psi_background_cartesian`, superposes the two Bowen-York vector potentials
with `VU_cart_two_punctures`, constructs the conformal tracefree extrinsic
curvature with `ADD_conf_cartesian`, contracts it with
`ADD_times_AUU_conf_cartesian`, and maps Cartesian expressions to the requested
reference-metric coordinates through `replace_cart_coord_by_xx`.

Validation covers both sides of the split. `ConformallyFlat_RHSs.py` generates
trusted RHS and residual dictionaries for Cartesian, SinhCartesian, Spherical,
SinhSpherical, Cylindrical, SinhCylindrical, SymTP, and SinhSymTP. The source
term module generates representative `psi_background_<Coord>` and
`ADD_times_AUU_<Coord>` trusted values for Spherical, SinhSpherical,
Cartesian, SinhCylindrical, and SinhSymTP.

All RHS trusted variants use `enable_rfm_precompute=False`; the implemented
precompute branch is not selected by this module entry point. Trusted matching
checks sampled symbolic values only. It does not show that a generated solver
builds, reaches a steady state, satisfies boundary conditions, or attains any
residual or physical-accuracy threshold. Performance and accuracy claims in the
NRPyElliptic paper describe its reported solver configuration, not evidence for
every current backend or option in this repository.

BHaH `TwoPunctures` generated-C registration and persistent-library lifecycle
belong to [BHaH GR Application Wiring](../infrastructures/bhah/gr-application-wiring.md).

## Sources

- [ConformallyFlat_RHSs.py](../../nrpy/equations/nrpyelliptic/ConformallyFlat_RHSs.py) - `HyperbolicRelaxationCurvilinearRHSs`, `NRPyElliptic_RHSs_varname_to_expr_dict`
- [ConformallyFlat_SourceTerms.py](../../nrpy/equations/nrpyelliptic/ConformallyFlat_SourceTerms.py) - `compute_psi_background_and_ADD_times_AUU`, `psi_background_cartesian`
- [ConformallyFlat_SourceTerms.py](../../nrpy/equations/nrpyelliptic/ConformallyFlat_SourceTerms.py) - `VU_cart_two_punctures`, `ADD_conf_cartesian`, `ADD_times_AUU_conf_cartesian`, `replace_cart_coord_by_xx`
- [ConformallyFlat_RHSs_Cartesian.py](../../nrpy/equations/nrpyelliptic/tests/ConformallyFlat_RHSs_Cartesian.py) - `trusted_dict`
- [ConformallyFlat_RHSs_Spherical.py](../../nrpy/equations/nrpyelliptic/tests/ConformallyFlat_RHSs_Spherical.py) - `trusted_dict`
- [ConformallyFlat_SourceTerms_Cartesian.py](../../nrpy/equations/nrpyelliptic/tests/ConformallyFlat_SourceTerms_Cartesian.py) - `trusted_dict`
- [ConformallyFlat_SourceTerms_Spherical.py](../../nrpy/equations/nrpyelliptic/tests/ConformallyFlat_SourceTerms_Spherical.py) - `trusted_dict`
- [NRPyElliptic paper](https://arxiv.org/abs/2111.02424) - defining scientific background for hyperbolic relaxation and conformally flat binary-puncture initial data

## See Also

- [Equations](index.md)
- [Wave Equation](wave-equation.md)
- [Initial Data](general-relativity/initial-data.md)
- [Geometry And Special-Function Support](geometry-and-special-function-support.md)
- [Trusted Expression Pipeline](trusted-expression-pipeline.md)
- [Reference Metrics](../core/reference-metrics.md)
