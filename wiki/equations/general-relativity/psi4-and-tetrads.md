# Psi4 And Tetrads

> Map the Psi4 radiation scalar and tetrad construction modules. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [General Relativity](index.md)

## Summary

`Psi4` builds the real and imaginary parts of the Weyl scalar `psi4` from BSSN
variables converted to ADM quantities, metric derivatives, extrinsic-curvature
derivatives, and a null tetrad. `Psi4Tetrads` supplies the implemented
quasi-Kinnersley tetrad path; `Psi4` can also leave the tetrad vectors symbolic
for validation and code-generation workflows that provide them externally.

## Detail

`Psi4` starts from `BSSN_to_ADM`, constructs the 3-Riemann tensor, the
extrinsic-curvature terms, and the contractions that produce `psi4_re` and
`psi4_im`. To keep generated expressions manageable, it declares
`gammaDDdDD`, `GammaUDD`, and `KDDdD` as intermediate symbols inside the main
contractions, then records `metric_derivs_varname_list`,
`metric_derivs_varname_arr_list`, `metric_derivs_expr_list`, and related string
helpers for the generated C functions that compute those metric derivative
inputs separately.

The `tetrad` argument controls how the null tetrad enters. With
`tetrad="leave_symbolic"`, `Psi4` declares symbolic `mre4U`, `mim4U`, and
`n4U` vectors. With any other value, it delegates to `Psi4Tetrads`; that class
currently accepts only `tetrad="quasiKinnersley"` and raises `ValueError` for
unsupported tetrad choices.

`Psi4Tetrads` uses the requested reference metric, converts BSSN variables to
ADM variables, builds two Cartesian seed vectors and a third vector from the
Levi-Civita construction, transforms the seed vectors into the active
coordinate basis, and orthonormalizes them with the ADM 3-metric. It then
stores the four-vectors `l4U`, `n4U`, `mre4U`, and `mim4U`. By default the
timelike unit-normal component is simplified to `u4U[0] = 1`; the optional
`use_metric_to_construct_unit_normal` branch instead uses `alpha` and `betaU`
from `BSSN_quantities`.

The cited validation evidence covers representative `Spherical` trusted files
for quasi-Kinnersley Psi4, symbolic-tetrad Psi4, and direct tetrad construction,
plus one quasi-Kinnersley `SinhSpherical_rfm_precompute` representative. The
page does not enumerate every coordinate variant as direct source evidence.

## Sources

- [psi4.py](../../../nrpy/equations/general_relativity/psi4.py) - `Psi4`, `psi4_re`, `psi4_im`, `metric_derivs_varname_list`
- [psi4_tetrads.py](../../../nrpy/equations/general_relativity/psi4_tetrads.py) - `Psi4Tetrads`, `l4U`, `n4U`, `mre4U`, `mim4U`
- [psi4_quasiKinnersley_Spherical.py](../../../nrpy/equations/general_relativity/tests/psi4_quasiKinnersley_Spherical.py) - `trusted_dict`
- [psi4_leave_symbolic_Spherical.py](../../../nrpy/equations/general_relativity/tests/psi4_leave_symbolic_Spherical.py) - `trusted_dict`
- [psi4_tetrads_Spherical.py](../../../nrpy/equations/general_relativity/tests/psi4_tetrads_Spherical.py) - `trusted_dict`
- [psi4_quasiKinnersley_SinhSpherical_rfm_precompute.py](../../../nrpy/equations/general_relativity/tests/psi4_quasiKinnersley_SinhSpherical_rfm_precompute.py) - `trusted_dict`

## See Also

- [General Relativity](index.md)
- [BSSN Family](bssn-family.md)
- [Metric Conversions And Matter](metric-conversions-and-matter.md)
- [Geodesics](geodesics.md)
- [Reference Metrics](../../core/reference-metrics.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
