# Psi4 And Tetrads

> Map the Psi4 radiation scalar and tetrad construction modules. · Status: confirmed
> Up: [General Relativity](index.md)

## Summary

`Psi4` builds a vacuum-derived reconstruction of the real and imaginary parts
of the Weyl scalar `psi4` from BSSN variables converted to ADM quantities,
metric derivatives, extrinsic-curvature derivatives, and a null tetrad. It
provides no general matter or spacetime-Ricci corrections for nonvacuum use. `Psi4Tetrads` supplies the implemented
Baker-Campanelli-Lousto arXiv:gr-qc/0104063v3 Sec. V.A step (a) tetrad path;
`Psi4` can also leave the tetrad vectors symbolic for validation and
code-generation workflows that provide them externally.

## Detail

`Psi4` starts from `BSSN_to_ADM`, constructs the 3-Riemann tensor, the
extrinsic-curvature terms, and the contractions that produce `psi4_re` and
`psi4_im`. To keep generated expressions manageable, it declares
`gammaDDdDD`, `GammaUDD`, and `KDDdD` as intermediate symbols inside the main
contractions, then records `metric_derivs_varname_list`,
`metric_derivs_varname_arr_list`, `metric_derivs_expr_list`, and related string
helpers for the generated C functions that compute those metric derivative
inputs separately. The curvature reconstruction uses
`R_ij + K*K_ij - K_ip*K^p_j` in its rank-two term. General nonvacuum spacetime
curvature requires additional information absent from this interface; the
result is not guaranteed to be the Weyl scalar for arbitrary matter data.

Evidence (descriptive behavior): `general_relativity/psi4.py::Psi4.__init__`, `rank2term3DD` and curvature contractions.

The `tetrad` argument controls how the null tetrad enters. With
`tetrad="leave_symbolic"`, `Psi4` declares symbolic `mre4U`, `mim4U`, and
`n4U` vectors. With any other value, it delegates to `Psi4Tetrads`; that class
currently accepts only
`tetrad="BCL_arXiv_gr_qc_0104063v3_Eq_5p6_tetrad"` and raises `ValueError`
for unsupported tetrad choices.

`Psi4Tetrads` implements Baker-Campanelli-Lousto v3 Sec. V.A step (a): Eq.
(5.6) forms the null tetrad from the timelike unit normal and an orthonormal
spatial frame; Eq. (5.7) supplies the spatial seed vectors; and the unnumbered
Gram-Schmidt procedure immediately following Eq. (5.7) orthonormalizes them.
The later Kerr-alignment rotation in Eq. (5.9) is intentionally absent.

Evidence (scientific/interface requirement): `arXiv:gr-qc/0104063v3`, Sec. V.A step (a), Eqs. (5.6)-(5.7), and the following unnumbered Gram-Schmidt procedure; corroboration: `general_relativity/tests/psi4_BCL_arXiv_gr_qc_0104063v3_Eq_5p6_tetrad_Spherical.py::trusted_dict`.

The implementation uses the requested reference metric, converts BSSN
variables to ADM variables, builds the Eq. (5.7) Cartesian seed vectors,
transforms them into the active coordinate basis, and applies the cited
Gram-Schmidt procedure with the ADM 3-metric. It then stores `l4U`, `n4U`,
`mre4U`, and `mim4U`. Their slots are components in the mixed
`{n_hat, partial_i}` frame used by the `Psi4` curvature projections: component
zero is along the hypersurface unit normal, while components one through three
are along the spatial coordinate basis. The unit normal therefore has frame
components `u4U = (1, 0, 0, 0)`. `Psi4Tetrads` always uses that representation
and exposes no option to replace it with coordinate-basis lapse and shift
components. The construction requires a nonsingular coordinate chart and
positive-definite spatial metric with nondegenerate Gram-Schmidt seeds. Its
Cartesian seed `(-y,x,0)` vanishes on the z-axis, including the origin, and
normalization divides by zero there; no alternate seed is supplied.
`leave_symbolic` permits a caller to supply a valid external tetrad separately.

Evidence (descriptive behavior): `general_relativity/psi4_tetrads.py::Psi4Tetrads.__init__`; `general_relativity/psi4.py::Psi4`.

Checked-in `trusted_dict` modules cover `Spherical` BCL v3 Eq. (5.6) tetrad, symbolic-tetrad `Psi4`, and direct tetrad construction, plus `SinhSpherical_rfm_precompute`.

## Sources

- [Baker, Campanelli, and Lousto, arXiv:gr-qc/0104063v3](https://arxiv.org/pdf/gr-qc/0104063v3) - Sec. V.A step (a), Eqs. (5.6)-(5.7), and unnumbered Gram-Schmidt procedure
- [Phys. Rev. D 65, 044001](https://doi.org/10.1103/PhysRevD.65.044001) - published Baker, Campanelli, and Lousto reference
- [psi4.py](../../../nrpy/equations/general_relativity/psi4.py) - `Psi4`, `psi4_re`, `psi4_im`, `metric_derivs_varname_list`
- [psi4_tetrads.py](../../../nrpy/equations/general_relativity/psi4_tetrads.py) - `Psi4Tetrads`, `l4U`, `n4U`, `mre4U`, `mim4U`
- [psi4_BCL_arXiv_gr_qc_0104063v3_Eq_5p6_tetrad_Spherical.py](../../../nrpy/equations/general_relativity/tests/psi4_BCL_arXiv_gr_qc_0104063v3_Eq_5p6_tetrad_Spherical.py) - `trusted_dict`
- [psi4_leave_symbolic_Spherical.py](../../../nrpy/equations/general_relativity/tests/psi4_leave_symbolic_Spherical.py) - `trusted_dict`
- [psi4_tetrads_Spherical.py](../../../nrpy/equations/general_relativity/tests/psi4_tetrads_Spherical.py) - `trusted_dict`
- [psi4_BCL_arXiv_gr_qc_0104063v3_Eq_5p6_tetrad_SinhSpherical_rfm_precompute.py](../../../nrpy/equations/general_relativity/tests/psi4_BCL_arXiv_gr_qc_0104063v3_Eq_5p6_tetrad_SinhSpherical_rfm_precompute.py) - `trusted_dict`

## See Also

- See also: [BSSN Family](bssn-family.md)
- Depends on: [Reference Metrics](../../core/reference-metrics.md)
- Validated by: [Expression Validation Helpers](../../validation/expression-validation-helpers.md)
