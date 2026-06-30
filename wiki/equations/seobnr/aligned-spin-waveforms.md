# Aligned-Spin Waveforms

> Map SEOBNRv5 aligned-spin factorized modes, flux, strain, and merger attachment quantities. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [SEOBNR And BOB](index.md)

## Summary

This page owns the aligned-spin SEOBNR waveform side: factorized-resummed mode
ingredients, flux construction, strain-mode construction, and the native
`(l=2,m=2)` merger-ringdown attachment quantities. It depends on conservative
Hamiltonian outputs and calibration/remnant quantities but does not own those
fits.

## Detail

`SEOBNRv5_aligned_spin_waveform_quantities` declares mass, spin, canonical,
phase, Hamiltonian, and orbital-frequency symbols, including `Hreal`, `Omega`,
and `Omega_circ`. It derives mass-ratio and spin combinations such as `nu`,
`delta`, `chi_A`, `chi_S`, `vomega`, `vphi`, `vh3`, and the effective-source
list. It also handles equal-mass limits through branchless `noneqcond`,
`eqcond`, `deltainvertible`, and `deltainv` expressions.

The waveform class stores factorized-mode ingredients in dictionaries and
matrices. The `rho` dictionary covers modes with `2 <= l <= 8` and
`1 <= m <= l`; `deltalm`, `fspin`, and `fspin_limit` cover the modeled
low-order aligned-spin modes; `pn_contribution_f` stores special PN amplitude
contributions for `(2,1)`, `(4,3)`, and `(5,5)`. Optional special amplitude
coefficients are controlled by `apply_special_amplitude_coefficients`, which
keeps `c_43`, `c_21`, and `c_55` symbolic when enabled and sets them to zero by
default.

`flux()` loops over modes through `l = 8`, combines Newtonian amplitudes,
effective sources, tail terms, and PN factors, and returns the factorized
resummed flux divided by `nu`. `strain()` returns an `hlms` dictionary for the
implemented strain modes `(2,2)`, `(2,1)`, `(3,3)`, `(3,2)`, `(4,4)`, `(4,3)`,
and `(5,5)`.

`SEOBNRv5_aligned_spin_merger_quantities` owns the native aligned-spin
merger-ringdown attachment surface for the `(2,2)` mode. It builds symbolic
amplitude and phase outputs `h` and `phi`, then exposes attachment quantities
`h_t_attach`, `hdot_t_attach`, `hddot_t_attach`, `w_t_attach`, and
`wdot_t_attach` from NR-fitted expressions and the quasi-normal-mode inputs.

Validation is split by module. The waveform-quantity script validates the class
state and flattens the mode dictionaries into keys such as `rho(2 , 2)`,
`fspin(2 , 1)`, `fspin_limit(5 , 5)`, `deltalm(4 , 3)`, and
`pn_contribution_f(5 , 5)`. The merger-quantity script validates the merger
class `__dict__`, including the attachment keys and the final `h` and `phi`
expressions.

## Sources

- [SEOBNRv5_aligned_spin_waveform_quantities.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_waveform_quantities.py) - `SEOBNRv5_aligned_spin_waveform_quantities`
- [SEOBNRv5_aligned_spin_waveform_quantities.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_waveform_quantities.py) - `flux`, `strain`
- [SEOBNRv5_aligned_spin_merger_quantities.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_merger_quantities.py) - `SEOBNRv5_aligned_spin_merger_quantities`
- [SEOBNRv5_aligned_spin_waveform_quantities.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_waveform_quantities.py) - `trusted_dict`
- [SEOBNRv5_aligned_spin_merger_quantities.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_merger_quantities.py) - `trusted_dict`

## See Also

- [SEOBNR And BOB](index.md)
- [Aligned-Spin Hamiltonian](aligned-spin-hamiltonian.md)
- [Aligned-Spin Calibration And Remnant](aligned-spin-calibration-and-remnant.md)
- [BOB Waveforms](bob-waveforms.md)
- [Precessing Rotations And Ringdown](precessing-rotations-and-ringdown.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
