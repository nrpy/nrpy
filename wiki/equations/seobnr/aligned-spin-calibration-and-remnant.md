# Aligned-Spin Calibration And Remnant

> Map SEOBNR aligned-spin calibration constants, remnant fits, and NR attachment data. · Status: confirmed · Last reconciled: 07-13-2026
> Up: [SEOBNR And BOB](index.md)

## Summary

This page owns the aligned-spin SEOBNR constants that sit between conservative
dynamics and waveform attachment. The implementation computes or exposes
calibration parameters, final remnant mass and spin, Kerr ISCO radius, stopping
radius, and NR-fitted mode amplitudes and frequencies at the waveform matching
time.

## Detail

`SEOBNR_aligned_spin_constants` accepts two mutually exclusive calibration
flags: `calibration_no_spin` and `calibration_spin`. In nonspinning calibration
mode it exposes `a6` and `Delta_t_NS` as symbols, sets `dSO` and `Delta_t_S` to
zero, and expects the calibration workflow to provide `chi1=chi2=0`. The class
does not enforce that spin condition. In spin calibration mode it computes the
nonspinning pieces first, then exposes
`dSO` and `Delta_t_S` as symbols. In the default post-calibration mode it calls
`compute_calibration_params()` and stores calibrated expressions such as
`pyseobnr_a6`, `pyseobnr_dSO`, `Delta_t_NS`, and `Delta_t_S`.

For all modes, `Delta_t` is defined as `Delta_t_NS + Delta_t_S`. The class then
computes remnant properties through `final_spin_non_precessing_HBR2016()` and
`final_mass_non_precessing_UIB2016()`, stores the remnant outputs as `a_f` and
`M_f`, evaluates `rISCO` with `Kerr_ISCO_radius(a_f)`, and forms `rstop` with
the branchless coordinate comparison helpers used elsewhere in NRPy equations.
For negative `Delta_t` the intended selected value is `-1`; for positive
`Delta_t` it is `0.98*rISCO`; at exactly zero both strict-mask terms vanish and
the expression returns zero.

The class also stores NR-fitted attachment data in dictionaries keyed by mode
strings. `hNR` and `omegaNR` cover `(2,2)`, `(3,3)`, `(2,1)`, `(4,4)`, `(4,3)`,
`(5,5)`, and `(3,2)`. The script validation flattens those dictionaries into
stable keys such as `hNR_22`, `omegaNR_22`, `hNR_55`, and `omegaNR_32` before
calling the trusted-expression pipeline.

Script validation instantiates only the default post-calibration mode. Neither
`calibration_no_spin=True` nor `calibration_spin=True` has a sibling trusted
variant here. The stored dictionary is sampled numerical evidence for the
current formulas, not an independent reproduction of the SEOBNRv5HM, HBR2016,
or UIB2016 calibration data and not a remnant-fit accuracy test.

## Sources

- [SEOBNRv5_aligned_spin_constants.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py) - `SEOBNR_aligned_spin_constants`
- [SEOBNRv5_aligned_spin_constants.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py) - `compute_calibration_params`, `Kerr_ISCO_radius`
- [SEOBNRv5_aligned_spin_constants.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py) - `final_spin_non_precessing_HBR2016`, `final_mass_non_precessing_UIB2016`
- [SEOBNRv5_aligned_spin_constants.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_constants.py) - `trusted_dict`
- [SEOBNRv5HM current latest paper page](https://arxiv.org/abs/2303.18039) - background orientation only; Equations 78-81 and calibration context are not yet audited to a pinned revision
- [HBR2016 current latest paper page](https://arxiv.org/abs/1605.01938) - background orientation only; final-spin fit mapping is not yet audited to a pinned revision
- [UIB2016 current latest paper page](https://arxiv.org/abs/1611.00332) - background orientation only; final-state fit and ancillary mapping are not yet audited to a pinned revision

## See Also

- [SEOBNR And BOB](index.md)
- [Aligned-Spin Hamiltonian](aligned-spin-hamiltonian.md)
- [Aligned-Spin Waveforms](aligned-spin-waveforms.md)
- [Precessing Rotations And Ringdown](precessing-rotations-and-ringdown.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
