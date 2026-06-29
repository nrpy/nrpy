# Aligned-Spin Calibration And Remnant

> Map SEOBNR aligned-spin calibration constants, remnant fits, and NR attachment data. · Status: confirmed · Last reconciled: 2026-06-29
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
zero, and expects the calibration workflow to provide nonspinning inputs. In
spin calibration mode it computes the nonspinning pieces first, then exposes
`dSO` and `Delta_t_S` as symbols. In the default post-calibration mode it calls
`compute_calibration_params()` and stores calibrated expressions such as
`pyseobnr_a6`, `pyseobnr_dSO`, `Delta_t_NS`, and `Delta_t_S`.

For all modes, `Delta_t` is defined as `Delta_t_NS + Delta_t_S`. The class then
computes remnant properties through `final_spin_non_precessing_HBR2016()` and
`final_mass_non_precessing_UIB2016()`, stores the remnant outputs as `a_f` and
`M_f`, evaluates `rISCO` with `Kerr_ISCO_radius(a_f)`, and forms `rstop` with
the branchless coordinate comparison helpers used elsewhere in NRPy equations.

The class also stores NR-fitted attachment data in dictionaries keyed by mode
strings. `hNR` and `omegaNR` cover `(2,2)`, `(3,3)`, `(2,1)`, `(4,4)`, `(4,3)`,
`(5,5)`, and `(3,2)`. The script validation flattens those dictionaries into
stable keys such as `hNR_22`, `omegaNR_22`, `hNR_55`, and `omegaNR_32` before
calling the trusted-expression pipeline.

## Sources

- [SEOBNRv5_aligned_spin_constants.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py) - `SEOBNR_aligned_spin_constants`
- [SEOBNRv5_aligned_spin_constants.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py) - `compute_calibration_params`, `Kerr_ISCO_radius`
- [SEOBNRv5_aligned_spin_constants.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py) - `final_spin_non_precessing_HBR2016`, `final_mass_non_precessing_UIB2016`
- [SEOBNRv5_aligned_spin_constants.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_constants.py) - `trusted_dict`

## See Also

- [SEOBNR And BOB](index.md)
- [Aligned-Spin Hamiltonian](aligned-spin-hamiltonian.md)
- [Aligned-Spin Waveforms](aligned-spin-waveforms.md)
- [Precessing Rotations And Ringdown](precessing-rotations-and-ringdown.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
