# Aligned-Spin Waveforms

> Map SEOBNRv5 aligned-spin factorized modes, flux, strain, and merger attachment quantities. · Status: confirmed
> Up: [SEOBNR And BOB](index.md)

## Summary

The aligned-spin SEOBNR waveform modules build factorized-resummed mode
ingredients, radiation-reaction flux, strain modes, and native `(l=2,m=2)`
merger-ringdown attachment quantities from conservative Hamiltonian outputs and
calibration/remnant fits.

These are current symbolic ingredients associated with the quasi-circular,
spinning, nonprecessing binary-black-hole SEOBNRv5HM model. The unversioned paper
landing page is model background; model-wide calibration and accuracy claims
neither decide this transcription's
equations nor validate every generated C caller.

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

The `(3,3)` mode has an imaginary spin correction, `fspinimag`, with the
equal-mass limit `fspinimag_limit`. Both `strain()` and `flux()` select between
them with `noneqcond` and `eqcond`, as for the real corrections. In `strain()`
the `(3,3)` PN factor is the real correction plus `I * fspinimag`. In `flux()`
the `(3,3)` term adds `fspinimag**2` to the square of the real correction,
which is the squared magnitude of that complex factor. No other mode has an
imaginary spin correction.

Evidence (descriptive behavior): `seobnr/SEOBNRv5_aligned_spin_waveform_quantities.py::flux`, `strain`; corroboration: `seobnr/tests/SEOBNRv5_aligned_spin_waveform_quantities.py::trusted_dict`, pins the sampled `fspinimag`/`fspinimag_limit` symbolic values.

`SEOBNRv5_aligned_spin_merger_quantities` constructs the native aligned-spin
merger-ringdown attachment expressions for the `(2,2)` mode. It builds symbolic
amplitude and phase outputs `h` and `phi`, then exposes attachment quantities
`h_t_attach`, `hdot_t_attach`, `hddot_t_attach`, `w_t_attach`, and
`wdot_t_attach` from NR-fitted expressions and the quasi-normal-mode inputs.

Validation is split by module. The waveform-quantity script validates the class
state and flattens the mode dictionaries into keys such as `rho(2 , 2)`,
`fspin(2 , 1)`, `fspin_limit(5 , 5)`, `deltalm(4 , 3)`, and
`pn_contribution_f(5 , 5)`. The merger-quantity script validates the merger
class `__dict__`, including the attachment keys and the final `h` and `phi`
expressions.

The waveform trusted entry point constructs
`SEOBNRv5_aligned_spin_waveform_quantities()` with
`apply_special_amplitude_coefficients=False`. It does not call `flux()` or
`strain()`, so their returned assembled flux and `hlms` dictionaries are not
direct trusted keys. The enabled special-coefficient branch is also absent from
that trusted variant. Trusted matching here proves only sampled values for the
exported default object state; it does not reproduce paper-level calibration,
build generated code, integrate dynamics, or measure waveform accuracy.

## Sources

- [SEOBNRv5_aligned_spin_waveform_quantities.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_waveform_quantities.py) - `SEOBNRv5_aligned_spin_waveform_quantities`; `flux`, `strain`
- [SEOBNRv5_aligned_spin_merger_quantities.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_merger_quantities.py) - `SEOBNRv5_aligned_spin_merger_quantities`
- [SEOBNRv5_aligned_spin_waveform_quantities.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_waveform_quantities.py) - `trusted_dict`
- [SEOBNRv5_aligned_spin_merger_quantities.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_merger_quantities.py) - `trusted_dict`
- [SEOBNRv5HM paper landing page](https://arxiv.org/abs/2303.18039) - background orientation only; aligned-spin waveform/calibration claims have no exact revision mapping here

## See Also

- See also: [Aligned-Spin Hamiltonian](aligned-spin-hamiltonian.md)
- See also: [Aligned-Spin Calibration And Remnant](aligned-spin-calibration-and-remnant.md)
- Validated by: [Expression Validation Helpers](../../validation/expression-validation-helpers.md)
