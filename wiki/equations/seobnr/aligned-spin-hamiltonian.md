# Aligned-Spin Hamiltonian

> Map the conservative SEOBNRv5 aligned-spin Hamiltonian quantities and derivative outputs. · Status: confirmed
> Up: [SEOBNR And BOB](index.md)

## Summary

The SEOBNRv5 aligned-spin Hamiltonian modules implement conservative
Effective-One-Body dynamics for quasi-circular, spinning, nonprecessing binary
black holes, building symbolic expressions for the tortoise map, effective and
real Hamiltonian pieces, Hamiltonian derivatives, circular-orbit derivatives,
and orbital frequencies without generalizing to eccentric, precessing, or matter
systems. The unversioned SEOBNRv5HM paper landing page provides model background
only until an immutable publication locator is registered.

## Detail

`SEOBNRv5_aligned_spin_Hamiltonian_quantities` declares the mass, aligned-spin,
radius, tortoise radial momentum, angular momentum, and calibration symbols
`m1`, `m2`, `chi1`, `chi2`, `r`, `prstar`, `pphi`, `a6`, and `dSO`. It stores
the symmetric mass ratio as `nu`, constructs the aligned-spin nonspinning,
spin-spin, and spin-orbit pieces, and exposes stable object outputs including
`xi`, `Aalign`, `Balignnp`, `Heven`, and `Hreal`.

The same class differentiates the real Hamiltonian into quantities used for
integration and initial-data construction. Public non-circular derivative
outputs are `dHreal_dr`, `dHreal_dprstar`, `dHreal_dpphi`, `dHreal_dr_dr`,
`dHreal_dr_dpphi`, `dQalign_dprstar`, and `dHoddbar_dr`.

Evidence (descriptive behavior): `seobnr/SEOBNRv5_aligned_spin_Hamiltonian.py::SEOBNRv5_aligned_spin_Hamiltonian_quantities.__init__`; corroboration: `seobnr/tests/SEOBNRv5_aligned_spin_Hamiltonian.py::trusted_dict`, pins the exact key set produced by the module's own default-constructor evaluation.

Circular counterparts are formed by substituting `prstar = 0` before
differentiating, with outputs `dHreal_dr_circ`, `dHreal_dpphi_circ`,
`dHreal_dr_dr_circ`, `dHreal_dr_dpphi_circ`, and
`dHreal_dpphi_dpphi_circ`. The circular first derivatives are divided by `nu`;
each circular second derivative differentiates one of those normalized first
derivatives without a second division by `nu`. The two first derivatives, with
the input orbital frequency subtracted from `dHreal_dpphi_circ`, form the
residual of the circular-orbit root solve, and the three second derivatives are
the entries of its Jacobian.

Evidence (descriptive behavior): `seobnr/SEOBNRv5_aligned_spin_Hamiltonian.py::SEOBNRv5_aligned_spin_Hamiltonian_quantities.__init__`; corroboration: `seobnr/tests/SEOBNRv5_aligned_spin_Hamiltonian.py::trusted_dict`.

Every first derivative of `Hreal` is divided by `nu`. The non-circular second
derivatives `dHreal_dr_dr` and `dHreal_dr_dpphi` differentiate the normalized
`dHreal_dr` and divide by `nu` again, so they carry `1/nu**2` where the circular
second derivatives carry `1/nu`; at `prstar = 0`, `dHreal_dr_dr` equals
`dHreal_dr_dr_circ / nu`. The BHaH radial-momentum initial condition uses only
their ratio, in which the extra factor cancels; see [SEOBNR BOB Generated Library](../../infrastructures/bhah/seobnr-bob-generated-library.md).

Evidence (descriptive behavior): `seobnr/SEOBNRv5_aligned_spin_Hamiltonian.py::SEOBNRv5_aligned_spin_Hamiltonian_quantities.__init__`.

Validation is module-local. Running the Hamiltonian source as a script executes
doctests, processes the class `__dict__` through the trusted-expression helper,
and compares the resulting expression dictionary against the sibling trusted
file. The trusted dictionary pins the main conservative outputs, derivative
outputs, and canonical input symbols that survive expression processing.
This is a default-constructor sampled numerical fingerprint, not a formal
derivation check, ODE integration, waveform validation, generated build, or
runtime-accuracy result.

## Sources

- [SEOBNRv5_aligned_spin_Hamiltonian.py](../../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_Hamiltonian.py) - `SEOBNRv5_aligned_spin_Hamiltonian_quantities`
- [SEOBNRv5_aligned_spin_Hamiltonian.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_Hamiltonian.py) - `trusted_dict`
- [SEOBNRv5HM paper landing page](https://arxiv.org/abs/2303.18039) - background orientation only; aligned-spin model equations have no exact revision mapping here

## See Also

- See also: [Aligned-Spin Calibration And Remnant](aligned-spin-calibration-and-remnant.md)
- See also: [Aligned-Spin Waveforms](aligned-spin-waveforms.md)
- Validated by: [Expression Validation Helpers](../../validation/expression-validation-helpers.md)
