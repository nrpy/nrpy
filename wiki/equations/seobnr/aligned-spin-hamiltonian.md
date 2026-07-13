# Aligned-Spin Hamiltonian

> Map the conservative SEOBNRv5 aligned-spin Hamiltonian quantities and derivative outputs. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [SEOBNR And BOB](index.md)

## Summary

The aligned-spin Hamiltonian page owns the conservative Effective-One-Body
dynamics implemented for SEOBNRv5 aligned-spin binaries. The source builds
symbolic expressions for the tortoise map, effective and real Hamiltonian
pieces, Hamiltonian derivatives, circular-orbit derivatives, and frequency-like
quantities used by the waveform and orbital-evolution consumers.

Scope is the source's current SEOBNRv5 aligned-spin implementation. The
SEOBNRv5HM paper defines the intended quasi-circular, spinning,
nonprecessing-binary-black-hole model; this page does not generalize these
expressions to eccentric, precessing, or matter systems.

## Detail

`SEOBNRv5_aligned_spin_Hamiltonian_quantities` declares the mass, aligned-spin,
radius, tortoise radial momentum, angular momentum, and calibration symbols
`m1`, `m2`, `chi1`, `chi2`, `r`, `prstar`, `pphi`, `a6`, and `dSO`. It stores
the symmetric mass ratio as `nu`, constructs the aligned-spin nonspinning,
spin-spin, and spin-orbit pieces, and exposes stable object outputs including
`xi`, `Aalign`, `Balignnp`, `Heven`, and `Hreal`.

The same class differentiates the real Hamiltonian into quantities used for
integration and initial-data construction. Public derivative outputs include
`dHreal_dr`, `dHreal_dprstar`, `dHreal_dpphi`, `dHreal_dr_dr`,
`dHreal_dr_dpphi`, `dQalign_dprstar`, `dHoddbar_dr`, and `r_dot`. Circular
counterparts are formed by substituting `prstar = 0` before differentiating, with
outputs such as `dHreal_dr_circ`, `dHreal_dpphi_circ`,
`dHreal_dr_dr_circ`, `dHreal_dr_dpphi_circ`, and
`dHreal_dpphi_dpphi_circ`.

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
- [validate_expressions.py](../../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `compare_or_generate_trusted_results`
- [SEOBNRv5_aligned_spin_Hamiltonian.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_aligned_spin_Hamiltonian.py) - `trusted_dict`
- [SEOBNRv5HM paper](https://arxiv.org/abs/2303.18039) - intended aligned-spin model scope and equations

## See Also

- [SEOBNR And BOB](index.md)
- [Aligned-Spin Calibration And Remnant](aligned-spin-calibration-and-remnant.md)
- [Aligned-Spin Waveforms](aligned-spin-waveforms.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
