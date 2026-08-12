# SEOBNR Precessing Dynamics

> Map NRPy's quasi-precessing SEOBNRv5 Hamiltonian and spin/orbital evolution equation builders. · Status: confirmed · Last reconciled: 07-13-2026
> Up: [SEOBNR And BOB](index.md)

## Summary

The precessing-dynamics SEOBNR modules build symbolic expressions for two
adjacent parts of the precessing binary-black-hole model. The quasi-precessing
Hamiltonian module owns conservative EOB quantities, Hamiltonian derivatives,
and circular-orbit derivatives. The spin-evolution module owns the
post-Newtonian right-hand sides for the spin vectors, Newtonian orbital-plane
unit vector, orbital angular momentum, and orbital frequency.

“Quasi-precessing” is material scope. Its relationship to the unversioned
current-latest precessing-dynamics paper page is unaudited model background, not
a versioned intended contract. The implementation orbit-averages in-plane spin
terms under circular-orbit assumptions and sets non-black-hole multipolar
coefficients to zero; audit a specific paper revision before asserting an exact
paper-to-code specialization mapping.

## Detail

`SEOBNRv5_quasi_precessing_spin_Hamiltonian_quantities` declares masses,
canonical variables `r`, `phi`, `prstar`, and `pphi`, vector spin components,
spin projections onto orbital angular momentum, orbit-averaged spin
projections onto Newtonian angular momentum, and calibration symbols `a6` and
`dSO`. It forms normalized mass and spin combinations, applies the
black-hole-binary assumption documented in the constructor, and builds the
quasi-precessing Hamiltonian through nonspinning potentials, even-in-spin
terms, odd-in-spin terms, and spin-orbit calibration terms.

The Hamiltonian object's stable outputs include `xi`, `Hreal`,
`dHreal_dr`, `dHreal_dprstar`, `dHreal_dpphi`, `dHreal_dr_dr`, and
`dHreal_dr_dpphi`. It also stores circular-orbit derivatives obtained by
substituting `prstar = 0`: `dHreal_dr_circ`, `dHreal_dpphi_circ`,
`dHreal_dr_dr_circ`, `dHreal_dr_dpphi_circ`, and
`dHreal_dpphi_dpphi_circ`. These outputs are the symbolic interface for ODE
integration, circular initial-data solves, and waveform or flux consumers that
need instantaneous Hamiltonian derivatives.

`SEOBNRv5_spin_evolution_equations` declares masses, orbital frequency
`omega`, spin-vector components, and the Newtonian angular-momentum unit vector
`ln_x`, `ln_y`, and `ln_z`. It builds SymPy matrices for the spin vectors and
`ln`, derives normalized mass combinations, constructs the PN velocity
`v = omega^(1/3)`, and implements the documented SEOBNRv5 spin and orbital
right-hand sides in the quasi-circular limit.

The spin-evolution object's outputs are componentwise and code-generation
friendly. It stores `chi1_dot_x`, `chi1_dot_y`, `chi1_dot_z`,
`chi2_dot_x`, `chi2_dot_y`, `chi2_dot_z`, `ln_dot_x`, `ln_dot_y`,
`ln_dot_z`, `L_x`, `L_y`, `L_z`, and `omega_dot`. The implementation computes
`omega_dot` by differentiating through `omega = v^3`, so the stored output is
the orbital-frequency derivative rather than only a PN velocity derivative.

Both modules validate by running doctests, processing an instantiated object's
`__dict__` with `process_dictionary_of_expressions`, and comparing the
deterministic result against the sibling trusted file. Representative trusted
keys pin the Hamiltonian outputs `Hreal`, `xi`, and the `dHreal_*` derivatives,
and the spin-evolution outputs `chi*_dot_*`, `ln_dot_*`, `L_*`, and
`omega_dot`.

Those trusted files compare default object state at one deterministic sampled
symbol assignment. They do not test preservation of spin or unit-vector
constraints under integration, compare a trajectory with pySEOBNR, build a
generated implementation, or establish waveform accuracy.

## Sources

- [SEOBNRv5_quasi_precessing_spin_Hamiltonian.py](../../../nrpy/equations/seobnr/SEOBNRv5_quasi_precessing_spin_Hamiltonian.py) - `SEOBNRv5_quasi_precessing_spin_Hamiltonian_quantities`, `Hreal`, `xi`, `dHreal_dr`, `dHreal_dpphi_circ`
- [SEOBNRv5_spin_evolution_equations.py](../../../nrpy/equations/seobnr/SEOBNRv5_spin_evolution_equations.py) - `SEOBNRv5_spin_evolution_equations`, `chi1_dot_x`, `chi2_dot_x`, `ln_dot_x`, `L_x`, `omega_dot`
- [SEOBNRv5_quasi_precessing_spin_Hamiltonian.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_quasi_precessing_spin_Hamiltonian.py) - `trusted_dict`
- [SEOBNRv5_spin_evolution_equations.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_spin_evolution_equations.py) - `trusted_dict`
- [SEOBNRv5 dynamics current latest paper page](https://arxiv.org/abs/2303.18143) - background orientation only; cited sections/equations are not yet audited to a pinned revision

## See Also

- [SEOBNR And BOB](index.md)
- [Aligned-Spin Hamiltonian](aligned-spin-hamiltonian.md)
- [Precessing Rotations And Ringdown](precessing-rotations-and-ringdown.md)
- [BOB Waveforms](bob-waveforms.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
