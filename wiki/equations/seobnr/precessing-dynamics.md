# SEOBNR Precessing Dynamics

> Map NRPy's quasi-precessing SEOBNRv5 Hamiltonian and spin/orbital evolution equation builders. · Status: confirmed
> Up: [SEOBNR And BOB](index.md)

## Summary

The precessing-dynamics SEOBNR modules build symbolic expressions for two
adjacent parts of the precessing binary-black-hole model. The quasi-precessing
Hamiltonian module owns conservative EOB quantities, Hamiltonian derivatives,
and circular-orbit derivatives. The spin-evolution module owns the
post-Newtonian right-hand sides for the spin vectors, Newtonian orbital-plane
unit vector, orbital angular momentum, and orbital frequency.

“Quasi-precessing” is material scope. The unversioned precessing-dynamics paper
landing page is model background, not a version-pinned statement of the
implemented equations. The
implementation orbit-averages in-plane spin
terms under circular-orbit assumptions and sets non-black-hole multipolar
coefficients to zero. The spin-evolution right-hand sides are mapped below to
the version-pinned arXiv:2303.18143v2; the Hamiltonian module has no exact
paper-revision mapping here, so audit a specific revision before asserting one.

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
`dHreal_dpphi_dpphi_circ`. As in the aligned-spin Hamiltonian, the circular first
derivatives are divided by `nu`, and each circular second derivative
differentiates a normalized first derivative without a second division by `nu`.
These outputs are the symbolic interface for ODE
integration, circular initial-data solves, and waveform or flux consumers that
need instantaneous Hamiltonian derivatives.

Claim evidence:
- Claim: `SEOBNRv5_quasi_precessing_spin_Hamiltonian_quantities`'s circular
  first Hamiltonian derivatives (`dHreal_dr_circ`, `dHreal_dpphi_circ`) are
  divided by `nu`, and each circular second derivative differentiates one of
  those already-normalized first derivatives without a second division by
  `nu`.
- Role: descriptive behavior
- Deciding authority: [SEOBNRv5_quasi_precessing_spin_Hamiltonian.py](../../../nrpy/equations/seobnr/SEOBNRv5_quasi_precessing_spin_Hamiltonian.py), `SEOBNRv5_quasi_precessing_spin_Hamiltonian_quantities.__init__`
- Corroboration: [SEOBNRv5_quasi_precessing_spin_Hamiltonian.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_quasi_precessing_spin_Hamiltonian.py), `trusted_dict`, pins the sampled circular-derivative outputs against the module's own default-constructor evaluation

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

In arXiv:2303.18143v2, the spin derivatives follow Eq. (66), `L` is Eq. (65)
divided by `M mu` (the module works in `M = 1` units, where `mu = nu`),
`ln_dot` follows Eqs. (71a)-(71d), and `omega_dot` is `3 v^2` times `v_dot`
from Eq. (69). The terms proportional to `Ctilde_1ES2 = C_1ES2 - 1` in
Eqs. (65f), (66b), (69f), and (71e) vanish for black holes and are omitted. The
radiation-reaction contribution to `ln_dot` in Eq. (77) is not included.

`ln_dot` is not orthogonal to `ln`. For a unit `ln`, the two `delta`
spin1-spin2 terms of Eq. (71c) give
`ln . ln_dot = nu ln . (S1 x S2) / (M^3 mu^2) [-(delta/4) v^9 + (delta/96)(35 nu - 54) v^11]`,
which is Eq. (72) of the paper. Integrating these right-hand sides therefore
does not keep `|ln| = 1` unless the integration renormalizes `ln`.

Claim evidence:
- Claim: `SEOBNRv5_spin_evolution_equations` implements Eqs. (65), (66), (69), and (71a)-(71d) of arXiv:2303.18143v2 for black holes, omitting the terms proportional to `Ctilde_1ES2`, and stores `L/(M mu)` and `omega_dot = 3 v^2 v_dot`; its `ln_dot` has the component parallel to a unit `ln` given by Eq. (72), so these right-hand sides do not preserve `|ln| = 1`.
- Role: descriptive behavior
- Deciding authority: [SEOBNRv5_spin_evolution_equations.py](../../../nrpy/equations/seobnr/SEOBNRv5_spin_evolution_equations.py), `SEOBNRv5_spin_evolution_equations`
- Corroboration: [Khalil et al., arXiv:2303.18143v2](https://arxiv.org/pdf/2303.18143v2), Eqs. (65), (66), (69), (71), and (72), and the paragraph before Eq. (72)

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
- [SEOBNRv5 dynamics paper landing page](https://arxiv.org/abs/2303.18143) - background orientation only; the version-pinned PDF below carries the equation mapping
- [Khalil et al., arXiv:2303.18143v2](https://arxiv.org/pdf/2303.18143v2) - Eqs. (65), (66), (69), (71), and (72), spin-evolution right-hand sides

## See Also

- [SEOBNR And BOB](index.md)
- [Aligned-Spin Hamiltonian](aligned-spin-hamiltonian.md)
- [Precessing Rotations And Ringdown](precessing-rotations-and-ringdown.md)
- [BOB Waveforms](bob-waveforms.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
