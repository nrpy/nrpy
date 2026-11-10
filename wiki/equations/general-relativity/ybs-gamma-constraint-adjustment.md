# YBS Gamma-Constraint Adjustment

> Define NRPy's opt-in Yo-Baumgarte-Shapiro Gamma-constraint addition to the BSSN conformal-connection RHS, its coefficient mapping, timestep effect, parameter registration, and gauge forwarding. · Status: confirmed
> Up: [General Relativity](index.md)

## Summary

`BSSNRHSs(..., enable_YBS_Gamma_constraint_adjustment=True)` adds the
lower-order term `-YBS_chi*C^i*Dbar_j beta^j` to the unrescaled BSSN
conformal-connection RHS, where `C^i` is the Lambda constraint. The term is an
increment to NRPy's existing equation, so `YBS_chi = 2*xi/3` in the
Yo-Lin-Cao parameterization, and the runtime default `YBS_chi = 0.0` removes it
exactly. It leaves the principal characteristic speeds unchanged, and its
isolated local factor attenuates `C^i` only where `YBS_chi*Dbar_j beta^j > 0`.
It is separate from the YBS-MOM momentum-gradient adjustment.

## Detail

`BSSNRHSs(..., enable_YBS_Gamma_constraint_adjustment=True)` adds

```text
-YBS_chi * (LambdabarU[i] - DGammaU[i]) * Dbarbetacontraction
```

to the unrescaled conformal-connection RHS, where
`Dbarbetacontraction = Dbar_j beta^j` and the parenthesized factor is the
Lambda constraint `C^i`. This is the incremental term for NRPy's existing
equation, not the full paper coefficient: the existing RHS already contains
`+(2/3) DGammaU[i] * Dbarbetacontraction`, which supplies the paper's baseline
`-(2/3) C^i Dbar_j beta^j` contribution when the connection equation is
rewritten in terms of `C^i`. The full Yo-Baumgarte-Shapiro coefficient is
`-(YBS_chi + 2/3) C^i Dbar_j beta^j`; adding that full expression again would
double-count the baseline. Yo-Lin-Cao write the same coefficient as
`-(2/3)(xi + 1)`, so NRPy's incremental parameter satisfies
`YBS_chi = 2*xi/3`. Their usual `xi=1` choice therefore maps to
`YBS_chi = 2/3`. The runtime default is `0.0`, which removes the term exactly;
a user who enables the option sets `YBS_chi` above zero in the parameter file.
The recommended maximum is `2/3` for BSSN. It is the `xi=1` value that Yo, Lin,
and Cao state is usually chosen in Eq. (47), a customary choice and not a
stability bound. The term is lower differential order and changes the local
growth rate of `C^i` by `-YBS_chi*Dbar_j beta^j`. It leaves the principal
characteristic speeds unchanged, but sufficiently large decay rates can
restrict explicit timesteps. For example, classical RK4 applied to
`dC/dt=-kappa*C` has amplification `R(-3)=1.375` when `kappa*Delta t=3`,
despite continuum decay. Larger parameter values have no supporting bound
here, and negative values amplify `C^i` wherever `Dbar_j beta^j > 0`.

Evidence (descriptive behavior): [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`; [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2), Eq. (47); corroboration: [Dendro rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`; [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS Gamma branch

Evidence (descriptive behavior): [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs.__init__` YBS Gamma addition; [rk_butcher_table_dictionary.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/rk_butcher_table_dictionary.py), `generate_Butcher_tables`, classical RK4 coefficients; corroboration: the classical RK4 polynomial `R(z)=1+z+z^2/2+z^3/6+z^4/24` gives the displayed amplification directly

The equation layer uses a plain real SymPy symbol for `YBS_chi`; constructing
enabled symbolic RHSs does not register a global `CodeParameter`. The BHaH
function `register_CFunction_rhs_eval` registers `YBS_chi` only when the option is
enabled, for either formulation, as a runtime `REAL`, using the
default `0.0` (term disabled at run time), with `commondata=True` and
`add_to_parfile=True`; there is no Python coefficient argument.
`BSSNRHSs_dict.get_rhs` receives the Boolean
separately from the unchanged coordinate/options string and keeps enabled
expressions in a separate internal cache. With the option false, neither the
term nor BHaH parameter registration is present, and the existing base cache
is used.

`BSSN_gauge_RHSs` forwards the Boolean to the shared BSSN RHS cache. The
derivative-based shift choices `GammaDriving2ndOrder_NoCovariant`,
`GammaDriving2ndOrder_Covariant`,
`GammaDriving2ndOrder_Covariant__Hatted`, and `NonAdvectingGammaDriving`
consume the adjusted Lambda RHS by default. Supplying
`evolved_connection_rhsU` replaces that input with the caller's unrescaled
connection RHS; forwarding the Boolean alone does not adjust a replacement
vector. `Frozen` and both first-order covariant variants do not consume that
RHS and remain unchanged.

Evidence (descriptive behavior): [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py), `BSSN_gauge_RHSs`; corroboration: [fCCZ4_gauge_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_gauge_RHSs.py), `fCCZ4_gauge_RHSs`, supplies the fCCZ4 connection RHS

With `Theta = Dbar_j beta^j`, the corresponding constraint propagation has
the lower-order addition `-YBS_chi*Theta*C^i`:

```text
partial_t C^i = beta^j Dhat_j C^i
                + 2 alpha exp(4 phi) M^i - YBS_chi Theta C^i.
```

This term attenuates its isolated local factor only where
`YBS_chi*Theta > 0`; it amplifies where that product is negative and is neutral
where it vanishes. Advection and the momentum-constraint source can dominate,
so a globally positive `YBS_chi` is not a uniform damping guarantee when
`Theta` changes sign. Because the addition is lower order, it does not change
the current second-order system's principal symbol or high-frequency
characteristic speeds; it can still change finite-frequency growth.

Evidence (public/scientific contract): [Yo, Baumgarte, and Shapiro, arXiv:gr-qc/0209066v2](https://arxiv.org/pdf/gr-qc/0209066v2), Eq. (45); [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2), Eq. (47); corroboration: [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py), `BSSNRHSs`; [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py), `BSSN_gauge_RHSs`; [BSSN_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_RHSs_Cartesian.py), representative jointly enabled `trusted_dict`

## Sources

- [BSSN_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs.__init__` YBS Gamma branch, `BSSNRHSs_dict.get_rhs`
- [BSSN_gauge_RHSs.py](../../../nrpy/equations/general_relativity/BSSN_gauge_RHSs.py) - `BSSN_gauge_RHSs`
- [fCCZ4_gauge_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_gauge_RHSs.py) - `fCCZ4_gauge_RHSs`
- [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py) - BHaH `register_CFunction_rhs_eval`, `YBS_chi` registration
- [Dendro rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`
- [rk_butcher_table_dictionary.py](../../../nrpy/infrastructures/BHaH/MoLtimestepping/rk_butcher_table_dictionary.py) - `generate_Butcher_tables`
- [BSSN_RHSs_Cartesian.py](../../../nrpy/equations/general_relativity/tests/BSSN_RHSs_Cartesian.py) - `trusted_dict`
- [Yo, Baumgarte, and Shapiro, arXiv:gr-qc/0209066v2](https://arxiv.org/pdf/gr-qc/0209066v2) - Eq. (45), adjusted Gamma-constraint coefficient
- [Yo, Lin, and Cao, arXiv:1205.5111v2](https://arxiv.org/pdf/1205.5111v2) - Eq. (47), `xi` parameterization of the adjusted coefficient

## See Also

- Parent: [General Relativity](index.md)
- Depends on: [BSSN Family](bssn-family.md)
- Contrasts with: [YBS-MOM Timestep-Scaled Momentum Adjustment](ybs-momentum-damping.md)
- See also: [Fully Covariant Conformal Z4](fccz4.md)
