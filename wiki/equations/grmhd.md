# GRMHD

> Map NRPy's magnetic extension of the GRHD symbolic equations and face-flux helpers. · Status: confirmed
> Up: [Equations](index.md)

## Summary

`GRMHD_Equations` inherits the GRHD equation builder. It adds magnetic
stress-energy to `compute_T4UU`. The inherited energy and momentum equations
then use the total tensor. Density, electron fraction, and entropy retain
their GRHD fluid-advection forms. Companion modules compute GRMHD signal-speed
bounds from a fast magnetosonic speed estimate and HLL fluid fluxes. The GRMHD
modules provide no evolved magnetic field, induction equation or flux,
magnetic-constraint treatment, or primitive recovery.

## Detail

Notation: `alpha` is the lapse, `beta^i` the shift with
`beta_i = gamma_ij beta^j`, `gamma_ij` the spatial metric, and `g_{mu nu}` the
spacetime metric with signature (-,+,+,+) in geometrized units (`G = c = 1`).
`u^mu` is the fluid four-velocity, `rho_b` the baryon density, `h` the specific
enthalpy, and `c_s^2` the squared sound speed. Greek indices run 0-3 and Latin
indices 1-3. Python lists such as `BmagU[i]` and `ReU[i]` use zero-based
spatial indices, so `BmagU[i]` holds `B^{i+1}` and pairs with `u4U[i+1]`.
In this page, `grhd.` and `grmhd.` abbreviate `nrpy.equations.grhd.` and
`nrpy.equations.grmhd.`.

`BmagU^i = B_Gaussian^i / sqrt(4*pi)` denotes the magnetic field measured by
Eulerian observers, who are normal to the spatial slices, in Heaviside-Lorentz
units. Callers of `GRMHD_Equations`,
`grmhd.characteristic_speeds.find_cmax_cmin`, and
`grmhd.HLL_fluxes.calculate_HLL_fluxes` supply magnetic components already
divided by `sqrt(4*pi)`; none of these functions applies that factor.

`GRMHD_Equations.BmagU` holds components in the coordinate basis of
`CoordSystem`. The class declares `rescaledBmagU0`, `rescaledBmagU1`, and
`rescaledBmagU2` and constructs `BmagU[i] = rescaledBmagU[i] * ReU[i]` during
initialization. Thus values substituted for those symbols must be
`BmagU^i / ReU[i]` in the evolved coordinate basis. To assign a field after
construction, set `BmagU[i]` before calling `compute_T4UU` or
`construct_all_equations`. Later changes to `rescaledBmagU[i]` do not change
`BmagU`. The scale factors `ReU` follow
[Jacques et al., Eq. (22)](https://arxiv.org/pdf/2412.03659v2).

The face functions `grmhd.characteristic_speeds.find_cmax_cmin` and
`grmhd.HLL_fluxes.calculate_HLL_fluxes` do not apply `ReU`. They contract the
face fields (`BmagU_r` and `BmagU_l`, or `BmagrU` and `BmaglU`) with the
supplied face metric, shift, and four-velocities, so these fields must use the
same basis as `gamma_faceDD`, `beta_faceU`, and the face four-velocities. For
reference-metric evolutions this is the rescaled basis that the GRHD face
functions use: supply `B^i / (sqrt(4*pi) * ReU[i])`, and the returned speeds and
fluxes are rescaled quantities. The [Fishbone-Moncrief
conversion](general-relativity/fishbone-moncrief.md) also requires a basis
transform when the evolved coordinates differ from its spherical coordinates.

Claim evidence:
- Claim: GRMHD functions take Eulerian magnetic components already divided by `sqrt(4*pi)`. `GRMHD_Equations` multiplies its `rescaledBmagU` symbols by `ReU` during initialization, so `GRMHD_Equations.BmagU` uses the `CoordSystem` coordinate basis. The face functions apply no `ReU` and require their magnetic fields in the basis of the supplied face metric, shift, and four-velocities. Fishbone-Moncrief components need conversion to the evolved basis and scaling before use.
- Role: descriptive behavior
- Deciding authority: [GRMHD_equations.py](../../nrpy/equations/grmhd/GRMHD_equations.py), `GRMHD_Equations.__init__` and `compute_smallb4U`; [characteristic_speeds.py](../../nrpy/equations/grmhd/characteristic_speeds.py), `find_cmax_cmin`; [HLL_fluxes.py](../../nrpy/equations/grmhd/HLL_fluxes.py), `calculate_HLL_fluxes`; [fishbone_moncrief.py](../../nrpy/equations/general_relativity/fishbone_moncrief/fishbone_moncrief.py), `FishboneMoncriefID._compute_initial_data` magnetic-field construction
- Corroboration: [GRMHD_equations_Spherical.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_Spherical.py), `trusted_dict` for the nonunit `ReU` mapping; [GRHD HLL_fluxes.py](../../nrpy/equations/grhd/HLL_fluxes.py), `calculate_Tmunu_and_contractions_from_equations` rescaled face inputs

`compute_T4UU` calls the module-level functions `compute_smallb4U` and
`compute_smallb2` and stores `smallb4U` and `smallb2`; these attributes do not
exist before that call. `compute_smallb4U` forms `u_i BmagU^i`, where
`u_i = g_{i mu} u^mu = beta_i u^0 + gamma_{ij} u^j`, and requires `BmagU`,
`gammaDD`, `betaU`, and `u4U` in one basis:

```text
b^0 = u_i BmagU^i / alpha
b^i = (BmagU^i + (u_j BmagU^j) u^i) / (alpha u^0)
b^2 = g_{mu nu} b^mu b^nu
```

The comoving-field components and Gaussian scaling follow [Duez et al., Eqs.
(23)-(24) and (31)](https://arxiv.org/pdf/astro-ph/0503420v2); their stress
tensor uses `b^2` in Eq. (32).

Claim evidence:
- Claim: `compute_smallb4U` and `compute_smallb2` compute the displayed comoving-field expressions.
- Role: public/scientific contract
- Deciding authority: [Duez et al. (2005)](https://arxiv.org/pdf/astro-ph/0503420v2), Eqs. (23)-(24) and (31)-(32)
- Corroboration: [GRMHD_equations.py](../../nrpy/equations/grmhd/GRMHD_equations.py), `compute_smallb4U` and `compute_smallb2`; [GRMHD_equations_Cartesian.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_Cartesian.py), `trusted_dict`

`compute_T4UU` starts with the inherited perfect-fluid tensor, then adds `b^2
u^mu u^nu + (b^2/2) g^{mu nu} - b^mu b^nu`. It updates the inherited
reference-metric-rescaled tensor with the same contribution. `compute_T4UD`
lowers the fluid part with the GRHD method and adds the magnetic part `b^2 u^mu
u_nu + (b^2/2) delta^mu_nu - b^mu b_nu` directly. `compute_tau_tilde_fluxU`
likewise applies the GRHD method to the fluid part and adds `alpha^2 e6phi
T_EM^{0j}`. These two methods avoid combining the magnetic terms over a common
denominator with `sp.together`, which the GRHD methods apply to their fluid
expressions. The other inherited methods form the energy and momentum Valencia
conserved quantities, fluxes, source terms, and reference-metric connection
terms from the total tensor. Density, electron fraction, and entropy branches
retain their GRHD advection expressions. No function in `nrpy/equations/grmhd/`
constructs an evolved magnetic conserved variable, induction flux, vector
potential, magnetic-constraint treatment, or primitive recovery. The magnetic
addition is [Duez et al., Eqs.
(32)-(33)](https://arxiv.org/pdf/astro-ph/0503420v2); its tensor-component
rescaling follows [Jacques et al., Eq.
(22)](https://arxiv.org/pdf/2412.03659v2).

Claim evidence:
- Claim: The GRMHD class adds the displayed magnetic stress-energy to the GRHD tensor, rescales it with `ReU`, lowers its index and forms the energy flux with GRHD methods for the fluid part plus directly added magnetic parts, and inherits the other energy and momentum equations using the total tensor; density, electron fraction, and entropy retain fluid advection. The GRMHD modules do not compute induction terms or recover primitive variables.
- Role: descriptive behavior
- Deciding authority: [GRMHD_equations.py](../../nrpy/equations/grmhd/GRMHD_equations.py), `GRMHD_Equations.compute_T4UU`, `compute_T4UD`, and `compute_tau_tilde_fluxU`; [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py), `construct_all_equations` and downstream methods
- Corroboration: [GRMHD_equations_Cartesian.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_Cartesian.py) and [GRMHD_equations_Spherical.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_Spherical.py), `trusted_dict`

`grhd.characteristic_speeds.find_cmax_cmin` and
`grmhd.characteristic_speeds.find_cmax_cmin` both take `flux_dirn`,
`gamma_faceDD`, `beta_faceU`, `alpha_face`, `u4U_r`, and `u4U_l` first. GRHD
then takes `v02_r` and `v02_l`; GRMHD instead takes `BmagU_r`, `BmagU_l`,
`rho_b_r`, `rho_b_l`, `h_r`, `h_l`, `cs2_r`, and `cs2_l`. The last two GRMHD
inputs are sound-speed squares. For each reconstructed face state,
`grmhd.characteristic_speeds.compute_v02` computes the squared fluid-frame
Alfvén speed `v_A^2 = b^2/(rho_b h + b^2)` and returns the fast magnetosonic
speed estimate `v_0^2 = v_A^2 + c_s^2(1-v_A^2)`. GRMHD `find_cmax_cmin` passes
both `v_0^2` values to GRHD `find_cmax_cmin` and returns nonnegative `(cmin,
cmax)` HLL bounds. The Alfvén and approximate signal speeds come from [Duez et
al., Eq. (50) and the text following
it](https://arxiv.org/pdf/astro-ph/0503420v2).
`grmhd.HLL_fluxes.calculate_HLL_fluxes` takes the GRHD HLL arguments with
`BmagrU` and `BmaglU` inserted after `u4lU`. Its argument order `(flux_dirn,
alpha_face, gamma_faceDD, beta_faceU, ...)` differs from the `find_cmax_cmin`
order `(flux_dirn, gamma_faceDD, beta_faceU, alpha_face, ...)`, as in GRHD. It
follows the GRHD function step by step: for each state, it creates a Cartesian
`GRMHD_Equations` object, assigns the face field to `BmagU`, and calls the GRHD
`calculate_Tmunu_and_contractions_from_equations` function to compute fluid
conserved variables and physical fluxes including magnetic stress-energy.
`HLL_solver` combines the right and left results. The HLL flux is [Duez et al.,
Eq. (48)](https://arxiv.org/pdf/astro-ph/0503420v2).

Claim evidence:
- Claim: `grmhd.characteristic_speeds.compute_v02` computes the displayed `v_A^2` and `v_0^2`; GRMHD `find_cmax_cmin` passes the two face-state values to GRHD `find_cmax_cmin` and returns `(cmin, cmax)`; GRMHD HLL fluid fluxes use the GRHD contraction function and `HLL_solver`.
- Role: descriptive behavior
- Deciding authority: [characteristic_speeds.py](../../nrpy/equations/grmhd/characteristic_speeds.py), `compute_v02` and `find_cmax_cmin`; [HLL_fluxes.py](../../nrpy/equations/grmhd/HLL_fluxes.py), `calculate_HLL_fluxes`
- Corroboration: [GRHD characteristic_speeds.py](../../nrpy/equations/grhd/characteristic_speeds.py), `find_cmax_cmin`; [GRMHD characteristic_speeds.py](../../nrpy/equations/grmhd/tests/characteristic_speeds.py) and [GRMHD HLL_fluxes.py](../../nrpy/equations/grmhd/tests/HLL_fluxes.py), `trusted_dict`

Each GRMHD validation path follows its GRHD counterpart. `GRMHD_equations.py`
first sets `BmagU = 0` in Cartesian coordinates and requires every GRHD
expression attribute, including source and connection terms, to equal the GRHD
result exactly. It then loops over the seven coordinate systems that
`GRHD_equations.py` validates: Spherical, SinhSpherical, SinhSpherical with
`enable_rfm_precompute`, Cartesian, SinhCartesian, SinhCylindrical, and
SinhSymTP. For each, it runs `construct_all_equations` with the symbolic
`rescaledBmagU` field and compares every stored expression, including magnetic
source terms, connection terms, and curvilinear fluxes, with
`GRMHD_equations_<CoordSystem>.py`.

The speed and HLL modules use the symbolic BSSN-derived Cartesian face states of
the GRHD tests, with symbolic `BmagrU` and `BmaglU` added. With zero field, the
speed module requires the returned bounds to equal the GRHD bounds exactly for
every `flux_dirn` from `0` to `2`, and the HLL module requires the returned
fluxes to equal the GRHD fluxes exactly for `flux_dirn = 0` on a flat face
metric. As in GRHD, the speed trusted dictionary records the left-state
`find_cp_cm` roots for `flux_dirn = 1` and `(cmin, cmax)` for `flux_dirn = 2`.
At the sample values, `cmax` there is zero and the left state sets `cmin`, so
the speed module also requires GRMHD `find_cmax_cmin` to equal GRHD
`find_cmax_cmin` built from the two face-state `compute_v02` values on a flat
face metric in every `flux_dirn`; this ties each face field, density, and
enthalpy to its own state. Also as in GRHD, the HLL trusted dictionary records
the right-state conserved variables and physical fluxes for `flux_dirn = 2` and
one `HLL_solver` combination with the `flux_dirn = 1` speed bounds. It records
all HLL fluxes for `flux_dirn = 1`, where both speed bounds are nonzero at the
sample values, so both face states and the magnetic speed bounds enter the
stored fluxes. Because the left state sets both bounds at the sample values, it
also records the `flux_dirn = 1` fluxes with the two face states exchanged,
which tests the right-state inputs to the speed bounds. `find_cmax_cmin` and
`calculate_HLL_fluxes` use the same expressions for every `flux_dirn`; only the
selected component index changes.

Before computing trusted speed bounds and HLL fluxes, the GRMHD validation paths
rebuild each expression once with `sp.Abs` in place of `nrpyAbs`, sharing
repeated subexpressions, because `.subs` on these large magnetic expressions is
slow. Trusted values evaluate the expressions at fixed sample values of their
free symbols; these samples need not form a normalized four-velocity or a
physically admissible state. These comparisons check symbolic construction and
sampled regression, not generated C, numerical evolution, or physical
admissibility of arbitrary caller-supplied states. As in GRHD, callers must
supply valid lapse, metric, four-velocity, densities, and nonzero denominators.

Claim evidence:
- Claim: The GRMHD equation class checks exact Cartesian B=0 equality with GRHD and, for the seven GRHD coordinate systems, compares all stored magnetic expressions, including source and connection terms and curvilinear fluxes, with trusted values. The speed module checks exact B=0 equality with GRHD for symbolic face states in all three flux directions, and the HLL module in flux direction 0 on a flat face metric; they compare left-state roots, `(cmin, cmax)`, right-state conserved variables and fluxes, and one direction of HLL fluxes, in which both speed bounds are nonzero, with and without exchanged face states, with trusted values. The speed module requires exact agreement with GRHD bounds built from the two face-state `compute_v02` values on a flat face metric in all three flux directions. Trusted samples need not be physical states.
- Role: descriptive behavior
- Deciding authority: [GRMHD_equations.py](../../nrpy/equations/grmhd/GRMHD_equations.py), [characteristic_speeds.py](../../nrpy/equations/grmhd/characteristic_speeds.py), and [HLL_fluxes.py](../../nrpy/equations/grmhd/HLL_fluxes.py), `__main__` validation paths
- Corroboration: `GRMHD_equations_<CoordSystem>.py` files such as [GRMHD_equations_Cartesian.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_Cartesian.py) and [GRMHD_equations_SinhSpherical_rfm_precompute.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_SinhSpherical_rfm_precompute.py), [GRMHD characteristic_speeds.py](../../nrpy/equations/grmhd/tests/characteristic_speeds.py), and [GRMHD HLL_fluxes.py](../../nrpy/equations/grmhd/tests/HLL_fluxes.py), `trusted_dict`

## Sources

- [Duez et al. (2005)](https://arxiv.org/pdf/astro-ph/0503420v2) - Eqs. (23)-(24), (31)-(33), (48), (50)
- [Jacques et al.](https://arxiv.org/pdf/2412.03659v2) - Eq. (22)
- [GRMHD_equations.py](../../nrpy/equations/grmhd/GRMHD_equations.py) - `compute_smallb4U`, `compute_smallb2`, `GRMHD_Equations.compute_T4UU`
- [characteristic_speeds.py](../../nrpy/equations/grmhd/characteristic_speeds.py) - `compute_v02`, `find_cmax_cmin`
- [HLL_fluxes.py](../../nrpy/equations/grmhd/HLL_fluxes.py) - `calculate_HLL_fluxes`
- [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py) - `GRHD_Equations.construct_all_equations`
- [GRHD characteristic_speeds.py](../../nrpy/equations/grhd/characteristic_speeds.py) - `find_cmax_cmin`
- [GRHD HLL_fluxes.py](../../nrpy/equations/grhd/HLL_fluxes.py) - `calculate_Tmunu_and_contractions_from_equations`, `HLL_solver`
- [fishbone_moncrief.py](../../nrpy/equations/general_relativity/fishbone_moncrief/fishbone_moncrief.py) - `FishboneMoncriefID._compute_initial_data`
- [GRMHD_equations_Spherical.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_Spherical.py) - `trusted_dict`
- [GRMHD_equations_SinhSpherical.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_SinhSpherical.py) - `trusted_dict`
- [GRMHD_equations_SinhSpherical_rfm_precompute.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_SinhSpherical_rfm_precompute.py) - `trusted_dict`
- [GRMHD_equations_Cartesian.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_Cartesian.py) - `trusted_dict`
- [GRMHD_equations_SinhCartesian.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_SinhCartesian.py) - `trusted_dict`
- [GRMHD_equations_SinhCylindrical.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_SinhCylindrical.py) - `trusted_dict`
- [GRMHD_equations_SinhSymTP.py](../../nrpy/equations/grmhd/tests/GRMHD_equations_SinhSymTP.py) - `trusted_dict`
- [GRMHD characteristic_speeds.py](../../nrpy/equations/grmhd/tests/characteristic_speeds.py) - `trusted_dict`
- [GRMHD HLL_fluxes.py](../../nrpy/equations/grmhd/tests/HLL_fluxes.py) - `trusted_dict`

## See Also

- Parent: [Equations](index.md)
- Depends on: [GRHD](grhd.md)
- Validated by: [Trusted Expression Pipeline](trusted-expression-pipeline.md)
- Depends on: [Metric Conversions And Matter](general-relativity/metric-conversions-and-matter.md)
- Depends on: [Fishbone-Moncrief](general-relativity/fishbone-moncrief.md)
- Depends on: [Reference Metrics](../core/reference-metrics.md)
