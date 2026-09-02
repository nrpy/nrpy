# SEOBNR Precessing Rotations And Ringdown

> Map co-precessing frame rotations, inertial polarizations, and precessing merger-ringdown frame quantities. · Status: confirmed · Last reconciled: 07-13-2026
> Up: [SEOBNR And BOB](index.md)

## Summary

The rotations and ringdown SEOBNR modules turn precessing-frame ingredients into
observer-frame and post-merger symbolic quantities. `SEOBNRv5_Coprecessing_Rotations`
owns the final-angular-momentum frame triad, inertial Euler-angle extraction,
Wigner small-d rotations, and inertial plus/cross polarizations.
`SEOBNRv5_MergerRingdown_Quantities` owns post-merger Euler angles, precession
frequency selection, and QNM-frequency conversion from the final-angular-momentum
frame to the co-precessing frame.

## Detail

`SEOBNRv5_Coprecessing_Rotations` defaults to the mode set `(2,2)`, `(2,1)`,
`(3,3)`, `(3,2)`, `(4,4)`, `(4,3)`, and `(5,5)`. Callers may pass another
nonempty list of `(l, m)` modes. The constructor directly rejects negative `m`
values because negative-mode symmetries are handled internally; subsequent
`wigner_d_small` calls reject `l < 0` and `|m| > l`. No duplicate-mode check is
performed.

The class first constructs the `J_f`-frame triad from symbols
`J_f_x`, `J_f_y`, and `J_f_z`. To keep generated code well behaved near
coordinate degeneracies, the construction normalizes with small lower bounds
and blends Gram-Schmidt seed directions before storing the matrix entries as
`R_Eq15_00` through `R_Eq15_22`.

The rotation module then combines observer-angle rotations, the transposed
`J_f` frame, and co-precessing Euler angles `alpha_JP`, `beta_JP`, and
`gamma_JP` to build an inertial transformation matrix. It stores `M33`,
clamps it before computing `beta_PI`, and exposes explicit Euler-angle
branches: `alpha_PI_generic`, `gamma_PI_generic`, `alpha_PI_pole_pos`,
`gamma_PI_pole_pos`, `alpha_PI_pole_neg`, and `gamma_PI_pole_neg`. The aliases
`alpha_PI` and `gamma_PI` point at the generic branch.

Wigner rotations are implemented by the cached
`wigner_d_small_template` and checked wrapper `wigner_d_small`. The projection
helper `_polarizations_from_angles` combines each co-precessing mode's real and
imaginary symbols with Wigner small-d factors and Euler phases. The constructor
stores branch-specific polarizations `h_plus_I_generic`, `h_cross_I_generic`,
`h_plus_I_pole_pos`, `h_cross_I_pole_pos`, `h_plus_I_pole_neg`, and
`h_cross_I_pole_neg`, with `h_plus_I` and `h_cross_I` as generic-branch
aliases.

`SEOBNRv5_MergerRingdown_Quantities` starts from final spin and final orbital
angular-momentum symbols, stores the matching angles `alpha_match`,
`beta_match`, and `gamma_match`, and stores times `t` and `t_match`.
It computes `chi_f_dot_L_f` and uses no-branch coordinate-bound helpers to
select the prograde or retrograde precession-frequency expression. The selected
frequency is stored as `omega_prec`.

At `chi_f_dot_L_f == 0`, both strict branch masks vanish, so the current
expression selects `omega_prec = 0`.

The merger-ringdown class advances post-merger Euler angles as
`alpha_merger_RD`, `beta_merger_RD`, and `gamma_merger_RD`. It declares
J-frame QNM symbols for `(2,2)`, `(2,1)`, `(3,3)`, `(3,2)`, `(4,4)`, `(4,3)`,
and `(5,5)`, then stores co-precessing-frame QNM outputs with names such as
`omega_QNM_P_l2_m2`, `omega_QNM_P_l3_m3`, and `omega_QNM_P_l5_m5`.

Both modules validate by running doctests, converting an instantiated object's
state through `process_dictionary_of_expressions`, and comparing against the
module's trusted dictionary. Representative trusted keys cover rotation matrix
entries, Euler-angle branches, inertial polarizations, post-merger Euler angles,
`omega_prec`, and J-frame/P-frame QNM frequencies.

The trusted rotation object uses only the default mode list, and the trusted
merger-ringdown object uses its fixed implemented mode list. Alternate caller
mode lists and Euler-pole branch selection in generated code have no separate
runtime evidence here. Sampled symbolic matching does not prove orthonormality
for all inputs, generated branch correctness, or waveform accuracy.

## Sources

- [SEOBNRv5_coprecessing_rotations_quantities.py](../../../nrpy/equations/seobnr/SEOBNRv5_coprecessing_rotations_quantities.py) - `SEOBNRv5_Coprecessing_Rotations`, `wigner_d_small_template`, `wigner_d_small`, `_polarizations_from_angles`
- [SEOBNRv5_coprecessing_rotations_quantities.py](../../../nrpy/equations/seobnr/SEOBNRv5_coprecessing_rotations_quantities.py) - `R_Eq15_00`, `R_Eq15_22`, `M33`, `beta_PI`, `alpha_PI_generic`, `h_plus_I_generic`, `h_cross_I_generic`
- [SEOBNRv5_merger_ringdown.py](../../../nrpy/equations/seobnr/SEOBNRv5_merger_ringdown.py) - `SEOBNRv5_MergerRingdown_Quantities`, `chi_f_dot_L_f`, `omega_prec`, `alpha_merger_RD`
- [SEOBNRv5_coprecessing_rotations_quantities.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_coprecessing_rotations_quantities.py) - `trusted_dict`
- [SEOBNRv5_merger_ringdown.py](../../../nrpy/equations/seobnr/tests/SEOBNRv5_merger_ringdown.py) - `trusted_dict`
- [Min_Max_and_Piecewise_Expressions.py](../../../nrpy/equations/grhd/Min_Max_and_Piecewise_Expressions.py) - `coord_greater_bound`, `coord_less_bound`
- [SEOBNRv5PHM current latest paper page](https://arxiv.org/abs/2303.18046) - background orientation only; cited rotation/ringdown equations are not yet audited to a pinned revision

## See Also

- [SEOBNR And BOB](index.md)
- [Aligned-Spin Waveforms](aligned-spin-waveforms.md)
- [Precessing Dynamics](precessing-dynamics.md)
- [BOB Waveforms](bob-waveforms.md)
- [Geometry And Special-Function Support](../geometry-and-special-function-support.md)
- [Trusted Expression Pipeline](../trusted-expression-pipeline.md)
