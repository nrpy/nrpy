# TOV Equations

> Map the symbolic Tolman-Oppenheimer-Volkoff ODE RHS module and validation hook. · Status: confirmed
> Up: [Equations](index.md)

## Summary

`TOV_Equations` is the equation-side source for the four ordinary differential
equations used by relativistic-star initial data. It stores the
Schwarzschild-radius pressure, metric-potential, enclosed-mass, and isotropic
radius RHS expressions as SymPy object attributes.

## Detail

The class declares the integration-facing symbols `r_Schw`,
`y[TOVOLA_R_ISO]`, `rho_energy`, `y[TOVOLA_PRESSURE]`,
`y[TOVOLA_MASS]`, and `M_PI`. Those names intentionally match the
integration-facing state-vector slots used by relativistic-star
initial-data solvers.

The stored RHS attributes are `dP_dr`, `dnu_dr`, `dM_dr`, and `dr_iso_dr`.
The metric-potential derivative depends on enclosed mass, radius, and pressure:
`dnu_dr = 2*(M + 4*M_PI*r_Schw**3*P)/(r_Schw*(r_Schw - 2*M))`.
Energy density enters the pressure derivative through
`dP_dr = -(rho_energy + P)*dnu_dr/2`; it is absent from `dnu_dr` itself.
`dM_dr` evolves the enclosed mass from `rho_energy`, and `dr_iso_dr` relates
the areal coordinate `r_Schw` to isotropic radius.

Evidence (descriptive behavior): `tov/TOV_equations.py::TOV_Equations.__init__`.

`TOV_Equations` defines only the symbolic stellar-equilibrium ODE RHSs and their
trusted-value checks; equation-of-state closures, radial integration steppers,
outer-boundary detection, interpolation, and `CFunction` registration modules live in
`nrpy/infrastructures/BHaH/general_relativity/TOVola/` (see [GR Initial Data, Diagnostics, And Psi4](../infrastructures/bhah/gr-initial-data-and-psi4.md#twopunctures-and-tovola-solvers)).

The expressions divide by `r_Schw`, by `1 - 2*M/r_Schw`, and, for the isotropic
radius equation, by its square root. The class performs no origin expansion,
horizon/domain check, or physical-state validation. A downstream integrator
must supply its regular-center treatment and remain in the intended stellar
domain.

Validation is compact: the module's `__main__` path runs doctests, processes
the class dictionary with `process_dictionary_of_expressions`, and compares it
against `nrpy/equations/tov/tests/TOV_equations.py`. The trusted dictionary
covers the declared integration symbols and all four RHS attributes.

That comparison is one deterministic sampled-symbol fingerprint. It does not
integrate a star, enforce an equation of state or zero-pressure boundary condition, build or
run generated code, or establish mass-radius or numerical-accuracy results.

## Sources

- [TOV_equations.py](../../nrpy/equations/tov/TOV_equations.py) - `TOV_Equations`, `dP_dr`, `dnu_dr`, `dM_dr`, `dr_iso_dr`
- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `compare_or_generate_trusted_results`
- [TOV_equations.py](../../nrpy/equations/tov/tests/TOV_equations.py) - `trusted_dict`
- [Oppenheimer and Volkoff, On Massive Neutron Cores](https://link.aps.org/doi/10.1103/PhysRev.55.374) - background for the relativistic stellar-equilibrium system

## See Also

- See also: [Conformally Flat Elliptic](conformally-flat-elliptic.md)
- See also: [GRHD](grhd.md)
- See also: [Expression Validation Helpers](../validation/expression-validation-helpers.md)
