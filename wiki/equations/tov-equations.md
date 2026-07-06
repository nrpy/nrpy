# TOV Equations

> Map the symbolic Tolman-Oppenheimer-Volkoff ODE RHS module and validation hook. · Status: confirmed · Last reconciled: 06-29-2026
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
integration-facing state-vector slots used by downstream relativistic-star
initial-data consumers.

The stored RHS attributes are `dP_dr`, `dnu_dr`, `dM_dr`, and `dr_iso_dr`.
`dP_dr` and `dnu_dr` use the same relativistic pressure-gravity factor built
from enclosed mass, radius, pressure, and total energy density; `dM_dr` evolves
the enclosed mass from `rho_energy`; and `dr_iso_dr` converts the
Schwarzschild-radius integration to isotropic radius.

The equation module does not choose an equation of state, integration stepper,
stellar surface rule, interpolation strategy, or generated C layout. Those are
infrastructure responsibilities; this page owns only the symbolic ODE RHSs and
their validation evidence.

Validation is compact: the module's `__main__` path runs doctests, processes
the class dictionary through the trusted-expression pipeline, and compares it
against `nrpy/equations/tov/tests/TOV_equations.py`. The trusted dictionary
covers the declared integration symbols and all four RHS attributes.

## Sources

- [TOV_equations.py](../../nrpy/equations/tov/TOV_equations.py) - `TOV_Equations`, `dP_dr`, `dnu_dr`, `dM_dr`, `dr_iso_dr`
- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `compare_or_generate_trusted_results`
- [TOV_equations.py](../../nrpy/equations/tov/tests/TOV_equations.py) - `trusted_dict`
- [Oppenheimer and Volkoff, On Massive Neutron Cores](https://link.aps.org/doi/10.1103/PhysRev.55.374) - background for the relativistic stellar-equilibrium system

## See Also

- [Equations](index.md)
- [Conformally Flat Elliptic](conformally-flat-elliptic.md)
- [GRHD](grhd.md)
- [Trusted Expression Pipeline](trusted-expression-pipeline.md)
- [Lifecycle And Project Assembly](../infrastructures/bhah/lifecycle-and-project-assembly.md)
