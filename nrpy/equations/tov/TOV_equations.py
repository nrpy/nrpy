"""
Construct the Tolman-Oppenheimer-Volkoff equations.

These equations are used to set up initial data for relativistic stars, like neutron stars for example.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Step 1.a: import all needed modules from NRPy+:
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends


class TOV_Equations:
    """Construct and store expressions for the TOV equations."""

    def __init__(
        self,
    ) -> None:
        """Set up the right-hand-sides for all four TOV ordinary differential equations, storing them within the class object."""
        self.r_Schw, self.r_iso, self.rho_energy, self.P, self.M, self.M_PI = (
            sp.symbols(
                "r_Schw y[TOVOLA_R_ISO] rho_energy y[TOVOLA_PRESSURE] y[TOVOLA_MASS] M_PI",
                real=True,
            )
        )
        r_Schw = self.r_Schw
        r_iso = self.r_iso
        rho_energy = self.rho_energy
        P = self.P
        M = self.M
        M_PI = self.M_PI
        self.dP_dr = -(
            (rho_energy + P) * ((2 * M) / (r_Schw) + 8 * M_PI * r_Schw * r_Schw * P)
        ) / (r_Schw * 2 * (1 - (2 * M) / (r_Schw)))
        self.dnu_dr = ((2 * M) / (r_Schw) + 8 * M_PI * r_Schw * r_Schw * P) / (
            r_Schw * (1 - (2 * M) / (r_Schw))
        )
        self.dM_dr = 4 * M_PI * r_Schw * r_Schw * rho_energy
        # r_iso == isotropic radius, sometimes called rbar.
        self.dr_iso_dr = (r_iso) / (r_Schw * sp.sqrt(1 - (2 * M) / r_Schw))


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    TOV_eqs = TOV_Equations()
    results_dict = ve.process_dictionary_of_expressions(
        TOV_eqs.__dict__, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        # File basename. If this is set to "trusted_module_test1", then
        #   trusted results_dict will be stored in tests/trusted_module_test1.py
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )
