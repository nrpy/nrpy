"""
Construct quantities for the SEOBNRv5 merger-ringdown model.

Authors:

Suchindram Dasgupta
sd00113 at mix dot wvu dot edu

Zachariah B. Etienne
zachetie at gmail *dot com

This module implements the key quantities needed to construct the
merger-ringdown waveform in a precessing binary black hole system. It provides:
1. The post-merger Euler rotation angles (alpha, beta, gamma).
2. The precession frequency (omega_prec) of the co-precessing frame.
3. The transformation of quasi-normal mode (QNM) frequencies from the
   remnant's Jf-frame to the co-precessing P-frame.

This implementation is based on Equations 18-22 of
https://arxiv.org/pdf/2303.18046.pdf (SEOBNRv5PHM paper), with the list
of generated modes constrained to match the available NRPy inspiral implementation.

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
import sympy as sp

from nrpy.equations.grhd.Min_Max_and_Piecewise_Expressions import (
    coord_greater_bound,
    coord_less_bound,
)

# The name of this module is given by __name__:
thismodule = __name__


class SEOBNRv5_MergerRingdown_Quantities:
    """
    Compute the SEOBNRv5 merger-ringdown quantities.

    This class sets up the symbolic expressions for the merger-ringdown phase of
    the SEOBNRv5 waveform model, including remnant dynamics and QNM frequencies.
    """

    def __init__(self) -> None:
        """
        Compute the symbolic expressions for the SEOBNRv5 merger-ringdown.

        Inputs (as symbolic variables):
        -----------------------------
        - chi_f_x, chi_f_y, chi_f_z: Cartesian components of the remnant spin vector.
        - L_f_x, L_f_y, L_f_z: Cartesian components of the remnant orbital angular momentum.
        - alpha_match, beta_match, gamma_match: Euler angles at the matching time.
        - t, t_match: Symbolic time and the matching time.

        Outputs (as class attributes):
        ------------------------------
        - self.chi_f_dot_L_f: Symbolic dot product of remnant spin and orbital angular momentum.
        - self.omega_QNM_l_2_m_2, self.omega_QNM_l_2_m_1, self.omega_QNM_l_2_m_n1, self.omega_QNM_l_2_m_n2:
          Input symbols for the l=2 QNM frequencies needed to construct omega_prec.
          These represent the (2,2), (2,1), (2,-1), and (2,-2) modes in the J-frame.
        - self.omega_prec: Precession frequency of the co-precessing frame (Eq. 21).
        - self.alpha_merger_RD, self.beta_merger_RD, self.gamma_merger_RD: Post-merger angles (Eqs. 18-20).
        - self.omega_QNM_J_...: Input symbols for J-frame frequencies (e.g., self.omega_QNM_J_l2_m2).
        - self.omega_QNM_P_...: Output expressions for co-precessing frequencies (e.g., self.omega_QNM_P_l2_m2) (Eq. 22).
        """
        # Step 1: Define fundamental input symbols
        chi_f_x, chi_f_y, chi_f_z = sp.symbols("chi_f_x chi_f_y chi_f_z", real=True)
        L_f_x, L_f_y, L_f_z = sp.symbols("L_f_x L_f_y L_f_z", real=True)
        self.alpha_match, self.beta_match, self.gamma_match = sp.symbols(
            "alpha_match beta_match gamma_match", real=True
        )
        self.t, self.t_match = sp.symbols("t t_match", real=True)

        # Step 2: Calculate the dot product explicitly, storing it as an attribute
        self.chi_f_dot_L_f = chi_f_x * L_f_x + chi_f_y * L_f_y + chi_f_z * L_f_z

        # Step 3: Implement precession frequency (Eq. 21)
        # Define symbols for all four l=2 QNM frequencies needed for omega_prec calculation.
        self.omega_QNM_l_2_m_2 = sp.Symbol("omega_QNM_l_2_m_2", real=True)
        self.omega_QNM_l_2_m_1 = sp.Symbol("omega_QNM_l_2_m_1", real=True)
        self.omega_QNM_l_2_m_n1 = sp.Symbol("omega_QNM_l_2_m_n1", real=True)
        self.omega_QNM_l_2_m_n2 = sp.Symbol("omega_QNM_l_2_m_n2", real=True)

        # Prograde branch: omega_prec = omega_22 - omega_21
        omega_prec_prograde = self.omega_QNM_l_2_m_2 - self.omega_QNM_l_2_m_1
        # Retrograde branch: omega_prec = omega_2,-1 - omega_2,-2
        omega_prec_retrograde = self.omega_QNM_l_2_m_n1 - self.omega_QNM_l_2_m_n2

        prograde_part = omega_prec_prograde * coord_greater_bound(
            self.chi_f_dot_L_f, 0
        ).subs(sp.Function("nrpyAbs"), sp.Abs)
        retrograde_part = omega_prec_retrograde * coord_less_bound(
            self.chi_f_dot_L_f, 0
        ).subs(sp.Function("nrpyAbs"), sp.Abs)
        self.omega_prec = prograde_part + retrograde_part

        # Step 4: Implement post-merger Euler angles (Eqs. 18-20)
        time_diff = self.t - self.t_match
        self.alpha_merger_RD = self.alpha_match + self.omega_prec * time_diff
        self.beta_merger_RD = self.beta_match
        self.gamma_merger_RD = self.gamma_match - self.omega_prec * time_diff * sp.cos(
            self.beta_match
        )

        # Step 5: Loop through modes to calculate co-precessing QNM frequencies (Eq. 22)
        modes = [
            (2, 2),
            (2, 1),
            (3, 3),
            (3, 2),
            (4, 4),
            (4, 3),
            (5, 5),
        ]

        for l, m in modes:
            omega_J_attr_name = f"omega_QNM_J_l{l}_m{m}"
            omega_P_attr_name = f"omega_QNM_P_l{l}_m{m}"

            omega_J_symbol = sp.Symbol(omega_J_attr_name, real=True)
            setattr(self, omega_J_attr_name, omega_J_symbol)

            omega_P_expr = (
                omega_J_symbol
                - m * (1 - sp.Abs(sp.cos(self.beta_match))) * self.omega_prec
            )
            setattr(self, omega_P_attr_name, omega_P_expr)


if __name__ == "__main__":
    import doctest
    import os
    import sys

    project_root = os.path.join(os.path.dirname(__file__), "..", "..", "..")
    sys.path.append(project_root)

    try:
        import nrpy.validate_expressions.validate_expressions as ve
    except ImportError as exc:
        raise ImportError(
            "Could not import NRPy validation module. Please ensure "
            "the script is in the correct directory structure "
            "and __init__.py files are present."
        ) from exc

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    print("Running NRPy expression validation...")
    results_dict = ve.process_dictionary_of_expressions(
        SEOBNRv5_MergerRingdown_Quantities().__dict__,
        fixed_mpfs_for_free_symbols=True,
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )
