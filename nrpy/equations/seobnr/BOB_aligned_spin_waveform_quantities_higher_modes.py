"""
Construct the BOB merger-ringdown (l,m) = (2,2), (3,3), (4,4), (5,5), (2,1), (3,2), (4,3) modes and NQC correction factors for aligned-spin binaries.

Authors: Siddharth Mahesh
sm0193 at mix dot wvu dot edu
Zachariah B. Etienne
zachetie at gmail *dot com

The Backwards-One Body (BOB) formalism is a first principles merger-ringdown model that
maps the properties of null congruences in the spacetime of the remnant black hole to
the merger-ringdown waveform of a binary black hole merger.

The Non Quasi-Circular corrections are performed by the SEOBNRv5 model
to ensure that the inspiral waveform matches the NR waveform
at the peak of the (l=2,m=2) strain mode up to second derivatives in amplitude and phase.
See Appendix B of https://arxiv.org/pdf/2303.18039 for the strain and angular frequency terms.
See Section IV B of Mahesh, McWilliams, and Etienne, "Spinning Effective-to-Backwards-One Body"
for the BOB-derived NQC corrections.

The modes are expressed in terms of the initial amplitudes and frequencies (h_lm, omega_lm),
the quasi-normal modes of the remnant black hole (omega_lm_qnm, tau_lm_qnm), and NQC attachment time (t_0).
(see Equations 24-31 of Mahesh, McWilliams, and Etienne, "Spinning Effective-to-Backwards-One Body"
for the full list of terms).

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
import sympy as sp

# The name of this module is given by __name__:
thismodule = __name__


class BOB_aligned_spin_waveform_quantities_higher_modes:
    """Class for computing the BOB aligned-spin gravitational-wave strain and NQC corrections."""

    def __init__(self) -> None:
        """
        Compute the BOB aligned-spin waveform.

        This constructor sets up the necessary symbolic variables and expressions
        used in the computation of the aligned-spin BOB strain. It
        initializes class variables like mass parameters, spin parameters, and
        various coefficients required for the waveforms's amplitude and phase.
        The key outputs of the BOB_aligned_spin_waveform_quantities_higher_modes class are:
            - 'h' : the amplitude of the merger-ringdown (l,m) mode.
            - 'phi' : the phase of the merger-ringdown (l,m) mode.
            - 'h_t_attach' : the NR-fitted strain amplitude of the (l,m) mode
                                at the NQC attachment time.
            - 'hddot_t_attach' : the BOB-derived second time derivative of the strain amplitude (l,m) mode
                                at the NQC attachment time.
            - 'w_t_attach' : the NR-fitted angular frequency of the (l,m) mode
                                at the NQC attachment time.
                                (Equation C29 of https://arxiv.org/pdf/2303.18039)
            - 'wdot_t_attach' : the BOB-derived first time derivative of the angular frequency of the (l=2,m=2) mode
                                at the NQC attachment time.

        :return None:
        """
        self.modes = [(2, 2), (3, 3), (4, 4), (5, 5), (2, 1), (3, 2), (4, 3)]
        self.h_lm = {}
        self.phi_lm = {}
        self.hdot_t_attach_lm = {}
        self.hddot_t_attach_lm = {}
        self.wdot_t_attach_lm = {}

        t, t_0 = sp.symbols("t t_0", real=True)
        self.modewise_input_symbols = {}
        for lm in self.modes:
            self.modewise_input_symbols[lm] = [
                sp.Symbol(f"h_{lm[0]}{lm[1]}_t_attach", real=True),
                sp.Symbol(f"omega_{lm[0]}{lm[1]}_t_attach", real=True),
                sp.Symbol(f"omega_{lm[0]}{lm[1]}_qnm", real=True),
                sp.Symbol(f"tau_{lm[0]}{lm[1]}_qnm", real=True),
            ]

        for lm in self.modes:
            m = lm[1]
            h_0 = self.modewise_input_symbols[lm][0]
            omega_0 = self.modewise_input_symbols[lm][1]
            omega_qnm = self.modewise_input_symbols[lm][2]
            tau_qnm = self.modewise_input_symbols[lm][3]
            Omega_0 = omega_0 / m
            Omega_qnm = omega_qnm / m
            t_p = t_0 - 2 * tau_qnm * sp.log(Omega_0 / Omega_qnm)
            A_p = h_0 * (omega_0**2) * sp.cosh((t_0 - t_p) / tau_qnm)
            k = (Omega_qnm**4 - Omega_0**4) / (1 - sp.tanh((t_0 - t_p) / tau_qnm))
            kappa_m = Omega_0 * Omega_0 / Omega_qnm
            kappa_p = Omega_qnm
            Omega = (
                Omega_0**4
                + k * (sp.tanh((t - t_p) / tau_qnm) - sp.tanh((t_0 - t_p) / tau_qnm))
            ) ** (1 / 4)
            Omega_0_over_kappa_p = Omega_0 / kappa_p
            Omega_0_over_kappa_m = Omega_0 / kappa_m
            Omega_over_kappa_p = Omega / kappa_p
            Omega_over_kappa_m = Omega / kappa_m
            arctanh_m = (
                0.5
                * kappa_m
                * tau_qnm
                * sp.log(
                    (1 + Omega_over_kappa_m)
                    * (1 - Omega_0_over_kappa_m)
                    / ((1 - Omega_over_kappa_m) * (1 + Omega_0_over_kappa_m))
                )
            )
            arctanh_p = (
                0.5
                * kappa_p
                * tau_qnm
                * sp.log(
                    (1 + Omega_over_kappa_p)
                    * (1 - Omega_0_over_kappa_p)
                    / ((1 - Omega_over_kappa_p) * (1 + Omega_0_over_kappa_p))
                )
            )
            arctan_m = (
                kappa_m
                * tau_qnm
                * (sp.atan2(Omega, kappa_m) - sp.atan2(Omega_0, kappa_m))
            )
            arctan_p = (
                kappa_p
                * tau_qnm
                * (sp.atan2(Omega, kappa_p) - sp.atan2(Omega_0, kappa_p))
            )
            self.h_lm[lm] = (A_p / m**2 / (Omega**2)) * (
                1 / sp.cosh((t - t_p) / tau_qnm)
            )
            Phi = arctan_p + arctanh_p - arctan_m - arctanh_m
            self.phi_lm[lm] = m * Phi
            self.hdot_t_attach_lm[lm] = sp.sympify(0)
            hddot = sp.diff(sp.diff(self.h_lm[lm], t), t)
            self.hddot_t_attach_lm[lm] = hddot.subs(t, t_0)
            Omegadot = sp.diff(Omega, t)
            self.wdot_t_attach_lm[lm] = m * Omegadot.subs(t, t_0)


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

    wf = BOB_aligned_spin_waveform_quantities_higher_modes()
    expressions_dict = {}
    for ell, em in wf.modes:
        label = f"{ell}{em}"
        ellm = (ell, em)
        expressions_dict[f"h_{label}"] = wf.h_lm[ellm]
        expressions_dict[f"phi_{label}"] = wf.phi_lm[ellm]
        expressions_dict[f"hdot_t_attach_{label}"] = wf.hdot_t_attach_lm[ellm]
        expressions_dict[f"hddot_t_attach_{label}"] = wf.hddot_t_attach_lm[ellm]
        expressions_dict[f"wdot_t_attach{label}"] = wf.wdot_t_attach_lm[ellm]

    results_dict = ve.process_dictionary_of_expressions(
        expressions_dict,
        fixed_mpfs_for_free_symbols=True,
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        # File basename. If this is set to "trusted_module_test1", then
        #   trusted results_dict will be stored in tests/trusted_module_test1.py
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )
