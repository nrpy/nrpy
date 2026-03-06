"""
Construct the BOBv2 merger-ringdown modes, as described in https://www.arxiv.org/abs/2510.25012, and NQC-relevant matching quantities for aligned-spin binaries.

Authors: Siddharth Mahesh
sm0193 at mix dot wvu dot edu
Zachariah B. Etienne
zachetie at gmail *dot com
Anuj Kankani
aak00009 at mix dot wvu dot edu

The Backwards-One Body (BOB) formalism is a first principles merger-ringdown model that
maps the properties of null congruences in the spacetime of the remnant black hole to
the merger-ringdown waveform of a binary black hole merger.

The modes are expressed in terms of the binary masses (m1, m2), spins (chi1, chi2),
and the quasi-normal modes of the remnant black hole (omega_qnm, tau_qnm).

We implement an uncalibrated version of the BOBv2 waveform. Because the peak of the BOB strain will not match the NR peak strain time (usually off by 1-2M),
calibrations are required to determine the ideal attachment point.

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
import sympy as sp

from nrpy.helpers.float_to_rational import f2r

# The name of this module is given by __name__:
thismodule = __name__


class BOB_v2_waveform_quantities:
    """Class for computing the BOBv2 gravitational-wave strain and NQC corrections."""

    def __init__(self) -> None:
        """
        Compute the BOBv2 waveform as described in https://www.arxiv.org/abs/2510.25012.

        This constructor sets up the necessary symbolic variables and expressions
        used in the computation of the BOBv2 waveform. It initializes
        class variables like mass parameters, spin parameters, and various
        coefficients required for the waveform's amplitude and phase.
        The waveform is currently calculated for the (2,2) mode for the sake of
        debugging/sanity checks. Higher mode information will be added in an
        upcoming commit.
        The key outputs of the BOB_v2_waveform_quantities class are:
            - 'strain_amp_deriv' : time derivative of the strain amplitude
            - 'h_complex' : complex merger-ringdown strain for the (2,2) mode
            - 'h_t_attach' : the strain amplitude of the BOB merger-ringdown strain for the (l=2,m=2) mode at t_attach
            - 'hdot_t_attach' : the first time derivative of the strain amplitude (l=2,m=2) mode at the t_attach.
            - 'hddot_t_attach' : the second time derivative of the strain amplitude (l=2,m=2) mode at the t_attach.
            - 'w_t_attach' : the angular frequency of the (2,2) mode at the t_attach
            - 'wdot_t_attach' : the first time derivative of the angular frequency of the (l=2,m=2) mode at the t_attach.
        :return None:
        """
        m1, m2, chi1, chi2, omega_qnm, tau_qnm, t_attach, t_p, t = sp.symbols(
            "m1 m2 chi1 chi2 omega_qnm tau_qnm t_attach t_p t", real=True
        )

        M_f, a_f = sp.symbols("M_f a_f", real=True)
        chif = a_f / M_f

        # NQC matching parameters
        M = m1 + m2
        nu = m1 * m2 / M**2

        chi_s = (chi1 + chi2) / 2
        chi_a = (chi1 - chi2) / 2
        delta = (m1 - m2) / M
        chi_eob = chi_s + chi_a * delta / (
            1 - 2 * nu
        )  # Eq (C7 in https://arxiv.org/pdf/2303.18039)

        # Fit for peak amplitude of the news (2,2) mode
        chip1 = chi_eob
        chip2 = chi_eob**2
        chip3 = chi_eob**3
        chip4 = chi_eob**4
        nup1 = nu
        nup2 = nu**2
        nup3 = nu**3
        nup4 = nu**4

        coeff1 = f2r(-0.1690964613229652)
        coeff2 = f2r(0.0804022444632999)
        coeff3 = f2r(-0.1636344956438827)
        coeff4 = f2r(0.0957036572344713)
        coeff5 = f2r(-8.2209626908806417)
        coeff6 = f2r(3.5325407860101454)
        coeff7 = f2r(-0.5562670282104044)
        coeff8 = f2r(0.1923923251183420)
        coeff9 = f2r(-1.5773993140208711)
        coeff10 = f2r(2.6183097662645647)
        coeff11 = f2r(0.9216415432031680)
        coeff12 = f2r(0.4715159651693307)
        coeff13 = f2r(0.4273652956250933)
        coeff14 = f2r(0.0327224221106050)

        term1 = coeff1 * chip3 * nup1
        term2 = coeff2 * chip3
        term3 = coeff3 * chip2 * nup2
        term4 = coeff4 * chip2
        term5 = coeff5 * chip1 * nup3
        term6 = coeff6 * chip1 * nup2
        term7 = coeff7 * chip1 * nup1
        term8 = coeff8 * chip1
        term9 = coeff9 * nup4
        term10 = coeff10 * nup3
        term11 = coeff11 * nup2
        term12 = coeff12 * nup1
        term13 = coeff13
        term14 = coeff14 * chip4

        news22NR_Ap = nu * sp.Abs(
            term1
            + term2
            + term3
            + term4
            + term5
            + term6
            + term7
            + term8
            + term9
            + term10
            + term11
            + term12
            + term13
            + term14
        )

        # Fit for Omega0 for the News (2,2) frequency
        A = f2r(0.33568227)
        B = f2r(0.03450997)
        C = f2r(-0.18763176)
        Omega_0 = A * M_f + B * chif + C

        # BOB computation begins here
        # We build BOB for the news according to equation (1) for the amplitude and (13) for the phase in https://www.arxiv.org/abs/2510.25012
        Omega_QNM = omega_qnm / 2
        T = (t - t_p) / tau_qnm

        # Eq. (10) and (11) in https://www.arxiv.org/abs/2510.25012 with t_0 = -infty
        Omega_minus = Omega_QNM**2 - Omega_0**2
        Omega_plus = Omega_QNM**2 + Omega_0**2
        Omega_orb = sp.sqrt(
            Omega_minus * sp.tanh(T) / sp.Integer(2) + Omega_plus / sp.Integer(2)
        )
        # Eq. (1) in https://www.arxiv.org/abs/2510.25012
        Ap = news22NR_Ap / sp.cosh(T)

        # Eq. (13) in https://www.arxiv.org/abs/2510.25012. Because we take t_0 = -infty, the phase has a closed form expression
        Omega_minus_Q = (
            Omega_QNM - Omega_orb
        )  # equivalent to |Omega_orb - Omega_QNM| since Omega_orb < Omega_QNM
        Omega_minus_0 = Omega_orb - Omega_0
        outer = tau_qnm / 2
        inner1 = sp.log((Omega_orb + Omega_QNM) / (Omega_minus_Q))
        inner2 = sp.log((Omega_orb + Omega_0) / (Omega_minus_0))
        Phi_orb = outer * (Omega_QNM * inner1 - Omega_0 * inner2)
        omega_news = 2 * Omega_orb
        phi_news = 2 * Phi_orb

        # Going from news to strain
        # The BOB strain is given as
        # h = H*exp(-i * m * Phi_orb)
        # where H is complex and given by the infinite series
        # H = Sum_n H_n
        # First the series is truncated at N = N_provisional

        # Prior standalone tests (Kankani & McWilliams) suggest N_provisional = 7 often minimizes
        # news→strain mismatch, but that study did not assess NQC calibrations and some configurations
        # (e.g., chi_eff <= 0) favor smaller N. We therefore default to 0 here until an NQC-optimized
        # choice is established.

        # Going beyond N_provisional=3 causes issues with this sympy implementation due to the complexity of the expressions
        # However, Kankani & McWilliams found that beyond N_provisional=3, the improvements are marginal while the expresssions increase dramatically in complexity
        N_provisional = 0
        i_times_omega = sp.I * omega_news
        H_n = 1.0 / i_times_omega
        Hbar_n = -1.0 / i_times_omega
        H = H_n
        Hbar = Hbar_n
        for _ in range(1, N_provisional + 1):
            H_n = -1 * sp.diff(H_n, t) / i_times_omega
            Hbar_n = sp.diff(Hbar_n, t) / i_times_omega
            H += H_n
            Hbar += Hbar_n
        H *= Ap
        Hbar *= Ap

        strain_amplitude = sp.sqrt(H * Hbar)
        # This is kept as a class vairable because it is used in BOB_v2_setup_peak_attachment to determine the attachment point
        self.strain_amp_deriv = sp.diff(strain_amplitude, t)
        strain_amp_dderiv = sp.diff(strain_amplitude, t, 2)
        # For h = H*exp(-i*phi_news) and defining omega := -d/dt arg(h) = -Im(hdot/h):
        # omega = phidot + i*(Hbar*Hdot - Hbardot*H)/(2*H*Hbar)
        strain_frequency = omega_news + sp.I * (
            Hbar * sp.diff(H, t) - sp.diff(Hbar, t) * H
        ) / (2 * H * Hbar)
        strain_freq_deriv = sp.diff(strain_frequency, t)
        self.h_complex = H * sp.exp(-sp.I * phi_news)

        self.h_t_attach = strain_amplitude.subs(t, t_attach)
        self.hdot_t_attach = self.strain_amp_deriv.subs(t, t_attach)
        self.hddot_t_attach = strain_amp_dderiv.subs(t, t_attach)
        self.w_t_attach = strain_frequency.subs(t, t_attach)
        self.wdot_t_attach = strain_freq_deriv.subs(t, t_attach)


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

    results_dict = ve.process_dictionary_of_expressions(
        BOB_v2_waveform_quantities().__dict__,
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
