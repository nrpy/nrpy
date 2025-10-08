"""
Construct the BOBv2 merger-ringdown modes and NQC correction factors for aligned-spin binaries.

Authors: Siddharth Mahesh
sm0193 at mix dot wvu dot edu
Zachariah B. Etienne
zachetie at gmail *dot com
Anuj A. Kankani
aak00009 at mix dot wvu dot edu

The Backwards-One Body (BOB) formalism is a first principles merger-ringdown model that
maps the properties of null congruences in the spacetime of the remnant black hole to
the merger-ringdown waveform of a binary black hole merger.

The Non Quasi-Circular corrections are performed by the SEOBNRv5 model
to ensure that the inspiral waveform matches the NR waveform
at the peak of the (l=2,m=2) strain mode up to second derivatives in amplitude and phase.
See *citation to Anuj's paper*.
See Section IV B of Mahesh, McWilliams, and Etienne, "Spinning Effective-to-Backwards-One Body"
for the BOB-derived NQC corrections.

The modes are expressed in terms of the binary masses (m1, m2), spins (chi1, chi2),
the quasi-normal modes of the remnant black hole (omega_qnm, tau_qnm), and NQC attachment time (t_0).

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
        Compute the BOBv2 waveform.

        This constructor sets up the necessary symbolic variables and expressions
        used in the computation of the BOBv2 waveform. It initializes
        class variables like mass parameters, spin parameters, and various
        coefficients required for the waveform's amplitude and phase.
        The waveform is currently calculated for the (2,2) mode for the sake of
        debugging/sanity checks. Higher mode information will be added in an
        upcoming commit.
        The key outputs of the BOB_v2_waveform_quantities class are:
            - 'h' : the amplitude of the merger-ringdown (l=2,m=2) mode.
            - 'phi' : the phase of the merger-ringdown (l=2,m=2) mode.
            - 'h_t_attach' : the NR-fitted strain amplitude of the (l=2,m=2) mode
                                at the NQC attachment time.
                                (Equation C8 of https://arxiv.org/pdf/2303.18039)
            - 'hddot_t_attach' : the BOB-derived second time derivative of the strain amplitude (l=2,m=2) mode
                                at the NQC attachment time.
            - 'w_t_attach' : the NR-fitted angular frequency of the (l=2,m=2) mode
                                at the NQC attachment time.
                                (Equation C29 of https://arxiv.org/pdf/2303.18039)
            - 'wdot_t_attach' : the BOB-derived first time derivative of the angular frequency of the (l=2,m=2) mode
                                at the NQC attachment time.
            - 't_p_condition' : the BOBv2 peak strain condition.
            - 'Omega_0_condition' : the BOBv2 reference orbital frequency condition.

        :return None:
        """
        (m1, m2, chi1, chi2, Omega_0, omega_qnm, tau_qnm, t_0, t_p, t) = sp.symbols(
            "m1 m2 chi1 chi2 Omega_0 omega_qnm tau_qnm t_0 t_p t", real=True
        )

        (Mf, chif) = sp.symbols("Mf chif", real=True)
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

        coeff1 = f2r(-0.1502782927808088)
        coeff2 = f2r(0.0772527566672910)
        coeff3 = f2r(-0.4915078035114737)
        coeff4 = f2r(0.0844626781097461)
        coeff5 = f2r(-9.9006874770940101)
        coeff6 = f2r(4.3195764159638648)
        coeff7 = f2r(-0.6726797586869773)
        coeff8 = f2r(-14310.0848419331578043)
        coeff9 = f2r(8.9275792268810878)
        coeff10 = f2r(-4.0137844777738332)
        coeff11 = f2r(2.4354075600224134)
        coeff12 = f2r(0.3226728252422290)
        coeff13 = f2r(0.4325654236226389)
        coeff14 = f2r(0.0305373070314258)
        coeff15 = f2r(14310.2825118014534382)

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
        term15 = coeff15 * chip1

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
            + term15
        )

        # Fit for Omega0 for the News (2,2) frequency
        A = f2r(0.33568227)
        B = f2r(0.03450997)
        C = f2r(-0.18763176)
        Omega0 = A * Mf + B * chif + C

        # BOB computation begins here
        Omega_QNM = omega_qnm / 2
        T = (t - t_p) / tau_qnm

        Omega_minus = Omega_QNM**2 - Omega0**2
        Omega_plus = Omega_QNM**2 + Omega0**2
        Omega_orb = sp.sqrt(
            Omega_minus * sp.tanh(T) / sp.Integer(2) + Omega_plus / sp.Integer(2)
        )

        Ap = news22NR_Ap / sp.cosh(T)

        Omega_minus_Q = sp.Abs(Omega_orb - Omega_QNM)
        Omega_minus_0 = sp.Abs(Omega_orb - Omega_0)
        outer = tau_qnm / 2
        inner1 = sp.log((Omega_orb + Omega_QNM) / (Omega_minus_Q))
        inner2 = sp.log((Omega_orb + Omega_0) / (Omega_minus_0))
        Phi_orb = outer * (Omega_QNM * inner1 - Omega_0 * inner2)
        omega_news = -2 * Omega_orb
        phi_news = -2 * Phi_orb

        # Going from news to strain
        # The BOB strain is given as
        # h = H*exp(i * m * Phi_orb)
        # where H is complex and given by the infinite series
        # H = Sum_n H_n
        # First the series is truncated at N = N_provisional

        N_provisional = 0
        i_times_omega = sp.I * omega_news
        H_n = Ap / i_times_omega
        Hbar_n = -Ap / i_times_omega
        H = H_n
        Hbar = Hbar_n
        for _ in range(1, N_provisional + 1):
            H_n = -1 * sp.diff(H_n, t) / i_times_omega
            Hbar_n = sp.diff(Hbar_n, t) / i_times_omega
            H += H_n
            Hbar += Hbar_n
        # Since the merger and NQC handling need amplitude and frequency
        # we use the below two identities
        # amp(z) = sqrt(z * zbar)
        strain_amplitude = sp.sqrt(H * Hbar)
        # frequency(z*exp(i * phi)) = phidot -i * (zbar * zdot - zbardot * z) / (2 * z * zbar)
        strain_frequency = omega_news - sp.I * (
            Hbar * sp.diff(H, t) - sp.diff(Hbar, t) * H
        ) / (2 * H * Hbar)

        self.h_complex = H * sp.exp(sp.I * phi_news)

        # BOB strain peak values will not match the NR values
        # Therefore we calculate the new BOB strain peak values and use those for NQCs
        # Ask Sid for more
        self.t_p_condition = sp.diff(strain_amplitude, t).subs(t, t_0)
        self.t_p_guess = t_0 + tau_qnm * sp.log(Omega_QNM / Omega_0)
        self.h_t_attach = strain_amplitude.subs(t, t_0)
        self.hdot_t_attach = sp.sympify(0)
        hddot = sp.diff(sp.diff(strain_amplitude, t), t)
        self.hddot_t_attach = hddot.subs(t, t_0)
        self.w_t_attach = strain_frequency.subs(t, t_0)
        self.wdot_t_attach = -sp.diff(strain_frequency, t).subs(t, t_0)


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
