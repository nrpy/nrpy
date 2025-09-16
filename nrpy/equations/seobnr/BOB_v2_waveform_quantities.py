"""
Construct the BOBv2 merger-ringdown modes and NQC correction factors for aligned-spin binaries.

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
        # NQC matching parameters
        M = m1 + m2
        nu = m1 * m2 / M**2
        Shat = (m1 * m1 * chi1 + m2 * m2 * chi2) / (m1 * m1 + m2 * m2)
        Shat2 = Shat * Shat
        Shat3 = Shat2 * Shat
        Shat4 = Shat3 * Shat
        nu2 = nu * nu
        nu3 = nu2 * nu
        nu4 = nu3 * nu
        h22NR = nu * sp.Abs(
            f2r(71.97969776036882194603) * nu4
            - f2r(13.35761402231352157344) * nu3 * Shat
            - f2r(46.87585958426210908101) * nu3
            + f2r(0.61988944517825661507) * nu2 * Shat2
            + f2r(7.19426416189229733789) * nu2 * Shat
            + f2r(12.44040490932310127903) * nu2
            + f2r(0.43014673069078152023) * nu * Shat3
            - f2r(1.74313546783413597652) * nu * Shat
            - f2r(0.86828935763242798274) * nu
            - f2r(0.08493901280736430859) * Shat3
            - f2r(0.02082621429567295401) * Shat2
            + f2r(0.18693991146784910695) * Shat
            + f2r(1.46709663479911811557)
        )
        omega22NR = -(
            f2r(5.89352329617707670906) * nu4
            + f2r(3.75145580491965446868) * nu3 * Shat
            - f2r(3.34930536209472551334) * nu3
            - f2r(0.97140932625194231775) * nu2 * Shat2
            - f2r(1.69734302394369973577) * nu2 * Shat
            + f2r(0.28539204856044564362) * nu2
            + f2r(0.2419483723662931296) * nu * Shat3
            + f2r(0.51801427018052081941) * nu * Shat2
            + f2r(0.25096450064948544467) * nu * Shat
            - f2r(0.31709602351033533418) * nu
            - f2r(0.01525897668158244028) * Shat4
            - f2r(0.06692658483513916345) * Shat3
            - f2r(0.08715176045684569495) * Shat2
            - f2r(0.09133931944098934441) * Shat
            - f2r(0.2685414392185025978)
        )
        # BOB computation begins here
        Omega_qnm = omega_qnm / 2
        T = (t - t_p) / tau_qnm
        T_0 = (t_0 - t_p) / tau_qnm
        # BOB news quantities
        # Note:
        # While the derivation below follows Kankani et al
        # they fix constants of integration/continuity by requiring
        # t_0 = - infinity
        # A_p = peak news amplitude from NR
        # Omega_0 = asymptotic initial orbital frequency
        # For IMR attachment it can be more optimal
        # to define the constants of integration/continuity
        # as follows
        # t_0 = time of peak strain
        # A_p = rescaled peak news amplitude such that
        #         BOB strain amplitude matches NR fitted quantity
        #         at t_0
        # Omega_0 = rescaled initial orbital frequency such that
        #         BOB phase frequency matches NR fitted quantity
        #         at t_0
        # The condition for A_p can be expressed as a closed form.
        # The definition of t_0 implies that
        # Omega_0 and T_0 have to be evaluated numerically
        # given the NQC matching parameters.
        # A_noap = The BOB news amplitude without the A_p scale factor
        A_noap = 1 / sp.cosh(T)
        k = (Omega_qnm**2 - Omega_0**2) / (1 - sp.tanh(T_0))
        Omega_orb = (Omega_0**2 + k * (sp.tanh(T) - sp.tanh(T_0))) ** sp.Rational(1, 2)
        Phi_orb = -tau_qnm * (
            Omega_qnm * sp.atanh(Omega_orb / Omega_qnm)
            - Omega_0 * sp.atanh(Omega_0 / Omega_orb)
        )
        omega_news = 2 * Omega_orb
        # Going from news to strain
        # The BOB strain is given as
        # h = H*exp(i * m * Phi_orb)
        # where H is complex and given by the infinite series
        # H = Sum_n H_n
        # First the series is truncated at N = N_provisional
        N_provisional = 2
        i_times_omega = sp.I * omega_news
        H_n = A_noap / i_times_omega
        Hbar_n = -A_noap / i_times_omega
        H = H_n
        Hbar = Hbar_n
        for _ in range(1, N_provisional + 1):
            H_n = -1 * sp.diff(H_n, t) / i_times_omega
            Hbar_n = sp.diff(Hbar_n, t) / i_times_omega
            H += H_n
            Hbar += Hbar_n
        # Since the merger and NQC handling need amplitude and phase
        # we use the below three identities
        # amp(z) = sqrt(z * zbar)
        strain_amplitude_noap = sp.sqrt(H * Hbar)
        # phase(z*exp(i * phi)) = phi + atan((z - zbar)/(i * (z + zbar)))
        strain_phase = 2 * Phi_orb + sp.atan((H - Hbar) / (sp.I * (H + Hbar)))
        # frequency(z*exp(i * phi)) = phidot -i * (zbar * zdot - zbardot * z) / (2 * z * zbar)
        strain_frequency = omega_news - sp.I * (
            Hbar * sp.diff(H, t) - sp.diff(Hbar, t) * H
        ) / (2 * H * Hbar)
        # define A_p based on continuity at strain peak
        A_p = h22NR / strain_amplitude_noap.subs(t, t_0)
        self.h = A_p * strain_amplitude_noap
        self.phi = strain_phase
        # mostly trivial
        # unlike in BOBv1 where t_0 - t_p has a closed form,
        # we will need to solve a non-linear equation to
        # get t_0 - t_p
        self.t_p_condition = sp.diff(strain_amplitude_noap, t).subs(t, t_0)
        self.Omega_0_condition = strain_frequency.subs(t, t_0) - omega22NR
        self.h_t_attach = h22NR
        self.hdot_t_attach = sp.sympify(0)
        hddot = sp.diff(sp.diff(self.h, t), t)
        self.hddot_t_attach = hddot.subs(t, t_0)
        self.w_t_attach = omega22NR
        self.wdot_t_attach = sp.diff(strain_frequency, t).subs(t, t_0)


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
