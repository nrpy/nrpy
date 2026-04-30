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
        The waveform construction below still uses only the (2,2) BOB mode.
        This code stores higher-mode peak-news fits in `self.newsNR` for future use, but only `(2,2)` is consumed downstream in this module.
        The key outputs of the BOB_v2_waveform_quantities class are:
            - `newsNR["(l,m)"]` for the stored peak-news fits
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
        # these values follow the SEOBNRv5 convention for m1, chi1, m2, chi2 assigment.
        # NQC matching parameters
        M = m1 + m2
        nu = m1 * m2 / M**2

        chi_s = (chi1 + chi2) / 2
        chi_a = (chi1 - chi2) / 2
        delta = (m1 - m2) / M
        chi_eob = chi_s + chi_a * delta / (
            1 - 2 * nu
        )  # Eq (C7 in https://arxiv.org/pdf/2303.18039)

        self.m1 = m1
        self.m2 = m2
        self.chi1 = chi1
        self.chi2 = chi2
        self.M = M
        self.nu = nu
        self.chi_s = chi_s
        self.chi_a = chi_a
        self.delta = delta
        self.chi_eob = chi_eob

        self.newsNR = {}
        modes = [(2, 2), (3, 3), (2, 1), (4, 4), (4, 3), (5, 5), (3, 2)]
        for mode in modes:
            l, m = mode
            self.newsNR.update({f"({l} , {m})": 0})

        self.news_Ap_22()
        self.news_Ap_33()
        self.news_Ap_21()
        self.news_Ap_44()
        self.news_Ap_43()
        self.news_Ap_55()
        self.news_Ap_32()

        # Only the (2,2) mode is implemented for BOBv2 right now.
        newsNR_Ap = self.newsNR["(2 , 2)"]

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
        Ap = newsNR_Ap / sp.cosh(T)

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

    # The following fits for the peak news amplitude for quasi-circular and nonprecessing systems using the V3 SXS catalog
    # For higher modes, we only consider cases where the highest and second-highest resolution levels agree for the peak news amplitude within 0.5%
    # For odd modes we only use cases where delta = (m1-m2)/M > 0.05
    # A list of simulations can be found in tests/BOB_v2_fit_sim_list.md
    def news_Ap_22(self) -> None:
        """Peak news amplitude fit for the (2,2) mode."""
        nu = self.nu
        chi_eob = self.chi_eob
        self.newsNR["(2 , 2)"] = nu * sp.Abs(
            f2r(-0.2099758179309237) * chi_eob**3 * nu
            + f2r(0.0898134166069989) * chi_eob**3
            + f2r(-0.5250725154653605) * chi_eob**2 * nu**2
            + f2r(0.0816049306512471) * chi_eob**2
            + f2r(-12.0901691532595219) * chi_eob * nu**3
            + f2r(5.4706706511122052) * chi_eob * nu**2
            + f2r(-0.8396428499197531) * chi_eob * nu
            + f2r(0.2023493896580952) * chi_eob
            + f2r(-2.2883509124802992) * nu**4
            + f2r(3.1136105962828000) * nu**3
            + f2r(0.8560341070995979) * nu**2
            + f2r(0.4641563023784367) * nu
            + f2r(0.4287100151978567)
            + f2r(0.0368066872502243) * chi_eob**4
        )

    def news_Ap_44(self) -> None:
        """Peak news amplitude fit for the (4,4) mode."""
        nu = self.nu
        chi_s = self.chi_s
        chi_a = self.chi_a
        delta = self.delta
        cad = chi_a * delta
        self.newsNR["(4 , 4)"] = nu * sp.Abs(
            f2r(0.0719107382817115)
            + f2r(4.2714251102339107) * nu
            + f2r(0.1396698884753976) * cad
            + f2r(0.0211970728942547) * chi_s
            + f2r(2903.1233303248659467) * nu**5
            + f2r(1.9754625362184948) * chi_s * nu
            + f2r(-2021.8010919869418558) * nu**4
            + f2r(39.1886106351795505) * chi_s * nu**3
            + f2r(0.0935679366275437) * chi_s**2
            + f2r(-0.1096472808214704) * cad**2 * nu**2
            + f2r(0.0889492950528227) * chi_s * cad
            + f2r(-17.7239654443563381) * chi_s * nu**2
            + f2r(-0.5704327212567502) * chi_s**2 * nu
            + f2r(4.1855358880221187) * chi_s**2 * nu**3
            + f2r(0.0038031784062489) * chi_s**3
            + f2r(0.0574617585831314) * chi_a**2
            + f2r(546.1255161933245290) * nu**3
            + f2r(-72.5837215029558536) * nu**2
            + f2r(0.0310815986621801) * chi_s**2 * cad
            + f2r(-0.2054689578210708) * chi_a**2 * nu
            + f2r(-0.3129620434936395) * cad * nu
        )

    def news_Ap_33(self) -> None:
        """Peak news amplitude fit for the (3,3) mode."""
        nu = self.nu
        chi_s = self.chi_s
        chi_a = self.chi_a
        delta = self.delta
        self.newsNR["(3 , 3)"] = nu * sp.Abs(
            f2r(0.6384426194909012) * delta * nu
            + f2r(0.2542952355821977) * delta
            + f2r(0.1706622629114520) * chi_s * delta
            + f2r(0.0753264566831321) * chi_a
            + f2r(-3.4232149668010052) * chi_a * nu**2
            + f2r(0.1055909587261154) * chi_s**2 * delta
            + f2r(0.1922377901407913) * chi_a**2
            + f2r(0.1820724696424940) * chi_s * chi_a * delta
            + f2r(4.8520577105691292) * delta * nu**3
            + f2r(-0.7024915881908390) * chi_a**2 * nu
            + f2r(0.6522794667725517) * chi_a * nu
            + f2r(-0.0532165822496429) * chi_s * chi_a * nu
            + f2r(1.4114556049561673) * chi_a**3 * delta
            + f2r(1.0418686518147002) * chi_s * delta * nu**2
            + f2r(-1.1769576670699471) * chi_a**3
            + f2r(4.5922498483270955) * chi_a**3 * nu
            + f2r(0.1462996887257039) * chi_s**2 * chi_a * delta
            + f2r(-4.0528225769501720) * chi_a**3 * delta * nu
            + f2r(-0.0580119688642456) * chi_a**2 * delta
            + f2r(-2.2787502829429074) * chi_s * delta * nu**3
        )

    def news_Ap_21(self) -> None:
        """Peak news amplitude fit for the (2,1) mode."""
        nu = self.nu
        chi_s = self.chi_s
        chi_a = self.chi_a
        delta = self.delta
        chi21a = chi_s * delta / (1 - f2r(1.3) * nu) + chi_a
        self.newsNR["(2 , 1)"] = nu * sp.Abs(
            f2r(0.1484235341204765) * delta
            + f2r(-0.4341177039139482) * chi21a * delta
            + f2r(-0.1563736030928236) * chi21a**2 * nu
            + f2r(-0.3945591030310062) * chi_a * nu
            + f2r(-0.7640682856955867) * chi21a**2 * delta
            + f2r(-9.0202567299503364) * chi_s * delta * nu**2
            + f2r(0.8038961066195870) * chi21a**3 * nu
            + f2r(-0.0836258572886352) * chi21a**3 * delta
            + f2r(0.0620966200296146) * chi_a**2 * nu
            + f2r(-0.0155719044439703) * chi_a**2 * delta
            + f2r(0.1838174243321941) * chi_s * delta * nu
            + f2r(-4.0913087391425815) * chi_a * nu**2
            + f2r(0.6661363746524982) * chi21a**2
            + f2r(0.3254466136701785) * chi21a
            + f2r(0.0277640055297656) * delta * nu
            + f2r(-0.0572251308416217) * chi_s**2 * delta * nu
            + f2r(-0.3500926254975283) * chi21a**4 * nu
            + f2r(-0.0014720860009937) * chi21a**4 * delta
            + f2r(-0.0161438168688514) * chi_a**3 * nu
            + f2r(-6.7691739719733510) * chi21a**2 * nu**2
            + f2r(-0.7602742077125189) * chi21a**5 * nu
        )

    def news_Ap_55(self) -> None:
        """Peak news amplitude fit for the (5,5) mode."""
        nu = self.nu
        chi_s = self.chi_s
        chi_a = self.chi_a
        delta = self.delta
        self.newsNR["(5 , 5)"] = nu * sp.Abs(
            f2r(0.1313926678) * delta
            + f2r(0.0693397679) * chi_s * delta
            + f2r(0.1177465005) * chi_a
            + f2r(-5.8928702047) * chi_a * nu**3
            + f2r(-1.0394896698) * chi_a**2 * nu**2
            + f2r(3.1890531795) * chi_s**2 * delta * nu**2
            + f2r(0.1078599864) * chi_s * chi_a * delta
            + f2r(-0.3067283949) * delta * nu
            + f2r(84.1512722160) * delta * nu**4
            + f2r(-32.9500900752) * delta * nu**3
            + f2r(0.0766497896) * chi_a**2
            + f2r(4.1298677252) * delta * nu**2
            + f2r(-0.0519665024) * chi_a**2 * delta
            + f2r(-0.8566261944) * chi_a * nu
            + f2r(3.1432371074) * chi_a * nu**2
            + f2r(-1.0569307007) * chi_s**2 * delta * nu
            + f2r(0.4958660056) * chi_s * delta * nu
            + f2r(18.7848714086) * chi_s * delta * nu**3
            + f2r(0.1068696584) * chi_s**2 * delta
            + f2r(-6.9871345213) * chi_s * delta * nu**2
            + f2r(-0.3879987253) * chi_s * chi_a * delta * nu
        )

    def news_Ap_32(self) -> None:
        """Peak news amplitude fit for the (3,2) mode."""
        nu = self.nu
        chi_s = self.chi_s
        chi_a = self.chi_a
        delta = self.delta
        cad = chi_a * delta
        self.newsNR["(3 , 2)"] = nu * sp.Abs(
            f2r(0.0577545572673071)
            + f2r(11.9406845223270235) * cad * nu**3
            + f2r(17.0580436575595620) * chi_s * nu**2
            + f2r(109.9438513270014681) * nu**4
            + f2r(-30.0228076416933121) * nu**3
            + f2r(-0.1443807511730615) * chi_s**2
            + f2r(-0.0453626058575554) * cad
            + f2r(-52.3889618079092614) * chi_s * nu**3
            + f2r(-1.4284634631896036) * chi_s * nu
            + f2r(2.8546941351622510) * cad**2 * nu**2
            + f2r(-0.2271693889602418) * chi_s * cad * nu
            + f2r(0.3441847817715949) * cad**2 * nu
            + f2r(0.0089555404141321) * chi_s**3
            + f2r(0.1472774555772416) * nu
            + f2r(0.1062807860491596) * chi_s * cad
            + f2r(0.0950917657779650) * chi_s**2 * cad
            + f2r(0.1741453741496004) * chi_s * cad**2
            + f2r(0.1146704472597599) * cad**3
            + f2r(-4.7385633203438493) * chi_s**2 * nu**2
            + f2r(1.8285010762094736) * chi_s**2 * nu
            + f2r(30.4155122363908923) * chi_s * nu**4
        )

    def news_Ap_43(self) -> None:
        """Peak news amplitude fit for the (4,3) mode."""
        nu = self.nu
        chi_s = self.chi_s
        chi_a = self.chi_a
        delta = self.delta
        self.newsNR["(4 , 3)"] = nu * sp.Abs(
            f2r(0.0299884197215471) * delta
            + f2r(0.0851154731131591) * chi_a**2 * nu
            + f2r(1.0281190429001357) * chi_s * delta * nu
            + f2r(17.0363543490113400) * delta * nu**3
            + f2r(-5.7379134076446610) * delta * nu**2
            + f2r(0.0035788067699685) * chi_s**2 * delta
            + f2r(0.3496042267046198) * chi_a**3 * delta * nu
            + f2r(0.0462878750795055) * chi_s * chi_a * delta
            + f2r(0.5170339896107963) * delta * nu
            + f2r(-0.0080442896806339) * chi_a**3
            + f2r(-0.0087657741673305) * chi_a**3 * delta
            + f2r(-2.7102795229228880) * chi_s * delta * nu**2
            + f2r(-0.0695995217823558) * chi_s * delta
            + f2r(0.0018629241397624) * chi_a
            + f2r(-0.0629034591329580) * chi_a * nu**2
            + f2r(0.0569224478624692) * chi_s * chi_a**2 * delta
            + f2r(0.4590840469943437) * chi_s**2 * delta * nu**2
            + f2r(0.0397201775916978) * chi_s**2 * chi_a * delta
        )


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
    obj = BOB_v2_waveform_quantities()
    test_dict = obj.__dict__.copy()
    if "newsNR" in test_dict:
        del test_dict["newsNR"]
    test_dict.update(
        {
            "news_Ap_22": obj.newsNR["(2 , 2)"],
            "news_Ap_33": obj.newsNR["(3 , 3)"],
            "news_Ap_21": obj.newsNR["(2 , 1)"],
            "news_Ap_44": obj.newsNR["(4 , 4)"],
            "news_Ap_43": obj.newsNR["(4 , 3)"],
            "news_Ap_55": obj.newsNR["(5 , 5)"],
            "news_Ap_32": obj.newsNR["(3 , 2)"],
        }
    )
    results_dict = ve.process_dictionary_of_expressions(
        test_dict,
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
