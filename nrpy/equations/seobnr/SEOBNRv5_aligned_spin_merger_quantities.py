"""
Construct the native SEOBNRv5 merger-ringdown (l=2,m=2) mode and NQC correction factors for aligned-spin binaries.

Authors: Siddharth Mahesh
sm0193 at mix dot wvu dot edu
Zachariah B. Etienne
zachetie at gmail *dot com

The phenomenological modes documented below are attached as the merger-ringdown component
to inspiral waveforms generated using the Effective-One-Body (EOB) formalism, which maps
the conservative part of the two-body problem to the motion of an effective particle in a deformed metric.
This approach combines post-Newtonian (PN) theory with numerical relativity (NR)
calibrations to model the inspiral phase of a binary black hole merger.

The Non Quasi-Circular corrections are performed by the SEOBNRv5 model
to ensure that the inspiral waveform matches the NR waveform
at the peak of the (l=2,m=2) strain mode up to second derivatives in amplitude and phase.
See Appendix B of https://arxiv.org/pdf/2303.18039 for the full list of NQC terms.

The modes are expressed in terms of the binary masses (m1, m2), spins (chi1, chi2),
the quasi-normal modes of the remnant black hole (omega_qnm, tau_qnm),
and initial conditions derived from the inspiral waveform (t_0, h_0, hdot_0, phi_0, phidot_0)
(see Equations 48-53 and Appendix B of https://arxiv.org/pdf/2303.18039 for the full list of terms).
These are used for documenting the native SEOBNRv5 model in NRPy
and as part of the analysis of the SEBOB formalism  described in Mahesh, McWilliams, and Etienne,
"Spinning Effective-to-Backwards-One Body".
License: BSD 2-Clause
"""

# Step P1: Import needed modules:
import sympy as sp

from nrpy.helpers.float_to_rational import f2r

# The name of this module ("WaveEquation") is given by __name__:
thismodule = __name__


class SEOBNRv5_aligned_spin_merger_quantities:
    """Class for computing the SEOBNRv5 aligned-spin merger gravitational-wave strain and NQC corrections."""

    def __init__(self) -> None:
        """
        Compute the SEONRv5 aligned-spin merger waveform.

        This constructor sets up the necessary symbolic variables and expressions
        used in the computation of the SEOBNRv5 aligned-spin merger strain. It
        initializes class variables like mass parameters, spin parameters, and
        various coefficients required for the waveforms's amplitude and phase.
        The key outputs of the SEOBNRv5_aligned_spin_merger_quantities class are:
            - 'h' : the amplitude of the merger-ringdown (l=2,m=2) mode.
            - 'phi' : the phase of the merger-ringdown (l=2,m=2) mode.
            - 'h_t_attach' : the NR-fitted strain amplitude of the (l=2,m=2) mode
                                at the NQC attachment time.
                                (Equation C8 of https://arxiv.org/pdf/2303.18039)
            - 'hddot_t_attach' : the NR-fitted second time derivative of the strain amplitude (l=2,m=2) mode
                                at the NQC attachment time.
                                (Equation C22 of https://arxiv.org/pdf/2303.18039)
            - 'w_t_attach' : the NR-fitted angular frequency of the (l=2,m=2) mode
                                at the NQC attachment time.
                                (Equation C29 of https://arxiv.org/pdf/2303.18039)
            - 'wdot_t_attach' : the NR-fitted first time derivative of the angular frequency of the (l=2,m=2) mode
                                at the NQC attachment time.
                                (Equation C36 of https://arxiv.org/pdf/2303.18039)

        :return None:
        """
        (
            m1,
            m2,
            chi1,
            chi2,
            omega_qnm,
            tau_qnm,
            t_0,
            h_0,
            hdot_0,
            phi_0,
            phidot_0,
            t,
        ) = sp.symbols(
            "m1 m2 chi1 chi2 omega_qnm tau_qnm t_0 h_0 hdot_0 phi_0 phidot_0 t",
            real=True,
        )
        M = m1 + m2
        q = m1 / m2
        nu = m1 * m2 / M**2
        dm = (q - 1) / (q + 1)
        chiAS = (chi1 - chi2) / 2
        chiS = (chi1 + chi2) / 2
        Shat = chiS + chiAS * dm / (1 - 2 * nu)
        Shat2 = Shat * Shat
        Shat3 = Shat2 * Shat
        Shat4 = Shat3 * Shat
        nu2 = nu * nu
        nu3 = nu2 * nu
        nu4 = nu3 * nu
        # NQC fits
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
        hddot22NR = nu * (
            -f2r(0.00335300225882774906) * nu2
            + f2r(0.00358851415965951012) * nu * Shat
            - f2r(0.00561520957901851664) * nu
            + f2r(0.00038615328462788281) * Shat2
            + f2r(0.00132564277249644174) * Shat
            - f2r(0.00245697909115198589)
        )
        omega22NR = (
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
        omegadot22NR = (
            -f2r(0.23712612963269574795) * nu3
            + f2r(0.07799016321986182443) * nu2 * Shat
            + f2r(0.09221479145462828375) * nu2
            - f2r(0.00839285104015016943) * nu * Shat2
            - f2r(0.02877175649350346628) * nu * Shat
            - f2r(0.03103970973029888253) * nu
            + f2r(0.00061394267373741083) * Shat3
            + f2r(0.0019481328417233967) * Shat2
            + f2r(0.0017051416772119448) * Shat
            - f2r(0.0054839528158373476)
        )
        # MR fits

        c1f = (
            -f2r(4.23824640099272276217) * nu4
            + f2r(1.86797765630645606905) * nu3 * Shat
            + f2r(2.04371156181773017124) * nu3
            + f2r(0.01316051994812337048) * nu2 * Shat2
            - f2r(0.70248812998796750229) * nu2 * Shat
            - f2r(0.40699224653253718298) * nu2
            + f2r(0.06284152952186422558) * nu * Shat3
            + f2r(0.04938791239730106614) * nu * Shat2
            + f2r(0.03388527758147390212) * nu * Shat
            + f2r(0.05358902726316702725) * nu
            - f2r(0.00177674004113037185) * Shat4
            - f2r(0.01890815113487190682) * Shat3
            - f2r(0.01931426101231131093) * Shat2
            - f2r(0.01161150126773277842) * Shat
            + f2r(0.08625435880606084627)
        )

        c2f = (
            -f2r(63.28645899089733006804) * nu4
            + f2r(2.00294725303467924249) * nu3 * Shat
            + f2r(44.33138899436394098075) * nu3
            - f2r(3.55617293922388588712) * nu2 * Shat2
            - f2r(5.58585057654383287939) * nu2 * Shat
            - f2r(9.5295732728313318205) * nu2
            + f2r(1.02187518454288950309) * nu * Shat3
            + f2r(1.97008188121834493245) * nu * Shat2
            + f2r(1.83772448389004638969) * nu * Shat
            + f2r(1.15569522525235401922) * nu
            - f2r(0.20348032514327910047) * Shat3
            - f2r(0.2642970192733161694) * Shat2
            - f2r(0.27076037187561419195) * Shat
            - f2r(0.52876279548305116229)
        )

        d1f = (
            -f2r(28.42370101139921700906) * nu4
            + f2r(4.11346289839689127632) * nu3 * Shat
            + f2r(20.71987362022024470321) * nu3
            + f2r(1.03335215030655280799) * nu2 * Shat2
            - f2r(1.65292430358775521704) * nu2 * Shat
            - f2r(6.07567868511363951001) * nu2
            + f2r(0.04730524488221983515) * nu * Shat3
            - f2r(0.254350860993373451) * nu * Shat2
            + f2r(0.09083410987717309426) * nu * Shat
            + f2r(0.78009259453928059269) * nu
            - f2r(0.01332056979451640664) * Shat4
            - f2r(0.0242033556149483034) * Shat3
            - f2r(0.00784682245346276369) * Shat2
            + f2r(0.1357578010277912528)
        )

        d2f = sp.exp(
            -f2r(352.24938296898454836992) * nu4
            + f2r(9.05730635731021394008) * nu3 * Shat
            + f2r(275.84349920209979245556) * nu3
            + f2r(23.975132253988164166) * nu2 * Shat2
            - f2r(5.26829618908132601973) * nu2 * Shat
            - f2r(81.48331396357356481985) * nu2
            - f2r(3.39885766491276442025) * nu * Shat3
            - f2r(10.06495407151063048445) * nu * Shat2
            + f2r(0.46455322692280037744) * nu * Shat
            + f2r(11.18457585889310479388) * nu
            - f2r(0.1631129108825159213) * Shat4
            + f2r(0.728816370357556087) * Shat3
            + f2r(1.2114999080794128794) * Shat2
            + f2r(0.56269034372727599891) * Shat
            + f2r(0.03570970180918431325)
        )
        # require amplitude and amplitude derivative to match values at maximum
        c1c = 1 / (c1f) * (hdot_0 + h_0 / tau_qnm) * sp.cosh(c2f) ** 2
        c2c = h_0 - 1 / (c1f) * (hdot_0 + h_0 / tau_qnm) * sp.cosh(c2f) * sp.sinh(c2f)

        # require frequency to match value at peak
        d1c = (phidot_0 + omega_qnm) * (1 + d2f) / (d1f * d2f)

        amp = c1c * sp.tanh(c1f * (t - t_0) + c2f) + c2c
        phase = phi_0 - d1c * sp.log((1 + d2f * sp.exp(-d1f * (t - t_0))) / (1 + d2f))
        self.h = amp * sp.exp(-(t - t_0) / tau_qnm)
        self.phi = phase - omega_qnm * (t - t_0)
        self.h_t_attach = h22NR
        self.hdot_t_attach = sp.sympify(0)
        self.hddot_t_attach = hddot22NR
        self.w_t_attach = omega22NR
        self.wdot_t_attach = omegadot22NR


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
        SEOBNRv5_aligned_spin_merger_quantities().__dict__,
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
