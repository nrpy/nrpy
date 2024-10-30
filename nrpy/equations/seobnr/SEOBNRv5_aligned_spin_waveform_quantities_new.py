"""
Construct symbolic expression for the SEOBNRv5 aligned-spin gravitational-wave strain and flux.

Author: Siddharth Mahesh

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
from typing import Any, Dict, List

import sympy as sp

# The name of this module ("WaveEquation") is given by __name__:
thismodule = __name__


def complex_mult(z1: List[Any], z2: List[Any]) -> List[Any]:
    """
    Multiply two complex numbers given as list of real and imaginary parts.

    This functions takes two lists containing the real and imaginary part of a complex number
    and returns a list with the real and imaginary part of the resulting multiple.

    :param z1: Complex number 1 as list [Real(z1),Imag(z1)]
    :param z2: Complex number 2 as list [Real(z2),Imag(z2)]
    :return: Complex number z1 x z2 as list [Real(z1*z2),Imag(z1*z2)]

    >>> z1 = [1,2]
    >>> z2 = [3,5]
    >>> complex_mult(z1,z2)
    [-7, 11]

    >>> import sympy as sp
    >>> x1 , y1 , x2 , y2 = sp.symbols('x1 y1 x2 y2',real = True)
    >>> z1 = [x1,y1]
    >>> z2 = [x2,y2]
    >>> complex_mult(z1,z2)
    [x1*x2 - y1*y2, x1*y2 + x2*y1]
    """
    # complex multiplication
    # z1 = x1 + I*y1
    # z2 = x2 + I*y2
    # z1*z2 = x1*x2 - y1*y2 + I*(x1*y2 + x2*y1)

    return [z1[0] * z2[0] - z1[1] * z2[1], z1[0] * z2[1] + z1[1] * z2[0]]


def f2r(input_float: float) -> sp.Rational:
    """
    Convert a floating-point number to a high-precision rational number.

    This function takes a floating-point number, converts it to a string,
    and appends 60 zeros to increase the precision of the conversion to a rational number.

    :param input_float: The floating-point number to convert.
    :return: A sympy Rational number with high precision.

    >>> f2r(0.1)
    1/10
    >>> f2r(1.5)
    3/2
    >>> f2r(2.0)
    2
    """
    # Convert the input float to a string
    float_as_string = str(input_float)

    # Ensure the string has a decimal point
    if "." not in float_as_string:
        float_as_string = f"{float_as_string}."

    # Append 60 zeros after the decimal of the floating point number to increase precision
    return sp.Rational(float_as_string + "0" * 60)


class SEOBNRv5_aligned_spin_waveform_quantities:
    """Class for computing the SEOBNRv5 aligned-spin gravitational-wave strain and flux."""

    def __init__(self) -> None:
        """
        Compute the SEOBNRv5 aligned-spin Hamiltonian.

        This constructor sets up the necessary symbolic variables and expressions
        used in the computation of the aligned-spin SEOBNRv5 flux and strain. It
        initializes class variables like mass parameters, spin parameters, and
        various coefficients required for the waveforms's multipole modes.

        Inputs: 'm1', 'm2', 'r', 'phi', 'prstar', 'pphi', 'chi1', 'chi2', 'Hreal', 'Omega' and 'Omega_circ'
        Outputs: 'flux' and 'hlms'
        """
        m1, m2, r, phi, prstar, pphi, chi1, chi2,Hreal,Omega,Omega_circ = sp.symbols('m1 m2 r phi prstar pphi chi1 chi2 Hreal Omega Omega_circ',real = True)
        gsl_gamma_wrapper = sp.Function("SEOBNRv5_aligned_spin_gamma_wrapper")
        M = m1+m2
        eta = m1*m2/(M**2)
        TINYDOUBLE = 1e-100
        delta = (m1 - m2)/M
        min_delta_1em14 = sp.Rational(1,2)*(delta - 1e-14 + 0 - sp.Abs(delta - 1e-14 - 0))
        max_delta_1em14 = sp.Rational(1,2)*(delta - 1e-14 + TINYDOUBLE + 0 + sp.Abs(delta - 1e-14 + TINYDOUBLE - 0))
        noneqcond = max_delta_1em14/(delta - 1e-14+TINYDOUBLE)
        eqcond = min_delta_1em14/(delta - 1e-14-TINYDOUBLE)
        deltainvertible = delta*noneqcond + 1*eqcond
        deltainv = 1/deltainvertible
        chi_S = sp.Rational(1,2)*(chi1 + chi2)
        chi_A = sp.Rational(1,2)*(chi1 - chi2)
        vomega = Omega**(sp.Rational(1,3))
        vphi = Omega*(Omega_circ**(-sp.Rational(2,3)))
        Heff = (Hreal**2 - 1)/(2*eta) + 1
        eulerlog5 = sp.EulerGamma + sp.log(2*5*vomega)
        eulerlog4 = sp.EulerGamma + sp.log(2*4*vomega)
        eulerlog3 = sp.EulerGamma + sp.log(2*3*vomega)
        eulerlog2v2 = sp.EulerGamma + sp.log(2*2*vomega**2)
        eulerlog2 = sp.EulerGamma + sp.log(2*2*vomega)
        eulerlog1v2 = sp.EulerGamma + sp.log(2*1*vomega**2)
        eulerlog1 = sp.EulerGamma + sp.log(2*1*vomega)
        khat1 = 1*Omega*Hreal
        khat2 = 2*Omega*Hreal
        khat3 = 3*Omega*Hreal
        khat4 = 4*Omega*Hreal
        khat5 = 5*Omega*Hreal
        khat6 = 6*Omega*Hreal
        khat7 = 7*Omega*Hreal
        khat8 = 8*Omega*Hreal
        delta55=((sp.Rational(96875,131250)+sp.Rational(857528,131250)*eta)*(Omega*Hreal)/(1-2*eta)+sp.Rational(3865,429)*sp.pi*(Omega*Hreal)**2+((sp.Rational(-7686949127,31783752)+sp.Rational(954500400,31783752)*sp.pi**2)*(Omega*Hreal)**3))
        delta43=(((sp.Rational(4961,810)*eta+sp.Rational(3,5))*(Omega*Hreal))/(1-2*eta)+sp.Rational(1571,385)*sp.pi*(Omega*Hreal)**2)
        delta44=((sp.Rational(112,120)+sp.Rational(219,120)*eta)*(Omega*Hreal)+sp.Rational(25136,3465)*(Omega*Hreal)**2+(sp.Rational(201088,10395)*sp.pi**2-sp.Rational(55144,375))*(Omega*Hreal)**3)
        delta32=(((sp.Rational(11,5)*eta+sp.Rational(2,3))*(Hreal*Omega))/(1-3*eta)+sp.Rational(52,21)*sp.pi*(Hreal*Omega)**2+((sp.Rational(208,63)*sp.pi**2)-sp.Rational(9112,405))*(Hreal*Omega)**3)
        delta33=(sp.Rational(13,10)*(Hreal*Omega)+sp.Rational(39,7)*(Hreal*Omega)**2+(-sp.Rational(227827,3000)+sp.Rational(78,7)*sp.pi**2)*(Hreal*Omega)**3-sp.Rational(80897,2430)*eta*vomega**5)
        delta21=(sp.Rational(2,3)*Omega*Hreal+(sp.Rational(107,105)*sp.pi)*(Omega*Hreal)**2+((sp.Rational(214,315)*sp.pi**2)-(sp.Rational(272,81)))*(Omega*Hreal)**3-(sp.Rational(25,2)*eta*vomega**5))
        delta22=(sp.Rational(7,3)*Omega*Hreal+(Omega*Hreal)**2*(((sp.Rational(8,3)*eta-sp.Rational(4,3))*chi_S)-(sp.Rational(4,3)*delta*chi_A)+(sp.Rational(428,105)*sp.pi))+(Omega*Hreal)**3*((sp.Rational(1712,315)*sp.pi**2)-(sp.Rational(2203,81)))-24*eta*vomega**5)
        fspin55_limit=(((sp.Rational(-70,3)*eta/((-1+2*eta))+sp.Rational(110,3)*eta**2/((-1+2*eta))+sp.Rational(10,3)/((-1+2*eta)))*(chi_A))*vomega**3+((-5.0/(-1.0+2.0*eta)+30.0*eta/(-1.0+2.0*eta)-40.0*eta**2/(-1.0+2.0*eta))*chi_A*chi_S)*vomega**4)
        fspin55=(((-70*eta/(3*(-1+2*eta))+110*eta**2/(3*(-1+2*eta))+10/(3*(-1+2*eta)))*(chi_A*deltainv)+(10/(3*(-1+2*eta))-10*eta/(-1+2*eta)+10*eta**2/(-1+2*eta))*chi_S)*vomega**3+(sp.Rational(5,2)*chi_S**2+(-5.0/(-1.0+2.0*eta)+30.0*eta/(-1.0+2.0*eta)-40.0*eta**2/(-1.0+2.0*eta))*chi_A*chi_S*deltainv+(-5.0/(2.0*(-1.0+2.0*eta))+15.0*eta/(-1.0+2.0*eta)-20.0*eta**2/(-1.0+2.0*eta))*chi_A**2)*vomega**4)
        fspin41_limit=(sp.Rational(5,2)*eta*vomega/(1-2*eta)*(-chi_A))
        fspin41=(sp.Rational(5,2)*eta*vomega/(1-2*eta)*(chi_S-chi_A*deltainv))
        fspin43_limit=(vomega/(1-2*eta)*(-sp.Rational(5,2)*eta*chi_A)+vomega**3/(1-2*eta)*((sp.Rational(887,44)*eta-sp.Rational(3143,132)*eta**2)*chi_A)+vomega**4/(1-2*eta)*((sp.Rational(137,6)*eta**2-18*eta+3)*chi_A*chi_S))
        fspin43=(vomega/(1-2*eta)*(sp.Rational(5,2)*eta*chi_S-sp.Rational(5,2)*eta*chi_A*deltainv)+vomega**3/(1-2*eta)*((sp.Rational(887,44)*eta-sp.Rational(3143,132)*eta**2)*chi_A*deltainv+(-sp.Rational(529,132)*eta**2-sp.Rational(667,44)*eta)*chi_S)+vomega**4/(1-2*eta)*((12*eta**2-sp.Rational(37,3)*eta+sp.Rational(3,2))*chi_A**2+(sp.Rational(137,6)*eta**2-18*eta+3)*chi_A*chi_S*deltainv+(sp.Rational(35,6)*eta**2+sp.Rational(1,3)*eta+sp.Rational(3,2))*chi_S**2))
        fspin31_limit = vomega**3 * (-chi_A * 5.0 / 8.0)
        fspin31=vomega**3*(chi_A*deltainv*(-4.0+11.0*eta)+chi_S*(-4.0+13.0*eta))/(2.0)
        fspin33_limit=vomega**3*((sp.Rational(19,2)*eta-2)*chi_A)+vomega**4*((3-12*eta)*(chi_A*chi_S))+vomega**5*((sp.Rational(407,30)*eta**2-sp.Rational(593,60)*eta+sp.Rational(2,3))*chi_A)+vomega**6*((44*eta**2-eta-sp.Rational(7,2))*(chi_A*chi_S))+sp.I*(Omega*Hreal)**2*(sp.Rational(7339,540)*eta-sp.Rational(81,20))*chi_A
        fspin33amp_limit=vomega**3*((sp.Rational(19,2)*eta-2)*chi_A)+vomega**4*((3-12*eta)*(chi_A*chi_S))+vomega**5*((sp.Rational(407,30)*eta**2-sp.Rational(593,60)*eta+sp.Rational(2,3))*chi_A)+vomega**6*((44*eta**2-eta-sp.Rational(7,2))*(chi_A*chi_S))
        fspin33=(vomega**3*(((sp.Rational(19,2)*eta-2)*chi_A*deltainv)+((sp.Rational(5,2)*eta-2)*chi_S))+vomega**4*((sp.Rational(3,2)-6*eta)*chi_A**2+(3-12*eta)*(chi_A*chi_S*deltainv)+sp.Rational(3,2)*chi_S**2)+vomega**5*(((sp.Rational(407,30)*eta**2-sp.Rational(593,60)*eta+sp.Rational(2,3))*chi_A*deltainv)+((sp.Rational(241,30)*eta**2+sp.Rational(11,20)*eta+sp.Rational(2,3))*chi_S))+vomega**6*((-12*eta**2+sp.Rational(11,2)*eta-sp.Rational(7,4))*chi_A**2+(44*eta**2-eta-sp.Rational(7,2))*(chi_A*chi_S*deltainv)+(6*eta**2-sp.Rational(27,2)*eta-sp.Rational(7,4))*chi_S**2)+sp.I*((Omega*Hreal)**2*(sp.Rational(7339,540)*eta-sp.Rational(81,20))*chi_A*deltainv+(sp.Rational(593,108)*eta-sp.Rational(81,20))*chi_S))
        fspin33amp=(vomega**3*(((sp.Rational(19,2)*eta-2)*chi_A*deltainv)+((sp.Rational(5,2)*eta-2)*chi_S))+vomega**4*((sp.Rational(3,2)-6*eta)*chi_A**2+(3-12*eta)*(chi_A*chi_S*deltainv)+sp.Rational(3,2)*chi_S**2)+vomega**5*(((sp.Rational(407,30)*eta**2-sp.Rational(593,60)*eta+sp.Rational(2,3))*chi_A*deltainv)+((sp.Rational(241,30)*eta**2+sp.Rational(11,20)*eta+sp.Rational(2,3))*chi_S))+vomega**6*((-12*eta**2+sp.Rational(11,2)*eta-sp.Rational(7,4))*chi_A**2+(44*eta**2-eta-sp.Rational(7,2))*(chi_A*chi_S*deltainv)+(6*eta**2-sp.Rational(27,2)*eta-sp.Rational(7,4))*chi_S**2))
        fspin21_limit=(-sp.Rational(3,2)*vomega*(chi_A)+vomega**3*((sp.Rational(131,84)*eta+sp.Rational(61,12))*(chi_A))+vomega**4*(+(sp.Rational(21,2)*eta-6)*(chi_A*chi_S))+vomega**5*((sp.Rational(-703,112)*eta**2+sp.Rational(8797,1008)*eta-sp.Rational(81,16))*(chi_A)+(sp.Rational(3,4)-3*eta)*(chi_A**3)+(sp.Rational(9,4)-6*eta)*(chi_A*chi_S**2))+vomega**6*((sp.Rational(9487,504)*eta**2-sp.Rational(1636,21)*eta+sp.Rational(4163,126))*(chi_A*chi_S)))
        fspin21=(-sp.Rational(3,2)*vomega*(chi_A*deltainv+chi_S)+vomega**3*((sp.Rational(131,84)*eta+sp.Rational(61,12))*(chi_A*deltainv)+(sp.Rational(79,84)*eta+sp.Rational(61,12))*chi_S)+vomega**4*((-2*eta-3)*chi_A**2+(sp.Rational(21,2)*eta-6)*(chi_A*chi_S*deltainv)+(sp.Rational(1,2)*eta-3)*chi_S**2)+vomega**5*((sp.Rational(-703,112)*eta**2+sp.Rational(8797,1008)*eta-sp.Rational(81,16))*(chi_A*deltainv)+(sp.Rational(613,1008)*eta**2+sp.Rational(1709,1008)*eta-sp.Rational(81,16))*chi_S+(sp.Rational(3,4)-3*eta)*(chi_A**3*deltainv)+(sp.Rational(9,4)-6*eta)*(chi_A*chi_S**2*deltainv)+(sp.Rational(9,4)-3*eta)*(chi_A**2*chi_S)+sp.Rational(3,4)*chi_S**3)+vomega**6*((sp.Rational(5,7)*eta**2-sp.Rational(9287,1008)*eta+sp.Rational(4163,252))*chi_A**2+(sp.Rational(139,72)*eta**2-sp.Rational(2633,1008)*eta+sp.Rational(4163,252))*chi_S**2+(sp.Rational(9487,504)*eta**2-sp.Rational(1636,21)*eta+sp.Rational(4163,126))*(chi_A*chi_S*deltainv)))
        rhoNS81=1+((20022-126451*eta+236922*eta**2-138430*eta**3+21640*eta**4)/(18240*(-1+6*eta-10*eta**2+4*eta**3)))*vomega**2
        rho82=1+((2462-17598*eta+37119*eta**2-22845*eta**3+3063*eta**4)/(2736*(-1+7*eta-14*eta**2+7*eta**3)))*vomega**2
        rhoNS83=1+((20598-131059*eta+249018*eta**2-149950*eta**3+24520*eta**4)/(18240*(-1+6*eta-10*eta**2+4*eta**3)))*vomega**2
        rho84=1+((2666-19434*eta+42627*eta**2-28965*eta**3+4899*eta**4)/(2736*(-1+7*eta-14*eta**2+7*eta**3)))*vomega**2
        rhoNS85=1+((4350-28055*eta+54642*eta**2-34598*eta**3+6056*eta**4)/(3648*(-1+6*eta-10*eta**2+4*eta**3)))*vomega**2
        rho86=1+((1002-7498*eta+17269*eta**2-13055*eta**3+2653*eta**4)/(912*(-1+7*eta-14*eta**2+7*eta**3)))*vomega**2
        rhoNS87=1+((23478-154099*eta+309498*eta**2-207550*eta**3+38920*eta**4)/(18240*(-1+6*eta-10*eta**2+4*eta**3)))*vomega**2
        rho88=1+((3482-26778*eta+64659*eta**2-53445*eta**3+12243*eta**4)/(2736*(-1+7*eta-14*eta**2+7*eta**3)))*vomega**2
        rhoNS71=1+((228*eta**3-2083*eta**2+2518*eta-618)/(714*(3*eta**2-4*eta+1)))*vomega**2
        rho72=1+((32760*eta**4-190239*eta**3+273924*eta**2-123489*eta+16832)/(14994*(7*eta**3-14*eta**2+7*eta-1)))*vomega**2
        rhoNS73=1+((420*eta**3-2563*eta**2+2806*eta-666)/(714*(3*eta**2-4*eta+1)))*vomega**2
        rho74=1+((41076*eta**4-217959*eta**3+298872*eta**2-131805*eta+17756)/(14994*(7*eta**3-14*eta**2+7*eta-1)))*vomega**2
        rhoNS75=1+((804*eta**3-3523*eta**2+3382*eta-762)/(714*(3*eta**2-4*eta+1)))*vomega**2
        rho76=1+((6104*eta**4-29351*eta**3+37828*eta**2-16185*eta+2144)/(1666*(7*eta**3-14*eta**2+7*eta-1)))*vomega**2
        rhoNS77=1+((1380*eta**3-4963*eta**2+4246*eta-906)/(714*(3*eta**2-4*eta+1)))*vomega**2
        rhoNS61=1+((124*eta**3-670*eta**2+694*eta-161)/(144*(3*eta**2-4*eta+1)))*vomega**2
        rho62=1+((49*eta**3-413*eta**2+378*eta-74)/(84*(5*eta**2-5*eta+1)))*vomega**2-sp.Rational(817991,3298680)*vomega**4
        rhoNS63=1+((156*eta**3-750*eta**2+742*eta-169)/(144*(3*eta**2-4*eta+1)))*vomega**2
        rho64=1+vomega**2*((133*eta**3-581*eta**2+462*eta-86)/(84*(5*eta**2-5*eta+1)))-sp.Rational(476887,659736)*vomega**4
        rhoNS65=1+((220*eta**3-910*eta**2+838*eta-185)/(144*(3*eta**2-4*eta+1)))*vomega**2
        rho66=1+((273*eta**3-861*eta**2+602*eta-106)/(84*(5*eta**2-5*eta+1)))*vomega**2-sp.Rational(1025435,659736)*vomega**4
        rhoNS51=1+((8*eta**2-626*eta+319)/(390*(2*eta-1)))*vomega**2-sp.Rational(31877,304200)*vomega**4
        rho52=1+((21980*eta**3-104930*eta**2+84679*eta-15828)/(13650*(5*eta**2-5*eta+1)))*vomega**2-sp.Rational(7187914,15526875)*vomega**4
        rhoNS53=1+((176*eta**2-850*eta+375)/(390*(2*eta-1)))*vomega**2-sp.Rational(410833,709800)*vomega**4
        rho54=1+((33320*eta**3-127610*eta**2+96019*eta-17448)/(13650*(5*eta**2-5*eta+1)))*vomega**2-sp.Rational(16213384,15526875)*vomega**4
        rhoNS55=1+vomega**2*(sp.Rational(487,390)/(-1+2*eta)-sp.Rational(649,195)*eta/(-1+2*eta)+sp.Rational(256,195)*eta**2/(-1+2*eta))-sp.Rational(3353747,2129400)*vomega**4+(sp.Rational(190606537999247,11957879934000)-sp.Rational(1546,429)*eulerlog5)*vomega**6+(-sp.Rational(1213641959949291437,118143853747920000)+sp.Rational(376451,83655)*eulerlog5)*vomega**8+(-sp.Rational(150082616449726042201261,4837990810977324000000)+sp.Rational(2592446431,456756300)*eulerlog5)*vomega**10+eta*(-2.61*vomega**4+1.25*vomega**6+-35.7*vomega**8)
        rhoNS41=1+((288*eta**2-1385*eta+602)/(528*(2*eta-1)))*vomega**2-(sp.Rational(7775491,21141120))*vomega**4+(sp.Rational(1227423222031,1758095539200)-sp.Rational(1571,6930)*eulerlog1)*vomega**6
        rho42=(1+vomega**2*((1146.0-3530.0*eta+285.0*eta**2)/(1320.0*(-1+3*eta)))+vomega**3*((chi_A*(10.0-21.0*eta)*delta+chi_S*(10.0-59.0*eta+78.0*eta**2))/(15.0*(-1+3*eta)))+vomega**4*((-114859044.0+295834536.0*eta+1204388696.0*eta**2-3047981160.0*eta**3-379526805.0*eta**4)/(317116800.0*(-1+3*eta)**2))+vomega**6*(848238724511.0/219761942400.-(3142.0/3465.0)*eulerlog2))
        rhoNS43=1+vomega**2/(1-2*eta)*(-sp.Rational(10,11)*eta**2+sp.Rational(547,176)*eta-sp.Rational(111,88))-sp.Rational(6894273,7047040)*vomega**4+vomega**6*(sp.Rational(1664224207351,195343948800)-sp.Rational(1571,770)*eulerlog3)+vomega**8*(-sp.Rational(2465107182496333,460490801971200)+sp.Rational(174381,67760)*eulerlog3)+eta*(-0.654*vomega**4+-3.69*vomega**6+18.5*vomega**8)
        rho44=(1+vomega**2*((1614.0-5870.0*eta+2625.0*eta**2)/(1320.0*(-1+3*eta)))+vomega**3*(chi_A*(10.0-39.0*eta)*delta+chi_S*(10.0-41.0*eta+42.0*eta**2))/(15.0*(-1+3*eta))+vomega**4*((-511573572.0+2338945704.0*eta-313857376.0*eta**2-6733146000.0*eta**3+1252563795.0*eta**4)/(317116800.0*(-1+3*eta)**2)+chi_S**2/2.0+delta*chi_S*chi_A+delta**2*chi_A**2/2.0)+vomega**5*(chi_A*delta*(-8280.0+42716.0*eta-57990.0*eta**2+8955*eta**3)/(6600.0*(-1+3*eta)**2)+chi_S*(-8280.0+66284.0*eta-176418.0*eta**2+128085.0*eta**3+88650*eta**2*eta**2)/(6600.0*(-1+3*eta)**2))+vomega**6*(sp.Rational(16600939332793,1098809712000)-sp.Rational(12568.0,3465.0)*eulerlog4)+vomega**8*(-sp.Rational(172066910136202271,19426955708160000)+sp.Rational(845198.0,190575.0)*eulerlog4)+vomega**10*(-sp.Rational(17154485653213713419357,568432724020761600000)+sp.Rational(22324502267,3815311500)*eulerlog4)+eta*(-3.56*vomega**6+15.6*vomega**8+-216*vomega**10))
        rhoNS31=(1-vomega**2*(sp.Rational(2,9)*eta+sp.Rational(13,18))+vomega**4*(-sp.Rational(829,1782)*eta**2-sp.Rational(1685,1782)*eta+sp.Rational(101,7128))+vomega**6*(sp.Rational(11706720301,6129723600)-sp.Rational(26,63)*eulerlog1)+vomega**8*(sp.Rational(169,567)*eulerlog1+sp.Rational(2606097992581,4854741091200)))
        rho32=1+vomega*((4*eta*chi_S)/(3*(1-3*eta)))+vomega**2*((sp.Rational(-32,27)*eta**2+sp.Rational(223,54)*eta-sp.Rational(164,135))/(1-3*eta)-(16*eta**2*chi_S**2)/(9*(1-3*eta)**2))+vomega**3*((sp.Rational(13,9)*eta+sp.Rational(2,9))*(delta*chi_A)/(1-3*eta)+(sp.Rational(607,81)*eta**3+sp.Rational(503,81)*eta**2-sp.Rational(1478,405)*eta+sp.Rational(2,9))*chi_S/(1-3*eta)**2+(320*eta**3*chi_S**3)/(81*(1-3*eta)**3))+vomega**4*((sp.Rational(77141,40095)*eta**4-sp.Rational(508474,40095)*eta**3-sp.Rational(945121,320760)*eta**2+sp.Rational(1610009,320760)*eta-sp.Rational(180566,200475))/(1-3*eta)**2+(4*eta**2-3*eta+sp.Rational(1,3))*(chi_A**2)/(1-3*eta)+(sp.Rational(-50,27)*eta**2-sp.Rational(88,27)*eta+sp.Rational(2,3))*(delta*chi_A*chi_S)/(1-3*eta)**2+(sp.Rational(-2452,243)*eta**4-sp.Rational(1997,243)*eta**3+sp.Rational(1435,243)*eta**2-sp.Rational(43,27)*eta+sp.Rational(1,3))*(chi_S**2)/((1-3*eta)**3))+vomega**5*((sp.Rational(-1184225,96228)*eta**5-sp.Rational(40204523,962280)*eta**4+sp.Rational(101706029,962280)*eta**3-sp.Rational(14103833,192456)*eta**2+sp.Rational(20471053,962280)*eta-sp.Rational(2788,1215))*chi_S/(1-3*eta)**3+(sp.Rational(608,81)*eta**3+sp.Rational(736,81)*eta**2-sp.Rational(16,9)*eta)*(delta*chi_A*chi_S**2)/(1-3*eta)**3+(sp.Rational(889673,106920)*eta**3-sp.Rational(75737,5346)*eta**2+sp.Rational(376177,35640)*eta-sp.Rational(2788,1215))*(delta*chi_A)/(1-3*eta)**2+(sp.Rational(96176,2187)*eta**5+sp.Rational(43528,2187)*eta**4-sp.Rational(40232,2187)*eta**3+sp.Rational(376,81)*eta**2-sp.Rational(8,9)*eta)*(chi_S**3)/(1-3*eta)**4+(sp.Rational(-32,3)*eta**3+8*eta**2-sp.Rational(8,9)*eta)*(chi_A**2*chi_S)/((1-3*eta)**2))+vomega**6*(sp.Rational(5849948554,940355325)-sp.Rational(104,63)*eulerlog2)+vomega**8*(sp.Rational(17056,8505)*eulerlog2-sp.Rational(10607269449358,3072140846775))+vomega**10*(-sp.Rational(1312549797426453052,176264081083715625)+sp.Rational(18778864,12629925)*eulerlog2)+eta*(+.333*vomega**6-6.5*vomega**8+98*vomega**10)
        rhoNS33=1+vomega**2*(sp.Rational(2,3)*eta-sp.Rational(7,6))+vomega**4*(-sp.Rational(6719,3960)-sp.Rational(1861,990)*eta+sp.Rational(149,330)*eta**2)+vomega**6*(sp.Rational(3203101567,227026800)+(-sp.Rational(129509,25740)+sp.Rational(41,192)*sp.pi**2)*eta-sp.Rational(274621,154440)*eta**2+sp.Rational(12011,46332)*eta**3-sp.Rational(26,7)*eulerlog3)+vomega**8*(-sp.Rational(57566572157,8562153600)+sp.Rational(13,3)*eulerlog3)+vomega**10*(-sp.Rational(903823148417327,30566888352000)+sp.Rational(87347,13860)*eulerlog3)+eta*(12*vomega**8+-215*vomega**10)
        rhoNS21=1+vomega**2*(sp.Rational(23,84)*eta-sp.Rational(59,56))+vomega**4*(sp.Rational(617,4704)*eta**2-sp.Rational(10993,14112)*eta-sp.Rational(47009,56448))+vomega**6*(sp.Rational(7613184941,2607897600)-sp.Rational(107,105)*eulerlog1)+vomega**8*(-sp.Rational(1168617463883,911303737344)+sp.Rational(6313,5880)*eulerlog1)+vomega**10*(-sp.Rational(63735873771463,16569158860800)+sp.Rational(5029963,5927040)*eulerlog1)+eta*(1.65*vomega**6+26.5*vomega**8+80*vomega**10)
        rho22=1+vomega**2*(sp.Rational(55,84)*eta-sp.Rational(43,42))+vomega**3*((-2.0*(chi_S+chi_A*delta-chi_S*eta))/3.0)+vomega**4*(sp.Rational(19583,42336)*eta**2-sp.Rational(33025,21168)*eta-sp.Rational(20555,10584)+(sp.Rational(1,2)-2*eta)*chi_A**2+delta*chi_A*chi_S+sp.Rational(1,2)*chi_S**2)+vomega**5*(delta*(-sp.Rational(19,42)*eta-sp.Rational(34,21))*chi_A+(sp.Rational(209,126)*eta**2+sp.Rational(49,18)*eta-sp.Rational(34,21))*chi_S)+vomega**6*(sp.Rational(10620745,39118464)*eta**3-sp.Rational(6292061,3259872)*eta**2+sp.Rational(41,192)*sp.pi**2*eta-sp.Rational(48993925,9779616)*eta-sp.Rational(428,105)*eulerlog2+sp.Rational(1556919113,122245200)+delta*(sp.Rational(89,126)-sp.Rational(781,252)*eta)*chi_A*chi_S+(-sp.Rational(27,14)*eta**2-sp.Rational(457,504)*eta+sp.Rational(89,252))*chi_A**2+(sp.Rational(10,9)*eta**2-sp.Rational(1817,504)*eta+sp.Rational(89,252))*chi_S**2)+vomega**7*(delta*(sp.Rational(97865,63504)*eta**2+sp.Rational(50140,3969)*eta+sp.Rational(18733,15876))*chi_A+(sp.Rational(50803,63504)*eta**3-sp.Rational(245717,63504)*eta**2+sp.Rational(74749,5292)*eta+sp.Rational(18733,15876))*chi_S+delta*chi_A**3*(sp.Rational(1,3)-sp.Rational(4,3)*eta)+delta*(2*eta+1)*chi_A*chi_S**2+(sp.Rational(-4,1)*eta**2-sp.Rational(3,1)*eta+sp.Rational(1,1))*chi_A**2*chi_S+(eta+sp.Rational(1,3))*chi_S**3)+vomega**8*(sp.Rational(9202,2205)*eulerlog2-sp.Rational(387216563023,160190110080))+vomega**10*(sp.Rational(439877,55566)*eulerlog2-sp.Rational(16094530514677,533967033600))+eta*(21.2*vomega**8+-411*vomega**10)
        f81amp = (rhoNS81)**8
        f82amp = (rho82)**8
        f83amp = (rhoNS83)**8
        f84amp = (rho84)**8
        f85amp = (rhoNS85)**8
        f86amp = (rho86)**8
        f87amp = (rhoNS87)**8
        f88amp = (rho88)**8
        f71amp = (rhoNS71)**7
        f72amp = (rho72)**7
        f73amp = (rhoNS73)**7
        f74amp = (rho74)**7
        f75amp = (rhoNS75)**7
        f76amp = (rho76)**7
        f77amp = (rhoNS77)**7
        f61amp = (rhoNS61)**6
        f62amp = (rho62)**6
        f63amp = (rhoNS63)**6
        f64amp = (rho64)**6
        f65amp = (rhoNS65)**6
        f66amp = (rho66)**6
        f51amp = ((rhoNS51)**5)
        f52amp = (rho52)**5
        f53amp = ((rhoNS53)**5)
        f54amp = (rho54)**5
        f55 = ((rhoNS55)**5 + fspin55)*noneqcond + fspin55_limit*eqcond
        f55amp = f55
        f41amp = ((rhoNS41)**4 + fspin41)*noneqcond + fspin41_limit*eqcond
        f42amp = (rho42)**4
        f43 = ((rhoNS43)**4 + fspin43)*noneqcond + fspin43_limit*eqcond
        f43amp = f43
        f44 = (rho44)**4
        f44amp = f44
        f31amp = ((rhoNS31)**3 + fspin31)*noneqcond + fspin31_limit*eqcond
        f32 = (rho32)**3
        f32amp = f32
        f33 = ((rhoNS33)**3 + fspin33)*noneqcond + fspin33_limit*eqcond
        f33amp = ((rhoNS33)**3 + fspin33amp)*noneqcond + fspin33amp_limit*eqcond
        f21 = ((rhoNS21)**2 + fspin21)*noneqcond + fspin21_limit*eqcond
        f21amp = f21
        f22 = (rho22)**2
        f22amp = f22
        Y82amp = 0.32254835519288305
        Y84amp = 0.3382915688890245
        Y86amp = 0.3764161087284946
        Y88amp = 0.5154289843972844
        Y71amp = 0.31937046138540076
        Y73amp = 0.331899519333737
        Y75amp = 0.3669287245764378
        Y77amp = 0.5000395635705508
        Y62amp = 0.32569524293385776
        Y64amp = 0.3567812628539981
        Y66amp = 0.48308411358006625
        Y51amp = 0.32028164857621516
        Y53amp = 0.34594371914684025
        Y55amp = 0.46413220344085826
        Y55 = 0.46413220344085826*sp.exp(-5*sp.I*phi)
        Y42amp = 0.33452327177864466
        Y44amp = 0.4425326924449826
        Y44 = 0.4425326924449826*sp.exp(-4*sp.I*phi)
        Y31amp = 0.3231801841141506
        Y33amp = 0.4172238236327842
        Y33 = 0.4172238236327842*sp.exp(-3*sp.I*phi)
        Y11amp = 0.3454941494713355
        Y11 = 0.3454941494713355*sp.exp(-sp.I*phi)
        Y22amp = 0.3862742020231896
        Y22 = 0.3862742020231896*sp.exp(-2*sp.I*phi)
        c9 = ((1 - delta)/2)**8 + ((-1)**9)*((1 + delta)/2)**8
        c8 = ((1 - delta)/2)**7 + ((-1)**8)*((1 + delta)/2)**7
        c7 = ((1 - delta)/2)**6 + ((-1)**7)*((1 + delta)/2)**6
        c6 = ((1 - delta)/2)**5 + ((-1)**6)*((1 + delta)/2)**5
        c5 = (((1 - delta)/2)**4 + ((-1)**5)*((1 + delta)/2)**4)*noneqcond + (-sp.Rational(1,2))*eqcond
        c4 = ((1 - delta)/2)**3 + ((-1)**4)*((1 + delta)/2)**3
        c3 = (((1 - delta)/2)**2 + ((-1)**3)*((1 + delta)/2)**2)*noneqcond + (-1)*eqcond
        c2 = ((1 - delta)/2)**1 + ((-1)**2)*((1 + delta)/2)**1
        N81 = (-16*sp.I*sp.pi*(sp.I*1)**8/sp.factorial2(2*8+1))*sp.sqrt( ( (2*8+1)*(8+2)*(8**2 - 1**2) ) / ( (2*8-1)*(8+1)*(8)*(8-1) ) )
        N82 = (8*sp.I*sp.pi*(sp.I*2)**8/sp.factorial2(2*8 + 1))*sp.sqrt( ( (8+1)*(8+2) ) / ( (8)*(8-1) ) )
        N83 = (-16*sp.I*sp.pi*(sp.I*3)**8/sp.factorial2(2*8+1))*sp.sqrt( ( (2*8+1)*(8+2)*(8**2 - 3**2) ) / ( (2*8-1)*(8+1)*(8)*(8-1) ) )
        N84 = (8*sp.I*sp.pi*(sp.I*4)**8/sp.factorial2(2*8 + 1))*sp.sqrt( ( (8+1)*(8+2) ) / ( (8)*(8-1) ) )
        N85 = (-16*sp.I*sp.pi*(sp.I*5)**8/sp.factorial2(2*8+1))*sp.sqrt( ( (2*8+1)*(8+2)*(8**2 - 5**2) ) / ( (2*8-1)*(8+1)*(8)*(8-1) ) )
        N86 = (8*sp.I*sp.pi*(sp.I*6)**8/sp.factorial2(2*8 + 1))*sp.sqrt( ( (8+1)*(8+2) ) / ( (8)*(8-1) ) )
        N87 = (-16*sp.I*sp.pi*(sp.I*7)**8/sp.factorial2(2*8+1))*sp.sqrt( ( (2*8+1)*(8+2)*(8**2 - 7**2) ) / ( (2*8-1)*(8+1)*(8)*(8-1) ) )
        N88 = (8*sp.I*sp.pi*(sp.I*8)**8/sp.factorial2(2*8 + 1))*sp.sqrt( ( (8+1)*(8+2) ) / ( (8)*(8-1) ) )
        N71 = (8*sp.I*sp.pi*(sp.I*1)**7/sp.factorial2(2*7 + 1))*sp.sqrt( ( (7+1)*(7+2) ) / ( (7)*(7-1) ) )
        N72 = (-16*sp.I*sp.pi*(sp.I*2)**7/sp.factorial2(2*7+1))*sp.sqrt( ( (2*7+1)*(7+2)*(7**2 - 2**2) ) / ( (2*7-1)*(7+1)*(7)*(7-1) ) )
        N73 = (8*sp.I*sp.pi*(sp.I*3)**7/sp.factorial2(2*7 + 1))*sp.sqrt( ( (7+1)*(7+2) ) / ( (7)*(7-1) ) )
        N74 = (-16*sp.I*sp.pi*(sp.I*4)**7/sp.factorial2(2*7+1))*sp.sqrt( ( (2*7+1)*(7+2)*(7**2 - 4**2) ) / ( (2*7-1)*(7+1)*(7)*(7-1) ) )
        N75 = (8*sp.I*sp.pi*(sp.I*5)**7/sp.factorial2(2*7 + 1))*sp.sqrt( ( (7+1)*(7+2) ) / ( (7)*(7-1) ) )
        N76 = (-16*sp.I*sp.pi*(sp.I*6)**7/sp.factorial2(2*7+1))*sp.sqrt( ( (2*7+1)*(7+2)*(7**2 - 6**2) ) / ( (2*7-1)*(7+1)*(7)*(7-1) ) )
        N77 = (8*sp.I*sp.pi*(sp.I*7)**7/sp.factorial2(2*7 + 1))*sp.sqrt( ( (7+1)*(7+2) ) / ( (7)*(7-1) ) )
        N61 = (-16*sp.I*sp.pi*(sp.I*1)**6/sp.factorial2(2*6+1))*sp.sqrt( ( (2*6+1)*(6+2)*(6**2 - 1**2) ) / ( (2*6-1)*(6+1)*(6)*(6-1) ) )
        N62 = (8*sp.I*sp.pi*(sp.I*2)**6/sp.factorial2(2*6 + 1))*sp.sqrt( ( (6+1)*(6+2) ) / ( (6)*(6-1) ) )
        N63 = (-16*sp.I*sp.pi*(sp.I*3)**6/sp.factorial2(2*6+1))*sp.sqrt( ( (2*6+1)*(6+2)*(6**2 - 3**2) ) / ( (2*6-1)*(6+1)*(6)*(6-1) ) )
        N64 = (8*sp.I*sp.pi*(sp.I*4)**6/sp.factorial2(2*6 + 1))*sp.sqrt( ( (6+1)*(6+2) ) / ( (6)*(6-1) ) )
        N65 = (-16*sp.I*sp.pi*(sp.I*5)**6/sp.factorial2(2*6+1))*sp.sqrt( ( (2*6+1)*(6+2)*(6**2 - 5**2) ) / ( (2*6-1)*(6+1)*(6)*(6-1) ) )
        N66 = (8*sp.I*sp.pi*(sp.I*6)**6/sp.factorial2(2*6 + 1))*sp.sqrt( ( (6+1)*(6+2) ) / ( (6)*(6-1) ) )
        N51 = (8*sp.I*sp.pi*(sp.I*1)**5/sp.factorial2(2*5 + 1))*sp.sqrt( ( (5+1)*(5+2) ) / ( (5)*(5-1) ) )
        N52 = (-16*sp.I*sp.pi*(sp.I*2)**5/sp.factorial2(2*5+1))*sp.sqrt( ( (2*5+1)*(5+2)*(5**2 - 2**2) ) / ( (2*5-1)*(5+1)*(5)*(5-1) ) )
        N53 = (8*sp.I*sp.pi*(sp.I*3)**5/sp.factorial2(2*5 + 1))*sp.sqrt( ( (5+1)*(5+2) ) / ( (5)*(5-1) ) )
        N54 = (-16*sp.I*sp.pi*(sp.I*4)**5/sp.factorial2(2*5+1))*sp.sqrt( ( (2*5+1)*(5+2)*(5**2 - 4**2) ) / ( (2*5-1)*(5+1)*(5)*(5-1) ) )
        N55 = (8*sp.I*sp.pi*(sp.I*5)**5/sp.factorial2(2*5 + 1))*sp.sqrt( ( (5+1)*(5+2) ) / ( (5)*(5-1) ) )
        N41 = (-16*sp.I*sp.pi*(sp.I*1)**4/sp.factorial2(2*4+1))*sp.sqrt( ( (2*4+1)*(4+2)*(4**2 - 1**2) ) / ( (2*4-1)*(4+1)*(4)*(4-1) ) )
        N42 = (8*sp.I*sp.pi*(sp.I*2)**4/sp.factorial2(2*4 + 1))*sp.sqrt( ( (4+1)*(4+2) ) / ( (4)*(4-1) ) )
        N43 = (-16*sp.I*sp.pi*(sp.I*3)**4/sp.factorial2(2*4+1))*sp.sqrt( ( (2*4+1)*(4+2)*(4**2 - 3**2) ) / ( (2*4-1)*(4+1)*(4)*(4-1) ) )
        N44 = (8*sp.I*sp.pi*(sp.I*4)**4/sp.factorial2(2*4 + 1))*sp.sqrt( ( (4+1)*(4+2) ) / ( (4)*(4-1) ) )
        N31 = (8*sp.I*sp.pi*(sp.I*1)**3/sp.factorial2(2*3 + 1))*sp.sqrt( ( (3+1)*(3+2) ) / ( (3)*(3-1) ) )
        N32 = (-16*sp.I*sp.pi*(sp.I*2)**3/sp.factorial2(2*3+1))*sp.sqrt( ( (2*3+1)*(3+2)*(3**2 - 2**2) ) / ( (2*3-1)*(3+1)*(3)*(3-1) ) )
        N33 = (8*sp.I*sp.pi*(sp.I*3)**3/sp.factorial2(2*3 + 1))*sp.sqrt( ( (3+1)*(3+2) ) / ( (3)*(3-1) ) )
        N21 = (-16*sp.I*sp.pi*(sp.I*1)**2/sp.factorial2(2*2+1))*sp.sqrt( ( (2*2+1)*(2+2)*(2**2 - 1**2) ) / ( (2*2-1)*(2+1)*(2)*(2-1) ) )
        N22 = (8*sp.I*sp.pi*(sp.I*2)**2/sp.factorial2(2*2 + 1))*sp.sqrt( ( (2+1)*(2+2) ) / ( (2)*(2-1) ) )
        N81amp = sp.Abs(N81)
        N82amp = sp.Abs(N82)
        N83amp = sp.Abs(N83)
        N84amp = sp.Abs(N84)
        N85amp = sp.Abs(N85)
        N86amp = sp.Abs(N86)
        N87amp = sp.Abs(N87)
        N88amp = sp.Abs(N88)
        N71amp = sp.Abs(N71)
        N72amp = sp.Abs(N72)
        N73amp = sp.Abs(N73)
        N74amp = sp.Abs(N74)
        N75amp = sp.Abs(N75)
        N76amp = sp.Abs(N76)
        N77amp = sp.Abs(N77)
        N61amp = sp.Abs(N61)
        N62amp = sp.Abs(N62)
        N63amp = sp.Abs(N63)
        N64amp = sp.Abs(N64)
        N65amp = sp.Abs(N65)
        N66amp = sp.Abs(N66)
        N51amp = sp.Abs(N51)
        N52amp = sp.Abs(N52)
        N53amp = sp.Abs(N53)
        N54amp = sp.Abs(N54)
        N55amp = sp.Abs(N55)
        N41amp = sp.Abs(N41)
        N42amp = sp.Abs(N42)
        N43amp = sp.Abs(N43)
        N44amp = sp.Abs(N44)
        N31amp = sp.Abs(N31)
        N32amp = sp.Abs(N32)
        N33amp = sp.Abs(N33)
        N21amp = sp.Abs(N21)
        N22amp = sp.Abs(N22)
        r0 = 2/sp.sqrt(sp.exp(1))
        b8 = -2*khat8
        b7 = -2*khat7
        b6 = -2*khat6
        b5 = -2*khat5
        b4 = -2*khat4
        b3 = -2*khat3
        b2 = -2*khat2
        b1 = -2*khat1
        T81prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)*(4**2 + b1**2)*(5**2 + b1**2)*(6**2 + b1**2)*(7**2 + b1**2)*(8**2 + b1**2)
        T82prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)*(4**2 + b2**2)*(5**2 + b2**2)*(6**2 + b2**2)*(7**2 + b2**2)*(8**2 + b2**2)
        T83prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)*(4**2 + b3**2)*(5**2 + b3**2)*(6**2 + b3**2)*(7**2 + b3**2)*(8**2 + b3**2)
        T84prodfac = (1**2 + b4**2)*(2**2 + b4**2)*(3**2 + b4**2)*(4**2 + b4**2)*(5**2 + b4**2)*(6**2 + b4**2)*(7**2 + b4**2)*(8**2 + b4**2)
        T85prodfac = (1**2 + b5**2)*(2**2 + b5**2)*(3**2 + b5**2)*(4**2 + b5**2)*(5**2 + b5**2)*(6**2 + b5**2)*(7**2 + b5**2)*(8**2 + b5**2)
        T86prodfac = (1**2 + b6**2)*(2**2 + b6**2)*(3**2 + b6**2)*(4**2 + b6**2)*(5**2 + b6**2)*(6**2 + b6**2)*(7**2 + b6**2)*(8**2 + b6**2)
        T87prodfac = (1**2 + b7**2)*(2**2 + b7**2)*(3**2 + b7**2)*(4**2 + b7**2)*(5**2 + b7**2)*(6**2 + b7**2)*(7**2 + b7**2)*(8**2 + b7**2)
        T88prodfac = (1**2 + b8**2)*(2**2 + b8**2)*(3**2 + b8**2)*(4**2 + b8**2)*(5**2 + b8**2)*(6**2 + b8**2)*(7**2 + b8**2)*(8**2 + b8**2)
        T71prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)*(4**2 + b1**2)*(5**2 + b1**2)*(6**2 + b1**2)*(7**2 + b1**2)
        T72prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)*(4**2 + b2**2)*(5**2 + b2**2)*(6**2 + b2**2)*(7**2 + b2**2)
        T73prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)*(4**2 + b3**2)*(5**2 + b3**2)*(6**2 + b3**2)*(7**2 + b3**2)
        T74prodfac = (1**2 + b4**2)*(2**2 + b4**2)*(3**2 + b4**2)*(4**2 + b4**2)*(5**2 + b4**2)*(6**2 + b4**2)*(7**2 + b4**2)
        T75prodfac = (1**2 + b5**2)*(2**2 + b5**2)*(3**2 + b5**2)*(4**2 + b5**2)*(5**2 + b5**2)*(6**2 + b5**2)*(7**2 + b5**2)
        T76prodfac = (1**2 + b6**2)*(2**2 + b6**2)*(3**2 + b6**2)*(4**2 + b6**2)*(5**2 + b6**2)*(6**2 + b6**2)*(7**2 + b6**2)
        T77prodfac = (1**2 + b7**2)*(2**2 + b7**2)*(3**2 + b7**2)*(4**2 + b7**2)*(5**2 + b7**2)*(6**2 + b7**2)*(7**2 + b7**2)
        T61prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)*(4**2 + b1**2)*(5**2 + b1**2)*(6**2 + b1**2)
        T62prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)*(4**2 + b2**2)*(5**2 + b2**2)*(6**2 + b2**2)
        T63prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)*(4**2 + b3**2)*(5**2 + b3**2)*(6**2 + b3**2)
        T64prodfac = (1**2 + b4**2)*(2**2 + b4**2)*(3**2 + b4**2)*(4**2 + b4**2)*(5**2 + b4**2)*(6**2 + b4**2)
        T65prodfac = (1**2 + b5**2)*(2**2 + b5**2)*(3**2 + b5**2)*(4**2 + b5**2)*(5**2 + b5**2)*(6**2 + b5**2)
        T66prodfac = (1**2 + b6**2)*(2**2 + b6**2)*(3**2 + b6**2)*(4**2 + b6**2)*(5**2 + b6**2)*(6**2 + b6**2)
        T51prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)*(4**2 + b1**2)*(5**2 + b1**2)
        T52prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)*(4**2 + b2**2)*(5**2 + b2**2)
        T53prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)*(4**2 + b3**2)*(5**2 + b3**2)
        T54prodfac = (1**2 + b4**2)*(2**2 + b4**2)*(3**2 + b4**2)*(4**2 + b4**2)*(5**2 + b4**2)
        T55prodfac = (1**2 + b5**2)*(2**2 + b5**2)*(3**2 + b5**2)*(4**2 + b5**2)*(5**2 + b5**2)
        T41prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)*(4**2 + b1**2)
        T42prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)*(4**2 + b2**2)
        T43prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)*(4**2 + b3**2)
        T44prodfac = (1**2 + b4**2)*(2**2 + b4**2)*(3**2 + b4**2)*(4**2 + b4**2)
        T31prodfac = (1**2 + b1**2)*(2**2 + b1**2)*(3**2 + b1**2)
        T32prodfac = (1**2 + b2**2)*(2**2 + b2**2)*(3**2 + b2**2)
        T33prodfac = (1**2 + b3**2)*(2**2 + b3**2)*(3**2 + b3**2)
        T21prodfac = (1**2 + b1**2)*(2**2 + b1**2)
        T22prodfac = (1**2 + b2**2)*(2**2 + b2**2)
        T8prefac = 2*sp.pi*b8/(sp.exp(sp.pi*b8) - sp.exp(-sp.pi*b8))
        T7prefac = 2*sp.pi*b7/(sp.exp(sp.pi*b7) - sp.exp(-sp.pi*b7))
        T6prefac = 2*sp.pi*b6/(sp.exp(sp.pi*b6) - sp.exp(-sp.pi*b6))
        T5prefac = 2*sp.pi*b5/(sp.exp(sp.pi*b5) - sp.exp(-sp.pi*b5))
        T4prefac = 2*sp.pi*b4/(sp.exp(sp.pi*b4) - sp.exp(-sp.pi*b4))
        T3prefac = 2*sp.pi*b3/(sp.exp(sp.pi*b3) - sp.exp(-sp.pi*b3))
        T2prefac = 2*sp.pi*b2/(sp.exp(sp.pi*b2) - sp.exp(-sp.pi*b2))
        T1prefac = 2*sp.pi*b1/(sp.exp(sp.pi*b1) - sp.exp(-sp.pi*b1))
        gamma_amp_81 = sp.sqrt(T1prefac*T81prodfac)
        gamma_amp_82 = sp.sqrt(T2prefac*T82prodfac)
        gamma_amp_83 = sp.sqrt(T3prefac*T83prodfac)
        gamma_amp_84 = sp.sqrt(T4prefac*T84prodfac)
        gamma_amp_85 = sp.sqrt(T5prefac*T85prodfac)
        gamma_amp_86 = sp.sqrt(T6prefac*T86prodfac)
        gamma_amp_87 = sp.sqrt(T7prefac*T87prodfac)
        gamma_amp_88 = sp.sqrt(T8prefac*T88prodfac)
        gamma_amp_71 = sp.sqrt(T1prefac*T71prodfac)
        gamma_amp_72 = sp.sqrt(T2prefac*T72prodfac)
        gamma_amp_73 = sp.sqrt(T3prefac*T73prodfac)
        gamma_amp_74 = sp.sqrt(T4prefac*T74prodfac)
        gamma_amp_75 = sp.sqrt(T5prefac*T75prodfac)
        gamma_amp_76 = sp.sqrt(T6prefac*T76prodfac)
        gamma_amp_77 = sp.sqrt(T7prefac*T77prodfac)
        gamma_amp_61 = sp.sqrt(T1prefac*T61prodfac)
        gamma_amp_62 = sp.sqrt(T2prefac*T62prodfac)
        gamma_amp_63 = sp.sqrt(T3prefac*T63prodfac)
        gamma_amp_64 = sp.sqrt(T4prefac*T64prodfac)
        gamma_amp_65 = sp.sqrt(T5prefac*T65prodfac)
        gamma_amp_66 = sp.sqrt(T6prefac*T66prodfac)
        gamma_amp_51 = sp.sqrt(T1prefac*T51prodfac)
        gamma_amp_52 = sp.sqrt(T2prefac*T52prodfac)
        gamma_amp_53 = sp.sqrt(T3prefac*T53prodfac)
        gamma_amp_54 = sp.sqrt(T4prefac*T54prodfac)
        gamma_amp_55 = sp.sqrt(T5prefac*T55prodfac)
        gamma_amp_41 = sp.sqrt(T1prefac*T41prodfac)
        gamma_amp_42 = sp.sqrt(T2prefac*T42prodfac)
        gamma_amp_43 = sp.sqrt(T3prefac*T43prodfac)
        gamma_amp_44 = sp.sqrt(T4prefac*T44prodfac)
        gamma_amp_31 = sp.sqrt(T1prefac*T31prodfac)
        gamma_amp_32 = sp.sqrt(T2prefac*T32prodfac)
        gamma_amp_33 = sp.sqrt(T3prefac*T33prodfac)
        gamma_amp_21 = sp.sqrt(T1prefac*T21prodfac)
        gamma_amp_22 = sp.sqrt(T2prefac*T22prodfac)
        T81amp = gamma_amp_81*sp.exp(sp.pi*khat1)/sp.factorial(8)
        T82amp = gamma_amp_82*sp.exp(sp.pi*khat2)/sp.factorial(8)
        T83amp = gamma_amp_83*sp.exp(sp.pi*khat3)/sp.factorial(8)
        T84amp = gamma_amp_84*sp.exp(sp.pi*khat4)/sp.factorial(8)
        T85amp = gamma_amp_85*sp.exp(sp.pi*khat5)/sp.factorial(8)
        T86amp = gamma_amp_86*sp.exp(sp.pi*khat6)/sp.factorial(8)
        T87amp = gamma_amp_87*sp.exp(sp.pi*khat7)/sp.factorial(8)
        T88amp = gamma_amp_88*sp.exp(sp.pi*khat8)/sp.factorial(8)
        T71amp = gamma_amp_71*sp.exp(sp.pi*khat1)/sp.factorial(7)
        T72amp = gamma_amp_72*sp.exp(sp.pi*khat2)/sp.factorial(7)
        T73amp = gamma_amp_73*sp.exp(sp.pi*khat3)/sp.factorial(7)
        T74amp = gamma_amp_74*sp.exp(sp.pi*khat4)/sp.factorial(7)
        T75amp = gamma_amp_75*sp.exp(sp.pi*khat5)/sp.factorial(7)
        T76amp = gamma_amp_76*sp.exp(sp.pi*khat6)/sp.factorial(7)
        T77amp = gamma_amp_77*sp.exp(sp.pi*khat7)/sp.factorial(7)
        T61amp = gamma_amp_61*sp.exp(sp.pi*khat1)/sp.factorial(6)
        T62amp = gamma_amp_62*sp.exp(sp.pi*khat2)/sp.factorial(6)
        T63amp = gamma_amp_63*sp.exp(sp.pi*khat3)/sp.factorial(6)
        T64amp = gamma_amp_64*sp.exp(sp.pi*khat4)/sp.factorial(6)
        T65amp = gamma_amp_65*sp.exp(sp.pi*khat5)/sp.factorial(6)
        T66amp = gamma_amp_66*sp.exp(sp.pi*khat6)/sp.factorial(6)
        T51amp = gamma_amp_51*sp.exp(sp.pi*khat1)/sp.factorial(5)
        T52amp = gamma_amp_52*sp.exp(sp.pi*khat2)/sp.factorial(5)
        T53amp = gamma_amp_53*sp.exp(sp.pi*khat3)/sp.factorial(5)
        T54amp = gamma_amp_54*sp.exp(sp.pi*khat4)/sp.factorial(5)
        T55amp = gamma_amp_55*sp.exp(sp.pi*khat5)/sp.factorial(5)
        T55 = sp.gamma(5 + 1 - 2*sp.I*khat5)*sp.exp(sp.pi*khat5)*(sp.exp(2*sp.I*khat5*sp.log(2*5*Omega*r0)))/sp.factorial(5)
        T41amp = gamma_amp_41*sp.exp(sp.pi*khat1)/sp.factorial(4)
        T42amp = gamma_amp_42*sp.exp(sp.pi*khat2)/sp.factorial(4)
        T43amp = gamma_amp_43*sp.exp(sp.pi*khat3)/sp.factorial(4)
        T43 = sp.gamma(4 + 1 - 2*sp.I*khat3)*sp.exp(sp.pi*khat3)*(sp.exp(2*sp.I*khat3*sp.log(2*3*Omega*r0)))/sp.factorial(4)
        T44amp = gamma_amp_44*sp.exp(sp.pi*khat4)/sp.factorial(4)
        T44 = sp.gamma(4 + 1 - 2*sp.I*khat4)*sp.exp(sp.pi*khat4)*(sp.exp(2*sp.I*khat4*sp.log(2*4*Omega*r0)))/sp.factorial(4)
        T31amp = gamma_amp_31*sp.exp(sp.pi*khat1)/sp.factorial(3)
        T32amp = gamma_amp_32*sp.exp(sp.pi*khat2)/sp.factorial(3)
        T32 = sp.gamma(3 + 1 - 2*sp.I*khat2)*sp.exp(sp.pi*khat2)*(sp.exp(2*sp.I*khat2*sp.log(2*2*Omega*r0)))/sp.factorial(3)
        T33amp = gamma_amp_33*sp.exp(sp.pi*khat3)/sp.factorial(3)
        T33 = sp.gamma(3 + 1 - 2*sp.I*khat3)*sp.exp(sp.pi*khat3)*(sp.exp(2*sp.I*khat3*sp.log(2*3*Omega*r0)))/sp.factorial(3)
        T21amp = gamma_amp_21*sp.exp(sp.pi*khat1)/sp.factorial(2)
        T21 = sp.gamma(2 + 1 - 2*sp.I*khat1)*sp.exp(sp.pi*khat1)*(sp.exp(2*sp.I*khat1*sp.log(2*1*Omega*r0)))/sp.factorial(2)
        T22amp = gamma_amp_22*sp.exp(sp.pi*khat2)/sp.factorial(2)
        T22 = sp.gamma(2 + 1 - 2*sp.I*khat2)*sp.exp(sp.pi*khat2)*(sp.exp(2*sp.I*khat2*sp.log(2*2*Omega*r0)))/sp.factorial(2)
        source_odd = vomega*pphi
        source_even = Heff
        hN81amp = eta*N81amp*sp.Abs(c9)*vphi**9*Y71amp
        hN82amp = eta*N82amp*sp.Abs(c8)*vphi**8*Y82amp
        hN83amp = eta*N83amp*sp.Abs(c9)*vphi**9*Y73amp
        hN84amp = eta*N84amp*sp.Abs(c8)*vphi**8*Y84amp
        hN85amp = eta*N85amp*sp.Abs(c9)*vphi**9*Y75amp
        hN86amp = eta*N86amp*sp.Abs(c8)*vphi**8*Y86amp
        hN87amp = eta*N87amp*sp.Abs(c9)*vphi**9*Y77amp
        hN88amp = eta*N88amp*sp.Abs(c8)*vphi**8*Y88amp
        hN71amp = eta*N71amp*sp.Abs(c7)*vphi**7*Y71amp
        hN72amp = eta*N72amp*sp.Abs(c8)*vphi**8*Y62amp
        hN73amp = eta*N73amp*sp.Abs(c7)*vphi**7*Y73amp
        hN74amp = eta*N74amp*sp.Abs(c8)*vphi**8*Y64amp
        hN75amp = eta*N75amp*sp.Abs(c7)*vphi**7*Y75amp
        hN76amp = eta*N76amp*sp.Abs(c8)*vphi**8*Y66amp
        hN77amp = eta*N77amp*sp.Abs(c7)*vphi**7*Y77amp
        hN61amp = eta*N61amp*sp.Abs(c7)*vphi**7*Y51amp
        hN62amp = eta*N62amp*sp.Abs(c6)*vphi**6*Y62amp
        hN63amp = eta*N63amp*sp.Abs(c7)*vphi**7*Y53amp
        hN64amp = eta*N64amp*sp.Abs(c6)*vphi**6*Y64amp
        hN65amp = eta*N65amp*sp.Abs(c7)*vphi**7*Y55amp
        hN66amp = eta*N66amp*sp.Abs(c6)*vphi**6*Y66amp
        hN51amp = eta*N51amp*sp.Abs(c5)*vphi**5*Y51amp
        hN52amp = eta*N52amp*sp.Abs(c6)*vphi**6*Y42amp
        hN53amp = eta*N53amp*sp.Abs(c5)*vphi**5*Y53amp
        hN54amp = eta*N54amp*sp.Abs(c6)*vphi**6*Y44amp
        hN55amp = eta*N55amp*sp.Abs(c5)*vphi**5*Y55amp
        hN55 = eta*N55*c5*vphi**5*Y55
        hN41amp = eta*N41amp*sp.Abs(c5)*vphi**5*Y31amp
        hN42amp = eta*N42amp*sp.Abs(c4)*vphi**4*Y42amp
        hN43amp = eta*N43amp*sp.Abs(c5)*vphi**5*Y33amp
        hN43 = eta*N43*c5*vphi**5*Y33
        hN44amp = eta*N44amp*sp.Abs(c4)*vphi**4*Y44amp
        hN44 = eta*N44*c4*vphi**4*Y44
        hN31amp = eta*N31amp*sp.Abs(c3)*vphi**3*Y31amp
        hN32amp = eta*N32amp*sp.Abs(c4)*vphi**4*Y22amp
        hN32 = eta*N32*c4*vphi**4*Y22
        hN33amp = eta*N33amp*sp.Abs(c3)*vphi**3*Y33amp
        hN33 = eta*N33*c3*vphi**3*Y33
        hN21amp = eta*N21amp*sp.Abs(c3)*vphi**3*Y11amp
        hN21 = eta*N21*c3*vphi**3*Y11
        hN22amp = eta*N22amp*sp.Abs(c2)*vphi**2*Y22amp
        hN22 = eta*N22*c2*vphi**2*Y22
        h81amp = hN81amp*source_odd*T81amp*f81amp
        h82amp = hN82amp*source_even*T82amp*f82amp
        h83amp = hN83amp*source_odd*T83amp*f83amp
        h84amp = hN84amp*source_even*T84amp*f84amp
        h85amp = hN85amp*source_odd*T85amp*f85amp
        h86amp = hN86amp*source_even*T86amp*f86amp
        h87amp = hN87amp*source_odd*T87amp*f87amp
        h88amp = hN88amp*source_even*T88amp*f88amp
        h71amp = hN71amp*source_even*T71amp*f71amp
        h72amp = hN72amp*source_odd*T72amp*f72amp
        h73amp = hN73amp*source_even*T73amp*f73amp
        h74amp = hN74amp*source_odd*T74amp*f74amp
        h75amp = hN75amp*source_even*T75amp*f75amp
        h76amp = hN76amp*source_odd*T76amp*f76amp
        h77amp = hN77amp*source_even*T77amp*f77amp
        h61amp = hN61amp*source_odd*T61amp*f61amp
        h62amp = hN62amp*source_even*T62amp*f62amp
        h63amp = hN63amp*source_odd*T63amp*f63amp
        h64amp = hN64amp*source_even*T64amp*f64amp
        h65amp = hN65amp*source_odd*T65amp*f65amp
        h66amp = hN66amp*source_even*T66amp*f66amp
        h51amp = hN51amp*source_even*T51amp*f51amp
        h52amp = hN52amp*source_odd*T52amp*f52amp
        h53amp = hN53amp*source_even*T53amp*f53amp
        h54amp = hN54amp*source_odd*T54amp*f54amp
        h55amp = hN55amp*source_even*T55amp*f55amp
        h55 = hN55*source_even*T55*f55*sp.exp(1j*delta55)
        h41amp = hN41amp*source_odd*T41amp*f41amp
        h42amp = hN42amp*source_even*T42amp*f42amp
        h43amp = hN43amp*source_odd*T43amp*f43amp
        h43 = hN43*source_odd*T43*f43*sp.exp(1j*delta43)
        h44amp = hN44amp*source_even*T44amp*f44amp
        h44 = hN44*source_even*T44*f44*sp.exp(1j*delta44)
        h31amp = hN31amp*source_even*T31amp*f31amp
        h32amp = hN32amp*source_odd*T32amp*f32amp
        h32 = hN32*source_odd*T32*f32*sp.exp(1j*delta32)
        h33amp = hN33amp*source_even*T33amp*f33amp
        h33 = hN33*source_even*T33*f33*sp.exp(1j*delta33)
        h21amp = hN21amp*source_odd*T21amp*f21amp
        h21 = hN21*source_odd*T21*f21*sp.exp(1j*delta21)
        h22amp = hN22amp*source_even*T22amp*f22amp
        self.h22 = hN22*source_even*T22*f22*sp.exp(1j*delta22)
        self.flux=-(Omega**2/8/sp.pi)*(2*2*h22amp**2+1*1*h21amp**2+3*3*h33amp**2+2*2*h32amp**2+1*1*h31amp**2+4*4*h44amp**2+3*3*h43amp**2+2*2*h42amp**2+1*1*h41amp**2+5*5*h55amp**2+4*4*h54amp**2+3*3*h53amp**2+2*2*h52amp**2+1*1*h51amp**2+6*6*h66amp**2+5*5*h65amp**2+4*4*h64amp**2+3*3*h63amp**2+2*2*h62amp**2+1*1*h61amp**2+7*7*h77amp**2+6*6*h76amp**2+5*5*h75amp**2+4*4*h74amp**2+3*3*h73amp**2+2*2*h72amp**2+1*1*h71amp**2+8*8*h88amp**2+7*7*h87amp**2+6*6*h86amp**2+5*5*h85amp**2+4*4*h84amp**2+3*3*h83amp**2+2*2*h82amp**2+1*1*h81amp**2) / eta


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
        SEOBNRv5_aligned_spin_waveform_quantities().__dict__,
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
