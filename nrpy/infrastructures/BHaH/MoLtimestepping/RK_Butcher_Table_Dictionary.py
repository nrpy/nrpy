"""
Construct a Python dictionary of Butcher tables for explicit RK methods.

Authors: Brandon Clark
         Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Gabriel M Steward
"""

# Step 1: Initialize needed Python/NRPy+ modules
from typing import Dict, List, Union, Tuple
import sympy as sp  # Import SymPy, Python's computer algebra system


# Step 2a: Generate a Dictionary of Butcher Tables for Explicit Runge Kutta Techniques
def generate_Butcher_tables(
    generate_adams_bashforth_method: bool = False,
    adams_bashforth_order: int = 7,
) -> Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]]:
    """
    Generate a dictionary of Butcher tables for Explicit Runge Kutta techniques.

    :param generate_adams_bashforth_method: If True, generate the Adams-Bashforth method. Default is False.
    :param adams_bashforth_order: The order of Adams-Bashforth method to generate. Default is 7.

    :return: A dictionary with method names as keys and Butcher tables as values.

    >>> Butcher_dict = generate_Butcher_tables(generate_adams_bashforth_method=True, adams_bashforth_order=9)
    >>> Butcher_table = Butcher_dict["AB"][0]
    >>> print(Butcher_table)
    [[1], [3/2, -1/2], [23/12, -4/3, 5/12], [55/24, -59/24, 37/24, -3/8], [1901/720, -1387/360, 109/30, -637/360, 251/720], [4277/1440, -2641/480, 4991/720, -3649/720, 959/480, -95/288], [198721/60480, -18637/2520, 235183/20160, -10754/945, 135713/20160, -5603/2520, 19087/60480], [16083/4480, -1152169/120960, 242653/13440, -296053/13440, 2102243/120960, -115747/13440, 32863/13440, -5257/17280], [14097247/3628800, -21562603/1814400, 47738393/1814400, -69927631/1814400, 862303/22680, -45586321/1814400, 19416743/1814400, -4832053/1814400, 1070017/3628800]]
    """
    # Initialize the dictionary Butcher_dict
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]] = {}

    # Disable black formatting for the Butcher tables
    # fmt: off

    # Step 2.a.i: Euler's Method
    Butcher_dict['Euler'] = (
    [[0],
     ["", 1]]
    , 1)

    # Step 2.a.ii: RK2 Heun's Method

    Butcher_dict['RK2 Heun'] = (
    [[0],
    [1, 1],
    ["", sp.Rational(1,2), sp.Rational(1,2)]]
    , 2)

    # Step 2.a.iii: RK2 Midpoint (MP) Method

    Butcher_dict['RK2 MP'] = (
    [[0],
    [sp.Rational(1,2), sp.Rational(1,2)],
    ["", 0, 1]]
    , 2)

    # Step 2.a.iv: RK2 Ralston's Method

    Butcher_dict['RK2 Ralston'] = (
    [[0],
    [sp.Rational(2,3), sp.Rational(2,3)],
    ["", sp.Rational(1,4), sp.Rational(3,4)]]
    , 2)

    # Step 2.a.v: Kutta's  Third-order Method

    Butcher_dict['RK3'] = (
    [[0],
    [sp.Rational(1,2), sp.Rational(1,2)],
    [1, -1, 2],
    ["", sp.Rational(1,6), sp.Rational(2,3), sp.Rational(1,6)]]
    , 3)

    # Step 2.a.vi: RK3 Heun's Method

    Butcher_dict['RK3 Heun'] = (
    [[0],
    [sp.Rational(1,3), sp.Rational(1,3)],
    [sp.Rational(2,3), 0, sp.Rational(2,3)],
    ["", sp.Rational(1,4), 0, sp.Rational(3,4)]]
    , 3)

    # Step 2.a.vii: RK3 Ralton's Method

    Butcher_dict['RK3 Ralston'] = (
    [[0],
    [sp.Rational(1,2), sp.Rational(1,2)],
    [sp.Rational(3,4), 0, sp.Rational(3,4)],
    ["", sp.Rational(2,9), sp.Rational(1,3), sp.Rational(4,9)]]
    , 3)

    # Step 2.a.viii: Strong Stability Preserving Runge-Kutta (SSPRK3) Method
    Butcher_dict['SSPRK3'] = (
    [[0],
    [1, 1],
    [sp.Rational(1,2), sp.Rational(1,4), sp.Rational(1,4)],
    ["", sp.Rational(1,6), sp.Rational(1,6), sp.Rational(2,3)]]
    , 3)

    # Step 2.a.ix: Classic RK4 Method

    Butcher_dict['RK4'] = (
    [[0],
    [sp.Rational(1,2), sp.Rational(1,2)],
    [sp.Rational(1,2), 0, sp.Rational(1,2)],
    [1, 0, 0, 1],
    ["", sp.Rational(1,6), sp.Rational(1,3), sp.Rational(1,3), sp.Rational(1,6)]]
    , 4)

    # Step 2.a.x:  RK5 Dormand-Prince Method

    Butcher_dict['DP5'] = (
    [[0],
    [sp.Rational(1,5), sp.Rational(1,5)],
    [sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
    [sp.Rational(4,5), sp.Rational(44,45), sp.Rational(-56,15), sp.Rational(32,9)],
    [sp.Rational(8,9), sp.Rational(19372,6561), sp.Rational(-25360,2187), sp.Rational(64448,6561), sp.Rational(-212,729)],
    [1, sp.Rational(9017,3168), sp.Rational(-355,33), sp.Rational(46732,5247), sp.Rational(49,176), sp.Rational(-5103,18656)],
    [1, sp.Rational(35,384), 0, sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84)],
    ["", sp.Rational(35,384), 0, sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84), 0]]
    , 5)

    # Step 2.a.xi:  RK5 Dormand-Prince Method Alternative

    Butcher_dict['DP5alt'] = (
    [[0],
    [sp.Rational(1,10), sp.Rational(1,10)],
    [sp.Rational(2,9), sp.Rational(-2, 81), sp.Rational(20, 81)],
    [sp.Rational(3,7), sp.Rational(615, 1372), sp.Rational(-270, 343), sp.Rational(1053, 1372)],
    [sp.Rational(3,5), sp.Rational(3243, 5500), sp.Rational(-54, 55), sp.Rational(50949, 71500), sp.Rational(4998, 17875)],
    [sp.Rational(4, 5), sp.Rational(-26492, 37125), sp.Rational(72, 55), sp.Rational(2808, 23375), sp.Rational(-24206, 37125), sp.Rational(338, 459)],
    [1, sp.Rational(5561, 2376), sp.Rational(-35, 11), sp.Rational(-24117, 31603), sp.Rational(899983, 200772), sp.Rational(-5225, 1836), sp.Rational(3925, 4056)],
    ["", sp.Rational(821, 10800), 0, sp.Rational(19683, 71825), sp.Rational(175273, 912600), sp.Rational(395, 3672), sp.Rational(785, 2704), sp.Rational(3, 50)]]
    , 5)

    # Step 2.a.xiii:  RK5 Cash-Karp Method

    Butcher_dict['CK5'] = (
    [[0],
    [sp.Rational(1,5), sp.Rational(1,5)],
    [sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
    [sp.Rational(3,5), sp.Rational(3,10), sp.Rational(-9,10), sp.Rational(6,5)],
    [1, sp.Rational(-11,54), sp.Rational(5,2), sp.Rational(-70,27), sp.Rational(35,27)],
    [sp.Rational(7,8), sp.Rational(1631,55296), sp.Rational(175,512), sp.Rational(575,13824), sp.Rational(44275,110592), sp.Rational(253,4096)],
    ["",sp.Rational(37,378), 0, sp.Rational(250,621), sp.Rational(125,594), 0, sp.Rational(512,1771)]]
    , 5)

    # Step 2.a.xiv:  RK6 Dormand-Prince Method

    Butcher_dict['DP6'] = (
    [[0],
    [sp.Rational(1,10), sp.Rational(1,10)],
    [sp.Rational(2,9), sp.Rational(-2, 81), sp.Rational(20, 81)],
    [sp.Rational(3,7), sp.Rational(615, 1372), sp.Rational(-270, 343), sp.Rational(1053, 1372)],
    [sp.Rational(3,5), sp.Rational(3243, 5500), sp.Rational(-54, 55), sp.Rational(50949, 71500), sp.Rational(4998, 17875)],
    [sp.Rational(4, 5), sp.Rational(-26492, 37125), sp.Rational(72, 55), sp.Rational(2808, 23375), sp.Rational(-24206, 37125), sp.Rational(338, 459)],
    [1, sp.Rational(5561, 2376), sp.Rational(-35, 11), sp.Rational(-24117, 31603), sp.Rational(899983, 200772), sp.Rational(-5225, 1836), sp.Rational(3925, 4056)],
    [1, sp.Rational(465467, 266112), sp.Rational(-2945, 1232), sp.Rational(-5610201, 14158144), sp.Rational(10513573, 3212352), sp.Rational(-424325, 205632), sp.Rational(376225, 454272), 0],
    ["", sp.Rational(61, 864), 0, sp.Rational(98415, 321776), sp.Rational(16807, 146016), sp.Rational(1375, 7344), sp.Rational(1375, 5408), sp.Rational(-37, 1120), sp.Rational(1,10)]]
    , 6)

    #  Step 2.a.xv:  RK6 Luther's Method

    q = sp.sqrt(21)
    Butcher_dict['L6'] = (
    [[0],
    [1, 1],
    [sp.Rational(1,2), sp.Rational(3,8), sp.Rational(1,8)],
    [sp.Rational(2,3), sp.Rational(8,27), sp.Rational(2,27), sp.Rational(8,27)],
    [(7 - q)/14, (-21 + 9*q)/392, (-56 + 8*q)/392, (336 -48*q)/392, (-63 + 3*q)/392],
    [(7 + q)/14, (-1155 - 255*q)/1960, (-280 -  40*q)/1960, (-320*q)/1960, (63 + 363*q)/1960, (2352 + 392*q)/1960],
    [1, ( 330 + 105*q)/180, sp.Rational(2,3), (-200 + 280*q)/180, (126 - 189*q)/180, (-686 - 126*q)/180, (490 -  70*q)/180],
    ["", sp.Rational(1, 20), 0, sp.Rational(16, 45), 0, sp.Rational(49, 180), sp.Rational(49, 180), sp.Rational(1, 20)]]
    , 6)

    # Step 2.a.xvi: RK8 Dormand-Prince Method

    Butcher_dict['DP8']=(
    [[0],
    [sp.Rational(1, 18), sp.Rational(1, 18)],
    [sp.Rational(1, 12), sp.Rational(1, 48), sp.Rational(1, 16)],
    [sp.Rational(1, 8), sp.Rational(1, 32), 0, sp.Rational(3, 32)],
    [sp.Rational(5, 16), sp.Rational(5, 16), 0, sp.Rational(-75, 64), sp.Rational(75, 64)],
    [sp.Rational(3, 8), sp.Rational(3, 80), 0, 0, sp.Rational(3, 16), sp.Rational(3, 20)],
    [sp.Rational(59, 400), sp.Rational(29443841, 614563906), 0, 0, sp.Rational(77736538, 692538347), sp.Rational(-28693883, 1125000000), sp.Rational(23124283, 1800000000)],
    [sp.Rational(93, 200), sp.Rational(16016141, 946692911), 0, 0, sp.Rational(61564180, 158732637), sp.Rational(22789713, 633445777), sp.Rational(545815736, 2771057229), sp.Rational(-180193667, 1043307555)],
    [sp.Rational(5490023248, 9719169821), sp.Rational(39632708, 573591083), 0, 0, sp.Rational(-433636366, 683701615), sp.Rational(-421739975, 2616292301), sp.Rational(100302831, 723423059), sp.Rational(790204164, 839813087), sp.Rational(800635310, 3783071287)],
    [sp.Rational(13, 20), sp.Rational(246121993, 1340847787), 0, 0, sp.Rational(-37695042795, 15268766246), sp.Rational(-309121744, 1061227803), sp.Rational(-12992083, 490766935), sp.Rational(6005943493, 2108947869), sp.Rational(393006217, 1396673457), sp.Rational(123872331, 1001029789)],
    [sp.Rational(1201146811, 1299019798), sp.Rational(-1028468189, 846180014), 0, 0, sp.Rational(8478235783, 508512852), sp.Rational(1311729495, 1432422823), sp.Rational(-10304129995, 1701304382), sp.Rational(-48777925059, 3047939560), sp.Rational(15336726248, 1032824649), sp.Rational(-45442868181, 3398467696), sp.Rational(3065993473, 597172653)],
    [1, sp.Rational(185892177, 718116043), 0, 0, sp.Rational(-3185094517, 667107341), sp.Rational(-477755414, 1098053517), sp.Rational(-703635378, 230739211), sp.Rational(5731566787, 1027545527), sp.Rational(5232866602, 850066563), sp.Rational(-4093664535, 808688257), sp.Rational(3962137247, 1805957418), sp.Rational(65686358, 487910083)],
    [1, sp.Rational(403863854, 491063109), 0, 0, sp.Rational(-5068492393, 434740067), sp.Rational(-411421997, 543043805), sp.Rational(652783627, 914296604), sp.Rational(11173962825, 925320556), sp.Rational(-13158990841, 6184727034), sp.Rational(3936647629, 1978049680), sp.Rational(-160528059, 685178525), sp.Rational(248638103, 1413531060), 0],
    ["", sp.Rational(14005451, 335480064), 0, 0, 0, 0, sp.Rational(-59238493, 1068277825), sp.Rational(181606767, 758867731), sp.Rational(561292985, 797845732), sp.Rational(-1041891430, 1371343529), sp.Rational(760417239, 1151165299), sp.Rational(118820643, 751138087), sp.Rational(-528747749, 2220607170), sp.Rational(1, 4)]]
    , 8)

    # Step 3.a:  Generating a Dictionary of Butcher Tables for Explicit Runge Kutta Techniques

    # Step 3.a.i: Adaptive Heun-Euler Method

    Butcher_dict['AHE'] = (
    [[0],
    [1, 1],
    ["", sp.Rational(1,2), sp.Rational(1,2)],
    ["", 1, 0]]
    , 2)

    # Step 3.a.ii: Adaptive Bogacki-Shampine Method

    Butcher_dict['ABS'] = (
    [[0],
    [sp.Rational(1,2), sp.Rational(1,2)],
    [sp.Rational(3,4), 0, sp.Rational(3,4)],
    [1, sp.Rational(2,9), sp.Rational(1,3), sp.Rational(4,9)],
    ["", sp.Rational(2,9), sp.Rational(1,3), sp.Rational(4,9), 0],
    ["", sp.Rational(7,24), sp.Rational(1,4), sp.Rational(1,3), sp.Rational(1,8)]]
    , 3)

    # Step 3.a.iii: Adaptive Runge-Kutta-Fehlberg

    Butcher_dict['ARKF'] = (
    [[0],
    [sp.Rational(1,4), sp.Rational(1,4)],
    [sp.Rational(3,8), sp.Rational(3,32), sp.Rational(9,32)],
    [sp.Rational(12,13), sp.Rational(1932,2197), sp.Rational(-7200,2197), sp.Rational(7296,2197)],
    [1, sp.Rational(439,216), -8, sp.Rational(3680,513), sp.Rational(-845,4104)],
    [sp.Rational(1,2), sp.Rational(-8,27), 2, sp.Rational(-3544,2565), sp.Rational(1859,4104), sp.Rational(-11,40)],
    ["", sp.Rational(16,135), 0, sp.Rational(6656,12825), sp.Rational(28561,56430), sp.Rational(-9,50), sp.Rational(2,55)],
    ["", sp.Rational(25,216), 0, sp.Rational(1408,2565), sp.Rational(2197,4104), sp.Rational(-1,5), 0]]
    , 5)

    # Step 3.a.iv: Adaptive Cash-Karp

    Butcher_dict['ACK'] = (
    [[0],
    [sp.Rational(1,5), sp.Rational(1,5)],
    [sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
    [sp.Rational(3,5), sp.Rational(3,10), sp.Rational(-9,10), sp.Rational(6,5)],
    [1, sp.Rational(-11,54), sp.Rational(5,2), sp.Rational(-70,27), sp.Rational(35,27)],
    [sp.Rational(7,8), sp.Rational(1631,55296), sp.Rational(175,512), sp.Rational(575,13824), sp.Rational(44275,110592), sp.Rational(253,4096)],
    ["",sp.Rational(37,378), 0, sp.Rational(250,621), sp.Rational(125,594), 0, sp.Rational(512,1771)],
    ["",sp.Rational(2825,27648), 0, sp.Rational(18575,48384), sp.Rational(13525,55296), sp.Rational(277,14336), sp.Rational(1,4)]]
    , 5)

    # Step 3.a.v: Adaptive Dormand-Prince 5(4)

    Butcher_dict['ADP5'] = (
    [[0],
    [sp.Rational(1,5), sp.Rational(1,5)],
    [sp.Rational(3,10),sp.Rational(3,40), sp.Rational(9,40)],
    [sp.Rational(4,5), sp.Rational(44,45), sp.Rational(-56,15), sp.Rational(32,9)],
    [sp.Rational(8,9), sp.Rational(19372,6561), sp.Rational(-25360,2187), sp.Rational(64448,6561), sp.Rational(-212,729)],
    [1, sp.Rational(9017,3168), sp.Rational(-355,33), sp.Rational(46732,5247), sp.Rational(49,176), sp.Rational(-5103,18656)],
    [1, sp.Rational(35,384), 0, sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84)],
    ["", sp.Rational(35,384), 0, sp.Rational(500,1113), sp.Rational(125,192), sp.Rational(-2187,6784), sp.Rational(11,84), 0],
    ["", sp.Rational(5179,57600), 0, sp.Rational(7571,16695), sp.Rational(393,640), sp.Rational(-92097,339200), sp.Rational(187,2100), sp.Rational(1,40)]]
    , 5)

    # Step 3.a.vi: Adaptive Dormand-Prince 8(7)

    Butcher_dict['ADP8']=(
    [[0],
    [sp.Rational(1, 18), sp.Rational(1, 18)],
    [sp.Rational(1, 12), sp.Rational(1, 48), sp.Rational(1, 16)],
    [sp.Rational(1, 8), sp.Rational(1, 32), 0, sp.Rational(3, 32)],
    [sp.Rational(5, 16), sp.Rational(5, 16), 0, sp.Rational(-75, 64), sp.Rational(75, 64)],
    [sp.Rational(3, 8), sp.Rational(3, 80), 0, 0, sp.Rational(3, 16), sp.Rational(3, 20)],
    [sp.Rational(59, 400), sp.Rational(29443841, 614563906), 0, 0, sp.Rational(77736538, 692538347), sp.Rational(-28693883, 1125000000), sp.Rational(23124283, 1800000000)],
    [sp.Rational(93, 200), sp.Rational(16016141, 946692911), 0, 0, sp.Rational(61564180, 158732637), sp.Rational(22789713, 633445777), sp.Rational(545815736, 2771057229), sp.Rational(-180193667, 1043307555)],
    [sp.Rational(5490023248, 9719169821), sp.Rational(39632708, 573591083), 0, 0, sp.Rational(-433636366, 683701615), sp.Rational(-421739975, 2616292301), sp.Rational(100302831, 723423059), sp.Rational(790204164, 839813087), sp.Rational(800635310, 3783071287)],
    [sp.Rational(13, 20), sp.Rational(246121993, 1340847787), 0, 0, sp.Rational(-37695042795, 15268766246), sp.Rational(-309121744, 1061227803), sp.Rational(-12992083, 490766935), sp.Rational(6005943493, 2108947869), sp.Rational(393006217, 1396673457), sp.Rational(123872331, 1001029789)],
    [sp.Rational(1201146811, 1299019798), sp.Rational(-1028468189, 846180014), 0, 0, sp.Rational(8478235783, 508512852), sp.Rational(1311729495, 1432422823), sp.Rational(-10304129995, 1701304382), sp.Rational(-48777925059, 3047939560), sp.Rational(15336726248, 1032824649), sp.Rational(-45442868181, 3398467696), sp.Rational(3065993473, 597172653)],
    [1, sp.Rational(185892177, 718116043), 0, 0, sp.Rational(-3185094517, 667107341), sp.Rational(-477755414, 1098053517), sp.Rational(-703635378, 230739211), sp.Rational(5731566787, 1027545527), sp.Rational(5232866602, 850066563), sp.Rational(-4093664535, 808688257), sp.Rational(3962137247, 1805957418), sp.Rational(65686358, 487910083)],
    [1, sp.Rational(403863854, 491063109), 0, 0, sp.Rational(-5068492393, 434740067), sp.Rational(-411421997, 543043805), sp.Rational(652783627, 914296604), sp.Rational(11173962825, 925320556), sp.Rational(-13158990841, 6184727034), sp.Rational(3936647629, 1978049680), sp.Rational(-160528059, 685178525), sp.Rational(248638103, 1413531060), 0],
    ["", sp.Rational(14005451, 335480064), 0, 0, 0, 0, sp.Rational(-59238493, 1068277825), sp.Rational(181606767, 758867731), sp.Rational(561292985, 797845732), sp.Rational(-1041891430, 1371343529), sp.Rational(760417239, 1151165299), sp.Rational(118820643, 751138087), sp.Rational(-528747749, 2220607170), sp.Rational(1, 4)],
    ["", sp.Rational(13451932, 455176623), 0, 0, 0, 0, sp.Rational(-808719846, 976000145), sp.Rational(1757004468, 5645159321), sp.Rational(656045339, 265891186), sp.Rational(-3867574721, 1518517206), sp.Rational(465885868, 322736535), sp.Rational(53011238, 667516719), sp.Rational(2, 45), 0]]
    , 8)

    # fmt: on
    # Check if the Adams-Bashforth method is required to be generated
    if generate_adams_bashforth_method:
        # Step 4: Generate coefficients for Adams-Bashforth method
        # Note: Adams-Bashforth methods are explicit multi-step methods
        # which predict the future value of the function using information from several preceding steps.
        # This contrasts with methods like Runge-Kutta, which are single-step and use only the current step's information.

        # Initial matrix setup. The order-1 method isn't used since it's just a single "1".
        pythonButcher: List[List[Union[sp.Basic, int, str]]]
        pythonButcher = [
            [0, 0],
            [0, 0],
        ]

        # Define the symbolic variable 'x' which will be used for mathematical manipulations.
        x = sp.symbols("x")

        # Loop through each order up to the desired Adams-Bashforth order
        for n in range(1, adams_bashforth_order + 1):
            # Compute coefficients for each step up to the current order 'n'
            for j in range(n):
                # Create an expression that's a product of (x+i) terms. This is a part of the Adams-Bashforth formula.
                # 'x' is the symbolic representation of time steps.
                expr = x * sp.prod(x + i for i in range(n) if i != j) / x

                # Integrate the expression over [0, 1]. This helps in evaluating coefficients for Adams-Bashforth.
                expr2 = sp.integrate(expr, (x, 0, 1))

                # Apply the combinatorial factor to the integral. This is part of the Adams-Bashforth formula.
                expr2 *= sp.Rational(
                    (-1) ** j, (sp.factorial(j) * sp.factorial(n - j - 1))
                )

                # If the current matrix isn't large enough to accommodate the current order, expand it.
                if len(pythonButcher) < n:
                    pythonButcher.append([sp.Rational(0) for _ in range(n)])

                # Insert the computed coefficient into the matrix.
                # Depending on the order and step, either append a new coefficient or replace an existing one.
                if len(pythonButcher[n - 1]) < j + 1:
                    pythonButcher[n - 1].append(expr2)
                else:
                    pythonButcher[n - 1][j] = expr2

        # Cleanup: Remove unnecessary zeros which were added initially in the matrix setup.
        pythonButcher[0] = [1]

        # Store the computed matrix and the order in a dictionary for future reference or use.
        Butcher_dict["AB"] = (pythonButcher, adams_bashforth_order)

    return Butcher_dict


# Validation
def validate(
    Butcher_dict: Dict[str, Tuple[List[List[Union[sp.Basic, int, str]]], int]],
    included_keys: List[str],
    verbose: bool = False,
) -> None:
    """
    Validate a Runge-Kutta method against the exact solution of an ODE.

    :param Butcher_dict: A dictionary of Butcher tables.
    :param included_keys: A list of keys to access the desired Butcher table in Butcher_dict.
    :param verbose: A boolean indicating whether to print verbose output.

    Doctests:
    >>> Butcher_dict = generate_Butcher_tables()
    >>> # validate(Butcher_dict, included_keys=["DP5", "DP5alt", "CK5", "DP6", "L6", "DP8", "AHE", "ABS", "ARKF", "ACK", "ADP5", "ADP8", "AB"], verbose=False)
    >>> validate(Butcher_dict, included_keys=["Euler", "RK2 Heun", "RK3", "RK3 Ralston", "RK4"], verbose=False)
    Euler: PASSED VALIDATION
    RK2 Heun: PASSED VALIDATION
    RK3: PASSED VALIDATION
    RK3 Ralston: PASSED VALIDATION
    RK4: PASSED VALIDATION
    """
    # Defining the right-hand side of the ODE using lambdas
    rhs_dict = {
        "ypt": lambda y, t: y + t,
        "eypt": lambda y, t: sp.exp(1.0 * (y + t)),
        "y": lambda y, _t: y,
        "tpoly6": lambda _y, t: (
            2 * t**6 - 389 * t**5 + 15 * t**4 - 22 * t**3 + 81 * t**2 - t + 42
        ),
        "tpoly5": lambda _y, t: t**5 + t**4 + t**3 + t**2 + t + 1,
        "ty2pt": lambda y, t: t * y**2 + t,
        "ymcost": lambda y, t: y - sp.cos(t),
        "ypsint": lambda y, t: y + sp.sin(t),
    }

    initial_dict = {
        "ypt": (1.0, 0.0),
        "eypt": (1.0, 0.0),
        "y": (1.0, 0.0),
        "tpoly6": (1.0, 0.0),
        "tpoly5": (1.0, 0.0),
        "ty2pt": (1.0, 0.0),
        "ymcost": (1.0, 0.0),
        "ypsint": (0.0, 0.0),
    }

    for Butcher_key, Butcher_value in Butcher_dict.items():
        if Butcher_key in included_keys:
            Butcher = Butcher_value[0]
            for rhs_dict_key, rhs_dict_value in rhs_dict.items():
                # print(f"Validating {Butcher_key}, {rhs_dict_key}...")
                y = sp.Function("y")
                t, dt = sp.symbols("t dt")
                yn, tn = initial_dict[rhs_dict_key]

                if verbose:
                    print(
                        f"RK METHOD: {Butcher_key}: Solving y'(t) = {rhs_dict_value(y(t), t)}, y({tn})={yn},"
                    )

                # 1. First we solve the ODE exactly
                sol = sp.dsolve(sp.Eq(y(t).diff(t), rhs_dict_value(y(t), t)), y(t)).rhs
                constants = sp.solve([sol.subs(t, tn) - yn])
                if isinstance(constants, list):
                    exact = sol.subs(constants[0].items())
                else:
                    exact = sol.subs(constants.items())

                # 2. Now we solve the ODE numerically using specified Butcher table

                # Access the requested Butcher table
                # Determine number of predictor-corrector steps
                L = len(Butcher) - 1
                # Set a temporary array for update values
                k = sp.zeros(L, 1)
                # Initialize the updated solution
                ynp1 = 0.0
                for i in range(L):
                    # Initialize and approximate update for solution
                    yhat = float(yn)
                    for j in range(i):
                        # Update yhat for solution using a_ij Butcher table coefficients
                        yhat += Butcher[i][j + 1] * k[j]
                        # if Butcher_key == "DP8" or Butcher_key == "L6":
                        # Otherwise the adding of fractions kills performance.
                        yhat = 1.0 * sp.N(yhat, 20)
                    # Determine the next corrector variable k_i using c_i Butcher table coefficients
                    k[i] = dt * rhs_dict[rhs_dict_key](yhat, tn + Butcher[i][0] * dt)
                    # Update the solution at the next iteration ynp1 using Butcher table coefficients
                    ynp1 = ynp1 + Butcher[L][i + 1] * k[i]

                # Finish determining the solution for the next iteration
                ynp1 = ynp1 + yn

                # Determine the order of the RK method
                local_truncation_order = Butcher_dict[Butcher_key][1] + 1
                # Produces Taylor series of exact solution at t=tn about t = t_0 with the specified order
                exact_series = sp.series(
                    exact.subs(t, dt),
                    dt,
                    initial_dict[rhs_dict_key][1],
                    local_truncation_order,
                ).removeO()
                num_series = sp.series(ynp1, dt, 0, local_truncation_order).removeO()
                diff_with_roundoff = exact_series - num_series

                def print_verbose_output(
                    local_truncation_order: int,
                    exact: sp.Basic,
                    exact_series: sp.Basic,
                    num_series: sp.Basic,
                    diff_with_roundoff: sp.Basic,
                ) -> str:
                    output_lines = [
                        f" the first nonzero term should have local truncation error proportional to O(dt^{local_truncation_order}) or a higher power of dt.",
                        f" EXACT: {exact}",
                        f" EXACT SERIES: {exact_series}",
                        f" NUMCL SERIES: {num_series}",
                        " Subtracting the numerical result from the exact Taylor expansion, we find a local truncation error of",
                        " (Coefficients of order 1e-15 or less may generally be ignored, as these are at roundoff error):",
                        sp.pretty(diff_with_roundoff)
                        # Using `pretty` instead of `pretty_print` since the latter prints directly to the console
                    ]
                    return "\n".join(output_lines)

                for n in range(local_truncation_order):
                    coeff = diff_with_roundoff.coeff(dt, n)

                    # Set coefficients to zero if their magnitude is less than 1e-14
                    if sp.Abs(coeff) < 1e-12:
                        coeff = 0

                    if coeff != 0:
                        print(
                            "BAD:",
                            Butcher_key,
                            rhs_dict_key,
                            n,
                            coeff,
                            local_truncation_order,
                        )
                        raise ValueError(
                            print_verbose_output(
                                local_truncation_order,
                                exact,
                                exact_series,
                                num_series,
                                diff_with_roundoff,
                            )
                        )

                if verbose:
                    print_verbose_output(
                        local_truncation_order,
                        exact,
                        exact_series,
                        num_series,
                        diff_with_roundoff,
                    )

            print(f"{Butcher_key}: PASSED VALIDATION")


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
