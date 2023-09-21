"""
This module implements the  Goldberg formula for
  computing spin-weighted spherical harmonics.
  (https://aip.scitation.org/doi/10.1063/1.1705135)
Wikipedia also has an article on Spin-Weighted Spherical Hamronics:
  (https://en.wikipedia.org/w/index.php?title=Spin-weighted_spherical_harmonics&oldid=853425244)

Authors: Brandon Clark
         Zachariah B. Etienne
         zachetie **at** gmail **dot* com
"""
# Step 1: Initialize needed Python/NRPy+ modules
from typing import cast, Dict
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends


# Step 2: Define the Goldberg formula for spin-weighted spherical harmonics
#         (https://aip.scitation.org/doi/10.1063/1.1705135);
#         referenced & described in Wikipedia Spin-weighted spherical harmonics article:
#         https://en.wikipedia.org/w/index.php?title=Spin-weighted_spherical_harmonics&oldid=853425244
def Y(
    s: int,
    l: int,
    m: int,
    th: sp.Symbol,
    ph: sp.Symbol,
    GenerateMathematicaCode: bool = False,
) -> sp.Expr:
    Sum: sp.Expr = sp.sympify(0)
    for r in range(l - s + 1):
        if GenerateMathematicaCode:
            # Mathematica needs expression to be in terms of cotangent, so that code validation below
            #    yields identity with existing Mathematica notebook on spin-weighted spherical harmonics.
            Sum += (
                sp.binomial(l - s, r)
                * sp.binomial(l + s, r + s - m)
                * (-1) ** (l - r - s)
                * sp.exp(sp.I * m * ph)
                * sp.cot(th / 2) ** (2 * r + s - m)
            )
        else:
            # SymPy C code generation cannot handle the cotangent function, so define cot(th/2) as 1/tan(th/2):
            Sum += (
                sp.binomial(l - s, r)
                * sp.binomial(l + s, r + s - m)
                * (-1) ** (l - r - s)
                * sp.exp(sp.I * m * ph)
                / sp.tan(th / 2) ** (2 * r + s - m)
            )

    return cast(
        sp.Expr,
        (-1) ** m
        * sp.simplify(
            sp.sqrt(
                sp.factorial(l + m)
                * sp.factorial(l - m)
                * (2 * l + 1)
                / (4 * sp.pi * sp.factorial(l + s) * sp.factorial(l - s))
            )
            * sp.sin(th / 2) ** (2 * l)
            * Sum
        ),
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

    expr_dict: Dict[str, sp.Expr] = {}
    _th, _ph = sp.symbols("th ph")
    for _l in range(2, 9):
        for _m in range(-_l, +_l):
            expr_dict[f"Y_{{s=-2, l={_l}, m={_m}}}"] = Y(-2, _l, _m, _th, _ph)
    results_dict = ve.process_dictionary_of_expressions(
        expr_dict, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        # File basename. If this is set to "trusted_module_test1", then
        #   trusted results_dict will be stored in tests/trusted_module_test1.py
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )


# def SpinWeight_minus2_SphHarmonics(
#     maximum_l=8,
#     filename=os.path.join(
#         "SpinWeight_minus2_SphHarmonics", "SpinWeight_minus2_SphHarmonics.h"
#     ),
# ):
#     # Step 3: (DISABLED FOR NOW; PASSES TEST).
#     #         Code Validation against Mathematica notebook:
#     #         https://demonstrations.wolfram.com/versions/source.jsp?id=SpinWeightedSphericalHarmonics&version=0012
#
#     # # For the l=0 case m=0, otherwise there is a divide-by-zero in the Y() function above.
#     # print("FullSimplify[Y[-2, 0, 0, th, ph]-"+str(sp.mathematica_code(sp.simplify(Y(-2, 0, 0, th, ph,GenerateMathematicaCode=True))))+"] \n") # Agrees with Mathematica notebook for l = 0
#
#     # # Check the other cases
#     # for l in range(1,9): # Agrees with Mathematica notebook for  l = 1, 2, 4, 5, 6, 7, 8;
#     #     print("FullSimplify[Y[-2, "+str(l)+", m, th, ph]-("+
#     #           str(sp.mathematica_code(sp.simplify(Y(-2, l, m, th, ph, GenerateMathematicaCode=True)))).replace("binomial","Binomial").replace("factorial","Factorial")+")] \n")
#
#     # Step 4: Generating C Code function for computing
#     #         s=-2 spin-weighted spherical harmonics,
#     #         using NRPy+'s outputC() function.
#
#     outCparams = "preindent=3,outCfileaccess=a,outCverbose=False,includebraces=True"
#
#     with open(filename, "w") as file:
#         file.write(
#             """
# void SpinWeight_minus2_SphHarmonics(const int l, const int m, const REAL th, const REAL ph,
#                                    REAL *reYlmswm2_l_m, REAL *imYlmswm2_l_m) {
# if(l<0 || l>"""
#             + str(maximum_l)
#             + """ || m<-l || m>+l) {
#     printf("ERROR: SpinWeight_minus2_SphHarmonics handles only l=[0,"""
#             + str(maximum_l)
#             + """] and only m=[-l,+l] is defined.\\n");
#     printf("       You chose l=%d and m=%d, which is out of these bounds.\\n",l,m);
#     exit(1);
# }\n"""
#         )
#
#         file.write("switch(l) {\n")
#         for l in range(maximum_l + 1):  # Output values up to and including l=8.
#             file.write("    case " + str(l) + ":\n")
#             file.write("        switch(m) {\n")
#             for m in range(-l, l + 1):
#                 file.write("            case " + str(m) + ":\n")
#                 Y_m2_lm = Y(-2, l, m, th, ph)
#                 Cstring = outputC(
#                     [sp.re(Y_m2_lm), sp.im(Y_m2_lm)],
#                     ["*reYlmswm2_l_m", "*imYlmswm2_l_m"],
#                     "returnstring",
#                     outCparams,
#                 )
#                 file.write(Cstring)
#                 file.write("                  return;\n")
#             file.write("        }  /* End switch(m) */\n")
#         file.write("    } /* End switch(l) */\n")
#         file.write("} /* End function SpinWeight_minus2_SphHarmonics() */\n")
