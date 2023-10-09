"""
Construct symbolic expressions for the right-hand-side of the wave equation, in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com

License: BSD 2-Clause

The wave equation in curvilinear coordinates can be
  written as
hat{g}^{mu nu} partial_{mu} partial_{nu} u - hat{Gamma}^{tau} partial_{tau} u ,
where hat{Gamma}^{tau} is the *contracted* Christoffel symbol.

For reference metrics supported by NRPy+,
  hat{g}_{t nu} = -delta_{t nu},
  where delta_{t nu} is the Kronecker delta, we have
partial_t^2 u = hat{g}^{i j} partial_{i} partial_{j} u - hat{Gamma}^i partial_i u ,
  where Latin indices imply the spatial components *only*.

We break up the above equation into two first-order-in-time PDEs,
  defining v = partial_t u:

partial_t u = v
partial_t v = hat{g}^{i j} partial_{i} partial_{j} u - hat{Gamma}^i partial_i u
"""

# Step P1: Import needed modules:
import sympy as sp  # For symbolic computations
import nrpy.grid as gri  # NRPy+: Functionality for handling numerical grids
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
from nrpy.equations.wave_equation.CommonParams import (
    wavespeed,
)  # NRPy+: Common parameters for all WaveEquation modules (defines wavespeed)

# The name of this module ("WaveEquationCurvilinear") is given by __name__:
thismodule: str = __name__


# Specify RHSs as class variables,
# to enable access outside this
# function (e.g., for C code output)
class WaveEquationCurvilinear_RHSs:
    """Class sets up and stores sympy expressions for wave equation RHSs in curvilinear coordinates."""

    def __init__(self, CoordSystem: str, enable_rfm_precompute: bool) -> None:
        """
        Compute the right-hand sides (RHSs) of the scalar wave equation in curvilinear coordinates.

        :param CoordSystem: The coordinate system being used.

        .. note::
            Class variables uu_rhs and vv_rhs will be set in this function.
        """
        # Step 1: Set up the reference metric and
        #         quantities derived from the
        #         reference metric.
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]

        # Step 2: Compute the contracted Christoffel symbols:
        contractedGammahatU = ixp.zerorank1()
        for k in range(3):
            for i in range(3):
                for j in range(3):
                    contractedGammahatU[k] += (
                        rfm.ghatUU[i][j] * rfm.GammahatUDD[k][i][j]
                    )

        # Step 3: Register gridfunctions that are needed as input
        #         to the scalar wave RHS expressions.
        if "uu" not in gri.glb_gridfcs_dict:
            _uu, vv = gri.register_gridfunctions(
                ["uu", "vv"], group="EVOL", f_infinity=[2.0, 0.0], wavespeed=[1.0, 1.0]
            )
        else:
            vv = sp.Symbol("vv", real=True)

        # Step 4a: Declare the rank-1 indexed expression \partial_{i} u,
        #          Derivative variables like these must have an underscore
        #          in them, so the finite difference module can parse the
        #          variable name properly.
        uu_dD = ixp.declarerank1("uu_dD")

        # Step 4b: Declare the rank-2 indexed expression \partial_{ij} u,
        #          which is symmetric about interchange of indices i and j
        #          Derivative variables like these must have an underscore
        #          in them, so the finite difference module can parse the
        #          variable name properly.
        uu_dDD = ixp.declarerank2("uu_dDD", symmetry="sym01")

        # Step 5: Define right-hand sides for the evolution.
        # Step 5a: uu_rhs = vv:
        self.uu_rhs = vv
        # Step 5b: The right-hand side of the \partial_t v equation
        #          is given by:
        #          \hat{g}^{ij} \partial_i \partial_j u - \hat{\Gamma}^i \partial_i u.
        #          ^^^^^^^^^^^^ PART 1 ^^^^^^^^^^^^^^^^ ^^^^^^^^^^ PART 2 ^^^^^^^^^^^
        self.vv_rhs = sp.sympify(0)
        for i in range(3):
            # PART 2:
            self.vv_rhs -= contractedGammahatU[i] * uu_dD[i]
            for j in range(3):
                # PART 1:
                self.vv_rhs += rfm.ghatUU[i][j] * uu_dDD[i][j]

        self.vv_rhs *= wavespeed * wavespeed

        # Step 6: Generate C code for scalarwave evolution equations,
        #         print output to the screen (standard out, or stdout).
        # fin.FD_outputC("stdout",
        #                [lhrh(lhs=gri.gfaccess("rhs_gfs","uu"),rhs=uu_rhs),
        #                 lhrh(lhs=gri.gfaccess("rhs_gfs","vv"),rhs=vv_rhs)])


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

    for Coord in [
        "Spherical",
        "SinhSpherical",
        "Cartesian",
        "SinhCylindrical",
        "SinhSymTP",
    ]:
        RHS = WaveEquationCurvilinear_RHSs(
            CoordSystem=Coord, enable_rfm_precompute=False
        )
        results_dict = ve.process_dictionary_of_expressions(
            RHS.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )
