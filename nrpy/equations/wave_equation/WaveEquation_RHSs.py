"""
Generating C code for the right-hand-side
 of the scalar wave equation, in
 ***Cartesian*** coordinates, in
 arbitrary spatial dimension
 (up to four spatial dimensions
  supported at the moment)

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com

License: BSD 2-Clause

partial_t^2 u  =  c^2 nabla^2 u

where u is a function of space and time
and nabla is the Laplacian differential operator.

We rewrite as set of two first-order-in-time PDEs,
defining partial_t u = v.

Then the above becomes

partial_t u = v
partial_t v = c^2 nabla^2 u

where u and v are functions of space and time,
and nabla is the Laplacian differential operator.
"""

# Step P1: Import needed modules:
import sympy as sp
import nrpy.grid as gri  # NRPy+: Functionality for handling numerical grids
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

# NRPy+: Common parameters for all WaveEquation modules (defines wavespeed)
from nrpy.equations.wave_equation.CommonParams import wavespeed

# The name of this module ("WaveEquation") is given by __name__:
thismodule = __name__


class WaveEquation_RHSs:
    """Class for storing wave equation RHS equations."""

    def __init__(self) -> None:
        """
        Main function for scalar wave RHS expressions. It defines and computes the right-hand sides
        for the scalar wave equation evolution.

        :return: None

        .. note::
            Class variables uu_rhs and vv_rhs will be set in this function.
        """

        # Step 1: Declare the rank-2 indexed expression \partial_{ij} u,
        #         which is symmetric about the interchange of indices i and j.
        #         Derivative variables like these must have an underscore
        #         in them, so the finite difference module can parse the
        #         variable name properly.
        uu_dDD = ixp.declarerank2("uu_dDD", symmetry="sym01")

        # Step 2: Specify RHSs as global variables,
        #         to enable access outside this function (e.g., for C code output).

        # Step 3: Define right-hand sides for the evolution.
        if "uu" not in gri.glb_gridfcs_dict:
            _uu, vv = gri.register_gridfunctions(
                ["uu", "vv"], group="EVOL", f_infinity=[2.0, 0.0], wavespeed=[1.0, 1.0]
            )
        else:
            vv = sp.Symbol("vv", real=True)
        self.uu_rhs = vv
        self.vv_rhs = sp.sympify(0)
        for i in range(3):
            self.vv_rhs += wavespeed * wavespeed * uu_dDD[i][i]


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

    for RHS_type in ["WaveEquation"]:
        RHS = WaveEquation_RHSs()
        results_dict = ve.process_dictionary_of_expressions(
            RHS.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{RHS_type}",
            results_dict,
        )
