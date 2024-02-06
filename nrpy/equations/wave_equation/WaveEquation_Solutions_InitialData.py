"""
Construct symbolic expressions for solution for the wave equation, as a function of time and Cartesian coordinates.

These functions are used to set up initial data and the exact solution at any later time.

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Thiago Assumpcao
         assumpcaothiago **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
from typing import Tuple
import sys
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.indexedexp as ixp
import nrpy.grid as gri

# NRPy+: Common parameters for all WaveEquation modules (defines wavespeed)
from nrpy.equations.wave_equation.CommonParams import wavespeed

# The name of this module ("InitialData") is given by __name__:
thismodule = __name__


class WaveEquation_solution_Cartesian:
    """
    Set up wave equation solution as a function of time & Cartesian coordinates, based on the given parameters.

    :param WaveType: The type of wave. Can be either "PlaneWave" or "SphericalGaussian"
    :param default_k0: The default value for k0. Default is 1
    :param default_k1: The default value for k1. Default is 1
    :param default_k2: The default value for k2. Default is 1
    :param default_sigma: The default value for sigma. Default is 3
    """

    def __init__(
        self,
        WaveType: str = "PlaneWave",
        default_k0: float = 1.0,
        default_k1: float = 1.0,
        default_k2: float = 1.0,
        default_sigma: float = 3.0,
    ):
        if "uu" not in gri.glb_gridfcs_dict:
            gri.register_gridfunctions(
                ["uu", "vv"], group="EVOL", f_infinity=[2.0, 0.0], wavespeed=[1.0, 1.0]
            )

        if WaveType == "PlaneWave":
            self.uu_exactsoln, self.vv_exactsoln = PlaneWave(
                default_k0=default_k0,
                default_k1=default_k1,
                default_k2=default_k2,
            )
        elif WaveType == "SphericalGaussian":
            (
                self.uu_exactsoln,
                self.vv_exactsoln,
                self.uu_exactsoln_r0,
                self.vv_exactsoln_r0,
            ) = SphericalGaussian(default_sigma=default_sigma)
        else:
            raise ValueError(
                f"Error: WaveEquation initial data WaveType={WaveType} not supported."
            )


# Set up spherically-symmetric Gaussian initial data in Cartesian coordinates
def SphericalGaussian(
    default_sigma: float = 3.0,
) -> Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """
    Set up initial data for a spherically-symmetric Gaussian wave.

    :param default_sigma: The default value for sigma, defaults to 3
    """
    # Step 1: Set up Cartesian coordinates (x, y, z) = (xCart0, xCart1, xCart2)
    xCart = ixp.declarerank1("xCart")

    # Step 2: Declare free parameters intrinsic to these initial data
    # provided as a C parameter by MoLtimestepping.MoL
    time = sp.symbols("time", real=True)
    sigma = par.register_CodeParameter(
        c_type_alias="REAL",
        module=thismodule,
        name="sigma",
        defaultvalue=default_sigma,
        commondata=True,
    )

    # Step 3: Compute r
    r = sp.sympify(0)
    for i in range(3):
        r += xCart[i] ** 2
    r = sp.sqrt(r)

    # Step 4: Set initial data for uu and vv, where vv_exactsoln = \partial_t uu_exactsoln.
    # uu_exactsoln = (r - wavespeed*time)/r * sp.exp(- (r - wavespeed*time)**2 / (2*sigma**2) )
    # By convention, limit(uu, r->infinity) = 2. This ensures that relative error is well defined.
    uu_exactsoln = (
        +((r - wavespeed * time) / r)
        * sp.exp(-((r - wavespeed * time) ** 2) / (2 * sigma**2))
        + ((r + wavespeed * time) / r)
        * sp.exp(-((r + wavespeed * time) ** 2) / (2 * sigma**2))
    ) + sp.sympify(
        2
    )  # Adding 2 ensures relative error is well defined (solution does not cross zero)
    vv_exactsoln = sp.diff(uu_exactsoln, time)

    # In the limit r->0, this equation is undefined, but we can use L'Hôpital's rule to
    # find the solution. The denominator is 1/r, which becomes 1. We evaluate the numerator
    # using sympy.
    rvar = sp.symbols("rvar", real=True)
    uu_exactsoln_r0 = sp.diff(
        +(rvar - wavespeed * time)
        * sp.exp(-((rvar - wavespeed * time) ** 2) / (2 * sigma**2))
        + (rvar + wavespeed * time)
        * sp.exp(-((rvar + wavespeed * time) ** 2) / (2 * sigma**2)),
        rvar,
    ) + sp.sympify(
        2
    )  # Adding 2 ensures relative error is well defined (solution does not cross zero)
    uu_exactsoln_r0 = uu_exactsoln_r0.subs(rvar, r)
    vv_exactsoln_r0 = sp.diff(uu_exactsoln_r0, time)

    return uu_exactsoln, vv_exactsoln, uu_exactsoln_r0, vv_exactsoln_r0


# Set up monochromatic plane-wave initial data
def PlaneWave(
    default_k0: float = 1,
    default_k1: float = 1,
    default_k2: float = 1,
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Set up initial data for a monochromatic plane wave.

    :param default_k0: The default value for k0, defaults to 1
    :param default_k1: The default value for k1, defaults to 1
    :param default_k2: The default value for k2, defaults to 1
    """
    # Step 1: Set up Cartesian coordinates (x, y, z) = (xCart0, xCart1, xCart2)
    xCart = ixp.declarerank1("xCart")

    # Step 2: Declare free parameters intrinsic to these initial data
    time = sp.symbols(
        "time", real=True
    )  # provided as a C parameter by MoLtimestepping.MoL
    kk = par.register_CodeParameters(
        c_type_alias="REAL",
        module=thismodule,
        names=["kk0", "kk1", "kk2"],
        defaultvalues=[default_k0, default_k1, default_k2],
        commondata=True,
    )

    # Step 3: Normalize the k vector
    kk_norm_factor = sp.sqrt(kk[0] ** 2 + kk[1] ** 2 + kk[2] ** 2)

    # Step 4: Compute k_norm.x
    dot_product = sp.sympify(0)
    for i in range(3):
        dot_product += kk[i] * xCart[i]
    dot_product /= kk_norm_factor

    # Step 5: Set initial data for uu and vv, where vv_exactsoln = \partial_t uu_exactsoln.
    # By convention, we set uu such that it is never zero. This ensures that relative error in uu is well defined.
    uu_exactsoln = sp.sin(dot_product - wavespeed * time) + 2
    vv_exactsoln = sp.diff(uu_exactsoln, time)

    return uu_exactsoln, vv_exactsoln


if __name__ == "__main__":
    import doctest
    import os
    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    for exact_WaveType in ["PlaneWave", "SphericalGaussian"]:
        ID = WaveEquation_solution_Cartesian(WaveType=exact_WaveType)
        results_dict = ve.process_dictionary_of_expressions(
            ID.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{exact_WaveType}",
            results_dict,
        )
