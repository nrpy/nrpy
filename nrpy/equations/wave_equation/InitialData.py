"""
Generating C code for plane wave initial
 data for the scalar wave equation in
 ***Cartesian*** coordinates, in up to
 *three* spatial dimensions

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Thiago Assumpcao
         assumpcaothiago **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
from typing import List, cast, Tuple
import sys
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.indexedexp as ixp
import nrpy.grid as gri

# NRPy+: Common parameters for all WaveEquation modules (defines wavespeed)
from nrpy.equations.wave_equation.CommonParams import wavespeed

# The name of this module ("InitialData") is given by __name__:
thismodule = __name__


class InitialData:
    """
    Initialize wave data based on the given parameters.

    :param WaveType: The type of wave. Can be either "PlaneWave" or "SphericalGaussian"
    :param CoordSystem: The coordinate system. Default is "Cartesian"
    :param default_k0: The default value for k0. Default is 1
    :param default_k1: The default value for k1. Default is 1
    :param default_k2: The default value for k2. Default is 1
    :param default_sigma: The default value for sigma. Default is 3
    """

    def __init__(
        self,
        WaveType: str = "PlaneWave",
        CoordSystem: str = "Cartesian",
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
            self.uu_ID, self.vv_ID = PlaneWave(
                CoordSystem=CoordSystem,
                default_k0=default_k0,
                default_k1=default_k1,
                default_k2=default_k2,
            )
        elif WaveType == "SphericalGaussian":
            self.uu_ID, self.vv_ID = SphericalGaussian(
                CoordSystem=CoordSystem, default_sigma=default_sigma
            )
        else:
            raise ValueError(
                f"Error: WaveEquation initial data WaveType={WaveType} not supported."
            )


# Set up spherically-symmetric Gaussian initial data
def SphericalGaussian(
    CoordSystem: str = "Cartesian", default_sigma: float = 3.0
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Set up initial data for a spherically-symmetric Gaussian wave.

    :param CoordSystem: The coordinate system, defaults to "Cartesian"
    :param default_sigma: The default value for sigma, defaults to 3
    """

    # Step 1: Set up Cartesian coordinates in terms of the native CoordSystem we have chosen.
    #         E.g., if CoordSystem="Cartesian", then xx_to_Cart = [xx[0],xx[1],xx[2]]
    #         or if CoordSystem="Spherical", then xx_to_Cart = [xx[0]*sp.sin(xx[1])*sp.cos(xx[2]),
    #                                                           xx[0]*sp.sin(xx[1])*sp.sin(xx[2]),
    #
    #                                                           xx[0]*sp.cos(xx[1])]
    if any("reference_metric" in key for key in sys.modules):
        # pylint: disable=C0415
        import nrpy.reference_metric as refmetric

        rfm = refmetric.reference_metric[CoordSystem]
        xx_to_Cart = rfm.xx_to_Cart
    else:
        xx_to_Cart = cast(List[sp.Symbol], ixp.declarerank1("xx"))

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
        r += xx_to_Cart[i] ** 2
    r = sp.sqrt(r)

    # Step 4: Set initial data for uu and vv, where vv_ID = \partial_t uu_ID.
    # uu_ID = (r - wavespeed*time)/r * sp.exp(- (r - wavespeed*time)**2 / (2*sigma**2) )
    # By convention, limit(uu, r->infinity) = 2. This ensures that relative error is well defined.
    uu_ID = (
        +((r - wavespeed * time) / r)
        * sp.exp(-((r - wavespeed * time) ** 2) / (2 * sigma**2))
        + ((r + wavespeed * time) / r)
        * sp.exp(-((r + wavespeed * time) ** 2) / (2 * sigma**2))
    ) + sp.sympify(
        2
    )  # Adding 2 ensures relative error is well defined (solution does not cross zero)
    vv_ID = sp.diff(uu_ID, time)

    return uu_ID, vv_ID


# Set up monochromatic plane-wave initial data
def PlaneWave(
    CoordSystem: str = "Cartesian",
    default_k0: float = 1,
    default_k1: float = 1,
    default_k2: float = 1,
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Set up initial data for a monochromatic plane wave.

    :param CoordSystem: The coordinate system, defaults to "Cartesian"
    :param default_k0: The default value for k0, defaults to 1
    :param default_k1: The default value for k1, defaults to 1
    :param default_k2: The default value for k2, defaults to 1
    """
    # Step 1: Set up Cartesian coordinates in terms of the native CoordSystem we have chosen.
    #         E.g., if CoordSystem="Cartesian", then xx_to_Cart = [xx[0],xx[1],xx[2]]
    #         or if CoordSystem="Spherical", then xx_to_Cart = [xx[0]*sp.sin(xx[1])*sp.cos(xx[2]),
    #                                                       xx[0]*sp.sin(xx[1])*sp.sin(xx[2]),
    #                                                       xx[0]*sp.cos(xx[1])]
    if any("reference_metric" in key for key in sys.modules):
        # pylint: disable=C0415
        import nrpy.reference_metric as refmetric

        rfm = refmetric.reference_metric[CoordSystem]
        xx_to_Cart = rfm.xx_to_Cart
    else:
        xx_to_Cart = cast(List[sp.Symbol], ixp.declarerank1("xx"))

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
        dot_product += kk[i] * xx_to_Cart[i]
    dot_product /= kk_norm_factor

    # Step 5: Set initial data for uu and vv, where vv_ID = \partial_t uu_ID.
    # By convention, we set uu such that it is never zero. This ensures that relative error in uu is well defined.
    uu_ID = sp.sin(dot_product - wavespeed * time) + 2
    vv_ID = sp.diff(uu_ID, time)

    return uu_ID, vv_ID


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

    for ID_type in ["PlaneWave", "SphericalGaussian"]:
        ID = InitialData(ID_type)
        results_dict = ve.process_dictionary_of_expressions(
            ID.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{ID_type}",
            results_dict,
        )
