"""
Register necessary gridfunctions to BH@H infrastructure, for GRHD evolution using GRoovy.

Consolidated into separate python file for simplicity.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

import nrpy.grid as gri  # NRPy+: Functionality for handling numerical grids
import nrpy.params as par
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities

# ==============================================================================
# Code Parameters Registration
# ==============================================================================
# Register a parameter to control the frequency of Conservative-to-Primitive (C2P)
# diagnostic output. This is useful for debugging C2P failures or monitoring
# the health of the primitive recovery routine.
_ = par.register_CodeParameter(
    "int",
    __name__,
    "C2P_diagnostics_every",
    2,
    commondata=True,
    add_to_parfile=True,
)


def register_all_grhd_gridfunctions(
    CoordSystem: str,
    evolving_spacetime: bool = True,
    evolving_temperature: bool = False,
    evolving_entropy: bool = False,
    evolving_neutrinos: bool = False,
) -> None:
    """
    Python function to register all gridfunctions needed for GRHD evolution. Best to run this first before generating C code.

    :param CoordSystem: The coordinate system (e.g., "Cartesian", "Spherical").
    :param evolving_spacetime: If True, metric variables (alpha, beta, gamma, K) are EVOL variables (BSSN). If False, they are stored as AUXEVOL variables.
    :param evolving_temperature: If True, register variables for temperature and electron fraction, for evolution with a tabulated equation of state.
    :param evolving_entropy: If True, register variables for entropy evolution.
    :param evolving_neutrinos: If True, register variables for the NRPyLeakage neutrino transport scheme (optical depths, opacities).
    """
    # ==========================================================================
    # Spacetime / Metric Gridfunctions
    # ==========================================================================
    if evolving_spacetime:
        # If the metric evolves, we register the full BSSN set (alpha, trK, cf, aDD, etc.)
        # as EVOL gridfunctions via the standard BSSN module.
        _ = BSSN_quantities[CoordSystem]

    else:
        # If the spacetime is static (Cowling approximation) or fixed background,
        # we register metric quantities as AUXEVOL (Auxiliary variables, not evolved via MoL).

        # Lapse function alpha
        _ = gri.register_gridfunctions(
            "alpha", group="AUXEVOL", gf_array_name="auxevol_gfs"
        )
        # Conformal factor cf
        _ = gri.register_gridfunctions(
            "cf", group="AUXEVOL", gf_array_name="auxevol_gfs"
        )
        # Trace of extrinsic curvature trK
        _ = gri.register_gridfunctions(
            "trK", group="AUXEVOL", gf_array_name="auxevol_gfs"
        )
        # hDD
        _ = gri.register_gridfunctions_for_single_rank2(
            "hDD",
            symmetry="sym01",
            dimension=3,
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
        )
        # aDD
        _ = gri.register_gridfunctions_for_single_rank2(
            "aDD",
            symmetry="sym01",
            dimension=3,
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
        )
        # Shift vector vetU
        _ = gri.register_gridfunctions_for_single_rank1(
            "vetU", dimension=3, group="AUXEVOL", gf_array_name="auxevol_gfs"
        )

    # ==========================================================================
    # Face-Interpolated Metric Quantities
    # ==========================================================================
    # These variables store metric values interpolated to cell faces (i+1/2),
    # which are required for computing fluxes at interfaces.
    _ = gri.register_gridfunctions(
        "alpha_face", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )
    _ = gri.register_gridfunctions(
        "cf_face", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )
    _ = gri.register_gridfunctions_for_single_rank2(
        "h_faceDD",
        symmetry="sym01",
        dimension=3,
        group="AUXEVOL",
        gf_array_name="auxevol_gfs",
    )
    _ = gri.register_gridfunctions_for_single_rank1(
        "vet_faceU", dimension=3, group="AUXEVOL", gf_array_name="auxevol_gfs"
    )

    # ==========================================================================
    # Relativistic Velocity Quantities
    # ==========================================================================
    # u4Ut: Time component of the 4-velocity (u^0)
    # u4rUt / u4lUt: Reconstructed values at Right/Left interfaces
    # Testing found we don't need to sync this across Charm++ boundaries.
    _ = gri.register_gridfunctions("u4Ut", group="AUXEVOL", gf_array_name="auxevol_gfs")
    _ = gri.register_gridfunctions(
        "u4rUt", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )
    _ = gri.register_gridfunctions(
        "u4lUt", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )

    # Reference metric rescaled three velocity rescaledvU
    # Includes Right/Left reconstructed states.
    _ = gri.register_gridfunctions_for_single_rank1(
        "rescaledvU",
        dimension=3,
        group="AUXEVOL",
        gf_array_name="auxevol_gfs",
        sync_gf_in_superB=True,  # Needs synchronization for boundaries
    )
    _ = gri.register_gridfunctions_for_single_rank1(
        "rescaledvrU", dimension=3, group="AUXEVOL", gf_array_name="auxevol_gfs"
    )
    _ = gri.register_gridfunctions_for_single_rank1(
        "rescaledvlU", dimension=3, group="AUXEVOL", gf_array_name="auxevol_gfs"
    )

    # ==========================================================================
    # Hydrodynamic Conserved Variables (EVOL)
    # ==========================================================================
    # These are the variables actually evolved by the Method of Lines (MoL).

    # Reference metric rescaled densitized conserved momentum rescaledstildeD
    _ = gri.register_gridfunctions_for_single_rank1(
        "rescaledstildeD", dimension=3, group="EVOL"
    )
    # Densitized conserved baryonic density rho_star
    _ = gri.register_gridfunctions("rho_star", group="EVOL")
    # Densitized conserved energy tau_tilde
    _ = gri.register_gridfunctions("tau_tilde", group="EVOL")

    # ==========================================================================
    # Flux Terms (Riemann Solver Outputs)
    # ==========================================================================
    # Stores the fluxes computed by the HLLE approximate Riemann solver
    # before they are used to update the conserved variables.
    _ = gri.register_gridfunctions_for_single_rank1(
        "rescaledStilde_flux_HLLD",
        dimension=3,
        group="AUXEVOL",
        gf_array_name="auxevol_gfs",
    )
    _ = gri.register_gridfunctions(
        "rho_star_HLL_flux", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )
    _ = gri.register_gridfunctions(
        "tau_tilde_HLL_flux", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )

    # ==========================================================================
    # Hydrodynamic Primitive Variables
    # ==========================================================================
    # Variables used for reconstruction and EOS calls.
    # Note: Primitives are AUXEVOL because they are computed from Conserved vars via C2P.

    # Baryon rest mass density rhob
    # Includes Right (_r) and Left (_l) reconstructed states.
    _ = gri.register_gridfunctions(
        "rhob", group="AUXEVOL", gf_array_name="auxevol_gfs", sync_gf_in_superB=True
    )
    _ = gri.register_gridfunctions(
        "rhob_r", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )
    _ = gri.register_gridfunctions(
        "rhob_l", group="AUXEVOL", gf_array_name="auxevol_gfs"
    )

    # Hydrodynamic pressure P
    _ = gri.register_gridfunctions(
        "P", group="AUXEVOL", gf_array_name="auxevol_gfs", sync_gf_in_superB=True
    )
    _ = gri.register_gridfunctions("P_r", group="AUXEVOL", gf_array_name="auxevol_gfs")
    _ = gri.register_gridfunctions("P_l", group="AUXEVOL", gf_array_name="auxevol_gfs")

    # ==========================================================================
    # Temperature and Electron Fraction (Ye) Evolution
    # ==========================================================================
    if evolving_temperature:
        # Ye: Electron fraction
        # Used for tabulated equations of state
        _ = gri.register_gridfunctions(
            "Ye", group="AUXEVOL", gf_array_name="auxevol_gfs", sync_gf_in_superB=True
        )
        _ = gri.register_gridfunctions(
            "Ye_r", group="AUXEVOL", gf_array_name="auxevol_gfs"
        )
        _ = gri.register_gridfunctions(
            "Ye_l", group="AUXEVOL", gf_array_name="auxevol_gfs"
        )
        # Flux for the electron fraction evolution
        _ = gri.register_gridfunctions(
            "Ye_star_HLL_flux", group="AUXEVOL", gf_array_name="auxevol_gfs"
        )

        # Temperature
        _ = gri.register_gridfunctions(
            "temperature",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "temperature_r", group="AUXEVOL", gf_array_name="auxevol_gfs"
        )
        _ = gri.register_gridfunctions(
            "temperature_l", group="AUXEVOL", gf_array_name="auxevol_gfs"
        )

        # Densitized conserved electron fraction Ye_star
        _ = gri.register_gridfunctions("Ye_star", group="EVOL")

    # ==========================================================================
    # Entropy Evolution
    # ==========================================================================
    if evolving_entropy:
        # Specific entropy S
        _ = gri.register_gridfunctions(
            "S", group="AUXEVOL", gf_array_name="auxevol_gfs", sync_gf_in_superB=True
        )
        _ = gri.register_gridfunctions(
            "S_r", group="AUXEVOL", gf_array_name="auxevol_gfs"
        )
        _ = gri.register_gridfunctions(
            "S_l", group="AUXEVOL", gf_array_name="auxevol_gfs"
        )
        # Flux for entropy evolution
        _ = gri.register_gridfunctions(
            "S_star_HLL_flux", group="AUXEVOL", gf_array_name="auxevol_gfs"
        )

        # Densitized conserved entropy S_star
        _ = gri.register_gridfunctions("S_star", group="EVOL")

    # ==========================================================================
    # Neutrino Leakage Quantities
    # ==========================================================================
    if evolving_neutrinos:
        # NOTE: Gridfunctions below support 3 species:
        # 1. nue  : Electron neutrinos
        # 2. anue : Electron anti-neutrinos
        # 3. nux  : Heavy lepton neutrinos (mu, tau)

        # ----------------------------------------------------------------------
        # Optical depths tau
        # ----------------------------------------------------------------------
        _ = gri.register_gridfunctions(
            "tau_0_nue",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "tau_1_nue",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "tau_0_anue",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "tau_1_anue",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "tau_0_nux",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "tau_1_nux",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )

        # ----------------------------------------------------------------------
        # Opacities kappa
        # ----------------------------------------------------------------------
        _ = gri.register_gridfunctions(
            "kappa_0_nue",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "kappa_1_nue",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "kappa_0_anue",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "kappa_1_anue",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "kappa_0_nux",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "kappa_1_nux",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )

        # ----------------------------------------------------------------------
        # Luminosities lum_
        # ----------------------------------------------------------------------
        _ = gri.register_gridfunctions(
            "lum_nue",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "lum_anue",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
        _ = gri.register_gridfunctions(
            "lum_nux",
            group="AUXEVOL",
            gf_array_name="auxevol_gfs",
            sync_gf_in_superB=True,
        )
