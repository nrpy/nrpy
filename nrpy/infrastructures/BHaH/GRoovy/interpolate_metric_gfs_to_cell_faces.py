"""
C function to interpolate metric data to cell-faces, or interfaces, from cell-centers, to 3rd order accuracy.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_interpolate_metric_gfs_to_cell_faces(
    evolving_spacetime: bool = True,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function for interpolating BSSN quantities to cell-faces.

    This function generates a C routine that interpolates metric quantities
    (lapse, shift, conformal factor, and metric tensor) from cell centers
    to cell faces along a specified flux direction using a 3rd-order stencil.
    This is typically required for reconstructing fluxes in GRHD/GRMHD solvers.

    :param evolving_spacetime: whether we're evolving the spacetime. If so, metric variables are just in_gfs, if not, auxevol.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Set up C function headers, signature, and parameters
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    cfunc_type = "void"
    desc = "Interpolate metric gridfunctions to cell faces"
    name = "interpolate_metric_gfs_to_cell_faces"
    params = "const commondata_struct *restrict commondata, const params_struct *params, const int flux_dirn, const REAL *restrict in_gfs, REAL *auxevol_gfs"

    # Step 3: Define C preprocessor macros for interpolation
    #         Defines a 4-point stencil for 3rd-order accurate interpolation
    #         to the interface at i - 1/2.
    prefunc = r"""
    // Side note: the following values could be used for cell averaged gfs:
    //     am2=-1.0/12.0, am1=7.0/12.0, a0=7.0/12.0, a1=-1.0/12.0
    // However, since the metric gfs store the grid point values instead of the cell average,
    //     the following coefficients should be used:
    //     am2 = -1/16, am1 = 9/16, a0 = 9/16, a1 = -1/16
    // This will yield the third-order-accurate face values at m-1/2,
    //      using values specified at {m-2,m-1,m,m+1}
    #define AM2 -0.0625
    #define AM1  0.5625
    #define A0   0.5625
    #define A1  -0.0625
    #define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))
    """

    # Step 4: Define the Loop Body
    #         1. Define lists of input (cell-center) and output (cell-face) gridfunctions.
    #         2. Loop over grid points and interpolate.
    loop_body = r"""
    const int metric_gfs_list[11] = {HDD00GF,
                                    HDD01GF,
                                    HDD02GF,
                                    HDD11GF,
                                    HDD12GF,
                                    HDD22GF,
                                    VETU0GF,
                                    VETU1GF,
                                    VETU2GF,
                                    ALPHAGF,
                                    CFGF};

    const int metric_gfs_face_list[11] = {H_FACEDD00GF,
                                        H_FACEDD01GF,
                                        H_FACEDD02GF,
                                        H_FACEDD11GF,
                                        H_FACEDD12GF,
                                        H_FACEDD22GF,
                                        VET_FACEU0GF,
                                        VET_FACEU1GF,
                                        VET_FACEU2GF,
                                        ALPHA_FACEGF,
                                        CF_FACEGF};

    const int num_metric_gfs = 11;

    """

    loop_body += r"""
        int in_gf,out_gf;
        REAL Qm2,Qm1,Qp0,Qp1;

        // Determine direction offsets based on flux_dirn (0=x, 1=y, 2=z)
        const int dirn0 = (flux_dirn == 0);
        const int dirn1 = (flux_dirn == 1);
        const int dirn2 = (flux_dirn == 2);

        // Loop over all metric gridfunctions to interpolate
        #pragma omp parallel
        for(int gf = 0;gf < num_metric_gfs;gf++) {
            in_gf  = metric_gfs_list[gf];
            out_gf = metric_gfs_face_list[gf];
            
            // Parallel loop over the grid, respecting ghost zones
            #pragma omp for schedule(static)
            for (int i2 = 2;i2 < Nxx_plus_2NGHOSTS2-1;i2++) {
                for (int i1 = 2;i1 < Nxx_plus_2NGHOSTS1-1;i1++) {
                    for (int i0 = 2;i0 < Nxx_plus_2NGHOSTS0-1;i0++) {
                        // Read values from stencil {i-2, i-1, i, i+1} along flux_dirn
                        Qm2 = in_gfs[IDX4(in_gf,i0-2*dirn0,i1-2*dirn1,i2-2*dirn2)];
                        Qm1 = in_gfs[IDX4(in_gf,i0-  dirn0,i1-  dirn1,i2-  dirn2)];
                        Qp0 = in_gfs[IDX4(in_gf,i0,        i1,        i2        )];
                        Qp1 = in_gfs[IDX4(in_gf,i0+  dirn0,i1+  dirn1,i2+  dirn2)];
                        
                        // Compute and store interpolated value at face i-1/2
                        auxevol_gfs[IDX4(out_gf,i0,i1,i2)] = COMPUTE_FCVAL(Qm2,Qm1,Qp0,Qp1);

                    }
                }
            }
        }
    """

    # Step 5: Handle static vs dynamic spacetime
    #         If the spacetime is not evolving, metric data resides in auxevol_gfs
    if not evolving_spacetime:
        loop_body = loop_body.replace("in_gfs", "auxevol_gfs")

    # Step 6: Register the final C function
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        prefunc=prefunc,
        body=loop_body,
    )
    return pcg.NRPyEnv()
