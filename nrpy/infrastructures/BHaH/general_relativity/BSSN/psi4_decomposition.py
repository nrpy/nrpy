"""
Generate Psi4 decomposition on spherical like grids using spin-weighted spherical harmonics.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc


def register_CFunction_psi4_spinweightm2_decomposition_on_sphlike_grids() -> None:
    """Register C function for decomposing psi4 into spin-weighted spherical harmonics."""
    prefunc = r"""
static void lowlevel_decompose_psi4_into_swm2_modes(const int Nxx_plus_2NGHOSTS1,const int Nxx_plus_2NGHOSTS2,
                                                    const REAL dxx1, const REAL dxx2,
                                                    const int swm2sh_maximum_l_mode_to_compute, const REAL curr_time, const REAL R_ext,
                                                    const REAL *restrict th_array, const REAL *restrict sinth_array, const REAL *restrict ph_array,
                                                    const REAL *restrict psi4r_at_R_ext, const REAL *restrict psi4i_at_R_ext) {
  char filename[100];  FILE *outpsi4_l_m;
  // Output header at t=0:
  if(curr_time==0) {
    for(int l=2;l<=swm2sh_maximum_l_mode_to_compute;l++) {
      sprintf(filename,"Rpsi4_l%d-r%06.1f.txt",l,(double)R_ext);
      outpsi4_l_m = fopen(filename, "w");
      fprintf(outpsi4_l_m, "# column 1: t-R_ext = [retarded time]\n");
      int col=2;
      for(int m=-l;m<=l;m++) {
        fprintf(outpsi4_l_m, "# column %d: Re(psi4_{l=%d,m=%d}) * R_ext\n", col,l,m);  col++;
        fprintf(outpsi4_l_m, "# column %d: Im(psi4_{l=%d,m=%d}) * R_ext\n", col,l,m);  col++;
      }
      fclose(outpsi4_l_m);
    }
  }

  // Output one file per l mode; each column represents a unique complex component of l,m
  for(int l=2;l<=swm2sh_maximum_l_mode_to_compute;l++) {
    sprintf(filename,"Rpsi4_l%d-r%06.1f.txt",l,(double)R_ext);
    outpsi4_l_m = fopen(filename, "a");
    char oneline[10000];
    sprintf(oneline, "%e", (double)(curr_time - R_ext));
    for(int m=-l;m<=l;m++) {
      // Parallelize the integration loop:
      REAL psi4r_l_m = 0.0;
      REAL psi4i_l_m = 0.0;
#pragma omp parallel for reduction(+:psi4r_l_m,psi4i_l_m)
      for(int i1=0;i1<Nxx_plus_2NGHOSTS1-2*NGHOSTS;i1++) {
        const REAL th    = th_array[i1];
        const REAL sinth = sinth_array[i1];
        for(int i2=0;i2<Nxx_plus_2NGHOSTS2-2*NGHOSTS;i2++) {
          const REAL ph = ph_array[i2];
          // Construct integrand for psi4 spin-weight s=-2 spherical harmonic
          REAL ReY_sm2_l_m,ImY_sm2_l_m;
          spin_weight_minus2_sph_harmonics(l,m, th,ph,  &ReY_sm2_l_m,&ImY_sm2_l_m);

          const int idx2d = i1*(Nxx_plus_2NGHOSTS2-2*NGHOSTS)+i2;
          const REAL a = psi4r_at_R_ext[idx2d];
          const REAL b = psi4i_at_R_ext[idx2d];
          const REAL c = ReY_sm2_l_m;
          const REAL d = ImY_sm2_l_m;
          psi4r_l_m += (a*c + b*d) * dxx2  * sinth*dxx1;
          psi4i_l_m += (b*c - a*d) * dxx2  * sinth*dxx1;
        }
      }
      sprintf(oneline + strlen(oneline), " %.15e %.15e", (double)(R_ext*psi4r_l_m), (double)(R_ext*psi4i_l_m));
    }
    fprintf(outpsi4_l_m, "%s\n", oneline);
    fclose(outpsi4_l_m);
  }
}
"""

    desc = "Decompose psi4 across all l,m modes from l=2 up to and including L_MAX (global variable)"
    name = "psi4_spinweightm2_decomposition_on_sphlike_grids"
    params = r"""const commondata_struct *restrict commondata,const params_struct *restrict params,
    REAL *restrict diagnostic_output_gfs,
    const REAL *restrict list_of_R_exts, const int num_of_R_exts,
    const int psi4_spinweightm2_sph_harmonics_max_l, REAL *restrict xx[3]"""

    body = r"""  // Step 1: Allocate memory for 2D arrays used to store psi4, theta, sin(theta), and phi.
  const int sizeof_2Darray = sizeof(REAL)*(Nxx_plus_2NGHOSTS1-2*NGHOSTS)*(Nxx_plus_2NGHOSTS2-2*NGHOSTS);
  REAL *restrict psi4r_at_R_ext = (REAL *restrict)malloc(sizeof_2Darray);
  REAL *restrict psi4i_at_R_ext = (REAL *restrict)malloc(sizeof_2Darray);
  //         ... also store theta, sin(theta), and phi to corresponding 1D arrays.
  REAL *restrict sinth_array = (REAL *restrict)malloc(sizeof(REAL)*(Nxx_plus_2NGHOSTS1-2*NGHOSTS));
  REAL *restrict th_array    = (REAL *restrict)malloc(sizeof(REAL)*(Nxx_plus_2NGHOSTS1-2*NGHOSTS));
  REAL *restrict ph_array    = (REAL *restrict)malloc(sizeof(REAL)*(Nxx_plus_2NGHOSTS2-2*NGHOSTS));

  const int NinterpGHOSTS = MIN(2, NGHOSTS-1);
  const int N0 = 2*NinterpGHOSTS; // Interp stencil is 2*NinterpGHOSTS+1 in size;
  //                                 reaches NinterpGHOSTS to the left & right of
  //                                 central point.
  const REAL pow_dxx0__N0 = pow(params->dxx0, N0);

  // Step 2: Loop over all extraction indices:
  for(int which_R_ext=0;which_R_ext<num_of_R_exts;which_R_ext++) {
    // Step 2.a: Set the extraction radius R_ext based on the radial index R_ext_idx
    const REAL R_ext = list_of_R_exts[which_R_ext];
    const REAL xCart_R_ext[3] = { R_ext, 0.0, 0.0 }; // just put a point on the x-axis.

    int Cart_to_i0i1i2[3]; REAL closest_xx[3];
    Cart_to_xx_and_nearest_i0i1i2(params, xCart_R_ext, closest_xx, Cart_to_i0i1i2);

    const int closest_i0=Cart_to_i0i1i2[0];

    // We want a src grid point inside the source grid (duh) with
    //  mask=+0, as all mask=+0 points will have at least
    //  NGHOSTS>=NinterpGHOSTS of filled neighbor pts.
    if(IS_IN_GRID_INTERIOR(Cart_to_i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NGHOSTS)) {

      // Step 2.a.i: Set radial interpolation coefficients for r=R_ext.
      REAL l0i__times__w0i_inv[2*NinterpGHOSTS+1];
      {
        for(int i=0;i<=N0;i++) {
          REAL prod_numer_i = 1.0;
          int  prod_denom_i = 1;
          for(int l=0;  l<i;  l++) { prod_denom_i *= i-l; prod_numer_i *= closest_xx[0] - xx[0][closest_i0 - N0/2 + l]; }
          for(int l=i+1;l<=N0;l++) { prod_denom_i *= i-l; prod_numer_i *= closest_xx[0] - xx[0][closest_i0 - N0/2 + l]; }
          l0i__times__w0i_inv[i] = prod_numer_i / ( pow_dxx0__N0 * (REAL)prod_denom_i );
        }
      }

      // Step 2.b: Compute psi_4 at this extraction radius and store to a local 2D array.
#pragma omp parallel for
      for(int i1=NGHOSTS;i1<Nxx_plus_2NGHOSTS1-NGHOSTS;i1++) {
        th_array[i1-NGHOSTS]    =     xx[1][i1];
        sinth_array[i1-NGHOSTS] = sin(xx[1][i1]);
        for(int i2=NGHOSTS;i2<Nxx_plus_2NGHOSTS2-NGHOSTS;i2++) {
          ph_array[i2-NGHOSTS] = xx[2][i2];

          REAL sum_psi4r=0;
          REAL sum_psi4i=0;
          // Perform radial interpolation to get psi4 at desired extraction radius R_ext.
          for(int i=0;i<=N0;i++) {
            // psi4r and psi4i in fixed frame have been stored to diagnostic_output_gfs.
            //  Here we interpolate to specific radius.
            sum_psi4r += diagnostic_output_gfs[IDX4(PSI4_REGF, i + closest_i0 - N0/2, i1, i2)] * l0i__times__w0i_inv[i];
            sum_psi4i += diagnostic_output_gfs[IDX4(PSI4_IMGF, i + closest_i0 - N0/2, i1, i2)] * l0i__times__w0i_inv[i];
          }
          // Store result to "2D" array (actually 1D array with 2D storage):
          const int idx2d = (i1-NGHOSTS)*(Nxx_plus_2NGHOSTS2-2*NGHOSTS)+(i2-NGHOSTS);
          psi4r_at_R_ext[idx2d] = sum_psi4r;
          psi4i_at_R_ext[idx2d] = sum_psi4i;
        }
      }
      // Step 3: Perform integrations across all l,m modes from l=2 up to and including L_MAX (global variable):
      lowlevel_decompose_psi4_into_swm2_modes(Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2,
                                              dxx1,dxx2, swm2sh_maximum_l_mode_to_compute,
                                              time, R_ext, th_array, sinth_array, ph_array,
                                              psi4r_at_R_ext,psi4i_at_R_ext);
    }
  }
  // Step 4: Free all allocated memory:
  free(psi4r_at_R_ext); free(psi4i_at_R_ext);
  free(sinth_array); free(th_array); free(ph_array);
"""

    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        prefunc=prefunc,
        desc=desc,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
