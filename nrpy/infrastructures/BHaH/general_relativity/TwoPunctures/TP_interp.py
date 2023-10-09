"""
Register TwoPunctures code TwoPunctures_interp.c.

TwoPunctures creates initial for two puncture black holes using a
single domain spectral method.  This method is described in
Marcus Ansorg, Bernd Brügmann, Wolfgang Tichy, "A single-domain
spectral method for black hole puncture data", PRD 70, 064011 (2004),
arXiv:gr-qc/0404056.

License: Lesser GNU Public License, version 2.0+

Authors: Marcus Ansorg
         Erik Schnetter
         Frank Löffler
         Zachariah B. Etienne (pasted original code into Python strings)
         zachetie **at** gmail **dot* com
"""
import nrpy.c_function as cfc


def register_CFunction_TP_Interp() -> None:
    """
    Register C function TP_Interp().

    Provides spectral interpolator to provide data at arbitrary point x,y,z in Cartesian basis.
    """
    includes = ["assert.h", "TP_utilities.h", "TwoPunctures.h"]
    desc = "Spectral interpolation from TwoPunctures grids onto an arbitrary point xCart[3] = {x,y,z} in the Cartesian basis."
    prefunc = f"// {desc}\n\n"
    prefunc += """
//#define STANDALONE_UNIT_TEST

/* Swap two variables */
static inline void swap(REAL *restrict const a, REAL *restrict const b) {
  REAL const t = *a;
  *a = *b;
  *b = t;
}

#undef SWAP
#define SWAP(a, b) (swap(&(a), &(b)))
"""
    name = "TP_Interp"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
     const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data"""
    body = r"""  const REAL mp = ID_persist->mp;
  const REAL mm = ID_persist->mm;

  const REAL x = xCart[0];
  const REAL y = xCart[1];
  const REAL z = xCart[2];

  int const nvar = 1, n1 = ID_persist->npoints_A, n2 = ID_persist->npoints_B, n3 = ID_persist->npoints_phi;
  int const ntotal = n1 * n2 * n3 * nvar;

  int antisymmetric_lapse, averaged_lapse, pmn_lapse, brownsville_lapse;

  enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
  enum GRID_SETUP_METHOD gsm;

  if (CCTK_EQUALS(ID_persist->grid_setup_method, "Taylor expansion")) {
    gsm = GSM_Taylor_expansion;
  } else if (CCTK_EQUALS(ID_persist->grid_setup_method, "evaluation")) {
    gsm = GSM_evaluation;
  } else {
    fprintf(stderr, "internal error\n");
    exit(1);
  }

  antisymmetric_lapse = CCTK_EQUALS(ID_persist->initial_lapse, "twopunctures-antisymmetric");
  averaged_lapse = CCTK_EQUALS(ID_persist->initial_lapse, "twopunctures-averaged");
  pmn_lapse = CCTK_EQUALS(ID_persist->initial_lapse, "psi^n");
  /*
    if (pmn_lapse)
    fprintf(stderr,"Setting initial lapse to psi^%f profile.\n",(double)ID_persist->initial_lapse_psi_exponent);
  */
  brownsville_lapse = CCTK_EQUALS(ID_persist->initial_lapse, "brownsville");
  /*
    if (brownsville_lapse)
    fprintf(stderr,"Setting initial lapse to a Brownsville-style profile with exp %f.\n",(double)ID_persist->initial_lapse_psi_exponent);
  */

  // printf("Interpolating result\n");
  /* ZACH DISABLED. ADM VARIABLE INPUT SHOULD SUFFICE.
     if (CCTK_EQUALS(metric_type, "static conformal")) {
     if (CCTK_EQUALS(conformal_storage, "factor")) {
     *conformal_state = 1;
     } else if (CCTK_EQUALS(conformal_storage, "factor+derivs")) {
     *conformal_state = 2;
     } else if (CCTK_EQUALS(conformal_storage, "factor+derivs+2nd derivs")) {
     *conformal_state = 3;
     }
     } else {
     *conformal_state = 0;
     }
  */
  const REAL conformal_state = 0;

  // MAIN OUTPUT LOOP

  /*
    char filename[100];
    sprintf(filename,"outTwoPunctures.txt");
    FILE *outascii_file = fopen(filename, "w");
    const REAL dx = 0.1;
    for(REAL x=-10.0;x<10.0;x+=dx) {
    REAL y=0.0,z=0.0;
  */
  REAL alp_out = -1.0e100; // Set to crazy value in case parameters set incorrectly
  REAL gxx_out, gxy_out, gxz_out, gyy_out, gyz_out, gzz_out;
  // REAL puncture_u_out;
  REAL Kxx_out, Kxy_out, Kxz_out, Kyy_out, Kyz_out, Kzz_out;

  // rotate xy to zx via: x->z , y->x, z->y.
  // z_dest = x_src; x_dest = y_src; y_dest = z_src
  REAL xx, yy, zz;
  {
    const REAL x_dest = x - ID_persist->center_offset[0];
    const REAL y_dest = y - ID_persist->center_offset[1];
    const REAL z_dest = z - ID_persist->center_offset[2];
    xx = z_dest;
    yy = x_dest;
    zz = y_dest;
  }

  // We implement swapping the x and z coordinates as follows.
  //    The bulk of the code that performs the actual calculations
  //    is unchanged.  This code looks only at local variables.
  //    Before the bulk --i.e., here-- we swap all x and z tensor
  //    components, and after the code --i.e., at the end of this
  //    main loop-- we swap everything back.
  if (ID_persist->swap_xz) {
    // Swap the x and z coordinates
    SWAP(xx, zz);
  }

  REAL r_plus = sqrt(pow(xx - ID_persist->par_b, 2) + pow(yy, 2) + pow(zz, 2));
  REAL r_minus = sqrt(pow(xx + ID_persist->par_b, 2) + pow(yy, 2) + pow(zz, 2));

  REAL U;
  switch (gsm) {
  case GSM_Taylor_expansion:
    U = PunctTaylorExpandAtArbitPosition(*ID_persist, 0, nvar, n1, n2, n3, ID_persist->v, xx, yy, zz);
    break;
  case GSM_evaluation:
    U = PunctIntPolAtArbitPositionFast(*ID_persist, 0, nvar, n1, n2, n3, ID_persist->cf_v, xx, yy, zz);
    break;
  default:
    assert(0);
  }
  r_plus = pow(pow(r_plus, 4) + pow(ID_persist->TP_epsilon, 4), 0.25);
  r_minus = pow(pow(r_minus, 4) + pow(ID_persist->TP_epsilon, 4), 0.25);
  if (r_plus < ID_persist->TP_Tiny)
    r_plus = ID_persist->TP_Tiny;
  if (r_minus < ID_persist->TP_Tiny)
    r_minus = ID_persist->TP_Tiny;
  REAL psi1 = 1 + 0.5 * mp / r_plus + 0.5 * mm / r_minus + U;
#define EXTEND(M, r) (M * (3. / 8 * pow(r, 4) / pow(ID_persist->TP_Extend_Radius, 5) - 5. / 4 * pow(r, 2) / pow(ID_persist->TP_Extend_Radius, 3) + 15. / 8 / ID_persist->TP_Extend_Radius))
  if (r_plus < ID_persist->TP_Extend_Radius) {
    psi1 = 1 + 0.5 * EXTEND(mp, r_plus) + 0.5 * mm / r_minus + U;
  }
  if (r_minus < ID_persist->TP_Extend_Radius) {
    psi1 = 1 + 0.5 * EXTEND(mm, r_minus) + 0.5 * mp / r_plus + U;
  }
  REAL static_psi = 1;

  REAL Aij[3][3];
  BY_Aijofxyz(*ID_persist, xx, yy, zz, Aij);

  REAL old_alp = 1.0;
  if (ID_persist->multiply_old_lapse)
    old_alp = alp_out;

  if ((conformal_state > 0) || (pmn_lapse) || (brownsville_lapse)) {

    REAL xp, yp, zp, rp, ir;
    REAL s1, s3, s5;
    REAL p, px, py, pz, pxx, pxy, pxz, pyy, pyz, pzz;
    p = 1.0;
    px = py = pz = 0.0;
    pxx = pxy = pxz = 0.0;
    pyy = pyz = pzz = 0.0;

    /* first puncture */
    xp = xx - ID_persist->par_b;
    yp = yy;
    zp = zz;
    rp = sqrt(xp * xp + yp * yp + zp * zp);
    rp = pow(pow(rp, 4) + pow(ID_persist->TP_epsilon, 4), 0.25);
    if (rp < ID_persist->TP_Tiny)
      rp = ID_persist->TP_Tiny;
    ir = 1.0 / rp;

    if (rp < ID_persist->TP_Extend_Radius) {
      ir = EXTEND(1., rp);
    }

    s1 = 0.5 * mp * ir;
    s3 = -s1 * ir * ir;
    s5 = -3.0 * s3 * ir * ir;

    p += s1;

    px += xp * s3;
    py += yp * s3;
    pz += zp * s3;

    pxx += xp * xp * s5 + s3;
    pxy += xp * yp * s5;
    pxz += xp * zp * s5;
    pyy += yp * yp * s5 + s3;
    pyz += yp * zp * s5;
    pzz += zp * zp * s5 + s3;

    /* second puncture */
    xp = xx + ID_persist->par_b;
    yp = yy;
    zp = zz;
    rp = sqrt(xp * xp + yp * yp + zp * zp);
    rp = pow(pow(rp, 4) + pow(ID_persist->TP_epsilon, 4), 0.25);
    if (rp < ID_persist->TP_Tiny)
      rp = ID_persist->TP_Tiny;
    ir = 1.0 / rp;

    if (rp < ID_persist->TP_Extend_Radius) {
      ir = EXTEND(1., rp);
    }

    s1 = 0.5 * mm * ir;
    s3 = -s1 * ir * ir;
    s5 = -3.0 * s3 * ir * ir;

    p += s1;

    px += xp * s3;
    py += yp * s3;
    pz += zp * s3;

    pxx += xp * xp * s5 + s3;
    pxy += xp * yp * s5;
    pxz += xp * zp * s5;
    pyy += yp * yp * s5 + s3;
    pyz += yp * zp * s5;
    pzz += zp * zp * s5 + s3;

    /* ZACH DISABLED:
       if (*conformal_state >= 1) {
       static_psi = p;
       psi_out = static_psi;
       }
       if (*conformal_state >= 2) {
       psix_out = px / static_psi;
       psiy_out = py / static_psi;
       psiz_out = pz / static_psi;
       }
       if (*conformal_state >= 3) {
       psixx_out = pxx / static_psi;
       psixy_out = pxy / static_psi;
       psixz_out = pxz / static_psi;
       psiyy_out = pyy / static_psi;
       psiyz_out = pyz / static_psi;
       psizz_out = pzz / static_psi;
       }
    */

    if (pmn_lapse)
      alp_out = pow(p, ID_persist->initial_lapse_psi_exponent);
    if (brownsville_lapse)
      alp_out = 2.0 / (1.0 + pow(p, ID_persist->initial_lapse_psi_exponent));

  } /* if conformal-state > 0 */

  // puncture_u_out = U;

  gxx_out = pow(psi1 / static_psi, 4);
  gxy_out = 0;
  gxz_out = 0;
  gyy_out = pow(psi1 / static_psi, 4);
  gyz_out = 0;
  gzz_out = pow(psi1 / static_psi, 4);

  Kxx_out = Aij[0][0] / pow(psi1, 2);
  Kxy_out = Aij[0][1] / pow(psi1, 2);
  Kxz_out = Aij[0][2] / pow(psi1, 2);
  Kyy_out = Aij[1][1] / pow(psi1, 2);
  Kyz_out = Aij[1][2] / pow(psi1, 2);
  Kzz_out = Aij[2][2] / pow(psi1, 2);

  if (antisymmetric_lapse || averaged_lapse) {
    alp_out = ((1.0 - 0.5 * mp / r_plus - 0.5 * mm / r_minus) / (1.0 + 0.5 * mp / r_plus + 0.5 * mm / r_minus));

    if (r_plus < ID_persist->TP_Extend_Radius) {
      alp_out = ((1.0 - 0.5 * EXTEND(mp, r_plus) - 0.5 * mm / r_minus) / (1.0 + 0.5 * EXTEND(mp, r_plus) + 0.5 * mm / r_minus));
    }
    if (r_minus < ID_persist->TP_Extend_Radius) {
      alp_out = ((1.0 - 0.5 * EXTEND(mm, r_minus) - 0.5 * mp / r_plus) / (1.0 + 0.5 * EXTEND(mp, r_minus) + 0.5 * mp / r_plus));
    }

    if (averaged_lapse) {
      alp_out = 0.5 * (1.0 + alp_out);
    }
  }
  if (ID_persist->multiply_old_lapse)
    alp_out *= old_alp;

  if (ID_persist->swap_xz) {
    /*
    // Swap the x and z components of all tensors
    if (*conformal_state >= 2) {
    SWAP (psix_out, psiz_out);
    }
    if (*conformal_state >= 3) {
    SWAP (psixx_out, psizz_out);
    SWAP (psixy_out, psiyz_out);
    }
    */
    SWAP(gxx_out, gzz_out);
    SWAP(gxy_out, gyz_out);
    SWAP(Kxx_out, Kzz_out);
    SWAP(Kxy_out, Kyz_out);
  } /* if swap_xz */

  /* ZACH DISABLED: DO NOT USE MATTER SOURCE TERMS.
     if (ID_persist->use_sources && ID_persist->rescale_sources)
     {
     Rescale_Sources(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2],
     x, y, z,
     (conformal_state > 0) ? psi : NULL,
     gxx, gyy, gzz,
     gxy, gxz, gyz);
     }
     fprintf(outascii_file, "%e %e %e | gij's: %e %e %e %e %e %e | Kij's: %e %e %e %e %e %e | alp: %e | TwoPunctures OUTPUT\n",xx,yy,zz,
     gxx_out,gxy_out,gxz_out,gyy_out,gyz_out,gzz_out,
     Kxx_out,Kxy_out,Kxz_out,Kyy_out,Kyz_out,Kzz_out,
     alp_out);
     }
     fclose(outascii_file);
  */
  // New z = old x
  // New x = old y
  // New y = old z
  initial_data->gammaSphorCartDD22 = gxx_out;
  initial_data->gammaSphorCartDD02 = gxy_out;  // Technically gammaSphorCartDD20 = gammaSphorCartDD02
  initial_data->gammaSphorCartDD12 = gxz_out;  // Technically gammaSphorCartDD21 = gammaSphorCartDD12
  initial_data->gammaSphorCartDD00 = gyy_out;
  initial_data->gammaSphorCartDD01 = gyz_out;
  initial_data->gammaSphorCartDD11 = gzz_out;

  initial_data->KSphorCartDD22 = Kxx_out;
  initial_data->KSphorCartDD02 = Kxy_out;  // Technically KSphorCartDD20 = KSphorCartDD02
  initial_data->KSphorCartDD12 = Kxz_out;  // Technically KSphorCartDD21 = KSphorCartDD12
  initial_data->KSphorCartDD00 = Kyy_out;
  initial_data->KSphorCartDD01 = Kyz_out;
  initial_data->KSphorCartDD11 = Kzz_out;

  initial_data->alpha = alp_out;
  initial_data->betaSphorCartU0 = 0;
  initial_data->betaSphorCartU1 = 0;
  initial_data->betaSphorCartU2 = 0;
  initial_data->BSphorCartU0 = 0;
  initial_data->BSphorCartU1 = 0;
  initial_data->BSphorCartU2 = 0;

  ///////////////////////////////////
#ifdef STANDALONE_UNIT_TEST

  int main() {

    ID_persist_struct par;
    char ID_type[100];
    snprintf(ID_type, 100, "QC0");
    TwoPunctures_initialize_parameters(ID_type, &par);

    TwoPunctures_solve(&par);

    TwoPunctures_Interp(&par);

    free_derivs(&par.v, par.npoints_A * par.npoints_B * par.npoints_phi * 1);
    free_derivs(&par.cf_v, par.npoints_A * par.npoints_B * par.npoints_phi * 1);

    return 0;
  }

#endif
"""

    cfc.register_CFunction(
        subdirectory="TwoPunctures",
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )
