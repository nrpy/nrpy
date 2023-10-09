"""
Register TwoPunctures code TwoPunctures_solve.c.

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


def register_CFunction_TP_solve() -> None:
    """Register C function TP_solve(), main driver function for pseudospectral solve."""
    includes = ["stdio.h", "TP_utilities.h", "TwoPunctures.h"]
    desc = "Driver routine for pseudospectral solve."
    name = "TP_solve"
    params = "ID_persist_struct *par"
    body = r"""  par->mp = par->par_m_plus;
  par->mm = par->par_m_minus;

  int const nvar = 1, n1 = par->npoints_A, n2 = par->npoints_B, n3 = par->npoints_phi;

  int imin[3], imax[3];
  int const ntotal = n1 * n2 * n3 * nvar;
#if 0
  int percent10 = 0;
#endif
  static REAL *F = NULL;
  static derivs u; //, v, cf_v;
  REAL admMass;

  if (!F) {
    REAL up, um;
    /* Solve only when called for the first time */
    F = dvector(0, ntotal - 1);
    allocate_derivs(&u, ntotal);
    allocate_derivs(&par->v, ntotal);
    allocate_derivs(&par->cf_v, ntotal);

    if (par->use_sources) {
      fprintf(stderr, "Solving puncture equation for BH-NS/NS-NS system\n");
    } else {
      fprintf(stderr, "Solving puncture equation for BH-BH system\n");
    }
    fprintf(stderr, "b = %g\n", par->par_b);

    /* initialise to 0 */
    for (int j = 0; j < ntotal; j++) {
      par->cf_v.d0[j] = 0.0;
      par->cf_v.d1[j] = 0.0;
      par->cf_v.d2[j] = 0.0;
      par->cf_v.d3[j] = 0.0;
      par->cf_v.d11[j] = 0.0;
      par->cf_v.d12[j] = 0.0;
      par->cf_v.d13[j] = 0.0;
      par->cf_v.d22[j] = 0.0;
      par->cf_v.d23[j] = 0.0;
      par->cf_v.d33[j] = 0.0;
      par->v.d0[j] = 0.0;
      par->v.d1[j] = 0.0;
      par->v.d2[j] = 0.0;
      par->v.d3[j] = 0.0;
      par->v.d11[j] = 0.0;
      par->v.d12[j] = 0.0;
      par->v.d13[j] = 0.0;
      par->v.d22[j] = 0.0;
      par->v.d23[j] = 0.0;
      par->v.d33[j] = 0.0;
    }
    /* call for external initial guess */
    // ZACH SAYS: DISABLED.
    /*
      if (par->use_external_initial_guess)
      {
      set_initial_guess(v);
      }
    */

    /* If bare masses are not given, iteratively solve for them given the
       target ADM masses target_M_plus and target_M_minus and with initial
       guesses given by par->par_m_plus and par->par_m_minus. */
    if (!(par->give_bare_mass)) {
      REAL tmp, mp_adm_err, mm_adm_err;
      char valbuf[100];

      REAL M_p = par->target_M_plus;
      REAL M_m = par->target_M_minus;

      fprintf(stderr, "Attempting to find bare masses.\n");
      fprintf(stderr, "Target ADM masses: M_p=%g and M_m=%g\n", (double)M_p, (double)M_m);
      fprintf(stderr, "ADM mass tolerance: %g\n", (double)par->adm_tol);

      /* Loop until both ADM masses are within adm_tol of their target */
      do {
        fprintf(stderr, "Bare masses: mp=%.15g, mm=%.15g\n", (double)par->mp, (double)par->mm);
        TP_Newton(*par, nvar, n1, n2, n3, par->v, par->Newton_tol, 1);

        F_of_v(*par, nvar, n1, n2, n3, par->v, F, u);

        up = PunctIntPolAtArbitPosition(*par, 0, nvar, n1, n2, n3, par->v, par->par_b, 0., 0.);
        um = PunctIntPolAtArbitPosition(*par, 0, nvar, n1, n2, n3, par->v, -par->par_b, 0., 0.);

        /* Calculate the ADM masses from the current bare mass guess */
        par->mp_adm = (1 + up) * par->mp + par->mp * par->mm / (4. * par->par_b);
        par->mm_adm = (1 + um) * par->mm + par->mp * par->mm / (4. * par->par_b);

        /* Check how far the current ADM masses are from the target */
        mp_adm_err = fabs(M_p - par->mp_adm);
        mm_adm_err = fabs(M_m - par->mm_adm);
        fprintf(stderr, "ADM mass error: M_p_err=%.15g, M_m_err=%.15g\n", (double)mp_adm_err, (double)mm_adm_err);

        /* Invert the ADM mass equation and update the bare mass guess so that
           it gives the correct target ADM masses */
        tmp = -4 * par->par_b * (1 + um + up + um * up) + sqrt(16 * par->par_b * M_m * (1 + um) * (1 + up) + pow(-M_m + M_p + 4 * par->par_b * (1 + um) * (1 + up), 2));
        par->mp = (tmp + M_p - M_m) / (2. * (1 + up));
        par->mm = (tmp - M_p + M_m) / (2. * (1 + um));

        /* Set the par->par_m_plus and par->par_m_minus parameters */
        par->par_m_plus = par->mp;
        par->par_m_minus = par->mm;

      } while ((mp_adm_err > par->adm_tol) || (mm_adm_err > par->adm_tol));

      fprintf(stderr, "Found bare masses.\n");
    }

    TP_Newton(*par, nvar, n1, n2, n3, par->v, par->Newton_tol, par->Newton_maxit);

    F_of_v(*par, nvar, n1, n2, n3, par->v, F, u);

    SpecCoef(n1, n2, n3, 0, par->v.d0, par->cf_v.d0);

    fprintf(stderr, "The two puncture masses are mp=%.17g and mm=%.17g\n", (double)par->mp, (double)par->mm);

    up = PunctIntPolAtArbitPosition(*par, 0, nvar, n1, n2, n3, par->v, par->par_b, 0., 0.);
    um = PunctIntPolAtArbitPosition(*par, 0, nvar, n1, n2, n3, par->v, -par->par_b, 0., 0.);

    /* Calculate the ADM masses from the current bare mass guess */
    par->mp_adm = (1 + up) * par->mp + par->mp * par->mm / (4. * par->par_b);
    par->mm_adm = (1 + um) * par->mm + par->mp * par->mm / (4. * par->par_b);

    fprintf(stderr, "Puncture 1 ADM mass is %g\n", (double)par->mp_adm);
    fprintf(stderr, "Puncture 2 ADM mass is %g\n", (double)par->mm_adm);

    /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*bpar->v for A=1,B=0,phi=0 */
    admMass = (par->mp + par->mm - 4 * par->par_b * PunctEvalAtArbitPosition(par->v.d0, 0, 1, 0, 0, nvar, n1, n2, n3));
    fprintf(stderr, "The total ADM mass is %g\n", (double)admMass);
    par->E = admMass;

    /*
      Run this in Mathematica (version 8 or later) with
      math -script <file>

      Needs["SymbolicC`"];
      co = Table["center_offset[" <> ToString[i] <> "]", {i, 0, 2}];
      r1 = co + {"par->par_b", 0, 0};
      r2 = co + {-"par->par_b", 0, 0};
      {p1, p2} = Table["par->par_P_" <> bh <> "[" <> ToString[i] <> "]", {bh, {"plus", "minus"}}, {i, 0, 2}];
      {s1, s2} = Table["par->par_S_" <> bh <> "[" <> ToString[i] <> "]", {bh, {"plus", "minus"}}, {i, 0, 2}];

      J = Cross[r1, p1] + Cross[r2, p2] + s1 + s2;

      JVar = Table["J" <> ToString[i], {i, 1, 3}];
      Print[OutputForm@StringReplace[
      ToCCodeString@MapThread[CAssign[#1, CExpression[#2]] &, {JVar, J}],
      "\"" -> ""]];
    */

    par->J1 = -(par->center_offset[2] * par->par_P_minus[1]) + par->center_offset[1] * par->par_P_minus[2] - par->center_offset[2] * par->par_P_plus[1] + par->center_offset[1] * par->par_P_plus[2] +
              par->par_S_minus[0] + par->par_S_plus[0];
    par->J2 = par->center_offset[2] * par->par_P_minus[0] - par->center_offset[0] * par->par_P_minus[2] + par->par_b * par->par_P_minus[2] + par->center_offset[2] * par->par_P_plus[0] -
              par->center_offset[0] * par->par_P_plus[2] - par->par_b * par->par_P_plus[2] + par->par_S_minus[1] + par->par_S_plus[1];
    par->J3 = -(par->center_offset[1] * par->par_P_minus[0]) + par->center_offset[0] * par->par_P_minus[1] - par->par_b * par->par_P_minus[1] - par->center_offset[1] * par->par_P_plus[0] +
              par->center_offset[0] * par->par_P_plus[1] + par->par_b * par->par_P_plus[1] + par->par_S_minus[2] + par->par_S_plus[2];
  }

  free_dvector(F, 0, ntotal - 1);
  free_derivs(&u, ntotal);
"""
    cfc.register_CFunction(
        subdirectory="TwoPunctures",
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )
