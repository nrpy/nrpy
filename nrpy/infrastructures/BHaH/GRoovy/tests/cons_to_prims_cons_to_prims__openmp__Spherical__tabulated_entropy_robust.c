#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Recover primitive variables from conservative GRHD variables using the TabulatedEntropy GRHayLHD strategy with the robust tabulated-entropy
 * fallback.
 *
 * @param[in] commondata Common simulation data and diagnostics settings.
 * @param[in] params Grid-local runtime parameters.
 * @param[in] ghl_params GRHayL hydrodynamics parameters.
 * @param[in] eos GRHayL equation-of-state parameters.
 * @param[in] xx Reference-metric coordinate arrays.
 * @param[in,out] evol_gfs Conservative gridfunctions.
 * @param[in,out] auxevol_gfs Primitive and auxiliary gridfunctions.
 */
void cons_to_prims__rfm__Spherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                   const ghl_parameters *restrict ghl_params, const ghl_eos_parameters *restrict eos, REAL *restrict xx[3],
                                   REAL *restrict evol_gfs, REAL *restrict auxevol_gfs) {
#include "set_CodeParameters.h"
  // Step 1: Set up the interior-loop bounds and recovery diagnostics.
  const int imin = NGHOSTS;
  const int imax = Nxx_plus_2NGHOSTS0 - NGHOSTS;
  const int jmin = NGHOSTS;
  const int jmax = Nxx_plus_2NGHOSTS1 - NGHOSTS;
  const int kmin = NGHOSTS;
  const int kmax = Nxx_plus_2NGHOSTS2 - NGHOSTS;
  const int pointcount = (imax - imin) * (jmax - jmin) * (kmax - kmin);
  const char *recovery_name = "TabulatedEntropyRobust";

  int failures = 0;
  int failures_inhoriz = 0;
  int pointcount_inhoriz = 0;
  int vel_limited_ptcount = 0;
  int rho_star_fix_applied = 0;
  int backup0 = 0;
  int backup1 = 0;
  int backup2 = 0;
  int n_iter = 0;
  int pointcount_avg = 0;

#pragma omp parallel for reduction(+ : backup0, backup1, backup2, vel_limited_ptcount, rho_star_fix_applied, failures, failures_inhoriz,             \
                                       pointcount_inhoriz, n_iter, pointcount_avg) schedule(static)
  for (int k = kmin; k < kmax; k++) {
    for (int j = jmin; j < jmax; j++) {
      for (int i = imin; i < imax; i++) {
        const int index = IDX3(i, j, k);

        ghl_primitive_quantities prims;
        ghl_conservative_quantities cons, cons_orig;
        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);

        ghl_metric_quantities ADM_metric;
        basis_transform_rfm_basis_to_Cartesian(commondata, params, &prims, &cons, &ADM_metric, i, j, k, xx, auxevol_gfs, evol_gfs);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        const int in_horizon = (ADM_metric.sqrt_detgamma > ghl_params->psi6threshold);
        pointcount_inhoriz += in_horizon;
        cons_orig = cons;
        ghl_tabulated_enforce_bounds_rho_Ye_P(eos, &prims.rho, &prims.Y_e, &prims.press);
        ghl_tabulated_compute_eps_T_from_P(eos, prims.rho, prims.Y_e, prims.press, &prims.eps, &prims.temperature);

        ghl_error_codes_t error = ghl_success;
        if (cons.rho > 0.0) {
          ghl_primitive_quantities prims1, prims2, prims3, prims4;
          ghl_conservative_quantities cons_undens1, cons_undens2, cons_undens3, cons_undens4;
          prims1 = prims;
          prims2 = prims;
          prims3 = prims;
          prims4 = prims;

          ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons, &cons_undens1);
          cons_undens2 = cons_undens1;
          cons_undens3 = cons_undens1;
          cons_undens4 = cons_undens1;

          ghl_error_codes_t error1 =
              ghl_tabulated_Palenzuela1D_energy(ghl_params, eos, &ADM_metric, &metric_aux, &cons_undens1, &prims1, &diagnostics);
          ghl_error_codes_t error2 = ghl_tabulated_Newman1D_energy(ghl_params, eos, &ADM_metric, &metric_aux, &cons_undens2, &prims2, &diagnostics);
          ghl_error_codes_t error3 =
              ghl_tabulated_Palenzuela1D_entropy(ghl_params, eos, &ADM_metric, &metric_aux, &cons_undens3, &prims3, &diagnostics);
          ghl_error_codes_t error4 = ghl_tabulated_Newman1D_entropy(ghl_params, eos, &ADM_metric, &metric_aux, &cons_undens4, &prims4, &diagnostics);
          if (error1 == ghl_success && isnan(prims1.rho * prims1.press * prims1.eps * prims1.vU[0] * prims1.vU[1] * prims1.vU[2] * prims1.entropy *
                                             prims1.Y_e * prims1.temperature))
            error1 = ghl_error_c2p_singular;
          if (error2 == ghl_success && isnan(prims2.rho * prims2.press * prims2.eps * prims2.vU[0] * prims2.vU[1] * prims2.vU[2] * prims2.entropy *
                                             prims2.Y_e * prims2.temperature))
            error2 = ghl_error_c2p_singular;
          if (error3 == ghl_success && isnan(prims3.rho * prims3.press * prims3.eps * prims3.vU[0] * prims3.vU[1] * prims3.vU[2] * prims3.entropy *
                                             prims3.Y_e * prims3.temperature))
            error3 = ghl_error_c2p_singular;
          if (error4 == ghl_success && isnan(prims4.rho * prims4.press * prims4.eps * prims4.vU[0] * prims4.vU[1] * prims4.vU[2] * prims4.entropy *
                                             prims4.Y_e * prims4.temperature))
            error4 = ghl_error_c2p_singular;

          REAL err1 = 1e300;
          REAL err2 = 1e300;
          REAL err3 = 1e300;
          REAL err4 = 1e300;
          bool speed_limited_dummy = false;

          // Note: in the following lines we compare the conserved energy and electron
          // fraction between routines. We do not include entropy, since it is not
          // conserved at shocks. We ignore the diagnostic from
          // ghl_enforce_primitive_limits_and_compute_u0 here; a candidate that fails
          // primitive limiting is unlikely to be selected.

          if (error1 == ghl_success) {
            ghl_conservative_quantities cons_temp;
            ghl_enforce_primitive_limits_and_compute_u0(ghl_params, eos, &ADM_metric, &prims1, &speed_limited_dummy);
            ghl_compute_conservs(&ADM_metric, &metric_aux, &prims1, &cons_temp);
            err1 =
                fabs((cons_temp.Y_e - cons_orig.Y_e) / (cons_orig.Y_e + 1e-100)) + fabs((cons_temp.tau - cons_orig.tau) / (cons_orig.tau + 1e-100));
          } // END IF: Palenzuela energy recovery succeeded

          if (error2 == ghl_success) {
            ghl_conservative_quantities cons_temp;
            ghl_enforce_primitive_limits_and_compute_u0(ghl_params, eos, &ADM_metric, &prims2, &speed_limited_dummy);
            ghl_compute_conservs(&ADM_metric, &metric_aux, &prims2, &cons_temp);
            err2 =
                fabs((cons_temp.Y_e - cons_orig.Y_e) / (cons_orig.Y_e + 1e-100)) + fabs((cons_temp.tau - cons_orig.tau) / (cons_orig.tau + 1e-100));
          } // END IF: Newman energy recovery succeeded

          if (error3 == ghl_success) {
            ghl_conservative_quantities cons_temp;
            ghl_enforce_primitive_limits_and_compute_u0(ghl_params, eos, &ADM_metric, &prims3, &speed_limited_dummy);
            ghl_compute_conservs(&ADM_metric, &metric_aux, &prims3, &cons_temp);
            err3 =
                fabs((cons_temp.Y_e - cons_orig.Y_e) / (cons_orig.Y_e + 1e-100)) + fabs((cons_temp.tau - cons_orig.tau) / (cons_orig.tau + 1e-100));
          } // END IF: Palenzuela entropy recovery succeeded

          if (error4 == ghl_success) {
            ghl_conservative_quantities cons_temp;
            ghl_enforce_primitive_limits_and_compute_u0(ghl_params, eos, &ADM_metric, &prims4, &speed_limited_dummy);
            ghl_compute_conservs(&ADM_metric, &metric_aux, &prims4, &cons_temp);
            err4 =
                fabs((cons_temp.Y_e - cons_orig.Y_e) / (cons_orig.Y_e + 1e-100)) + fabs((cons_temp.tau - cons_orig.tau) / (cons_orig.tau + 1e-100));
          } // END IF: Newman entropy recovery succeeded

          REAL min_err = err1;
          int best_method_index = (error1 == ghl_success) ? 1 : 0;
          if (err2 < min_err) {
            min_err = err2;
            best_method_index = 2;
          } // END IF: Newman energy gives a smaller mismatch
          if (err3 < min_err) {
            min_err = err3;
            best_method_index = 3;
          } // END IF: Palenzuela entropy gives a smaller mismatch
          if (err4 < min_err) {
            min_err = err4;
            best_method_index = 4;
          } // END IF: Newman entropy gives a smaller mismatch

          if (best_method_index == 0) {
            pointcount_avg++;
            ghl_conservative_quantities cons_neigh_avg, cons_avg;
            cons_neigh_avg.rho = 0.0;
            cons_neigh_avg.tau = 0.0;
            cons_neigh_avg.SD[0] = 0.0;
            cons_neigh_avg.SD[1] = 0.0;
            cons_neigh_avg.SD[2] = 0.0;
            cons_neigh_avg.entropy = 0.0;
            cons_neigh_avg.Y_e = 0.0;

            const int iavg_min = MAX(imin, i - 1);
            const int javg_min = MAX(jmin, j - 1);
            const int kavg_min = MAX(kmin, k - 1);
            const int iavg_max = MIN(imax, i + 2);
            const int javg_max = MIN(jmax, j + 2);
            const int kavg_max = MIN(kmax, k + 2);

            int n_avg = 0;
            for (int kavg = kavg_min; kavg < kavg_max; kavg++) {
              for (int javg = javg_min; javg < javg_max; javg++) {
                for (int iavg = iavg_min; iavg < iavg_max; iavg++) {
                  const int index_avg = IDX3(iavg, javg, kavg);
                  if (index_avg == index)
                    continue;

                  ghl_conservative_quantities cons_avg_loop;
                  ghl_metric_quantities ADM_metric_avg;
                  basis_transform_rfm_basis_to_Cartesian__read_cons_only(commondata, params, &cons_avg_loop, &ADM_metric_avg, iavg, javg, kavg, xx,
                                                                         auxevol_gfs, evol_gfs);
                  cons_neigh_avg.rho += cons_avg_loop.rho;
                  cons_neigh_avg.tau += cons_avg_loop.tau;
                  cons_neigh_avg.SD[0] += cons_avg_loop.SD[0];
                  cons_neigh_avg.SD[1] += cons_avg_loop.SD[1];
                  cons_neigh_avg.SD[2] += cons_avg_loop.SD[2];
                  cons_neigh_avg.entropy += cons_avg_loop.entropy;
                  cons_neigh_avg.Y_e += cons_avg_loop.Y_e;
                  n_avg++;
                } // END LOOP: for iavg over neighboring x indices
              } // END LOOP: for javg over neighboring y indices
            } // END LOOP: for kavg over neighboring z indices

            for (int avg_weight = 1; avg_weight <= 4 && best_method_index == 0; avg_weight++) {
              const REAL w_neigh = ((REAL)avg_weight) / 4.0;
              const REAL w_self = 1.0 - w_neigh;
              const REAL inv_n_avg = 1.0 / (REAL)n_avg;

              cons_avg.rho = (w_neigh * cons_neigh_avg.rho * inv_n_avg + w_self * cons_orig.rho);
              cons_avg.tau = (w_neigh * cons_neigh_avg.tau * inv_n_avg + w_self * cons_orig.tau);
              cons_avg.SD[0] = (w_neigh * cons_neigh_avg.SD[0] * inv_n_avg + w_self * cons_orig.SD[0]);
              cons_avg.SD[1] = (w_neigh * cons_neigh_avg.SD[1] * inv_n_avg + w_self * cons_orig.SD[1]);
              cons_avg.SD[2] = (w_neigh * cons_neigh_avg.SD[2] * inv_n_avg + w_self * cons_orig.SD[2]);
              cons_avg.entropy = (w_neigh * cons_neigh_avg.entropy * inv_n_avg + w_self * cons_orig.entropy);
              cons_avg.Y_e = (w_neigh * cons_neigh_avg.Y_e * inv_n_avg + w_self * cons_orig.Y_e);

              prims1 = prims;
              prims2 = prims;
              prims3 = prims;
              prims4 = prims;
              ghl_undensitize_conservatives(ADM_metric.sqrt_detgamma, &cons_avg, &cons_undens1);
              cons_undens2 = cons_undens1;
              cons_undens3 = cons_undens1;
              cons_undens4 = cons_undens1;

              error1 = ghl_tabulated_Palenzuela1D_energy(ghl_params, eos, &ADM_metric, &metric_aux, &cons_undens1, &prims1, &diagnostics);
              error2 = ghl_tabulated_Newman1D_energy(ghl_params, eos, &ADM_metric, &metric_aux, &cons_undens2, &prims2, &diagnostics);
              error3 = ghl_tabulated_Palenzuela1D_entropy(ghl_params, eos, &ADM_metric, &metric_aux, &cons_undens3, &prims3, &diagnostics);
              error4 = ghl_tabulated_Newman1D_entropy(ghl_params, eos, &ADM_metric, &metric_aux, &cons_undens4, &prims4, &diagnostics);
              if (error1 == ghl_success && isnan(prims1.rho * prims1.press * prims1.eps * prims1.vU[0] * prims1.vU[1] * prims1.vU[2] *
                                                 prims1.entropy * prims1.Y_e * prims1.temperature))
                error1 = ghl_error_c2p_singular;
              if (error2 == ghl_success && isnan(prims2.rho * prims2.press * prims2.eps * prims2.vU[0] * prims2.vU[1] * prims2.vU[2] *
                                                 prims2.entropy * prims2.Y_e * prims2.temperature))
                error2 = ghl_error_c2p_singular;
              if (error3 == ghl_success && isnan(prims3.rho * prims3.press * prims3.eps * prims3.vU[0] * prims3.vU[1] * prims3.vU[2] *
                                                 prims3.entropy * prims3.Y_e * prims3.temperature))
                error3 = ghl_error_c2p_singular;
              if (error4 == ghl_success && isnan(prims4.rho * prims4.press * prims4.eps * prims4.vU[0] * prims4.vU[1] * prims4.vU[2] *
                                                 prims4.entropy * prims4.Y_e * prims4.temperature))
                error4 = ghl_error_c2p_singular;

              err1 = 1e300;
              err2 = 1e300;
              err3 = 1e300;
              err4 = 1e300;

              if (error1 == ghl_success) {
                ghl_conservative_quantities cons_temp;
                ghl_enforce_primitive_limits_and_compute_u0(ghl_params, eos, &ADM_metric, &prims1, &speed_limited_dummy);
                ghl_compute_conservs(&ADM_metric, &metric_aux, &prims1, &cons_temp);
                err1 =
                    fabs((cons_temp.Y_e - cons_avg.Y_e) / (cons_avg.Y_e + 1e-100)) + fabs((cons_temp.tau - cons_avg.tau) / (cons_avg.tau + 1e-100));
              } // END IF: averaged Palenzuela energy recovery succeeded

              if (error2 == ghl_success) {
                ghl_conservative_quantities cons_temp;
                ghl_enforce_primitive_limits_and_compute_u0(ghl_params, eos, &ADM_metric, &prims2, &speed_limited_dummy);
                ghl_compute_conservs(&ADM_metric, &metric_aux, &prims2, &cons_temp);
                err2 =
                    fabs((cons_temp.Y_e - cons_avg.Y_e) / (cons_avg.Y_e + 1e-100)) + fabs((cons_temp.tau - cons_avg.tau) / (cons_avg.tau + 1e-100));
              } // END IF: averaged Newman energy recovery succeeded

              if (error3 == ghl_success) {
                ghl_conservative_quantities cons_temp;
                ghl_enforce_primitive_limits_and_compute_u0(ghl_params, eos, &ADM_metric, &prims3, &speed_limited_dummy);
                ghl_compute_conservs(&ADM_metric, &metric_aux, &prims3, &cons_temp);
                err3 =
                    fabs((cons_temp.Y_e - cons_avg.Y_e) / (cons_avg.Y_e + 1e-100)) + fabs((cons_temp.tau - cons_avg.tau) / (cons_avg.tau + 1e-100));
              } // END IF: averaged Palenzuela entropy recovery succeeded

              if (error4 == ghl_success) {
                ghl_conservative_quantities cons_temp;
                ghl_enforce_primitive_limits_and_compute_u0(ghl_params, eos, &ADM_metric, &prims4, &speed_limited_dummy);
                ghl_compute_conservs(&ADM_metric, &metric_aux, &prims4, &cons_temp);
                err4 =
                    fabs((cons_temp.Y_e - cons_avg.Y_e) / (cons_avg.Y_e + 1e-100)) + fabs((cons_temp.tau - cons_avg.tau) / (cons_avg.tau + 1e-100));
              } // END IF: averaged Newman entropy recovery succeeded

              min_err = err1;
              best_method_index = (error1 == ghl_success) ? 1 : 0;
              if (err2 < min_err) {
                min_err = err2;
                best_method_index = 2;
              } // END IF: averaged Newman energy gives a smaller mismatch
              if (err3 < min_err) {
                min_err = err3;
                best_method_index = 3;
              } // END IF: averaged Palenzuela entropy gives a smaller mismatch
              if (err4 < min_err) {
                min_err = err4;
                best_method_index = 4;
              } // END IF: averaged Newman entropy gives a smaller mismatch
            } // END LOOP: for avg_weight over weighted averaging attempts
          } // END IF: all direct tabulated-entropy methods failed

          if (best_method_index == 1) {
            prims = prims1;
          } else if (best_method_index == 2) {
            prims = prims2;
          } else if (best_method_index == 3) {
            prims = prims3;
          } else if (best_method_index == 4) {
            prims = prims4;
          } else {
            ghl_set_prims_to_constant_atm(eos, &prims);
            failures++;
            failures_inhoriz += in_horizon;
          } // END IF: robust tabulated-entropy method selection completed
        } else {
          ghl_set_prims_to_constant_atm(eos, &prims);
          rho_star_fix_applied++;
        } // END IF: conservative density is positive

        ghl_enforce_primitive_limits_and_compute_u0(ghl_params, eos, &ADM_metric, &prims, &diagnostics.speed_limited);
        if (error != ghl_success) {
          ghl_set_prims_to_constant_atm(eos, &prims);
          failures++;
          failures_inhoriz += in_horizon;
          ghl_enforce_primitive_limits_and_compute_u0(ghl_params, eos, &ADM_metric, &prims, &diagnostics.speed_limited);
        } // END IF: primitive post-processing failed and atmosphere fallback was applied

        if (diagnostics.speed_limited)
          vel_limited_ptcount++;
        backup0 += diagnostics.backup[0];
        backup1 += diagnostics.backup[1];
        backup2 += diagnostics.backup[2];
        n_iter += diagnostics.n_iter;
        // Store the recovered primitive state without touching evol_gfs yet.
        // Conserved variables are rewritten in a later pass once all points have
        // finished their primitive recovery.
        const REAL xx0 = xx[0][i];
        const REAL xx1 = xx[1][j];
        const REAL xx2 = xx[2][k];

        auxevol_gfs[IDX4pt(RHOBGF, index)] = prims.rho;
        auxevol_gfs[IDX4pt(PGF, index)] = prims.press;

        auxevol_gfs[IDX4pt(YEGF, index)] = prims.Y_e;
        auxevol_gfs[IDX4pt(TEMPERATUREGF, index)] = prims.temperature;

        auxevol_gfs[IDX4pt(SGF, index)] = prims.entropy;

        const REAL vx = prims.vU[0];
        const REAL vy = prims.vU[1];
        const REAL vz = prims.vU[2];
        /*
         *  Original SymPy expressions:
         *  "[auxevol_gfs[IDX4pt(RESCALEDVU0GF, index)] = (vx*sin(xx1)*cos(xx2) + vy*sin(xx1)*sin(xx2) + vz*(sin(xx2)**2 +
         * cos(xx2)**2)*cos(xx1))/(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + sin(xx2)**2*cos(xx1)**2 + cos(xx1)**2*cos(xx2)**2)]"
         *  "[auxevol_gfs[IDX4pt(RESCALEDVU1GF, index)] = (vx*cos(xx1)*cos(xx2) + vy*sin(xx2)*cos(xx1) + vz*(-sin(xx2)**2 -
         * cos(xx2)**2)*sin(xx1))/(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + sin(xx2)**2*cos(xx1)**2 + cos(xx1)**2*cos(xx2)**2)]"
         *  "[auxevol_gfs[IDX4pt(RESCALEDVU2GF, index)] = (vx*(-sin(xx1)**2 - cos(xx1)**2)*sin(xx2) + vy*(sin(xx1)**2 +
         * cos(xx1)**2)*cos(xx2))/(sin(xx1)**2*sin(xx2)**2 + sin(xx1)**2*cos(xx2)**2 + sin(xx2)**2*cos(xx1)**2 + cos(xx1)**2*cos(xx2)**2)]"
         */
        {
          const REAL tmp0 = sin(xx1);
          const REAL tmp1 = cos(xx2);
          const REAL tmp3 = sin(xx2);
          const REAL tmp5 = cos(xx1);
          const REAL tmp6 = ((tmp3) * (tmp3));
          const REAL tmp7 = ((tmp1) * (tmp1));
          const REAL tmp9 = ((tmp0) * (tmp0));
          const REAL tmp10 = ((tmp5) * (tmp5));
          const REAL tmp11 = (1.0 / (tmp10 * tmp6 + tmp10 * tmp7 + tmp6 * tmp9 + tmp7 * tmp9));
          auxevol_gfs[IDX4pt(RESCALEDVU0GF, index)] = tmp11 * (tmp0 * tmp1 * vx + tmp0 * tmp3 * vy + tmp5 * vz * (tmp6 + tmp7));
          auxevol_gfs[IDX4pt(RESCALEDVU1GF, index)] = tmp11 * (-tmp0 * vz * (tmp6 + tmp7) + tmp1 * tmp5 * vx + tmp3 * tmp5 * vy);
          auxevol_gfs[IDX4pt(RESCALEDVU2GF, index)] = tmp11 * (tmp1 * vy * (tmp10 + tmp9) - tmp3 * vx * (tmp10 + tmp9));
        }

      } // END LOOP: for i over interior x indices
    } // END LOOP: for j over interior y indices
  } // END LOOP: for k over interior z indices

  // Step 2: Recompute conservatives from the recovered primitives in a separate pass.

  REAL error_rho_numer = 0.0;
  REAL error_tau_numer = 0.0;
  REAL error_Sx_numer = 0.0;
  REAL error_Sy_numer = 0.0;
  REAL error_Sz_numer = 0.0;

  REAL error_rho_denom = 0.0;
  REAL error_tau_denom = 0.0;
  REAL error_Sx_denom = 0.0;
  REAL error_Sy_denom = 0.0;
  REAL error_Sz_denom = 0.0;

  REAL error_ent_numer = 0.0;
  REAL error_ent_denom = 0.0;

  REAL error_Ye_numer = 0.0;
  REAL error_Ye_denom = 0.0;

#pragma omp parallel for reduction(+ : error_rho_numer, error_tau_numer, error_Sx_numer, error_Sy_numer, error_Sz_numer, error_rho_denom,            \
                                       error_tau_denom, error_Sx_denom, error_Sy_denom, error_Sz_denom, error_ent_numer, error_ent_denom,            \
                                       error_Ye_numer, error_Ye_denom) schedule(static)
  for (int k = kmin; k < kmax; k++) {
    for (int j = jmin; j < jmax; j++) {
      for (int i = imin; i < imax; i++) {
        const int index = IDX3(i, j, k);

        ghl_primitive_quantities prims;
        ghl_conservative_quantities cons, cons_orig;
        ghl_metric_quantities ADM_metric;
        basis_transform_rfm_basis_to_Cartesian(commondata, params, &prims, &cons, &ADM_metric, i, j, k, xx, auxevol_gfs, evol_gfs);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        cons_orig = cons;
        ghl_compute_conservs(&ADM_metric, &metric_aux, &prims, &cons);
        basis_transform_Cartesian_to_rfm_basis(commondata, params, &prims, &cons, i, j, k, xx, auxevol_gfs, evol_gfs);
        auxevol_gfs[IDX4pt(U4UTGF, index)] = prims.u0;

        error_rho_numer += fabs(cons.rho - cons_orig.rho);
        error_tau_numer += fabs(cons.tau - cons_orig.tau);
        error_Sx_numer += fabs(cons.SD[0] - cons_orig.SD[0]);
        error_Sy_numer += fabs(cons.SD[1] - cons_orig.SD[1]);
        error_Sz_numer += fabs(cons.SD[2] - cons_orig.SD[2]);

        error_rho_denom += cons_orig.rho;
        error_tau_denom += cons_orig.tau;
        error_Sx_denom += fabs(cons_orig.SD[0]);
        error_Sy_denom += fabs(cons_orig.SD[1]);
        error_Sz_denom += fabs(cons_orig.SD[2]);

        error_ent_numer += fabs(cons.entropy - cons_orig.entropy);
        error_ent_denom += cons_orig.entropy;

        error_Ye_numer += fabs(cons.Y_e - cons_orig.Y_e);
        error_Ye_denom += cons_orig.Y_e;

      } // END LOOP: for i over interior x indices
    } // END LOOP: for j over interior y indices
  } // END LOOP: for k over interior z indices

  // Step 3: Emit periodic conservative-to-primitive diagnostics.
  if (commondata->C2P_diagnostics_every > 0 && (commondata->nn % commondata->C2P_diagnostics_every == 0)) {
    const REAL rho_error = (error_rho_denom == 0.0) ? error_rho_numer : error_rho_numer / error_rho_denom;
    const REAL tau_error = (error_tau_denom == 0.0) ? error_tau_numer : error_tau_numer / error_tau_denom;
    const REAL Sx_error = (error_Sx_denom == 0.0) ? error_Sx_numer : error_Sx_numer / error_Sx_denom;
    const REAL Sy_error = (error_Sy_denom == 0.0) ? error_Sy_numer : error_Sy_numer / error_Sy_denom;
    const REAL Sz_error = (error_Sz_denom == 0.0) ? error_Sz_numer : error_Sz_numer / error_Sz_denom;

    const REAL ent_error = (error_ent_denom == 0.0) ? error_ent_numer : error_ent_numer / error_ent_denom;

    const REAL Ye_error = (error_Ye_denom == 0.0) ? error_Ye_numer : error_Ye_numer / error_Ye_denom;

    const REAL avg_iters = (pointcount == 0) ? 0.0 : ((REAL)n_iter) / ((REAL)pointcount);

    printf("C2P[%s]: NumPts= %d | Backups: %d %d %d | Fixes: VL= %d rho*= %d\n"
           "         Averaged= %d | Failures: %d InHoriz= %d / %d | %.2f iters/gridpt\n"
           "         Error, Sum: rho %.3e, %.3e | tau %.3e, %.3e | entropy %.3e, %.3e | Y_e %.3e, %.3e\n"

           "                     Sx %.3e, %.3e | Sy %.3e, %.3e | Sz %.3e, %.3e\n",
           recovery_name, pointcount, backup0, backup1, backup2, vel_limited_ptcount, rho_star_fix_applied,

           pointcount_avg, failures, failures_inhoriz, pointcount_inhoriz, avg_iters,

           rho_error, error_rho_denom, tau_error, error_tau_denom,

           ent_error, error_ent_denom,

           Ye_error, error_Ye_denom,

           Sx_error, error_Sx_denom, Sy_error, error_Sy_denom, Sz_error, error_Sz_denom);
  } // END IF: conservative-to-primitive diagnostics are enabled this timestep
} // END FUNCTION: cons_to_prims__rfm__Spherical
