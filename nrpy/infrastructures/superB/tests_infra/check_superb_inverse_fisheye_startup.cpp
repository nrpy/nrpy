#include <cassert>
#include <chrono>

#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

// Compile and run with check_superb_inverse_fisheye_startup.py PROJECT_DIR.
static double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
static void freebc(bc_struct &b) {
    free(b.inner_bc_array);
    free(b.inner_bc_array_nonlocal);
    for (int j = 0; j < 3 * NGHOSTS; j++)
        free(b.pure_outer_bc_array[j]);
}
/**
 * Compare local boundary geometry, distributed CFL, and coordinate components.
 *
 * @return Zero when all comparisons pass.
 */
int main() {
    commondata_struct c = {};
    params_struct p = {};
    p.inverse_fisheye_a1 = 2.0;
    p.inverse_fisheye_R1 = 4.640939;
    p.inverse_fisheye_Rz1 = 1.5;
    p.inverse_fisheye_w1 = 0.887533;
    p.inverse_fisheye_a2 = 4.0;
    p.inverse_fisheye_R2 = 5.33048;
    p.inverse_fisheye_Rz2 = 2.5;
    p.inverse_fisheye_w2 = 1.107081;
    p.inverse_fisheye_a3 = 8.0;
    p.inverse_fisheye_R3 = 8.717449;
    p.inverse_fisheye_Rz3 = 8.717449;
    p.inverse_fisheye_w3 = 1.903712;
    p.inverse_fisheye_a4 = 16.0;
    p.inverse_fisheye_R4 = 15.437891;
    p.inverse_fisheye_Rz4 = 15.437891;
    p.inverse_fisheye_w4 = 3.802808;
    p.inverse_fisheye_a5 = 32.0;
    p.inverse_fisheye_R5 = 32.381844;
    p.inverse_fisheye_Rz5 = 32.381844;
    p.inverse_fisheye_w5 = 11.451039;
    p.inverse_fisheye_a6 = 64.0;
    p.inverse_fisheye_R6 = 57.301167;
    p.inverse_fisheye_Rz6 = 57.301167;
    p.inverse_fisheye_w6 = 19.1984;
    p.inverse_fisheye_a7 = 128.0;
    p.inverse_fisheye_R7 = 109.597767;
    p.inverse_fisheye_Rz7 = 109.597767;
    p.inverse_fisheye_w7 = 38.894618;
    p.inverse_fisheye_a8 = 256.0;
    p.inverse_fisheye_R8 = 196.835023;
    p.inverse_fisheye_Rz8 = 196.835023;
    p.inverse_fisheye_w8 = 74.89408;
    p.CoordSystem_hash = GENERALRFM_INVERSE_FISHEYEN8;
    p.grid_physical_size = 400;
    inverse_fisheye_normalize__N8(&p);
    p.xxmin0 = p.xmin;
    p.xxmin1 = p.ymin;
    p.xxmin2 = p.zmin;
    p.Nxx0 = 24;
    p.Nxx1 = 36;
    p.Nxx2 = 24;
    p.Nxx_plus_2NGHOSTS0 = p.Nxx0 + 2 * NGHOSTS;
    p.Nxx_plus_2NGHOSTS1 = p.Nxx1 + 2 * NGHOSTS;
    p.Nxx_plus_2NGHOSTS2 = p.Nxx2 + 2 * NGHOSTS;
    p.dxx0 = 2 * p.xmax / p.Nxx0;
    p.dxx1 = 2 * p.ymax / p.Nxx1;
    p.dxx2 = 2 * p.zmax / p.Nxx2;
    REAL *xx[3];
    const int nt[3] = {p.Nxx_plus_2NGHOSTS0, p.Nxx_plus_2NGHOSTS1, p.Nxx_plus_2NGHOSTS2};
    const REAL d[3] = {p.dxx0, p.dxx1, p.dxx2}, lo[3] = {p.xmin, p.ymin, p.zmin};
    for (int k = 0; k < 3; k++) {
        xx[k] = (REAL *)malloc(nt[k] * sizeof(REAL));
        for (int i = 0; i < nt[k]; i++)
            xx[k][i] = lo[k] + (i - NGHOSTS + .5) * d[k];
    }
    bc_struct global = {};
    double start = now();
    bcstruct_set_up__rfm__GeneralRFM_inverse_fisheyeN8(&c, &p, xx, &global);
    printf("baseline global boundary setup %.6f seconds\n", now() - start);
    fflush(stdout);
    assert(global.bc_info.num_inner_boundary_points == 0);
    for (int mode = 0; mode < 3; mode++) {
        c.Nchare0 = mode == 2 ? 3 : mode ? 2 : 1;
        c.Nchare1 = mode ? 3 : 1;
        c.Nchare2 = mode == 2 ? 3 : mode ? 2 : 1;
        params_struct local = p;
        local.Nxx0 = p.Nxx0 / c.Nchare0;
        local.Nxx1 = p.Nxx1 / c.Nchare1;
        local.Nxx2 = p.Nxx2 / c.Nchare2;
        local.Nxx_plus_2NGHOSTS0 = local.Nxx0 + 2 * NGHOSTS;
        local.Nxx_plus_2NGHOSTS1 = local.Nxx1 + 2 * NGHOSTS;
        local.Nxx_plus_2NGHOSTS2 = local.Nxx2 + 2 * NGHOSTS;
        double fast_seconds = 0;
        start = now();
        for (int z = 0; z < c.Nchare2; z++)
            for (int y = 0; y < c.Nchare1; y++)
                for (int x = 0; x < c.Nchare0; x++) {
                    const int ci[3] = {x, y, z};
                    bc_struct fast = {};
                    nonlocalinnerbc_struct nonlocal = {};
                    double fast_start = now();
                    bcstruct_chare_set_up__rfm__GeneralRFM_inverse_fisheyeN8(&c, &p, &local, xx, NULL, &fast, &nonlocal,
                                                                             ci);
                    fast_seconds += now() - fast_start;
                    assert(!fast.bc_info.num_inner_boundary_points && !nonlocal.tot_num_src_chares &&
                           !nonlocal.tot_num_dst_chares);
                    for (int gz = 0; gz < NGHOSTS; gz++)
                        for (int dir = 0; dir < 3; dir++) {
                            int count = 0;
                            for (int j = 0; j < global.bc_info.num_pure_outer_boundary_points[gz][dir]; j++) {
                                auto b = global.pure_outer_bc_array[3 * gz + dir][j];
                                const int64_t idx = IDX3GENERAL_64(b.i0, b.i1, b.i2, nt[0], nt[1]);
                                int owner = globalidx3pt_to_chareidx3(idx, nt[0], nt[1], local.Nxx0, local.Nxx1,
                                                                      local.Nxx2, c.Nchare0, c.Nchare1, c.Nchare2);
                                if (owner != x + c.Nchare0 * (y + c.Nchare1 * z))
                                    continue;
                                assert(count < fast.bc_info.num_pure_outer_boundary_points[gz][dir]);
                                auto a = fast.pure_outer_bc_array[3 * gz + dir][count++];
                                assert(a.i0 == b.i0 - x * local.Nxx0 && a.i1 == b.i1 - y * local.Nxx1 &&
                                       a.i2 == b.i2 - z * local.Nxx2);
                                assert(a.FACEX0 == b.FACEX0 && a.FACEX1 == b.FACEX1 && a.FACEX2 == b.FACEX2);
                                assert(isfinite(a.r) && fabs(a.r - b.r) <= 1e-11 * fmax(1., fabs(b.r)));
                                assert(isfinite(a.partial_x0_partial_r) &&
                                       fabs(a.partial_x0_partial_r - b.partial_x0_partial_r) <=
                                           1e-11 * fmax(1., fabs(b.partial_x0_partial_r)));
                                assert(isfinite(a.partial_x1_partial_r) &&
                                       fabs(a.partial_x1_partial_r - b.partial_x1_partial_r) <=
                                           1e-11 * fmax(1., fabs(b.partial_x1_partial_r)));
                                assert(isfinite(a.partial_x2_partial_r) &&
                                       fabs(a.partial_x2_partial_r - b.partial_x2_partial_r) <=
                                           1e-11 * fmax(1., fabs(b.partial_x2_partial_r)));
                                assert(isfinite(a.r_int) && fabs(a.r_int - b.r_int) <= 1e-11 * fmax(1., fabs(b.r_int)));
                                assert(isfinite(a.partial_x0_partial_r_int) &&
                                       fabs(a.partial_x0_partial_r_int - b.partial_x0_partial_r_int) <=
                                           1e-11 * fmax(1., fabs(b.partial_x0_partial_r_int)));
                                assert(isfinite(a.partial_x1_partial_r_int) &&
                                       fabs(a.partial_x1_partial_r_int - b.partial_x1_partial_r_int) <=
                                           1e-11 * fmax(1., fabs(b.partial_x1_partial_r_int)));
                                assert(isfinite(a.partial_x2_partial_r_int) &&
                                       fabs(a.partial_x2_partial_r_int - b.partial_x2_partial_r_int) <=
                                           1e-11 * fmax(1., fabs(b.partial_x2_partial_r_int)));
                            } // END LOOP: compare boundary points
                            assert(count == fast.bc_info.num_pure_outer_boundary_points[gz][dir]);
                        } // END LOOP: compare boundary faces
                    freebc(fast);
                    free(nonlocal.idx3chare_to_src_chare_id);
                    free(nonlocal.idx3chare_to_dst_chare_id);
                } // END LOOP: compare local chares
        printf("local boundary setup total %.9fs\n", fast_seconds);
        printf("PASS boundary ordering/ownership, %dx%dx%d chares (%.6fs including "
               "oracle)\n",
               c.Nchare0, c.Nchare1, c.Nchare2, now() - start);
    } // END LOOP: compare boundary decompositions
    c.CFL_FACTOR = .5;
    c.dt = 1e30;
    cfl_limited_timestep(&c, &p, xx);
    const double scanned_dt = c.dt;
    double minimum = 1e30;
    start = now();
    for (int z = 0; z < nt[2]; z++)
        for (int y = 0; y < nt[1]; y++)
            for (int x = 0; x < nt[0]; x++) {
                REAL ds;
                ds_min_single_pt__rfm__GeneralRFM_inverse_fisheyeN8(&p, xx[0][x], xx[1][y], xx[2][z], &ds);
                minimum = fmin(minimum, ds);
            } // END LOOP: scan global grid
    const double full_scan_seconds = now() - start;
    assert(scanned_dt == c.CFL_FACTOR * minimum);
    printf("PASS full CFL scan dt %.15g; exhaustive oracle %.6fs\n", scanned_dt, full_scan_seconds);
    // Local coordinate arrays are slices of the same global arrays. Their union,
    // including the outer ghosts, must produce precisely the same minimum.
    for (int mode = 0; mode < 3; mode++) {
        const int counts[3] = {mode == 0 ? 1 : mode == 1 ? 2 : 3, mode == 0 ? 1 : 3, mode == 0 ? 1 : mode == 1 ? 2 : 3};
        params_struct local = p;
        local.Nxx0 = p.Nxx0 / counts[0];
        local.Nxx1 = p.Nxx1 / counts[1];
        local.Nxx2 = p.Nxx2 / counts[2];
        local.Nxx_plus_2NGHOSTS0 = local.Nxx0 + 2 * NGHOSTS;
        local.Nxx_plus_2NGHOSTS1 = local.Nxx1 + 2 * NGHOSTS;
        local.Nxx_plus_2NGHOSTS2 = local.Nxx2 + 2 * NGHOSTS;
        double reduced_dt = 1e30;
        start = now();
        for (int z = 0; z < counts[2]; z++)
            for (int y = 0; y < counts[1]; y++)
                for (int x = 0; x < counts[0]; x++) {
                    REAL *local_xx[3] = {xx[0] + x * local.Nxx0, xx[1] + y * local.Nxx1, xx[2] + z * local.Nxx2};
                    c.dt = 1e30;
                    cfl_limited_timestep(&c, &local, local_xx);
                    reduced_dt = fmin(reduced_dt, c.dt);
                } // END LOOP: scan local chares
        assert(reduced_dt == scanned_dt);
        const double local_seconds = now() - start;
        printf("PASS distributed CFL %dx%dx%d: dt %.15g, all local scans %.6fs\n", counts[0], counts[1], counts[2],
               reduced_dt, local_seconds);
        if (mode == 2) {
            start = now();
            for (int chare = 0; chare < counts[0] * counts[1] * counts[2]; chare++) {
                c.dt = 1e30;
                cfl_limited_timestep(&c, &p, xx);
                assert(c.dt == reduced_dt);
            } // END LOOP: measure global scans
            const double repeated_seconds = now() - start;
            printf("Measured 27-chare scan work: repeated global %.6fs, local %.6fs "
                   "(%.2fx less CPU time)\n",
                   repeated_seconds, local_seconds, repeated_seconds / local_seconds);
        } // END IF: measure repeated global scan
    } // END LOOP: compare CFL decompositions

    // Compare the component cache with independent solves and Jacobians. Include
    // a changed point on every iteration, so this measures reuse across
    // components.
    double sums[2] = {0, 0};
    for (int cached = 0; cached < 2; cached++) {
        start = now();
        for (int pt = 0; pt < 1000; pt++) {
            const REAL q[3] = {.03 * pt + 1., 2., 3.};
            for (int kind = 0; kind < 3; kind++)
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < (kind ? 3 : 1); j++) {
                        if (cached)
                            sums[cached] += inverse_fisheye_component__N8(&p, q[0], q[1], q[2], kind, i, j);
                        else {
                            REAL x[3], J[3][3], K[3][3], g[3][3], dg[3][3][3], ddg[3][3][3][3];
                            assert(inverse_fisheye_solve__N8(&p, q, x) == 0);
                            if (kind)
                                inverse_fisheye_geometry__N8(&p, x, false, J, K, g, dg, ddg);
                            sums[cached] += kind == 0 ? x[i] : kind == 1 ? K[i][j] : J[i][j];
                        } // END ELSE: compare direct component
                    } // END LOOP: compare component indices
        } // END LOOP: compare coordinate points
        printf("%s component sequence %.6fs\n", cached ? "cached" : "uncached", now() - start);
    } // END LOOP: compare coordinate components
    assert(fabs(sums[0] - sums[1]) <= 1e-12 * fabs(sums[0]));
    freebc(global);
    for (int d = 0; d < 3; d++)
        free(xx[d]);
} // END FUNCTION: main
