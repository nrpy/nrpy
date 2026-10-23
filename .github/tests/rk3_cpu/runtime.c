/* Numerical CI for generated projected RK3, with scalar and SIMD CPU kernels.
 */
#include "BHaH_function_prototypes.h"

_Static_assert(NUM_EVOL_GFS == 24, "This CPU comparison evolves BSSN.");
_Static_assert(sizeof(MoL_gridfunctions_struct) == 4 * sizeof(REAL *),
               "RK3 has y_n, two intermediate arrays, and auxevol.");

static void make_boundary(const params_struct *params, REAL *xx[3],
                          bc_struct *bc, int mixed) {
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  int lo[3] = {NGHOSTS, NGHOSTS, NGHOSTS};
  int hi[3] = {Nxx_plus_2NGHOSTS0 - NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
               Nxx_plus_2NGHOSTS2 - NGHOSTS};
  const int dims[3] = {Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1,
                       Nxx_plus_2NGHOSTS2};
  memset(bc, 0, sizeof(*bc));
  for (int layer = 0; layer < NGHOSTS; layer++)
    for (int dir = 0; dir < 3; dir++) {
      int count = 0;
      outerpt_bc_struct *array =
          calloc(2 * dims[(dir + 1) % 3] * dims[(dir + 2) % 3], sizeof(*array));
      bc->pure_outer_bc_array[3 * layer + dir] = array;
      for (int face = 0; face < 2; face++) {
        int start[3] = {lo[0], lo[1], lo[2]}, end[3] = {hi[0], hi[1], hi[2]};
        start[dir] = face ? hi[dir] : lo[dir] - 1;
        end[dir] = start[dir] + 1;
        for (int k = start[2]; k < end[2]; k++)
          for (int j = start[1]; j < end[1]; j++)
            for (int i = start[0]; i < end[0]; i++) {
              if (mixed &&
                  ((i == 0 && j == 0 && k == 0) ||
                   (i == dims[0] - 1 && j == dims[1] - 1 && k == dims[2] - 1)))
                continue;
              outerpt_bc_struct *b = &array[count++];
              int ijk[3] = {i, j, k}, src[3] = {i, j, k};
              src[dir] += face ? -1 : 1;
              b->i0 = i;
              b->i1 = j;
              b->i2 = k;
              b->FACEX0 = dir == 0 ? (face ? -1 : 1) : 0;
              b->FACEX1 = dir == 1 ? (face ? -1 : 1) : 0;
              b->FACEX2 = dir == 2 ? (face ? -1 : 1) : 0;
              REAL x[3] = {xx[0][ijk[0]], xx[1][ijk[1]], xx[2][ijk[2]]};
              REAL y[3] = {xx[0][src[0]], xx[1][src[1]], xx[2][src[2]]};
              b->r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
              b->r_int = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
              b->partial_x0_partial_r = x[0] / b->r;
              b->partial_x1_partial_r = x[1] / b->r;
              b->partial_x2_partial_r = x[2] / b->r;
              b->partial_x0_partial_r_int = y[0] / b->r_int;
              b->partial_x1_partial_r_int = y[1] / b->r_int;
              b->partial_x2_partial_r_int = y[2] / b->r_int;
            }
      }
      lo[dir]--;
      hi[dir]++;
      bc->bc_info.num_pure_outer_boundary_points[layer][dir] = count;
    }
  if (mixed) {
    bc->bc_info.num_inner_boundary_points = 2;
    bc->inner_bc_array = calloc(2, sizeof(*bc->inner_bc_array));
    bc->inner_bc_array[0].dstpt = IDX3(0, 0, 0);
    bc->inner_bc_array[0].srcpt = IDX3(4, 4, 4);
    bc->inner_bc_array[1].dstpt = IDX3(dims[0] - 1, dims[1] - 1, dims[2] - 1);
    bc->inner_bc_array[1].srcpt = IDX3(2, 2, 0);
    for (int i = 0; i < 10; i++) {
      bc->inner_bc_array[0].parity[i] = i % 2 ? -1 : 1;
      bc->inner_bc_array[1].parity[i] = i % 3 ? 1 : -1;
    }
  }
}

/* Count actual linked evaluations; dense-reference calls are excluded. */
static int counting, ricci_calls, rhs_calls, projection_calls;
static REAL rhs_times[3];
static void count_rhs(const commondata_struct *commondata) {
  if (counting) {
    if (rhs_calls < 3)
      rhs_times[rhs_calls] = commondata->time;
    rhs_calls++;
  }
}
#define REAL_FUNCTION(name) extern __typeof__(name) __real_##name
REAL_FUNCTION(Ricci_eval__rfm__Cartesian);
REAL_FUNCTION(Ricci_eval_tile__rfm__Cartesian);
REAL_FUNCTION(rhs_eval__rfm__Cartesian);
REAL_FUNCTION(rhs_eval_tile__rfm__Cartesian);
REAL_FUNCTION(rhs_eval_accumulate__rfm__Cartesian);
REAL_FUNCTION(enforce_detgbar_equals_detghat_trAzero__rfm__Cartesian);

void __wrap_Ricci_eval__rfm__Cartesian(const params_struct *p,
                                       const rfm_struct *r, const REAL *in,
                                       REAL *out) {
  if (counting)
    ricci_calls++;
  __real_Ricci_eval__rfm__Cartesian(p, r, in, out);
}
void __wrap_Ricci_eval_tile__rfm__Cartesian(const params_struct *p,
                                            const rfm_struct *r, const REAL *in,
                                            REAL *out, int lo1, int hi1,
                                            int lo2, int hi2) {
  if (counting)
    ricci_calls++;
  __real_Ricci_eval_tile__rfm__Cartesian(p, r, in, out, lo1, hi1, lo2, hi2);
}
void __wrap_rhs_eval__rfm__Cartesian(const commondata_struct *c,
                                     const params_struct *p,
                                     const rfm_struct *r, const REAL *aux,
                                     const REAL *in, REAL *out) {
  count_rhs(c);
  __real_rhs_eval__rfm__Cartesian(c, p, r, aux, in, out);
}
void __wrap_rhs_eval_tile__rfm__Cartesian(const commondata_struct *c,
                                          const params_struct *p,
                                          const rfm_struct *r, const REAL *aux,
                                          const REAL *in, REAL *out, int lo1,
                                          int hi1, int lo2, int hi2) {
  count_rhs(c);
  __real_rhs_eval_tile__rfm__Cartesian(c, p, r, aux, in, out, lo1, hi1, lo2,
                                       hi2);
}
void __wrap_rhs_eval_accumulate__rfm__Cartesian(
    const commondata_struct *c, const params_struct *p, const rfm_struct *r,
    const REAL *aux, const REAL *in, REAL *stage, REAL *accum, REAL a, REAL b) {
  count_rhs(c);
  __real_rhs_eval_accumulate__rfm__Cartesian(c, p, r, aux, in, stage, accum, a,
                                             b);
}
void __wrap_enforce_detgbar_equals_detghat_trAzero__rfm__Cartesian(
    const params_struct *p, const rfm_struct *r, REAL *in, const REAL *aux) {
  if (counting)
    projection_calls++;
  __real_enforce_detgbar_equals_detghat_trAzero__rfm__Cartesian(p, r, in, aux);
}

static int is_interior(size_t p, const params_struct *params) {
  const int n0 = params->Nxx_plus_2NGHOSTS0, n1 = params->Nxx_plus_2NGHOSTS1;
  int i = p % n0, j = (p / n0) % n1, k = p / (n0 * n1);
  return i >= NGHOSTS && i < n0 - NGHOSTS && j >= NGHOSTS && j < n1 - NGHOSTS &&
         k >= NGHOSTS && k < params->Nxx_plus_2NGHOSTS2 - NGHOSTS;
}

static void prepare_state(commondata_struct *common, griddata_struct *grid,
                          REAL *state) {
  if (strcmp(common->outer_bc_type, "extrapolation") == 0)
    apply_bcs_outerextrap_and_inner(common, &grid->params, &grid->bcstruct,
                                    state);
  else
    apply_bcs_inner_only(common, &grid->params, &grid->bcstruct, state);
  enforce_detgbar_equals_detghat_trAzero__rfm__Cartesian(
      &grid->params, grid->rfmstruct, state, grid->gridfuncs.auxevol_gfs);
}

static void raw_rhs(commondata_struct *common, griddata_struct *grid,
                    const REAL *state, REAL *rhs, size_t size) {
  memset(rhs, 0, size * sizeof(REAL));
  Ricci_eval__rfm__Cartesian(&grid->params, grid->rfmstruct, state,
                             grid->gridfuncs.auxevol_gfs);
  rhs_eval__rfm__Cartesian(common, &grid->params, grid->rfmstruct,
                           grid->gridfuncs.auxevol_gfs, state, rhs);
  if (strcmp(common->outer_bc_type, "radiation") == 0)
    apply_bcs_outerradiation_and_inner__rfm__Cartesian(
        common, &grid->params, &grid->bcstruct, grid->xx,
        gridfunctions_wavespeed, gridfunctions_f_infinity, state, rhs);
  else {
    const size_t ntot = size / NUM_EVOL_GFS;
    for (size_t p = 0; p < size; p++)
      if (!is_interior(p % ntot, &grid->params))
        rhs[p] = 0;
  }
}

/* Independent dense tableau: retain all three derivatives across projections.
 */
static void dense_step(commondata_struct *common, griddata_struct *grid,
                       REAL *y, int heun, size_t size) {
  REAL *stage = calloc(size, sizeof(REAL)), *k1 = calloc(size, sizeof(REAL));
  REAL *k2 = calloc(size, sizeof(REAL)), *k3 = calloc(size, sizeof(REAL));
  if (!stage || !k1 || !k2 || !k3)
    exit(2);
  const REAL h = common->dt, t = common->time;
  const REAL a21 = heun ? 1.0 / 3 : .5, a32 = heun ? 2.0 / 3 : .75;
  const REAL b1 = heun ? .25 : 2.0 / 9, b2 = heun ? 0 : 1.0 / 3,
             b3 = heun ? .75 : 4.0 / 9;
  raw_rhs(common, grid, y, k1, size);
  for (size_t p = 0; p < size; p++)
    stage[p] = y[p] + a21 * h * k1[p];
  prepare_state(common, grid, stage);
  common->time = t + a21 * h;
  raw_rhs(common, grid, stage, k2, size);
  for (size_t p = 0; p < size; p++)
    stage[p] = y[p] + a32 * h * k2[p];
  prepare_state(common, grid, stage);
  common->time = t + a32 * h;
  raw_rhs(common, grid, stage, k3, size);
  for (size_t p = 0; p < size; p++)
    y[p] += h * (b1 * k1[p] + b2 * k2[p] + b3 * k3[p]);
  prepare_state(common, grid, y);
  common->time = common->t_0 + (common->nn - common->nn_0 + 1) * h;
  common->nn++;
  free(stage);
  free(k1);
  free(k2);
  free(k3);
}

static REAL compare(const REAL *actual, const REAL *expected, size_t size) {
  REAL maxerr = 0;
  for (size_t p = 0; p < size; p++) {
    REAL err = fabs(actual[p] - expected[p]);
    if (!isfinite(actual[p]) || !isfinite(expected[p]) ||
        err > 5e-12 * (1 + fabs(expected[p]))) {
      fprintf(stderr, "Numerical mismatch at %zu: %.17g versus %.17g\n", p,
              actual[p], expected[p]);
      exit(1);
    }
    if (err > maxerr)
      maxerr = err;
  }
  return maxerr;
}

static void check_accumulation(commondata_struct *common, griddata_struct *grid,
                               const REAL *input, size_t size, REAL *maxerr) {
  REAL *saved = malloc(size * sizeof(REAL)), *raw = calloc(size, sizeof(REAL));
  REAL *stage = malloc(size * sizeof(REAL)),
       *accum = malloc(size * sizeof(REAL));
  REAL *original = malloc(size * sizeof(REAL)),
       *oldaccum = malloc(size * sizeof(REAL)),
       *expected = malloc(size * sizeof(REAL));
  if (!saved || !raw || !stage || !accum || !original || !oldaccum || !expected)
    exit(2);
  memcpy(saved, input, size * sizeof(REAL));
  for (size_t p = 0; p < size; p++) {
    original[p] = .15 * sin(.019 * p) + .2 * cos(.011 * p);
    oldaccum[p] = -.31 * cos(.017 * p) + .13 * sin(.023 * p);
  }
  raw_rhs(common, grid, input, raw, size);
  const REAL coefficients[][2] = {
      {.02775, .037 / 3}, {.037 * 2 / 3, 0}, {-.013, .021}, {0, 0}};
  const int radiation = strcmp(common->outer_bc_type, "radiation") == 0;
  for (int c = 0; c < 4; c++) {
    REAL a = coefficients[c][0], b = coefficients[c][1];
    memcpy(stage, original, size * sizeof(REAL));
    memcpy(accum, oldaccum, size * sizeof(REAL));
    if (radiation)
      apply_bcs_outerradiation_and_inner_accumulate__rfm__Cartesian(
          common, &grid->params, &grid->bcstruct, grid->xx,
          gridfunctions_wavespeed, gridfunctions_f_infinity, input, stage,
          accum, a, b, 1);
    rhs_eval_accumulate__rfm__Cartesian(common, &grid->params, grid->rfmstruct,
                                        grid->gridfuncs.auxevol_gfs, input,
                                        stage, accum, a, b);
    if (radiation)
      apply_bcs_outerradiation_and_inner_accumulate__rfm__Cartesian(
          common, &grid->params, &grid->bcstruct, grid->xx,
          gridfunctions_wavespeed, gridfunctions_f_infinity, input, stage,
          accum, a, b, 0);
    for (size_t p = 0; p < size; p++)
      expected[p] = original[p] + a * raw[p];
    REAL err = compare(stage, expected, size);
    if (err > *maxerr)
      *maxerr = err;
    for (size_t p = 0; p < size; p++)
      expected[p] = oldaccum[p] + b * raw[p];
    err = compare(accum, expected, size);
    if (err > *maxerr)
      *maxerr = err;
    if (memcmp(input, saved, size * sizeof(REAL)) ||
        (a == 0 && memcmp(stage, original, size * sizeof(REAL))) ||
        (b == 0 && memcmp(accum, oldaccum, size * sizeof(REAL)))) {
      fputs("Input or zero-weight destination changed\n", stderr);
      exit(1);
    }
  }
  free(saved);
  free(raw);
  free(stage);
  free(accum);
  free(original);
  free(oldaccum);
  free(expected);
}

int main(void) {
  REAL max_stage_error = 0, max_step_error = 0;
  int stage_cases = 0, step_cases = 0;
  for (int n0 = 16; n0 <= 18; n0 += 2)
    for (int mixed = 0; mixed < 2; mixed++) {
      const int Nxx_plus_2NGHOSTS0 = n0 + 2 * NGHOSTS, Nxx_plus_2NGHOSTS1 = 12,
                Nxx_plus_2NGHOSTS2 = 11;
      const int ntot =
          Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
      const size_t size = NUM_EVOL_GFS * (size_t)ntot;
      griddata_struct grid = {0};
      grid.params.Nxx_plus_2NGHOSTS0 = Nxx_plus_2NGHOSTS0;
      grid.params.Nxx_plus_2NGHOSTS1 = Nxx_plus_2NGHOSTS1;
      grid.params.Nxx_plus_2NGHOSTS2 = Nxx_plus_2NGHOSTS2;
      grid.params.invdxx0 = grid.params.invdxx1 = grid.params.invdxx2 = 5;
      rfm_struct rfm = {0};
      grid.rfmstruct = &rfm;
      commondata_struct initial = {.NUMGRIDS = 1,
                                   .KreissOliger_strength_gauge = .3,
                                   .KreissOliger_strength_nongauge = .2,
                                   .eta = 1.2,
                                   .dt = .01,
                                   .t_0 = .125,
                                   .time = .125};
      int dims[3] = {Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1,
                     Nxx_plus_2NGHOSTS2};
      for (int d = 0; d < 3; d++) {
        grid.xx[d] = calloc(dims[d], sizeof(REAL));
        if (!grid.xx[d])
          return 2;
        for (int i = 0; i < dims[d]; i++)
          grid.xx[d][i] = 1.1 + .2 * i;
      }
      make_boundary(&grid.params, grid.xx, &grid.bcstruct, mixed);
      grid.gridfuncs.auxevol_gfs =
          calloc(NUM_AUXEVOL_GFS * (size_t)ntot, sizeof(REAL));
      grid.gridfuncs.y_n_gfs = calloc(size, sizeof(REAL));
      MoL_malloc_intermediate_stage_gfs(&initial, &grid.params,
                                        &grid.gridfuncs);
      REAL *original = calloc(size, sizeof(REAL)),
           *dense = calloc(size, sizeof(REAL));
      if (!grid.gridfuncs.auxevol_gfs || !grid.gridfuncs.y_n_gfs || !original ||
          !dense)
        return 2;
      for (int gf = 0; gf < NUM_EVOL_GFS; gf++)
        for (int k = 0; k < dims[2]; k++)
          for (int j = 0; j < dims[1]; j++)
            for (int i = 0; i < dims[0]; i++) {
              int p = IDX4(gf, i, j, k);
              REAL arg = .09 * (gf + 1) + .37 * grid.xx[0][i] -
                         .28 * grid.xx[1][j] + .19 * grid.xx[2][k];
              original[p] = .008 * sin(arg) + .003 * cos(.7 * arg);
              if (gf == ALPHAGF || gf == CFGF)
                original[p] += 1;
            }
      for (int threads = 1; threads <= 4; threads += 3)
        for (int radiation = 0; radiation < 2; radiation++) {
          omp_set_num_threads(threads);
          commondata_struct common = initial;
          strcpy(common.outer_bc_type,
                 radiation ? "radiation" : "extrapolation");
          check_accumulation(&common, &grid, original, size, &max_stage_error);
          stage_cases += 4;
          for (int heun = 0; heun < 2; heun++) {
            common = initial;
            strcpy(common.outer_bc_type,
                   radiation ? "radiation" : "extrapolation");
            memcpy(grid.gridfuncs.y_n_gfs, original, size * sizeof(REAL));
            prepare_state(&common, &grid, grid.gridfuncs.y_n_gfs);
            memcpy(dense, grid.gridfuncs.y_n_gfs, size * sizeof(REAL));
            commondata_struct reference = common;
            REAL *const allocation = grid.gridfuncs.y_n_gfs;
            for (int step = 0; step < 3; step++) {
              REAL time = common.time;
              ricci_calls = rhs_calls = projection_calls = 0;
              counting = 1;
              if (heun)
                MoL_Heun(&common, &grid);
              else
                MoL_Ralston(&common, &grid);
              counting = 0;
              const REAL a21 = heun ? 1.0 / 3 : .5, a32 = heun ? 2.0 / 3 : .75;
              if (ricci_calls != 3 || rhs_calls != 3 || projection_calls != 3 ||
                  fabs(rhs_times[0] - time) > 1e-15 ||
                  fabs(rhs_times[1] - (time + a21 * common.dt)) > 1e-15 ||
                  fabs(rhs_times[2] - (time + a32 * common.dt)) > 1e-15 ||
                  grid.gridfuncs.y_n_gfs != allocation) {
                fprintf(
                    stderr,
                    "RK accounting failed: Ricci=%d RHS=%d projections=%d\n",
                    ricci_calls, rhs_calls, projection_calls);
                return 1;
              }
              dense_step(&reference, &grid, dense, heun, size);
              REAL error = compare(grid.gridfuncs.y_n_gfs, dense, size);
              if (error > max_step_error)
                max_step_error = error;
              if (common.nn != reference.nn || common.time != reference.time)
                return 1;
              step_cases++;
            }
          }
        }
      MoL_free_intermediate_stage_gfs(&grid.gridfuncs);
      free(grid.gridfuncs.y_n_gfs);
      free(grid.gridfuncs.auxevol_gfs);
      free(original);
      free(dense);
      for (int d = 0; d < 3; d++)
        free(grid.xx[d]);
      for (int i = 0; i < 3 * NGHOSTS; i++)
        free(grid.bcstruct.pure_outer_bc_array[i]);
      free(grid.bcstruct.inner_bc_array);
    }
  printf("PASS accumulated RHS: %d cases, max abs error %.17g\n", stage_cases,
         max_stage_error);
  printf("PASS complete projected RK3: %d steps, max abs error %.17g; three "
         "arrays, three Ricci/RHS calls per step\n",
         step_cases, max_step_error);
  return 0;
}
