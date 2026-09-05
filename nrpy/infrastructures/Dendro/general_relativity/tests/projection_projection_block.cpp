const DendroScalar* const in_aDD00 = state_gfs[0] + geom.component_offset;
const DendroScalar* const in_aDD01 = state_gfs[1] + geom.component_offset;
const DendroScalar* const in_aDD02 = state_gfs[2] + geom.component_offset;
const DendroScalar* const in_aDD11 = state_gfs[3] + geom.component_offset;
const DendroScalar* const in_aDD12 = state_gfs[4] + geom.component_offset;
const DendroScalar* const in_aDD22 = state_gfs[5] + geom.component_offset;
const DendroScalar* const in_hDD00 = state_gfs[11] + geom.component_offset;
const DendroScalar* const in_hDD01 = state_gfs[12] + geom.component_offset;
const DendroScalar* const in_hDD02 = state_gfs[13] + geom.component_offset;
const DendroScalar* const in_hDD11 = state_gfs[14] + geom.component_offset;
const DendroScalar* const in_hDD12 = state_gfs[15] + geom.component_offset;
const DendroScalar* const in_hDD22 = state_gfs[16] + geom.component_offset;
DendroScalar* const out_hDD00 = state_gfs[11] + geom.component_offset;
DendroScalar* const out_hDD01 = state_gfs[12] + geom.component_offset;
DendroScalar* const out_hDD02 = state_gfs[13] + geom.component_offset;
DendroScalar* const out_hDD11 = state_gfs[14] + geom.component_offset;
DendroScalar* const out_hDD12 = state_gfs[15] + geom.component_offset;
DendroScalar* const out_hDD22 = state_gfs[16] + geom.component_offset;
DendroScalar* const out_aDD00 = state_gfs[0] + geom.component_offset;
DendroScalar* const out_aDD01 = state_gfs[1] + geom.component_offset;
DendroScalar* const out_aDD02 = state_gfs[2] + geom.component_offset;
DendroScalar* const out_aDD11 = state_gfs[3] + geom.component_offset;
DendroScalar* const out_aDD12 = state_gfs[4] + geom.component_offset;
DendroScalar* const out_aDD22 = state_gfs[5] + geom.component_offset;
const std::ptrdiff_t nx = static_cast<std::ptrdiff_t>(geom.nx);
const std::ptrdiff_t ny = static_cast<std::ptrdiff_t>(geom.ny);
const std::ptrdiff_t nz = static_cast<std::ptrdiff_t>(geom.nz);
[[maybe_unused]] const std::ptrdiff_t nxy = nx * ny;
const std::ptrdiff_t padding = static_cast<std::ptrdiff_t>(0);
[[maybe_unused]] const DendroScalar invdxx0 = static_cast<DendroScalar>(1) / geom.dx[0];
[[maybe_unused]] const DendroScalar invdxx1 = static_cast<DendroScalar>(1) / geom.dx[1];
[[maybe_unused]] const DendroScalar invdxx2 = static_cast<DendroScalar>(1) / geom.dx[2];
for (int i2 = static_cast<int>(padding); i2 < static_cast<int>(nz - padding); i2++) {
for (int i1 = static_cast<int>(padding); i1 < static_cast<int>(ny - padding); i1++) {
for (int i0 = static_cast<int>(padding); i0 < static_cast<int>(nx - padding); i0++) {
const std::ptrdiff_t pp = i0 + nx * (i1 + ny * i2);
[[maybe_unused]] const DendroScalar xx0 = geom.pmin_padded[0] + static_cast<DendroScalar>(i0) * geom.dx[0];
[[maybe_unused]] const DendroScalar xx1 = geom.pmin_padded[1] + static_cast<DendroScalar>(i1) * geom.dx[1];
[[maybe_unused]] const DendroScalar xx2 = geom.pmin_padded[2] + static_cast<DendroScalar>(i2) * geom.dx[2];
/*
 * NRPy-Generated GF Access/FD Code, Step 1 of 2:
 * Read gridfunction(s) from main memory and compute FD stencils as needed.
 */
const DendroScalar aDD00 = in_aDD00[pp];
const DendroScalar aDD01 = in_aDD01[pp];
const DendroScalar aDD02 = in_aDD02[pp];
const DendroScalar aDD11 = in_aDD11[pp];
const DendroScalar aDD12 = in_aDD12[pp];
const DendroScalar aDD22 = in_aDD22[pp];
const DendroScalar hDD00 = in_hDD00[pp];
const DendroScalar hDD01 = in_hDD01[pp];
const DendroScalar hDD02 = in_hDD02[pp];
const DendroScalar hDD11 = in_hDD11[pp];
const DendroScalar hDD12 = in_hDD12[pp];
const DendroScalar hDD22 = in_hDD22[pp];

/*
 * NRPy-Generated GF Access/FD Code, Step 2 of 2:
 * Evaluate SymPy expressions and write to main memory.
 */
const DendroScalar FDPart3tmp1 = hDD22 + 1;
const DendroScalar FDPart3tmp3 = hDD11 + 1;
const DendroScalar FDPart3tmp5 = hDD00 + 1;
const DendroScalar FDPart3tmp6 = FDPart3tmp1*FDPart3tmp3*FDPart3tmp5 - FDPart3tmp1*((hDD01)*(hDD01)) - FDPart3tmp3*((hDD02)*(hDD02)) - FDPart3tmp5*((hDD12)*(hDD12)) + 2*hDD01*hDD02*hDD12;
const DendroScalar FDPart3tmp7 = (1.0/(FDPart3tmp6));
const DendroScalar FDPart3tmp10 = (1.0/cbrt(FDPart3tmp6));
const DendroScalar FDPart3tmp8 = 2*FDPart3tmp7;
const DendroScalar FDPart3tmp9 = FDPart3tmp7*aDD00*(FDPart3tmp1*FDPart3tmp3 - ((hDD12)*(hDD12))) + FDPart3tmp7*aDD11*(FDPart3tmp1*FDPart3tmp5 - ((hDD02)*(hDD02))) + FDPart3tmp7*aDD22*(FDPart3tmp3*FDPart3tmp5 - ((hDD01)*(hDD01))) + FDPart3tmp8*aDD01*(-FDPart3tmp1*hDD01 + hDD02*hDD12) + FDPart3tmp8*aDD02*(-FDPart3tmp3*hDD02 + hDD01*hDD12) + FDPart3tmp8*aDD12*(-FDPart3tmp5*hDD12 + hDD01*hDD02);
const DendroScalar FDPart3tmp11 = (1.0/3.0)*FDPart3tmp9;
const DendroScalar det_ratio = FDPart3tmp6;
const DendroScalar trace_residual = FDPart3tmp9;
const DendroScalar projected_hDD00 = FDPart3tmp10*FDPart3tmp5 - 1;
const DendroScalar projected_hDD01 = FDPart3tmp10*hDD01;
const DendroScalar projected_hDD02 = FDPart3tmp10*hDD02;
const DendroScalar projected_hDD11 = FDPart3tmp10*FDPart3tmp3 - 1;
const DendroScalar projected_hDD12 = FDPart3tmp10*hDD12;
const DendroScalar projected_hDD22 = FDPart3tmp1*FDPart3tmp10 - 1;
const DendroScalar projected_aDD00 = -FDPart3tmp9*((1.0/3.0)*hDD00 + 1.0/3.0) + aDD00;
const DendroScalar projected_aDD01 = -FDPart3tmp11*hDD01 + aDD01;
const DendroScalar projected_aDD02 = -FDPart3tmp11*hDD02 + aDD02;
const DendroScalar projected_aDD11 = -FDPart3tmp9*((1.0/3.0)*hDD11 + 1.0/3.0) + aDD11;
const DendroScalar projected_aDD12 = -FDPart3tmp11*hDD12 + aDD12;
const DendroScalar projected_aDD22 = -FDPart3tmp9*((1.0/3.0)*hDD22 + 1.0/3.0) + aDD22;
if (!(det_ratio > 0) || !std::isfinite(det_ratio) ||
    !std::isfinite(trace_residual)) {
  if (!std::isfinite(det_ratio) || !std::isfinite(trace_residual)) {
    status->nonfinite_points += 1;
  } // END IF: nonfinite determinant or trace
  if (status->failed_points == 0) {
    status->first_failing_index = static_cast<long long>(pp);
    if (status->first_failing_field < 0 && !std::isfinite(in_aDD00[pp])) status->first_failing_field = 0;
    if (status->first_failing_field < 0 && !std::isfinite(in_aDD01[pp])) status->first_failing_field = 1;
    if (status->first_failing_field < 0 && !std::isfinite(in_aDD02[pp])) status->first_failing_field = 2;
    if (status->first_failing_field < 0 && !std::isfinite(in_aDD11[pp])) status->first_failing_field = 3;
    if (status->first_failing_field < 0 && !std::isfinite(in_aDD12[pp])) status->first_failing_field = 4;
    if (status->first_failing_field < 0 && !std::isfinite(in_aDD22[pp])) status->first_failing_field = 5;
    if (status->first_failing_field < 0 && !std::isfinite(in_hDD00[pp])) status->first_failing_field = 11;
    if (status->first_failing_field < 0 && !std::isfinite(in_hDD01[pp])) status->first_failing_field = 12;
    if (status->first_failing_field < 0 && !std::isfinite(in_hDD02[pp])) status->first_failing_field = 13;
    if (status->first_failing_field < 0 && !std::isfinite(in_hDD11[pp])) status->first_failing_field = 14;
    if (status->first_failing_field < 0 && !std::isfinite(in_hDD12[pp])) status->first_failing_field = 15;
    if (status->first_failing_field < 0 && !std::isfinite(in_hDD22[pp])) status->first_failing_field = 16;
  } // END IF: first refused point
  status->failed_points += 1;
  continue;
} // END IF: nonpositive or nonfinite determinant
status->projected_points += 1;
status->max_abs_det_minus_one = std::fmax(
    status->max_abs_det_minus_one,
    static_cast<double>(std::fabs(det_ratio - 1)));
status->max_abs_trace_residual = std::fmax(
    status->max_abs_trace_residual,
    static_cast<double>(std::fabs(trace_residual)));
out_hDD00[pp] = projected_hDD00;
out_hDD01[pp] = projected_hDD01;
out_hDD02[pp] = projected_hDD02;
out_hDD11[pp] = projected_hDD11;
out_hDD12[pp] = projected_hDD12;
out_hDD22[pp] = projected_hDD22;
out_aDD00[pp] = projected_aDD00;
out_aDD01[pp] = projected_aDD01;
out_aDD02[pp] = projected_aDD02;
out_aDD11[pp] = projected_aDD11;
out_aDD12[pp] = projected_aDD12;
out_aDD22[pp] = projected_aDD22;

} // END LOOP: for i0 over [static_cast<int>(padding), static_cast<int>(nx - padding))
} // END LOOP: for i1 over [static_cast<int>(padding), static_cast<int>(ny - padding))
} // END LOOP: for i2 over [static_cast<int>(padding), static_cast<int>(nz - padding))
