// Define points for output along the xy-plane in SinhSymTP coordinates.
const int numpts_i0 = Nxx0, numpts_i1 = 1, numpts_i2 = Nxx2;
int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];
#pragma omp parallel for
for (int i0 = NGHOSTS; i0 < Nxx0 + NGHOSTS; i0++)
  i0_pts[i0 - NGHOSTS] = i0;
i1_pts[0] = (int)((1.0 / 2.0) * Nxx_plus_2NGHOSTS1);
#pragma omp parallel for
for (int i2 = NGHOSTS; i2 < Nxx2 + NGHOSTS; i2++)
  i2_pts[i2 - NGHOSTS] = i2;
// Main loop to store data points in the specified plane
LOOP_NOOMP(i0_pt, 0, numpts_i0, i1_pt, 0, numpts_i1, i2_pt, 0, numpts_i2) {
  const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
  const int idx3 = IDX3(i0, i1, i2);
  REAL xCart[3];
  xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
  {
    // Collect diagnostic data
    const REAL log10HL = log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16));
    fprintf(outfile, "%.15e %.15e %.15e\n", xCart[0], xCart[1], log10HL);
  }
}
