const DendroScalar* const in_betaU0 = adm_gfs[0] + geom.component_offset;
const DendroScalar* const in_betaU1 = adm_gfs[1] + geom.component_offset;
const DendroScalar* const in_betaU2 = adm_gfs[2] + geom.component_offset;
const DendroScalar* const in_BU0 = adm_gfs[3] + geom.component_offset;
const DendroScalar* const in_BU1 = adm_gfs[4] + geom.component_offset;
const DendroScalar* const in_BU2 = adm_gfs[5] + geom.component_offset;
const DendroScalar* const in_gammaDD00 = adm_gfs[6] + geom.component_offset;
const DendroScalar* const in_gammaDD01 = adm_gfs[7] + geom.component_offset;
const DendroScalar* const in_gammaDD02 = adm_gfs[8] + geom.component_offset;
const DendroScalar* const in_gammaDD11 = adm_gfs[9] + geom.component_offset;
const DendroScalar* const in_gammaDD12 = adm_gfs[10] + geom.component_offset;
const DendroScalar* const in_gammaDD22 = adm_gfs[11] + geom.component_offset;
const DendroScalar* const in_KDD00 = adm_gfs[12] + geom.component_offset;
const DendroScalar* const in_KDD01 = adm_gfs[13] + geom.component_offset;
const DendroScalar* const in_KDD02 = adm_gfs[14] + geom.component_offset;
const DendroScalar* const in_KDD11 = adm_gfs[15] + geom.component_offset;
const DendroScalar* const in_KDD12 = adm_gfs[16] + geom.component_offset;
const DendroScalar* const in_KDD22 = adm_gfs[17] + geom.component_offset;
DendroScalar* const out_aDD00 = out_gfs[0] + geom.component_offset;
DendroScalar* const out_aDD01 = out_gfs[1] + geom.component_offset;
DendroScalar* const out_aDD02 = out_gfs[2] + geom.component_offset;
DendroScalar* const out_aDD11 = out_gfs[3] + geom.component_offset;
DendroScalar* const out_aDD12 = out_gfs[4] + geom.component_offset;
DendroScalar* const out_aDD22 = out_gfs[5] + geom.component_offset;
DendroScalar* const out_alpha = out_gfs[6] + geom.component_offset;
DendroScalar* const out_betU0 = out_gfs[7] + geom.component_offset;
DendroScalar* const out_betU1 = out_gfs[8] + geom.component_offset;
DendroScalar* const out_betU2 = out_gfs[9] + geom.component_offset;
DendroScalar* const out_cf = out_gfs[10] + geom.component_offset;
DendroScalar* const out_hDD00 = out_gfs[11] + geom.component_offset;
DendroScalar* const out_hDD01 = out_gfs[12] + geom.component_offset;
DendroScalar* const out_hDD02 = out_gfs[13] + geom.component_offset;
DendroScalar* const out_hDD11 = out_gfs[14] + geom.component_offset;
DendroScalar* const out_hDD12 = out_gfs[15] + geom.component_offset;
DendroScalar* const out_hDD22 = out_gfs[16] + geom.component_offset;
DendroScalar* const out_Theta_fCCZ4 = out_gfs[20] + geom.component_offset;
DendroScalar* const out_trK = out_gfs[21] + geom.component_offset;
DendroScalar* const out_vetU0 = out_gfs[22] + geom.component_offset;
DendroScalar* const out_vetU1 = out_gfs[23] + geom.component_offset;
DendroScalar* const out_vetU2 = out_gfs[24] + geom.component_offset;
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
const DendroScalar BU0 = in_BU0[pp];
const DendroScalar BU1 = in_BU1[pp];
const DendroScalar BU2 = in_BU2[pp];
const DendroScalar KDD00 = in_KDD00[pp];
const DendroScalar KDD01 = in_KDD01[pp];
const DendroScalar KDD02 = in_KDD02[pp];
const DendroScalar KDD11 = in_KDD11[pp];
const DendroScalar KDD12 = in_KDD12[pp];
const DendroScalar KDD22 = in_KDD22[pp];
const DendroScalar betaU0 = in_betaU0[pp];
const DendroScalar betaU1 = in_betaU1[pp];
const DendroScalar betaU2 = in_betaU2[pp];
const DendroScalar gammaDD00 = in_gammaDD00[pp];
const DendroScalar gammaDD01 = in_gammaDD01[pp];
const DendroScalar gammaDD02 = in_gammaDD02[pp];
const DendroScalar gammaDD11 = in_gammaDD11[pp];
const DendroScalar gammaDD12 = in_gammaDD12[pp];
const DendroScalar gammaDD22 = in_gammaDD22[pp];

/*
 * NRPy-Generated GF Access/FD Code, Step 2 of 2:
 * Evaluate SymPy expressions and write to main memory.
 */
const DendroScalar FDPart3tmp1 = gammaDD00*((gammaDD12)*(gammaDD12));
const DendroScalar FDPart3tmp3 = ((gammaDD01)*(gammaDD01))*gammaDD22;
const DendroScalar FDPart3tmp5 = ((gammaDD02)*(gammaDD02))*gammaDD11;
const DendroScalar FDPart3tmp6 = -FDPart3tmp1 - FDPart3tmp3 - FDPart3tmp5 + gammaDD00*gammaDD11*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12;
const DendroScalar FDPart3tmp7 = (1.0/(FDPart3tmp6));
const DendroScalar FDPart3tmp8 = cbrt(FDPart3tmp7);
const DendroScalar FDPart3tmp9 = 2*FDPart3tmp7;
const DendroScalar FDPart3tmp10 = FDPart3tmp7*KDD00*(gammaDD11*gammaDD22 - ((gammaDD12)*(gammaDD12))) + FDPart3tmp7*KDD11*(gammaDD00*gammaDD22 - ((gammaDD02)*(gammaDD02))) + FDPart3tmp7*KDD22*(gammaDD00*gammaDD11 - ((gammaDD01)*(gammaDD01))) + FDPart3tmp9*KDD01*(-gammaDD01*gammaDD22 + gammaDD02*gammaDD12) + FDPart3tmp9*KDD02*(gammaDD01*gammaDD12 - gammaDD02*gammaDD11) + FDPart3tmp9*KDD12*(-gammaDD00*gammaDD12 + gammaDD01*gammaDD02);
const DendroScalar FDPart3tmp11 = (1.0/3.0)*FDPart3tmp10;
out_aDD00[pp] = FDPart3tmp8*(-FDPart3tmp11*gammaDD00 + KDD00);
out_aDD01[pp] = FDPart3tmp8*(-FDPart3tmp11*gammaDD01 + KDD01);
out_aDD02[pp] = FDPart3tmp8*(-FDPart3tmp11*gammaDD02 + KDD02);
out_aDD11[pp] = FDPart3tmp8*(-FDPart3tmp11*gammaDD11 + KDD11);
out_aDD12[pp] = FDPart3tmp8*(-FDPart3tmp11*gammaDD12 + KDD12);
out_aDD22[pp] = FDPart3tmp8*(-FDPart3tmp11*gammaDD22 + KDD22);
out_alpha[pp] = 1.0;
out_betU0[pp] = BU0;
out_betU1[pp] = BU1;
out_betU2[pp] = BU2;
out_cf[pp] = (1.0/cbrt(FDPart3tmp6/(-FDPart3tmp1*FDPart3tmp7 - FDPart3tmp3*FDPart3tmp7 - FDPart3tmp5*FDPart3tmp7 + FDPart3tmp7*gammaDD00*gammaDD11*gammaDD22 + 2*FDPart3tmp7*gammaDD01*gammaDD02*gammaDD12)));
out_hDD00[pp] = FDPart3tmp8*gammaDD00 - 1;
out_hDD01[pp] = FDPart3tmp8*gammaDD01;
out_hDD02[pp] = FDPart3tmp8*gammaDD02;
out_hDD11[pp] = FDPart3tmp8*gammaDD11 - 1;
out_hDD12[pp] = FDPart3tmp8*gammaDD12;
out_hDD22[pp] = FDPart3tmp8*gammaDD22 - 1;
out_Theta_fCCZ4[pp] = 0;
out_trK[pp] = FDPart3tmp10;
out_vetU0[pp] = betaU0;
out_vetU1[pp] = betaU1;
out_vetU2[pp] = betaU2;

} // END LOOP: for i0 over [static_cast<int>(padding), static_cast<int>(nx - padding))
} // END LOOP: for i1 over [static_cast<int>(padding), static_cast<int>(ny - padding))
} // END LOOP: for i2 over [static_cast<int>(padding), static_cast<int>(nz - padding))
