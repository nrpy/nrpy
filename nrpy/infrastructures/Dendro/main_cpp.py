# MIT License
# Copied numerical routine: Dendro-GR/BSSN_GR/src/grUtils.cpp,
# punctureDataPhysicalCoord. The function signature and BH parameter access
# were adapted for the generated NRPy solver; its numerical body is unchanged.
#
# Copyright (c) 2018 DendroGR
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is furnished
# to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Emit the executable entry point for a generated Dendro application.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.params as par
from nrpy.infrastructures.Dendro.CodeParameters import output_toml_bindings
from nrpy.infrastructures.Dendro.state_h import BSSN_EVOLVED_GRIDFUNCTIONS

DENDRO_LICENSE = """// MIT License
// Source: Dendro-GR/BSSN_GR/src/grUtils.cpp, punctureDataPhysicalCoord.
// Adapted signature and parameter access; numerical body unchanged.
//
// Copyright (c) 2018 DendroGR
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is furnished
// to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

"""


DENDRO_PUNCTURE_SEED = r"""
// Dendro-GR BSSN_GR/src/grUtils.cpp: punctureDataPhysicalCoord.
// Numerical expressions and branch conditions are retained from Dendro-GR.
namespace nrpy_dendro_seed {
struct PunctureParameters {
    double mass, x, y, z, vx, vy, vz, spin, spin_theta, spin_phi;
};
namespace VAR {
enum {
    U_ALPHA, U_CHI, U_K, U_GT0, U_GT1, U_GT2,
    U_BETA0, U_BETA1, U_BETA2, U_B0, U_B1, U_B2,
    U_SYMGT0, U_SYMGT1, U_SYMGT2, U_SYMGT3, U_SYMGT4, U_SYMGT5,
    U_SYMAT0, U_SYMAT1, U_SYMAT2, U_SYMAT3, U_SYMAT4, U_SYMAT5
};
}
void punctureDataPhysicalCoord(const double xx, const double yy,
                               const double zz, double* var,
                               const PunctureParameters& BH1,
                               const PunctureParameters& BH2,
                               const double CHI_FLOOR) {
    /* Define the Levi-Cevita pseudo-tensor and Kroneckar delta */
    double epijk[3][3][3];
    int i, j, k;
    for (k = 0; k < 3; k++) {
        for (j = 0; j < 3; j++) {
            for (i = 0; i < 3; i++) {
                epijk[k][j][i] = 0.0;
            }
        }
    }
    epijk[0][1][2] = 1.0;
    epijk[1][2][0] = 1.0;
    epijk[2][0][1] = 1.0;
    epijk[0][2][1] = -1.0;
    epijk[2][1][0] = -1.0;
    epijk[1][0][2] = -1.0;

    double deltaij[3][3];
    for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++) {
            deltaij[j][i] = 0.0;
        }
    }

    deltaij[0][0] = 1.0;
    deltaij[1][1] = 1.0;
    deltaij[2][2] = 1.0;

    double x1, y1, z1, rv1;
    double x2, y2, z2, rv2;
    double vn1[3], vn2[3];

    double vpsibl;
    double v_u_corr, amp_capj, amp_capr, l_r, u0_j, u2_j, mu_j, p2_mu_j, v_u_j1;
    double v1, v2, v3, v4, vt1, vt2;

    int i1, i2, i3, i4;
    double amp_capp, u0_p, u2_p, mu_p, p2_mu_p;
    double v_u_p1, v_u_c1, v_u_j2, v_u_p2;
    double v_u_c2, vpsibl_u, vpsibl_u2;

    // bh 1
    double mass1 = BH1.mass;
    double bh1x  = BH1.x;
    double bh1y  = BH1.y;
    double bh1z  = BH1.z;

    double vp1[3];
    vp1[0]          = BH1.vx;
    vp1[1]          = BH1.vy;
    vp1[2]          = BH1.vz;

    double vp1tot   = sqrt(vp1[0] * vp1[0] + vp1[1] * vp1[1] + vp1[2] * vp1[2]);
    double spin1    = BH1.spin;
    double spin1_th = BH1.spin_theta;
    double spin1_phi = BH1.spin_phi;
    double vs1[3];

    vs1[0]       = spin1 * sin(spin1_th) * cos(spin1_phi);
    vs1[1]       = spin1 * sin(spin1_th) * sin(spin1_phi);
    vs1[2]       = spin1 * cos(spin1_th);

    // bh 2
    double mass2 = BH2.mass;
    double bh2x  = BH2.x;
    double bh2y  = BH2.y;
    double bh2z  = BH2.z;

    double vp2[3];
    vp2[0]          = BH2.vx;
    vp2[1]          = BH2.vy;
    vp2[2]          = BH2.vz;

    double vp2tot   = sqrt(vp2[0] * vp2[0] + vp2[1] * vp2[1] + vp2[2] * vp2[2]);
    double spin2    = BH2.spin;
    double spin2_th = BH2.spin_theta;
    double spin2_phi = BH2.spin_phi;

    double vs2[3];
    vs2[0]   = spin2 * sin(spin2_th) * cos(spin2_phi);
    vs2[1]   = spin2 * sin(spin2_th) * sin(spin2_phi);
    vs2[2]   = spin2 * cos(spin2_th);

    // coordinates with respect to center of bh1
    x1       = xx - bh1x;
    y1       = yy - bh1y;
    z1       = zz - bh1z;

    // locating as a radial form
    rv1      = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    vn1[0]   = x1 / rv1;
    vn1[1]   = y1 / rv1;
    vn1[2]   = z1 / rv1;

    // same as BH2
    x2       = xx - bh2x;
    y2       = yy - bh2y;
    z2       = zz - bh2z;

    rv2      = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    vn2[0]   = x2 / rv2;
    vn2[1]   = y2 / rv2;
    vn2[2]   = z2 / rv2;

    // Initial data is related with the paper: http://arxiv.org/abs/0711.1165
    // Brill-Lindquist conformal factor
    vpsibl   = 1.0 + mass1 / (2.0 * rv1);
    vpsibl   = vpsibl + mass2 / (2.0 * rv2);

    v_u_corr = 0.0;
    // bh 1

    // For spinning puncture
    if (fabs(spin1) > 1.e-6) {
        amp_capj = 4.0 * spin1 / (mass1 * mass1);
        amp_capr = 2.0 * rv1 / mass1;
        l_r      = 1.0 / (1.0 + amp_capr);
        u0_j =
            (l_r + l_r * l_r + l_r * l_r * l_r - 4.0 * l_r * l_r * l_r * l_r +
             2.0 * l_r * l_r * l_r * l_r * l_r) /
            40.0;
        u2_j    = -pow(l_r, 5) / 20.0;
        mu_j    = vn1[0] * vs1[0];
        mu_j    = mu_j + vn1[1] * vs1[1];
        mu_j    = (mu_j + vn1[2] * vs1[2]) / fabs(spin1);
        p2_mu_j = (3.0 * mu_j * mu_j - 1.0) / 2.0;
        v_u_j1 =
            amp_capj * amp_capj * (u0_j + u2_j * amp_capr * amp_capr * p2_mu_j);
        v_u_corr = v_u_corr + v_u_j1;
    }
    // For boosting puncture
    if (vp1tot > 1.e-6) {
        amp_capp = 2.0 * vp1tot / mass1;
        amp_capr = 2.0 * rv1 / mass1;
        l_r      = 1.0 / (1.0 + amp_capr);
        u0_p     = l_r - 2.0 * l_r * l_r + 2.0 * pow(l_r, 3);
        u0_p     = (u0_p - pow(l_r, 4) + 0.20 * pow(l_r, 5)) * (5.0 / 32.0);
        u2_p     = 15.0 * l_r + 132.0 * l_r * l_r + 53.0 * pow(l_r, 3);
        u2_p     = u2_p + 96.0 * pow(l_r, 4) + 82.0 * pow(l_r, 5);
        u2_p = u2_p + (84.0 / amp_capr) * (pow(l_r, 5) + log(l_r) / amp_capr);
        u2_p = (u2_p) / (80.0 * amp_capr);
        mu_p = vn1[0] * vp1[0] / vp1tot;
        mu_p = mu_p + vn1[1] * vp1[1] / vp1tot;
        mu_p = mu_p + vn1[2] * vp1[2] / vp1tot;
        p2_mu_p  = (3.0 * pow(mu_p, 2) - 1.0) / 2.0;
        v_u_p1   = pow(amp_capp, 2) * (u0_p + u2_p * p2_mu_p);
        v_u_corr = v_u_corr + v_u_p1;
    }
    // For spinning boosted pucture
    if (vp1tot > 1.e-6 && fabs(spin1) > 1.e-6) {
        v1       = (vp1[1] * vs1[2] - vp1[2] * vs1[1]) * vn1[0];
        v1       = v1 + (vp1[2] * vs1[0] - vp1[0] * vs1[2]) * vn1[1];
        v1       = v1 + (vp1[0] * vs1[1] - vp1[1] * vs1[0]) * vn1[2];
        v1       = v1 * (16.0 / pow(mass1, 4)) * rv1;

        amp_capr = 2.0 * rv1 / mass1;
        l_r      = 1.0 / (1.0 + amp_capr);

        v2       = 1.0 + 5.0 * amp_capr + 10.0 * pow(amp_capr, 2);

        v_u_c1   = (v1 * v2 * pow(l_r, 5)) / 80.0;
        v_u_corr = v_u_corr + v_u_c1;
    }
    // bh 2 same puncture as bh 1
    if (fabs(spin2) > 1.e-6) {
        amp_capj = 4.0 * spin2 / (mass2 * mass2);
        amp_capr = 2.0 * rv2 / mass2;
        l_r      = 1.0 / (1.0 + amp_capr);
        u0_j =
            (l_r + l_r * l_r + l_r * l_r * l_r - 4.0 * l_r * l_r * l_r * l_r +
             2.0 * l_r * l_r * l_r * l_r * l_r) /
            40.0;
        u2_j    = -pow(l_r, 5) / 20.0;
        mu_j    = vn2[0] * vs2[0];
        mu_j    = mu_j + vn2[1] * vs2[1];
        mu_j    = (mu_j + vn2[2] * vs2[2]) / fabs(spin2);
        p2_mu_j = (3.0 * mu_j * mu_j - 1.0) / 2.0;
        v_u_j2 =
            amp_capj * amp_capj * (u0_j + u2_j * amp_capr * amp_capr * p2_mu_j);
        v_u_corr = v_u_corr + v_u_j2;
    }

    if (vp2tot > 1.e-6) {
        amp_capp = 2.0 * vp2tot / mass2;
        amp_capr = 2.0 * rv2 / mass2;
        l_r      = 1.0 / (1.0 + amp_capr);
        u0_p     = l_r - 2.0 * l_r * l_r + 2.0 * pow(l_r, 3);
        u0_p     = (u0_p - pow(l_r, 4) + 0.20 * pow(l_r, 5)) * (5.0 / 32.0);
        u2_p     = 15.0 * l_r + 132.0 * l_r * l_r + 53.0 * pow(l_r, 3);
        u2_p     = u2_p + 96.0 * pow(l_r, 4) + 82.0 * pow(l_r, 5);
        u2_p = u2_p + (84.0 / amp_capr) * (pow(l_r, 5) + log(l_r) / amp_capr);
        u2_p = (u2_p) / (80.0 * amp_capr);
        mu_p = vn2[0] * vp2[0] / vp2tot;
        mu_p = mu_p + vn2[1] * vp2[1] / vp2tot;
        mu_p = mu_p + vn2[2] * vp2[2] / vp2tot;
        p2_mu_p  = (3.0 * pow(mu_p, 2) - 1.0) / 2.0;
        v_u_p2   = pow(amp_capp, 2) * (u0_p + u2_p * p2_mu_p);
        v_u_corr = v_u_corr + v_u_p2;
    }

    if (vp2tot > 1.e-6 && fabs(spin2) > 1.e-6) {
        v1       = (vp2[1] * vs2[2] - vp2[2] * vs2[1]) * vn2[0];
        v1       = v1 + (vp2[2] * vs2[0] - vp2[0] * vs2[2]) * vn2[1];
        v1       = v1 + (vp2[0] * vs2[1] - vp2[1] * vs2[0]) * vn2[2];
        v1       = v1 * (16.0 / pow(mass2, 4)) * rv2;

        amp_capr = 2.0 * rv2 / mass2;
        l_r      = 1.0 / (1.0 + amp_capr);

        v2       = 1.0 + 5.0 * amp_capr + 10.0 * pow(amp_capr, 2);

        v_u_c2   = (v1 * v2 * pow(l_r, 5)) / 80.0;
        v_u_corr = v_u_corr + v_u_c2;
    }

    // vpsibl_u will be used for the conformal factor,
    vpsibl_u          = vpsibl + v_u_corr;
    // vpsibl_u2 is for the Aij terms...
    // ! since the corrections are first order...
    // ! adding half of the correction seems to give the best results...
    // ! update - do a fit for spin = 0.6...
    vpsibl_u2         = vpsibl + v_u_corr;

    var[VAR::U_ALPHA] = 1.0 / (vpsibl_u * vpsibl_u);
    // std::cout<<"Alpha: "<<u[U_ALPHA]<<" vpsibl_u: "<< vpsibl_u<<std::endl;
    var[VAR::U_ALPHA] = std::max(var[VAR::U_ALPHA], CHI_FLOOR);

    v2                = 1.0 / pow(vpsibl_u, 4);
    var[VAR::U_CHI]   = v2;

    if (var[VAR::U_CHI] < CHI_FLOOR) var[VAR::U_CHI] = CHI_FLOOR;

    var[VAR::U_K]      = 0.0;

    var[VAR::U_BETA0]  = 0.0;
    var[VAR::U_BETA1]  = 0.0;
    var[VAR::U_BETA2]  = 0.0;

    var[VAR::U_GT0]    = 0.0;
    var[VAR::U_GT1]    = 0.0;
    var[VAR::U_GT2]    = 0.0;

    var[VAR::U_B0]     = 0.0;
    var[VAR::U_B1]     = 0.0;
    var[VAR::U_B2]     = 0.0;

    var[VAR::U_SYMGT0] = 1.0;  // XX
    var[VAR::U_SYMGT1] = 0.0;  // XY
    var[VAR::U_SYMGT2] = 0.0;  // XZ
    var[VAR::U_SYMGT3] = 1.0;  // YY
    var[VAR::U_SYMGT4] = 0.0;  // YZ
    var[VAR::U_SYMGT5] = 1.0;  // ZZ

    for (i1 = 0; i1 < 3; i1++) {
        for (i2 = 0; i2 < 3; i2++) {
            // first BH
            v2 = 0.0;
            for (i3 = 0; i3 < 3; i3++) {
                for (i4 = 0; i4 < 3; i4++) {
                    vt1 = epijk[i1][i3][i4] * vs1[i3] * vn1[i4] * vn1[i2];
                    vt2 = epijk[i2][i3][i4] * vs1[i3] * vn1[i4] * vn1[i1];
                    v2  = v2 + vt1 + vt2;
                }
            }

            v3  = vp1[i1] * vn1[i2] + vp1[i2] * vn1[i1];
            vt1 = 0.0;
            for (i3 = 0; i3 < 3; i3++) {
                vt1 = vt1 + vp1[i3] * vn1[i3];
            }
            vt1 = vt1 * (vn1[i1] * vn1[i2] - deltaij[i1][i2]);
            v3  = v3 + vt1;

            v1  = 3.0 / (pow(vpsibl_u2, 6) * pow(rv1, 3));
            v4  = v1 * (v2 + (rv1 / 2.0) * v3);

            // second BH
            v2  = 0.0;
            for (i3 = 0; i3 < 3; i3++) {
                for (i4 = 0; i4 < 3; i4++) {
                    vt1 = epijk[i1][i3][i4] * vs2[i3] * vn2[i4] * vn2[i2];
                    vt2 = epijk[i2][i3][i4] * vs2[i3] * vn2[i4] * vn2[i1];
                    v2  = v2 + vt1 + vt2;
                }
            }

            v3  = vp2[i1] * vn2[i2] + vp2[i2] * vn2[i1];
            vt1 = 0.0;
            for (i3 = 0; i3 < 3; i3++) {
                vt1 = vt1 + vp2[i3] * vn2[i3];
            }
            vt1 = vt1 * (vn2[i1] * vn2[i2] - deltaij[i1][i2]);
            v3  = v3 + vt1;

            v1  = 3.0 / (pow(vpsibl_u2, 6) * pow(rv2, 3));
            v4  = v4 + v1 * (v2 + (rv2 / 2.0) * v3);

            if (i1 == 0 && i2 == 0) {
                var[VAR::U_SYMAT0] = v4;  // XX
            } else if (i1 == 0 && i2 == 1) {
                var[VAR::U_SYMAT1] = v4;  // XY
            } else if (i1 == 0 && i2 == 2) {
                var[VAR::U_SYMAT2] = v4;  // XZ
            } else if (i1 == 1 && i2 == 1) {
                var[VAR::U_SYMAT3] = v4;  // YY
            } else if (i1 == 1 && i2 == 2) {
                var[VAR::U_SYMAT4] = v4;  // YZ
            } else if (i1 == 2 && i2 == 2) {
                var[VAR::U_SYMAT5] = v4;  // ZZ
            }
        }
    }
}
}  // namespace nrpy_dendro_seed

"""


def output_main_cpp(
    solver_stem: str,
    solver_namespace: str,
    executable_name: str,
    profile_name: str,
) -> str:
    """
    Emit the complete TwoPunctures, mesh, and RK driver.

    :param solver_stem: Lowercase formulation name used in generated files.
    :param solver_namespace: C++ namespace owning the generated application.
    :param executable_name: Executable name printed in usage diagnostics.
    :param profile_name: Generated numerical-profile identity.
    :return: Complete generated C++ source.
    """
    horizon_fields = (
        "cf",
        "trK",
        "aDD00",
        "aDD01",
        "aDD02",
        "aDD11",
        "aDD12",
        "aDD22",
        "hDD00",
        "hDD01",
        "hDD02",
        "hDD11",
        "hDD12",
        "hDD22",
    )
    horizon_indices = ", ".join(
        str(BSSN_EVOLVED_GRIDFUNCTIONS.index(name)) for name in horizon_fields
    )
    inverse_chi_expression = (
        "1.0 / values[0]"
        if par.parval_from_str("EvolvedConformalFactor_cf") == "chi"
        else "1.0 / (values[0] * values[0])"
    )
    return (
        DENDRO_LICENSE + "// GENERATED FILE - DO NOT EDIT\n"
        "// AUTOMATICALLY GENERATED BY NRPy\n"
        f'#include "{solver_stem}Ctx.h"\n'
        r"""#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "ets.h"
#include "meshUtils.h"
#include "octUtils.h"
#include <toml.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

void allocate_derivs(derivs* derivatives, int count);

"""
        + DENDRO_PUNCTURE_SEED
        + r"""
int main(int argc, char** argv) {
  const bool generate_tpid =
      argc == 3 && std::strcmp(argv[1], "--tpid") == 0;
  if (argc != 2 && !generate_tpid) {
    std::cerr << "usage: """
        + executable_name
        + r""" [--tpid] PARAM_FILE\n";
    return 2;
  }
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    int mpi_tasks = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
    if (generate_tpid && mpi_tasks != 1)
      throw std::runtime_error("--tpid requires exactly one MPI task");
    const toml::value document = toml::parse(argv[generate_tpid ? 2 : 1]);
    const unsigned element_order =
        toml::find_or<unsigned>(document, "BSSN_ELE_ORDER", """
        + solver_namespace
        + r"""::generated::FD_ORDER);
    if (element_order != 4 && element_order != 6 && element_order != 8)
      throw std::runtime_error("BSSN_ELE_ORDER must be 4, 6, or 8");
    const unsigned minimum_depth =
        toml::find_or<unsigned>(document, "BSSN_MINDEPTH", 4);
    const unsigned maximum_depth =
        toml::find_or<unsigned>(document, "BSSN_MAXDEPTH", 14);
    if (minimum_depth > maximum_depth || maximum_depth >= 31)
      throw std::runtime_error("invalid octree depth range");
    const DendroScalar wavelet_tolerance =
        toml::find_or<DendroScalar>(document, "BSSN_WAVELET_TOL", 1.0e-5);
    const unsigned refinement_mode =
        toml::find_or<unsigned>(document, "BSSN_REFINEMENT_MODE", 4);
    const unsigned wavelet_tolerance_mode =
        toml::find_or<unsigned>(document, "BSSN_USE_WAVELET_TOL_FUNCTION", 0);
    const DendroScalar maximum_wavelet_tolerance = toml::find_or<DendroScalar>(
        document, "BSSN_WAVELET_TOL_MAX", wavelet_tolerance);
    const DendroScalar gravitational_wave_tolerance = toml::find_or<DendroScalar>(
        document, "BSSN_GW_REFINE_WTOL", wavelet_tolerance);
    const DendroScalar amr_coarsening_factor = toml::find_or<DendroScalar>(
        document, "BSSN_DENDRO_AMR_FAC", 0.1);
    const DendroScalar postmerger_amr_coarsening_factor = toml::find_or<DendroScalar>(
        document, "BSSN_DENDRO_AMR_FAC_POST_MERGER", 0.0);
    std::vector<unsigned> refinement_variables;
    for (unsigned field = 0; field <
         """
        + solver_namespace
        + r"""::generated::NUM_EVOL_GFS; ++field)
      refinement_variables.push_back(field);
    refinement_variables = toml::find_or<std::vector<unsigned>>(
        document, "BSSN_REFINE_VARIABLE_INDICES", refinement_variables);
    const unsigned number_refinement_variables = toml::find_or<unsigned>(
        document, "BSSN_NUM_REFINE_VARS",
        static_cast<unsigned>(refinement_variables.size()));
    if (number_refinement_variables == 0 ||
        number_refinement_variables > refinement_variables.size())
      throw std::runtime_error("invalid BSSN_NUM_REFINE_VARS");
    refinement_variables.resize(number_refinement_variables);
    for (const unsigned field : refinement_variables)
      if (field >= """
        + solver_namespace
        + r"""::generated::NUM_EVOL_GFS)
        throw std::runtime_error("BSSN_REFINE_VARIABLE_INDICES out of range");
    if (refinement_mode != 4 || (wavelet_tolerance_mode != 0 &&
                                  wavelet_tolerance_mode != 6))
      throw std::runtime_error(
          "generated Dendro solver supports BH_WAMR (mode 4) with "
          "constant or causal wavelet tolerance (mode 0 or 6)");
    const DendroScalar cfl =
        toml::find_or<DendroScalar>(document, "BSSN_CFL_FACTOR", 0.25);
    const DendroScalar time_begin =
        toml::find_or<DendroScalar>(document, "BSSN_RK_TIME_BEGIN", 0.0);
    const DendroScalar time_end =
        toml::find_or<DendroScalar>(document, "BSSN_RK_TIME_END", 700.0);
    const unsigned maximum_iterations = toml::find_or<unsigned>(
        document, "BSSN_MAX_ITERATIONS", std::numeric_limits<unsigned>::max());
    const unsigned remesh_frequency =
        toml::find_or<unsigned>(document, "BSSN_REMESH_TEST_FREQ", 50);
    const unsigned postmerger_remesh_frequency = toml::find_or<unsigned>(
        document, "BSSN_REMESH_TEST_FREQ_AFTER_MERGER", 10);
    const unsigned initial_grid_iterations =
        toml::find_or<unsigned>(document, "BSSN_INIT_GRID_ITER", 10);
    const bool use_refinement_mode_for_initial_grid = toml::find_or<bool>(
        document, "BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE", true);
    if (initial_grid_iterations > 0 && !use_refinement_mode_for_initial_grid)
      throw std::runtime_error(
          "generated BH_WAMR solver requires "
          "BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE=true");
    const unsigned grain_size =
        toml::find_or<unsigned>(document, "BSSN_DENDRO_GRAIN_SZ", 1000);
    const DendroScalar load_imbalance_tolerance =
        toml::find_or<DendroScalar>(document, "BSSN_LOAD_IMB_TOL", 0.1);
    const unsigned split_fix =
        toml::find_or<unsigned>(document, "BSSN_SPLIT_FIX", 2);
    const unsigned diagnostic_frequency =
        toml::find_or<unsigned>(document, "BSSN_TIME_STEP_OUTPUT_FREQ", 25);
    const unsigned vtu_frequency =
        toml::find_or<unsigned>(document, "BSSN_IO_OUTPUT_FREQ", 8);
    const unsigned checkpoint_frequency =
        toml::find_or<unsigned>(document, "BSSN_CHECKPT_FREQ", 100);
    const unsigned apparent_horizon_frequency =
        toml::find_or<unsigned>(document, "AEH_SOLVER_FREQ", 0);
    const unsigned gravitational_wave_frequency =
        toml::find_or<unsigned>(document, "BSSN_GW_EXTRACT_FREQ", 0);
    const std::vector<DendroScalar> gravitational_wave_radii =
        toml::find_or<std::vector<DendroScalar>>(
            document, "BSSN_GW_RADAII", std::vector<DendroScalar>{{50.0}});
    const unsigned gravitational_wave_num_radii = toml::find_or<unsigned>(
        document, "BSSN_GW_NUM_RADAII",
        static_cast<unsigned>(gravitational_wave_radii.size()));
    const std::vector<unsigned> gravitational_wave_l_modes =
        toml::find_or<std::vector<unsigned>>(
            document, "BSSN_GW_L_MODES", std::vector<unsigned>{{2}});
    const unsigned gravitational_wave_num_l_modes = toml::find_or<unsigned>(
        document, "BSSN_GW_NUM_LMODES",
        static_cast<unsigned>(gravitational_wave_l_modes.size()));
    if (gravitational_wave_num_radii != gravitational_wave_radii.size() ||
        gravitational_wave_num_l_modes != gravitational_wave_l_modes.size() ||
        gravitational_wave_l_modes.empty() ||
        gravitational_wave_radii.empty())
      throw std::runtime_error("inconsistent gravitational-wave extraction arrays");
    if (wavelet_tolerance_mode == 6 &&
        (!(gravitational_wave_radii.front() > 8.0) ||
         !(gravitational_wave_radii.back() >= gravitational_wave_radii.front())))
      throw std::runtime_error("invalid mode-6 wavelet radial interval");
    const unsigned gravitational_wave_maximum_l =
        *std::max_element(gravitational_wave_l_modes.begin(),
                          gravitational_wave_l_modes.end());
    const unsigned nyquist_mode =
        toml::find_or<unsigned>(document, "BSSN_NYQUIST_M", 0);
    if (gravitational_wave_maximum_l < 2 || gravitational_wave_maximum_l > 8)
      throw std::runtime_error("BSSN_GW_L_MODES must lie in [2, 8]");
    const bool restore_solver =
        toml::find_or<unsigned>(document, "BSSN_RESTORE_SOLVER", 0) != 0;
    const std::string checkpoint_prefix = toml::find_or<std::string>(
        document, "BSSN_CHKPT_FILE_PREFIX", """
        + f'"cp/{profile_name}"'
        + r""");
    int checkpoint_index = -1;
    std::filesystem::file_time_type newest_checkpoint_time{};
    if (restore_solver) {
      for (unsigned index = 0; index < 2; ++index) {
        const std::filesystem::path metadata =
            checkpoint_prefix + "_" + std::to_string(index) + "_step.cp";
        if (!std::filesystem::exists(metadata)) continue;
        const auto write_time = std::filesystem::last_write_time(metadata);
        if (checkpoint_index < 0 || write_time > newest_checkpoint_time) {
          checkpoint_index = static_cast<int>(index);
          newest_checkpoint_time = write_time;
        }
      }
    }
    const DendroScalar grid_min_x =
        toml::find_or<DendroScalar>(document, "BSSN_GRID_MIN_X", -400.0);
    const DendroScalar grid_min_y =
        toml::find_or<DendroScalar>(document, "BSSN_GRID_MIN_Y", -400.0);
    const DendroScalar grid_min_z =
        toml::find_or<DendroScalar>(document, "BSSN_GRID_MIN_Z", -400.0);
    const DendroScalar grid_max_x =
        toml::find_or<DendroScalar>(document, "BSSN_GRID_MAX_X", 400.0);
    const DendroScalar grid_max_y =
        toml::find_or<DendroScalar>(document, "BSSN_GRID_MAX_Y", 400.0);
    const DendroScalar grid_max_z =
        toml::find_or<DendroScalar>(document, "BSSN_GRID_MAX_Z", 400.0);
    if (!(wavelet_tolerance > 0.0) || !std::isfinite(wavelet_tolerance) ||
        !(maximum_wavelet_tolerance > 0.0) ||
        !std::isfinite(maximum_wavelet_tolerance) ||
        !(gravitational_wave_tolerance > 0.0) ||
        !std::isfinite(gravitational_wave_tolerance) ||
        !(amr_coarsening_factor > 0.0) ||
        !(amr_coarsening_factor <= 1.0) ||
        !std::isfinite(amr_coarsening_factor) || !(cfl > 0.0) ||
        !(postmerger_amr_coarsening_factor >= 0.0) ||
        !(postmerger_amr_coarsening_factor <= 1.0) ||
        !std::isfinite(postmerger_amr_coarsening_factor) ||
        !(time_end > time_begin) || !(grid_max_x > grid_min_x) ||
        !(grid_max_y > grid_min_y) || !(grid_max_z > grid_min_z))
      throw std::runtime_error("invalid runtime parameter range");

    commondata_struct commondata{};
    commondata.NUMGRIDS = 1;
    const DendroScalar mass_1 =
        toml::find_or<DendroScalar>(document, "BSSN_BH1", "MASS", 0.5);
    const DendroScalar mass_2 =
        toml::find_or<DendroScalar>(document, "BSSN_BH2", "MASS", 0.5);
    const std::array<DendroScalar, 2> black_hole_masses{{mass_1, mass_2}};
    const std::array<DendroScalar, 2> black_hole_amr_radii{{
        toml::find_or<DendroScalar>(document, "BSSN_BH1_AMR_R", 2.0),
        toml::find_or<DendroScalar>(document, "BSSN_BH2_AMR_R", 2.0)}};
    const std::array<unsigned, 2> black_hole_maximum_levels{{
        toml::find_or<unsigned>(document, "BSSN_BH1_MAX_LEV", maximum_depth),
        toml::find_or<unsigned>(document, "BSSN_BH2_MAX_LEV", maximum_depth)}};
    const DendroScalar black_hole_amr_ratio =
        toml::find_or<DendroScalar>(document, "BSSN_AMR_R_RATIO", 2.0);
    if (grain_size == 0 || !(load_imbalance_tolerance >= 0.0) ||
        !std::isfinite(load_imbalance_tolerance) ||
        !(black_hole_amr_radii[0] > 0.0) ||
        !(black_hole_amr_radii[1] > 0.0) ||
        !std::isfinite(black_hole_amr_radii[0]) ||
        !std::isfinite(black_hole_amr_radii[1]) ||
        !(black_hole_amr_ratio > 1.0) ||
        !std::isfinite(black_hole_amr_ratio) ||
        black_hole_maximum_levels[0] < MAXDEAPTH_LEVEL_DIFF + 2 ||
        black_hole_maximum_levels[1] < MAXDEAPTH_LEVEL_DIFF + 2 ||
        black_hole_maximum_levels[0] > maximum_depth ||
        black_hole_maximum_levels[1] > maximum_depth)
      throw std::runtime_error("invalid mesh-adaptation parameters");
    const std::array<Point, 2> excision_centers{{
        Point(toml::find_or<DendroScalar>(document, "BSSN_BH1", "X", 4.0),
              toml::find_or<DendroScalar>(document, "BSSN_BH1", "Y", 0.0),
              toml::find_or<DendroScalar>(document, "BSSN_BH1", "Z", 0.0)),
        Point(toml::find_or<DendroScalar>(document, "BSSN_BH2", "X", -4.0),
              toml::find_or<DendroScalar>(document, "BSSN_BH2", "Y", 0.0),
              toml::find_or<DendroScalar>(document, "BSSN_BH2", "Z", 0.0))}};
    const std::array<Point, 2> initial_black_hole_velocities{{
        Point(toml::find_or<DendroScalar>(document, "BSSN_BH1", "V_X", 0.0),
              toml::find_or<DendroScalar>(document, "BSSN_BH1", "V_Y", 0.0),
              toml::find_or<DendroScalar>(document, "BSSN_BH1", "V_Z", 0.0)),
        Point(toml::find_or<DendroScalar>(document, "BSSN_BH2", "V_X", 0.0),
              toml::find_or<DendroScalar>(document, "BSSN_BH2", "V_Y", 0.0),
              toml::find_or<DendroScalar>(document, "BSSN_BH2", "V_Z", 0.0))}};
    const std::array<DendroScalar, 2> excision_radii{{
        toml::find_or<DendroScalar>(document, "BSSN_BH1_CONSTRAINT_R", 1.0),
        toml::find_or<DendroScalar>(document, "BSSN_BH2_CONSTRAINT_R", 1.0)}};
    commondata.mass_ratio = std::max(mass_1, mass_2) / std::min(mass_1, mass_2);
    commondata.initial_sep =
        2.0 * toml::find_or<DendroScalar>(document, "TPID_PAR_B", 4.0);
    commondata.initial_p_t = std::abs(toml::find_or<DendroScalar>(
        document, "BSSN_BH1", "V_Y", 0.11284523509709575));
    commondata.initial_p_r = std::abs(toml::find_or<DendroScalar>(
        document, "BSSN_BH1", "V_X", -0.002284343811437988));
    commondata.bbhxy_BH_m_chix = commondata.bbhxy_BH_m_chiy =
        commondata.bbhxy_BH_m_chiz = 0.0;
    commondata.bbhxy_BH_M_chix = commondata.bbhxy_BH_M_chiy =
        commondata.bbhxy_BH_M_chiz = 0.0;
    commondata.TP_npoints_A =
        toml::find_or<int>(document, "TPID_NPOINTS_A", 65);
    commondata.TP_npoints_B =
        toml::find_or<int>(document, "TPID_NPOINTS_B", 78);
    commondata.TP_npoints_phi =
        toml::find_or<int>(document, "TPID_NPOINTS_PHI", 10);
    if (toml::find_or<int>(document, "TPID_GIVE_BARE_MASS", 1) != 0) {
      commondata.TP_bare_mass_M = toml::find_or<DendroScalar>(
          document, "TPID_TARGET_M_PLUS", 0.48236442246752931);
      commondata.TP_bare_mass_m = toml::find_or<DendroScalar>(
          document, "TPID_TARGET_M_MINUS", 0.48236442246752931);
    } else {
      commondata.TP_bare_mass_M = commondata.TP_bare_mass_m = -1.0;
    }
    std::snprintf(commondata.TP_BBH_description,
                  sizeof(commondata.TP_BBH_description),
                  "q=%g TwoPunctures", commondata.mass_ratio);
    params_struct tp_params{};
    ID_persist_struct punctures{};
    initialize_ID_persist_struct(&commondata, &punctures);
    punctures.Newton_tol = toml::find_or<DendroScalar>(
        document, "TPID_NEWTON_TOL", punctures.Newton_tol);
    punctures.adm_tol = toml::find_or<DendroScalar>(
        document, "TPID_ADM_TOL", punctures.adm_tol);
    punctures.initial_lapse_psi_exponent = -2.0;
    std::snprintf(punctures.initial_lapse, sizeof(punctures.initial_lapse),
                  "W");

    // Compare only inputs to the puncture solve, before target-mass iteration.
    const std::array<REAL, 33> tpid_inputs{{
        static_cast<REAL>(punctures.npoints_A),
        static_cast<REAL>(punctures.npoints_B),
        static_cast<REAL>(punctures.npoints_phi),
        static_cast<REAL>(punctures.Newton_maxit),
        punctures.adm_tol, punctures.Newton_tol, punctures.TP_epsilon,
        punctures.TP_Tiny, punctures.TP_Extend_Radius, punctures.par_b,
        punctures.par_m_plus, punctures.par_m_minus,
        punctures.target_M_plus, punctures.target_M_minus,
        punctures.par_P_plus[0], punctures.par_P_plus[1],
        punctures.par_P_plus[2], punctures.par_P_minus[0],
        punctures.par_P_minus[1], punctures.par_P_minus[2],
        punctures.par_S_plus[0], punctures.par_S_plus[1],
        punctures.par_S_plus[2], punctures.par_S_minus[0],
        punctures.par_S_minus[1], punctures.par_S_minus[2],
        punctures.center_offset[0], punctures.center_offset[1],
        punctures.center_offset[2],
        static_cast<REAL>(punctures.give_bare_mass),
        static_cast<REAL>(punctures.use_sources),
        static_cast<REAL>(punctures.rescale_sources),
        static_cast<REAL>(punctures.solve_momentum_constraint)}};
    const std::string tpid_prefix =
        toml::find_or<std::string>(document, "TPID_FILEPREFIX", "tp");
    if (tpid_prefix.empty())
      throw std::runtime_error("TPID_FILEPREFIX must not be empty");
    const std::filesystem::path tpid_file =
        tpid_prefix + "_nrpy_tpid_sol.bin";
    if (generate_tpid || checkpoint_index < 0) {
      const std::filesystem::path output_file =
          generate_tpid ? std::filesystem::path(tpid_file.string() + ".tmp")
                        : tpid_file;
      std::fstream coefficients(
          output_file, std::ios::binary |
                           (generate_tpid ? std::ios::out | std::ios::trunc
                                          : std::ios::in));
      if (!coefficients)
        throw std::runtime_error("cannot open NRPy TwoPunctures file: " +
                                 output_file.string());
      constexpr std::array<char, 16> file_tag{{'N', 'R', 'P', 'y', '-', 'T',
                                               'P', 'I', 'D', '-', '1'}};
      if (generate_tpid) {
        TP_solve(&punctures);
        coefficients.write(file_tag.data(), file_tag.size());
        coefficients.write(reinterpret_cast<const char*>(tpid_inputs.data()),
                           sizeof(tpid_inputs));
        const std::array<REAL, 8> results{{
            punctures.mp, punctures.mm, punctures.mp_adm, punctures.mm_adm,
            punctures.E, punctures.J1, punctures.J2, punctures.J3}};
        coefficients.write(reinterpret_cast<const char*>(results.data()),
                           sizeof(results));
      } // END IF: write puncture metadata
      else {
        std::array<char, 16> stored_tag{};
        std::array<REAL, 33> stored_inputs{};
        std::array<REAL, 8> results{};
        coefficients.read(stored_tag.data(), stored_tag.size());
        coefficients.read(reinterpret_cast<char*>(stored_inputs.data()),
                          sizeof(stored_inputs));
        coefficients.read(reinterpret_cast<char*>(results.data()),
                          sizeof(results));
        if (!coefficients || stored_tag != file_tag ||
            stored_inputs != tpid_inputs)
          throw std::runtime_error(
              "NRPy TwoPunctures file does not match the input parameters; "
              "run the solver with --tpid first: " + tpid_file.string());
        punctures.mp = results[0];
        punctures.mm = results[1];
        punctures.mp_adm = results[2];
        punctures.mm_adm = results[3];
        punctures.E = results[4];
        punctures.J1 = results[5];
        punctures.J2 = results[6];
        punctures.J3 = results[7];
        punctures.par_m_plus = punctures.mp;
        punctures.par_m_minus = punctures.mm;
        const std::size_t coefficient_count =
            static_cast<std::size_t>(punctures.npoints_A) *
            static_cast<std::size_t>(punctures.npoints_B) *
            static_cast<std::size_t>(punctures.npoints_phi);
        if (coefficient_count > static_cast<std::size_t>(
                                    std::numeric_limits<int>::max()))
          throw std::runtime_error("TwoPunctures coefficient count overflows int");
        allocate_derivs(&punctures.v, static_cast<int>(coefficient_count));
        allocate_derivs(&punctures.cf_v, static_cast<int>(coefficient_count));
      } // END ELSE: restore puncture metadata
      const std::streamsize coefficient_bytes =
          static_cast<std::streamsize>(punctures.npoints_A) *
          punctures.npoints_B * punctures.npoints_phi * sizeof(REAL);
      for (derivs* derivatives : {&punctures.v, &punctures.cf_v}) {
        for (REAL* values : {derivatives->d0, derivatives->d1,
                             derivatives->d2, derivatives->d3,
                             derivatives->d11, derivatives->d12,
                             derivatives->d13, derivatives->d22,
                             derivatives->d23, derivatives->d33}) {
          if (generate_tpid)
            coefficients.write(reinterpret_cast<const char*>(values),
                               coefficient_bytes);
          else
            coefficients.read(reinterpret_cast<char*>(values),
                              coefficient_bytes);
        } // END LOOP: puncture derivative arrays
      } // END LOOP: puncture coefficient families
      if (!coefficients)
        throw std::runtime_error("incomplete NRPy TwoPunctures file: " +
                                 output_file.string());
      if (generate_tpid) {
        coefficients.close();
        if (!coefficients)
          throw std::runtime_error("cannot finish NRPy TwoPunctures file");
        std::filesystem::rename(output_file, tpid_file);
        std::cout << "wrote NRPy TwoPunctures data: " << tpid_file << '\n';
      } // END IF: finalize puncture file
      else {
        if (coefficients.peek() != std::char_traits<char>::eof())
          throw std::runtime_error("NRPy TwoPunctures file has extra data");
        if (rank == 0)
          std::cout << "loaded NRPy TwoPunctures data: " << tpid_file << '\n';
      } // END ELSE: validate loaded puncture file
    } // END IF: produce or load punctures
    if (generate_tpid) {
      MPI_Finalize();
      return 0;
    } // END IF: exit after puncture solve

    const Point domain_minimum(grid_min_x, grid_min_y, grid_min_z);
    const Point domain_maximum(grid_max_x, grid_max_y, grid_max_z);
    m_uiMaxDepth = maximum_depth;
    _InitializeHcurve(m_uiDim);
    std::vector<ot::TreeNode> octree;
    const DendroScalar octree_coordinate_scale =
        std::ldexp(1.0, -static_cast<int>(maximum_depth));
    const double initial_mesh_chi_floor = toml::find_or<DendroScalar>(
        document, "CHI_FLOOR", 0.1);
    const nrpy_dendro_seed::PunctureParameters seed_black_hole_1{
        mass_1, excision_centers[0].x(), excision_centers[0].y(),
        excision_centers[0].z(), initial_black_hole_velocities[0].x(),
        initial_black_hole_velocities[0].y(),
        initial_black_hole_velocities[0].z(),
        toml::find_or<DendroScalar>(document, "BSSN_BH1", "SPIN", 0.0),
        toml::find_or<DendroScalar>(document, "BSSN_BH1", "SPIN_THETA", 0.0),
        toml::find_or<DendroScalar>(document, "BSSN_BH1", "SPIN_PHI", 0.0)};
    const nrpy_dendro_seed::PunctureParameters seed_black_hole_2{
        mass_2, excision_centers[1].x(), excision_centers[1].y(),
        excision_centers[1].z(), initial_black_hole_velocities[1].x(),
        initial_black_hole_velocities[1].y(),
        initial_black_hole_velocities[1].z(),
        toml::find_or<DendroScalar>(document, "BSSN_BH2", "SPIN", 0.0),
        toml::find_or<DendroScalar>(document, "BSSN_BH2", "SPIN_THETA", 0.0),
        toml::find_or<DendroScalar>(document, "BSSN_BH2", "SPIN_PHI", 0.0)};
    std::function<void(double, double, double, double*)> initial_bssn_fields =
        [&](double x, double y, double z, double* fields) {
          const double physical_x = grid_min_x + x * octree_coordinate_scale *
                                                   (grid_max_x - grid_min_x);
          const double physical_y = grid_min_y + y * octree_coordinate_scale *
                                                   (grid_max_y - grid_min_y);
          const double physical_z = grid_min_z + z * octree_coordinate_scale *
                                                   (grid_max_z - grid_min_z);
          nrpy_dendro_seed::punctureDataPhysicalCoord(
              physical_x, physical_y, physical_z, fields, seed_black_hole_1,
              seed_black_hole_2, initial_mesh_chi_floor);
        };
    std::array<unsigned, 24> initial_field_indices{};
    for (unsigned field = 0; field < initial_field_indices.size(); ++field)
      initial_field_indices[field] = field;
    const unsigned initial_refinement_depth =
        std::min(black_hole_maximum_levels[0],
                 black_hole_maximum_levels[1]) - MAXDEAPTH_LEVEL_DIFF - 2;
    if (checkpoint_index >= 0)
      createRegularOctree(octree, minimum_depth, m_uiDim, maximum_depth,
                          MPI_COMM_WORLD);
    else
      function2Octree(initial_bssn_fields,
                      initial_field_indices.size(),
                      initial_field_indices.data(),
                      initial_field_indices.size(), octree,
                      initial_refinement_depth, wavelet_tolerance,
                      element_order, MPI_COMM_WORLD);
    ot::Mesh* mesh = ot::createMesh(
        octree.data(), octree.size(), element_order, MPI_COMM_WORLD, 1,
        ot::SM_TYPE::FDM, grain_size, load_imbalance_tolerance, split_fix);
    if (!mesh) throw std::runtime_error("mesh construction failed");
    mesh->setDomainBounds(domain_minimum, domain_maximum);
    unsigned local_minimum_depth = 0, local_maximum_depth = 0;
    mesh->computeMinMaxLevel(local_minimum_depth, local_maximum_depth);
    DendroScalar minimum_dx = std::ldexp(
        (grid_max_x - grid_min_x) / element_order,
        -static_cast<int>(local_maximum_depth));
    DendroScalar time_step = cfl * minimum_dx;

    """
        + solver_namespace
        + r"""::generated::params_struct params{};
    """
        + solver_stem
        + r"""_params_struct_set_to_default(params);
    for (const auto& item : document.as_table()) {
"""
        + output_toml_bindings()
        + r"""
    }
    """
        + solver_stem
        + r"""_params_validate(params);
    const std::string output_prefix = toml::find_or<std::string>(
        document, "BSSN_PROFILE_FILE_PREFIX", """
        + f'"{profile_name}"'
        + r""");
    const std::filesystem::path output_path(output_prefix);
    if (!output_path.parent_path().empty())
      std::filesystem::create_directories(output_path.parent_path());
    const std::string vtu_prefix = toml::find_or<std::string>(
        document, "BSSN_VTU_FILE_PREFIX", """
        + f'"vtu/{profile_name}"'
        + r""");
    const std::filesystem::path vtu_path(vtu_prefix);
    if (!vtu_path.parent_path().empty())
      std::filesystem::create_directories(vtu_path.parent_path());
    const std::filesystem::path checkpoint_path(checkpoint_prefix);
    if (!checkpoint_path.parent_path().empty())
      std::filesystem::create_directories(checkpoint_path.parent_path());
    std::unique_ptr<dendro_aeh::AEH_BHaHAHA> apparent_horizon_finder;
    if (apparent_horizon_frequency > 0) {
      const std::string horizon_directory = toml::find_or<std::string>(
          document, "AEH_PARAMS", "AEH_SAVE_DIR", "bah");
      std::filesystem::create_directories(horizon_directory);
      const std::vector<double> initial_x{
          excision_centers[0].x(), excision_centers[1].x(), 0.0};
      const std::vector<double> initial_y{
          excision_centers[0].y(), excision_centers[1].y(), 0.0};
      const std::vector<double> initial_z{
          excision_centers[0].z(), excision_centers[1].z(), 0.0};
      const std::vector<dendro_aeh::SimpleBlackHoleData> black_holes{
          dendro_aeh::SimpleBlackHoleData(
              excision_centers[0].x(), excision_centers[0].y(),
              excision_centers[0].z(), mass_1),
          dendro_aeh::SimpleBlackHoleData(
              excision_centers[1].x(), excision_centers[1].y(),
              excision_centers[1].z(), mass_2)};
      const std::vector<int> horizon_indices{"""
        + horizon_indices
        + r"""};
      const auto bssn_to_adm = [](const std::vector<double>& values) {
        const double inverse_chi = """
        + inverse_chi_expression
        + r""";
        const double one_third_trK = values[1] / 3.0;
        std::vector<double> adm(12);
        adm[0] = (values[8] + 1.0) * inverse_chi;
        adm[1] = values[9] * inverse_chi;
        adm[2] = values[10] * inverse_chi;
        adm[3] = (values[11] + 1.0) * inverse_chi;
        adm[4] = values[12] * inverse_chi;
        adm[5] = (values[13] + 1.0) * inverse_chi;
        for (unsigned component = 0; component < 6; ++component)
          adm[6 + component] = values[2 + component] * inverse_chi +
                               adm[component] * one_third_trK;
        return adm;
      };
      const Point grid_limits[2]{
          Point(0.0, 0.0, 0.0),
          Point(1u << maximum_depth, 1u << maximum_depth,
                1u << maximum_depth)};
      const Point domain_limits[2]{domain_minimum, domain_maximum};
      apparent_horizon_finder = std::make_unique<dendro_aeh::AEH_BHaHAHA>(
          3, true, initial_x, initial_y, initial_z, 3,
          std::vector<double>{0.0, 0.0, 0.0},
          toml::find_or<std::vector<double>>(
              document, "AEH_PARAMS", "CFL_FACTOR",
              std::vector<double>{1.0, 1.0, 1.0}),
          std::vector<int>{10000, 10000, 10000},
          toml::find_or<std::vector<double>>(
              document, "AEH_PARAMS", "THETA_L2_M_TOL",
              std::vector<double>{1.0e-5, 1.0e-5, 1.0e-5}),
          toml::find_or<std::vector<double>>(
              document, "AEH_PARAMS", "THETA_LINF_M_TOL",
              std::vector<double>{1.0e-2, 1.0e-2, 1.0e-2}),
          std::vector<double>{7.0, 7.0, 7.0},
          std::vector<double>{0.0, 0.0, 0.0},
          toml::find_or<std::vector<double>>(
              document, "AEH_PARAMS", "MAX_SEARCH_RADIUS",
              std::vector<double>{1.5, 1.5, 1.5}),
          toml::find_or<std::vector<int>>(
              document, "AEH_PARAMS", "NR_INTERP_MAX",
              std::vector<int>{48, 48, 48}),
          32, 64, horizon_directory, black_holes, horizon_indices,
          bssn_to_adm, grid_limits, domain_limits, vtu_frequency, 3,
          std::vector<int>{8, 16, 32}, std::vector<int>{16, 32, 64}, 0,
          toml::find_or<int>(document, "AEH_PARAMS", "VERBOSITY_LEVEL", 1));
    }
    {
    """
        + solver_namespace
        + r"""::Ctx context(
        mesh, domain_minimum, domain_maximum, time_step, wavelet_tolerance,
        wavelet_tolerance_mode, maximum_wavelet_tolerance,
        gravitational_wave_tolerance, amr_coarsening_factor,
        postmerger_amr_coarsening_factor, refinement_variables,
        remesh_frequency, postmerger_remesh_frequency,
        diagnostic_frequency, vtu_frequency,
        checkpoint_frequency, output_prefix, vtu_prefix, checkpoint_prefix,
        excision_centers, excision_radii, black_hole_masses,
        black_hole_amr_radii, black_hole_maximum_levels,
        black_hole_amr_ratio, minimum_depth, apparent_horizon_frequency,
        apparent_horizon_finder.get(), gravitational_wave_frequency,
        gravitational_wave_radii, gravitational_wave_maximum_l, nyquist_mode,
        initial_black_hole_velocities, time_begin,
        Point(0.0, 0.0, 0.0));
    context.params = params;
    ts::TSInfo time_info{};
    time_info._m_uiStep = 0;
    time_info._m_uiT = time_begin;
    time_info._m_uiTh = time_step;
    context.set_ts_info(time_info);
    bool restored = false;
    if (checkpoint_index >= 0) {
      if (context.restore_checkpt(static_cast<unsigned>(checkpoint_index)) != 0)
        throw std::runtime_error("checkpoint restore failed");
      mesh = context.get_mesh();
      restored = true;
    }
    if (!restored) {
      if (context.initialize(commondata, tp_params, punctures) != 0)
        throw std::runtime_error("initial-data construction failed");
      const unsigned initial_grid_remesh_passes =
          initial_grid_iterations > 1 ? initial_grid_iterations - 1
                                      : initial_grid_iterations;
      std::array<unsigned long long, 2> local_grid_counts{
          context.get_mesh()->getNumLocalMeshElements(),
          context.get_mesh()->getNumLocalMeshNodes()};
      std::array<unsigned long long, 2> global_grid_counts{};
      MPI_Allreduce(local_grid_counts.data(), global_grid_counts.data(),
                    2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      bool initial_grid_remeshed = false;
      for (unsigned pass = 0; pass < initial_grid_remesh_passes; ++pass) {
        if (!context.is_remesh(true)) break;
        context.remesh_and_gridtransfer(
            grain_size, load_imbalance_tolerance, split_fix);
        initial_grid_remeshed = true;
        const auto old_grid_counts = global_grid_counts;
        local_grid_counts = {
            context.get_mesh()->getNumLocalMeshElements(),
            context.get_mesh()->getNumLocalMeshNodes()};
        MPI_Allreduce(local_grid_counts.data(), global_grid_counts.data(),
                      2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (rank == 0)
          std::cout << "initial-grid remesh pass=" << pass + 1
                    << " elements=" << old_grid_counts[0] << "->"
                    << global_grid_counts[0] << " nodes="
                    << old_grid_counts[1] << "->" << global_grid_counts[1]
                    << std::endl;
        if (global_grid_counts == old_grid_counts) break;
      }
      if (initial_grid_remeshed &&
          context.initialize(commondata, tp_params, punctures) != 0)
        throw std::runtime_error("initial-grid reconstruction failed");
      unsigned converged_minimum_depth = 0, converged_maximum_depth = 0;
      context.get_mesh()->computeMinMaxLevel(converged_minimum_depth,
                                             converged_maximum_depth);
      minimum_dx = std::ldexp(
          (grid_max_x - grid_min_x) / element_order,
          -static_cast<int>(converged_maximum_depth));
      time_step = cfl * minimum_dx;
      ts::TSInfo converged_time_info = context.get_ts_info();
      converged_time_info._m_uiTh = time_step;
      context.set_ts_info(converged_time_info);
    } else {
      unsigned restored_minimum_depth = 0, restored_maximum_depth = 0;
      context.get_mesh()->computeMinMaxLevel(restored_minimum_depth,
                                             restored_maximum_depth);
      minimum_dx = std::ldexp(
          (grid_max_x - grid_min_x) / element_order,
          -static_cast<int>(restored_maximum_depth));
      time_step = context.get_ts_info()._m_uiTh;
    }
    ts::ETS<DendroScalar, """
        + solver_namespace
        + r"""::Ctx> time_stepper(&context);
    time_stepper.set_ets_coefficients(ts::ETSType::RK4);
    time_stepper.set_evolve_vars(context.get_evolution_vars());
    time_stepper.init();
    if (context.diagnostic_output() != 0)
      throw std::runtime_error("initial diagnostic output failed");
    if (context.gravitational_wave_output() != 0)
      throw std::runtime_error("initial gravitational-wave output failed");
    if (context.adm_output() != 0)
      throw std::runtime_error("initial ADM output failed");
    if (context.write_vtu() != 0)
      throw std::runtime_error("initial VTU output failed");
    if (context.apparent_horizon_output() != 0)
      throw std::runtime_error("initial apparent-horizon output failed");
    if (rank == 0)
      std::cout << """
        + f'"{executable_name}: profile={profile_name} FD="'
        + r""" << element_order
                << " KO=" << element_order - 2
                << " dx_min=" << minimum_dx << " dt=" << time_step << '\n';
    while (time_stepper.curr_time() < time_end &&
           time_stepper.curr_step() < maximum_iterations) {
      context.terminal_output();
      time_stepper.evolve();
      if (context.is_remesh()) {
        context.remesh_and_gridtransfer(
            grain_size, load_imbalance_tolerance, split_fix);
        context.post_timestep(context.get_evolution_vars());
        time_stepper.sync_with_mesh();
        unsigned remeshed_minimum_depth = 0, remeshed_maximum_depth = 0;
        context.get_mesh()->computeMinMaxLevel(remeshed_minimum_depth,
                                               remeshed_maximum_depth);
        const DendroScalar remeshed_minimum_dx = std::ldexp(
            (grid_max_x - grid_min_x) / element_order,
            -static_cast<int>(remeshed_maximum_depth));
        ts::TSInfo remeshed_time_info = context.get_ts_info();
        remeshed_time_info._m_uiTh = cfl * remeshed_minimum_dx;
        context.set_ts_info(remeshed_time_info);
      }
      if (context.evolve_excision_centers() != 0)
        throw std::runtime_error("puncture-center evolution failed");
      if (context.diagnostic_output() != 0)
        throw std::runtime_error("diagnostic output failed");
      if (context.gravitational_wave_output() != 0)
        throw std::runtime_error("gravitational-wave output failed");
      if (context.adm_output() != 0)
        throw std::runtime_error("ADM output failed");
      if (context.write_vtu() != 0)
        throw std::runtime_error("VTU output failed");
      if (context.apparent_horizon_output() != 0)
        throw std::runtime_error("apparent-horizon output failed");
      if (context.write_checkpt() != 0)
        throw std::runtime_error("checkpoint write failed");
    }
    context.terminal_output();
    mesh = context.get_mesh();
    }
    delete mesh;
  } catch (const std::exception& error) {
    if (rank == 0) std::cerr << """
        + f'"{executable_name}: "'
        + r""" << error.what() << '\n';
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  MPI_Finalize();
  return 0;
}
"""
    )
