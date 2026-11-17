"""
Generate and run a BHaH project checking spheroidal fisheye inversion.

Author: Nishita Jadoo
"""

import argparse
import contextlib
import doctest
import io
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict, List

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures.BHaH.BHaH_defines_h import output_BHaH_defines_h
from nrpy.infrastructures.BHaH.generalrfm_cart_to_xx import (
    register_CFunction_generalrfm_Cart_to_xx,
)
from nrpy.infrastructures.BHaH.Makefile_helpers import (
    output_CFunctions_function_prototypes_and_construct_Makefile,
)
from nrpy.infrastructures.BHaH.xx_tofrom_Cart import register_CFunction_xx_to_Cart

PROJECT_NAME = "spheroidal_fisheye_roundtrip"
VECTOR_CHECK = r"""
static int vectors_match(const REAL expected[3], const REAL actual[3]) {
  const REAL tolerance = sizeof(REAL) == sizeof(float) ? (REAL)1.0e-4 : (REAL)1.0e-10;
  for (int i = 0; i < 3; i++) {
    const REAL scale = fmax((REAL)1.0, fabs(expected[i]));
    if (!isfinite(actual[i]) || fabs(expected[i] - actual[i]) > tolerance * scale)
      return 0;
  } // END LOOP: vector components
  return 1;
} // END FUNCTION: vectors_match

static int check_roundtrip(const params_struct *params, const REAL xx[3],
                           void (*forward)(const params_struct *, const REAL[3], REAL[3]),
                           int (*inverse)(const params_struct *, const REAL[3], REAL[3])) {
  REAL Cart[3], xx_back[3];
  forward(params, xx, Cart);
  if (inverse(params, Cart, xx_back) != 0)
    return 0;
  return vectors_match(xx, xx_back);
} // END FUNCTION: check_roundtrip

static int check_small_radius_inverse(const params_struct *params,
                                     int (*inverse)(const params_struct *, const REAL[3], REAL[3])) {
  const int exponents[] = {8, 12, 16, 20, 30, 200};
  const REAL directions[3][3] = {{1, 0, 0}, {0, 0, 1}, {0.6, -0.8, 1}};
  const REAL tolerance = sizeof(REAL) == sizeof(float) ? (REAL)1.0e-4 : (REAL)1.0e-10;
  const int num_radii = sizeof(REAL) == sizeof(float) ? 5 : 6;
  for (int n = 0; n < num_radii; n++) {
    const REAL radius = (REAL)pow(10.0, -exponents[n]);
    for (int direction = 0; direction < 3; direction++) {
      REAL Cart[3], expected[3], actual[3];
      for (int i = 0; i < 3; i++) {
        Cart[i] = radius * directions[direction][i];
        const REAL scale = i == 2 ? params->fisheye_z_c : params->fisheye_xy_c;
        expected[i] = Cart[i] / (params->fisheye_a0 * scale);
      } // END LOOP: limiting inverse components
      const int status = inverse(params, Cart, actual);
      const REAL expected_scale = fmax(fmax(fabs(expected[0]), fabs(expected[1])), fabs(expected[2]));
      for (int i = 0; i < 3; i++) {
        if (status != 0 || !isfinite(actual[i]) ||
            fabs((actual[i] - expected[i]) / expected_scale) > tolerance) {
          fprintf(stderr, "Small-radius inverse failed: exponent=%d direction=%d status=%d component=%d expected=%.17g actual=%.17g\n",
                  exponents[n], direction, status, i, (double)expected[i], status ? 0.0 : (double)actual[i]);
          return 0;
        } // END IF: limiting inverse mismatch
      } // END LOOP: check relative errors
    } // END LOOP: Cartesian directions
  } // END LOOP: small nonzero radii
  return 1;
} // END FUNCTION: check_small_radius_inverse
"""

N8_PARAMETERS: Dict[str, List[float]] = {
    "fisheye_a": [1, 2, 4, 8, 16, 32, 64, 128, 256],
    "fisheye_xy_R": [
        171.06625518081691,
        184.44748588094367,
        215.80453713804999,
        246.99533671450504,
        286.33984889444974,
        315.21952800621472,
        345.56945122994739,
        371.01972035005042,
    ],
    "fisheye_xy_s": [
        11.383190952826789,
        6.3566576581017591,
        6.144116855070628,
        6.1420873714876416,
        9.2010079885841947,
        7.6806497869980452,
        7.657345408705031,
        7.6860700359661864,
    ],
    "fisheye_z_R": [
        205.49362806002722,
        220.78389866486555,
        223.95116266565509,
        248.54941568121561,
        279.15394051173627,
        307.36671974291738,
        333.5953870005842,
        355.59171204875031,
    ],
    "fisheye_z_s": [
        13.477846089992086,
        5.4664273952348825,
        5.3509689956296613,
        4.778580730986703,
        7.1649852613864402,
        6.633355484239388,
        6.6181317532559625,
        6.6429682390042704,
    ],
    "fisheye_xy_c": [0.027126574950438669],
    "fisheye_z_c": [0.022597979833417235],
}


def _run(command: List[str], cwd: Path, log_path: Path, timeout: int) -> None:
    """
    Run a generated-project command with bounded resources and output.

    :param command: Command and arguments.
    :param cwd: Generated project directory.
    :param log_path: File receiving command output.
    :param timeout: Maximum duration in seconds.
    :raises RuntimeError: If the command fails.
    """
    with log_path.open("w", encoding="utf-8") as log_file:
        result = subprocess.run(
            command,
            cwd=cwd,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    if result.returncode != 0:
        with log_path.open("rb") as log_file:
            log_file.seek(0, 2)
            log_file.seek(max(0, log_file.tell() - 4000))
            tail = log_file.read().decode("utf-8", errors="replace")
        raise RuntimeError(
            f"{' '.join(command)} failed with exit status {result.returncode}:\n{tail}"
        )


def main() -> None:
    """Build and run the generated BHaH coordinate check for one precision."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", choices=("float", "double"), required=True)
    precision = parser.parse_args().precision

    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("parallelization", "openmp")
    par.set_parval_from_str("fp_type", precision)
    cfc.CFunction_dict.clear()
    with contextlib.redirect_stdout(io.StringIO()):
        for transitions in (1, 8):
            coord = f"GeneralRFM_spheroidal_fisheyeN{transitions}"
            register_CFunction_xx_to_Cart(coord)
            register_CFunction_generalrfm_Cart_to_xx(coord)
    assignment_lines = []
    for prefix, entries in N8_PARAMETERS.items():
        for index, value in enumerate(entries):
            if prefix.endswith("_c"):
                parameter_name = prefix
            elif prefix == "fisheye_a":
                parameter_name = f"{prefix}{index}"
            else:
                parameter_name = f"{prefix}{index + 1}"
            assignment_lines.append(f"params.{parameter_name} = (REAL){value!r};")
    n8_assignments = "\n".join(assignment_lines)
    desc = """Check spheroidal fisheye inverse coordinates.

    @return 0 if all coordinate round trips pass; nonzero test-case code on failure.
    """
    cfunc_type = "int"
    name = "main"
    params = "void"
    body = (
        r"""
params_struct params = {0};
params.fisheye_a0 = (REAL)1.0;
params.fisheye_a1 = (REAL)256.0;
params.fisheye_xy_R1 = params.fisheye_z_R1 = (REAL)9.99;
params.fisheye_xy_s1 = params.fisheye_z_s1 = (REAL)0.01;
params.fisheye_xy_c = params.fisheye_z_c = (REAL)0.7866685970142996;
const REAL strongly_scaled_xx[3] = {(REAL)9.0, (REAL)9.0, (REAL)9.0};
if (!check_roundtrip(&params, strongly_scaled_xx,
                     xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN1,
                     generalrfm_Cart_to_xx__GeneralRFM_spheroidal_fisheyeN1))
  return 1;

params.fisheye_a1 = (REAL)16.0;
params.fisheye_xy_R1 = (REAL)4.0;
params.fisheye_xy_s1 = (REAL)0.5;
params.fisheye_xy_c = (REAL)0.8;
params.fisheye_z_R1 = (REAL)8.0;
params.fisheye_z_s1 = (REAL)1.2;
params.fisheye_z_c = (REAL)1.1;
if (!check_small_radius_inverse(&params, generalrfm_Cart_to_xx__GeneralRFM_spheroidal_fisheyeN1))
  return 10;
const REAL anisotropic_xx[3] = {(REAL)3.0, (REAL)-5.0, (REAL)7.0};
if (!check_roundtrip(&params, anisotropic_xx,
                     xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN1,
                     generalrfm_Cart_to_xx__GeneralRFM_spheroidal_fisheyeN1))
  return 2;

if (sizeof(REAL) == sizeof(double)) {
  params.fisheye_a0 = params.fisheye_a1 = (REAL)1.0;
  params.fisheye_xy_R1 = params.fisheye_z_R1 = (REAL)1.0;
  params.fisheye_xy_s1 = params.fisheye_z_s1 = (REAL)0.5;
  params.fisheye_xy_c = params.fisheye_z_c = (REAL)1.0;
  const REAL tiny_xx[3] = {(REAL)1.0e-15, (REAL)0.0, (REAL)0.0};
  REAL tiny_Cart[3], tiny_xx_back[3];
  xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN1(&params, tiny_xx, tiny_Cart);
  if (generalrfm_Cart_to_xx__GeneralRFM_spheroidal_fisheyeN1(
          &params, tiny_Cart, tiny_xx_back) != 0 ||
      fabs(tiny_xx_back[0] / tiny_xx[0] - (REAL)1.0) > (REAL)1.0e-10)
    return 9;
} // END IF: tiny identity radius

params.fisheye_a0 = params.fisheye_a1 = (REAL)1.0;
params.fisheye_xy_R1 = params.fisheye_z_R1 = (REAL)1.0;
params.fisheye_xy_s1 = params.fisheye_z_s1 = (REAL)0.5;
params.fisheye_xy_c = params.fisheye_z_c = (REAL)1.0e-30;
const REAL small_scale_xx[3] = {(REAL)1.0, (REAL)0.0, (REAL)0.0};
if (!check_roundtrip(&params, small_scale_xx,
                     xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN1,
                     generalrfm_Cart_to_xx__GeneralRFM_spheroidal_fisheyeN1))
  return 12;

if (sizeof(REAL) == sizeof(double)) {
  params.fisheye_a0 = (REAL)1.0e10;
  params.fisheye_a1 = (REAL)1.0;
  params.fisheye_xy_R1 = params.fisheye_z_R1 = (REAL)1.0;
  params.fisheye_xy_s1 = params.fisheye_z_s1 = (REAL)0.1;
  params.fisheye_xy_c = params.fisheye_z_c = (REAL)1.0;
  const REAL poorly_conditioned_xx[3] = {(REAL)100.0, (REAL)0.0, (REAL)0.0};
  REAL Cart[3], xx_back[3];
  xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN1(
      &params, poorly_conditioned_xx, Cart);
  const int status = generalrfm_Cart_to_xx__GeneralRFM_spheroidal_fisheyeN1(
      &params, Cart, xx_back);
  if (status == 0 && !vectors_match(poorly_conditioned_xx, xx_back))
    return 3;
} // END IF: double precision conditioning

"""
        + n8_assignments
        + r"""
if (!check_small_radius_inverse(&params, generalrfm_Cart_to_xx__GeneralRFM_spheroidal_fisheyeN8))
  return 11;
const REAL shipped_n8_xx[3] = {
    (REAL)-192.3581085205078,
    (REAL)-94.86515808105469,
    (REAL)127.530517578125,
};
if (!check_roundtrip(&params, shipped_n8_xx,
                     xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN8,
                     generalrfm_Cart_to_xx__GeneralRFM_spheroidal_fisheyeN8))
  return 4;
if (sizeof(REAL) == sizeof(double)) {
  const REAL stalled_n8_xx[3] = {
      (REAL)58.90804597701151,
      (REAL)-61.781609195402325,
      (REAL)234.375,
  };
  if (!check_roundtrip(&params, stalled_n8_xx,
                       xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN8,
                       generalrfm_Cart_to_xx__GeneralRFM_spheroidal_fisheyeN8))
    return 5;
} // END IF: double precision N8 inversion
if (sizeof(REAL) == sizeof(float)) {
  const REAL stalled_n8_xx[3] = {
      (REAL)2.0114874839782715,
      (REAL)-125.57471466064453,
      (REAL)189.93052673339844,
  };
  if (!check_roundtrip(&params, stalled_n8_xx,
                       xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN8,
                       generalrfm_Cart_to_xx__GeneralRFM_spheroidal_fisheyeN8))
    return 6;
} // END IF: float precision N8 inversion
if (sizeof(REAL) == sizeof(double)) {
  const REAL pi = (REAL)3.14159265358979323846;
  const REAL dtheta = pi / (REAL)100.0;
  const REAL dphi = (REAL)2.0 * pi / (REAL)200.0;
  for (int i_ph = 0; i_ph < 200; i_ph++) {
    const REAL phi = -pi + ((REAL)i_ph + (REAL)0.5) * dphi;
    for (int i_th = 0; i_th < 100; i_th++) {
      const REAL theta = ((REAL)i_th + (REAL)0.5) * dtheta;
      const REAL xCart[3] = {(REAL)11.5 * sin(theta) * cos(phi),
                             (REAL)11.5 * sin(theta) * sin(phi),
                             (REAL)11.5 * cos(theta)};
      REAL xx[3], xCart_back[3];
      if (generalrfm_Cart_to_xx__GeneralRFM_spheroidal_fisheyeN8(
              &params, xCart, xx) != 0)
        return 7;
      xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN8(&params, xx, xCart_back);
      if (!vectors_match(xCart, xCart_back))
        return 8;
    } // END LOOP: shell theta points
  } // END LOOP: shell phi points
} // END IF: double precision shell
return 0;
"""
    )
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        prefunc=VECTOR_CHECK,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )

    with TemporaryDirectory(prefix="nrpy-spheroidal-roundtrip-") as directory:
        project_dir = Path(directory)
        output_BHaH_defines_h(str(project_dir), enable_rfm_precompute=False)
        output_CFunctions_function_prototypes_and_construct_Makefile(
            project_dir=str(project_dir),
            project_name=PROJECT_NAME,
            exec_or_library_name=PROJECT_NAME,
            CC="gcc",
            use_openmp=False,
        )
        _run(["make", "-j2"], project_dir, project_dir / "build.log", 120)
        _run(
            [str(project_dir / PROJECT_NAME)], project_dir, project_dir / "run.log", 30
        )
    print(
        f"PASS [generated-product numerical]: {precision} N1 and N8 coordinate round trips"
    )


if __name__ == "__main__":
    if len(sys.argv) == 1:
        if doctest.testmod().failed:
            sys.exit(1)
    else:
        try:
            main()
        except (RuntimeError, subprocess.TimeoutExpired) as error:
            print(error, file=sys.stderr)
            sys.exit(1)
