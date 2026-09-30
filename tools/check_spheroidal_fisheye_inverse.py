#!/usr/bin/env python3
"""
Compile and run a small generated-project inverse check for spheroidal fisheye.

This avoids launching the full Charm++/MPI evolution. It compiles only the
project defaults, fisheye parameter converters, Cart->xx inverse wrapper, and
xx->Cart forward map, then round-trips known Cartesian failure fixtures.
"""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

DEFAULT_PROJECT_DIR = Path(
    "project/superB_blackhole_spectroscopy_8Mseparation_spheroidal_fisheye"
)
DEFAULT_POINTS = (
    (3.562131984779170, -2.151107620260850e-02, 3.973390482844859e-01),
    (3.999591595248627, -2.006363886584629e-05, 8.323295468376436e-03),
)


def _format_point(point: tuple[float, float, float]) -> str:
    return "{" + ", ".join(f"{value:.17e}" for value in point) + "}"


def _parse_point(raw: str) -> tuple[float, float, float]:
    values = [float(item) for item in raw.replace(",", " ").split()]
    if len(values) != 3:
        raise argparse.ArgumentTypeError(
            f"point must contain exactly three numbers; got {len(values)}"
        )
    return values[0], values[1], values[2]


def _harness_source(
    points: list[tuple[float, float, float]], *, sweep_stride: int | None
) -> str:
    point_rows = ",\n      ".join(_format_point(point) for point in points)
    sweep_body = ""
    if sweep_stride is not None:
        sweep_body = f"""
  {{
    const int stride = {sweep_stride};
    int checked = 0;
    const int i0_values[] = {{
        0, 1, 2, NGHOSTS, griddata[0].params.Nxx_plus_2NGHOSTS0 / 2,
        griddata[0].params.Nxx_plus_2NGHOSTS0 - NGHOSTS - 1,
        griddata[0].params.Nxx_plus_2NGHOSTS0 - 3,
        griddata[0].params.Nxx_plus_2NGHOSTS0 - 2,
        griddata[0].params.Nxx_plus_2NGHOSTS0 - 1,
    }};
    const int i1_values[] = {{
        0, 1, 2, NGHOSTS, griddata[0].params.Nxx_plus_2NGHOSTS1 / 2,
        griddata[0].params.Nxx_plus_2NGHOSTS1 - NGHOSTS - 1,
        griddata[0].params.Nxx_plus_2NGHOSTS1 - 3,
        griddata[0].params.Nxx_plus_2NGHOSTS1 - 2,
        griddata[0].params.Nxx_plus_2NGHOSTS1 - 1,
    }};
    const int i2_values[] = {{
        0, 1, 2, NGHOSTS, griddata[0].params.Nxx_plus_2NGHOSTS2 / 2,
        griddata[0].params.Nxx_plus_2NGHOSTS2 - NGHOSTS - 1,
        griddata[0].params.Nxx_plus_2NGHOSTS2 - 3,
        griddata[0].params.Nxx_plus_2NGHOSTS2 - 2,
        griddata[0].params.Nxx_plus_2NGHOSTS2 - 1,
    }};

    for (int i0 = 0; i0 < griddata[0].params.Nxx_plus_2NGHOSTS0; i0 += stride) {{
      for (int i1 = 0; i1 < griddata[0].params.Nxx_plus_2NGHOSTS1; i1 += stride) {{
        for (int i2 = 0; i2 < griddata[0].params.Nxx_plus_2NGHOSTS2; i2 += stride) {{
          REAL raw[3] = {{xx_arrays[0][i0], xx_arrays[1][i1], xx_arrays[2][i2]}};
          REAL cart[3];
          REAL raw_back[3];
          xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN8(
              &griddata[0].params, raw, cart);
          Cart_to_xx_and_nearest_i0i1i2_assume_valid__rfm__GeneralRFM_spheroidal_fisheyeN8(
              &griddata[0].params, cart, raw_back, NULL);
          for (int d = 0; d < 3; d++) {{
            const REAL scale = fmax((REAL)1.0, fabs(raw[d]));
            if (fabs(raw[d] - raw_back[d]) > (REAL)1.0e-8 * scale) {{
              fprintf(stderr,
                      "sweep mismatch i=(%d,%d,%d) component %d raw=%.17g back=%.17g cart=(%.17g %.17g %.17g)\\n",
                      i0, i1, i2, d, (double)raw[d], (double)raw_back[d],
                      (double)cart[0], (double)cart[1], (double)cart[2]);
              return 4;
            }}
          }}
          checked++;
        }}
      }}
    }}

    for (int a = 0; a < (int)(sizeof(i0_values) / sizeof(i0_values[0])); a++) {{
      for (int b = 0; b < (int)(sizeof(i1_values) / sizeof(i1_values[0])); b++) {{
        for (int c = 0; c < (int)(sizeof(i2_values) / sizeof(i2_values[0])); c++) {{
          const int i0 = i0_values[a];
          const int i1 = i1_values[b];
          const int i2 = i2_values[c];
          REAL raw[3] = {{xx_arrays[0][i0], xx_arrays[1][i1], xx_arrays[2][i2]}};
          REAL cart[3];
          REAL raw_back[3];
          xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN8(
              &griddata[0].params, raw, cart);
          Cart_to_xx_and_nearest_i0i1i2_assume_valid__rfm__GeneralRFM_spheroidal_fisheyeN8(
              &griddata[0].params, cart, raw_back, NULL);
          for (int d = 0; d < 3; d++) {{
            const REAL scale = fmax((REAL)1.0, fabs(raw[d]));
            if (fabs(raw[d] - raw_back[d]) > (REAL)1.0e-8 * scale) {{
              fprintf(stderr,
                      "boundary mismatch i=(%d,%d,%d) component %d raw=%.17g back=%.17g cart=(%.17g %.17g %.17g)\\n",
                      i0, i1, i2, d, (double)raw[d], (double)raw_back[d],
                      (double)cart[0], (double)cart[1], (double)cart[2]);
              return 5;
            }}
          }}
          checked++;
        }}
      }}
    }}
    fprintf(stderr, "sweep checked %d generated grid points with stride %d plus boundary set\\n",
            checked, stride);
  }}
"""
    return f"""\
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

extern "C" void CmiAbort(const char *msg, ...) {{
  (void)msg;
  __builtin_trap();
}}

int main() {{
  commondata_struct commondata;
  griddata_struct griddata[MAXNUMGRIDS];
  REAL *xx_arrays[3];
  const int default_Nx[3] = {{-1, -1, -1}};
  commondata_struct_set_to_default(&commondata);
  params_struct_set_to_default(&commondata, griddata);

  for (int grid = 0; grid < commondata.NUMGRIDS; grid++) {{
    griddata[grid].params.fisheye_a0 = commondata.fisheye_phys_a0;
    griddata[grid].params.fisheye_a1 = commondata.fisheye_phys_a1;
    griddata[grid].params.fisheye_a2 = commondata.fisheye_phys_a2;
    griddata[grid].params.fisheye_a3 = commondata.fisheye_phys_a3;
    griddata[grid].params.fisheye_a4 = commondata.fisheye_phys_a4;
    griddata[grid].params.fisheye_a5 = commondata.fisheye_phys_a5;
    griddata[grid].params.fisheye_a6 = commondata.fisheye_phys_a6;
    griddata[grid].params.fisheye_a7 = commondata.fisheye_phys_a7;
    griddata[grid].params.fisheye_a8 = commondata.fisheye_phys_a8;
    if (spheroidal_fisheye_params_from_physical_N8(
            &commondata, &griddata[grid].params) != 0) {{
      fprintf(stderr, "spheroidal_fisheye_params_from_physical_N8 failed\\n");
      return 2;
    }}
  }}
  numerical_grid_params_Nxx_dxx_xx__rfm__GeneralRFM_spheroidal_fisheyeN8(
      &commondata, &griddata[0].params, xx_arrays, default_Nx, true);

  const REAL points[][3] = {{
      {point_rows},
  }};
  const int num_points = (int)(sizeof(points) / sizeof(points[0]));
  for (int n = 0; n < num_points; n++) {{
    REAL xx[3];
    Cart_to_xx_and_nearest_i0i1i2_assume_valid__rfm__GeneralRFM_spheroidal_fisheyeN8(
        &griddata[0].params, points[n], xx, NULL);
    REAL back[3];
    xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN8(
        &griddata[0].params, xx, back);
    fprintf(stderr,
            "point %d xx=(%.17g %.17g %.17g) back=(%.17g %.17g %.17g)\\n",
            n, (double)xx[0], (double)xx[1], (double)xx[2],
            (double)back[0], (double)back[1], (double)back[2]);
    for (int i = 0; i < 3; i++) {{
      const REAL scale = fmax((REAL)1.0, fabs(points[n][i]));
      if (fabs(points[n][i] - back[i]) > (REAL)1.0e-10 * scale) {{
        fprintf(stderr, "roundtrip mismatch point %d component %d\\n", n, i);
        return 3;
      }}
    }}
  }}
{sweep_body}
  return 0;
}}
"""


def main() -> int:
    """Compile and run the generated spheroidal-fisheye inverse-map checks."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-dir",
        type=Path,
        default=DEFAULT_PROJECT_DIR,
        help="generated project directory to test",
    )
    parser.add_argument(
        "--point",
        action="append",
        type=_parse_point,
        help="Cartesian point as 'x y z' or 'x,y,z'; may be repeated",
    )
    parser.add_argument(
        "--sweep-stride",
        type=int,
        default=None,
        help=(
            "also sweep generated grid cell centers using this stride; "
            "use 1 for an exhaustive full-grid scan"
        ),
    )
    args = parser.parse_args()

    project_dir = args.project_dir.resolve()
    points = args.point if args.point else list(DEFAULT_POINTS)
    if args.sweep_stride is not None and args.sweep_stride < 1:
        raise ValueError("--sweep-stride must be >= 1")

    sources = [
        "commondata_struct_set_to_default.cpp",
        "params_struct_set_to_default.cpp",
        "fisheye/fisheye_params_from_physical_N8.cpp",
        "fisheye/spheroidal_fisheye_params_from_physical_N8.cpp",
        "GeneralRFM_spheroidal_fisheyeN8/numerical_grid_params_Nxx_dxx_xx__rfm__GeneralRFM_spheroidal_fisheyeN8.cpp",
        "GeneralRFM_spheroidal_fisheyeN8/Cart_to_xx_and_nearest_i0i1i2_assume_valid__rfm__GeneralRFM_spheroidal_fisheyeN8.cpp",
        "GeneralRFM_spheroidal_fisheyeN8/xx_to_Cart__rfm__GeneralRFM_spheroidal_fisheyeN8.cpp",
    ]
    missing = [source for source in sources if not (project_dir / source).exists()]
    if missing:
        raise FileNotFoundError(
            f"{project_dir} is missing required generated sources: {missing}"
        )

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        harness = tmpdir / "check_inverse.cpp"
        executable = tmpdir / "check_inverse"
        harness.write_text(
            _harness_source(points, sweep_stride=args.sweep_stride),
            encoding="utf-8",
        )
        command = [
            "g++",
            "-std=gnu++17",
            "-O0",
            "-I.",
            str(harness),
            *sources,
            "-lm",
            "-o",
            str(executable),
        ]
        subprocess.run(command, cwd=project_dir, check=True)
        subprocess.run([str(executable)], cwd=project_dir, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
