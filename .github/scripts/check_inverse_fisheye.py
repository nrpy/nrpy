"""Compile generated inverse spheroidal fisheye C and check numerical contracts."""

import ctypes as ct
import math
import os
import subprocess
import sys
import tempfile
from functools import partial
from pathlib import Path
from typing import Any

import mpmath as mp  # type: ignore

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.BHaH.CurviBoundaryConditions.apply_bcs_outerradiation_and_inner import (
    register_CFunction_apply_bcs_outerradiation_and_inner,
)
from nrpy.infrastructures.BHaH.generalrfm_cart_to_xx import (
    register_CFunction_generalrfm_Cart_to_xx,
)
from nrpy.infrastructures.BHaH.inverse_fisheye import (
    register_CFunctions,
    register_precompute,
)
from nrpy.infrastructures.BHaH.numerical_grids_and_timestep import (
    register_CFunction_ds_min_single_pt,
)
from nrpy.infrastructures.BHaH.xx_tofrom_Cart import register_CFunction_xx_to_Cart


def reference_scale(params: ct.Structure, n: int, *position: mp.mpf) -> mp.mpf:
    """
    Evaluate equations (26)--(28) independently at high precision.

    :param params: Physical transition parameters and normalization.
    :param n: Transition count.
    :param position: Three physical Cartesian coordinates.
    :return: Computational radius divided by physical radius.
    """
    result = mp.mpf(1) / getattr(params, f"inverse_fisheye_a{n}")
    previous = mp.mpf(1)
    for j in range(1, n + 1):
        a = mp.mpf(getattr(params, f"inverse_fisheye_a{j}"))
        R = mp.mpf(getattr(params, f"inverse_fisheye_R{j}"))
        Rz = mp.mpf(getattr(params, f"inverse_fisheye_Rz{j}"))
        w = mp.mpf(getattr(params, f"inverse_fisheye_w{j}"))
        rho = mp.sqrt(position[0] ** 2 + position[1] ** 2 + (R / Rz * position[2]) ** 2)
        ratio = (
            1
            if not rho
            else w
            * mp.log(mp.cosh(2 * (rho + R) / w) / mp.cosh(2 * (rho - R) / w))
            / (4 * mp.tanh(2 * R / w) * rho)
        )
        result += (1 / previous - 1 / a) * ratio
        previous = a
    return result * getattr(params, "inverse_fisheye_c")


def main() -> None:
    """Check generated float/double maps, normalization, Jacobians, and metric jets."""
    par.set_parval_from_str("parallelization", "openmp")
    mp.mp.dps = 70
    for real_type, scalar, tolerance in [
        ("double", ct.c_double, 3e-12),
        ("float", ct.c_float, 3e-5),
    ]:
        cfc.CFunction_dict.clear()
        par.set_parval_from_str("fp_type", real_type)
        for n in (1, 2):
            register_CFunctions(n)
            coord = f"GeneralRFM_inverse_fisheyeN{n}"
            register_CFunction_xx_to_Cart(coord)
            register_CFunction_generalrfm_Cart_to_xx(coord)
            register_precompute(coord, n, enable_xx_to_Cart_precompute=n == 2)
            register_CFunction_ds_min_single_pt(coord)
        fields = [
            f"inverse_fisheye_{name}{j}"
            for name in ("a", "R", "Rz", "w")
            for j in (1, 2)
        ]
        fields += [
            "inverse_fisheye_c",
            "grid_physical_size",
            "xmin",
            "xmax",
            "ymin",
            "ymax",
            "zmin",
            "zmax",
            "Cart_originx",
            "Cart_originy",
            "Cart_originz",
        ]

        fields += ["dxx0", "dxx1", "dxx2"]
        int_fields = [f"Nxx_plus_2NGHOSTS{i}" for i in range(3)]
        gf_names = sorted(gri.glb_gridfcs_dict)

        class Params(ct.Structure):
            """Exact minimal params layout for the generated functions under test."""

            _fields_ = [(name, scalar) for name in fields] + [
                (name, ct.c_int) for name in int_fields
            ]

        # ctypes arrays return Python numbers; typeshed models these as C wrappers.
        vector: Any = scalar * 3
        matrix = vector * 3
        rank3 = matrix * 3
        rank4 = rank3 * 3
        source = f"""
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#define REAL {real_type}
#define DOUBLE double
#define NRPYMIN fmin
#define NRPYMAX fmax
typedef struct {{ {''.join(f'REAL {name};' for name in fields)} {''.join(f'int {name};' for name in int_fields)} }} params_struct;
typedef struct {{ int unused; }} commondata_struct;
#define IDX4(g,i,j,k) ((g)*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2+(i)+Nxx_plus_2NGHOSTS0*((j)+Nxx_plus_2NGHOSTS1*(k)))
"""
        source += (
            "\n".join(
                f"#define {name.upper()}GF {i}" for i, name in enumerate(gf_names)
            )
            + "\n"
        )
        source += "\n".join(
            f.cfunc_type + " " + f.name + "(" + f.params + ");"
            for f in cfc.CFunction_dict.values()
        )
        source += "\n".join(
            "\n".join(
                line
                for line in f.full_function.splitlines()
                if not line.startswith("#include")
            )
            for f in cfc.CFunction_dict.values()
        )
        with tempfile.TemporaryDirectory(prefix="nrpy-inverse-fisheye-") as tmp:
            c_path = Path(tmp) / "check.c"
            so_path = Path(tmp) / "check.so"
            c_path.write_text(source, encoding="utf-8")
            subprocess.run(
                [
                    os.environ.get("CC", "cc"),
                    "-std=c99",
                    "-O2",
                    "-shared",
                    "-fPIC",
                    str(c_path),
                    "-lm",
                    "-o",
                    str(so_path),
                ],
                check=True,
                timeout=120,
            )
            library = ct.CDLL(str(so_path))
            for n in (1, 2):
                for polar_ratio in (0.5, 1.0, 2.0):
                    params = Params(
                        inverse_fisheye_a1=2,
                        inverse_fisheye_a2=5,
                        inverse_fisheye_R1=2,
                        inverse_fisheye_R2=5,
                        inverse_fisheye_Rz1=2 * polar_ratio,
                        inverse_fisheye_Rz2=5 / polar_ratio,
                        inverse_fisheye_w1=0.5,
                        inverse_fisheye_w2=1.3,
                        grid_physical_size=12,
                    )
                    normalize = getattr(library, f"inverse_fisheye_normalize__N{n}")
                    scale = getattr(library, f"inverse_fisheye_scale__N{n}")
                    solve = getattr(library, f"inverse_fisheye_solve__N{n}")
                    geometry = getattr(library, f"inverse_fisheye_geometry__N{n}")
                    normalize(ct.byref(params))

                    reference = partial(reference_scale, params, n)

                    assert abs(reference(mp.mpf(12), 0, 0) - 1) < tolerance
                    assert (
                        abs(12 * reference(0, 0, mp.mpf(12)) - params.zmax)
                        < 20 * tolerance
                    )
                    if polar_ratio == 1:
                        assert abs(params.zmax - 12) < 20 * tolerance
                    points = [
                        (0, 0, 0),
                        (1e-12, 2e-12, -3e-12),
                        (0.003, 0.002, 0.001),
                        (0.02, 0, 0),
                        (1.3, -0.4, 0.6),
                        (3, 2, -1),
                        (1000, 200, -300),
                    ]
                    for point in points:
                        x = vector(*point)
                        physical = tuple(mp.mpf(float(v)) for v in x)
                        value, grad, hess, third = scalar(), vector(), matrix(), rank3()
                        scale(
                            ct.byref(params),
                            x,
                            3,
                            ct.byref(value),
                            grad,
                            hess,
                            third,
                        )
                        # The inexpensive paths must match full derivatives, including
                        # origin/Taylor cases. Unrequested arrays may be null.
                        fast_value, fast_grad = scalar(), vector()
                        scale(
                            ct.byref(params),
                            x,
                            0,
                            ct.byref(fast_value),
                            None,
                            None,
                            None,
                        )
                        assert fast_value.value == value.value
                        scale(
                            ct.byref(params),
                            x,
                            1,
                            ct.byref(fast_value),
                            fast_grad,
                            None,
                            None,
                        )
                        assert fast_value.value == value.value
                        assert list(fast_grad) == list(grad)
                        expected = reference(*physical)
                        assert abs(value.value - expected) < tolerance * max(
                            1, abs(expected)
                        )
                        q, recovered = vector(*(value.value * v for v in x)), vector()
                        assert solve(ct.byref(params), q, recovered) == 0
                        assert max(
                            abs(recovered[i] - x[i]) for i in range(3)
                        ) < 5 * tolerance * max(1, *map(abs, x))
                        J, K, g, dg, ddg = (
                            matrix(),
                            matrix(),
                            matrix(),
                            rank3(),
                            rank4(),
                        )
                        geometry(ct.byref(params), x, True, J, K, g, dg, ddg)
                        fast_J, fast_K = matrix(), matrix()
                        geometry(
                            ct.byref(params), x, False, fast_J, fast_K, None, None, None
                        )
                        assert bytes(fast_J) == bytes(J)
                        assert bytes(fast_K) == bytes(K)
                        for i in range(3):
                            for j in range(3):
                                assert (
                                    abs(
                                        sum(J[i][a] * K[a][j] for a in range(3))
                                        - (i == j)
                                    )
                                    < 20 * tolerance
                                )
                                assert (
                                    abs(
                                        g[i][j]
                                        - sum(K[a][i] * K[a][j] for a in range(3))
                                    )
                                    < 20 * tolerance
                                )
                        if point == (1.3, -0.4, 0.6) and real_type == "double":
                            for i in range(3):
                                order = [0, 0, 0]
                                order[i] += 1
                                assert (
                                    abs(
                                        grad[i]
                                        - mp.diff(reference, physical, tuple(order))
                                    )
                                    < 1e-11
                                )
                                for j in range(3):
                                    order[j] += 1
                                    assert (
                                        abs(
                                            hess[i][j]
                                            - mp.diff(reference, physical, tuple(order))
                                        )
                                        < 1e-10
                                    )
                                    for k in range(3):
                                        order[k] += 1
                                        assert (
                                            abs(
                                                third[i][j][k]
                                                - mp.diff(
                                                    reference,
                                                    physical,
                                                    tuple(order),
                                                )
                                            )
                                            < 1e-9
                                        )
                                        order[k] -= 1
                                    order[j] -= 1
                            # Differentiate g and dg along computational coordinates,
                            # exercising the inverse-function chain rule independently.
                            for k in range(3):
                                metrics, gradients = [], []
                                step = 1e-4
                                for sign in (-1, 1):
                                    shifted = vector(*q)
                                    shifted[k] += sign * step
                                    assert (
                                        solve(ct.byref(params), shifted, recovered) == 0
                                    )
                                    gm, dgm, ddgm = matrix(), rank3(), rank4()
                                    geometry(
                                        ct.byref(params),
                                        recovered,
                                        True,
                                        J,
                                        K,
                                        gm,
                                        dgm,
                                        ddgm,
                                    )
                                    metrics.append(gm)
                                    gradients.append(dgm)
                                for i in range(3):
                                    for j in range(3):
                                        estimate = (
                                            metrics[1][i][j] - metrics[0][i][j]
                                        ) / (2 * step)
                                        assert abs(estimate - dg[i][j][k]) < 2e-6
                                        for l in range(3):
                                            estimate = (
                                                gradients[1][i][j][l]
                                                - gradients[0][i][j][l]
                                            ) / (2 * step)
                                            assert (
                                                abs(estimate - ddg[i][j][l][k]) < 2e-5
                                            )
                        assert all(
                            math.isfinite(ddg[i][j][k][l])
                            for i in range(3)
                            for j in range(3)
                            for k in range(3)
                            for l in range(3)
                        )
                    # Exercise actual AUXEVOL storage and spacing consumers at one gridpoint.
                    for i in range(3):
                        setattr(params, f"Nxx_plus_2NGHOSTS{i}", 1)
                        setattr(params, f"dxx{i}", 0.1 * (i + 1))
                    q = vector(1, 2, 3)
                    recovered = vector()
                    assert solve(ct.byref(params), q, recovered) == 0
                    J, K, g, dg, ddg = matrix(), matrix(), matrix(), rank3(), rank4()
                    geometry(ct.byref(params), recovered, True, J, K, g, dg, ddg)
                    axes = [vector(q[i], 0, 0) for i in range(3)]
                    xx = (ct.POINTER(scalar) * 3)(*axes)
                    storage = (scalar * len(gf_names))()
                    precompute = getattr(
                        library,
                        f"generalrfm_precompute__rfm__GeneralRFM_inverse_fisheyeN{n}",
                    )
                    for axis, offset in zip("xyz", (1, -2, 3)):
                        setattr(params, f"Cart_origin{axis}", offset)
                    precompute(None, ct.byref(params), xx, storage)
                    if n == 2:
                        for i, offset in enumerate((1, -2, 3)):
                            assert (
                                abs(
                                    storage[gf_names.index(f"xCart{i}")]
                                    - recovered[i]
                                    - offset
                                )
                                < 20 * tolerance
                            )
                    for i in range(3):
                        for j in range(i, 3):
                            assert (
                                abs(storage[gf_names.index(f"ghatDD{i}{j}")] - g[i][j])
                                < 20 * tolerance
                            )
                            for k in range(3):
                                assert (
                                    abs(
                                        storage[gf_names.index(f"ghatDDdD{i}{j}{k}")]
                                        - dg[i][j][k]
                                    )
                                    < 20 * tolerance
                                )
                                for l in range(k, 3):
                                    assert (
                                        abs(
                                            storage[
                                                gf_names.index(f"ghatDDdDD{i}{j}{k}{l}")
                                            ]
                                            - ddg[i][j][k][l]
                                        )
                                        < 20 * tolerance
                                    )
                    spacing = getattr(
                        library,
                        f"ds_min_single_pt__rfm__GeneralRFM_inverse_fisheyeN{n}",
                    )
                    result_spacing = scalar()
                    spacing(
                        ct.byref(params),
                        scalar(q[0]),
                        scalar(q[1]),
                        scalar(q[2]),
                        ct.byref(result_spacing),
                    )
                    expected_spacing = min(
                        math.sqrt(g[i][i]) * getattr(params, f"dxx{i}")
                        for i in range(3)
                    )
                    assert abs(result_spacing.value - expected_spacing) < 20 * tolerance
                    # Public coordinate wrapper includes translation; the direct helper is patch-local.
                    for axis, offset in zip("xyz", [1, -2, 3]):
                        setattr(params, "Cart_origin" + axis, offset)
                    coord = f"GeneralRFM_inverse_fisheyeN{n}"
                    forward = getattr(library, f"xx_to_Cart__rfm__{coord}")
                    inverse = getattr(library, f"generalrfm_Cart_to_xx__{coord}")
                    q, cart, result = vector(1, 2, 3), vector(), vector()
                    forward(ct.byref(params), q, cart)
                    local = vector(cart[0] - 1, cart[1] + 2, cart[2] - 3)
                    assert inverse(ct.byref(params), local, result) == 0
                    assert max(abs(result[i] - q[i]) for i in range(3)) < 20 * tolerance
                    assert (
                        inverse(ct.byref(params), vector(float("nan"), 0, 0), result)
                        != 0
                    )
                    # Repeated components share a cache, but every map parameter
                    # must invalidate it even when the params pointer is unchanged.
                    component = getattr(library, f"inverse_fisheye_component__N{n}")
                    component.restype = scalar
                    component.argtypes = [
                        ct.POINTER(Params),
                        scalar,
                        scalar,
                        scalar,
                        ct.c_int,
                        ct.c_int,
                        ct.c_int,
                    ]
                    map_fields = ["inverse_fisheye_c"] + [
                        f"inverse_fisheye_{key}{j}"
                        for j in range(1, n + 1)
                        for key in ("a", "R", "Rz", "w")
                    ]
                    for field in map_fields:
                        original = getattr(params, field)
                        for factor in (1.01, 1.0):
                            setattr(params, field, original * factor)
                            assert solve(ct.byref(params), q, result) == 0
                            geometry(ct.byref(params), result, False, J, K, g, dg, ddg)
                            for i in range(3):
                                actual = component(ct.byref(params), *q, 0, i, 0)
                                assert abs(actual - result[i]) < 20 * tolerance
                                for j in range(3):
                                    for kind, expected in ((1, K), (2, J)):
                                        actual = component(
                                            ct.byref(params), *q, kind, i, j
                                        )
                                        assert (
                                            abs(actual - expected[i][j])
                                            < 20 * tolerance
                                        )
                    # Identity stretches normalize exactly to Cartesian coordinates.
                    for j in (1, 2):
                        setattr(params, f"inverse_fisheye_a{j}", 1)
                    normalize(ct.byref(params))
                    assert params.inverse_fisheye_c == 1 and params.zmax == 12
                    print(f"PASS {real_type} N={n} polar ratio={polar_ratio}")


def check_wave_project(enable_xx_to_Cart_precompute: bool = False) -> None:
    """
    Build and evolve a small oblate Gaussian-wave project with runtime inputs.

    :param enable_xx_to_Cart_precompute: Exercise radiation BCs using cached
        Cartesian positions in AUXEVOL, including the outer ghost layers.
    """
    root = Path(__file__).resolve().parents[2]
    source = (root / "nrpy/examples/wave_equation_curvilinear.py").read_text(
        encoding="utf-8"
    )
    source = source.replace(
        'CoordSystem = "SinhCylindrical"',
        'CoordSystem = "GeneralRFM_inverse_fisheyeN1"',
    )
    source = source.replace(
        '    "GeneralRFM_fisheyeN1": [64, 64, 64],',
        '    "GeneralRFM_inverse_fisheyeN1": [8, 8, 8],',
    )
    source = source.replace(
        "enable_parallel_codegen = True", "enable_parallel_codegen = False"
    )
    if enable_xx_to_Cart_precompute:
        source = source.replace(
            "BHaH.generalrfm_precompute.register_CFunction_generalrfm_precompute(\n",
            "BHaH.generalrfm_precompute.register_CFunction_generalrfm_precompute(\n"
            "                enable_xx_to_Cart_precompute=True,\n",
        )
        source = source.replace(
            "radiation_BC_fd_order=radiation_BC_fd_order,",
            "radiation_BC_fd_order=radiation_BC_fd_order,\n"
            "    enable_xx_to_Cart_precompute=True,",
        )
        source = source.replace(
            "\n                                     RK_INPUT_GFS, RK_OUTPUT_GFS);",
            "\n                                     RK_INPUT_GFS, RK_OUTPUT_GFS, auxevol_gfs);",
        )
    env = dict(os.environ, PYTHONPATH=str(root), OMP_NUM_THREADS="2")
    with tempfile.TemporaryDirectory(prefix="nrpy-inverse-wave-") as temporary:
        directory = Path(temporary)
        generator = directory / "generate.py"
        generator.write_text(source, encoding="utf-8")
        project = directory / "project/wave_equation_curvilinear"
        commands = [
            ([sys.executable, str(generator), "--disable_intrinsics"], directory),
            (["make", "-j2"], project),
        ]
        for command, cwd in commands:
            with (directory / "build.log").open("w", encoding="utf-8") as log:
                subprocess.run(
                    command,
                    cwd=cwd,
                    env=env,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                    timeout=180,
                )
        if enable_xx_to_Cart_precompute:
            coord = "GeneralRFM_inverse_fisheyeN1"
            register_CFunction_apply_bcs_outerradiation_and_inner(coord)
            name = f"apply_bcs_outerradiation_and_inner__rfm__{coord}"
            reference = directory / "reference_boundary.c"
            reference.write_text(
                cfc.CFunction_dict[name].full_function.replace(
                    name, "reference_boundary"
                ),
                encoding="utf-8",
            )
            driver = directory / "check_boundary.c"
            driver.write_text(
                r"""
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include <assert.h>
void reference_boundary(const commondata_struct *, const params_struct *, const bc_struct *,
                        REAL *restrict [3], const REAL *, const REAL *, REAL *, REAL *);
int main(int argc, const char *argv[]) {
  (void)argc; (void)argv;
  commondata_struct commondata;
  commondata_struct_set_to_default(&commondata);
  commondata.NUMGRIDS=1;
  commondata.inverse_fisheye_phys_a1=2;
  commondata.inverse_fisheye_phys_R1=2;
  commondata.inverse_fisheye_phys_Rz1=1;
  commondata.inverse_fisheye_phys_w1=0.5;
  griddata_struct *griddata=calloc(MAXNUMGRIDS,sizeof(griddata_struct));
  params_struct_set_to_default(&commondata,griddata);
  params_struct *params=&griddata[0].params;
  params->CoordSystem_hash=GENERALRFM_INVERSE_FISHEYEN1;
  params->grid_physical_size=4;
  params->Cart_originx=1; params->Cart_originy=-2; params->Cart_originz=3;
  const int Nx[3]={8,10,12};
  REAL *xx[3];
  commondata.inverse_fisheye_phys_Rz1=0;
  numerical_grid_params_Nxx_dxx_xx(&commondata,params,xx,Nx,true);
  assert(params->inverse_fisheye_Rz1 == params->inverse_fisheye_R1);
  assert(fabs(params->xmax-params->zmax)<1e-12);
  for(int i=0;i<3;i++) free(xx[i]);
  commondata.inverse_fisheye_phys_Rz1=1;
  numerical_grid_params_Nxx_dxx_xx(&commondata,params,xx,Nx,true);
  const int Nxx_plus_2NGHOSTS0=params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1=params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2=params->Nxx_plus_2NGHOSTS2;
  const int count=Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;
  REAL *auxevol_gfs=calloc(NUM_AUXEVOL_GFS*count,sizeof(REAL));
  generalrfm_precompute(&commondata,params,(const REAL *restrict *)xx,auxevol_gfs);
  for(int i2=0;i2<Nxx_plus_2NGHOSTS2;i2++)
    for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++)
      for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
        const REAL q[3]={xx[0][i0],xx[1][i1],xx[2][i2]};
        REAL x[3]; xx_to_Cart(params,q,x);
        assert(fabs(auxevol_gfs[IDX4(XCART0GF,i0,i1,i2)]-x[0])<1e-11);
        assert(fabs(auxevol_gfs[IDX4(XCART1GF,i0,i1,i2)]-x[1])<1e-11);
        assert(fabs(auxevol_gfs[IDX4(XCART2GF,i0,i1,i2)]-x[2])<1e-11);
      }
  bc_struct bc={0}; bcstruct_set_up(&commondata,params,xx,&bc);
  REAL *gfs=malloc(NUM_EVOL_GFS*count*sizeof(REAL));
  REAL *expected=malloc(NUM_EVOL_GFS*count*sizeof(REAL));
  REAL *actual=malloc(NUM_EVOL_GFS*count*sizeof(REAL));
  REAL speed[NUM_EVOL_GFS], infinity[NUM_EVOL_GFS];
  for(int gf=0;gf<NUM_EVOL_GFS;gf++) { speed[gf]=1+0.1*gf; infinity[gf]=0.02*gf; }
  for(int i=0;i<NUM_EVOL_GFS*count;i++) {
    gfs[i]=sin(0.03*i)+0.2*cos(0.011*i);
    expected[i]=actual[i]=0.1+0.003*(i%97);
  }
  reference_boundary(&commondata,params,&bc,xx,speed,infinity,gfs,expected);
  apply_bcs_outerradiation_and_inner(&commondata,params,&bc,xx,speed,infinity,gfs,actual,auxevol_gfs);
  REAL error=0;
  for(int i=0;i<NUM_EVOL_GFS*count;i++) {
    assert(isfinite(expected[i]) && isfinite(actual[i]));
    error=fmax(error,fabs(expected[i]-actual[i]));
  }
  assert(error<1e-10);
  printf("PASS cached coordinates including ghosts and translated origin; boundary RHS max error %.3e\n",error);
  for(int j=0;j<3;j++) free(xx[j]);
  for(int j=0;j<3*NGHOSTS;j++) free(bc.pure_outer_bc_array[j]);
  free(bc.inner_bc_array); free(auxevol_gfs); free(griddata);
  free(gfs); free(expected); free(actual);
  return 0;
}
""",
                encoding="utf-8",
            )
            executable = directory / "check_boundary"
            objects = sorted(
                path for path in project.rglob("*.o") if path.name != "main.o"
            )
            with (directory / "boundary_build.log").open("w", encoding="utf-8") as log:
                subprocess.run(
                    [
                        os.environ.get("CC", "cc"),
                        "-O2",
                        f"-I{project}",
                        str(driver),
                        str(reference),
                        *(str(path) for path in objects),
                        "-lm",
                        "-fopenmp",
                        "-o",
                        str(executable),
                    ],
                    env=env,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                    timeout=120,
                )
            subprocess.run([str(executable)], env=env, check=True, timeout=60)
        parfile = project / "wave_equation_curvilinear.par"
        text = parfile.read_text(encoding="utf-8")
        # Public shared physical inputs must be present and parsed by the generated application.
        for key in ("a", "R", "Rz", "w"):
            assert f"inverse_fisheye_phys_{key}1 =" in text
        text = text.replace("t_final = 8.0", "t_final = 0.2")
        text = text.replace("CFL_FACTOR = 0.5", "CFL_FACTOR = 0.01")
        text = text.replace(
            "diagnostics_output_every = 0.5", "diagnostics_output_every = 0.05"
        )
        text = text.replace(
            "inverse_fisheye_phys_a1 = 1.0", "inverse_fisheye_phys_a1 = 2.0"
        )
        text = text.replace(
            "inverse_fisheye_phys_R1 = 1.0", "inverse_fisheye_phys_R1 = 2.0"
        )
        text = text.replace(
            "inverse_fisheye_phys_Rz1 = 0.0", "inverse_fisheye_phys_Rz1 = 1.0"
        )
        parfile.write_text(text, encoding="utf-8")
        with (directory / "run.log").open("w", encoding="utf-8") as log:
            subprocess.run(
                [str(project / "wave_equation_curvilinear"), str(parfile)],
                cwd=project,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                timeout=120,
            )
        diagnostic = (
            project / "out0d-grid00-GeneralRFM_inverse_fisheyeN1-conv_factor-1.00.txt"
        )
        rows = [
            [float(value) for value in line.split()]
            for line in diagnostic.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")
        ]
        assert rows[-1][0] >= 0.14
        assert all(math.isfinite(value) for row in rows for value in row)
        assert max(abs(row[1]) for row in rows) < 0.001
        assert rows[-1][3] != rows[0][3]
        print(
            f"PASS oblate wave project (Cartesian cache={enable_xx_to_Cart_precompute}): "
            "t >= 0.14, sampled relative uu error < 0.001"
        )


if __name__ == "__main__":
    if "--wave-project" in sys.argv:
        check_wave_project()
        check_wave_project(enable_xx_to_Cart_precompute=True)
    else:
        main()
