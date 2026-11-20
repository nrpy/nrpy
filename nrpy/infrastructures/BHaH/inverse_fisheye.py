"""
Generate BHaH spherical or spheroidal inverse fisheye maps, normalization, and metric jets.

Geometry helpers use patch-local physical coordinates; coordinate wrappers own
rotations. The metric uses inverse-function derivatives through third order,
so its two derivatives are with respect to computational coordinates.

Scale evaluation requests derivative order 0 for coordinate maps, 1 for
Jacobians, and 3 for metric jets. Unrequested derivative arrays are untouched
and may be NULL. Jacobian-only calls do not evaluate the second or third
scale derivatives.

Author: OpenAI
"""

from typing import List, Set, Tuple, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
from nrpy.equations.generalrfm.inverse_fisheye import (
    GeneralRFMInverseFisheye,
)


def register_CFunctions(num_transitions: int) -> None:
    """
    Register provider-local functions needed by coordinate and grid wrappers.

    :param num_transitions: Number of physical spheroidal transitions.
    :raises ValueError: If CUDA generation is requested.
    """
    if par.parval_from_str("parallelization") == "cuda":
        raise ValueError("Inverse spheroidal fisheye does not support CUDA.")
    suffix = f"__N{num_transitions}"
    if "inverse_fisheye_scale" + suffix in cfc.CFunction_dict:
        return
    model = GeneralRFMInverseFisheye(num_transitions)
    # BHaH parameter files expose commondata only. These are physical inputs;
    # grid setup copies them verbatim into the per-grid provider parameters.
    for key, values in [
        ("a", [1.0] * num_transitions),
        ("R", [float(j + 1) for j in range(num_transitions)]),
        ("Rz", [0.0] * num_transitions),
        ("w", [0.5] * num_transitions),
    ]:
        for j, value in enumerate(values):
            par.register_CodeParameter(
                "REAL",
                __name__,
                f"inverse_fisheye_phys_{key}{j+1}",
                value,
                commondata=True,
            )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    calc_type = "DOUBLE" if par.parval_from_str("fp_type") == "float" else "REAL"
    radius = sp.Symbol("rho", real=True, nonnegative=True)
    coeffs = model.kernel_coefficients(radius, 0)
    R_symbol, w_symbol = model.R_list[0], model.w_list[0]
    generic_R, generic_w = sp.Symbol("R", positive=True), sp.Symbol("w", positive=True)
    replacements = {R_symbol: generic_R, w_symbol: generic_w}
    branch_bodies: List[str] = []
    branch_conditions: List[str] = []
    # Keep each kernel expression inside its Taylor, interior, or exterior C branch.
    for branch in range(3):
        condition = cast(Tuple[sp.Expr, sp.Expr], coeffs[0].args[branch])[1]
        branch_conditions.append(sp.ccode(condition.xreplace(replacements)))
        expressions = [
            cast(Tuple[sp.Expr, sp.Expr], coeff.args[branch])[0].xreplace(replacements)
            for coeff in coeffs
        ]
        options = {
            "fp_type": (
                "double" if calc_type == "DOUBLE" else par.parval_from_str("fp_type")
            ),
            "fp_type_alias": calc_type,
            "include_braces": False,
            "verbose": False,
        }
        value_code = ccg.c_codegen(expressions[0], "f[0]", **options)
        gradient_code = ccg.c_codegen(expressions[1], "f[1]", **options)
        metric_code = ccg.c_codegen(expressions[2:], ["f[2]", "f[3]"], **options)
        branch_bodies.append(
            f"{value_code}\n"
            f"if (derivative_order >= 1) {{\n{gradient_code}\n"
            f"if (derivative_order == 3) {{\n{metric_code}\n}}\n}}"
        )
    prefunc = f"""
static inline void inverse_fisheye_layer{suffix}(
    const {calc_type} R, const {calc_type} w, const {calc_type} rho,
    const int derivative_order, {calc_type} f[4]) {{
  if ({branch_conditions[0]}) {{
    {branch_bodies[0]}
  }} else if ({branch_conditions[1]}) {{
    {branch_bodies[1]}
  }} else {{
    {branch_bodies[2]}
  }}
}}
"""
    arrays = {
        key: ", ".join(
            f"params->inverse_fisheye_{key}{j}" for j in range(1, num_transitions + 1)
        )
        for key in ("a", "R", "Rz", "w")
    }
    body = f"""
const {calc_type} a[{num_transitions}] = {{{arrays['a']}}};
const {calc_type} R[{num_transitions}] = {{{arrays['R']}}};
const {calc_type} Rz[{num_transitions}] = {{{arrays['Rz']}}};
const {calc_type} w[{num_transitions}] = {{{arrays['w']}}};
{calc_type} lam = ({calc_type})params->inverse_fisheye_c / a[{num_transitions - 1}];
{calc_type} grad[3] = {{0}}, hess[3][3] = {{{{0}}}}, third[3][3][3] = {{{{{{0}}}}}};
{calc_type} previous_inverse_a = 1;
for (int layer=0; layer<{num_transitions}; layer++) {{
  const {calc_type} inverse_a = ({calc_type})1/a[layer];
  const {calc_type} delta = ({calc_type})params->inverse_fisheye_c*(previous_inverse_a-inverse_a);
  previous_inverse_a = inverse_a;
  if (delta == 0) continue;
  const {calc_type} axis_ratio = R[layer]/Rz[layer];
  const {calc_type} b[3] = {{1, 1, axis_ratio*axis_ratio}};
  const {calc_type} h[3] = {{x[0], x[1], b[2]*x[2]}};
  const {calc_type} rho = sqrt(({calc_type})x[0]*x[0]+({calc_type})x[1]*x[1]+b[2]*x[2]*x[2]);
  {calc_type} f[4];
  inverse_fisheye_layer{suffix}(R[layer], w[layer], rho, derivative_order, f);
  lam += delta*f[0];
  if (derivative_order >= 1) {{
    for (int i=0; i<3; i++) grad[i] += delta*2*f[1]*h[i];
    if (derivative_order == 3) {{
      for (int i=0; i<3; i++) for (int k=0; k<3; k++) {{
        hess[i][k] += delta*(4*f[2]*h[i]*h[k]+2*f[1]*b[i]*(i==k));
        for (int l=0; l<3; l++)
          third[i][k][l] += delta*(8*f[3]*h[i]*h[k]*h[l]
            +4*f[2]*((i==k)*b[i]*h[l]+(i==l)*b[i]*h[k]+(k==l)*b[k]*h[i]));
      }}
    }}
  }}
}}
"""
    body += """
*lambda = lam;
if (derivative_order >= 1) for (int i=0; i<3; i++) v[i] = grad[i];
if (derivative_order == 3)
  for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
    W[i][j] = hess[i][j];
    for (int k=0; k<3; k++) U[i][j][k] = third[i][j][k];
  }
"""
    cfc.register_CFunction(
        includes=includes,
        desc="Evaluate scale with derivative_order 0 (value), 1 (gradient), or 3 (full metric jet). Unrequested derivative arrays are untouched and may be NULL.",
        cfunc_type="void",
        name="inverse_fisheye_scale" + suffix,
        params="const params_struct *restrict params, const REAL x[3], const int derivative_order, REAL *lambda, REAL v[3], REAL W[3][3], REAL U[3][3][3]",
        include_CodeParameters_h=False,
        prefunc=prefunc,
        body=body,
    )

    body = f"""
const REAL qnorm = hypot(hypot(q[0], q[1]), q[2]);
if (!isfinite(qnorm) || !(params->inverse_fisheye_c > 0)) return 1;
if (qnorm == 0) {{ x[0]=x[1]=x[2]=0; return 0; }}
const REAL n[3] = {{q[0]/qnorm, q[1]/qnorm, q[2]/qnorm}};
REAL lo=0, hi=qnorm*params->inverse_fisheye_a{num_transitions}/params->inverse_fisheye_c;
if (!isfinite(hi) || hi <= 0) return 1;
// lambda >= c/aN gives an analytic upper bracket, for every ray.
const REAL tol = (sizeof(REAL)==sizeof(float)) ? 2e-6 : 2e-14;
for (int iteration=0; iteration<100; iteration++) {{
  const REAL r=lo+(hi-lo)/2;
  for (int i=0; i<3; i++) x[i]=r*n[i];
  REAL lambda;
  inverse_fisheye_scale{suffix}(params, x, 0, &lambda, NULL, NULL, NULL);
  const REAL residual=r*lambda-qnorm;
  if (!isfinite(residual)) return 1;
  if (fabs(residual) <= tol*qnorm) return 0;
  if (residual > 0) hi=r; else lo=r;
}}
return 1;
"""
    cfc.register_CFunction(
        includes=includes,
        desc="Solve the monotone scalar inverse along a physical ray.",
        cfunc_type="int",
        name="inverse_fisheye_solve" + suffix,
        params="const params_struct *restrict params, const REAL q[3], REAL x[3]",
        include_CodeParameters_h=False,
        body=body,
    )

    body = f"""
REAL lambda, v[3], W[3][3], U[3][3][3];
inverse_fisheye_scale{suffix}(params, x, derivatives ? 3 : 1, &lambda, v, W, U);
const REAL radial=lambda+x[0]*v[0]+x[1]*v[1]+x[2]*v[2];
for (int a=0; a<3; a++) for (int i=0; i<3; i++) {{
  J[a][i]=lambda*(a==i)+x[a]*v[i];
  K[a][i]=((a==i)-x[a]*v[i]/radial)/lambda;
}}
if (!derivatives) return;
REAL F2[3][3][3], F3[3][3][3][3];
for (int a=0; a<3; a++) for (int b=0; b<3; b++) for (int c=0; c<3; c++) {{
  F2[a][b][c]=(a==b)*v[c]+(a==c)*v[b]+x[a]*W[b][c];
  for (int d=0; d<3; d++)
    F3[a][b][c][d]=(a==b)*W[c][d]+(a==c)*W[b][d]+(a==d)*W[b][c]+x[a]*U[b][c][d];
}}
// H = d^2 x / dq^2, T = d^3 x / dq^3 by implicit differentiation.
REAL H[3][3][3]={{0}}, T[3][3][3][3]={{0}};
for (int a=0; a<3; a++) for (int i=0; i<3; i++) for (int j=0; j<3; j++)
  for (int b=0; b<3; b++) for (int c=0; c<3; c++) for (int d=0; d<3; d++)
    H[a][i][j]-=K[a][b]*F2[b][c][d]*K[c][i]*K[d][j];
for (int a=0; a<3; a++) for (int i=0; i<3; i++) for (int j=0; j<3; j++) for (int k=0; k<3; k++)
  for (int b=0; b<3; b++) for (int c=0; c<3; c++) for (int d=0; d<3; d++) {{
    T[a][i][j][k]-=K[a][b]*F2[b][c][d]*(H[c][i][j]*K[d][k]+H[c][i][k]*K[d][j]+H[c][j][k]*K[d][i]);
    for (int e=0; e<3; e++)
      T[a][i][j][k]-=K[a][b]*F3[b][c][d][e]*K[c][i]*K[d][j]*K[e][k];
  }}
for (int i=0; i<3; i++) for (int j=0; j<3; j++) {{
  g[i][j]=0;
  for (int a=0; a<3; a++) g[i][j]+=K[a][i]*K[a][j];
  for (int k=0; k<3; k++) {{
    dg[i][j][k]=0;
    for (int a=0; a<3; a++) dg[i][j][k]+=H[a][i][k]*K[a][j]+K[a][i]*H[a][j][k];
    for (int l=0; l<3; l++) {{
      ddg[i][j][k][l]=0;
      for (int a=0; a<3; a++) ddg[i][j][k][l]+=T[a][i][k][l]*K[a][j]
        +H[a][i][k]*H[a][j][l]+H[a][i][l]*H[a][j][k]+K[a][i]*T[a][j][k][l];
    }}
  }}
}}
"""
    cfc.register_CFunction(
        includes=includes,
        desc="Evaluate inverse-map Jacobians and computational metric derivatives at physical x.",
        cfunc_type="void",
        name="inverse_fisheye_geometry" + suffix,
        params="const params_struct *restrict params, const REAL x[3], const bool derivatives, REAL J[3][3], REAL K[3][3], REAL g[3][3], REAL dg[3][3][3], REAL ddg[3][3][3][3]",
        include_CodeParameters_h=False,
        body=body,
    )

    # Bounded thread-local cache: generic symbolic consumers request distinct
    # components at the same point. Key every map parameter, not a params pointer,
    # so switching chares or changing parameters cannot reuse stale geometry.
    cache_keys = ["q0", "q1", "q2", "params->inverse_fisheye_c"] + [
        f"params->inverse_fisheye_{key}{j}"
        for j in range(1, num_transitions + 1)
        for key in ("a", "R", "Rz", "w")
    ]
    key_count = len(cache_keys)
    body = f"""
#ifdef __cplusplus
#define INVERSE_FISHEYE_THREAD_LOCAL thread_local
#else
#define INVERSE_FISHEYE_THREAD_LOCAL _Thread_local
#endif
static INVERSE_FISHEYE_THREAD_LOCAL bool valid=false, jacobian_valid=false;
static INVERSE_FISHEYE_THREAD_LOCAL REAL saved_key[{key_count}], x[3], J[3][3], K[3][3];
#undef INVERSE_FISHEYE_THREAD_LOCAL
const REAL key[{key_count}]={{{", ".join(cache_keys)}}};
bool same=valid;
for (int k=0; k<{key_count}; k++) if (key[k]!=saved_key[k]) same=false;
if (!same) {{
  const REAL q[3]={{q0,q1,q2}};
  if (inverse_fisheye_solve{suffix}(params, q, x)) {{
    fprintf(stderr, "Inverse spheroidal fisheye solve failed.\\n"); exit(1);
  }}
  for (int k=0; k<{key_count}; k++) saved_key[k]=key[k];
  valid=true;
  jacobian_valid=false;
}}
if (kind==0) return x[i];
if (!jacobian_valid) {{
  REAL g[3][3], dg[3][3][3], ddg[3][3][3][3];
  inverse_fisheye_geometry{suffix}(params, x, false, J, K, g, dg, ddg);
  jacobian_valid=true;
}}
return (kind==1) ? K[i][j] : J[i][j];
"""
    cfc.register_CFunction(
        includes=includes,
        desc="Evaluate symbolic map/Jacobian components for generic BHaH consumers.",
        cfunc_type="REAL",
        name="inverse_fisheye_component" + suffix,
        params="const params_struct *restrict params, const REAL q0, const REAL q1, const REAL q2, const int kind, const int i, const int j",
        include_CodeParameters_h=False,
        body=body,
    )

    checks = """
if (!(isfinite(params->grid_physical_size) && params->grid_physical_size > 0)) {
  fprintf(stderr, "Inverse fisheye requires positive finite physical extent.\\n"); exit(1);
}
"""
    for layer in range(num_transitions):
        j = layer + 1
        prev = "1" if layer == 0 else f"params->inverse_fisheye_a{layer}"
        checks += f"""
if (!(isfinite(params->inverse_fisheye_a{j}) && params->inverse_fisheye_a{j} >= {prev}
   && isfinite(params->inverse_fisheye_R{j}) && params->inverse_fisheye_R{j} > 0
   && isfinite(params->inverse_fisheye_Rz{j}) && params->inverse_fisheye_Rz{j} > 0
   && isfinite(params->inverse_fisheye_w{j}) && params->inverse_fisheye_w{j} > 0)) {{
  fprintf(stderr, "Invalid inverse fisheye physical parameters at layer {j}.\\n"); exit(1);
}}
"""
    body = checks + f"""
const REAL L=params->grid_physical_size;
const REAL equator[3]={{L,0,0}}, pole[3]={{0,0,L}};
REAL e, p;
params->inverse_fisheye_c=1;
inverse_fisheye_scale{suffix}(params, equator, 0, &e, NULL, NULL, NULL);
inverse_fisheye_scale{suffix}(params, pole, 0, &p, NULL, NULL, NULL);
params->inverse_fisheye_c=1/e;
params->xmin=params->ymin=-L;
params->xmax=params->ymax=L;
params->zmax=L*p/e;
params->zmin=-params->zmax;
if (!(isfinite(params->zmax) && params->zmax > 0 && isfinite(params->inverse_fisheye_c))) {{
  fprintf(stderr, "Inverse fisheye normalization failed.\\n"); exit(1);
}}
"""
    cfc.register_CFunction(
        includes=includes,
        desc="Validate physical parameters and normalize equal physical axis extents.",
        cfunc_type="void",
        name="inverse_fisheye_normalize" + suffix,
        params="params_struct *restrict params",
        include_CodeParameters_h=False,
        body=body,
    )


def register_precompute(
    CoordSystem: str,
    num_transitions: int,
) -> None:
    """
    Register the grid loop storing the metric and two computational derivatives.

    :param CoordSystem: Selected inverse spheroidal coordinate system.
    :param num_transitions: Number of physical transitions.
    """
    register_CFunctions(num_transitions)
    suffix = f"__N{num_transitions}"
    # ReferenceMetric owns the covariant metric registrations.
    if "ghatUU00" not in gri.glb_gridfcs_dict:
        gri.register_gridfunctions_for_single_rankN(
            "ghatUU", rank=2, symmetry="sym01", group="AUXEVOL"
        )
    if "detgammahat" not in gri.glb_gridfcs_dict:
        gri.register_gridfunctions(["detgammahat"], group="AUXEVOL")
    assignments = []
    for i in range(3):
        for j in range(i, 3):
            assignments.extend(
                [
                    (f"ghatDD{i}{j}", f"g[{i}][{j}]"),
                    (
                        f"ghatUU{i}{j}",
                        f"J[{i}][0]*J[{j}][0]+J[{i}][1]*J[{j}][1]+J[{i}][2]*J[{j}][2]",
                    ),
                ]
            )
            for k in range(3):
                assignments.append((f"ghatDDdD{i}{j}{k}", f"dg[{i}][{j}][{k}]"))
                for l in range(k, 3):
                    assignments.append(
                        (f"ghatDDdDD{i}{j}{k}{l}", f"ddg[{i}][{j}][{k}][{l}]")
                    )
    assignments.append(("detgammahat", "detK*detK"))
    stores = ""
    for name, value in assignments:
        access = gri.BHaHGridFunction.access_gf(
            name, 0, 0, 0, gf_array_name="auxevol_gfs"
        )
        stores += f"{access} = {value};\n"
    body = f"""
const int Nxx_plus_2NGHOSTS0=params->Nxx_plus_2NGHOSTS0;
const int Nxx_plus_2NGHOSTS1=params->Nxx_plus_2NGHOSTS1;
const int Nxx_plus_2NGHOSTS2=params->Nxx_plus_2NGHOSTS2;
for (int i2=0; i2<Nxx_plus_2NGHOSTS2; i2++)
  for (int i1=0; i1<params->Nxx_plus_2NGHOSTS1; i1++)
    for (int i0=0; i0<params->Nxx_plus_2NGHOSTS0; i0++) {{
      const REAL q[3]={{xx[0][i0],xx[1][i1],xx[2][i2]}};
      REAL x[3], J[3][3], K[3][3], g[3][3], dg[3][3][3], ddg[3][3][3][3];
      if (inverse_fisheye_solve{suffix}(params, q, x)) {{
        fprintf(stderr, "Inverse fisheye precompute solve failed.\\n"); exit(1);
      }}
      inverse_fisheye_geometry{suffix}(params, x, true, J, K, g, dg, ddg);
      const REAL detK=K[0][0]*(K[1][1]*K[2][2]-K[1][2]*K[2][1])
                    -K[0][1]*(K[1][0]*K[2][2]-K[1][2]*K[2][0])
                    +K[0][2]*(K[1][0]*K[2][1]-K[1][1]*K[2][0]);
      {stores}
    }}
"""
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc="Precompute inverse spheroidal reference geometry, including ghost zones.",
        cfunc_type="void",
        name="generalrfm_precompute",
        CoordSystem_for_wrapper_func=CoordSystem,
        params="const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL *restrict xx[3], REAL *restrict auxevol_gfs",
        include_CodeParameters_h=False,
        body=body,
    )


def spacing_body(num_transitions: int, output: str) -> str:
    """
    Return code for minimum physical spacing along computational coordinate axes.

    :param num_transitions: Number of physical transitions.
    :param output: C pointer receiving the minimum coordinate-line spacing.
    :return: Generated function body, with no coordinate-orthogonality assumption.
    """
    register_CFunctions(num_transitions)
    suffix = f"__N{num_transitions}"
    return f"""
const REAL q[3]={{xx0,xx1,xx2}};
REAL x[3], J[3][3], K[3][3], g[3][3], dg[3][3][3], ddg[3][3][3][3];
if (inverse_fisheye_solve{suffix}(params, q, x)) {{
  fprintf(stderr, "Inverse fisheye spacing solve failed.\\n"); exit(1);
}}
inverse_fisheye_geometry{suffix}(params, x, false, J, K, g, dg, ddg);
const REAL d[3]={{params->dxx0,params->dxx1,params->dxx2}};
REAL minimum=INFINITY;
for (int i=0; i<3; i++) {{
  const REAL spacing=fabs(d[i])*hypot(hypot(K[0][i], K[1][i]), K[2][i]);
  minimum=NRPYMIN(minimum,spacing);
}}
*{output}=minimum;
"""


def register_CFunction_cfl_limited_timestep(set_of_CoordSystems: Set[str]) -> None:
    r"""
    Register an opt-in conservative CFL spacing bound without a grid scan.

    Write q = c*x*(1/aN + sum_j delta_j*f_j(rho_j)), with
    delta_j = 1/a_(j-1) - 1/a_j >= 0, f_j = G_j/rho_j, and
    B_j = diag(1, 1, (R_j/Rz_j)^2). The transition derivative G'_j decreases
    from 1 to 0, so 0 <= G'_j <= f_j <= 1. For rho_j > 0, define
    P_j = x*(B_j*x)^T/rho_j^2. Since P_j^2 = P_j, both P_j and I-P_j
    have Euclidean operator norm t_j(x) = |x|*|B_j*x|/(x^T B_j x), bounded
    by t_j = (k_j + 1/k_j)/2, where k_j = R_j/Rz_j.

    Each layer Jacobian is G'_j*I + (f_j-G'_j)*(I-P_j), whose norm is
    at most t_j. Thus ||dq/dx||_2 <= c*(1/aN + sum_j delta_j*t_j) = C.
    For K = dx/dq, ||K*e_i|| >= 1/C, so min_i(|dxx_i|)/C is a lower
    bound on every coordinate-line spacing used by the sampled CFL rule.
    The origin follows by continuity. The bound covers ghost zones too.
    It is exact for the identity map and uses only O(number of layers) work.
    It can reduce the timestep compared with the sampled minimum; it does
    not establish stability beyond the existing CFL rule.

    :param set_of_CoordSystems: Inverse spheroidal coordinate systems in the project.
    :raises ValueError: If a coordinate system is not inverse spheroidal fisheye.
    """
    body = "REAL bound;\nswitch (params->CoordSystem_hash) {\n"
    prefix = "GeneralRFM_inverse_fisheyeN"
    for coord in sorted(set_of_CoordSystems):
        if not coord.startswith(prefix) or not coord[len(prefix) :].isdigit():
            raise ValueError("Analytic CFL bound requires inverse spheroidal fisheye.")
        n = int(coord[len(prefix) :])
        body += (
            f"case {coord.upper()}: {{\nREAL norm = 1/params->inverse_fisheye_a{n};\n"
        )
        for j in range(1, n + 1):
            previous = "1" if j == 1 else f"1/params->inverse_fisheye_a{j-1}"
            body += f"""{{
const REAL k = params->inverse_fisheye_R{j}/params->inverse_fisheye_Rz{j};
norm += ({previous}-1/params->inverse_fisheye_a{j})*(k+1/k)/2;
}}
"""
        body += "bound = params->inverse_fisheye_c*norm;\nbreak;\n}\n"
    body += r"""
default: fprintf(stderr, "Unsupported coordinate system in inverse fisheye CFL bound.\n"); exit(1);
}
const REAL ds = NRPYMIN(fabs(params->dxx0), NRPYMIN(fabs(params->dxx1), fabs(params->dxx2)))/bound;
commondata->dt = NRPYMIN(commondata->dt, commondata->CFL_FACTOR*ds);
"""
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc="Conservative inverse-fisheye CFL bound, without global grid traversal.",
        cfunc_type="void",
        name="cfl_limited_timestep",
        params="commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3]",
        include_CodeParameters_h=False,
        body=body,
    )
