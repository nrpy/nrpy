# nrpy/infrastructures/BHaH/xx_tofrom_Cart.py
"""
C function registration for converting between uniform-grid coordinates (xx0, xx1, xx2) and Cartesian coordinates (x, y, z).

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Set, Tuple, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
import nrpy.reference_metric as refmetric


def _prepare_sympy_exprs_for_codegen(
    sympy_exprs: List[sp.Expr], local_vars: Set[str]
) -> List[sp.Expr]:
    """
    Substitute symbolic parameters with 'params->' prefix for C code generation.

    This function iterates through a list of SymPy expressions, identifies all unique
    free symbols that are not specified as local variables, and prepends 'params->'
    to their names. This is a common requirement for generating C code that
    accesses parameters from a 'params' struct.

    This implementation is efficient, creating a single substitution dictionary
    for all symbols across all expressions before applying it.

    :param sympy_exprs: A list of SymPy expressions to process.
    :param local_vars: A set of strings, where each string is the name of a
                       symbol that should be treated as a local variable (and thus
                       not be prefixed).
    :return: A new list of SymPy expressions with parameters appropriately substituted.

    Doctests:
    >>> import sympy as sp
    >>> # Define input symbols and expressions
    >>> a, b, c, xx0, xx1 = sp.symbols("a b c xx0 xx1")
    >>> input_exprs = [a + b*xx0, c - xx1]
    >>> local_vars = {"xx0", "xx1"}
    >>>
    >>> # Manually construct the expected SymPy objects for a robust comparison.
    >>> # The symbols "params->a", etc., are valid SymPy symbols.
    >>> params_a, params_b, params_c = sp.symbols("params->a params->b params->c")
    >>> expected_exprs_set = {xx0*params_b + params_a, -xx1 + params_c}
    >>>
    >>> # Generate the expressions and convert to a set for order-independent comparison.
    >>> generated_exprs_set = set(_prepare_sympy_exprs_for_codegen(input_exprs, local_vars))
    >>>
    >>> # The sets of expressions must be mathematically identical.
    >>> generated_exprs_set == expected_exprs_set
    True
    """
    # Create a single substitution dictionary for all unique parameters.
    all_symbols = set().union(*(expr.free_symbols for expr in sympy_exprs))
    substitutions: Dict[sp.Basic, sp.Basic] = {}

    for sym in all_symbols:
        # We explicitly cast Basic -> Symbol to access .name safely
        symbol = cast(sp.Symbol, sym)

        if symbol.name not in local_vars:
            # We explicitly type the replacement as a Symbol (which is an Expr)
            substitutions[symbol] = sp.symbols(f"params->{symbol.name}")

    # Apply exact symbol replacements without triggering algebraic rewriting.
    return [expr.xreplace(substitutions) for expr in sympy_exprs]


def _generate_bracketed_radial_inverse_body(
    radius_symbol: sp.Symbol,
    rbar_expr: sp.Expr,
    drbar_dr_expr: sp.Expr,
    asymptotic_scale_expr: sp.Expr,
    cart_components: Tuple[str, str, str],
    origin_body: str,
    success_body: str,
    failure_body: str,
) -> str:
    """
    Generate a bracketed radial inverse for monotone fisheye maps.

    :param radius_symbol: Symbol representing the raw radial coordinate.
    :param rbar_expr: Provider-owned scaled fisheye radius map.
    :param drbar_dr_expr: Provider-owned radial derivative of the scaled radius map.
    :param asymptotic_scale_expr: Provider-owned asymptotic scale for the far-field map.
    :param cart_components: Cartesian component expressions `(x, y, z)`.
    :param origin_body: C statements used when `rCart` is effectively zero.
    :param success_body: C statements used after a successful inverse solve.
    :param failure_body: C statements used when the inverse solve fails.
    :return: C code string for the inverse solve.
    """
    cartx, carty, cartz = cart_components

    def emit_codegen(
        radius_var: str,
        outputs: Tuple[Tuple[str, sp.Expr], ...],
        cse_varprefix: str,
    ) -> str:
        local_vars = {"rCart", "high", "radial_seed", "trial_seed", radius_var}
        local_radius = sp.Symbol(radius_var, real=True, nonnegative=True)
        processed_exprs = _prepare_sympy_exprs_for_codegen(
            [expr.xreplace({radius_symbol: local_radius}) for _, expr in outputs],
            local_vars,
        )
        return ccg.c_codegen(
            processed_exprs,
            [name for name, _ in outputs],
            include_braces=False,
            verbose=False,
            cse_varprefix=cse_varprefix,
        )

    asymptotic_scale_codegen = emit_codegen(
        "radial_seed",
        (("asymptotic_scale", asymptotic_scale_expr),),
        "asymptotic_",
    )
    high_map_codegen = emit_codegen(
        "high",
        (("high_map", rbar_expr),),
        "high_",
    )
    radial_map_codegen = emit_codegen(
        "radial_seed",
        (("radial_map", rbar_expr), ("radial_map_prime", drbar_dr_expr)),
        "radial_",
    )
    trial_map_codegen = emit_codegen(
        "trial_seed",
        (("trial_map", rbar_expr),),
        "trial_",
    )
    fallback_map_codegen = emit_codegen(
        "trial_seed",
        (("trial_map_fallback", rbar_expr),),
        "fallback_",
    )

    return f"""
  const REAL rCart = sqrt(({cartx}) * ({cartx}) + ({carty}) * ({carty}) + ({cartz}) * ({cartz}));
  if (!(isfinite(rCart))) {{
{failure_body}
  }} else if (rCart <= (REAL)1.0e-15) {{
{origin_body}
  }} else {{
    const REAL residual_tolerance = (REAL)1.0e-12 * NRPYMAX((REAL)1.0, rCart);
    const REAL bracket_tolerance = (REAL)1.0e-12 * NRPYMAX((REAL)1.0, rCart);
    REAL asymptotic_scale;
{asymptotic_scale_codegen}
    const REAL inv_asymptotic_scale =
        (fabs(asymptotic_scale) > (REAL)1.0e-15) ? (REAL)1.0 / asymptotic_scale : (REAL)1.0;
    REAL low = (REAL)0.0;
    REAL high = NRPYMAX(rCart * inv_asymptotic_scale, (REAL)1.0e-15);
    REAL radial_seed = (REAL)0.5 * high;
    int bracket_found = 0;
    int converged = 0;
    for (int expand = 0; expand < 80; expand++) {{
      REAL high_map;
{high_map_codegen}
      const REAL high_residual = high_map - rCart;
      if (isfinite(high_residual) && high_residual >= (REAL)0.0) {{
        bracket_found = 1;
        break;
      }}
      high = NRPYMAX(high * (REAL)2.0, (REAL)1.0);
    }} // END LOOP: for expand over bracket expansions
    if (!bracket_found) {{
{failure_body}
    }}
    radial_seed = (REAL)0.5 * (low + high);
    for (int iter = 0; iter < 80; iter++) {{
      REAL radial_map;
      REAL radial_map_prime;
{radial_map_codegen}
      const REAL radial_residual = radial_map - rCart;
      REAL trial_seed = (REAL)0.5 * (low + high);
      if (isfinite(radial_map_prime) && fabs(radial_map_prime) > (REAL)1.0e-14) {{
        const REAL newton_seed = radial_seed - radial_residual / radial_map_prime;
        if (isfinite(newton_seed) && newton_seed > low && newton_seed < high) {{
          trial_seed = newton_seed;
        }}
      }}
      REAL trial_map;
{trial_map_codegen}
      REAL trial_residual = trial_map - rCart;
      if (!isfinite(trial_residual)) {{
        trial_seed = (REAL)0.5 * (low + high);
        REAL trial_map_fallback;
{fallback_map_codegen}
        trial_residual = trial_map_fallback - rCart;
        if (!isfinite(trial_residual)) {{
{failure_body}
        }}
      }}
      if (trial_residual >= (REAL)0.0) {{
        high = trial_seed;
      }} else {{
        low = trial_seed;
      }}
      if (fabs(trial_residual) < residual_tolerance &&
          (fabs(high - low) < bracket_tolerance || fabs(trial_seed - radial_seed) < bracket_tolerance)) {{
        radial_seed = trial_seed;
        converged = 1;
        break;
      }}
      radial_seed = trial_seed;
    }} // END LOOP: for iter over bracketed Newton iterations
    if (!converged || !isfinite(radial_seed) || radial_seed < (REAL)0.0) {{
{failure_body}
    }}
{success_body}
  }} // END ELSE: invert fisheye radius away from the origin
"""


def register_CFunction__Cart_to_xx_and_nearest_i0i1i2(
    CoordSystem: str,
    relative_to: str = "local_grid_center",
    gridding_approach: str = "independent grid(s)",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register a C function that maps Cartesian coordinates to xx and finds the nearest grid indices.

    This function generates a C function which, given Cartesian coordinates (x, y, z),
    computes the corresponding (xx0, xx1, xx2) coordinates and determines the "closest"
    grid indices (i0, i1, i2) for the specified coordinate system. The C function is
    then registered for later use.

    :param CoordSystem: The coordinate system for the local grid patch.
    :param relative_to: Whether the computation is relative to the "local_grid_center"
                        (default) or "global_grid_center".
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :raises ValueError: When the value of `gridding_approach` is not "independent grid(s)"
                        or "multipatch".
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings, clang_format
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import supported_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "Cart_to_xx_and_nearest_i0i1i2"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in supported_CoordSystems:
    ...       if parallelization == "cuda" and CoordSystem.startswith("GeneralRFM"):
    ...          continue
    ...       cfc.CFunction_dict.clear()
    ...       _ = register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
    ...       generated_str = clang_format(cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function)
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}"
    ...       _ = validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    Setting up reference_metric[Spherical]...
    Setting up reference_metric[SinhSpherical]...
    Setting up reference_metric[SinhSphericalv2n2]...
    Setting up reference_metric[Cartesian]...
    Setting up reference_metric[SinhCartesian]...
    Setting up reference_metric[Cylindrical]...
    Setting up reference_metric[SinhCylindrical]...
    Setting up reference_metric[SinhCylindricalv2n2]...
    Setting up reference_metric[SymTP]...
    Setting up reference_metric[SinhSymTP]...
    Setting up reference_metric[LWedgeHSinhSph]...
    Setting up reference_metric[UWedgeHSinhSph]...
    Setting up reference_metric[RingHoleySinhSpherical]...
    Setting up reference_metric[HoleySinhSpherical]...
    Setting up reference_metric[GeneralRFM_fisheyeN2]...
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 1: Basic setup and parameter validation.
    if gridding_approach not in {"independent grid(s)", "multipatch"}:
        raise ValueError(
            "Invalid value for 'gridding_approach'. Must be 'independent grid(s)' or 'multipatch'."
        )

    parallelization = par.parval_from_str("parallelization")

    rfm = refmetric.reference_metric[CoordSystem]
    rfm_obj = rfm
    is_generalrfm = CoordSystem.startswith("GeneralRFM")
    provider_name = getattr(rfm, "general_rfm_provider_name", "")
    provider = getattr(rfm, "general_rfm_provider", None)
    is_fisheye_provider = is_generalrfm and provider_name == "fisheye"
    if is_generalrfm and not is_fisheye_provider:
        raise ValueError(
            f"GeneralRFM provider '{provider_name}' for {CoordSystem} is not yet supported in Cart_to_xx_and_nearest_i0i1i2."
        )
    if parallelization == "cuda" and is_generalrfm:
        raise ValueError(
            "GeneralRFM Cart_to_xx_and_nearest_i0i1i2 does not support CUDA parallelization."
        )
    local_C_vars = {"xx0", "xx1", "xx2", "Cartx", "Carty", "Cartz"} | (
        {"r", "rCart"} if is_fisheye_provider else set()
    )

    namesuffix = f"_{relative_to}" if relative_to == "global_grid_center" else ""
    name = f"Cart_to_xx_and_nearest_i0i1i2{namesuffix}"
    params = "const params_struct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]"
    cfunc_decorators = "__host__ __device__" if parallelization == "cuda" else ""

    desc = "Given Cartesian point (x,y,z), this function "
    if gridding_approach == "multipatch":
        desc += "assumes any required multipatch preprocessing has already been applied, then "
    desc += """unshifts the grid back to the origin to output the corresponding
            (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid"""

    # Step 2: Generate the core C-code for the coordinate transformation.
    core_body_list: List[str] = []
    if is_fisheye_provider:
        # ---------------------------------------------------------------------
        # Fisheye inverse map (Cart -> xx):
        #
        #   Cart^i = (rbar(r)/r) * xx^i,   with r = ||xx|| and rbar(r) monotone.
        #
        # Taking norms gives:
        #   rCart = ||Cart|| = rbar(r)
        #
        # So we solve the 1D equation rbar(r) - rCart = 0 for r using
        # safeguarded bracketed Newton iterations,
        # then recover:
        #   xx^i = (r / rCart) * Cart^i    (with the rCart=0 limit giving xx=0).
        # ---------------------------------------------------------------------
        if provider is None:
            raise ValueError(f"GeneralRFM provider object missing for {CoordSystem}.")
        r_local = sp.Symbol("r", real=True, nonnegative=True)
        rbar_expr, drbar_dr_expr, _, _ = provider.radius_map_and_derivs_for_inverse(
            r_local
        )
        asymptotic_scale_expr = provider.c * provider.a_list[-1]
        origin_body = """    xx[0] = (REAL)0.0;
    xx[1] = (REAL)0.0;
    xx[2] = (REAL)0.0;"""
        success_body = """    const REAL scale = radial_seed / rCart;
    xx[0] = scale * Cartx;
    xx[1] = scale * Carty;
    xx[2] = scale * Cartz;"""
        failure_body = f"""      fprintf(stderr, "ERROR: bracketed inverse failed for {CoordSystem} (fisheye): rCart, x,y,z = %.15e %.15e %.15e %.15e\\n",
              (double)rCart, (double)Cartx, (double)Carty, (double)Cartz);
      exit(1);"""
        fisheye_body = _generate_bracketed_radial_inverse_body(
            r_local,
            rbar_expr,
            drbar_dr_expr,
            asymptotic_scale_expr,
            ("Cartx", "Carty", "Cartz"),
            origin_body,
            success_body,
            failure_body,
        )
        core_body_list.append(fisheye_body)

    elif rfm_obj.requires_NewtonRaphson_for_Cart_to_xx:
        # Step 2.a: Handle mixed analytical and Newton-Raphson inversions.
        analytic_exprs: List[sp.Expr] = []
        analytic_names: List[str] = []
        for i in range(3):
            if rfm_obj.NewtonRaphson_f_of_xx[i] == sp.sympify(0):
                analytic_exprs.append(rfm_obj.Cart_to_xx[i])
                analytic_names.append(f"xx[{i}]")

        if analytic_exprs:
            core_body_list.append(
                "  // First compute analytical coordinate inversions:\n"
            )
            processed_analytic_exprs = _prepare_sympy_exprs_for_codegen(
                analytic_exprs, local_C_vars
            )
            core_body_list.append(
                ccg.c_codegen(
                    processed_analytic_exprs, analytic_names, include_braces=False
                )
            )

        core_body_list.append("""
  // Next perform Newton-Raphson iterations as needed:
  const REAL XX_TOLERANCE = 1e-12;  // that's 1 part in 1e12 dxxi.
  const REAL F_OF_XX_TOLERANCE = 1e-12;  // tolerance of function for which we're finding the root.
  const int ITER_MAX = 100;
  int iter;
""")
        for i in range(3):
            if rfm_obj.NewtonRaphson_f_of_xx[i] != sp.sympify(0):
                nr_input_exprs = [
                    rfm_obj.NewtonRaphson_f_of_xx[i],
                    sp.diff(
                        rfm_obj.NewtonRaphson_f_of_xx[i],
                        rfm_obj.xx[i],
                    ),
                ]
                nr_processed_exprs = _prepare_sympy_exprs_for_codegen(
                    nr_input_exprs, local_C_vars
                )
                nr_codegen_output = ccg.c_codegen(
                    nr_processed_exprs,
                    [f"f_of_xx{i}", f"fprime_of_xx{i}"],
                    include_braces=True,
                    verbose=False,
                )
                core_body_list.append(f"""
  {{
  int tolerance_has_been_met = 0;
  iter = 0;
  REAL xx{i} = (REAL)0.5 * (params->xxmin{i} + params->xxmax{i});
  while(iter < ITER_MAX && !tolerance_has_been_met) {{
    REAL f_of_xx{i}, fprime_of_xx{i};

{nr_codegen_output}
    if(fprime_of_xx{i} == (REAL)0.0) {{
      break;
    }}
    const REAL xx{i}_np1 = xx{i} - f_of_xx{i} / fprime_of_xx{i};

    if( fabs(xx{i} - xx{i}_np1) <= XX_TOLERANCE * params->dxx{i} && fabs(f_of_xx{i}) <= F_OF_XX_TOLERANCE ) {{
      tolerance_has_been_met = 1;
    }}
    xx{i} = xx{i}_np1;
    iter++;
  }} // END Newton-Raphson iterations to compute xx{i}
  if(iter >= ITER_MAX || !tolerance_has_been_met) {{
#ifdef __CUDA_ARCH__
    printf("ERROR: Newton-Raphson failed for {CoordSystem}: xx{i}=%.15e, x,y,z = %.15e %.15e %.15e\\n",
           (double)xx{i}, (double)Cartx, (double)Carty, (double)Cartz);
    asm("trap;");
#else
    fprintf(stderr, "ERROR: Newton-Raphson failed for {CoordSystem}: xx{i}=%.15e, x,y,z = %.15e %.15e %.15e\\n",
            (double)xx{i}, (double)Cartx, (double)Carty, (double)Cartz);
    exit(1);
#endif
  }}
  xx[{i}] = xx{i};
  }}
""")
    else:
        # Step 2.b: Handle purely analytical inversions.
        processed_exprs = _prepare_sympy_exprs_for_codegen(
            list(rfm_obj.Cart_to_xx), local_C_vars
        )
        core_body_list.append(
            ccg.c_codegen(
                processed_exprs,
                ["xx[0]", "xx[1]", "xx[2]"],
                include_braces=False,
            )
        )

    # Step 2.c: Add logic to find the nearest grid point index.
    core_body_list.append("""
      // Find the nearest grid indices (i0, i1, i2) for the given Cartesian coordinates (x, y, z).
      // Assuming a cell-centered grid, which follows the pattern:
      //   xx0[i0] = params->xxmin0 + ((REAL)(i0 - NGHOSTS) + 0.5) * params->dxx0
      // The index i0 can be derived as:
      //   i0 = (xx0[i0] - params->xxmin0) / params->dxx0 - 0.5 + NGHOSTS
      // Now, including typecasts:
      //   i0 = (int)((xx[0] - params->xxmin0) / params->dxx0 - 0.5 + (REAL)NGHOSTS)
      // C float-to-int conversion truncates toward zero; for nonnegative inputs this matches floor().
      // Assuming (xx - xxmin)/dxx + NGHOSTS is nonnegative (typical for valid interior points), this is safe.
      //   i0 = (int)((xx[0] - params->xxmin0) / params->dxx0 - 0.5 + (REAL)NGHOSTS + 0.5)
      // The 0.5 values cancel out:
      //   i0 =           (int)( ( xx[0] - params->xxmin0 ) / params->dxx0 + (REAL)NGHOSTS )
      Cart_to_i0i1i2[0] = (int)( ( xx[0] - params->xxmin0 ) / params->dxx0 + (REAL)NGHOSTS );
      Cart_to_i0i1i2[1] = (int)( ( xx[1] - params->xxmin1 ) / params->dxx1 + (REAL)NGHOSTS );
      Cart_to_i0i1i2[2] = (int)( ( xx[2] - params->xxmin2 ) / params->dxx2 + (REAL)NGHOSTS );
""")
    core_body = "".join(core_body_list)

    # Step 3: Assemble the full C function body.
    body_parts = ["""
  // Set (Cartx, Carty, Cartz) relative to the global (as opposed to local) grid.
  //   This local grid may be offset from the origin by adjusting
  //   (Cart_originx, Cart_originy, Cart_originz) to nonzero values.
  REAL Cartx = xCart[0];
  REAL Carty = xCart[1];
  REAL Cartz = xCart[2];
"""]
    if relative_to == "local_grid_center":
        body_parts.append("""
  // Set the origin, (Cartx, Carty, Cartz) = (0, 0, 0), to the center of the local grid patch.
  Cartx -= params->Cart_originx;
  Carty -= params->Cart_originy;
  Cartz -= params->Cart_originz;
  {
""")
        body_parts.append(core_body)
        body_parts.append("  }\n")
    else:
        body_parts.append(core_body)

    # Step 4: Register the C function.
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc=desc,
        cfunc_type="void",
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body="".join(body_parts),
        cfunc_decorators=cfunc_decorators,
    )
    return pcg.NRPyEnv()


def register_CFunction_xx_to_Cart(
    CoordSystem: str,
    gridding_approach: str = "independent grid(s)",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Convert uniform-grid coordinate (xx[0], xx[1], xx[2]) to the corresponding Cartesian coordinate.

    :param CoordSystem: The coordinate system name as a string.
    :param gridding_approach: Choices: "independent grid(s)" (default) or "multipatch".
    :raises ValueError: If an invalid gridding_approach is provided.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings, clang_format
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import supported_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "xx_to_Cart"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in supported_CoordSystems:
    ...       if parallelization == "cuda" and CoordSystem.startswith("GeneralRFM"):
    ...          continue
    ...       cfc.CFunction_dict.clear()
    ...       _ = register_CFunction_xx_to_Cart(CoordSystem)
    ...       generated_str = clang_format(cfc.CFunction_dict[f'{name}__rfm__{CoordSystem}'].full_function)
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}"
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 1: Basic setup and parameter validation.
    if gridding_approach not in {"independent grid(s)", "multipatch"}:
        raise ValueError(
            "Invalid value for 'gridding_approach'. Must be 'independent grid(s)' or 'multipatch'."
        )

    parallelization = par.parval_from_str("parallelization")

    rfm = refmetric.reference_metric[CoordSystem]
    is_generalrfm = CoordSystem.startswith("GeneralRFM")
    provider_name = getattr(rfm, "general_rfm_provider_name", "")
    provider = getattr(rfm, "general_rfm_provider", None)
    is_fisheye_provider = is_generalrfm and provider_name == "fisheye"
    if is_generalrfm and not is_fisheye_provider:
        raise ValueError(
            f"GeneralRFM provider '{provider_name}' for {CoordSystem} is not yet supported in xx_to_Cart."
        )
    if parallelization == "cuda" and is_generalrfm:
        raise ValueError("GeneralRFM xx_to_Cart does not support CUDA parallelization.")

    local_C_vars = {"xx0", "xx1", "xx2"} | ({"r"} if is_fisheye_provider else set())

    # Step 2: Prepare SymPy expressions for C code generation.
    if is_fisheye_provider:
        if provider is None:
            raise ValueError(f"GeneralRFM provider object missing for {CoordSystem}.")
        r_local = sp.Symbol("r", real=True, nonnegative=True)
        rbar_unscaled, _, _, _ = provider.radius_map_unscaled_and_derivs_closed_form(
            r_local
        )
        rbar_expr = provider.c * rbar_unscaled
        processed_expr = _prepare_sympy_exprs_for_codegen(
            [rbar_expr], local_C_vars | {"r"}
        )
        rbar_codegen = ccg.c_codegen(
            processed_expr, ["rbar"], include_braces=False, verbose=False
        )

        body = f"""
const REAL xx0 = xx[0];
const REAL xx1 = xx[1];
const REAL xx2 = xx[2];
const REAL r2 = xx0*xx0 + xx1*xx1 + xx2*xx2;

if(r2 <= (REAL)0.0) {{
  // Fisheye map sends the origin to the origin (plus any patch offset).
  xCart[0] = params->Cart_originx;
  xCart[1] = params->Cart_originy;
  xCart[2] = params->Cart_originz;
}} else {{
  const REAL r = sqrt(r2);
  REAL rbar;
{rbar_codegen}

  const REAL lam = rbar / r;

  xCart[0] = lam*xx0 + params->Cart_originx;
  xCart[1] = lam*xx1 + params->Cart_originy;
  xCart[2] = lam*xx2 + params->Cart_originz;
}}
"""
    else:
        raw_xx_to_Cart_exprs = [
            rfm.xx_to_Cart[i] + gri.Cart_origin[i] for i in range(3)
        ]
        processed_exprs = _prepare_sympy_exprs_for_codegen(
            raw_xx_to_Cart_exprs, local_C_vars
        )

        # Step 3: Construct the full C function body.
        codegen_results = ccg.c_codegen(
            processed_exprs,
            ["xCart[0]", "xCart[1]", "xCart[2]"],
        )
        body = """
const REAL xx0 = xx[0];
const REAL xx1 = xx[1];
const REAL xx2 = xx[2];
""" + codegen_results

    # Step 4: Register the C function.
    cfc.register_CFunction(
        includes=["BHaH_defines.h"],
        desc="""Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} from the
local coordinate vector {xx[0], xx[1], xx[2]} = {xx0, xx1, xx2},
taking into account the possibility that the origin of this grid is off-center.""",
        cfunc_type="void",
        CoordSystem_for_wrapper_func=CoordSystem,
        name="xx_to_Cart",
        params="const params_struct *restrict params, const REAL xx[3], REAL xCart[3]",
        include_CodeParameters_h=False,
        body=body,
        cfunc_decorators="__host__ __device__" if parallelization == "cuda" else "",
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
