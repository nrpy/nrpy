"""
Register a C function that exports Cartesian-output-basis raytracing data.

This module generates a BHaH diagnostics helper that, on each scheduled output
step, evaluates the Cartesian covariant four-metric from the native BSSN
evolution state using the transformed symbolic recipe in
``nrpy.equations.general_relativity.geodesics.geodesics``, and writes a stable
binary payload for later raytracing.

The payload stores the final Cartesian metric directly. The raytracer derives
metric derivatives and Christoffel symbols from the stored metric data.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Sequence, Tuple, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)
from nrpy.infrastructures import BHaH

SUPPORTED_RAYTRACING_DATA_MODES = ("g4DD", "g4DD_d0", "GammaUDD", "all")

METRIC_COMPONENTS: List[Tuple[Tuple[int, int], str]] = [
    ((0, 0), "g4DD00"),
    ((0, 1), "g4DD01"),
    ((0, 2), "g4DD02"),
    ((0, 3), "g4DD03"),
    ((1, 1), "g4DD11"),
    ((1, 2), "g4DD12"),
    ((1, 3), "g4DD13"),
    ((2, 2), "g4DD22"),
    ((2, 3), "g4DD23"),
    ((3, 3), "g4DD33"),
]

METRIC_D0_COMPONENTS: List[Tuple[Tuple[int, int], str]] = [
    ((mu, nu), f"g4DDdt{mu}{nu}") for (mu, nu), _ in METRIC_COMPONENTS
]

CHRISTOFFEL_COMPONENTS: List[Tuple[Tuple[int, int, int], str]] = [
    ((alpha, mu, nu), f"Gamma4UDD{alpha}{mu}{nu}")
    for alpha in range(4)
    for mu in range(4)
    for nu in range(mu, 4)
]


def _validate_fixed_width_string(text: str, field_bytes: int, label: str) -> None:
    """
    Validate that a fixed-width on-disk string fits its serialized field.

    The binary writer stores null-padded strings and rejects any text whose
    encoded length would fill or exceed the field. Checking this at Python
    code-generation time is cheaper than discovering the mismatch during a run.

    :param text: String that will be serialized into the fixed-width field.
    :param field_bytes: Exact field width in bytes.
    :param label: Human-readable label for error reporting.
    :raises ValueError: If the string cannot fit in the fixed-width field.
    """
    if len(text.encode("ascii")) >= field_bytes:
        raise ValueError(
            f"String '{text}' for {label} exceeds fixed-width field size "
            f"{field_bytes}."
        )


def _build_file_writer_block(
    mode_name: str,
    magic: str,
    header_label: str,
    format_name: str,
    CoordSystem: str,
    cf_convention: str,
    record_names: Sequence[str],
    metric_names: Sequence[str],
    secondary_names: Sequence[str],
    include_christoffel_count: bool,
) -> str:
    """
    Generate one robust all-mode file-writing block.

    The block uses the same fixed header and atomic-install helpers as the
    single-mode writer. It is a code-generation helper, not a runtime output
    channel abstraction.

    :param mode_name: Method directory and generated C variable suffix.
    :param magic: Existing slice-file magic.
    :param header_label: Fixed header label.
    :param format_name: Fixed format description.
    :param CoordSystem: Source coordinate-system name.
    :param cf_convention: Evolved conformal-factor convention.
    :param record_names: Names stored in the record-component table.
    :param metric_names: Names stored in the metric-component table.
    :param secondary_names: Names stored in the secondary-component table.
    :param include_christoffel_count: Whether the header has a Gamma count.
    :return: Generated C source for one all-mode writer block.
    """
    block = rf"""
  {{
    const char magic_{mode_name}[16] = "{magic}";
    char final_filename_{mode_name}[256];
    char temporary_filename_{mode_name}[256];
    if (mkdir("raytracing_slices", (mode_t)0777) != 0 && errno != EEXIST)
      raytracing_data_abort_with_message(
          "Error: output_raytracing_data could not create raytracing_slices.");
    if (mkdir("raytracing_slices/{mode_name}", (mode_t)0777) != 0 && errno != EEXIST)
      raytracing_data_abort_with_message(
          "Error: output_raytracing_data could not create an all-mode method directory.");
    snprintf(final_filename_{mode_name}, sizeof(final_filename_{mode_name}),
             "./raytracing_slices/{mode_name}/raytracing_data_t%08d.bin",
             time_output_index);
    snprintf(temporary_filename_{mode_name}, sizeof(temporary_filename_{mode_name}),
             "./raytracing_slices/{mode_name}/raytracing_data_t%08d.bin.tmp.XXXXXX",
             time_output_index);
    const int temporary_fd_{mode_name} = mkstemp(temporary_filename_{mode_name});
    if (temporary_fd_{mode_name} == -1)
      raytracing_data_abort_with_message(
          "Error: output_raytracing_data could not create an all-mode temporary file.");
    FILE *restrict fp_{mode_name} = fdopen(temporary_fd_{mode_name}, "wb");
    if (fp_{mode_name} == NULL) {{
      close(temporary_fd_{mode_name});
      remove(temporary_filename_{mode_name});
      raytracing_data_abort_with_message(
          "Error: output_raytracing_data could not open an all-mode temporary file.");
    }} // END IF: all-mode temporary file could not be opened
    raytracing_data_write_or_abort(fp_{mode_name}, magic_{mode_name}, sizeof(char), 16, "magic");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, header_size_{mode_name}, "header_size");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, output_index, "output_index");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, num_grids, "num_grids");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, serialized_real_bytes, "serialized_real_bytes");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, record_component_count_{mode_name}, "record_component_count");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, metric_component_count_{mode_name}, "metric_component_count");
"""
    if include_christoffel_count:
        block += rf"""    raytracing_data_write_u32_or_abort(
        fp_{mode_name}, christoffel_component_count_{mode_name}, "christoffel_component_count");
"""
    block += rf"""    raytracing_data_write_u32_or_abort(fp_{mode_name}, point_record_real_count_{mode_name}, "point_record_real_count");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, point_record_bytes_{mode_name}, "point_record_bytes");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, payload_includes_ghost_zones, "payload_includes_ghost_zones");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, file_is_little_endian, "file_is_little_endian");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, time_variable_is_f64, "time_variable_is_f64");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, reserved_u32, "reserved_u32");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, Nxx0, "Nxx0");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, Nxx1, "Nxx1");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, Nxx2, "Nxx2");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, Nxx_plus_2NGHOSTS0_u32, "Nxx_plus_2NGHOSTS0");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, Nxx_plus_2NGHOSTS1_u32, "Nxx_plus_2NGHOSTS1");
    raytracing_data_write_u32_or_abort(fp_{mode_name}, Nxx_plus_2NGHOSTS2_u32, "Nxx_plus_2NGHOSTS2");
    raytracing_data_write_u64_or_abort(fp_{mode_name}, point_record_count, "point_record_count");
    raytracing_data_write_u64_or_abort(fp_{mode_name}, (uint64_t)header_size_{mode_name}, "simulation_time_offset");
    raytracing_data_write_u64_or_abort(fp_{mode_name}, (uint64_t)header_size_{mode_name} + 8ULL, "point_records_offset");
    raytracing_data_write_u64_or_abort(fp_{mode_name}, point_records_bytes_{mode_name}, "point_records_bytes");
    raytracing_data_write_u64_or_abort(fp_{mode_name}, (uint64_t)header_size_{mode_name} + 8ULL + point_records_bytes_{mode_name}, "total_file_bytes");
    raytracing_data_write_u64_or_abort(fp_{mode_name}, (uint64_t)NGHOSTS, "NGHOSTS");
    raytracing_data_write_u64_or_abort(fp_{mode_name}, (uint64_t)Nxx_plus_2NGHOSTS0_u32, "payload_i0_count");
    raytracing_data_write_u64_or_abort(fp_{mode_name}, (uint64_t)Nxx_plus_2NGHOSTS1_u32, "payload_i1_count");
    raytracing_data_write_u64_or_abort(fp_{mode_name}, (uint64_t)Nxx_plus_2NGHOSTS2_u32, "payload_i2_count");
    raytracing_data_write_i64_or_abort(fp_{mode_name}, payload_i0_start, "payload_i0_start");
    raytracing_data_write_i64_or_abort(fp_{mode_name}, payload_i1_start, "payload_i1_start");
    raytracing_data_write_i64_or_abort(fp_{mode_name}, payload_i2_start, "payload_i2_start");
    raytracing_data_write_i64_or_abort(fp_{mode_name}, payload_i0_end, "payload_i0_end");
    raytracing_data_write_i64_or_abort(fp_{mode_name}, payload_i1_end, "payload_i1_end");
    raytracing_data_write_i64_or_abort(fp_{mode_name}, payload_i2_end, "payload_i2_end");
    for (int i = 0; i < 3; i++) raytracing_data_write_f64_or_abort(fp_{mode_name}, dxx[i], "dxx");
    for (int i = 0; i < 3; i++) raytracing_data_write_f64_or_abort(fp_{mode_name}, invdxx[i], "invdxx");
    for (int i = 0; i < 3; i++) raytracing_data_write_f64_or_abort(fp_{mode_name}, xxmin[i], "xxmin");
    for (int i = 0; i < 3; i++) raytracing_data_write_f64_or_abort(fp_{mode_name}, xxmax[i], "xxmax");
    for (int i = 0; i < 3; i++) raytracing_data_write_f64_or_abort(fp_{mode_name}, cart_origin[i], "cart_origin");
    raytracing_data_write_fixed_length_string_or_abort(fp_{mode_name}, "{header_label}", 32, "header_label");
    raytracing_data_write_fixed_length_string_or_abort(fp_{mode_name}, "{format_name}", 32, "format_name");
    raytracing_data_write_fixed_length_string_or_abort(fp_{mode_name}, "Cartesian", 16, "target_basis");
    raytracing_data_write_fixed_length_string_or_abort(fp_{mode_name}, "{CoordSystem}", 32, "source_coord_system");
    raytracing_data_write_fixed_length_string_or_abort(fp_{mode_name}, "i2maj_i0fast", 16, "loop_order");
    raytracing_data_write_fixed_length_string_or_abort(fp_{mode_name}, "{cf_convention}", 32, "cf_convention");
"""
    for component_name in record_names:
        block += (
            "    raytracing_data_write_fixed_length_string_or_abort("
            f'fp_{mode_name}, "{component_name}", 24, "record_component_names");\n'
        )
    for component_name in metric_names:
        block += (
            "    raytracing_data_write_fixed_length_string_or_abort("
            f'fp_{mode_name}, "{component_name}", 16, "metric_component_names");\n'
        )
    if include_christoffel_count:
        for component_name in secondary_names:
            block += (
                "    raytracing_data_write_fixed_length_string_or_abort("
                f'fp_{mode_name}, "{component_name}", 16, "christoffel_component_names");\n'
            )
    block += rf"""    raytracing_data_write_f64_or_abort(fp_{mode_name}, (double)commondata->time, "simulation_time");
    raytracing_data_write_f64_array_or_abort(
        fp_{mode_name}, payload_{mode_name}, payload_value_count_{mode_name}, "point-record payload");
    if (fflush(fp_{mode_name}) != 0) {{
      fclose(fp_{mode_name});
      remove(temporary_filename_{mode_name});
      raytracing_data_abort_with_message(
          "Error: output_raytracing_data failed while flushing an all-mode file.");
    }} // END IF: all-mode file flush failed
    if (fclose(fp_{mode_name}) != 0) {{
      remove(temporary_filename_{mode_name});
      raytracing_data_abort_with_message(
          "Error: output_raytracing_data failed while closing an all-mode file.");
    }} // END IF: all-mode file close failed
    BHAH_FREE(payload_{mode_name});
    if (link(temporary_filename_{mode_name}, final_filename_{mode_name}) != 0) {{
      const int saved_errno = errno;
      if (remove(temporary_filename_{mode_name}) != 0 && saved_errno != EEXIST)
        raytracing_data_abort_with_message(
            "Error: output_raytracing_data could not clean up an all-mode file.");
      if (saved_errno != EEXIST)
        raytracing_data_abort_with_message(
            "Error: output_raytracing_data could not install an all-mode file.");
    }} else if (remove(temporary_filename_{mode_name}) != 0) {{
      raytracing_data_abort_with_message(
          "Error: output_raytracing_data could not remove an all-mode temporary file.");
    }} // END IF: all-mode file installation failed or succeeded
  }} // END BLOCK: write {mode_name} all-mode slice
"""
    return block


def register_CFunction_output_raytracing_data(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    raytracing_data_mode: str = "g4DD",
    enable_RbarDD_gridfunctions: bool = False,
    enable_static_christoffels: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register the Cartesian-output-basis raytracing binary exporter.

    The generated C helper writes one binary file per diagnostics output in a
    method-specific directory. Each file contains
    a fixed-size self-describing header, the physical simulation time, and one
    point record per full logical-grid point, including ghost zones. Each
    point record stores:

    1. Cartesian coordinates ``(x, y, z)``.
    2. The 10 unique covariant four-metric components ``g4DD`` in a fixed
       upper-triangular order.

    The exporter evaluates Cartesian metric expressions only on interior
    points. Cartesian coordinates are evaluated on the full logical grid.
    Ghost-zone tensor payload records are then filled using the existing BHaH
    boundary ordering: pure outer ghost zones are extrapolated first, and
    inner-boundary records are copied from their mapped source records in
    ``bcstruct``. The payload-local ordering remains zero-based storage order,
    while the serialized payload window exposes the corresponding logical-grid
    half-open interval ``[-NGHOSTS, Nxx + NGHOSTS)``.

    :param CoordSystem: Coordinate system used by the evolved BSSN state.
    :param enable_rfm_precompute: Whether the generated project uses
        reference-metric precomputation.
    :param raytracing_data_mode: Data written by the generated exporter.
    :param enable_RbarDD_gridfunctions: Whether the generated project registers
        the Ricci path needed by derivative/Christoffel output.
    :param enable_static_christoffels: Whether the final qualifying GammaUDD
        output may use static-spacetime Christoffels. The diagnostics scheduler
        supplies the per-call selection flag.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If ``CoordSystem`` is not a string.
    :raises ValueError: If ``enable_rfm_precompute`` is not a bool.
    :raises ValueError: If ``raytracing_data_mode`` is unsupported.
    :raises ValueError: If a method requiring RHS data is used without the
        required Ricci/RHS setup.
    :raises ValueError: If the generated project uses CUDA parallelization.
    """
    if not isinstance(CoordSystem, str):
        raise ValueError(
            f"CoordSystem must be a string, got {type(CoordSystem).__name__}"
        )
    if not isinstance(enable_rfm_precompute, bool):
        raise ValueError(
            "enable_rfm_precompute must be a bool, "
            f"got {type(enable_rfm_precompute).__name__}"
        )
    if raytracing_data_mode not in SUPPORTED_RAYTRACING_DATA_MODES:
        raise ValueError(
            f"raytracing_data_mode must be one of {SUPPORTED_RAYTRACING_DATA_MODES}, "
            f"got '{raytracing_data_mode}'."
        )
    if not isinstance(enable_RbarDD_gridfunctions, bool):
        raise ValueError(
            "enable_RbarDD_gridfunctions must be a bool, "
            f"got {type(enable_RbarDD_gridfunctions).__name__}"
        )
    if not isinstance(enable_static_christoffels, bool):
        raise ValueError(
            "enable_static_christoffels must be a bool, "
            f"got {type(enable_static_christoffels).__name__}"
        )
    if raytracing_data_mode in ("g4DD_d0", "GammaUDD", "all"):
        if not enable_rfm_precompute:
            raise ValueError(
                f"{raytracing_data_mode} output requires enable_rfm_precompute=True."
            )
        if not enable_RbarDD_gridfunctions:
            raise ValueError(
                f"{raytracing_data_mode} output requires "
                "enable_RbarDD_gridfunctions=True."
            )
    if enable_static_christoffels and raytracing_data_mode not in (
        "GammaUDD",
        "all",
    ):
        raise ValueError(
            "enable_static_christoffels is valid only for GammaUDD or all output."
        )
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    parallelization = cast(str, par.parval_from_str("parallelization"))
    if parallelization == "cuda":
        raise ValueError(
            "Raytracing binary exporters currently support only host/OpenMP builds."
        )

    # Import the BSSN quantities module so the conformal-factor parameter is registered.
    # pylint: disable-next=import-outside-toplevel
    __import__("nrpy.equations.general_relativity.BSSN_quantities")

    # pylint: disable-next=import-outside-toplevel
    from nrpy.equations.general_relativity.geodesics.geodesics import (
        symbolic_christoffel_recipe_from_bssn_grid_basis,
        symbolic_g4DD_dt_recipe_from_bssn_grid_basis,
        symbolic_g4DD_recipe_from_bssn_grid_basis,
    )

    # Step 1: Build the transformed symbolic recipes required by the selected mode.
    g4DD = symbolic_g4DD_recipe_from_bssn_grid_basis(
        bssn_coord_system=CoordSystem,
        target_basis="Cartesian",
        enable_bssn_rfm_precompute=enable_rfm_precompute,
    )
    g4DD_d0 = None
    Gamma4UDD = None
    Gamma4UDD_static = None
    if raytracing_data_mode in ("g4DD_d0", "all"):
        g4DD_d0 = symbolic_g4DD_dt_recipe_from_bssn_grid_basis(
            bssn_coord_system=CoordSystem,
            target_basis="Cartesian",
            enable_bssn_rfm_precompute=enable_rfm_precompute,
        )
    if raytracing_data_mode in ("GammaUDD", "all"):
        Gamma4UDD = symbolic_christoffel_recipe_from_bssn_grid_basis(
            bssn_coord_system=CoordSystem,
            target_basis="Cartesian",
            enable_bssn_rfm_precompute=enable_rfm_precompute,
        )
        if enable_static_christoffels:
            Gamma4UDD_static = symbolic_christoffel_recipe_from_bssn_grid_basis(
                bssn_coord_system=CoordSystem,
                target_basis="Cartesian",
                enable_bssn_rfm_precompute=enable_rfm_precompute,
                use_static_time_derivatives=True,
            )

    metric_components = METRIC_COMPONENTS
    secondary_components: Sequence[
        Union[Tuple[Tuple[int, int], str], Tuple[Tuple[int, int, int], str]]
    ] = ()
    if raytracing_data_mode == "g4DD_d0":
        secondary_components = METRIC_D0_COMPONENTS
    elif raytracing_data_mode == "GammaUDD":
        secondary_components = CHRISTOFFEL_COMPONENTS

    # Step 2: Flatten only the symbolic expressions required by the selected
    #         output path. The all-mode path must not build a second unused
    #         single-dataset interpolation kernel.
    expr_list: List[sp.Expr] = []
    lhs_list: List[str] = []
    all_codegen_expr_list: List[sp.Expr] = []
    all_codegen_lhs_list: List[str] = []
    static_gamma_expr_list: List[sp.Expr] = []
    static_gamma_lhs_list: List[str] = []
    if raytracing_data_mode in ("GammaUDD", "all") and Gamma4UDD_static is not None:
        static_gamma_expr_list = [
            Gamma4UDD_static[alpha][mu][nu]
            for (alpha, mu, nu), _ in CHRISTOFFEL_COMPONENTS
        ]
        static_gamma_lhs_list = [
            f"REAL Gamma4UDD_static{alpha}{mu}{nu}"
            for (alpha, mu, nu), _ in CHRISTOFFEL_COMPONENTS
        ]
    if raytracing_data_mode == "all":
        all_codegen_expr_list = [g4DD[mu][nu] for (mu, nu), _ in METRIC_COMPONENTS]
        all_codegen_lhs_list = [f"REAL {name}" for _, name in METRIC_COMPONENTS]
        assert g4DD_d0 is not None
        all_codegen_expr_list.extend(
            g4DD_d0[mu][nu] for (mu, nu), _ in METRIC_D0_COMPONENTS
        )
        all_codegen_lhs_list.extend(f"REAL {name}" for _, name in METRIC_D0_COMPONENTS)
        assert Gamma4UDD is not None
        all_codegen_expr_list.extend(
            Gamma4UDD[alpha][mu][nu] for (alpha, mu, nu), _ in CHRISTOFFEL_COMPONENTS
        )
        all_codegen_lhs_list.extend(
            f"REAL {name}" for _, name in CHRISTOFFEL_COMPONENTS
        )
        all_codegen_expr_list.extend(static_gamma_expr_list)
        all_codegen_lhs_list.extend(static_gamma_lhs_list)
    else:
        expr_list = [g4DD[mu][nu] for (mu, nu), _ in metric_components]
        lhs_list = [f"REAL {name}" for _, name in metric_components]
        if raytracing_data_mode == "g4DD_d0":
            assert g4DD_d0 is not None
            expr_list.extend(g4DD_d0[mu][nu] for (mu, nu), _ in METRIC_D0_COMPONENTS)
            lhs_list.extend(f"REAL {name}" for _, name in METRIC_D0_COMPONENTS)
        elif raytracing_data_mode == "GammaUDD":
            assert Gamma4UDD is not None
            expr_list.extend(
                Gamma4UDD[alpha][mu][nu]
                for (alpha, mu, nu), _ in CHRISTOFFEL_COMPONENTS
            )
            lhs_list.extend(f"REAL {name}" for _, name in CHRISTOFFEL_COMPONENTS)
            expr_list.extend(static_gamma_expr_list)
            lhs_list.extend(static_gamma_lhs_list)
    expression_list_for_definitions = (
        all_codegen_expr_list if raytracing_data_mode == "all" else expr_list
    )

    # Step 4: Ask the shared expression helpers which runtime parameters the
    #         symbolic recipes need, then generate the corresponding C locals.
    params_symbols, commondata_symbols = get_params_commondata_symbols_from_expr_list(
        expression_list_for_definitions
    )
    params_definitions = generate_definition_header(
        params_symbols,
        enable_intrinsics=False,
        var_access="params->",
    )
    commondata_definitions = generate_definition_header(
        commondata_symbols,
        enable_intrinsics=False,
        var_access="commondata->",
    )

    # Step 5: Define the fixed-width header schema.
    record_component_names = ["x", "y", "z"]
    record_component_names.extend(name for _, name in metric_components)
    record_component_names.extend(name for _, name in secondary_components)
    header_label = {
        "g4DD": "NRPy RT metric data",
        "g4DD_d0": "NRPy RT metric and dt data",
        "GammaUDD": "NRPy RT spacetime data",
        "all": "NRPy RT spacetime data",
    }[raytracing_data_mode]
    record_component_name_bytes = 24
    header_label_bytes = 32
    format_name = {
        "g4DD": "Cartesian g4DD",
        "g4DD_d0": "Cartesian g4DD+g4DDdt",
        "GammaUDD": "Cartesian g4DD+Gamma4UDD",
        "all": "Cartesian g4DD+Gamma4UDD",
    }[raytracing_data_mode]
    format_name_bytes = 32
    target_basis_name = "Cartesian"
    basis_name_bytes = 16
    coord_name_bytes = 32
    loop_order_name = "i2maj_i0fast"
    loop_order_name_bytes = 16
    metric_component_name_bytes = 16
    num_metric_components = len(metric_components)
    num_secondary_components = len(secondary_components)
    num_record_components = len(record_component_names)
    point_record_real_count = 3 + num_metric_components + num_secondary_components
    point_record_bytes = 8 * point_record_real_count
    num_u32_header_fields = 19 if raytracing_data_mode == "GammaUDD" else 18
    num_u64_header_fields = 15
    num_f64_header_fields = 15
    cf_convention_bytes = format_name_bytes
    header_size_bytes = (
        16
        + num_u32_header_fields * 4
        + num_u64_header_fields * 8
        + num_f64_header_fields * 8
        + header_label_bytes
        + format_name_bytes
        + basis_name_bytes
        + coord_name_bytes
        + loop_order_name_bytes
        + cf_convention_bytes
        + num_record_components * record_component_name_bytes
        + num_metric_components * metric_component_name_bytes
        + (
            num_secondary_components * metric_component_name_bytes
            if raytracing_data_mode == "GammaUDD"
            else 0
        )
    )
    cf_convention = cast(str, par.parval_from_str("EvolvedConformalFactor_cf"))

    # Step 5.a: Validate every fixed-width string before generating the C exporter.
    _validate_fixed_width_string(header_label, header_label_bytes, "header_label")
    _validate_fixed_width_string(format_name, format_name_bytes, "format_name")
    _validate_fixed_width_string(target_basis_name, basis_name_bytes, "target_basis")
    _validate_fixed_width_string(CoordSystem, coord_name_bytes, "source_coord_system")
    _validate_fixed_width_string(loop_order_name, loop_order_name_bytes, "loop_order")
    _validate_fixed_width_string(cf_convention, cf_convention_bytes, "cf_convention")
    for component_name in record_component_names:
        _validate_fixed_width_string(
            component_name, record_component_name_bytes, "record_component_names"
        )
    for _, component_name in metric_components:
        _validate_fixed_width_string(
            component_name, metric_component_name_bytes, "metric_component_names"
        )
    if raytracing_data_mode == "GammaUDD":
        for _, component_name in secondary_components:
            _validate_fixed_width_string(
                component_name,
                metric_component_name_bytes,
                "christoffel_component_names",
            )

    # Step 6: Build the interior-point evaluation kernel for the selected
    #         output path. The all-mode path gets its own combined kernel below.
    rhs_bindings = ""
    all_rhs_bindings = ""
    if raytracing_data_mode in ("g4DD_d0", "GammaUDD"):
        rhs_bindings = """
  const REAL alpha_d0 = rhs_gfs[IDX4pt(ALPHAGF, idx3)];
  const REAL vetU_d00 = rhs_gfs[IDX4pt(VETU0GF, idx3)];
  const REAL vetU_d01 = rhs_gfs[IDX4pt(VETU1GF, idx3)];
  const REAL vetU_d02 = rhs_gfs[IDX4pt(VETU2GF, idx3)];
"""
    if raytracing_data_mode == "all":
        all_rhs_bindings = """
  const REAL alpha_d0 = rhs_gfs[IDX4pt(ALPHAGF, idx3)];
  const REAL vetU_d00 = rhs_gfs[IDX4pt(VETU0GF, idx3)];
  const REAL vetU_d01 = rhs_gfs[IDX4pt(VETU1GF, idx3)];
  const REAL vetU_d02 = rhs_gfs[IDX4pt(VETU2GF, idx3)];
"""

    loop = ""
    if raytracing_data_mode != "all":
        loop_body = rf"""
{commondata_definitions}
{params_definitions}
  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;
  MAYBE_UNUSED const REAL xx0 = xx[0][i0];
  MAYBE_UNUSED const REAL xx1 = xx[1][i1];
  MAYBE_UNUSED const REAL xx2 = xx[2][i2];
  const int idx3 = IDX3(i0, i1, i2);
  const uint64_t point_index = raytracing_data_point_index_from_logical_indices(
      i0 - NGHOSTS, i1 - NGHOSTS, i2 - NGHOSTS, payload_i0_start,
      payload_i1_start, payload_i2_start, (uint64_t)Nxx_plus_2NGHOSTS0_u32,
      (uint64_t)Nxx_plus_2NGHOSTS1_u32, (uint64_t)Nxx_plus_2NGHOSTS2_u32);
  const uint64_t base_index = point_index * (uint64_t)point_record_real_count;
{rhs_bindings}
  const REAL *restrict in_gfs = y_n_gfs;

{ccg.c_codegen(
        expr_list,
        lhs_list,
        enable_fd_codegen=True,
        enable_simd=False,
        enable_fd_functions=False,
        include_braces=False,
    )}
"""
        component_index = 3
        for _, name in metric_components:
            loop_body += (
                f"  payload_buffer[base_index + {component_index}] = (double){name};\n"
            )
            component_index += 1
        for _, name in secondary_components:
            value_name = (
                f"(use_static_christoffels ? Gamma4UDD_static{name[len('Gamma4UDD') :]} : {name})"
                if raytracing_data_mode == "GammaUDD" and enable_static_christoffels
                else name
            )
            loop_body += (
                f"  payload_buffer[base_index + {component_index}] = "
                f"(double){value_name};\n"
            )
            component_index += 1
        # Step 7: Wrap the pointwise kernel in the standard BHaH interior loop.
        loop = BHaH.simple_loop.simple_loop(
            loop_body=loop_body,
            loop_region="interior",
            enable_intrinsics=False,
            CoordSystem=CoordSystem,
            enable_rfm_precompute=enable_rfm_precompute,
            read_xxs=not enable_rfm_precompute,
            OMP_collapse=1,
        )

    all_loop = ""

    if raytracing_data_mode == "all":
        all_loop_body = rf"""
{commondata_definitions}
{params_definitions}
  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;
  MAYBE_UNUSED const REAL xx0 = xx[0][i0];
  MAYBE_UNUSED const REAL xx1 = xx[1][i1];
  MAYBE_UNUSED const REAL xx2 = xx[2][i2];
  const int idx3 = IDX3(i0, i1, i2);
  const uint64_t point_index = raytracing_data_point_index_from_logical_indices(
      i0 - NGHOSTS, i1 - NGHOSTS, i2 - NGHOSTS, payload_i0_start,
      payload_i1_start, payload_i2_start, (uint64_t)Nxx_plus_2NGHOSTS0_u32,
      (uint64_t)Nxx_plus_2NGHOSTS1_u32, (uint64_t)Nxx_plus_2NGHOSTS2_u32);
  const uint64_t base_index_g4DD =
      point_index * (uint64_t)point_record_real_count_g4DD;
  const uint64_t base_index_g4DD_d0 =
      point_index * (uint64_t)point_record_real_count_g4DD_d0;
  const uint64_t base_index_GammaUDD =
      point_index * (uint64_t)point_record_real_count_GammaUDD;
{all_rhs_bindings}
  const REAL *restrict in_gfs = y_n_gfs;

{ccg.c_codegen(
        all_codegen_expr_list,
        all_codegen_lhs_list,
        enable_fd_codegen=True,
        enable_simd=False,
        enable_fd_functions=False,
        include_braces=False,
    )}
"""
        for component_index, (_, name) in enumerate(METRIC_COMPONENTS, start=3):
            all_loop_body += (
                f"  payload_g4DD[base_index_g4DD + {component_index}] = "
                f"(double){name};\n"
            )
            all_loop_body += (
                f"  payload_g4DD_d0[base_index_g4DD_d0 + {component_index}] = "
                f"(double){name};\n"
            )
            all_loop_body += (
                f"  payload_GammaUDD[base_index_GammaUDD + {component_index}] = "
                f"(double){name};\n"
            )
        for component_index, (_, name) in enumerate(METRIC_D0_COMPONENTS, start=13):
            all_loop_body += (
                f"  payload_g4DD_d0[base_index_g4DD_d0 + {component_index}] = "
                f"(double){name};\n"
            )
        for component_index, (_, name) in enumerate(CHRISTOFFEL_COMPONENTS, start=13):
            gamma_value = (
                f"(use_static_christoffels ? "
                f"Gamma4UDD_static{name[len('Gamma4UDD') :]} : {name})"
                if enable_static_christoffels
                else name
            )
            all_loop_body += (
                f"  payload_GammaUDD[base_index_GammaUDD + {component_index}] = "
                f"(double){gamma_value};\n"
            )
        all_loop = BHaH.simple_loop.simple_loop(
            loop_body=all_loop_body,
            loop_region="interior",
            enable_intrinsics=False,
            CoordSystem=CoordSystem,
            enable_rfm_precompute=enable_rfm_precompute,
            read_xxs=not enable_rfm_precompute,
            OMP_collapse=1,
        )

    # Step 8: Register the generated C exporter and its helper functions.
    includes = [
        "<errno.h>",
        "<float.h>",
        "<inttypes.h>",
        "<math.h>",
        "<stdint.h>",
        "<stdio.h>",
        "<stdlib.h>",
        "<string.h>",
        "<sys/stat.h>",
        "<unistd.h>",
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = """
Export final Cartesian-basis raytracing data for the selected storage mode.

For a single mode, this function writes one binary file per diagnostics output
to ``./raytracing_slices/<mode>/``. In ``all`` mode it writes one file in each
of the three method directories during the same grid traversal.

The established slice magic and payload layout are retained for each method:
g4DD stores 13 doubles per point, g4DD_d0 stores 23, and GammaUDD stores 53.

where the eight-digit index is the diagnostics output index. Each file is
written to a unique temporary sibling path and then installed at the final
path without overwriting an existing file. The payload stores the physical
simulation time and one point record per full logical-grid point, including
ghost zones.

Each point record starts with Cartesian coordinates x, y, z and the 10 unique
covariant four-metric components. The selected method then appends either the
10 coordinate-time metric derivatives or the 40 four-Christoffel components.
The Cartesian metric expressions are evaluated directly from the transformed
symbolic recipes in geodesics.py. The exporter does not perform a second
runtime tensor basis transformation after importing those recipes.

The exporter evaluates Cartesian coordinates on the full logical grid and
evaluates the Cartesian metric only on interior logical-grid points.
It then fills pure outer ghost-zone tensor fields using the existing BHaH
2nd-order extrapolation ordering and copies inner-boundary ghost-zone tensor
fields from their mapped source records in bcstruct. Point records are
first assembled into a deterministic in-memory payload buffer using the
documented logical-grid ordering, then written to disk after the ghost-fill
phase completes. The serialized payload metadata exposes the logical-grid
half-open interval [-NGHOSTS, Nxx + NGHOSTS) while preserving the
underlying zero-based storage order.

@param[in] commondata  Global simulation metadata used in output naming and headers
@param[in] griddata  Host-side per-grid runtime state to export
@param time_output_index  Diagnostics output index used in the emitted filename
"""
    cfunc_type = "void"
    name = "output_raytracing_data"
    params = (
        "const commondata_struct *restrict commondata, "
        "const griddata_struct *restrict griddata, "
        "const int time_output_index"
    )
    if enable_static_christoffels:
        params += ", const int use_static_christoffels"

    prefunc = r"""
/**
 * Abort the program with a message.
 *
 * @param[in] message  Error message to print before aborting.
 */
static void raytracing_data_abort_with_message(const char *restrict message) {
  fprintf(stderr, "%s\n", message);
  exit(1);
} // END FUNCTION: raytracing_data_abort_with_message

/**
 * Cast a nonnegative integer to uint32_t, aborting on invalid input.
 *
 * @param value  Integer value to convert.
 * @param[in] label  Label for the error message.
 * @return Safely casted uint32_t value.
 */
static uint32_t raytracing_data_u32_from_nonnegative_int_or_abort(
    const int value,
    const char *restrict label
) {
  if (value < 0) {
    fprintf(stderr,
            "Error: output_raytracing_data received negative %s=%d.\n",
            label,
            value);
    exit(1);
  } // END IF: negative integers cannot be serialized as uint32_t
  if ((uint64_t)value > (uint64_t)UINT32_MAX) {
    fprintf(stderr,
            "Error: output_raytracing_data %s=%d does not fit in uint32_t.\n",
            label,
            value);
    exit(1);
  } // END IF: value exceeds the documented header format
  return (uint32_t)value;
} // END FUNCTION: raytracing_data_u32_from_nonnegative_int_or_abort

/**
 * Add two uint64_t values, aborting on overflow.
 *
 * @param a  First value.
 * @param b  Second value.
 * @param[in] label  Label for the error message.
 * @return The sum a + b.
 */
static uint64_t raytracing_data_add_u64_or_abort(
    const uint64_t a,
    const uint64_t b,
    const char *restrict label
) {
  if (UINT64_MAX - a < b) {
    fprintf(stderr,
            "Error: output_raytracing_data uint64 overflow while computing %s.\n",
            label);
    exit(1);
  } // END IF: addition would overflow
  return a + b;
} // END FUNCTION: raytracing_data_add_u64_or_abort

/**
 * Multiply two uint64_t values, aborting on overflow.
 *
 * @param a  First value.
 * @param b  Second value.
 * @param[in] label  Label for the error message.
 * @return The product a * b.
 */
static uint64_t raytracing_data_mul_u64_or_abort(
    const uint64_t a,
    const uint64_t b,
    const char *restrict label
) {
  if (a != 0ULL && b > UINT64_MAX / a) {
    fprintf(stderr,
            "Error: output_raytracing_data uint64 overflow while computing %s.\n",
            label);
    exit(1);
  } // END IF: multiplication would overflow
  return a * b;
} // END FUNCTION: raytracing_data_mul_u64_or_abort

/**
 * Cast a uint64_t value to size_t, aborting on overflow.
 *
 * @param value  Integer value to convert.
 * @param[in] label  Label for the error message.
 * @return Safely casted size_t value.
 */
static size_t raytracing_data_size_t_from_u64_or_abort(
    const uint64_t value,
    const char *restrict label
) {
  if (value > (uint64_t)SIZE_MAX) {
    fprintf(stderr,
            "Error: output_raytracing_data %s=%" PRIu64 " does not fit in size_t.\n",
            label,
            value);
    exit(1);
  } // END IF: value exceeds the local size_t range
  return (size_t)value;
} // END FUNCTION: raytracing_data_size_t_from_u64_or_abort

/**
 * Write elements to a file, aborting on failure.
 *
 * @param[in,out] fp  File pointer.
 * @param[in] data  Data to write.
 * @param element_size  Size of each element.
 * @param count  Number of elements.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_or_abort(
    FILE *restrict fp,
    const void *restrict data,
    const size_t element_size,
    const size_t count,
    const char *restrict label
) {
  if (count != 0 && element_size > SIZE_MAX / count) {
    fprintf(stderr,
            "Error: output_raytracing_data size_t overflow before writing %s.\n",
            label);
    fclose(fp);
    exit(1);
  } // END IF: fwrite byte count would overflow
  if (fwrite(data, element_size, count, fp) != count) {
    fprintf(stderr,
            "Error: output_raytracing_data failed while writing %s.\n",
            label);
    fclose(fp);
    exit(1);
  } // END IF: fwrite failed
} // END FUNCTION: raytracing_data_write_or_abort

/**
 * Write a uint32_t value in little-endian byte order.
 *
 * @param[in,out] fp  File pointer.
 * @param value  Value to write.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_u32_or_abort(
    FILE *restrict fp,
    const uint32_t value,
    const char *restrict label
) {
  uint8_t buffer[4];
  buffer[0] = (uint8_t)(value & 0xffU);
  buffer[1] = (uint8_t)((value >> 8U) & 0xffU);
  buffer[2] = (uint8_t)((value >> 16U) & 0xffU);
  buffer[3] = (uint8_t)((value >> 24U) & 0xffU);
  raytracing_data_write_or_abort(fp, buffer, sizeof(uint8_t), 4, label);
} // END FUNCTION: raytracing_data_write_u32_or_abort

/**
 * Write a uint64_t value in little-endian byte order.
 *
 * @param[in,out] fp  File pointer.
 * @param value  Value to write.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_u64_or_abort(
    FILE *restrict fp,
    const uint64_t value,
    const char *restrict label
) {
  uint8_t buffer[8];
  buffer[0] = (uint8_t)(value & 0xffULL);
  buffer[1] = (uint8_t)((value >> 8U) & 0xffULL);
  buffer[2] = (uint8_t)((value >> 16U) & 0xffULL);
  buffer[3] = (uint8_t)((value >> 24U) & 0xffULL);
  buffer[4] = (uint8_t)((value >> 32U) & 0xffULL);
  buffer[5] = (uint8_t)((value >> 40U) & 0xffULL);
  buffer[6] = (uint8_t)((value >> 48U) & 0xffULL);
  buffer[7] = (uint8_t)((value >> 56U) & 0xffULL);
  raytracing_data_write_or_abort(fp, buffer, sizeof(uint8_t), 8, label);
} // END FUNCTION: raytracing_data_write_u64_or_abort

/**
 * Write an int64_t value in little-endian byte order.
 *
 * @param[in,out] fp  File pointer.
 * @param value  Value to write.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_i64_or_abort(
    FILE *restrict fp,
    const int64_t value,
    const char *restrict label
) {
  raytracing_data_write_u64_or_abort(fp, (uint64_t)value, label);
} // END FUNCTION: raytracing_data_write_i64_or_abort

/**
 * Write a binary64 value in little-endian byte order.
 *
 * @param[in,out] fp  File pointer.
 * @param value  Value to write.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_f64_or_abort(
    FILE *restrict fp,
    const double value,
    const char *restrict label
) {
  uint64_t bits = 0ULL;
  memcpy(&bits, &value, sizeof(bits));
  raytracing_data_write_u64_or_abort(fp, bits, label);
} // END FUNCTION: raytracing_data_write_f64_or_abort

/**
 * Write an array of binary64 values in little-endian byte order.
 *
 * @param[in,out] fp  File pointer.
 * @param[in] values  Contiguous binary64 payload values.
 * @param count  Number of binary64 values to write.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_f64_array_or_abort(
    FILE *restrict fp,
    const double *restrict values,
    const size_t count,
    const char *restrict label
) {
  const uint16_t endianness_probe = 1U;
  if (((const uint8_t *restrict)&endianness_probe)[0] == 1U) {
    raytracing_data_write_or_abort(fp, values, sizeof(double), count, label);
    return;
  } // END IF: host memory already matches the documented little-endian payload format

  for (size_t i = 0; i < count; i++) {
    raytracing_data_write_f64_or_abort(fp, values[i], label);
  } // END LOOP: for i over payload binary64 values on non-little-endian hosts
} // END FUNCTION: raytracing_data_write_f64_array_or_abort

/**
 * Validate the binary64 output assumptions before serialization.
 */
static void raytracing_data_validate_binary64_output_or_abort(void) {
#if FLT_RADIX != 2 || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024 || DBL_MIN_EXP != -1021
  raytracing_data_abort_with_message(
      "Error: raytracing output requires C double to match IEEE-754 binary64 semantics.");
#endif
  if (sizeof(double) != 8)
    raytracing_data_abort_with_message(
        "Error: raytracing output requires 8-byte C double for binary64 serialization.");
  if (sizeof(REAL) > sizeof(double))
    raytracing_data_abort_with_message(
        "Error: raytracing output refuses to silently down-convert REAL values wider than binary64.");
} // END FUNCTION: raytracing_data_validate_binary64_output_or_abort

/**
 * Write a fixed-width string, null-padded to the requested size.
 *
 * @param[in,out] fp  File pointer.
 * @param[in] text  Text to write.
 * @param field_bytes  Exact field width in bytes.
 * @param[in] label  Label for the error message.
 */
static void raytracing_data_write_fixed_length_string_or_abort(
    FILE *restrict fp,
    const char *restrict text,
    const size_t field_bytes,
    const char *restrict label
) {
  char buffer[128];
  const size_t text_len = strlen(text);
  if (field_bytes > sizeof(buffer))
    raytracing_data_abort_with_message(
        "Error: raytracing header string field exceeds the internal serialization buffer.");
  if (field_bytes == 0 || text_len >= field_bytes) {
    fprintf(stderr,
            "Error: output_raytracing_data string field %s is too long for its fixed header field.\n",
            label);
    fclose(fp);
    exit(1);
  } // END IF: fixed-width string would be truncated
  memset(buffer, 0, sizeof(buffer));
  memcpy(buffer, text, text_len);
  raytracing_data_write_or_abort(fp, buffer, sizeof(char), field_bytes, label);
} // END FUNCTION: raytracing_data_write_fixed_length_string_or_abort

/**
 * Map logical-grid indices to the payload-local point-record order.
 *
 * @param i0  Signed logical x0 index in the payload window.
 * @param i1  Signed logical x1 index in the payload window.
 * @param i2  Signed logical x2 index in the payload window.
 * @param payload_i0_start  Logical x0 start index of the half-open payload window.
 * @param payload_i1_start  Logical x1 start index of the half-open payload window.
 * @param payload_i2_start  Logical x2 start index of the half-open payload window.
 * @param payload_i0_count  Full logical-grid point count in the x0 direction.
 * @param payload_i1_count  Full logical-grid point count in the x1 direction.
 * @param payload_i2_count  Full logical-grid point count in the x2 direction.
 * @return The zero-based payload-local point index in i2-major, i0-fast order.
 */
static uint64_t raytracing_data_point_index_from_logical_indices(
    const int i0,
    const int i1,
    const int i2,
    const int64_t payload_i0_start,
    const int64_t payload_i1_start,
    const int64_t payload_i2_start,
    const uint64_t payload_i0_count,
    const uint64_t payload_i1_count,
    const uint64_t payload_i2_count
) {
  const int64_t j0_signed = (int64_t)i0 - payload_i0_start;
  const int64_t j1_signed = (int64_t)i1 - payload_i1_start;
  const int64_t j2_signed = (int64_t)i2 - payload_i2_start;
  if (j0_signed < 0 || j1_signed < 0 || j2_signed < 0) {
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data logical indices fell below the payload window.");
  } // END IF: logical indices must not fall below the documented payload window
  const uint64_t j0 = (uint64_t)j0_signed;
  const uint64_t j1 = (uint64_t)j1_signed;
  const uint64_t j2 = (uint64_t)j2_signed;
  if (j0 >= payload_i0_count || j1 >= payload_i1_count ||
      j2 >= payload_i2_count) {
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data logical indices exceeded payload bounds.");
  } // END IF: logical indices must stay within the documented payload
  return j0 + payload_i0_count * (j1 + payload_i1_count * j2);
} // END FUNCTION: raytracing_data_point_index_from_logical_indices

"""

    if raytracing_data_mode == "all":
        all_header_prefix_bytes = (
            16
            + 18 * 4
            + 15 * 8
            + 15 * 8
            + header_label_bytes
            + format_name_bytes
            + basis_name_bytes
            + coord_name_bytes
            + loop_order_name_bytes
            + cf_convention_bytes
        )
        all_channel_specs = (
            (
                "g4DD",
                "NRPYRTMETRIC1",
                "NRPy RT metric data",
                "Cartesian g4DD",
                18,
                13,
                METRIC_COMPONENTS,
                (),
                False,
            ),
            (
                "g4DD_d0",
                "NRPYRTMETRIC2",
                "NRPy RT metric and dt data",
                "Cartesian g4DD+g4DDdt",
                18,
                23,
                METRIC_COMPONENTS,
                METRIC_D0_COMPONENTS,
                False,
            ),
            (
                "GammaUDD",
                "NRPYRTDATA4D",
                "NRPy RT spacetime data",
                "Cartesian g4DD+Gamma4UDD",
                19,
                53,
                METRIC_COMPONENTS,
                CHRISTOFFEL_COMPONENTS,
                True,
            ),
        )
        all_body = rf"""
  if (commondata->NUMGRIDS != 1)
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data currently requires commondata->NUMGRIDS == 1.");
  if (time_output_index < 0)
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data received a negative output index.");

  raytracing_data_validate_binary64_output_or_abort();

  const params_struct *restrict params = &griddata[0].params;
  const bc_struct *restrict bcstruct = &griddata[0].bcstruct;
  const rfm_struct *restrict rfmstruct = griddata[0].rfmstruct;
  const REAL *restrict y_n_gfs = griddata[0].gridfuncs.y_n_gfs;
  REAL *restrict auxevol_gfs = griddata[0].gridfuncs.auxevol_gfs;
  REAL *restrict xx[3] = {{
      griddata[0].xx[0],
      griddata[0].xx[1],
      griddata[0].xx[2],
  }};
  const uint32_t output_index =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          time_output_index, "time_output_index");
  const uint32_t num_grids = (uint32_t)commondata->NUMGRIDS;
  const uint32_t serialized_real_bytes = 8U;
  const uint32_t payload_includes_ghost_zones = 1U;
  const uint32_t file_is_little_endian = 1U;
  const uint32_t time_variable_is_f64 = 1U;
  const uint32_t reserved_u32 = 0U;
  const uint32_t Nxx0 =
      raytracing_data_u32_from_nonnegative_int_or_abort(params->Nxx0, "Nxx0");
  const uint32_t Nxx1 =
      raytracing_data_u32_from_nonnegative_int_or_abort(params->Nxx1, "Nxx1");
  const uint32_t Nxx2 =
      raytracing_data_u32_from_nonnegative_int_or_abort(params->Nxx2, "Nxx2");
  const uint32_t Nxx_plus_2NGHOSTS0_u32 =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS0, "Nxx_plus_2NGHOSTS0");
  const uint32_t Nxx_plus_2NGHOSTS1_u32 =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS1, "Nxx_plus_2NGHOSTS1");
  const uint32_t Nxx_plus_2NGHOSTS2_u32 =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS2, "Nxx_plus_2NGHOSTS2");
  const uint64_t point_record_count = raytracing_data_mul_u64_or_abort(
      raytracing_data_mul_u64_or_abort(
          (uint64_t)Nxx_plus_2NGHOSTS0_u32,
          (uint64_t)Nxx_plus_2NGHOSTS1_u32,
          "Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1"),
      (uint64_t)Nxx_plus_2NGHOSTS2_u32,
      "full logical-grid point count");
  const int64_t payload_i0_start = -(int64_t)NGHOSTS;
  const int64_t payload_i1_start = -(int64_t)NGHOSTS;
  const int64_t payload_i2_start = -(int64_t)NGHOSTS;
  const int64_t payload_i0_end = (int64_t)Nxx0 + (int64_t)NGHOSTS;
  const int64_t payload_i1_end = (int64_t)Nxx1 + (int64_t)NGHOSTS;
  const int64_t payload_i2_end = (int64_t)Nxx2 + (int64_t)NGHOSTS;
  SET_NXX_PLUS_2NGHOSTS_VARS(0);
  const double dxx[3] = {{
      (double)params->dxx0, (double)params->dxx1, (double)params->dxx2
  }};
  const double invdxx[3] = {{
      (double)params->invdxx0, (double)params->invdxx1, (double)params->invdxx2
  }};
  const double xxmin[3] = {{
      (double)params->xxmin0, (double)params->xxmin1, (double)params->xxmin2
  }};
  const double xxmax[3] = {{
      (double)params->xxmax0, (double)params->xxmax1, (double)params->xxmax2
  }};
  const double cart_origin[3] = {{
      (double)params->Cart_originx,
      (double)params->Cart_originy,
      (double)params->Cart_originz
  }};
"""
        for (
            mode_name,
            _,
            _,
            _,
            header_fields,
            record_count,
            metric_names,
            secondary_names,
            has_gamma_count,
        ) in all_channel_specs:
            metric_count = len(metric_names)
            secondary_count = len(secondary_names)
            header_size = (
                all_header_prefix_bytes
                + (header_fields - 18) * 4
                + record_count * record_component_name_bytes
                + metric_count * metric_component_name_bytes
                + (
                    secondary_count * metric_component_name_bytes
                    if has_gamma_count
                    else 0
                )
            )
            all_body += rf"""
  const uint32_t header_size_{mode_name} = {header_size}U;
  const uint32_t record_component_count_{mode_name} = {record_count}U;
  const uint32_t metric_component_count_{mode_name} = {metric_count}U;
  {f"const uint32_t christoffel_component_count_{mode_name} = {secondary_count}U;" if has_gamma_count else ""}
  const uint32_t point_record_real_count_{mode_name} = {record_count}U;
  const uint32_t point_record_bytes_{mode_name} =
      8U * point_record_real_count_{mode_name};
  const uint64_t point_records_bytes_{mode_name} =
      raytracing_data_mul_u64_or_abort(
          point_record_count, (uint64_t)point_record_bytes_{mode_name},
          "point-record payload size");
  const uint64_t payload_value_count_u64_{mode_name} =
      raytracing_data_mul_u64_or_abort(
          point_record_count, (uint64_t)point_record_real_count_{mode_name},
          "point-record payload value count");
  const size_t point_records_bytes_size_t_{mode_name} =
      raytracing_data_size_t_from_u64_or_abort(
          point_records_bytes_{mode_name}, "point-record payload size");
  const size_t payload_value_count_{mode_name} =
      raytracing_data_size_t_from_u64_or_abort(
          payload_value_count_u64_{mode_name}, "point-record payload value count");
  double *restrict payload_{mode_name};
  BHAH_MALLOC(payload_{mode_name}, point_records_bytes_size_t_{mode_name});
  if (payload_{mode_name} == NULL)
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data could not allocate an all-mode payload buffer.");
"""
        all_body += r"""
  const uint64_t Nxx_plus_2NGHOSTS_tot = raytracing_data_mul_u64_or_abort(
      raytracing_data_mul_u64_or_abort(
          (uint64_t)params->Nxx_plus_2NGHOSTS0,
          (uint64_t)params->Nxx_plus_2NGHOSTS1,
          "Nxx_plus_2NGHOSTS_tot"),
      (uint64_t)params->Nxx_plus_2NGHOSTS2,
      "Nxx_plus_2NGHOSTS_tot");
  const uint64_t rhs_gfs_bytes_u64 = raytracing_data_mul_u64_or_abort(
      raytracing_data_mul_u64_or_abort(
          (uint64_t)NUM_EVOL_GFS, Nxx_plus_2NGHOSTS_tot, "rhs_gfs value count"),
      (uint64_t)sizeof(REAL), "rhs_gfs byte count");
  const size_t rhs_gfs_bytes = raytracing_data_size_t_from_u64_or_abort(
      rhs_gfs_bytes_u64, "rhs_gfs byte count");
  REAL *restrict rhs_gfs;
  BHAH_MALLOC(rhs_gfs, rhs_gfs_bytes);
  if (rhs_gfs == NULL)
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data could not allocate rhs_gfs.");
  Ricci_eval(params, rfmstruct, y_n_gfs, auxevol_gfs);
  rhs_eval(commondata, params, rfmstruct, auxevol_gfs, y_n_gfs, rhs_gfs);

  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
      for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
        const uint64_t point_index = raytracing_data_point_index_from_logical_indices(
            i0 - NGHOSTS, i1 - NGHOSTS, i2 - NGHOSTS, payload_i0_start,
            payload_i1_start, payload_i2_start,
            (uint64_t)Nxx_plus_2NGHOSTS0_u32,
            (uint64_t)Nxx_plus_2NGHOSTS1_u32,
            (uint64_t)Nxx_plus_2NGHOSTS2_u32);
        const REAL xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
        REAL xCart[3];
        xx_to_Cart(params, xOrig, xCart);
"""
        for mode_name, _, _, _, _, _, _, _, _ in all_channel_specs:
            all_body += rf"""
        const uint64_t base_index_{mode_name} =
            point_index * (uint64_t)point_record_real_count_{mode_name};
        payload_{mode_name}[base_index_{mode_name} + 0] = (double)xCart[0];
        payload_{mode_name}[base_index_{mode_name} + 1] = (double)xCart[1];
        payload_{mode_name}[base_index_{mode_name} + 2] = (double)xCart[2];
"""
        all_body += r"""
      } // END LOOP: for i0 over all-mode coordinate points
    } // END LOOP: for i1 over all-mode coordinate points
  } // END LOOP: for i2 over all-mode coordinate points

"""
        all_body += all_loop
        all_body += r"""
  const bc_info_struct *restrict bc_info = &bcstruct->bc_info;
  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
    for (int dirn = 0; dirn < 3; dirn++) {
      const int num_pure_outer_boundary_points =
          bc_info->num_pure_outer_boundary_points[which_gz][dirn];
      if (num_pure_outer_boundary_points <= 0)
        continue;
      const size_t gz_idx = (size_t)dirn + 3U * (size_t)which_gz;
      const outerpt_bc_struct *restrict pure_outer_bc_array =
          bcstruct->pure_outer_bc_array[gz_idx];
#pragma omp parallel for schedule(static)
      for (int idx2d = 0; idx2d < num_pure_outer_boundary_points; idx2d++) {
        const short i0 = pure_outer_bc_array[idx2d].i0;
        const short i1 = pure_outer_bc_array[idx2d].i1;
        const short i2 = pure_outer_bc_array[idx2d].i2;
        const short face_x0 = pure_outer_bc_array[idx2d].FACEX0;
        const short face_x1 = pure_outer_bc_array[idx2d].FACEX1;
        const short face_x2 = pure_outer_bc_array[idx2d].FACEX2;
        const int idx_offset0 = IDX3(i0, i1, i2);
        const int idx_offset1 = IDX3(i0 + face_x0, i1 + face_x1, i2 + face_x2);
        const int idx_offset2 = IDX3(i0 + 2 * face_x0, i1 + 2 * face_x1, i2 + 2 * face_x2);
        const int idx_offset3 = IDX3(i0 + 3 * face_x0, i1 + 3 * face_x1, i2 + 3 * face_x2);
"""
        for mode_name, _, _, _, _, _, _, _, _ in all_channel_specs:
            all_body += rf"""
        const size_t base_offset0_{mode_name} =
            (size_t)idx_offset0 * (size_t)point_record_real_count_{mode_name};
        const size_t base_offset1_{mode_name} =
            (size_t)idx_offset1 * (size_t)point_record_real_count_{mode_name};
        const size_t base_offset2_{mode_name} =
            (size_t)idx_offset2 * (size_t)point_record_real_count_{mode_name};
        const size_t base_offset3_{mode_name} =
            (size_t)idx_offset3 * (size_t)point_record_real_count_{mode_name};
        for (size_t comp = 3; comp < (size_t)point_record_real_count_{mode_name}; comp++) {{
          payload_{mode_name}[base_offset0_{mode_name} + comp] =
              3.0 * payload_{mode_name}[base_offset1_{mode_name} + comp]
              - 3.0 * payload_{mode_name}[base_offset2_{mode_name} + comp]
              + payload_{mode_name}[base_offset3_{mode_name} + comp];
        }} // END LOOP: for comp over {mode_name} tensor payload components
"""
        all_body += r"""
      } // END LOOP: for idx2d over all-mode outer boundary points
    } // END LOOP: for dirn over all-mode outer boundary directions
  } // END LOOP: for which_gz over all-mode ghost layers

#pragma omp parallel for schedule(static)
  for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
    const innerpt_bc_struct *restrict bc = &bcstruct->inner_bc_array[pt];
"""
        for mode_name, _, _, _, _, _, _, _, _ in all_channel_specs:
            all_body += rf"""
    const size_t dst_base_offset_{mode_name} =
        (size_t)bc->dstpt * (size_t)point_record_real_count_{mode_name};
    const size_t src_base_offset_{mode_name} =
        (size_t)bc->srcpt * (size_t)point_record_real_count_{mode_name};
    for (size_t comp = 3; comp < (size_t)point_record_real_count_{mode_name}; comp++) {{
      payload_{mode_name}[dst_base_offset_{mode_name} + comp] =
          payload_{mode_name}[src_base_offset_{mode_name} + comp];
    }} // END LOOP: for comp over {mode_name} inner-boundary payload components
"""
        all_body += r"""
  } // END LOOP: for pt over all-mode inner boundary points

"""
        for (
            mode_name,
            magic,
            channel_label,
            channel_format,
            _header_fields,
            _record_count,
            metric_names,
            secondary_names,
            has_gamma_count,
        ) in all_channel_specs:
            record_names = (
                ["x", "y", "z"]
                + [name for _, name in metric_names]
                + [name for _, name in secondary_names]
            )
            all_body += _build_file_writer_block(
                mode_name=mode_name,
                magic=magic,
                header_label=channel_label,
                format_name=channel_format,
                CoordSystem=CoordSystem,
                cf_convention=cf_convention,
                record_names=record_names,
                metric_names=[name for _, name in metric_names],
                secondary_names=[name for _, name in secondary_names],
                include_christoffel_count=has_gamma_count,
            )
        all_body += "  BHAH_FREE(rhs_gfs);\n"
        body = all_body
    else:
        single_mode_magic = {
            "g4DD": "NRPYRTMETRIC1",
            "g4DD_d0": "NRPYRTMETRIC2",
            "GammaUDD": "NRPYRTDATA4D",
        }[raytracing_data_mode]
        rhs_gfs_cleanup = (
            "  BHAH_FREE(rhs_gfs);\n"
            if raytracing_data_mode in ("g4DD_d0", "GammaUDD")
            else ""
        )
        christoffel_count_declaration = (
            f"  const uint32_t christoffel_component_count = "
            f"{num_secondary_components}U;\n"
            if raytracing_data_mode == "GammaUDD"
            else ""
        )
        body = rf"""
  if (commondata->NUMGRIDS != 1)
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data currently requires commondata->NUMGRIDS == 1.");
  if (time_output_index < 0)
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data received a negative output index.");

  raytracing_data_validate_binary64_output_or_abort();

  const params_struct *restrict params = &griddata[0].params;
  const bc_struct *restrict bcstruct = &griddata[0].bcstruct;
  const rfm_struct *restrict rfmstruct = griddata[0].rfmstruct;
  const REAL *restrict y_n_gfs = griddata[0].gridfuncs.y_n_gfs;
{"  REAL *restrict auxevol_gfs = griddata[0].gridfuncs.auxevol_gfs;" if raytracing_data_mode in ("g4DD_d0", "GammaUDD") else ""}
  REAL *restrict xx[3] = {{
      griddata[0].xx[0],
      griddata[0].xx[1],
      griddata[0].xx[2],
  }};
  const uint32_t header_size = {header_size_bytes}U;
  const uint32_t output_index =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          time_output_index, "time_output_index");
  const uint32_t num_grids = (uint32_t)commondata->NUMGRIDS;
  const uint32_t serialized_real_bytes = 8U;
  const uint32_t record_component_count = {num_record_components}U;
  const uint32_t metric_component_count = {num_metric_components}U;
{christoffel_count_declaration}
  const uint32_t point_record_real_count = {point_record_real_count}U;
  const uint32_t point_record_bytes = {point_record_bytes}U;
  const uint32_t payload_includes_ghost_zones = 1U;
  const uint32_t file_is_little_endian = 1U;
  const uint32_t time_variable_is_f64 = 1U;
  const uint32_t reserved_u32 = 0U;
  const uint32_t Nxx0 =
      raytracing_data_u32_from_nonnegative_int_or_abort(params->Nxx0, "Nxx0");
  const uint32_t Nxx1 =
      raytracing_data_u32_from_nonnegative_int_or_abort(params->Nxx1, "Nxx1");
  const uint32_t Nxx2 =
      raytracing_data_u32_from_nonnegative_int_or_abort(params->Nxx2, "Nxx2");
  const uint32_t Nxx_plus_2NGHOSTS0_u32 =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS0, "Nxx_plus_2NGHOSTS0");
  const uint32_t Nxx_plus_2NGHOSTS1_u32 =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS1, "Nxx_plus_2NGHOSTS1");
  const uint32_t Nxx_plus_2NGHOSTS2_u32 =
      raytracing_data_u32_from_nonnegative_int_or_abort(
          params->Nxx_plus_2NGHOSTS2, "Nxx_plus_2NGHOSTS2");
  const uint64_t point_record_count = raytracing_data_mul_u64_or_abort(
      raytracing_data_mul_u64_or_abort(
          (uint64_t)Nxx_plus_2NGHOSTS0_u32,
          (uint64_t)Nxx_plus_2NGHOSTS1_u32,
          "Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1"),
      (uint64_t)Nxx_plus_2NGHOSTS2_u32,
      "full logical-grid point count");
  const uint64_t point_records_bytes = raytracing_data_mul_u64_or_abort(
      point_record_count,
      (uint64_t)point_record_bytes,
      "point-record payload size");
  const uint64_t payload_value_count_u64 = raytracing_data_mul_u64_or_abort(
      point_record_count,
      (uint64_t)point_record_real_count,
      "point-record payload value count");
  const uint64_t simulation_time_offset = (uint64_t)header_size;
  const uint64_t point_records_offset = raytracing_data_add_u64_or_abort(
      simulation_time_offset, 8ULL, "point-record payload offset");
  const uint64_t total_file_bytes = raytracing_data_add_u64_or_abort(
      point_records_offset, point_records_bytes, "total file bytes");
  const size_t point_records_bytes_size_t =
      raytracing_data_size_t_from_u64_or_abort(
          point_records_bytes, "point-record payload size");
  const size_t payload_value_count =
      raytracing_data_size_t_from_u64_or_abort(
          payload_value_count_u64, "point-record payload value count");
  const int64_t payload_i0_start = -(int64_t)NGHOSTS;
  const int64_t payload_i1_start = -(int64_t)NGHOSTS;
  const int64_t payload_i2_start = -(int64_t)NGHOSTS;
  const int64_t payload_i0_end = (int64_t)Nxx0 + (int64_t)NGHOSTS;
  const int64_t payload_i1_end = (int64_t)Nxx1 + (int64_t)NGHOSTS;
  const int64_t payload_i2_end = (int64_t)Nxx2 + (int64_t)NGHOSTS;
  SET_NXX_PLUS_2NGHOSTS_VARS(0);
  const double dxx[3] = {{
      (double)params->dxx0,
      (double)params->dxx1,
      (double)params->dxx2,
  }};
  const double invdxx[3] = {{
      (double)params->invdxx0,
      (double)params->invdxx1,
      (double)params->invdxx2,
  }};
  const double xxmin[3] = {{
      (double)params->xxmin0,
      (double)params->xxmin1,
      (double)params->xxmin2,
  }};
  const double xxmax[3] = {{
      (double)params->xxmax0,
      (double)params->xxmax1,
      (double)params->xxmax2,
  }};
  const double cart_origin[3] = {{
      (double)params->Cart_originx,
      (double)params->Cart_originy,
      (double)params->Cart_originz,
  }};
  const char magic[16] = "{single_mode_magic}";
  char final_filename[256];
  char temporary_filename[256];
  if (mkdir("raytracing_slices", (mode_t)0777) != 0 && errno != EEXIST)
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data could not create raytracing_slices.");
  if (mkdir("raytracing_slices/{raytracing_data_mode}", (mode_t)0777) != 0 &&
      errno != EEXIST)
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data could not create the method directory.");
  snprintf(final_filename, sizeof(final_filename),
           "./raytracing_slices/{raytracing_data_mode}/raytracing_data_t%08d.bin",
           time_output_index);
  snprintf(
      temporary_filename,
      sizeof(temporary_filename),
      "./raytracing_slices/{raytracing_data_mode}/raytracing_data_t%08d.bin.tmp.XXXXXX",
      time_output_index);

  const int temporary_fd = mkstemp(temporary_filename);
  if (temporary_fd == -1) {{
    fprintf(stderr,
            "Error: output_raytracing_data could not create temporary file %s. errno=%d\n",
            temporary_filename,
            errno);
    exit(1);
  }} // END IF: temporary output file could not be created

  const mode_t current_umask = umask((mode_t)0);
  umask(current_umask);
  if (fchmod(temporary_fd, (mode_t)(0666 & ~current_umask)) != 0) {{
    const int saved_errno = errno;
    const int close_result = close(temporary_fd);
    const int close_errno = close_result == 0 ? 0 : errno;
    const int cleanup_result = remove(temporary_filename);
    const int cleanup_errno = cleanup_result == 0 ? 0 : errno;
    if (close_result != 0 || cleanup_result != 0) {{
      fprintf(stderr,
              "Error: output_raytracing_data could not adjust permissions on %s and cleanup was incomplete. fchmod errno=%d close result=%d close errno=%d cleanup result=%d cleanup errno=%d\n",
              temporary_filename,
              saved_errno,
              close_result,
              close_errno,
              cleanup_result,
              cleanup_errno);
      exit(1);
    }} // END IF: cleanup after fchmod failure failed
    fprintf(stderr,
            "Error: output_raytracing_data could not adjust permissions on %s. errno=%d\n",
            temporary_filename,
            saved_errno);
    exit(1);
  }} // END IF: temporary output file permissions could not be adjusted

  FILE *restrict fp = fdopen(temporary_fd, "wb");
  if (fp == NULL) {{
    const int saved_errno = errno;
    const int close_result = close(temporary_fd);
    const int close_errno = close_result == 0 ? 0 : errno;
    const int cleanup_result = remove(temporary_filename);
    const int cleanup_errno = cleanup_result == 0 ? 0 : errno;
    if (close_result != 0 || cleanup_result != 0) {{
      fprintf(stderr,
              "Error: output_raytracing_data could not open %s for writing and cleanup was incomplete. fdopen errno=%d close result=%d close errno=%d cleanup result=%d cleanup errno=%d\n",
              temporary_filename,
              saved_errno,
              close_result,
              close_errno,
              cleanup_result,
              cleanup_errno);
      exit(1);
    }} // END IF: cleanup after fdopen failure failed
    fprintf(stderr,
            "Error: output_raytracing_data could not open %s for writing. errno=%d\n",
            temporary_filename,
            saved_errno);
    exit(1);
  }} // END IF: fdopen failed

{"  const uint64_t Nxx_plus_2NGHOSTS_tot = raytracing_data_mul_u64_or_abort(\n      raytracing_data_mul_u64_or_abort(\n          (uint64_t)params->Nxx_plus_2NGHOSTS0,\n          (uint64_t)params->Nxx_plus_2NGHOSTS1,\n          \"Nxx_plus_2NGHOSTS_tot\"),\n      (uint64_t)params->Nxx_plus_2NGHOSTS2,\n      \"Nxx_plus_2NGHOSTS_tot\");\n  const uint64_t rhs_gfs_bytes_u64 = raytracing_data_mul_u64_or_abort(\n      raytracing_data_mul_u64_or_abort(\n          (uint64_t)NUM_EVOL_GFS, Nxx_plus_2NGHOSTS_tot, \"rhs_gfs value count\"),\n      (uint64_t)sizeof(REAL), \"rhs_gfs byte count\");\n  const size_t rhs_gfs_bytes = raytracing_data_size_t_from_u64_or_abort(\n      rhs_gfs_bytes_u64, \"rhs_gfs byte count\");\n  REAL *restrict rhs_gfs;\n  BHAH_MALLOC(rhs_gfs, rhs_gfs_bytes);\n  if (rhs_gfs == NULL)\n    raytracing_data_abort_with_message(\n        \"Error: output_raytracing_data could not allocate rhs_gfs.\");\n  Ricci_eval(params, rfmstruct, y_n_gfs, auxevol_gfs);\n  rhs_eval(commondata, params, rfmstruct, auxevol_gfs, y_n_gfs, rhs_gfs);\n" if raytracing_data_mode in ("g4DD_d0", "GammaUDD") else ""}
  double *restrict payload_buffer;
  BHAH_MALLOC(payload_buffer, point_records_bytes_size_t);
  if (payload_buffer == NULL) {{
    fclose(fp);
    remove(temporary_filename);
    raytracing_data_abort_with_message(
        "Error: output_raytracing_data could not allocate the raytracing payload buffer.");
  }} // END IF: payload buffer allocation failed

  // Step 1: Write the fixed-size file header in documented little-endian field order.
  raytracing_data_write_or_abort(fp, magic, sizeof(char), 16, "magic");
  raytracing_data_write_u32_or_abort(fp, header_size, "header_size");
  raytracing_data_write_u32_or_abort(fp, output_index, "output_index");
  raytracing_data_write_u32_or_abort(fp, num_grids, "num_grids");
  raytracing_data_write_u32_or_abort(fp, serialized_real_bytes, "serialized_real_bytes");
  raytracing_data_write_u32_or_abort(fp, record_component_count, "record_component_count");
  raytracing_data_write_u32_or_abort(fp, metric_component_count, "metric_component_count");
{"  raytracing_data_write_u32_or_abort(\n      fp, christoffel_component_count, \"christoffel_component_count\");" if raytracing_data_mode == "GammaUDD" else ""}
  raytracing_data_write_u32_or_abort(fp, point_record_real_count, "point_record_real_count");
  raytracing_data_write_u32_or_abort(fp, point_record_bytes, "point_record_bytes");
  raytracing_data_write_u32_or_abort(
      fp, payload_includes_ghost_zones, "payload_includes_ghost_zones");
  raytracing_data_write_u32_or_abort(
      fp, file_is_little_endian, "file_is_little_endian");
  raytracing_data_write_u32_or_abort(
      fp, time_variable_is_f64, "time_variable_is_f64");
  raytracing_data_write_u32_or_abort(fp, reserved_u32, "reserved_u32");
  raytracing_data_write_u32_or_abort(fp, Nxx0, "Nxx0");
  raytracing_data_write_u32_or_abort(fp, Nxx1, "Nxx1");
  raytracing_data_write_u32_or_abort(fp, Nxx2, "Nxx2");
  raytracing_data_write_u32_or_abort(fp, Nxx_plus_2NGHOSTS0_u32, "Nxx_plus_2NGHOSTS0");
  raytracing_data_write_u32_or_abort(fp, Nxx_plus_2NGHOSTS1_u32, "Nxx_plus_2NGHOSTS1");
  raytracing_data_write_u32_or_abort(fp, Nxx_plus_2NGHOSTS2_u32, "Nxx_plus_2NGHOSTS2");
  raytracing_data_write_u64_or_abort(fp, point_record_count, "point_record_count");
  raytracing_data_write_u64_or_abort(fp, simulation_time_offset, "simulation_time_offset");
  raytracing_data_write_u64_or_abort(fp, point_records_offset, "point_records_offset");
  raytracing_data_write_u64_or_abort(fp, point_records_bytes, "point_records_bytes");
  raytracing_data_write_u64_or_abort(fp, total_file_bytes, "total_file_bytes");
  raytracing_data_write_u64_or_abort(fp, (uint64_t)NGHOSTS, "NGHOSTS");
  raytracing_data_write_u64_or_abort(
      fp, (uint64_t)Nxx_plus_2NGHOSTS0_u32, "payload_i0_count");
  raytracing_data_write_u64_or_abort(
      fp, (uint64_t)Nxx_plus_2NGHOSTS1_u32, "payload_i1_count");
  raytracing_data_write_u64_or_abort(
      fp, (uint64_t)Nxx_plus_2NGHOSTS2_u32, "payload_i2_count");
  raytracing_data_write_i64_or_abort(fp, payload_i0_start, "payload_i0_start");
  raytracing_data_write_i64_or_abort(fp, payload_i1_start, "payload_i1_start");
  raytracing_data_write_i64_or_abort(fp, payload_i2_start, "payload_i2_start");
  raytracing_data_write_i64_or_abort(fp, payload_i0_end, "payload_i0_end");
  raytracing_data_write_i64_or_abort(fp, payload_i1_end, "payload_i1_end");
  raytracing_data_write_i64_or_abort(fp, payload_i2_end, "payload_i2_end");

  for (int i = 0; i < 3; i++)
    raytracing_data_write_f64_or_abort(fp, dxx[i], "dxx");
  for (int i = 0; i < 3; i++)
    raytracing_data_write_f64_or_abort(fp, invdxx[i], "invdxx");
  for (int i = 0; i < 3; i++)
    raytracing_data_write_f64_or_abort(fp, xxmin[i], "xxmin");
  for (int i = 0; i < 3; i++)
    raytracing_data_write_f64_or_abort(fp, xxmax[i], "xxmax");
  for (int i = 0; i < 3; i++)
    raytracing_data_write_f64_or_abort(fp, cart_origin[i], "cart_origin");

"""
        body += (
            "  raytracing_data_write_fixed_length_string_or_abort("
            f'fp, "{header_label}", {header_label_bytes}, "header_label");\n'
        )
        body += (
            "  raytracing_data_write_fixed_length_string_or_abort("
            f'fp, "{format_name}", {format_name_bytes}, "format_name");\n'
        )
        body += (
            "  raytracing_data_write_fixed_length_string_or_abort("
            f'fp, "{target_basis_name}", {basis_name_bytes}, "target_basis");\n'
        )
        body += (
            f'  raytracing_data_write_fixed_length_string_or_abort(fp, "{CoordSystem}",'
            f' {coord_name_bytes}, "source_coord_system");\n'
        )
        body += (
            "  raytracing_data_write_fixed_length_string_or_abort("
            f'fp, "{loop_order_name}", {loop_order_name_bytes}, "loop_order");\n'
        )
        body += (
            "  raytracing_data_write_fixed_length_string_or_abort("
            f'fp, "{cf_convention}", {cf_convention_bytes}, "cf_convention");\n'
        )
        for component_name in record_component_names:
            body += (
                "  raytracing_data_write_fixed_length_string_or_abort("
                f'fp, "{component_name}", {record_component_name_bytes}, "record_component_names");\n'
            )
        for _, component_name in metric_components:
            body += (
                "  raytracing_data_write_fixed_length_string_or_abort("
                f'fp, "{component_name}", {metric_component_name_bytes}, "metric_component_names");\n'
            )
        if raytracing_data_mode == "GammaUDD":
            for _, component_name in secondary_components:
                body += (
                    "  raytracing_data_write_fixed_length_string_or_abort("
                    f'fp, "{component_name}", {metric_component_name_bytes}, "christoffel_component_names");\n'
                )
        body += r"""

      // Step 2: Write the physical simulation time.
      raytracing_data_write_f64_or_abort(fp, (double)commondata->time, "simulation_time");

      // Step 3: Evaluate Cartesian coordinates on the full logical grid. This
      //         loop is handwritten instead of generated by simple_loop()
      //         because it covers the full ghost-aware storage window and
      //         initializes payload coordinates before the interior tensor loop.
      for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
          for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
            const uint64_t point_index = raytracing_data_point_index_from_logical_indices(
                i0 - NGHOSTS, i1 - NGHOSTS, i2 - NGHOSTS, payload_i0_start,
                payload_i1_start, payload_i2_start, (uint64_t)Nxx_plus_2NGHOSTS0_u32,
                (uint64_t)Nxx_plus_2NGHOSTS1_u32, (uint64_t)Nxx_plus_2NGHOSTS2_u32);
            const uint64_t base_index =
                point_index * (uint64_t)point_record_real_count;
            const REAL xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
            REAL xCart[3];
            xx_to_Cart(params, xOrig, xCart);

            payload_buffer[base_index + 0] = (double)xCart[0];
            payload_buffer[base_index + 1] = (double)xCart[1];
            payload_buffer[base_index + 2] = (double)xCart[2];
          } // END LOOP: for i0 over full logical-grid x0 indices
        } // END LOOP: for i1 over full logical-grid x1 indices
      } // END LOOP: for i2 over full logical-grid x2 indices

      // Step 4: Evaluate the Cartesian metric into interior tensor payload records.
    """
        body += loop
        body += r"""

      // Step 5: Fill tensor fields for pure outer and inner ghost-zone payload records.
      //         Inner-boundary source records may map to outer-boundary records, so
      //         this phase preserves the canonical BHaH outer-then-inner ordering.
      const bc_info_struct *restrict bc_info = &bcstruct->bc_info;
      for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
        for (int dirn = 0; dirn < 3; dirn++) {
          const int num_pure_outer_boundary_points =
              bc_info->num_pure_outer_boundary_points[which_gz][dirn];
          if (num_pure_outer_boundary_points <= 0)
            continue;

          const size_t gz_idx = (size_t)dirn + 3U * (size_t)which_gz;
          const outerpt_bc_struct *restrict pure_outer_bc_array =
              bcstruct->pure_outer_bc_array[gz_idx];

    #pragma omp parallel for schedule(static)
          for (int idx2d = 0; idx2d < num_pure_outer_boundary_points; idx2d++) {
            const short i0 = pure_outer_bc_array[idx2d].i0;
            const short i1 = pure_outer_bc_array[idx2d].i1;
            const short i2 = pure_outer_bc_array[idx2d].i2;
            const short face_x0 = pure_outer_bc_array[idx2d].FACEX0;
            const short face_x1 = pure_outer_bc_array[idx2d].FACEX1;
            const short face_x2 = pure_outer_bc_array[idx2d].FACEX2;
            const int idx_offset0 = IDX3(i0, i1, i2);
            const int idx_offset1 =
                IDX3(i0 + face_x0, i1 + face_x1, i2 + face_x2);
            const int idx_offset2 =
                IDX3(i0 + 2 * face_x0, i1 + 2 * face_x1, i2 + 2 * face_x2);
            const int idx_offset3 =
                IDX3(i0 + 3 * face_x0, i1 + 3 * face_x1, i2 + 3 * face_x2);
            const size_t base_offset0 =
                (size_t)idx_offset0 * (size_t)point_record_real_count;
            const size_t base_offset1 =
                (size_t)idx_offset1 * (size_t)point_record_real_count;
            const size_t base_offset2 =
                (size_t)idx_offset2 * (size_t)point_record_real_count;
            const size_t base_offset3 =
                (size_t)idx_offset3 * (size_t)point_record_real_count;

            for (size_t comp = 3; comp < (size_t)point_record_real_count; comp++) {
              payload_buffer[base_offset0 + comp] =
                  +3.0 * payload_buffer[base_offset1 + comp]
                  - 3.0 * payload_buffer[base_offset2 + comp]
                  + 1.0 * payload_buffer[base_offset3 + comp];
            } // END LOOP: for comp over serialized tensor payload components
          } // END LOOP: for idx2d over pure outer boundary points on this face/layer
        } // END LOOP: for dirn over pure outer boundary directions
      } // END LOOP: for which_gz over pure outer ghost-zone layers

    #pragma omp parallel for schedule(static)
      for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
        const innerpt_bc_struct *restrict bc = &bcstruct->inner_bc_array[pt];
        const size_t dst_base_offset =
            (size_t)bc->dstpt * (size_t)point_record_real_count;
        const size_t src_base_offset =
            (size_t)bc->srcpt * (size_t)point_record_real_count;

        for (size_t comp = 3; comp < (size_t)point_record_real_count; comp++) {
          payload_buffer[dst_base_offset + comp] =
              payload_buffer[src_base_offset + comp];
        } // END LOOP: for comp over serialized tensor payload components
      } // END LOOP: for pt over inner boundary ghost-zone points

      // Step 6: Serialize the payload buffer using the documented little-endian binary64 layout.
      raytracing_data_write_f64_array_or_abort(
          fp, payload_buffer, payload_value_count, "point-record payload");

      // Step 7: Finalize the file and install it without overwriting an existing output.
      BHAH_FREE(payload_buffer);
"""
        body += rhs_gfs_cleanup
        body += r"""

      if (fflush(fp) != 0) {
        const int saved_errno = errno;
        fclose(fp);
        remove(temporary_filename);
        fprintf(stderr,
                "Error: output_raytracing_data failed while flushing %s. errno=%d\n",
                temporary_filename,
                saved_errno);
        exit(1);
      } // END IF: fflush failed
      if (fclose(fp) != 0) {
        const int saved_errno = errno;
        remove(temporary_filename);
        fprintf(stderr,
                "Error: output_raytracing_data failed while closing %s. errno=%d\n",
                temporary_filename,
                saved_errno);
        exit(1);
      } // END IF: fclose failed
      if (link(temporary_filename, final_filename) != 0) {
        const int saved_errno = errno;
        if (remove(temporary_filename) != 0) {
          fprintf(stderr,
                  "Error: output_raytracing_data could not install %s at %s and could not remove %s. link errno=%d cleanup errno=%d\n",
                  temporary_filename,
                  final_filename,
                  temporary_filename,
                  saved_errno,
                  errno);
          exit(1);
        } // END IF: cleanup failed after installation failure
        if (saved_errno == EEXIST)
          return;
        fprintf(stderr,
                "Error: output_raytracing_data could not install %s at %s. errno=%d\n",
                temporary_filename,
                final_filename,
                saved_errno);
        exit(1);
      } // END IF: final output file could not be installed without overwrite
      if (remove(temporary_filename) != 0) {
        fprintf(stderr,
                "Error: output_raytracing_data installed %s but could not remove %s. errno=%d\n",
                final_filename,
                temporary_filename,
                errno);
        exit(1);
      } // END IF: temporary output file could not be removed after installation
    """

    cfc.register_CFunction(
        subdirectory="diagnostics",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
