"""
Emit the numerical-photon axial-angular-momentum evaluator.

The generated CPU/OpenMP function evaluates ``L_z`` for one photon chunk.  It
uses the complete four-metric to lower direct four-momentum, or reconstructs
covariant momentum from the normalized variables ``u`` and ``Pi_i``.  One
diagnostic structure is filled for each active photon.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.indexedexp as ixp
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par
from nrpy.equations.general_relativity.geodesics.geodesic_diagnostics.conserved_quantities import (
    axial_angular_momentum_z_cartesian,
    photon_axial_angular_momentum_z_normalized,
)


def axial_angular_momentum(normalized_eom: bool = False) -> None:
    r"""
    Register the grouped-photon axial-angular-momentum evaluator.

    Direct evolution reads Cartesian position, contravariant four-momentum,
    and the interpolated covariant four-metric.  Normalized evolution reads
    Cartesian position, ``u``, ``Pi_1``, and ``Pi_2``; it does not require the
    metric.  Both variants fill one ``axial_angular_momentum_t`` structure per
    active photon.

    :param normalized_eom: Whether the photon state uses normalized
        coordinate-time variables instead of direct four-momentum.
    :raises ValueError: If ``normalized_eom`` is not Boolean.
    :raises ValueError: If the configured parallelization mode is not OpenMP.
    """
    if not isinstance(normalized_eom, bool):
        raise ValueError(
            f"normalized_eom must be a bool, got {type(normalized_eom).__name__}"
        )

    parallelization = par.parval_from_str("parallelization")
    if parallelization != "openmp":
        raise ValueError(
            "axial_angular_momentum currently supports only "
            "parallelization='openmp'."
        )

    # Step 1: Register the single-scalar output used by the batch and, later,
    # single-photon numerical integrators.
    angular_momentum_struct = r"""
    //==========================================
    // AXIAL ANGULAR MOMENTUM STRUCTURE
    //==========================================
    // Stores the symmetry-axis angular momentum evaluated along one trajectory.
    typedef struct {
        double Lz; // Axial angular momentum $L_z = x p_y - y p_x$.
    } axial_angular_momentum_t; // END STRUCT: axial_angular_momentum_t
    """
    Bdefines_h.register_BHaH_defines("axial_angular_momentum", angular_momentum_struct)

    # Step 2: Build the applicable symbolic expression from the equations
    # module. Both state representations use Cartesian x and y coordinates.
    x = sp.Symbol("x", real=True)
    y = sp.Symbol("y", real=True)
    if normalized_eom:
        u = sp.Symbol("u", real=True)
        PiD = ixp.declarerank1("PiD", dimension=3)
        angular_momentum_expr = photon_axial_angular_momentum_z_normalized(x, y, u, PiD)
    else:
        pU = ixp.declarerank1("pU", dimension=4)
        g4DD = ixp.declarerank2("metric_g4DD", symmetry="sym01", dimension=4)
        angular_momentum_expr = axial_angular_momentum_z_cartesian(x, y, pU, g4DD)

    math_kernel = ccg.c_codegen(
        [angular_momentum_expr],
        ["d_angular_momentum_bundle[c].Lz"],
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )
    used_symbol_names = {str(symbol) for symbol in angular_momentum_expr.free_symbols}

    # Step 3: Map symbols to the nine-component photon Structure of Arrays.
    preamble_lines = [
        "    //==========================================",
        "    // PHOTON STATE HYDRATION",
        "    //==========================================",
        "    const double x = d_f_bundle[IDX_STATE(1, c)];",
        "    const double y = d_f_bundle[IDX_STATE(2, c)];",
    ]
    if normalized_eom:
        preamble_lines.extend(
            [
                "    const double u = d_f_bundle[IDX_STATE(4, c)];",
                "    const double PiD0 = d_f_bundle[IDX_STATE(5, c)];",
                "    const double PiD1 = d_f_bundle[IDX_STATE(6, c)];",
            ]
        )
    else:
        for momentum_index in range(4):
            preamble_lines.append(
                "    const double pU"
                f"{momentum_index} = "
                f"d_f_bundle[IDX_STATE({momentum_index + 4}, c)];"
            )

        # Direct evolution stores p^mu, so load each metric component needed
        # to construct p_x and p_y. Metric storage follows upper-triangle order.
        metric_component_index = 0
        for mu in range(4):
            for nu in range(mu, 4):
                component_name = f"metric_g4DD{mu}{nu}"
                if component_name in used_symbol_names:
                    preamble_lines.append(
                        f"    const double {component_name} = "
                        "d_metric_bundle["
                        f"IDX_METRIC({metric_component_index}, c)];"
                    )
                metric_component_index += 1
    preamble = "\n".join(preamble_lines)

    # Step 4: Generate one OpenMP iteration per active photon.
    kernel_body = f"""
    #define IDX_STATE(component, ray_id) ((component) * BUNDLE_CAPACITY + (ray_id))
    #define IDX_METRIC(component, ray_id) ((component) * BUNDLE_CAPACITY + (ray_id))

    #pragma omp parallel for
    for (long int c = 0; c < current_chunk_size; ++c) {{
{preamble}

        // Evaluate $L_z$ using the selected photon state representation.
        {math_kernel}
    }} // END LOOP: for c over current_chunk_size

    #undef IDX_STATE
    #undef IDX_METRIC
    """

    arg_dict = {"d_f_bundle": "const double *restrict"}
    if not normalized_eom:
        arg_dict["d_metric_bundle"] = "const double *restrict"
    arg_dict["d_angular_momentum_bundle"] = "axial_angular_momentum_t *restrict"
    arg_dict["current_chunk_size"] = "const long int"

    prefunc, launch_body = parallel_utils.generate_kernel_and_launch_code(
        kernel_name="axial_angular_momentum_z_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization=parallelization,
        launch_dict=None,
        thread_tiling_macro_suffix="DEFAULT",
        cfunc_decorators="",
    )

    metric_description = (
        ""
        if normalized_eom
        else r"""
        @param[in] d_metric_bundle Interpolated covariant four-metric bundle
                                    $g_{\mu\nu}$."""
    )
    metric_parameter = (
        "" if normalized_eom else "const double *restrict d_metric_bundle, "
    )

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h"],
        desc=rf"""Evaluate axial angular momentum for one photon chunk.

        @param[in] d_f_bundle Photon state vectors in Structure-of-Arrays layout.{metric_description}
        @param[out] d_angular_momentum_bundle Per-photon $L_z$ structures.
        @param current_chunk_size Number of active photons in the chunk.
        @param stream_idx CPU compatibility argument.
        """,
        cfunc_type="void",
        name="axial_angular_momentum_z",
        params=(
            "const double *restrict d_f_bundle, "
            f"{metric_parameter}"
            "axial_angular_momentum_t *restrict d_angular_momentum_bundle, "
            "const long int current_chunk_size, "
            "const int stream_idx"
        ),
        include_CodeParameters_h=False,
        body=f"(void)stream_idx;\n{launch_body}",
    )

    writer_body = r"""
    if (filename == NULL || filename[0] == '\0' ||
        initial_angular_momentum == NULL || final_angular_momentum == NULL ||
        num_rays < 0) {
        fprintf(stderr, "ERROR: Invalid axial-angular-momentum output arguments.\n");
        return 1;
    }

    FILE *angular_momentum_file = fopen(filename, "wb");
    if (angular_momentum_file == NULL) {
        fprintf(
            stderr,
            "ERROR: Could not open axial-angular-momentum file '%s' for writing.\n",
            filename);
        return 1;
    }

    for (long int photon_index = 0; photon_index < num_rays; ++photon_index) {
        const uint64_t serialized_photon_index = (uint64_t)photon_index;
        const double initial_Lz = initial_angular_momentum[photon_index].Lz;
        const double final_Lz = final_angular_momentum[photon_index].Lz;
        if (fwrite(
                &serialized_photon_index,
                sizeof(serialized_photon_index),
                1,
                angular_momentum_file) != 1 ||
            fwrite(&initial_Lz, sizeof(initial_Lz), 1, angular_momentum_file) != 1 ||
            fwrite(&final_Lz, sizeof(final_Lz), 1, angular_momentum_file) != 1) {
            fprintf(
                stderr,
                "ERROR: Could not write axial angular momentum for photon %ld "
                "to '%s'.\n",
                photon_index,
                filename);
            fclose(angular_momentum_file);
            return 1;
        }
    }

    if (fclose(angular_momentum_file) != 0) {
        fprintf(
            stderr,
            "ERROR: Could not close axial-angular-momentum file '%s'.\n",
            filename);
        return 1;
    }
    return 0;
    """
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "stdint.h", "stdio.h"],
        desc=r"""Write initial and final axial angular momentum for every photon.

        Each binary record contains one ``uint64_t`` photon index followed by
        initial and final ``double`` values of $L_z$. A final value is NaN when
        terminal evaluation was unavailable.

        @param[in] filename Binary output filename.
        @param[in] initial_angular_momentum Initial $L_z$ in photon-index order.
        @param[in] final_angular_momentum Final $L_z$ in photon-index order.
        @param num_rays Number of photon records to write.
        @return Zero after successful output; one after invalid input or file failure.
        """,
        cfunc_type="int",
        name="write_axial_angular_momentum_z",
        params=(
            "const char *restrict filename, "
            "const axial_angular_momentum_t *restrict initial_angular_momentum, "
            "const axial_angular_momentum_t *restrict final_angular_momentum, "
            "const long int num_rays"
        ),
        include_CodeParameters_h=False,
        body=writer_body,
    )
