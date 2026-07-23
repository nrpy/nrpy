"""
Small regression tests for geodesic artifact and diagnostic contracts.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import os
import struct
import tempfile
import unittest
from contextlib import redirect_stdout
from importlib import import_module
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import numpy as np
import sympy as sp

import nrpy.c_function as cfc
from nrpy.equations.general_relativity.geodesics.geodesic_diagnostics.conserved_quantities import (
    GeodesicDiagnostics,
)
from nrpy.equations.general_relativity.geodesics.geodesics import GeodesicEquations
from nrpy.examples.geodesic_visualizations import blueprint_config_and_schema as cfg
from nrpy.examples.geodesic_visualizations import blueprint_io
from nrpy.examples.geodesic_visualizations.render_lensed_image import (
    pixel_indices_from_ordinal,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.massive.calculate_ode_rhs_massive import (
    calculate_ode_rhs_massive,
)

trajectory_visualizer = import_module(
    "nrpy.examples.geodesic_visualizations.visualize_trajectory"
)
blueprint_analysis = import_module(
    "nrpy.examples.geodesic_visualizations.blueprint_analysis"
)


class GeodesicFixTests(unittest.TestCase):
    """Exercise the narrow contracts changed by the geodesic fixes."""

    def test_blueprint_header_and_record_layout(self) -> None:
        """Require the Python artifact layout to be explicit and native-sized."""
        self.assertEqual(cfg.BLUEPRINT_HEADER_FORMAT, "=8sIIIIIIIIQ")
        self.assertEqual(cfg.BLUEPRINT_HEADER_SIZE, 48)
        self.assertEqual(cfg.BLUEPRINT_DTYPE.itemsize, 84)
        fields = cfg.BLUEPRINT_DTYPE.fields
        assert fields is not None
        self.assertEqual(fields["termination_type"][1], 0)
        self.assertEqual(fields["y_w"][1], 4)
        self.assertEqual(fields["t_f"][1], 76)

    def test_blueprint_reader_round_trip_and_terminal_statuses(self) -> None:
        """Read complete native records and reject nonterminal statuses."""
        header_values = (
            cfg.BLUEPRINT_MAGIC,
            cfg.BLUEPRINT_SCHEMA_VERSION,
            cfg.BLUEPRINT_HEADER_SIZE,
            cfg.BLUEPRINT_RECORD_SIZE,
            0,
            0,
            1,
            1,
            2,
            4,
        )
        records = np.zeros(4, dtype=cfg.BLUEPRINT_DTYPE)
        records["termination_type"] = cfg.TERM_SOURCE_PLANE
        with tempfile.TemporaryDirectory() as temp_dir:
            filename = os.path.join(temp_dir, "tile.bin")
            with open(filename, "wb") as blueprint_file:
                blueprint_file.write(
                    struct.pack(cfg.BLUEPRINT_HEADER_FORMAT, *header_values)
                )
                blueprint_file.write(records.tobytes())

            header = blueprint_io.read_blueprint_header(filename, expected_tile=(0, 0))
            chunks = list(blueprint_io.iter_blueprint_chunks(filename, 3))
            self.assertEqual(header.record_count, 4)
            self.assertEqual([start for _, start, _ in chunks], [0, 3])
            self.assertEqual(sum(len(chunk) for _, _, chunk in chunks), 4)

            records["termination_type"][0] = cfg.TERM_ACTIVE
            with open(filename, "wb") as blueprint_file:
                blueprint_file.write(
                    struct.pack(cfg.BLUEPRINT_HEADER_FORMAT, *header_values)
                )
                blueprint_file.write(records.tobytes())
            with self.assertRaises(ValueError):
                list(blueprint_io.iter_blueprint_chunks(filename, 4))

    def test_blueprint_analysis_prints_named_termination_statuses(self) -> None:
        """Require diagnostics to explain observed and configured enum values."""
        output = StringIO()
        with redirect_stdout(output):
            blueprint_analysis._print_termination_diagnostics({0: 10, 2: 5}, 15)

        diagnostics = output.getvalue()
        self.assertIn("Raw Enum  0 (TERM_COORD_RADIUS_EXCEEDED)", diagnostics)
        self.assertIn("Raw Enum  2 (TERM_EVOLUTION_MEASURE_EXCEEDED)", diagnostics)
        self.assertIn(
            "Enum  3 (TERM_RKF45_REJECTION_LIMIT): "
            "RKF45 rejected too many consecutive steps",
            diagnostics,
        )
        self.assertIn(
            "Configured evolution-measure status = 2 "
            "(TERM_EVOLUTION_MEASURE_EXCEEDED)",
            diagnostics,
        )

    def test_blueprint_analysis_uses_log10_photon_count_axis_and_named_key(
        self,
    ) -> None:
        """Require normalization histograms to use named groups and log counts."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        termination_types = np.array([0, 0, 2, 2, 2], dtype=np.int32)
        norm_abs_values = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)
        with patch.object(plt, "show"):
            blueprint_analysis.plot_norm_abs_log_histogram(
                termination_types, norm_abs_values
            )

        figure = plt.gcf()
        try:
            axis = figure.axes[0]
            self.assertEqual(axis.get_yscale(), "log")
            self.assertIn(r"\log_{10}", axis.get_ylabel())
            legend_labels = [text.get_text() for text in axis.get_legend().get_texts()]
            self.assertEqual(
                legend_labels,
                [
                    "Type 0: TERM_COORD_RADIUS_EXCEEDED",
                    "Type 2: TERM_EVOLUTION_MEASURE_EXCEEDED",
                ],
            )
        finally:
            plt.close(figure)

    def test_ordinal_mapping_is_center_preserving(self) -> None:
        """Map every dense source ray to an in-bounds output pixel."""
        pixels = [
            pixel_indices_from_ordinal(ordinal, 0, 0, 2, 1, 1, 5, 4)
            for ordinal in range(4)
        ]
        self.assertEqual(pixels, [(1, 2), (3, 2), (1, 0), (3, 0)])

    def test_carter_constant_is_finite_on_schwarzschild_axis(self) -> None:
        """The regular rotation axis must not create a removable 0/0."""
        diagnostics = GeodesicDiagnostics("KerrSchild_Cartesian", "massive")
        q_expr = diagnostics.Q_expr
        assert q_expr is not None
        substitutions = {
            diagnostics.xx[1]: 0,
            diagnostics.xx[2]: 0,
            diagnostics.xx[3]: 2,
            sp.Symbol("a_spin", real=True): 0,
            sp.Symbol("M_scale", real=True): 1,
            sp.Symbol("p0", real=True): 1,
            sp.Symbol("p1", real=True): 2,
            sp.Symbol("p2", real=True): 3,
            sp.Symbol("p3", real=True): 4,
        }
        value = q_expr.subs(substitutions).evalf()
        self.assertTrue(np.isfinite(float(value)))
        self.assertEqual(float(value), 52.0)

    def test_camera_basis_fallback_contract_is_consistent(self) -> None:
        """Require the camera basis implementations to share one fallback rule."""
        camera_sources = [
            Path(
                "nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/main_batch.py"
            ),
            Path(
                "nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/set_initial_conditions_kernel.py"
            ),
            Path(
                "nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/handle_window_plane_intersection.py"
            ),
        ]
        for source_path in camera_sources:
            source = source_path.read_text(encoding="utf-8")
            self.assertRegex(source, r"mag_[nw]_x < 1e-9")
            self.assertRegex(source, r"[nw]_z\[1\].*> 0\.999")
            self.assertRegex(
                source,
                r"double (?:alt_up|alternative_up)\[3\] = \{\{?0\.0, 1\.0, 0\.0\}\}?;",
            )

    def test_massive_rhs_does_not_require_metric_buffer(self) -> None:
        """Require the massive RHS interface to consume connections only."""
        cfc.CFunction_dict.pop("calculate_ode_rhs_massive", None)
        coordinates = list(sp.symbols("t x y z", real=True))
        velocity = list(sp.symbols("uU0 uU1 uU2 uU3", real=True))
        connection = sp.Symbol("conn_Gamma4UDD000", real=True)
        rhs = velocity + [connection * velocity[0] ** 2] * 4

        with patch(
            "nrpy.helpers.cached_functions.user_cache_dir",
            return_value="/tmp/nrpy",
        ):
            calculate_ode_rhs_massive(rhs, coordinates)
        generated = cfc.CFunction_dict["calculate_ode_rhs_massive"].full_function

        self.assertNotIn("metric_local", generated)
        self.assertIn("Gamma_local", generated)

    def test_conserved_diagnostics_expose_only_lz(self) -> None:
        """Require public symbolic diagnostics to retain Lz but not Lx/Ly."""
        for particle_type in ("massive", "photon"):
            diagnostics = GeodesicDiagnostics("KerrSchild_Cartesian", particle_type)
            self.assertFalse(hasattr(diagnostics, "L_exprs"))
            self.assertIsNotNone(diagnostics.Lz_expr)

    def test_unused_generic_christoffel_recipe_is_not_public(self) -> None:
        """Require the unused generic recipe to stay out of the public API."""
        self.assertFalse(
            hasattr(GeodesicEquations, "symbolic_christoffel_recipe_from_grid_basis")
        )

    def test_visualizer_keeps_single_row_trajectory_two_dimensional(self) -> None:
        """Require one-row trajectory files to reach the plotter as 2D data."""
        captured_data = []
        with tempfile.TemporaryDirectory() as temp_dir:
            trajectory_path = Path(temp_dir) / "trajectory.txt"
            trajectory_path.write_text("0.0 1.0 2.0 3.0 4.0\n", encoding="utf-8")

            with patch.object(
                trajectory_visualizer,
                "plot_trajectory",
                side_effect=lambda data, **_: captured_data.append(data),
            ):
                trajectory_visualizer.visualize_trajectory(str(trajectory_path))

        self.assertEqual(len(captured_data), 1)
        self.assertEqual(captured_data[0].shape, (1, 5))

    def test_visualizer_rejects_trajectory_without_spatial_columns(self) -> None:
        """Require malformed trajectory files to stop before plotting."""
        with tempfile.TemporaryDirectory() as temp_dir:
            trajectory_path = Path(temp_dir) / "trajectory.txt"
            trajectory_path.write_text("0.0 1.0 2.0 3.0\n", encoding="utf-8")

            with patch.object(trajectory_visualizer, "plot_trajectory") as plotter:
                trajectory_visualizer.visualize_trajectory(str(trajectory_path))

        plotter.assert_not_called()

    def test_photon_stage_six_computes_rhs_without_stage_update(self) -> None:
        """Require every photon integrator to keep stage 6 RHS and skip its update."""
        stage_sources = {
            Path(
                "nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_analytical.py"
            ): 2,
            Path(
                "nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py"
            ): 2,
            Path(
                "nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_analytical.py"
            ): 1,
            Path(
                "nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_numerical.py"
            ): 1,
        }
        for source_path, expected_guard_count in stage_sources.items():
            source = source_path.read_text(encoding="utf-8")
            self.assertEqual(source.count("if (stage < 6)"), expected_guard_count)

    def test_numerical_batch_initializes_result_records(self) -> None:
        """Require direct numerical-batch callers to receive deterministic records."""
        source = Path(
            "nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py"
        ).read_text(encoding="utf-8")
        self.assertIn(
            "results_buffer[i] = (blueprint_data_t){{0}};",
            source,
        )
        self.assertIn(
            "results_buffer[i].termination_type = TERMINATION_TYPE_FAILURE;",
            source,
        )


if __name__ == "__main__":
    unittest.main()
