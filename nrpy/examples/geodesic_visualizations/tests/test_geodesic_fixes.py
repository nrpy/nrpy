"""
Small regression tests for geodesic artifact and diagnostic contracts.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import os
import struct
import tempfile
import unittest

import numpy as np
import sympy as sp

from nrpy.equations.general_relativity.geodesics.geodesic_diagnostics.conserved_quantities import (
    GeodesicDiagnostics,
)
from nrpy.examples.geodesic_visualizations import blueprint_config_and_schema as cfg
from nrpy.examples.geodesic_visualizations import blueprint_io
from nrpy.examples.geodesic_visualizations.render_lensed_image import (
    pixel_indices_from_ordinal,
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


if __name__ == "__main__":
    unittest.main()
