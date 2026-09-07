# nrpy/infrastructures/Dendro/types_h.py
"""
Emit the generated scalar-contract header for a Dendro solver.

The header fixes the generated scalar alias against the registered ``fp_type``
and carries the structured status record the constraint enforcement reports.
Neither names a field nor carries a physics default, so both belong to the
generic scalar contract rather than to a formulation module.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.header_guards import header_guard

BANNER = generated_file_banner()


def output_types_h(solver_stem: str, solver_namespace: str) -> str:
    """
    Emit the generated scalar-contract header.

    :param solver_stem: Lowercase formulation stem for the emitted header name.
    :param solver_namespace: Solver namespace, following Dendro's lowercase
        formulation habit (``namespace bssn``).
    :return: The complete C++ header text.

    Doctests:
    >>> par.set_parval_from_str("fp_type", "double")
    >>> header = output_types_h("bssn", "bssn")
    >>> "using DendroScalar = double;" in header
    True
    >>> "#ifndef BSSN_TYPES_H" in header
    True
    >>> header.rstrip().endswith("#endif  // BSSN_TYPES_H")
    True
    >>> "namespace bssn::generated {" in header
    True
    >>> from nrpy.helpers.generic import clang_format
    >>> "}  // END NAMESPACE: bssn::generated" in clang_format(header)
    True
    """
    scalar_type = gri.DENDRO_SCALAR_TYPE
    fp_type = str(par.parval_from_str("fp_type"))
    opening, closing = header_guard(f"{solver_stem}_types.h")
    return BANNER + f"""{opening}

#include <type_traits>

#ifndef DENDRO_SCALAR_DEFINED
#define DENDRO_SCALAR_DEFINED
using {scalar_type} = {fp_type};
#endif

namespace {solver_namespace}::generated {{

using NRPyArithmetic = {fp_type};
using TargetScalar = {scalar_type};

static_assert(sizeof(TargetScalar) == sizeof(NRPyArithmetic),
              "generated scalar contract: width mismatch");
static_assert(std::is_same_v<TargetScalar, NRPyArithmetic>,
              "generated scalar contract: alias must be the registered fp_type");

// Structured status of one det/trace enforcement pass.  The record names no
// field and carries no physics default, so it belongs to the scalar contract.
// The kernel never calls exit(): a rank-local failure is reported here and
// the host owns the global reduction.
struct detgtrazero_status_struct {{
  // Largest |det(gammabar)/det(gammahat) - 1| seen before enforcement.
  double max_abs_det_minus_one = 0.0;
  // Largest |gammabar^ij Atilde_ij| seen before enforcement.
  double max_abs_trace_residual = 0.0;
  // Points projected, points refused, and points with nonfinite diagnostics.
  unsigned long long projected_points = 0;
  unsigned long long failed_points = 0;
  unsigned long long nonfinite_points = 0;
  // Padded-block index of the first refused point, or -1 when none.  The
  // index is block-local: under the all-block entry point it locates the
  // point within its block, not within the whole local vector.
  long long first_failing_index = -1;
  // Registry position of the first nonfinite input field there, or -1.
  int first_failing_field = -1;
}};  // END STRUCT: detgtrazero_status_struct

// clang-format off
}}  // END NAMESPACE: {solver_namespace}::generated
// clang-format on

{closing}
"""


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
