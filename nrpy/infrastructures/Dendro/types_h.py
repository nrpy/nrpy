# nrpy/infrastructures/Dendro/types_h.py
"""
Emit the generated scalar-contract header for a Dendro solver.

The header fixes the generated scalar alias against the registered ``fp_type``.
Application-owned declarations are supplied explicitly by the assembly recipe;
the generic emitter does not infer a physics contract.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.header_guards import header_guard

BANNER = generated_file_banner()


def output_types_h(
    solver_stem: str, solver_namespace: str, additional_declarations: str
) -> str:
    """
    Emit the generated scalar-contract header.

    :param solver_stem: Lowercase formulation stem for the emitted header name.
    :param solver_namespace: Solver namespace, following Dendro's lowercase
        formulation habit (``namespace bssn``).
    :param additional_declarations: Application-owned declarations placed in
        the generated namespace before dependent prototypes.
    :return: The complete C++ header text.

    Doctests:
    >>> par.set_parval_from_str("fp_type", "double")
    >>> header = output_types_h("wave", "wave", "struct RunStatus {};")
    >>> "using DendroScalar = double;" in header
    True
    >>> "#ifndef WAVE_TYPES_H" in header
    True
    >>> header.rstrip().endswith("#endif  // WAVE_TYPES_H")
    True
    >>> "namespace wave::generated {" in header
    True
    >>> from nrpy.helpers.generic import clang_format
    >>> "struct RunStatus {};" in header
    True
    >>> "}  // END NAMESPACE: wave::generated" in clang_format(header)
    True
    """
    scalar_type = gri.DENDRO_SCALAR_TYPE
    fp_type = str(par.parval_from_str("fp_type"))
    opening, closing = header_guard(f"{solver_stem}_types.h")
    return BANNER + f"""{opening}

#include <type_traits>

#if !defined(DENDRO_SCALAR_DEFINED) && !defined(DendroScalar)
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

{additional_declarations.rstrip()}

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
