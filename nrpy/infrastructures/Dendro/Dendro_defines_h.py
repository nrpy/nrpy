# nrpy/infrastructures/Dendro/Dendro_defines_h.py
"""
Emit the ``<stem>_defines.h`` header every generated CFunction source includes.

It pulls in the host header and generated scalar/state/parameter/constant
headers plus generated CFunction declarations.  It carries no field name,
finite-difference coefficient, or numerical loop.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.header_guards import header_guard

BANNER = generated_file_banner()


def output_Dendro_defines_h(solver_stem: str) -> str:
    """
    Emit ``<stem>_defines.h``, the one header every generated CFunction includes.

    :param solver_stem: Lowercase formulation stem for the generated header
        names, following Dendro's habit of naming solver files for the
        formulation (``bssnCtx.h``).
    :return: The complete C++ header text.

    Doctests:
    >>> header = output_Dendro_defines_h("bssn")
    >>> "#ifndef BSSN_DEFINES_H" in header
    True
    >>> header.rstrip().endswith("#endif  // BSSN_DEFINES_H")
    True
    >>> header.startswith("// GENERATED FILE - DO NOT EDIT")
    True
    >>> [line for line in header.splitlines() if line.startswith('#include "bssn')]
    ['#include "bssn_types.h"', '#include "bssn_constants.h"', '#include "bssn_parameters.h"', '#include "bssn_state.h"', '#include "bssn_function_prototypes.h"']
    >>> "UPWIND_ALG" in header
    False
    >>> '#include "block_geometry.h"' in header
    True
    """
    opening, closing = header_guard(f"{solver_stem}_defines.h")
    return BANNER + f"""{opening}

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string_view>

// clang-format off
// The include order is load-bearing: the prototypes header
// declares functions taking params_struct and the state enums, so the types,
// constants, parameters and state headers must precede it.
#if defined(NRPY_DENDRO_STANDALONE_HOST)
#include "dendro_standalone_host.h"  // NRPy-supplied host declarations
#else
#include "dendro.h"
#include "block_geometry.h"
#endif
#include "{solver_stem}_types.h"
#include "{solver_stem}_constants.h"
#include "{solver_stem}_parameters.h"
#include "{solver_stem}_state.h"
#include "{solver_stem}_function_prototypes.h"
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
