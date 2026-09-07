# nrpy/infrastructures/Dendro/Dendro_defines_h.py
"""
Emit the ``<stem>_defines.h`` header every generated CFunction source includes.

It pulls in the host header, the generated scalar/state/parameter/constant
headers and the generated CFunction declarations, and defines the upwind
selection macro NRPy owns.  It carries no field name, no finite-difference
coefficient and no numerical loop.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner
from nrpy.infrastructures.Dendro.header_guards import header_guard

BANNER = generated_file_banner()


def output_Dendro_defines_h(solver_stem: str, solver_prefix: str) -> str:
    """
    Emit ``<stem>_defines.h``, the one header every generated CFunction includes.

    :param solver_stem: Lowercase formulation stem for the generated header
        names, following Dendro's habit of naming solver files for the
        formulation (``bssnCtx.h``).
    :param solver_prefix: Bare formulation prefix, used to name the CMake
        option in the real-host diagnostic.
    :return: The complete C++ header text.

    Doctests:
    >>> header = output_Dendro_defines_h("bssn", "BSSN")
    >>> "#ifndef BSSN_DEFINES_H" in header
    True
    >>> header.rstrip().endswith("#endif  // BSSN_DEFINES_H")
    True
    >>> header.startswith("// GENERATED FILE - DO NOT EDIT")
    True
    >>> [line for line in header.splitlines() if line.startswith('#include "bssn')]
    ['#include "bssn_types.h"', '#include "bssn_constants.h"', '#include "bssn_parameters.h"', '#include "bssn_state.h"', '#include "bssn_function_prototypes.h"']
    >>> "#define UPWIND_ALG(UpwindVecU) ((UpwindVecU) > 0.0 ? 1.0 : 0.0)" in header
    True
    >>> "-DBSSN_STANDALONE_HOST=ON" in header
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
// Two things are protected from the formatter here.  The diagnostic below is
// one string literal that no column limit can break without a line
// continuation, and the include order is load-bearing: the prototypes header
// declares functions taking params_struct and the state enums, so the types,
// constants, parameters and state headers must precede it.
#if defined(NRPY_DENDRO_STANDALONE_HOST)
#include "dendro_standalone_host.h"  // NRPy-supplied host declarations
#else
#error "Not generated for a real Dendro-GR host yet; configure with -D{solver_prefix}_STANDALONE_HOST=ON, or define NRPY_DENDRO_STANDALONE_HOST when compiling a generated source by hand."
#endif
#include "{solver_stem}_types.h"
#include "{solver_stem}_constants.h"
#include "{solver_stem}_parameters.h"
#include "{solver_stem}_state.h"
#include "{solver_stem}_function_prototypes.h"
// clang-format on

// NRPy owns the upwind selection in the canonical backend, so the generated
// definition must win over any host definition -- matching
// nrpy/helpers/simd_intrinsics.h and cuda_intrinsics.h, which also #undef
// first.  A host macro with the opposite orientation would silently invert
// every advection term.
#ifdef UPWIND_ALG
#undef UPWIND_ALG
#endif
#define UPWIND_ALG(UpwindVecU) ((UpwindVecU) > 0.0 ? 1.0 : 0.0)

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
