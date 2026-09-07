# nrpy/infrastructures/Dendro/general_relativity/trusted_capture.py
"""
Shared support for the trusted-baseline sweeps in this subpackage.

Four modules capture trusted baselines in their ``__main__`` blocks -- the
generated-source ``.cpp`` files for the small emitted kernels and the trusted
expression dictionaries for the two large ones -- and all four need the same two
things: the profiles the shipped examples generate, and a way to clear one
profile out of the registries before building the next.  Keeping both here
means the four sweeps cannot drift apart in what they pin or in what they
forget to clear.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Tuple

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs
from nrpy.equations.general_relativity.fCCZ4_constraints import fCCZ4_constraints
from nrpy.equations.general_relativity.fCCZ4_RHSs import fCCZ4_RHSs

# The profiles the two examples ship, as (enable_fCCZ4, EvolvedConformalFactor_cf)
# pairs: nrpy/examples/dendro_fccz4.py generates chi and dendro_bssn.py generates
# W.  A baseline has to be the text its application emits, so the conformal
# factor is not a sweep axis -- it is a property of the profile, recorded here
# once rather than in each sweep.
SHIPPED_PROFILES: Tuple[Tuple[bool, str], ...] = ((True, "chi"), (False, "W"))

# The gauge both applications ship, spelled into the trusted-value file names as
# nrpy/infrastructures/CarpetX/general_relativity/rhs_eval.py spells its own, so
# a gauge change cannot silently reuse a trusted file.
SHIPPED_GAUGE = "OnePlusLog_GammaDriving2ndOrder_Covariant__Hatted"


def reset_generation_state() -> None:
    """
    Clear one generated profile out of the registries before building the next.

    The examples generate a single profile per process, so nothing in
    production needs this; a sweep that builds both profiles does.  The
    CFunction registry rejects a repeated name, the gridfunction registry and
    the Dendro role sidecar still describe the first profile, and -- the part
    that is easy to miss -- the equations layer's factories register the evolved
    state in their constructors, so clearing the gridfunction registry without
    clearing the memos leaves the next build with nothing registered.  The
    factories do rebuild on a changed conformal factor, so that is not the
    reason they are cleared.

    Doctests:
    >>> import nrpy.c_function as _cfc
    >>> import nrpy.grid as _gri
    >>> _ = _cfc.register_CFunction(
    ...     name="stale", desc="d", cfunc_type="void", params="", body="return;"
    ... )
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> _ = _gri.register_gridfunctions("staleGF", group="EVOL")
    >>> par.glb_extras_dict.setdefault("Dendro", {})["required_padding"] = 3
    >>> reset_generation_state()
    >>> "stale" in _cfc.CFunction_dict, "staleGF" in _gri.glb_gridfcs_dict
    (False, False)
    >>> par.glb_extras_dict.get("Dendro")
    """
    cfc.CFunction_dict.clear()
    gri.glb_gridfcs_dict.clear()
    par.glb_extras_dict.pop("Dendro", None)
    for factory in (
        BSSN_quantities,
        BSSN_RHSs,
        BSSN_constraints,
        fCCZ4_RHSs,
        fCCZ4_constraints,
    ):
        factory.clear()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
