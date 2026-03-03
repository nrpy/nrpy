"""
Utilities for fast parameter symbol detection in psi4 code generation.
"""

import re
from typing import List

import nrpy.params as par


def fast_generalrfm_fisheye_param_symbols_from_c_code(
    CoordSystem: str,
    generated_c_code: str,
) -> List[str]:
    """
    Return non-commondata fisheye parameter symbols used in generated C code.

    This avoids expensive SymPy free-symbol traversals over large expression trees.
    """
    match = re.match(r"GeneralRFM_fisheyeN(\d+)$", CoordSystem)
    if not match:
        return []

    ntransitions = int(match.group(1))
    candidates = (
        ["fisheye_c"]
        + [f"fisheye_a{i}" for i in range(ntransitions + 1)]
        + [f"fisheye_R{i + 1}" for i in range(ntransitions)]
        + [f"fisheye_s{i + 1}" for i in range(ntransitions)]
    )
    available = [
        sym
        for sym in candidates
        if sym in par.glb_code_params_dict and not par.glb_code_params_dict[sym].commondata
    ]
    return sorted(
        [
            sym
            for sym in available
            if re.search(rf"\b{re.escape(sym)}\b", generated_c_code) is not None
        ]
    )
