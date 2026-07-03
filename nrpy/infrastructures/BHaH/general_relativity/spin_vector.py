"""
Helpers for BHaH initial-data spin-vector parameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Tuple

import nrpy.c_function as cfc
import nrpy.params as par


def register_spin_vector_CodeParameters(
    spin_alignment_vector_params: Tuple[str, str, str],
    default_chiU: Tuple[float, float, float] = (0.0, 0.0, 0.99),
) -> None:
    """
    Register public commondata spin-vector parameters.

    :param spin_alignment_vector_params: Names of the public spin-vector components.
    :param default_chiU: Default spin-vector components.
    """
    par.register_CodeParameters(
        "REAL",
        __name__,
        list(spin_alignment_vector_params),
        list(default_chiU),
        commondata=True,
    )


def register_CFunction_validate_and_set_UIUC_spin_vector(
    spin_alignment_vector_params: Tuple[str, str, str],
) -> None:
    """
    Register runtime validation and scalar-spin derivation for UIUC vector spin.

    :param spin_alignment_vector_params: Names of the public spin-vector components.
    """
    chi_x_name, chi_y_name, chi_z_name = spin_alignment_vector_params
    includes = ["BHaH_defines.h"]
    desc = r"""
Validate public UIUC spin-vector parameters and set hidden scalar chi.

The symbolic UIUC equations remain aligned with +z spin and consume
``commondata->chi``. This helper derives that scalar from public parfile
components after parsing and before initial-data setup.

@param[in,out] commondata Commondata structure containing spin components and chi.
"""
    cfunc_type = "void"
    name = "validate_and_set_UIUC_spin_vector"
    params = "commondata_struct *restrict commondata"
    body = rf"""
  const REAL chi_x = commondata->{chi_x_name};
  const REAL chi_y = commondata->{chi_y_name};
  const REAL chi_z = commondata->{chi_z_name};

  if (!isfinite((double)chi_x) || !isfinite((double)chi_y) || !isfinite((double)chi_z)) {{
    fprintf(stderr,
            "ERROR: UIUC spin vector components must be finite: {chi_x_name}=%.15e {chi_y_name}=%.15e {chi_z_name}=%.15e\n",
            (double)chi_x, (double)chi_y, (double)chi_z);
    exit(1);
  }}

  const REAL chi2 = chi_x * chi_x + chi_y * chi_y + chi_z * chi_z;
  const REAL chi_norm = sqrt(chi2);
  if (!isfinite((double)chi_norm)) {{
    fprintf(stderr,
            "ERROR: UIUC spin magnitude is not finite: {chi_x_name}=%.15e {chi_y_name}=%.15e {chi_z_name}=%.15e |chi|=%.15e\n",
            (double)chi_x, (double)chi_y, (double)chi_z, (double)chi_norm);
    exit(1);
  }}
  if (!(chi2 < 1.0)) {{
    fprintf(stderr,
            "ERROR: UIUC spin vector must satisfy |chi| < 1: {chi_x_name}=%.15e {chi_y_name}=%.15e {chi_z_name}=%.15e |chi|=%.15e\n",
            (double)chi_x, (double)chi_y, (double)chi_z, (double)chi_norm);
    exit(1);
  }}

  commondata->chi = chi_norm;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
