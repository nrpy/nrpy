"""
C function registrations for staging GRHD primitive variable reconstruction routines.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot** com
"""

import textwrap
from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Sequence, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg

_INCLUDES = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
_PARAMS = "const commondata_struct *restrict commondata, const params_struct *restrict params, const ghl_parameters *restrict ghl_params, "
_PARAMS += "const int flux_dir, const ghl_eos_parameters *restrict eos, "
_PARAMS += "REAL *restrict auxevol_gfs"

_LOOP_PREAMBLE = r"""

  // Bounds are determined by the stencil, which requires a ghostzone of at least
  // NGHOSTS, but upper index includes first ghostzone point (stencil is only NGHOSTS - 1 on upper end)
  // This limit only applies to the direction of the stencil, hence the == logic below.

  const int xdir = (flux_dir == 0);
  const int ydir = (flux_dir == 1);
  const int zdir = (flux_dir == 2);

  const int imin = NGHOSTS;
  const int imax = Nxx_plus_2NGHOSTS0 - (NGHOSTS - xdir);
  const int jmin = NGHOSTS;
  const int jmax = Nxx_plus_2NGHOSTS1 - (NGHOSTS - ydir);
  const int kmin = NGHOSTS;
  const int kmax = Nxx_plus_2NGHOSTS2 - (NGHOSTS - zdir);
"""

_TEMPERATURE_POSTPROCESS = r"""
ghl_tabulated_enforce_bounds_rho_Ye_T(eos,
                                      &auxevol_gfs[IDX4pt(RHOB_RGF, index)],
                                      &auxevol_gfs[IDX4pt(YE_RGF, index)],
                                      &auxevol_gfs[IDX4pt(TEMPERATURE_RGF, index)]);
ghl_tabulated_compute_P_from_T(eos, auxevol_gfs[IDX4pt(RHOB_RGF, index)],
                               auxevol_gfs[IDX4pt(YE_RGF, index)],
                               auxevol_gfs[IDX4pt(TEMPERATURE_RGF, index)],
                               &auxevol_gfs[IDX4pt(P_RGF, index)]);

ghl_tabulated_enforce_bounds_rho_Ye_T(eos,
                                      &auxevol_gfs[IDX4pt(RHOB_LGF, index)],
                                      &auxevol_gfs[IDX4pt(YE_LGF, index)],
                                      &auxevol_gfs[IDX4pt(TEMPERATURE_LGF, index)]);
ghl_tabulated_compute_P_from_T(eos, auxevol_gfs[IDX4pt(RHOB_LGF, index)],
                               auxevol_gfs[IDX4pt(YE_LGF, index)],
                               auxevol_gfs[IDX4pt(TEMPERATURE_LGF, index)],
                               &auxevol_gfs[IDX4pt(P_LGF, index)]);
"""

_ENTROPY_POSTPROCESS = r"""
ghl_tabulated_enforce_bounds_rho_Ye_S(eos,
                                      &auxevol_gfs[IDX4pt(RHOB_RGF, index)],
                                      &auxevol_gfs[IDX4pt(YE_RGF, index)],
                                      &auxevol_gfs[IDX4pt(S_RGF, index)]);
ghl_tabulated_enforce_bounds_rho_Ye_S(eos,
                                      &auxevol_gfs[IDX4pt(RHOB_LGF, index)],
                                      &auxevol_gfs[IDX4pt(YE_LGF, index)],
                                      &auxevol_gfs[IDX4pt(S_LGF, index)]);
"""


def _scheme_settings(reconstruction_scheme: str) -> Tuple[str, int, int]:
    """
    Return reconstruction call details for a named scheme.

    :param reconstruction_scheme: Reconstruction scheme name (`"wenoz"` or `"mc"`).

    :return: Tuple of (reconstruction_call, stencil_size, stencil_shift).

    :raises ValueError: If `reconstruction_scheme` is not a supported scheme.
    """
    if reconstruction_scheme == "wenoz":
        return "ghl_wenoz_reconstruction", 6, 3
    if reconstruction_scheme == "mc":
        return "ghl_mc_reconstruction", 4, 2
    raise ValueError(f"Unknown reconstruction scheme: {reconstruction_scheme}")


def _build_reconstruction_body(
    reconstruction_scheme: str,
    V_GFs: Sequence[str],
    VR_GFs: Sequence[str],
    VL_GFs: Sequence[str],
    post_reconstruction: str = "",
) -> str:
    """
    Generate a reconstruction loop body for a given scheme and set of gridfunctions.

    :param reconstruction_scheme: Reconstruction scheme name (`"wenoz"` or `"mc"`).
    :param V_GFs: Volume-centered gridfunctions to reconstruct.
    :param VR_GFs: Right-state interface gridfunctions.
    :param VL_GFs: Left-state interface gridfunctions.
    :param post_reconstruction: Optional C code block applied after reconstruction per point.

    :return: Full C body string for the reconstruction loop.

    :raises ValueError: If provided gridfunction sequences are empty or have mismatched lengths.
    """
    if len(V_GFs) != len(VR_GFs) or len(V_GFs) != len(VL_GFs):
        raise ValueError("V_GFs, VR_GFs, and VL_GFs must have identical lengths.")
    if not V_GFs:
        raise ValueError("At least one gridfunction must be provided.")

    reconstruction_call, stencil_size, stencil_shift = _scheme_settings(
        reconstruction_scheme
    )
    stencil_upper = stencil_size - stencil_shift - 1
    num_vars = len(V_GFs)

    post_block = ""
    if post_reconstruction.strip():
        post_block = (
            textwrap.indent(textwrap.dedent(post_reconstruction).strip(), "        ")
            + "\n"
        )

    return f"""
{_LOOP_PREAMBLE}
  const int num_vars = {num_vars};
  const int V_GFs[{num_vars}] = {{{", ".join(V_GFs)}}};
  const int VR_GFs[{num_vars}] = {{{", ".join(VR_GFs)}}};
  const int VL_GFs[{num_vars}] = {{{", ".join(VL_GFs)}}};

#pragma omp parallel for
  for (int k = kmin; k < kmax; k++)
    for (int j = jmin; j < jmax; j++)
      for (int i = imin; i < imax; i++) {{
        const int index = IDX3(i, j, k);
        REAL var_data[{num_vars}][{stencil_size}];

        for (int gf = 0; gf < num_vars; gf++) {{
          for (int ind = 0; ind < {stencil_size}; ind++) {{
            // Stencil from -{stencil_shift} to +{stencil_upper} reconstructs to e.g. i-1/2
            const int stencil = IDX3(i + xdir * (ind - {stencil_shift}),
                                     j + ydir * (ind - {stencil_shift}),
                                     k + zdir * (ind - {stencil_shift}));
            var_data[gf][ind] = auxevol_gfs[IDX4pt(V_GFs[gf], stencil)];
          }}
          {reconstruction_call}(var_data[gf], &auxevol_gfs[IDX4pt(VR_GFs[gf], index)],
                                &auxevol_gfs[IDX4pt(VL_GFs[gf], index)]);
        }}
{post_block}      }}
"""


def _register_reconstruction_cfunc(
    CoordSystem: str,
    name: str,
    desc: str,
    body: str,
) -> None:
    """
    Register a reconstruction C function.

    :param CoordSystem: Coordinate system for wrapper generation.
    :param name: C function name to register.
    :param desc: Description string for function metadata.
    :param body: C function body string.
    """
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=_INCLUDES,
        desc=desc,
        cfunc_type="void",
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=_PARAMS,
        body=body,
    )


def _register_reconstruction_loops(
    CoordSystem: str,
    reconstruction_scheme: str,
    evolving_temperature: bool,
    evolving_entropy: bool,
) -> pcg.NRPyEnv_type:
    """
    Register primitive reconstruction loop and, optionally, entropy reconstruction loop.

    :param CoordSystem: Coordinate system for wrapper generation.
    :param reconstruction_scheme: Reconstruction scheme name (`"wenoz"` or `"mc"`).
    :param evolving_temperature: Whether primitives include temperature-based evolution fields.
    :param evolving_entropy: Whether to also register `reconstruction_entropy_loop`.

    :return: Updated NRPy environment.

    :raises ValueError: If `reconstruction_scheme` is not a supported scheme.
    """
    scheme_label = reconstruction_scheme.upper()
    primitive_V: Sequence[str]
    primitive_VR: Sequence[str]
    primitive_VL: Sequence[str]
    primitive_postprocess: str

    if reconstruction_scheme == "wenoz":
        if evolving_temperature:
            primitive_V = (
                "RHOBGF",
                "RESCALEDVU0GF",
                "RESCALEDVU1GF",
                "RESCALEDVU2GF",
                "YEGF",
                "TEMPERATUREGF",
            )
            primitive_VR = (
                "RHOB_RGF",
                "RESCALEDVRU0GF",
                "RESCALEDVRU1GF",
                "RESCALEDVRU2GF",
                "YE_RGF",
                "TEMPERATURE_RGF",
            )
            primitive_VL = (
                "RHOB_LGF",
                "RESCALEDVLU0GF",
                "RESCALEDVLU1GF",
                "RESCALEDVLU2GF",
                "YE_LGF",
                "TEMPERATURE_LGF",
            )
            primitive_postprocess = _TEMPERATURE_POSTPROCESS
        else:
            primitive_V = (
                "RHOBGF",
                "PGF",
                "RESCALEDVU0GF",
                "RESCALEDVU1GF",
                "RESCALEDVU2GF",
            )
            primitive_VR = (
                "RHOB_RGF",
                "P_RGF",
                "RESCALEDVRU0GF",
                "RESCALEDVRU1GF",
                "RESCALEDVRU2GF",
            )
            primitive_VL = (
                "RHOB_LGF",
                "P_LGF",
                "RESCALEDVLU0GF",
                "RESCALEDVLU1GF",
                "RESCALEDVLU2GF",
            )
            primitive_postprocess = ""
    elif reconstruction_scheme == "mc":
        if evolving_temperature:
            primitive_V = (
                "RHOBGF",
                "RESCALEDVU0GF",
                "RESCALEDVU1GF",
                "RESCALEDVU2GF",
                "YEGF",
                "TEMPERATUREGF",
            )
            primitive_VR = (
                "RHOB_RGF",
                "RESCALEDVRU0GF",
                "RESCALEDVRU1GF",
                "RESCALEDVRU2GF",
                "YE_RGF",
                "TEMPERATURE_RGF",
            )
            primitive_VL = (
                "RHOB_LGF",
                "RESCALEDVLU0GF",
                "RESCALEDVLU1GF",
                "RESCALEDVLU2GF",
                "YE_LGF",
                "TEMPERATURE_LGF",
            )
            primitive_postprocess = _TEMPERATURE_POSTPROCESS
        else:
            primitive_V = (
                "RHOBGF",
                "PGF",
                "RESCALEDVU0GF",
                "RESCALEDVU1GF",
                "RESCALEDVU2GF",
            )
            primitive_VR = (
                "RHOB_RGF",
                "P_RGF",
                "RESCALEDVRU0GF",
                "RESCALEDVRU1GF",
                "RESCALEDVRU2GF",
            )
            primitive_VL = (
                "RHOB_LGF",
                "P_LGF",
                "RESCALEDVLU0GF",
                "RESCALEDVLU1GF",
                "RESCALEDVLU2GF",
            )
            primitive_postprocess = ""
    else:
        raise ValueError(f"Unknown reconstruction scheme: {reconstruction_scheme}")

    primitive_body = _build_reconstruction_body(
        reconstruction_scheme,
        primitive_V,
        primitive_VR,
        primitive_VL,
        primitive_postprocess,
    )
    _register_reconstruction_cfunc(
        CoordSystem=CoordSystem,
        name="reconstruction_loop",
        desc=f"Reconstruct primitives using {scheme_label}",
        body=primitive_body,
    )

    if evolving_entropy:
        entropy_body = _build_reconstruction_body(
            reconstruction_scheme,
            ("SGF",),
            ("S_RGF",),
            ("S_LGF",),
            _ENTROPY_POSTPROCESS if evolving_temperature else "",
        )
        _register_reconstruction_cfunc(
            CoordSystem=CoordSystem,
            name="reconstruction_entropy_loop",
            desc=f"Reconstruct entropy using {scheme_label}",
            body=entropy_body,
        )

    return pcg.NRPyEnv()


def register_CFunction_reconstruction_loop_wenoz(
    CoordSystem: str,
    evolving_temperature: bool = False,
    evolving_entropy: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register WENOZ reconstruction loop(s).

    :param CoordSystem: The coordinate system.
    :param evolving_temperature: Whether primitives include evolving temperature variables.
    :param evolving_entropy: Whether to also register entropy reconstruction.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    return _register_reconstruction_loops(
        CoordSystem=CoordSystem,
        reconstruction_scheme="wenoz",
        evolving_temperature=evolving_temperature,
        evolving_entropy=evolving_entropy,
    )


def register_CFunction_reconstruction_loop_mc(
    CoordSystem: str,
    evolving_temperature: bool = False,
    evolving_entropy: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register MC reconstruction loop(s).

    :param CoordSystem: The coordinate system.
    :param evolving_temperature: Whether primitives include evolving temperature variables.
    :param evolving_entropy: Whether to also register entropy reconstruction.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    return _register_reconstruction_loops(
        CoordSystem=CoordSystem,
        reconstruction_scheme="mc",
        evolving_temperature=evolving_temperature,
        evolving_entropy=evolving_entropy,
    )
