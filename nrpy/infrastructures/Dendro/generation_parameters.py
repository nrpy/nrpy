# nrpy/infrastructures/Dendro/generation_parameters.py
"""
Dendro target generation parameters.

These are ordinary NRPy parameters registered through ``par.register_param``.
They hold the *non-scientific* generation choices for lowering a frozen NRPy
environment into a Dendro-GR module.  No parallel configuration object is kept:
the values live only in the NRPy parameter registry and are read back from it.

Author: NRPy Dendro fCCZ4 infrastructure (PR 1/2/3)
"""

from types import MappingProxyType
from typing import Dict, Mapping, Tuple

import nrpy.params as par

# Register the Dendro target generation choices.  register_param() is
# idempotent, so importing this module (or calling register again) is safe.
par.register_param(
    str,
    __name__,
    "Dendro_module_name",
    "FCCZ4_GR",
    description="Name of the Dendro-GR sibling module to generate.",
)
par.register_param(
    str,
    __name__,
    "Dendro_scalar_type",
    "DendroScalar",
    description="C/C++ scalar type alias used for all generated numerical code.",
)
par.register_param(
    str,
    __name__,
    "Dendro_derivative_backend",
    "full_stencil",
    description="Derivative backend; only 'full_stencil' (direct NRPy FD) is qualified.",
)
par.register_param(
    bool,
    __name__,
    "Dendro_enable_KO",
    True,
    description="Enable NRPy-owned Kreiss-Oliger dissipation in the canonical backend.",
)
par.register_param(
    str,
    __name__,
    "Dendro_fccz4_CoordSystem",
    "Cartesian",
    description="fCCZ4 reference-metric coordinate system.",
)
par.register_param(
    str,
    __name__,
    "Dendro_fccz4_LapseEvolutionOption",
    "OnePlusLog",
    description="fCCZ4 lapse evolution option.",
)
par.register_param(
    str,
    __name__,
    "Dendro_fccz4_ShiftEvolutionOption",
    "GammaDriving2ndOrder_Covariant__Hatted",
    description="fCCZ4 shift evolution option.",
)


_ALLOWED_VALUES: Dict[str, Tuple[object, ...]] = {
    "Dendro_derivative_backend": ("full_stencil",),
    "Dendro_fccz4_CoordSystem": ("Cartesian",),
    "Dendro_fccz4_LapseEvolutionOption": ("OnePlusLog",),
    "Dendro_fccz4_ShiftEvolutionOption": ("GammaDriving2ndOrder_Covariant__Hatted",),
    # The registered `f_infinity` of `cf` is 1 regardless of the selected
    # representation, so only representations whose Minkowski value is 1 may
    # be generated: `phi` has Minkowski value 0 and would make the generated
    # initial data wrong (whitepaper section 14.1).
    "EvolvedConformalFactor_cf": ("chi", "W"),
}

# Authoritative shared NRPy parameters that define the generation profile.
# The view copies their final values from the NRPy parameter registry
# (section 9.3); they are not Dendro-owned, so they are read-only here.
_SHARED_PROFILE_PARAMS: Tuple[str, ...] = (
    "Infrastructure",
    "fp_type",
    "parallelization",
    "fd_order",
    "EvolvedConformalFactor_cf",
    "detgbarOverdetghat_equals_one",
)


def validate_generation_parameters() -> None:
    """
    Validate the Dendro generation parameters against allowed values.

    Whitepaper section 9.3 requires allowed-value validation.  Only the
    qualified profiles are accepted; anything else fails generation
    instead of silently producing an unqualified configuration.

    :raises ValueError: If any parameter holds a disallowed value.
    """
    view = dict(get_generation_parameter_view())
    scalar_type = view.get("Dendro_scalar_type")
    if not isinstance(scalar_type, str) or not scalar_type.isidentifier():
        raise ValueError(
            f"Invalid Dendro_scalar_type: {scalar_type!r} (sections 4.1/6.6)."
        )
    for name, allowed in _ALLOWED_VALUES.items():
        # Shared NRPy parameters (e.g. EvolvedConformalFactor_cf) are
        # registered by the equation modules, so a profile that builds no
        # equations simply has nothing to validate for them.
        if name not in par.glb_params_dict:
            continue
        value = par.parval_from_str(name)
        if value not in allowed:
            raise ValueError(f"Unsupported {name}={value!r}; allowed: {list(allowed)}.")
    enable_ko = par.parval_from_str("Dendro_enable_KO")
    if not isinstance(enable_ko, bool):
        raise ValueError(f"Invalid Dendro_enable_KO: {enable_ko!r}; expected bool.")


def registered_Dendro_generation_parameter_names() -> Tuple[str, ...]:
    """
    Return the Dendro and shared generation parameter names, in order.

    The shared parameters (FD order, scalar type, parallelization,
    conformal-factor representation) are registered by other NRPy modules;
    they define the profile and are included so the frozen view describes
    the complete generation profile from one place (section 9.3).

    :return: Tuple of generation parameter names (Dendro parameters first,
        then the registered shared profile parameters).
    """
    module = __name__
    dendro_names = tuple(
        name for name, param in par.glb_params_dict.items() if param.module == module
    )
    shared = tuple(
        name for name in _SHARED_PROFILE_PARAMS if name in par.glb_params_dict
    )
    return dendro_names + shared


def get_generation_parameter_view() -> Mapping[str, object]:
    """
    Return an immutable read-only view of the Dendro generation parameters.

    The values are copied from the NRPy parameter registry after registration.
    This view is for display and hashing only; it is not accepted as an
    authoring input by any CFunction builder.

    :return: Immutable mapping of Dendro generation parameter names to values.
    """
    view: Dict[str, object] = {}
    for name in registered_Dendro_generation_parameter_names():
        view[name] = par.parval_from_str(name)
    return MappingProxyType(dict(sorted(view.items())))


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
