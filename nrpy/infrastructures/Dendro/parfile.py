# nrpy/infrastructures/Dendro/parfile.py
"""
Emit the generated solver's sample parameter file.

The profile block carries the generated finite-difference order, padding and
dissipation switch; the parameter table itself comes from the registered
CodeParameters. The real host binds this table through TOML. The standalone
host continues to reject supplied parameter files.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List

import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.Dendro import CodeParameters
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner

BANNER = generated_file_banner("#")


def output_parfile_sample() -> str:
    """
    Emit the sample parameter-file section for the runtime parameters.

    Only parameters used by the real host's block-RHS CFunction and marked
    ``add_to_parfile`` appear; every value is the registered default.

    :return: Parameter-file text with a trailing newline.
    :raises ValueError: If a registered ``cparam_type`` has no TOML mapping.
    """
    lines: List[str] = ["[params]"]
    for cp_name in CodeParameters.runtime_parameter_names():
        code_param = par.glb_code_params_dict[cp_name]
        cparam_type = code_param.cparam_type
        default_value = code_param.defaultvalue
        base_type, _size, is_array = par.parse_cparam_type(cparam_type)
        if is_array and base_type != "char":
            # Core broadcasts a scalar default across the length, so the
            # registered default is a list; TOML renders it as an array.
            elements = default_value
            renderer = (
                (lambda element: str(int(element)))
                if base_type == "int"
                else (lambda element: repr(float(element)))
            )
            literal = "[" + ", ".join(renderer(element) for element in elements) + "]"
        elif cparam_type == "bool":
            literal = "true" if bool(default_value) else "false"
        elif cparam_type == "int":
            literal = str(int(default_value))
        elif cparam_type in ("REAL", "float", "double", gri.DENDRO_SCALAR_TYPE):
            literal = repr(float(default_value))
        else:
            raise ValueError(f"No TOML mapping for cparam_type {cparam_type!r}.")
        lines.append(f"{cp_name} = {literal}")
    return "\n".join(lines) + "\n"


def generate_default_parfile(
    solver_stem: str,
    profile_name: str,
    fd_order: int,
    ko_fd_order: int,
    ko_effective_difference_order: int,
    required_padding: int,
    enable_KreissOliger_dissipation: bool,
) -> str:
    """
    Emit the sample parameter file for one generated profile.

    :param solver_stem: Lowercase formulation stem, used as the table prefix.
    :param profile_name: Name of the generation profile.
    :param fd_order: Centered finite-difference order.
    :param ko_fd_order: Base order supplied to NRPy's ``dKOD`` construction.
    :param ko_effective_difference_order: Actual even KO difference order.
    :param required_padding: Ghost points required on every axis.
    :param enable_KreissOliger_dissipation: Whether the generated profile emits
        Kreiss-Oliger dissipation, passed by the caller that built the kernels
        rather than read from a global parameter.
    :return: The complete parameter-file text.
    :raises ValueError: If the finite-difference order, KO order, effective KO
        difference order, and required padding do not define a supported
        Dendro profile.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.infrastructures.Dendro import CFunction_roles as roles
    >>> _saved_functions = dict(cfc.CFunction_dict)
    >>> _saved_extras = dict(par.glb_extras_dict)
    >>> try:
    ...     cfc.CFunction_dict.clear()
    ...     _ = par.glb_extras_dict.pop("Dendro", None)
    ...     cfc.register_CFunction(desc="test RHS", name="fixture_rhs", body="(void)0;")
    ...     roles.set_CFunction_role("fixture_rhs", "rhs_eval_block")
    ...     text = generate_default_parfile("bssn", "vacuum", 4, 2, 4, 2, False)
    ... finally:
    ...     cfc.CFunction_dict.clear()
    ...     cfc.CFunction_dict.update(_saved_functions)
    ...     _ = par.glb_extras_dict.pop("Dendro", None)
    ...     _ = (
    ...         par.glb_extras_dict.setdefault("Dendro", _saved_extras["Dendro"])
    ...         if "Dendro" in _saved_extras
    ...         else None
    ...     )
    >>> [line for line in text.splitlines() if not line.startswith("#")][1:8]
    ['[bssn.profile]', 'name = "vacuum"', 'fd_order = 4', 'ko_fd_order = 2', 'ko_effective_difference_order = 4', 'required_padding = 2', 'ko_enabled = false']
    """
    fd_order = int(fd_order)
    ko_order = int(ko_fd_order)
    ko_effective_order = int(ko_effective_difference_order)
    if (
        fd_order not in (4, 6, 8)
        or ko_order != fd_order - 2
        or ko_effective_order != fd_order
        or int(required_padding) != fd_order // 2
    ):
        raise ValueError(
            "Dendro profile requires fd_order in (4, 6, 8), "
            "ko_fd_order=fd_order-2, effective KO order equal to fd_order, "
            "and padding equal to fd_order/2."
        )
    enable_ko = bool(enable_KreissOliger_dissipation)
    sample = output_parfile_sample()
    return BANNER + f"""#
# The standalone-host entry point takes the block count, the block extent, the
# spacing and the output and refinement selections as command-line arguments.
# The timestep is fixed at 0.5*dx and the step count at 100.  The table below
# comes from the registered CodeParameter defaults.

[{solver_stem}.profile]
name = "{profile_name}"
fd_order = {fd_order}
ko_fd_order = {ko_order}
ko_effective_difference_order = {ko_effective_order}
required_padding = {int(required_padding)}
ko_enabled = {"true" if enable_ko else "false"}

# Runtime parameters apply with -t FILE in the real Dendro host build.
# The standalone host still refuses parameter files.
""" + sample


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
