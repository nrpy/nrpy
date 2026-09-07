# nrpy/infrastructures/Dendro/cmdline_input_and_parfiles.py
"""
Emit the generated solver's sample parameter file.

The profile block carries the generated finite-difference order, padding and
dissipation switch; the parameter table itself comes from the registered
CodeParameters.  This profile has no parameter-file binding yet -- the
generated parser refuses a supplied file -- so the table is emitted commented
out.  Shipping a live table would invite edits that silently do nothing.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List

import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.Dendro.generated_file_banner import generated_file_banner

BANNER = generated_file_banner("#")


def output_parfile_sample() -> str:
    """
    Emit the sample parameter-file section for the runtime parameters.

    Only parameters with ``add_to_parfile`` appear; every value is the
    registered default.

    :return: Parameter-file text with a trailing newline.
    :raises ValueError: If a registered ``cparam_type`` has no TOML mapping.
    """
    lines: List[str] = ["[params]"]
    for cp_name, code_param in sorted(par.glb_code_params_dict.items()):
        if not code_param.add_to_parfile or code_param.cparam_type == "#define":
            continue
        cparam_type = code_param.cparam_type
        default_value = code_param.defaultvalue
        base_type, size, is_array = par.parse_cparam_type(cparam_type)
        if is_array and base_type != "char":
            # Core broadcasts a scalar default across the length, so the
            # registered default is a list; TOML renders it as an array.
            elements = (
                default_value
                if isinstance(default_value, (list, tuple))
                else [default_value] * int(size or 0)
            )
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
    required_padding: int,
    enable_KreissOliger_dissipation: bool,
) -> str:
    """
    Emit the sample parameter file for one generated profile.

    :param solver_stem: Lowercase formulation stem, used as the table prefix.
    :param profile_name: Name of the generation profile.
    :param required_padding: Ghost points required on every axis.
    :param enable_KreissOliger_dissipation: Whether the generated profile emits
        Kreiss-Oliger dissipation, passed by the caller that built the kernels
        rather than read from a global parameter.
    :return: The complete parameter-file text.

    Doctests:
    >>> import nrpy.finite_difference  # noqa: F401
    >>> par.set_parval_from_str("fd_order", 4)
    >>> text = generate_default_parfile("bssn", "vacuum", 3, False)
    >>> [line for line in text.splitlines() if not line.startswith("#")][1:6]
    ['[bssn.profile]', 'name = "vacuum"', 'fd_order = 4', 'required_padding = 3', 'ko_enabled = false']
    """
    fd_order = int(par.parval_from_str("fd_order"))
    enable_ko = bool(enable_KreissOliger_dissipation)
    sample = output_parfile_sample()
    commented = "\n".join(
        f"# {line}" if line.strip() else "#"
        for line in sample.rstrip("\n").splitlines()
    )
    return BANNER + f"""#
# The standalone-host entry point takes the block count, the block extent, the
# spacing and the output and refinement selections as command-line arguments.
# The timestep is fixed at 0.5*dx and the step count at 100.  The table below
# comes from the registered CodeParameter defaults.

[{solver_stem}.profile]
name = "{profile_name}"
fd_order = {fd_order}
required_padding = {int(required_padding)}
ko_enabled = {"true" if enable_ko else "false"}

# The parameter table below is generated from the registered CodeParameters
# and is shown for reference only: this profile has no parameter-file binding
# yet, so the solver refuses a -t argument rather than appear to apply these
# values.  The effective values are printed at startup: the generated defaults,
# except smooth_perturbation_wavelength, which is a length and so is derived
# from the block extent and spacing the host was given.
""" + commented + "\n"


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
