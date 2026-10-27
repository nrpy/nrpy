"""Compile and check startup helpers against a generated N8 superB project."""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> None:
    """
    Check local boundaries and distributed CFL scans without evolving.

    :raises RuntimeError: If charmc is unavailable for locating Charm++ headers.
    """
    project = Path(sys.argv[1]).resolve()
    coord = "GeneralRFM_inverse_spheroidal_fisheyeN8"
    charmc = shutil.which("charmc")
    if charmc is None:
        raise RuntimeError("charmc must be on PATH to locate Charm++ headers.")
    includes = Path(charmc).resolve().parent.parent / "include"
    sources = [
        project / coord / f"bcstruct_set_up__rfm__{coord}.cpp",
        project / coord / f"bcstruct_chare_set_up__rfm__{coord}.cpp",
        project / coord / f"ds_min_single_pt__rfm__{coord}.cpp",
        project / f"generalrfm_Cart_to_xx__{coord}.cpp",
        project / "cfl_limited_timestep.cpp",
        project / "ds_min_single_pt.cpp",
    ] + [
        project / f"inverse_fisheye_{name}__N8.cpp"
        for name in ("scale", "solve", "geometry", "component", "normalize")
    ]
    with tempfile.TemporaryDirectory(prefix="nrpy-superb-startup-") as tmp:
        executable = Path(tmp) / "check"
        subprocess.run(
            [
                os.environ.get("CXX", "mpicxx"),
                "-O2",
                f"-I{includes}",
                f"-I{project}",
                str(Path(__file__).with_suffix(".cpp")),
                *(str(source) for source in sources),
                "-o",
                str(executable),
            ],
            check=True,
            timeout=180,
        )
        subprocess.run([str(executable)], check=True, timeout=180)


if __name__ == "__main__":
    main()
