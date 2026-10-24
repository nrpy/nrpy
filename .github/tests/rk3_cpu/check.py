"""
Compile and check complete projected RK3 CPU programs against a dense tableau.

Run from a checkout with the normal NRPy generation dependencies, GCC, and an
AVX2 CPU. Generated C, compiler output, and caches use a temporary directory
removed on success or failure. No network access or installation is performed.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot** com
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List


def run(command: List[str], cwd: Path, env: Dict[str, str], timeout: int) -> str:
    """
    Run one bounded generation, compilation, or numerical check.

    :param command: Argument vector, without shell expansion.
    :param cwd: Owned working directory.
    :param env: Environment with owned cache and temporary paths.
    :param timeout: Maximum subprocess duration in seconds.
    :raises RuntimeError: If the command fails or times out.
    :return: Bounded tail of stdout and stderr.
    """
    log = cwd / "command.log"
    try:
        with log.open("wb") as output:
            completed = subprocess.run(
                command,
                cwd=cwd,
                env=env,
                stdout=output,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"Command exceeded {timeout}s: {command[0]}") from exc
    with log.open("rb") as output:
        output.seek(max(0, log.stat().st_size - 12000))
        tail = output.read().decode("utf-8", errors="replace")
    if completed.returncode:
        raise RuntimeError(
            f"Command failed ({completed.returncode}): {command[0]}\n{tail}"
        )
    return tail


def main() -> None:
    """Check generated scalar and AVX2 programs, with one and four OpenMP threads."""
    source = Path(__file__).resolve().parent
    repository = source.parents[2]
    with tempfile.TemporaryDirectory(prefix="nrpy-rk3-cpu-") as temporary:
        root = Path(temporary)
        env = dict(os.environ)
        env.update(
            PYTHONPATH=str(repository),
            PYTHONDONTWRITEBYTECODE="1",
            XDG_CACHE_HOME=str(root / "cache"),
            TMPDIR=str(root),
            OMP_NUM_THREADS="4",
            OMP_THREAD_LIMIT="4",
            OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
        )
        print(run(["gcc", "--version"], root, env, 10).splitlines()[0], flush=True)
        print(sys.version.splitlines()[0], flush=True)
        for simd in (False, True):
            output = root / ("simd" if simd else "scalar")
            output.mkdir()
            generation = [sys.executable, str(source / "generate.py"), str(output)]
            if simd:
                generation.append("--simd")
            run(generation, output, env, 180)
            wrapped = (
                "Ricci_eval__rfm__Cartesian",
                "Ricci_eval_tile__rfm__Cartesian",
                "rhs_eval__rfm__Cartesian",
                "rhs_eval_tile__rfm__Cartesian",
                "rhs_eval_accumulate__rfm__Cartesian",
                "enforce_detgbar_equals_detghat_trAzero__rfm__Cartesian",
            )
            command = ["gcc", "-O2", "-std=gnu11", "-fopenmp", "-I", str(output)]
            if simd:
                command += ["-mavx2", "-mfma"]
            command += [str(path) for path in sorted(output.glob("*.c"))]
            command += [str(source / "runtime.c")]
            command += ["-Wl,--wrap=" + name for name in wrapped]
            command += ["-lm", "-o", str(output / "rk3_cpu")]
            run(command, output, env, 180)
            print("AVX2:" if simd else "Scalar:", flush=True)
            print(run([str(output / "rk3_cpu")], output, env, 60), end="", flush=True)


if __name__ == "__main__":
    main()
