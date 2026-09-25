r"""
Generate, build, run, and check the complete NRPy Dendro applications.

For one formulation (BSSN or fCCZ4), this helper generates the W and chi
applications twice each, builds them standalone against the Dendrolib commit
pinned by the generator, runs short MPI TwoPunctures evolutions, and checks
the output against exact identities, values computed independently of the
evolution, and comparisons between runs of the same binary. Each check prints
its identifier, the validation layer it proves, the measured value, and the
tolerance. The process exits nonzero if any check fails.

Usage (from the repository root):
    python nrpy/examples/tests/dendro_application_check.py \
        --formulation {bssn,fccz4} --work-dir DIR \
        [--launcher "mpiexec --oversubscribe --bind-to none"] [--ranks 4] \
        [--build-jobs 4]

The helper refuses rank counts that would start a 3-rank run, whose
horizon-finder checkpoint Dendrolib misroutes, and horizon finding with
BSSN_IO_OUTPUT_FREQ = 0, which aborts in Dendrolib.

Run without arguments, the helper runs its doctests. Configuring each
generated project downloads the Dendrolib and toml11 revisions pinned by the
generated CMakeLists.txt, so the run needs network access.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import argparse
import importlib.metadata
import json
import math
import os
import shlex
import shutil
import signal
import struct
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, cast

REPO_ROOT = Path(__file__).resolve().parents[3]

TIMEOUT_GENERATE = 600
TIMEOUT_CONFIGURE = 900
TIMEOUT_BUILD = 1800
TIMEOUT_RUN = 900
TIMEOUT_NEGATIVE = 120
TIMEOUT_VERSION = 60
LOG_TAIL_LINES = 60
NODE_CEILING = 2_500_000

SOLVER = {
    "bssn": ("NRPy_BSSN_GR", "nrpyBssnSolver", "bssn"),
    "fccz4": ("NRPy_fCCZ4_GR", "nrpyFccz4Solver", "fccz4"),
}

COMMON_OVERRIDES: Dict[str, str] = {
    "BSSN_MINDEPTH": "3",
    "BSSN_BH1_AMR_R": "1.0",
    "BSSN_BH2_AMR_R": "1.0",
    "BSSN_AMR_R_RATIO": "1.618033988749895",
    "BSSN_DENDRO_GRAIN_SZ": "50",
    "BSSN_WAVELET_TOL": "0.001",
    "BSSN_GW_EXTRACT_FREQ": "4",
    "BSSN_GW_RADAII": "[50.0, 100.0]",
    "BSSN_GW_NUM_RADAII": "2",
    "BSSN_GW_L_MODES": "[2, 3, 4]",
    "BSSN_GW_NUM_LMODES": "3",
    "TPID_NPOINTS_A": "24",
    "TPID_NPOINTS_B": "24",
    "TPID_NPOINTS_PHI": "6",
    "TPID_FILEPREFIX": '"tp_ci"',
    "BSSN_PROFILE_FILE_PREFIX": '"dat/dgr"',
    "BSSN_CHKPT_FILE_PREFIX": '"cp/ci_cp"',
    "BSSN_VTU_FILE_PREFIX": '"vtu/ci"',
}

# Profile P: FD6, horizons, a forced remesh at step 4, and checkpoints.
PROFILE_P: Dict[str, str] = {
    "BSSN_ELE_ORDER": "6",
    "BSSN_MAXDEPTH": "13",
    "BSSN_INIT_GRID_ITER": "1",
    "BSSN_MAX_ITERATIONS": "8",
    "BSSN_TIME_STEP_OUTPUT_FREQ": "4",
    "BSSN_REMESH_TEST_FREQ": "4",
    "BSSN_IO_OUTPUT_FREQ": "4",
    "BSSN_CHECKPT_FREQ": "4",
    "AEH_SOLVER_FREQ": "4",
}

# Profile O: one converged octree for every finite-difference order.
PROFILE_O: Dict[str, str] = {
    "BSSN_MAXDEPTH": "10",
    "BSSN_INIT_GRID_ITER": "10",
    "BSSN_MAX_ITERATIONS": "4",
    "BSSN_TIME_STEP_OUTPUT_FREQ": "1",
    "BSSN_REMESH_TEST_FREQ": "0",
    "AEH_SOLVER_FREQ": "0",
    "BSSN_IO_OUTPUT_FREQ": "0",
    "BSSN_CHECKPT_FREQ": "0",
}

TPID_FILE = "tp_ci_nrpy_tpid_sol.bin"
TPID_TAG = b"NRPy-TPID-1"
TPID_INPUT_COUNT = 33
TPID_RESULT_COUNT = 8
TPID_INDEX_P_PLUS_Y = 15
TPID_INDEX_P_MINUS_Y = 18


class CheckError(RuntimeError):
    """Raise when a step cannot produce the output a check needs."""


class Report:
    """Record check results and print them as they are made."""

    def __init__(self) -> None:
        """Start with no results."""
        self.failures: List[str] = []
        self.count = 0

    def check(
        self, ident: str, layer: str, what: str, measured: str, bound: str, ok: bool
    ) -> None:
        """
        Record and print one check.

        :param ident: Check identifier, for example ``R1``.
        :param layer: Validation layer the check proves.
        :param what: Short description of the checked property.
        :param measured: Measured value, formatted.
        :param bound: Tolerance or criterion, formatted.
        :param ok: Whether the check passed.
        """
        self.count += 1
        status = "PASS" if ok else "FAIL"
        print(
            f"[{status}] {ident} ({layer}) {what}: measured {measured}; criterion {bound}",
            flush=True,
        )
        if not ok:
            self.failures.append(f"{ident}: {what}")


def tail(path: Path, lines: int = LOG_TAIL_LINES) -> str:
    r"""
    Return at most the last ``lines`` lines of a text file.

    :param path: File to read.
    :param lines: Maximum number of lines.
    :return: Tail text, or an empty string if the file is missing.

    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as d:
    ...     p = Path(d) / "log"
    ...     _ = p.write_text("a\nb\nc\n")
    ...     tail(p, 2)
    'b\nc'
    """
    if not path.exists():
        return ""
    return "\n".join(path.read_text(errors="replace").splitlines()[-lines:])


def run_logged(
    argv: Sequence[str],
    cwd: Path,
    log: Path,
    timeout: int,
    env: Optional[Dict[str, str]] = None,
) -> Tuple[int, bool]:
    """
    Run a command with its output in a log file.

    :param argv: Argument vector; no shell is used.
    :param cwd: Working directory.
    :param log: Log file receiving standard output and standard error.
    :param timeout: Time limit in seconds.
    :param env: Environment, or None to inherit.
    :return: Exit status and whether the time limit was reached. On timeout, or
        if the helper is interrupted, the command's process group is killed;
        MPI ranks, which the launcher places in their own process groups, then
        exit through the launcher's failure detection.
    :raises BaseException: Any interruption other than the timeout, such as
        KeyboardInterrupt, re-raised after the process group is killed.
    """
    with open(log, "w", encoding="utf-8") as stream:
        with subprocess.Popen(
            list(argv),
            cwd=cwd,
            stdout=stream,
            stderr=subprocess.STDOUT,
            env=env,
            shell=False,
            start_new_session=True,
        ) as process:
            try:
                return process.wait(timeout=timeout), False
            except BaseException as error:
                # Kill the launcher's process group (including build jobs) on a
                # timeout or any interruption; MPI ranks exit once it is gone.
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
                if isinstance(error, subprocess.TimeoutExpired):
                    return 124, True
                raise


def run_checked(
    argv: Sequence[str],
    cwd: Path,
    log: Path,
    timeout: int,
    env: Optional[Dict[str, str]] = None,
) -> None:
    """
    Run a command that must succeed; on failure print the log tail and raise.

    :param argv: Argument vector; no shell is used.
    :param cwd: Working directory.
    :param log: Log file receiving standard output and standard error.
    :param timeout: Time limit in seconds.
    :param env: Environment, or None to inherit.
    :raises CheckError: If the command fails or reaches its time limit.
    """
    status, timed_out = run_logged(argv, cwd, log, timeout, env)
    if timed_out or status != 0:
        print(tail(log), flush=True)
        reason = f"timed out after {timeout} s" if timed_out else f"exit {status}"
        raise CheckError(f"{' '.join(argv)} in {cwd}: {reason}")


def apply_overrides(base_text: str, overrides: Dict[str, str]) -> str:
    r"""
    Replace or append top-level TOML keys, before the first table header.

    :param base_text: Text of the generated parameter file.
    :param overrides: Keys and literal TOML values.
    :return: Updated parameter-file text.
    :raises CheckError: If an override key appears in a non-flat form.

    >>> apply_overrides('A = 1\nB = 2\n[T]\nC = 3\n', {"B": "5", "D": "6"})
    'A = 1\nB = 5\nD = 6\n[T]\nC = 3\n'
    >>> apply_overrides('A = 1\n', {"A": "2", "E": '"x"'})
    'A = 2\nE = "x"\n'
    >>> apply_overrides('[A]\nx = 1\n', {"A": "2"})
    Traceback (most recent call last):
    ...
    CheckError: override key A appears as a table in the parameter file
    >>> apply_overrides('A.b = 1\n', {"A": "2"})
    Traceback (most recent call last):
    ...
    CheckError: override key A appears in dotted form
    """
    lines = base_text.splitlines()
    result: List[str] = []
    seen = set()
    appended = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            table = stripped.strip("[]").strip()
            if table.split(".")[0].strip() in overrides:
                raise CheckError(
                    f"override key {table} appears as a table in the parameter file"
                )
            if not appended:
                result.extend(
                    f"{k} = {v}" for k, v in overrides.items() if k not in seen
                )
                appended = True
            result.append(line)
            continue
        key, equals, _ = stripped.partition("=")
        key = key.strip()
        if equals and not appended and key in overrides:
            result.append(f"{key} = {overrides[key]}")
            seen.add(key)
            continue
        root = key.split(".")[0].strip()
        if equals and "." in key and root in overrides:
            raise CheckError(f"override key {root} appears in dotted form")
        result.append(line)
    if not appended:
        result.extend(f"{k} = {v}" for k, v in overrides.items() if k not in seen)
    return "\n".join(result) + "\n"


def parse_table(path: Path) -> Tuple[List[str], List[List[float]]]:
    r"""
    Parse a whitespace-separated numeric table with an optional header line.

    :param path: File to parse.
    :return: Header names (empty if none) and numeric rows.

    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as d:
    ...     p = Path(d) / "t.dat"
    ...     _ = p.write_text("TimeStep\t time\t H\t\n0\t0\t1e-6\n4\t0.5\tnan\n")
    ...     h, rows = parse_table(p)
    ...     h, rows[0], math.isnan(rows[1][2])
    (['TimeStep', 'time', 'H'], [0.0, 0.0, 1e-06], True)
    """
    header: List[str] = []
    rows: List[List[float]] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        try:
            rows.append([float(x) for x in fields])
        except ValueError:
            if not header:
                header = fields
    return header, rows


def parse_modes(path: Path) -> List[Tuple[int, List[complex]]]:
    r"""
    Parse a Dendro-GR-layout Psi4 mode file into (step, per-radius values).

    :param path: ``dgr_GW_l<l>_m<m>.dat`` file.
    :return: One entry per row: the step and the complex value at each radius.

    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as d:
    ...     p = Path(d) / "m.dat"
    ...     _ = p.write_text("TimeStep\t t\tr0\t\n4\t1e-1\t(1e-7,-2e-7)\t\n")
    ...     parse_modes(p)
    [(4, [(1e-07-2e-07j)])]
    """
    result = []
    for line in path.read_text().splitlines()[1:]:
        fields = [f for f in line.split("\t") if f.strip()]
        if not fields:
            continue
        values = []
        for field in fields[2:]:
            real, imag = field.strip().strip("()").split(",")
            values.append(complex(float(real), float(imag)))
        result.append((int(float(fields[0])), values))
    return result


def relative(a: float, b: float) -> float:
    """
    Return |a - b| / max(|a|, |b|), and 0 when both vanish.

    :param a: First value.
    :param b: Second value.
    :return: Relative difference.

    >>> relative(1.0, 1.001) < 1.0e-3
    True
    >>> relative(0.0, 0.0)
    0.0
    """
    scale = max(abs(a), abs(b))
    return 0.0 if scale == 0.0 else abs(a - b) / scale


def column(header: List[str], rows: List[List[float]], name: str) -> List[float]:
    """
    Return one named column of a parsed table.

    :param header: Header names.
    :param rows: Numeric rows.
    :param name: Column name.
    :return: Column values.
    :raises CheckError: If the column is missing.
    """
    if name not in header:
        raise CheckError(f"column {name} missing from header {header}")
    index = header.index(name)
    return [row[index] for row in rows]


def row_at(rows: List[List[float]], step: int) -> List[float]:
    """
    Return the row whose first column equals ``step``.

    :param rows: Numeric rows.
    :param step: Step number.
    :return: Matching row.
    :raises CheckError: If no row matches.
    """
    for row in rows:
        if int(round(row[0])) == step:
            return row
    raise CheckError(f"no row for step {step}")


def trees_identical(first: Path, second: Path) -> Tuple[bool, str]:
    """
    Compare two directory trees file by file.

    :param first: First tree.
    :param second: Second tree.
    :return: Whether they are identical, and the first difference found.

    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as d:
    ...     a, b = Path(d) / "a", Path(d) / "b"
    ...     for x in (a, b):
    ...         x.mkdir()
    ...         _ = (x / "f").write_text("1")
    ...     same = trees_identical(a, b)
    ...     _ = (b / "f").write_text("2")
    ...     same, trees_identical(a, b)
    ((True, ''), (False, 'f differs'))
    """
    names_a = sorted(p.relative_to(first) for p in first.rglob("*") if p.is_file())
    names_b = sorted(p.relative_to(second) for p in second.rglob("*") if p.is_file())
    if names_a != names_b:
        missing = sorted(set(map(str, names_a)) ^ set(map(str, names_b)))
        return False, f"file lists differ: {missing[:5]}"
    for name in names_a:
        if (first / name).read_bytes() != (second / name).read_bytes():
            return False, f"{name} differs"
    return True, ""


class Leg:
    """Generate, build, run, and check one formulation's W and chi applications."""

    def __init__(self, args: argparse.Namespace, report: Report) -> None:
        """
        Store the command-line settings.

        :param args: Parsed command-line arguments.
        :param report: Result recorder.
        """
        self.formulation: str = args.formulation
        self.work = Path(args.work_dir).resolve()
        self.launcher: List[str] = shlex.split(args.launcher)
        self.ranks: int = args.ranks
        self.build_jobs: int = args.build_jobs
        self.report = report
        self.solver_dir, self.exe_name, self.stem = SOLVER[self.formulation]
        self.executables: Dict[str, Path] = {}
        self.base_pars: Dict[str, str] = {}

    def mpi(self, ranks: int, exe: Path, *extra: str) -> List[str]:
        """
        Return the launcher argument vector for one solver run.

        :param ranks: Number of MPI ranks.
        :param exe: Solver executable.
        :param extra: Solver arguments.
        :return: Argument vector.
        """
        return self.launcher + ["-n", str(ranks), str(exe)] + list(extra)

    def record_versions(self) -> None:
        """
        Print the resolved versions of the dependencies that are not pinned.

        The apt packages, compilers, and ``requirements.txt`` packages come from
        the runner image and package indexes, which the workflow does not pin, so
        a pass is evidence only for the versions printed here.
        """
        commands = [
            [sys.executable, "--version"],
            ["clang-format", "--version"],
            ["cmake", "--version"],
            ["c++", "--version"],
            ["gfortran", "--version"],
            [self.launcher[0], "--version"],
            ["dpkg-query", "-W", "libopenmpi-dev", "libgsl-dev", "libopenblas-dev"],
        ]
        for requirement in (REPO_ROOT / "requirements.txt").read_text().split():
            try:
                version = importlib.metadata.version(requirement)
            except importlib.metadata.PackageNotFoundError:
                version = "not installed"
            print(f"version: {requirement}: {version}", flush=True)
        for argv in commands:
            try:
                completed = subprocess.run(
                    argv,
                    capture_output=True,
                    text=True,
                    timeout=TIMEOUT_VERSION,
                    check=False,
                    shell=False,
                )
                text = (completed.stdout or completed.stderr).strip().splitlines()
                print(
                    f"version: {' '.join(argv[:2])}: {' | '.join(text[:3])}", flush=True
                )
            except (OSError, subprocess.TimeoutExpired) as error:
                print(
                    f"version: {' '.join(argv[:2])}: unavailable ({error})", flush=True
                )

    def generate_and_build(self, conformal: str) -> None:
        """
        Generate one variant twice, compare the trees, and build the first.

        :param conformal: ``W`` or ``chi``.
        :raises CheckError: If the executable was not built.
        """
        trees = []
        for copy in (1, 2):
            project = self.work / f"gen-{conformal}-{copy}"
            cache = self.work / f"cache-{conformal}-{copy}"
            project.mkdir(parents=True)
            cache.mkdir(parents=True)
            env = dict(os.environ)
            env["XDG_CACHE_HOME"] = str(cache)
            env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
            run_checked(
                [
                    sys.executable,
                    "-m",
                    f"nrpy.examples.dendro_{self.formulation}",
                    "--project-dir",
                    str(project),
                    "--fd-order",
                    "6",
                    "--conformal-factor",
                    conformal,
                ],
                REPO_ROOT,
                self.work / f"generate-{conformal}-{copy}.log",
                TIMEOUT_GENERATE,
                env,
            )
            trees.append(project / self.solver_dir)
        same, detail = trees_identical(trees[0], trees[1])
        self.report.check(
            "G1",
            "generation determinism",
            f"{conformal}: two clean generations byte-identical",
            "identical" if same else detail,
            "identical",
            same,
        )
        build = self.work / f"build-{conformal}"
        run_checked(
            [
                "cmake",
                "-S",
                str(trees[0]),
                "-B",
                str(build),
                "-DCMAKE_BUILD_TYPE=Release",
                "-DCPU_ARCH=x86-64-v3",
                f"-DFETCHCONTENT_BASE_DIR={self.work / f'deps-{conformal}'}",
            ],
            self.work,
            self.work / f"configure-{conformal}.log",
            TIMEOUT_CONFIGURE,
        )
        run_checked(
            ["cmake", "--build", str(build), "--parallel", str(self.build_jobs)],
            self.work,
            self.work / f"build-{conformal}.log",
            TIMEOUT_BUILD,
        )
        exe = build / self.exe_name
        self.report.check(
            "B1",
            "build (compile/link compatibility)",
            f"{conformal}: standalone project builds {self.exe_name}",
            "built" if exe.is_file() else "missing",
            "executable exists",
            exe.is_file(),
        )
        if not exe.is_file():
            raise CheckError(f"{exe} was not built")
        self.executables[conformal] = exe
        self.base_pars[conformal] = (
            trees[0] / "pars" / f"{self.stem}.toml"
        ).read_text()

    def prepare(self, conformal: str, name: str, overrides: Dict[str, str]) -> Path:
        """
        Create a run directory with its parameter file and TwoPunctures file.

        :param conformal: ``W`` or ``chi``.
        :param name: Run name.
        :param overrides: Profile overrides on top of the common ones.
        :return: Run directory.
        :raises CheckError: If the horizon finder is enabled without VTU output,
            which aborts in Dendrolib.
        """
        merged = dict(COMMON_OVERRIDES)
        merged.update(overrides)
        if int(merged.get("AEH_SOLVER_FREQ", "0")) > 0 and (
            int(merged.get("BSSN_IO_OUTPUT_FREQ", "0")) == 0
        ):
            raise CheckError(
                "AEH_SOLVER_FREQ > 0 with BSSN_IO_OUTPUT_FREQ = 0 aborts in Dendrolib"
            )
        run_dir = self.work / conformal / name
        for sub in ("dat", "vtu", "cp", "bah"):
            (run_dir / sub).mkdir(parents=True, exist_ok=True)
        (run_dir / "ci.toml").write_text(
            apply_overrides(self.base_pars[conformal], merged)
        )
        tpid = self.work / conformal / TPID_FILE
        if tpid.exists() and not (run_dir / TPID_FILE).exists():
            os.link(tpid, run_dir / TPID_FILE)
        return run_dir

    def solve(
        self, conformal: str, name: str, overrides: Dict[str, str], ranks: int
    ) -> Path:
        """
        Run the solver in a prepared directory; it must succeed.

        :param conformal: ``W`` or ``chi``.
        :param name: Run name.
        :param overrides: Profile overrides.
        :param ranks: Number of MPI ranks.
        :return: Run directory.
        """
        run_dir = self.prepare(conformal, name, overrides)
        run_checked(
            self.mpi(ranks, self.executables[conformal], "ci.toml"),
            run_dir,
            run_dir / "run.log",
            TIMEOUT_RUN,
        )
        self.universal_checks(conformal, name, run_dir, overrides)
        return run_dir

    def universal_checks(
        self, conformal: str, name: str, run_dir: Path, overrides: Dict[str, str]
    ) -> None:
        """
        Check finiteness, cadence, and mesh size for every run.

        :param conformal: ``W`` or ``chi``.
        :param name: Run name.
        :param run_dir: Run directory.
        :param overrides: Profile overrides of the run.
        """
        tag = f"{conformal}/{name}"
        stdout = (run_dir / "run.log").read_text(errors="replace")
        alphas = [
            line.partition(" max_alpha=")[2].split()[0]
            for line in stdout.splitlines()
            if line.startswith("iteration=") and " max_alpha=" in line
        ]
        finite_alpha = bool(alphas) and all(math.isfinite(float(a)) for a in alphas)
        numbers: List[float] = []
        for path in sorted((run_dir / "dat").glob("*.dat")):
            if path.name.startswith("dgr_GW_"):
                for _, values in parse_modes(path):
                    numbers.extend(v for c in values for v in (c.real, c.imag))
            else:
                numbers.extend(v for row in parse_table(path)[1] for v in row)
        for path in sorted((run_dir / "bah").glob("BHaHAHA_diagnostics.ah*.gp")):
            numbers.extend(v for row in parse_table(path)[1] for v in row[:15])
        finite = finite_alpha and all(math.isfinite(v) for v in numbers)
        self.report.check(
            "U1",
            "runtime",
            f"{tag}: every output value finite",
            f"{len(numbers)} values, {len(alphas)} max_alpha lines",
            "all finite",
            finite,
        )
        header, rows = parse_table(run_dir / "dat" / "dgr_Constraints.dat")
        cadence = int(overrides["BSSN_TIME_STEP_OUTPUT_FREQ"])
        last = int(overrides["BSSN_MAX_ITERATIONS"])
        expected = list(range(0, last + 1, cadence))
        steps = [int(round(r[0])) for r in rows]
        times = [r[1] for r in rows]
        ordered = steps == expected and all(b > a for a, b in zip(times, times[1:]))
        self.report.check(
            "U2",
            "runtime",
            f"{tag}: constraint rows at the configured cadence",
            str(steps),
            str(expected),
            ordered,
        )
        nodes = column(header, rows, "unexcised_nodes")
        self.report.check(
            "U3",
            "runtime",
            f"{tag}: mesh below the node ceiling",
            f"max {int(max(nodes))}",
            f"<= {NODE_CEILING}",
            max(nodes) <= NODE_CEILING,
        )

    def run_variant(self, conformal: str) -> None:
        """
        Run and check one conformal-factor variant.

        :param conformal: ``W`` or ``chi``.
        :raises CheckError: If the TwoPunctures file lacks the NRPy tag.
        """
        exe = self.executables[conformal]
        tp_dir = self.prepare(conformal, "tpid", PROFILE_P)
        run_checked(
            self.mpi(1, exe, "--tpid", "ci.toml"),
            tp_dir,
            tp_dir / "run.log",
            TIMEOUT_RUN,
        )
        shutil.copyfile(tp_dir / TPID_FILE, self.work / conformal / TPID_FILE)
        data = (tp_dir / TPID_FILE).read_bytes()
        if not data.startswith(TPID_TAG):
            raise CheckError(f"{tp_dir / TPID_FILE} lacks the NRPy TwoPunctures tag")
        inputs = struct.unpack_from(f"<{TPID_INPUT_COUNT}d", data, 16)
        results = struct.unpack_from(
            f"<{TPID_RESULT_COUNT}d", data, 16 + 8 * TPID_INPUT_COUNT
        )

        run_a = self.solve(conformal, "A", PROFILE_P, self.ranks)
        self.check_run_a(conformal, run_a, results)

        stop4 = dict(PROFILE_P, BSSN_MAX_ITERATIONS="4")
        run_b = self.solve(conformal, "B", stop4, self.ranks)
        metadata = cast(
            Dict[str, Any],  # json.loads returns untyped JSON values
            json.loads((run_b / "cp" / "ci_cp_1_step.cp").read_text()),
        )
        self.check_checkpoint(conformal, metadata, inputs)
        restore = dict(PROFILE_P, BSSN_RESTORE_SOLVER="1")
        run_b = self.solve(conformal, "B", restore, self.ranks)
        for sub in ("dat", "bah", "vtu"):
            same, detail = trees_identical(run_a / sub, run_b / sub)
            self.report.check(
                "R1",
                "exact runtime identity",
                f"{conformal}: restart 0-4-8 equals uninterrupted run in {sub}/",
                "identical" if same else detail,
                "byte-identical",
                same,
            )

        run_c = self.solve(
            conformal, "C", dict(stop4, BSSN_CHECKPT_FREQ="0"), self.ranks // 2
        )
        self.compare_runs(
            f"{conformal}: C ({self.ranks // 2} ranks) vs A ({self.ranks} ranks)",
            run_a,
            run_c,
            [0, 4],
            True,
        )

        orders: Dict[int, Path] = {}
        for order in (4, 6, 8):
            orders[order] = self.solve(
                conformal,
                f"O{order}",
                dict(PROFILE_O, BSSN_ELE_ORDER=str(order)),
                self.ranks,
            )
        self.check_orders(conformal, orders)
        serial = self.solve(conformal, "O4s", dict(PROFILE_O, BSSN_ELE_ORDER="4"), 1)
        self.compare_runs(
            f"{conformal}: O4 1 rank vs {self.ranks} ranks",
            orders[4],
            serial,
            [0, 1, 2, 3, 4],
            False,
        )
        for name in ("A", "C", "O4", "O6", "O8", "O4s"):
            for sub in ("vtu", "cp"):
                shutil.rmtree(self.work / conformal / name / sub, ignore_errors=True)

    def check_run_a(
        self, conformal: str, run_a: Path, results: Sequence[float]
    ) -> None:
        """
        Check independent values, identities, and remeshing in run A.

        :param conformal: ``W`` or ``chi``.
        :param run_a: Run directory.
        :param results: TwoPunctures results ``mp, mm, mp_adm, mm_adm, E, J1, J2, J3``.
        """
        _, _, mp_adm, mm_adm, energy, _, _, j3 = results
        adm = parse_table(run_a / "dat" / "dgr_ADM.dat")[1]
        step0 = row_at(adm, 0)
        radius = step0[2]
        factor = 1.0 + energy / (2.0 * radius)
        e_expected = energy * factor**3
        j_expected = j3 / factor**2
        self.report.check(
            "X1",
            "regression/numerical, independent oracle",
            f"{conformal}: E({radius:g}) against E(1+E/2r)^3",
            f"{step0[3]:.6g} vs {e_expected:.6g}",
            "r = 100; relative <= 2e-4",
            radius == 100.0 and relative(step0[3], e_expected) <= 2e-4,
        )
        self.report.check(
            "X1",
            "regression/numerical, independent oracle",
            f"{conformal}: J_z({radius:g}) against J3/(1+E/2r)^2",
            f"{step0[9]:.6g} vs {j_expected:.6g}",
            "r = 100; relative <= 2e-4",
            radius == 100.0 and relative(step0[9], j_expected) <= 2e-4,
        )
        horizons = []
        for index, target in ((1, mp_adm), (2, mm_adm)):
            rows = parse_table(run_a / "bah" / f"BHaHAHA_diagnostics.ah{index}.gp")[1]
            found = sorted(int(round(r[0])) for r in rows)
            worst = max((relative(r[12], target) for r in rows), default=math.inf)
            self.report.check(
                "X2",
                "regression/numerical, independent oracle",
                f"{conformal}: horizon {index} irreducible mass against puncture ADM mass",
                f"steps {found}, max deviation {worst:.2e}",
                "steps [0, 4, 8]; relative <= 2e-3",
                found == [0, 4, 8] and worst <= 2e-3,
            )
            horizons.append({int(round(r[0])): r for r in rows})
        worst_mass = worst_x = 0.0
        for step in (0, 4, 8):
            if step in horizons[0] and step in horizons[1]:
                first, second = horizons[0][step], horizons[1][step]
                worst_mass = max(worst_mass, relative(first[12], second[12]))
                worst_x = max(worst_x, abs(first[2] + second[2]))
        # The horizon finder converges only to its own tolerance, so the two
        # horizons agree to that level, not to roundoff; the checkpoint check
        # I2 tests point reflection of the punctures at full precision.
        self.report.check(
            "X2",
            "numerical invariant",
            f"{conformal}: the two horizons are equal and point-reflected",
            f"mass {worst_mass:.1e}, centroid x {worst_x:.1e}",
            "mass <= 1e-5, x <= 1e-4 (horizon-finder tolerance)",
            worst_mass <= 1e-5 and worst_x <= 1e-4,
        )

        header, rows = parse_table(run_a / "dat" / "dgr_Constraints.dat")
        lam0 = column(header, rows, "LAMBDA_CONSTRAINT")[0]
        self.report.check(
            "I1",
            "numerical invariant",
            f"{conformal}: Lambda constraint at step 0",
            f"{lam0:.2e}",
            "<= 1e-12",
            lam0 <= 1e-12,
        )
        if self.formulation == "fccz4":
            h0 = column(header, rows, "H")[0]
            hz4 = column(header, rows, "H_Z4")[0]
            self.report.check(
                "I1",
                "numerical invariant",
                f"{conformal}: H_Z4 equals H at step 0",
                f"{hz4:.6g} vs {h0:.6g}",
                "relative <= 2e-5",
                relative(hz4, h0) <= 2e-5,
            )
            z4 = max(
                abs(column(header, rows, f"Z4constraintU{i}")[0]) for i in range(3)
            )
            self.report.check(
                "I1",
                "numerical invariant",
                f"{conformal}: Z4 vector at step 0",
                f"{z4:.2e}",
                "<= 1e-12",
                z4 <= 1e-12,
            )

        momentum = max(max(abs(r[4]), abs(r[5]), abs(r[6])) / abs(r[3]) for r in adm)
        spin = max(max(abs(r[7]), abs(r[8])) / abs(r[9]) for r in adm)
        self.report.check(
            "I2",
            "numerical invariant",
            f"{conformal}: P_ADM, J_x, J_y vanish",
            f"P/E {momentum:.1e}, J_xy/J_z {spin:.1e}",
            "<= 1e-10",
            momentum <= 1e-10 and spin <= 1e-10,
        )

        nodes = column(header, rows, "unexcised_nodes")
        step4 = row_at(adm, 4)
        continuity = relative(step4[3], step0[3])
        self.report.check(
            "M1",
            "runtime state plus numerical",
            f"{conformal}: remesh changes the mesh and preserves E",
            f"nodes {int(nodes[0])} -> {int(nodes[1])}, E change {continuity:.1e}",
            "nodes differ; relative E change <= 5e-3",
            nodes[0] != nodes[1] and continuity <= 5e-3,
        )

        modes: Dict[Tuple[int, int], List[Tuple[int, List[complex]]]] = {}
        for path in (run_a / "dat").glob("dgr_GW_l*_m*.dat"):
            ell, _, m = path.stem[len("dgr_GW_l") :].partition("_m")
            modes[(int(ell), int(m))] = parse_modes(path)
        expected_modes = {(l, m) for l in (2, 3, 4) for m in range(-l, l + 1)}
        radii = {len(values) for rows_lm in modes.values() for _, values in rows_lm}
        self.report.check(
            "I2",
            "runtime",
            f"{conformal}: Psi4 output covers l = 2-4, all m, two radii",
            f"{len(modes)} mode files, radii per row {sorted(radii)}",
            f"{len(expected_modes)} mode files, radii per row [2]",
            set(modes) == expected_modes and radii == {2},
        )
        if set(modes) != expected_modes:
            return  # the coverage check above has failed
        c22 = max(abs(v) for _, values in modes[(2, 2)] for v in values)
        odd = max(
            abs(v)
            for (l, m), rows_lm in modes.items()
            if m % 2
            for _, values in rows_lm
            for v in values
        )
        parity = max(
            abs(
                modes[(l, -m)][i][1][j] - (-1) ** l * modes[(l, m)][i][1][j].conjugate()
            )
            for (l, m) in modes
            if m > 0
            for i in range(len(modes[(l, m)]))
            for j in range(len(modes[(l, m)][i][1]))
        )
        self.report.check(
            "I2",
            "numerical invariant",
            f"{conformal}: odd-m Psi4 modes vanish (l = 2-4)",
            f"{odd / c22:.1e} of max|C22|",
            "<= 1e-8",
            odd / c22 <= 1e-8,
        )
        self.report.check(
            "I2",
            "numerical invariant",
            f"{conformal}: C(l,-m) = (-1)^l conj C(l,m) (l = 2-4)",
            f"{parity / c22:.1e} of max|C22|",
            "<= 1e-8",
            parity / c22 <= 1e-8,
        )

    def check_checkpoint(
        self,
        conformal: str,
        metadata: Dict[str, Any],  # untyped JSON values from json.loads
        inputs: Sequence[float],
    ) -> None:
        """
        Check the projected residual, reflection, and puncture motion at step 4.

        :param conformal: ``W`` or ``chi``.
        :param metadata: Step-4 checkpoint metadata.
        :param inputs: TwoPunctures solve inputs.
        """
        residual = float(str(metadata["NRPY_PROJECTED_ALGEBRAIC_RESIDUAL"]))
        self.report.check(
            "I1",
            "numerical invariant",
            f"{conformal}: projected algebraic residual at step 4",
            f"{residual:.2e}",
            "<= 1e-12",
            residual <= 1e-12,
        )
        centers = [float(x) for x in metadata["NRPY_EXCISION_CENTERS"]]
        x1, y1, z1, x2, y2, z2 = centers
        reflection = max(abs(x1 + x2), abs(y1 + y2), abs(z1), abs(z2))
        self.report.check(
            "I2",
            "numerical invariant",
            f"{conformal}: punctures point-reflected at step 4",
            f"{reflection:.1e}",
            "<= 1e-9",
            reflection <= 1e-9,
        )
        p_plus, p_minus = inputs[TPID_INDEX_P_PLUS_Y], inputs[TPID_INDEX_P_MINUS_Y]
        moved = y1 * p_plus > 0.0 and y2 * p_minus > 0.0
        self.report.check(
            "I3",
            "numerical sign check",
            f"{conformal}: punctures move along their momenta",
            f"y = ({y1:.2e}, {y2:.2e}), P_y = ({p_plus:.3g}, {p_minus:.3g})",
            "same signs",
            moved,
        )

    def compare_runs(
        self, label: str, first: Path, second: Path, steps: List[int], horizons: bool
    ) -> None:
        """
        Compare two runs of one binary at the given steps (check P1).

        :param label: Comparison label.
        :param first: Reference run directory.
        :param second: Run with a different rank count.
        :param steps: Steps to compare.
        :param horizons: Whether to compare horizon diagnostics.
        :raises CheckError: If a Psi4 file has no step in common with ``steps``.
        """
        used = 0.0
        nodes_equal = True
        for path in sorted((first / "dat").glob("*.dat")):
            other = second / "dat" / path.name
            if path.name.startswith("dgr_GW_"):
                mine = dict(parse_modes(path))
                theirs = dict(parse_modes(other))
                shared = [s for s in steps if s in mine and s in theirs]
                if not shared:
                    raise CheckError(f"{path.name}: no common steps in {steps}")
                pairs = [
                    (a, b)
                    for s in shared
                    for ca, cb in zip(mine[s], theirs[s])
                    for a, b in ((ca.real, cb.real), (ca.imag, cb.imag))
                ]
            else:
                header, rows_a = parse_table(path)
                rows_b = parse_table(other)[1]
                if header and "unexcised_nodes" in header:
                    index = header.index("unexcised_nodes")
                    nodes_equal = nodes_equal and all(
                        row_at(rows_a, s)[index] == row_at(rows_b, s)[index]
                        for s in steps
                    )
                pairs = [
                    (a, b)
                    for s in steps
                    for a, b in zip(row_at(rows_a, s), row_at(rows_b, s))
                ]
            for a, b in pairs:
                used = max(used, abs(a - b) / (1.0e-8 * max(abs(a), abs(b)) + 1.0e-12))
        if horizons:
            for index in (1, 2):
                name = f"BHaHAHA_diagnostics.ah{index}.gp"
                rows_a = parse_table(first / "bah" / name)[1]
                rows_b = parse_table(second / "bah" / name)[1]
                for s in steps:
                    for a, b in zip(row_at(rows_a, s)[:15], row_at(rows_b, s)[:15]):
                        used = max(
                            used, abs(a - b) / (1.0e-8 * max(abs(a), abs(b)) + 1.0e-12)
                        )
        self.report.check(
            "P1",
            "numerical invariant",
            f"{label} at steps {steps}",
            f"largest fraction of tolerance used {used:.1e}; nodes equal {nodes_equal}",
            "|a-b| <= 1e-8 max(|a|,|b|) + 1e-12; equal unexcised_nodes",
            used <= 1.0 and nodes_equal,
        )

    def check_orders(self, conformal: str, orders: Dict[int, Path]) -> None:
        """
        Check the finite-difference orders on one octree (check O1).

        :param conformal: ``W`` or ``chi``.
        :param orders: Run directory per order.
        """
        elements = {}
        h0 = {}
        for order, run_dir in orders.items():
            stdout = (run_dir / "run.log").read_text(errors="replace")
            passes = [
                line.partition(" elements=")[2].partition("->")[2].split()[0]
                for line in stdout.splitlines()
                if line.startswith("initial-grid remesh pass=")
            ]
            elements[order] = int(passes[-1]) if passes else -1
            reported = [
                line.partition(" FD=")[2].split()[0]
                for line in stdout.splitlines()
                if ": profile=" in line and " FD=" in line
            ]
            self.report.check(
                "O1",
                "runtime",
                f"{conformal}: FD{order} run reports its order",
                reported[0] if reported else "none",
                str(order),
                reported == [str(order)],
            )
            header, rows = parse_table(run_dir / "dat" / "dgr_Constraints.dat")
            h0[order] = column(header, rows, "H")[0]
            lam0 = column(header, rows, "LAMBDA_CONSTRAINT")[0]
            self.report.check(
                "I1",
                "numerical invariant",
                f"{conformal}: FD{order} Lambda constraint at step 0",
                f"{lam0:.2e}",
                "<= 1e-12",
                lam0 <= 1e-12,
            )
        same_octree = len(set(elements.values())) == 1 and -1 not in elements.values()
        self.report.check(
            "O1",
            "runtime",
            f"{conformal}: one octree for FD4, FD6, FD8",
            str(elements),
            "equal element counts",
            same_octree,
        )
        ordered = h0[6] <= h0[4] / 2.0 and h0[8] <= h0[6] / 2.0
        self.report.check(
            "O1",
            "numerical ordering (not convergence)",
            f"{conformal}: step-0 H decreases with order",
            f"FD4 {h0[4]:.3e}, FD6 {h0[6]:.3e}, FD8 {h0[8]:.3e}",
            "H6 <= H4/2 and H8 <= H6/2",
            same_octree and ordered,
        )

    def compare_variants(self) -> None:
        """Compare the W and chi builds of this formulation (check V1)."""
        run_w, run_chi = self.work / "W" / "A", self.work / "chi" / "A"
        adm_w = parse_table(run_w / "dat" / "dgr_ADM.dat")[1]
        adm_c = parse_table(run_chi / "dat" / "dgr_ADM.dat")[1]
        header_w, cons_w = parse_table(run_w / "dat" / "dgr_Constraints.dat")
        header_c, cons_c = parse_table(run_chi / "dat" / "dgr_Constraints.dat")
        nodes_equal = (
            column(header_w, cons_w, "unexcised_nodes")[0]
            == column(header_c, cons_c, "unexcised_nodes")[0]
        )
        e0 = relative(row_at(adm_w, 0)[3], row_at(adm_c, 0)[3])
        j0 = relative(row_at(adm_w, 0)[9], row_at(adm_c, 0)[9])
        c22 = max(
            abs(v)
            for _, values in parse_modes(run_w / "dat" / "dgr_GW_l2_m2.dat")
            for v in values
        )
        worst_mode = 0.0
        for path in sorted((run_w / "dat").glob("dgr_GW_l*_m*.dat")):
            first = dict(parse_modes(path))[0]
            second = dict(parse_modes(run_chi / "dat" / path.name))[0]
            worst_mode = max(worst_mode, max(abs(a - b) for a, b in zip(first, second)))
        self.report.check(
            "V1",
            "differential numerical",
            "W vs chi at step 0",
            f"nodes equal {nodes_equal}, E {e0:.1e}, J_z {j0:.1e}, Psi4 {worst_mode / c22:.1e}",
            "equal nodes; E, J_z relative <= 2e-5; Psi4 <= 1e-5 of max|C22|",
            nodes_equal and e0 <= 2e-5 and j0 <= 2e-5 and worst_mode <= 1e-5 * c22,
        )
        evolution = max(
            max(relative(a[3], b[3]), relative(a[9], b[9]))
            for a, b in zip(adm_w, adm_c)
        )
        ratios = []
        for name in ("H", "M_CONSTRAINT"):
            for a, b in zip(
                column(header_w, cons_w, name), column(header_c, cons_c, name)
            ):
                ratios.append(a / b if b != 0.0 else math.inf)
        self.report.check(
            "V1",
            "differential numerical",
            "W vs chi while evolving",
            f"E, J_z {evolution:.1e}; H, M ratios {min(ratios):.2f}-{max(ratios):.2f}",
            "E, J_z relative <= 1e-3; ratios in [0.5, 2] (gross-error check)",
            evolution <= 1e-3 and all(0.5 <= r <= 2.0 for r in ratios),
        )

    def negative(
        self, label: str, argv: List[str], run_dir: Path, expected: str
    ) -> None:
        """
        Run a case that must fail promptly with a named diagnostic (check N1).

        :param label: Case label.
        :param argv: Argument vector.
        :param run_dir: Working directory.
        :param expected: Text that must appear in the output; a trailing newline
            requires it to end a line.
        """
        safe = "".join(c if c.isalnum() else "_" for c in label)
        log = run_dir / f"negative-{safe}.log"
        status, timed_out = run_logged(argv, run_dir, log, TIMEOUT_NEGATIVE)
        text = log.read_text(errors="replace")
        ok = (not timed_out) and 0 < status < 124 and expected in text
        if not ok:
            print(tail(log), flush=True)
        self.report.check(
            "N1",
            "runtime error behavior",
            label,
            (
                "timed out"
                if timed_out
                else f"exit {status}, text {'found' if expected in text else 'missing'}"
            ),
            f"exit status 1-123 with '{expected.strip()}'",
            ok,
        )

    def run_negatives(self) -> None:
        """
        Run the invalid-input cases once for this formulation.

        :raises CheckError: If a copied checkpoint does not record its writer's
            formulation.
        """
        suffix = {"bssn": "BSSN", "fccz4": "fCCZ4"}[self.formulation]
        names = {"W": suffix, "chi": f"{suffix}_chi"}
        for writer, reader in (("W", "chi"), ("chi", "W")):
            target = self.work / reader / f"restore-from-{writer}"
            shutil.copytree(self.work / writer / "B", target)
            for slot in (0, 1):
                stored = json.loads(
                    (target / "cp" / f"ci_cp_{slot}_step.cp").read_text()
                ).get("NRPY_FORMULATION")
                if stored != names[writer]:
                    raise CheckError(
                        f"{writer} checkpoint slot {slot} records formulation {stored}"
                    )
            for sub in ("dat", "bah", "vtu"):
                shutil.rmtree(target / sub, ignore_errors=True)
                (target / sub).mkdir()
            (target / "ci.toml").write_text(
                apply_overrides(
                    self.base_pars[reader],
                    dict(COMMON_OVERRIDES, **PROFILE_P, BSSN_RESTORE_SOLVER="1"),
                )
            )
            self.negative(
                f"{reader} executable restoring the {writer} checkpoint",
                self.mpi(self.ranks, self.executables[reader], "ci.toml"),
                target,
                f"Checkpoint metadata does not match {names[reader]}\n",
            )
        exe = self.executables["W"]
        o4 = dict(PROFILE_O, BSSN_ELE_ORDER="4")
        cases = [
            (
                "--tpid on 2 ranks",
                o4,
                2,
                ["--tpid"],
                "--tpid requires exactly one MPI task",
                False,
            ),
            (
                "TPID_PAR_B changed after --tpid",
                dict(o4, TPID_PAR_B="4.5"),
                2,
                [],
                "NRPy TwoPunctures file does not match the input parameters",
                False,
            ),
            (
                "TwoPunctures file absent",
                o4,
                2,
                [],
                "cannot open NRPy TwoPunctures file",
                True,
            ),
            (
                "BSSN_ELE_ORDER = 5",
                dict(PROFILE_O, BSSN_ELE_ORDER="5"),
                2,
                [],
                "BSSN_ELE_ORDER must be 4, 6, or 8",
                False,
            ),
            (
                "BSSN_REFINEMENT_MODE = 0",
                dict(o4, BSSN_REFINEMENT_MODE="0"),
                2,
                [],
                "generated Dendro solver supports BH_WAMR",
                False,
            ),
            (
                "BSSN_CFL_FACTOR = 3.0",
                dict(o4, BSSN_CFL_FACTOR="3.0"),
                2,
                [],
                "constraint diagnostics found a nonfinite constraint",
                False,
            ),
        ]
        for index, (label, overrides, ranks, extra, expected, drop_tp) in enumerate(
            cases
        ):
            run_dir = self.prepare("W", f"N{index}", overrides)
            if drop_tp:
                (run_dir / TPID_FILE).unlink()
            self.negative(
                label, self.mpi(ranks, exe, *extra, "ci.toml"), run_dir, expected
            )
        self.negative(
            "no arguments", self.mpi(1, exe), self.work / "W", f"usage: {self.exe_name}"
        )

    def run(self) -> None:
        """
        Run the whole leg; the working directory is always removed.

        :raises CheckError: If the work directory already exists.
        """
        if self.work.exists():
            raise CheckError(f"work directory {self.work} already exists")
        self.work.mkdir(parents=True)
        try:
            self.record_versions()
            for conformal in ("W", "chi"):
                self.generate_and_build(conformal)
            for conformal in ("W", "chi"):
                self.run_variant(conformal)
            self.compare_variants()
            self.run_negatives()
        finally:
            shutil.rmtree(self.work, ignore_errors=True)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Parse arguments, run one formulation leg, and summarize.

    :param argv: Command-line arguments, or None for ``sys.argv``.
    :return: Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "--formulation",
        choices=sorted(SOLVER),
        required=True,
        help="application to generate and check",
    )
    parser.add_argument(
        "--work-dir", required=True, help="new directory, removed on exit"
    )
    parser.add_argument(
        "--launcher",
        default="mpiexec --oversubscribe --bind-to none",
        help="MPI launcher command, before its -n option",
    )
    parser.add_argument(
        "--ranks",
        type=int,
        default=4,
        help="MPI ranks of the main runs; at least 4, and neither it nor half of it 3",
    )
    parser.add_argument("--build-jobs", type=int, default=4, help="parallel build jobs")
    args = parser.parse_args(argv)
    if args.ranks < 4 or args.ranks // 2 == 3:
        parser.error(
            "--ranks must be at least 4, and 6 and 7 are refused because a 3-rank "
            "run misroutes the Dendrolib horizon-finder checkpoint"
        )
    report = Report()
    try:
        Leg(args, report).run()
    except (CheckError, OSError) as error:
        report.failures.append(f"aborted: {error}")
        print(f"[FAIL] aborted: {error}", flush=True)
    print(f"{report.count} checks, {len(report.failures)} failures", flush=True)
    for failure in report.failures:
        print(f"  FAILED {failure}", flush=True)
    return 1 if report.failures else 0


if __name__ == "__main__":
    if len(sys.argv) == 1:
        import doctest

        doctest_results = doctest.testmod()
        if doctest_results.failed > 0:
            print(
                f"Doctest failed: {doctest_results.failed} of "
                f"{doctest_results.attempted} test(s)"
            )
            sys.exit(1)
        print(f"Doctest passed: All {doctest_results.attempted} test(s) passed")
    else:
        sys.exit(main())
