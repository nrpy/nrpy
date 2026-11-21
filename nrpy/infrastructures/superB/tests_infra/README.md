# Inverse-fisheye superB startup check

This manual check compiles selected generated C++ functions against a generated
N=8 inverse-fisheye superB project and Charm++ headers. It requires a working
Charm++ installation, an MPI C++ compiler, and the generated project. The check
does not run an evolution.

From the NRPy repository, generate the project in a scratch directory so an
existing `project/` tree is left intact:

```bash
REPO=$PWD
mkdir -p /tmp/nrpy-inverse-fisheye-startup
cd /tmp/nrpy-inverse-fisheye-startup
PYTHONPATH="$REPO" python3 \
  "$REPO/nrpy/examples/superb_paper2/superB_blackhole_spectroscopy_8Mseparation.py" \
  --inverse-fisheye
```

Run the check with `charmc` on `PATH`. Set `CXX` if the MPI C++ compiler is not
named `mpicxx`:

```bash
PYTHONPATH="$REPO" python3 \
  "$REPO/nrpy/infrastructures/superB/tests_infra/check_superb_inverse_fisheye_startup.py" \
  /tmp/nrpy-inverse-fisheye-startup/project/superB_blackhole_spectroscopy_8Mseparation_inverse_fisheye
```

The program checks that the local outer-boundary lists have the same point
order, ownership, face flags, radius, and radial derivatives as the global
boundary builder. It also compares the distributed CFL minimum with an
exhaustive grid scan and compares cached coordinate components with direct
inverse solves and Jacobians. An assertion failure exits with a nonzero status.

Running the Python file without a project argument executes its doctests and
does not require Charm++ or a generated project. This allows the repository's
Python static-analysis job to inspect the file.
