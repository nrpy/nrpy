# Dendrolib capability mini-tests

Last real run: 09-06-2026, Ubuntu 24.04 x86_64, GCC 13.3.0, Open MPI 4.1.6,
CMake 3.28.3, against a Dendrolib clone at
`246043709e806021fcfc011fe657b8bf964cae4c` (`git -C dendrolib rev-parse HEAD`),
ending in `CAPABILITY_TESTS_OK`.  Whenever the pinned commit moves, record this
paragraph again -- platform, toolchain, the clone's `git rev-parse HEAD`, and the
closing `CAPABILITY_TESTS_OK` line.  A result without them is not evidence.

`dendrolib_capability_test.cpp` proves, against Dendrolib commit
`246043709e806021fcfc011fe657b8bf964cae4c`, the host contract the generated
Dendro solver assumes: the scalar ABI, the padded block dimensions, the padding
rule, the unzip offsets, the variable-major x-fastest layout, the padded origin,
and halo validity. Results are recorded on the Dendro validation page in the
NRPy knowledge base.

The harness is not built by NRPy and is not part of any generated project. It
is run by hand when the pin changes.

## Build

Clone and build the pinned Dendrolib, then compile the harness with that
build's own definitions and include paths:

```bash
git clone https://github.com/paralab/Dendro-5.01 dendrolib
git -C dendrolib checkout 246043709e806021fcfc011fe657b8bf964cae4c
cmake -S dendrolib -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target dendro5 -j "$(nproc)"

DEF=$(grep '^CXX_DEFINES' build/CMakeFiles/dendro5.dir/flags.make | cut -d= -f2-)
INC=$(grep '^CXX_INCLUDES' build/CMakeFiles/dendro5.dir/flags.make | cut -d= -f2-)
mpicxx -std=gnu++17 -O2 -fopenmp $DEF $INC \
  -o captest dendrolib_capability_test.cpp \
  -Lbuild -ldendro5 -Lbuild/_deps/spdlog-build -lspdlog \
  -lopenblas -lgfortran -lz -lm
```

## Run

```bash
mpirun -n 1 ./captest      # prints one PROVEN/FAILED line per axis
mpirun -n 2 ./captest      # same axes across a distributed block list
CAPTEST_ORDERS=10 mpirun -n 1 ./captest   # padding 5, the eighth-order probe
```

A clean run ends with `CAPABILITY_TESTS_OK` and exits zero.

## Exercising the checkers

A passing checker proves nothing until it has been shown to fail. Every axis
has a known-bad input, selected through `CAPTEST_INJECT`:

| Value | Defect introduced | Axis that must report `FAILED` |
| --- | --- | --- |
| `scalar` | expects a four-byte scalar | `scalar_abi` |
| `padding` | expects one more ghost point than the element order implies | `padding` |
| `dimensions` | expects one extra point per axis | `dimensions` |
| `offsets` | inflates the block extent past the unzip stride | `offsets` |
| `layout` | transposes the fastest and slowest indices | `layout` |
| `origin` | treats padded index zero as the block corner | `origin` |
| `ghost` | blanks the halo the host filled | `ghost_validity` |

Run each before quoting a clean result:

```bash
for ax in scalar padding dimensions offsets layout origin ghost; do
  CAPTEST_INJECT=$ax mpirun -n 1 ./captest | grep FAILED
done
```

## Generated real-host qualification

The generated fCCZ4 solver has a separate real-host build selected by
`FCCZ4_STANDALONE_HOST=OFF`. It uses actual blocks and vectors, synchronous
Dendrolib halo exchange, and the pinned `ts::ETS` RK4 lifecycle. This is a
fixed-mesh, serial-CPU-per-rank Minkowski qualification with analytic exterior
data. It does not qualify general physical boundaries, remeshing, LTS,
checkpoint/restart, output pipelines, GPU execution, or threaded kernels.

Use an isolated working directory. Set `NRPY_SOURCE` to the absolute path of
the source checkout containing this README. The following commands generate
owned products and modify only the isolated host clone:

```bash
export NRPY_SOURCE=/absolute/path/to/NRPy
mkdir nrpy-real-host
cd nrpy-real-host
git clone https://github.com/paralab/Dendro-GR host
git -C host checkout b3261e2a0d3457781b11d63ac5ab38375ffab93b
git clone https://github.com/paralab/Dendro-5.01 dendrolib
git -C dendrolib checkout 246043709e806021fcfc011fe657b8bf964cae4c
export XDG_CACHE_HOME="$PWD/cache"
PYTHONPATH="$NRPY_SOURCE" python -m nrpy.examples.dendro_fccz4 \
  --project-dir "$PWD/generated"
cat >> host/CMakeLists.txt <<'CMAKE'
add_subdirectory("${NRPY_GENERATED}/Dendro-GR/FCCZ4_GR" nrpy_fccz4)
add_executable(nrpy_runtime_test
  "${NRPY_SOURCE}/nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp")
target_link_libraries(nrpy_runtime_test PRIVATE fccz4_common MPI::MPI_CXX)
target_compile_definitions(nrpy_runtime_test PRIVATE
  RUNTIME_HEADER="fccz4Ctx.h" RUNTIME_NAMESPACE=fccz4)
CMAKE
cmake -S host -B build -DWITH_CUDA=OFF -DFCCZ4_STANDALONE_HOST=OFF \
  -DNRPY_SOURCE="$NRPY_SOURCE" -DNRPY_GENERATED="$PWD/generated" \
  -DDENDRO_dendrolib_DIR="$PWD/dendrolib" \
  -DDENDRO_dendrolib_GIT_TAG=246043709e806021fcfc011fe657b8bf964cae4c
cmake --build build --target fccz4Solver nrpy_runtime_test -j2
OMP_NUM_THREADS=1 timeout 300 mpiexec -n 2 build/nrpy_runtime_test
OMP_NUM_THREADS=1 timeout 300 mpiexec -n 2 build/nrpy_fccz4/fccz4Solver
```

The host's usual BLAS/LAPACK, MPI, C++17, and dependency-fetch requirements
apply. Cached `toml11` and `spdlog` source directories may be supplied through
CMake's `FETCHCONTENT_SOURCE_DIR_TOML11` and
`FETCHCONTENT_SOURCE_DIR_SPDLOG` variables. Verify both local checkout commits;
the generated tag comparison warns about a different configured tag but cannot
prove which commit a manually supplied source directory contains.

`runtime_integration_test.cpp` poisons all non-owned zipped slots and unzipped
storage before exchange. It checks every node in the host receive scatter map,
all in-domain padded block values against component-distinct affine data,
physical padded origins, nonzero block offsets, and nonconstant zip results.
The MPI mesh can reserve unused ghost slots; those are not promised receive
nodes. The test also changes `eta` from zero to one on a state with `B^0=0.01`
and checks the analytic RHS difference `-0.01`, with other components unchanged.
A passing run prints `REAL_TRANSPORT PASS` and `REAL_PARAMETER PASS`.

Exercise the negative cases before relying on these checkers:

```bash
for fault in offset halo nonfinite; do
  if OMP_NUM_THREADS=1 timeout 60 mpiexec -n 2 build/nrpy_runtime_test "$fault" \
      > "$fault.txt" 2>&1; then
    echo "ERROR: $fault escaped detection"
    exit 1
  else
    status=$?
    test "$status" -ne 124 || exit 1  # a hang is not successful detection
  fi
done
```

Inspect the diagnostics as well as exit status: offset/halo must produce a
transport mismatch or unfilled block diagnostic, and nonfinite must report
`nonfinite evolved state`. Each fault is injected on rank 1 only; failure must
terminate the entire MPI job without waiting for a timeout.

The solver defaults to 100 steps of `dt=0.001`; `--steps N` and `--dt T` change
those values. `-t FILE` reads TOML on rank 0, broadcasts its contents, binds
registered `[params]` values, and validates the optional generated profile.
Unknown keys, profile mismatch, malformed values, and nonfinite values fail
collectively. The sample file under `generated/Dendro-GR/FCCZ4_GR/pars/` is
usable by the real build. Parameters used by generated kernels, including
`eta`, `kappa1`, and `kappa2`, are forwarded through their registered signatures.
The qualification mesh/domain and Minkowski initial data remain fixed; geometry
and perturbation defaults in the shared registry do not select another problem.
The standalone executable still rejects `-t`.

The checked solver result requires finite output, the exact RK4 step/time and
projection schedule (one initial projection plus five per step), drift at most
`1e-11`, and projection residual at most `1e-13`. RHS and constraints must be at
most `256 * epsilon(double) / h_min^2`: coarse/fine interpolation introduces
roundoff in constant data, and second derivatives amplify it by inverse spacing
squared. This is a roundoff-scaled Minkowski check, not a convergence result.
`REAL_MINKOWSKI PASS` reports every measured maximum and the derivative bound.

Qualified on 09-07-2026 with Ubuntu 24.04 x86_64, GCC 13.3.0, CMake 3.28.3,
Open MPI 4.1.6, Python 3.12.3, and SymPy 1.14.0 at the two pins above, using
FD4, KO off, double precision, and two active MPI ranks. The transport run had
18 local blocks on one rank, 27 nonzero-offset blocks across ranks, 32,796
in-domain halo points, 3,259 receive nodes, and maximum affine/zip error
`5.9117155615240335e-12`. The 100-step run reached time
`0.10000000000000007`, RHS `2.6044445380985459e-11`, constraints
`2.6903897068313774e-11`, drift `2.7544910351229485e-13`, and projection
residual `9.9920072216264089e-16`; `h_min=0.020833333333333332` gives derivative
bound `1.3096723705530167e-10`. This was a local pinned-host qualification;
the repository workflow still runs the standalone vehicle.

The same qualification detected rank-local offset, halo, and NaN faults with
nonzero MPI exits before the 60-second timeout. A nondefault `eta=1.25` file
passed a two-step run with `dt=0.0005`; unknown parameter keys, NaN, wrong
value types, and a mismatched FD profile each failed collectively.
