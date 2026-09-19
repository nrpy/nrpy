# Dendrolib capability mini-tests

`dendrolib_capability_test.cpp` proves, against the selected Dendrolib checkout,
the host-interface requirements the generated Dendro solver assumes: the scalar ABI, the
padded block dimensions, the padding rule, the unzip offsets, the variable-major
x-fastest layout, the padded origin, and halo validity.

The test program is not built by NRPy and is not part of any generated project. It
is run by hand when the selected host changes.

## Build

Clone Dendrolib, select the revision under qualification, and build it. Then
compile the test program with that build's own definitions and include paths. Record
the selected revision with the active qualification evidence, not in this
durable procedure.

```bash
git clone https://github.com/paralab/Dendro-5.01 dendrolib
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
CAPTEST_ORDERS=4,6,8 mpirun -n 1 ./captest   # generated profile orders
```

A clean run completes every requested order, ends with `CAPABILITY_TESTS_OK`,
and exits zero on every rank. Rank zero selects `CAPTEST_ORDERS`: a nonempty,
comma-separated list of positive even element orders (for example `4,6` or
`4,6,8`). Empty, nonnumeric, or malformed overrides fail qualification.
Mesh/vector setup status is checked collectively before halo exchange, so a
local setup failure fails the MPI-wide result without stranding peer ranks.

## Exercising the checkers

A passing checker proves nothing until it has been shown to fail. Every axis
has a known-bad input, selected through `CAPTEST_INJECT`:

| Value | Defect introduced | Axis that must report `FAILED` |
| --- | --- | --- |
| `mesh` | fails mesh setup on rank 1 (rank 0 for a singleton) | `mesh_setup` |
| `zipped` | clears the zipped allocation on that rank | `vector_setup` |
| `unzipped` | clears the unzipped allocation on that rank | `vector_setup` |
| `scalar` | expects a four-byte scalar | `scalar_abi` |
| `padding` | expects one more ghost point than the element order implies | `padding` |
| `dimensions` | expects one extra point per axis | `dimensions` |
| `offsets` | inflates the block extent past the unzip stride | `offsets` |
| `layout` | transposes the fastest and slowest indices | `layout` |
| `origin` | treats padded index zero as the block corner | `origin` |
| `ghost` | blanks the halo the host filled | `ghost_validity` |

Run each before quoting a clean result:

```bash
for ax in mesh zipped unzipped scalar padding dimensions offsets layout origin ghost; do
  CAPTEST_INJECT=$ax mpirun -n 2 ./captest | grep FAILED
done
```

## Generated real-host qualification

Each generated formulation exposes a production library when its module is
added to a CMake tree containing `dendro5`. Enabling
`NRPY_DENDRO_BUILD_DRIVERS` also builds its real-host qualification executable.
Existing Dendro commands that define `BSSN_STANDALONE_HOST` or
`FCCZ4_STANDALONE_HOST` make the corresponding driver default on. When the
conventional target name is free, they may build the compatibility target
`bssnSolver` or `fccz4Solver`; that target copies the separately named NRPy
qualification executable into the conventional Dendro solver path. New
commands should use `NRPY_DENDRO_BUILD_DRIVERS` and the
`nrpy_*_dendro_qualify` target directly.
It uses actual blocks and vectors, synchronous Dendrolib halo exchange, and the
selected host's `ts::ETS` RK4 initialization, stage updates, and cleanup. This
is a fixed-mesh, serial-CPU-per-rank Minkowski qualification with generated
analytic exterior data. It does not qualify general physical boundaries, remeshing, LTS,
checkpoint/restart, output routines, GPU execution, or threaded kernels.

Use an isolated working directory. Set `NRPY_SOURCE` to the absolute path of
the source checkout containing this README. The following commands generate
solver files in the isolated working directory and modify only its `host` clone:

```bash
export NRPY_SOURCE=/absolute/path/to/NRPy
mkdir nrpy-real-host
cd nrpy-real-host
git clone https://github.com/paralab/Dendro-GR host
git clone https://github.com/paralab/Dendro-5.01 dendrolib
export XDG_CACHE_HOME="$PWD/cache"
PYTHONPATH="$NRPY_SOURCE" python -m nrpy.examples.dendro_fccz4 \
  --project-dir "$PWD/generated"
cat >> host/CMakeLists.txt <<'CMAKE'
add_subdirectory("${NRPY_GENERATED}/Dendro-GR/nrpy_fccz4" nrpy_fccz4)
add_executable(nrpy_runtime_test
  "${NRPY_SOURCE}/nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp")
target_link_libraries(nrpy_runtime_test PRIVATE nrpy_fccz4_dendro MPI::MPI_CXX)
target_compile_definitions(nrpy_runtime_test PRIVATE
  RUNTIME_HEADER="fccz4Ctx.h" RUNTIME_NAMESPACE=nrpy::fccz4)
CMAKE
cmake -S host -B build -DWITH_CUDA=OFF -DNRPY_DENDRO_BUILD_DRIVERS=ON \
  -DNRPY_DENDRO_BUILD_TESTS=OFF \
  -DNRPY_SOURCE="$NRPY_SOURCE" -DNRPY_GENERATED="$PWD/generated" \
  -DDENDRO_dendrolib_DIR="$PWD/dendrolib"
cmake --build build --target nrpy_fccz4_dendro_qualify nrpy_runtime_test -j2
OMP_NUM_THREADS=1 timeout 300 mpiexec -n 1 build/nrpy_runtime_test
OMP_NUM_THREADS=1 timeout 300 mpiexec -n 2 build/nrpy_runtime_test
OMP_NUM_THREADS=1 timeout 300 mpiexec -n 2 \
  build/nrpy_fccz4/nrpy_fccz4_dendro_qualify
```

For BSSN, generate with `nrpy.examples.dendro_bssn`, use
`nrpy_bssn_dendro`, `nrpy_bssn_dendro_qualify`, `bssnCtx.h`, and
`nrpy::bssn`. Both generated modules may be added to the same parent CMake
project because their production and qualification target names are distinct.

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
nodes. It calls `rhs_blkwise` on one block with a nonzero component offset,
proves all other storage retains a sentinel, and compares that result with
the same selected-block interior retained from whole-vector `rhs`. This
comparison uses `B^0=0.01` and `eta=1`, and requires a nonzero RHS magnitude,
so zero Minkowski data cannot hide a component-routing error. The test then
compares with zero-offset component-major storage passed to `rhs_blk`. It also
proves the three no-op block hooks preserve their input bytes and compares
`post_timestep_blk` with whole-vector algebraic projection. Finally, the test
changes `eta` from zero to one on a state with `B^0=0.01` and checks the
analytic RHS difference `-0.01`, with other components unchanged. A passing
run prints `REAL_TRANSPORT PASS`, `REAL_CALLBACKS PASS`, and
`REAL_PARAMETER PASS`.

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

for fault in blockwise_null_ids blockwise_bad_id blockwise_bad_dof \
             block_null block_bad_dof block_bad_id \
             projection_null projection_bad_dof projection_bad_id; do
  if OMP_NUM_THREADS=1 timeout 60 mpiexec -n 1 build/nrpy_runtime_test "$fault" \
      > "$fault.txt" 2>&1; then
    echo "ERROR: $fault escaped detection"
    exit 1
  else
    status=$?
    test "$status" -ne 124 || exit 1
  fi
done
```

Inspect the diagnostics as well as exit status: offset/halo must produce a
transport mismatch or unfilled block diagnostic, and nonfinite must report
`nonfinite evolved state`. Transport faults are injected on the last rank, so
the same source runs on one or two ranks. Callback faults check null pointers,
wrong field counts, and out-of-range local block identifiers. Every failure
must terminate the entire MPI job without waiting for a timeout.

The solver defaults to 100 steps of `dt=0.001`; `--steps N` and `--dt T` change
those values. `-t FILE` reads TOML on rank 0, broadcasts its contents, binds
registered `[params]` values, and validates the optional generated profile.
Unknown keys, profile mismatch, malformed values, and nonfinite values fail
collectively. The sample file under `generated/Dendro-GR/nrpy_fccz4/pars/` is
usable by the real build. Parameters used by the generated block RHS, including
`eta`, `kappa1`, and `kappa2`, are forwarded through its registered signature.
Parameters used only by standalone qualification kernels remain generated
struct members but are not real-host TOML keys. The qualification mesh/domain
and Minkowski initial data remain fixed; geometry and perturbation defaults in
the shared registry do not select another problem.
The standalone executable still rejects `-t`.

The checked solver result requires finite output, the exact RK4 step/time and
projection schedule (one initial projection plus five per step), drift at most
`1e-11`, and projection residual at most `1e-13`. RHS and constraints must be at
most `256 * epsilon(double) / h_min^2`: coarse/fine interpolation introduces
roundoff in constant data, and second derivatives amplify it by inverse spacing
squared. This is a roundoff-scaled Minkowski check, not a convergence result.
`REAL_MINKOWSKI PASS` reports every measured maximum and the derivative bound.
