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
