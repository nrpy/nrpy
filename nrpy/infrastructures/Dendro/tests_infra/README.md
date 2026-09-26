# Dendrolib capability mini-tests

`dendrolib_capability_test.cpp` proves, against the selected Dendrolib checkout,
the host-interface requirements the generated Dendro solver assumes: the scalar ABI, the
padded block dimensions, the padding rule, the unzip offsets, the variable-major
x-fastest layout, the padded origin, and halo validity.

The test program is not built by NRPy and is not part of any generated project. It
is run by hand when the selected host changes.

Validation of the complete generated applications (generation, build, MPI
evolution, restart, and the numerical checks) is the `dendro-validation`
GitHub Actions job, which runs
`nrpy/examples/tests/dendro_application_check.py`.

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
