# Einstein Toolkit Thorn Generators

> Explain NRPy's Einstein Toolkit thorn-generation examples and checked-in ET fixtures. · Status: confirmed · Last reconciled: 2026-06-30
> Up: [Examples](index.md)

## Summary

NRPy has four Einstein Toolkit thorn generators. The ETLegacy/Carpet examples
write classic Carpet thorns: `carpet_wavetoy_thorns.py` writes
`WaveToyNRPy`, `IDWaveToyNRPy`, and `diagWaveToyNRPy` under
`project/et_wavetoy/`, while `carpet_baikal_thorns.py` writes `Baikal` and
`BaikalVacuum` under `project/et_baikal/`. The CarpetX examples write CarpetX
counterparts: `carpetx_wavetoy_thorns.py` writes `WaveToyNRPyX`,
`IDWaveToyNRPyX`, and `diagWaveToyNRPyX` under `project/et_wavetoy/`, while
`carpetx_baikal_thorns.py` writes `BaikalX` and `BaikalVacuumX` under
`project/et_baikalx/`.

These generators are Python generation workflows, not standalone executables.
They create thorn directories and CCL/build metadata that must be copied or
linked into an Einstein Toolkit checkout before ET build or testsuite runs.

## Detail

The command shape follows the module names from the repository root:

```bash
python -m nrpy.examples.carpet_wavetoy_thorns
python -m nrpy.examples.carpet_baikal_thorns
python -m nrpy.examples.carpetx_wavetoy_thorns
python -m nrpy.examples.carpetx_baikal_thorns
```

Per-generator map:

| Generator | Output family | Thorn names | Validation route |
| --- | --- | --- | --- |
| `carpet_wavetoy_thorns.py` | ETLegacy/Carpet WaveToy thorns in `project/et_wavetoy/` | `WaveToyNRPy`, `IDWaveToyNRPy`, `diagWaveToyNRPy` | GitHub Actions `einsteintoolkit-validation` generates, links into ET, builds ET, and runs `WaveToyNRPy` testsuite. |
| `carpet_baikal_thorns.py` | ETLegacy/Carpet Baikal thorns in `project/et_baikal/` | `Baikal`, `BaikalVacuum` | GitHub Actions `einsteintoolkit-validation` generates, links into ET, builds ET, and runs `Baikal` and `BaikalVacuum` testsuites. |
| `carpetx_wavetoy_thorns.py` | CarpetX WaveToy thorns in `project/et_wavetoy/` | `WaveToyNRPyX`, `IDWaveToyNRPyX`, `diagWaveToyNRPyX` | Local full-CI script runs the generator and skips ET compile for `carpet*` scripts. |
| `carpetx_baikal_thorns.py` | CarpetX Baikal thorns in `project/et_baikalx/` | `BaikalX`, `BaikalVacuumX` | Local full-CI script runs the generator and skips ET compile for `carpet*` scripts. |

The ETLegacy/Carpet generators set `Infrastructure` to `ETLegacy`. Their
generated C functions target Cactus/Carpet schedule and storage conventions:
WaveToy RHS work is scheduled in `MoL_CalcRHS`, the evolution thorn inherits
`Boundary Grid MethodofLines`, and the Baikal thorns register Method of Lines,
symmetry, boundary, RHS-zeroing, BSSN conversion, Ricci, RHS, constraints,
slicing, lapse-floor, and determinant-enforcement functions. `Baikal` enables
stress-energy support through `T4DD_to_T4UU` and finite-difference orders
`2` and `4`; `BaikalVacuum` drops `T4UU` gridfunctions and emits vacuum
variants for finite-difference orders `4`, `6`, and `8`.

The CarpetX generators set `Infrastructure` to `CarpetX`. Their generated
functions target CarpetX loop and schedule conventions: WaveToy RHS work is
scheduled in `ODESolvers_RHS`, CarpetX functions include `loop_device.hxx`,
CarpetX Baikal writes `configuration.ccl`, and generated sources are C++ thorn
sources rather than ETLegacy C thorn sources. `BaikalX` is the matter-capable
thorn and `BaikalVacuumX` is the vacuum thorn; the finite-difference split
mirrors the Carpet Baikal generator. The CarpetX Baikal generator uses
`ADMBaseX`, `loop_device.hxx`, `REQUIRES Loop CarpetX`, CarpetX boundary and
zero-RHS helpers, and the same BSSN/Ricci/RHS/constraints/poststep-repair
registration family as the ETLegacy generator.

WaveToy is split into three thorns in both output families. The ID thorn writes
the spherical-Gaussian exact solution at `CCTK_INITIAL` and owns initial-data
parameters such as `sigma` and `wavespeed`. The evolution thorn owns `uu` and
`vv` evolution variables, registers RHS evaluation for the wave equation, and
shares `wavespeed` from the ID thorn. The diagnostics thorn writes exact
solution auxiliary variables during analysis and shares or reads the ID
parameters needed for comparison output.

Each generator creates CCL files and thorn-local build metadata through the
infrastructure writers, then writes registered C or C++ functions through
`make.code.defn`. These generated `interface.ccl`, `param.ccl`,
`schedule.ccl`, optional `configuration.ccl`, `src/make.code.defn`, and
`src/*` files are generated thorn products under `project/**`, not durable KB
sources. The checked-in fixture files under `nrpy/examples/et_WaveToyfiles/`
are different: `ThornList` is an ET component list including Cactus, Carpet,
Carpet IO, CarpetRegrid2, MoL, NewRad, and other thorns needed by the Carpet
WaveToy run; `WaveToyNRPy.par` is a larger Carpet WaveToy parfile; the
`test/` subtree contains a testsuite descriptor, a small eight-iteration
`WaveToyNRPy_test.par`, and `WaveToyNRPy_test/uuGF.x.asc` as checked-in fixture
evidence for the CarpetIOASCII `uuGF` x-line output. That ASCII file is not a
generated thorn product; it is registered directly in `raw/SOURCES.md`.

Python generation dependencies and ET build dependencies are separate. Running
the generator needs the NRPy Python environment and its symbolic/codegen stack.
Building or running the resulting thorns needs an Einstein Toolkit checkout and
external ET components. The repository README states this boundary directly:
the ET/CarpetX generators write thorns under `project/<name>/`, and an ET
checkout is only needed when compiling or running those thorns.

Validation is also split by family. GitHub Actions
`einsteintoolkit-validation` runs the Carpet WaveToy and Baikal generators in
an Apptainer ET image, symlinks `project/et_baikal/Baikal*` and
`project/et_wavetoy/*` into the ET arrangements tree, symlinks the checked-in
WaveToy tests into `WaveToyNRPy/`, builds ET, and runs the `Baikal`,
`BaikalVacuum`, and `WaveToyNRPy` testsuites. The local full-CI script also
runs all four generator scripts and deliberately skips compiling paths whose
script name contains `carpet`, so the CarpetX examples have a checked local
generation route but no local ET build route in that script.

## Sources

- [carpet_wavetoy_thorns.py](../../nrpy/examples/carpet_wavetoy_thorns.py) - `project_name`, `ID_thorn_name`, `diag_thorn_name`, `evol_thorn_name`, `register_CFunction_rhs_eval`
- [carpet_wavetoy_thorns.py](../../nrpy/examples/carpet_wavetoy_thorns.py) - `construct_schedule_ccl`, `construct_interface_ccl`, `construct_param_ccl`
- [carpet_baikal_thorns.py](../../nrpy/examples/carpet_baikal_thorns.py) - `project_name`, `thorn_names`, `enable_T4munu`, `fd_order_list`
- [carpet_baikal_thorns.py](../../nrpy/examples/carpet_baikal_thorns.py) - `register_CFunction_ADM_to_BSSN`, `register_CFunction_Ricci_eval`, `register_CFunction_rhs_eval`, `register_CFunction_BSSN_constraints`, `output_CFunctions_and_construct_make_code_defn`
- [carpetx_wavetoy_thorns.py](../../nrpy/examples/carpetx_wavetoy_thorns.py) - `project_name`, `ID_thorn_name`, `diag_thorn_name`, `evol_thorn_name`, `register_CFunction_rhs_eval`
- [carpetx_wavetoy_thorns.py](../../nrpy/examples/carpetx_wavetoy_thorns.py) - `ODESolvers_RHS`, `loop_device.hxx`, `DECLARE_CCTK_ARGUMENTSX`
- [carpetx_baikal_thorns.py](../../nrpy/examples/carpetx_baikal_thorns.py) - `project_name`, `thorn_names`, `configuration_ccl.construct_configuration_ccl`
- [carpetx_baikal_thorns.py](../../nrpy/examples/carpetx_baikal_thorns.py) - `ADMBaseX`, `USES INCLUDE: loop_device.hxx`, `register_CFunction_ADM_to_BSSN`, `register_CFunction_Ricci_eval`, `register_CFunction_rhs_eval`, `register_CFunction_BSSN_constraints`
- [ThornList](../../nrpy/examples/et_WaveToyfiles/ThornList) - Einstein Toolkit component list fixture
- [WaveToyNRPy.par](../../nrpy/examples/et_WaveToyfiles/WaveToyNRPy.par) - Carpet WaveToy parfile fixture
- [test.ccl](../../nrpy/examples/et_WaveToyfiles/test/test.ccl) - `TEST WaveToyNRPy_test`
- [WaveToyNRPy_test.par](../../nrpy/examples/et_WaveToyfiles/test/WaveToyNRPy_test.par) - eight-iteration Carpet WaveToy test parfile
- [uuGF.x.asc](../../nrpy/examples/et_WaveToyfiles/test/WaveToyNRPy_test/uuGF.x.asc) - checked-in fixture evidence for CarpetIOASCII `uuGF` x-line output
- [README.md](../../README.md) - `Einstein Toolkit and CarpetX Generators`
- [.github/workflows/main.yml](../../.github/workflows/main.yml) - `einsteintoolkit-validation`
- [.github/full_nrpy_local_ci.sh](../../.github/full_nrpy_local_ci.sh) - `example_scripts`

## See Also

- [Examples](index.md)
- [Generated Output Boundaries](../architecture/generated-output-boundaries.md)
- [Generated Project CI](../validation/generated-project-ci.md)
- [ETLegacy Thorn Assembly And CCL Files](../infrastructures/etlegacy/thorn-assembly-and-ccl-files.md)
- [CarpetX Thorn Assembly, Configuration, And CCL Files](../infrastructures/carpetx/thorn-assembly-configuration-and-ccl-files.md)
