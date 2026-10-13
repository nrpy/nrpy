# Generated Output Boundaries

> Distinguish handwritten source from generated projects, thorns, and other output files. · Status: confirmed
> Up: [Architecture](index.md)

## Summary

NRPy source files are the handwritten generators, symbolic equations, infrastructure modules, helpers, packaging files, and CI configuration. Generated C/CUDA projects, Einstein Toolkit thorns, Charm++ projects, and JAX projects are outputs of those generators. Selected generated text may be registered as frozen evidence; binaries, archives, images, and scratch outputs remain excluded from KB content.

## Detail

Most example generators write under `project/<name>/`. Standalone BHaH examples usually create a buildable project containing C source, headers, a Makefile, a parameter file, and an executable target for Cartesian or curvilinear-coordinate workflows. Waveform, geodesic, and some physics generators also produce standalone projects, often with a GSL build dependency.

Einstein Toolkit and CarpetX generators produce thorn directories rather than standalone executables. Those thorns are generated into `project/<name>/` and then copied or linked into an Einstein Toolkit checkout for build and run. `superB` generators produce Charm++ projects, `dendro_bssn.py` and `dendro_fccz4.py` produce standalone `Dendro_NRPy_BSSN` and `Dendro_NRPy_fCCZ4` CMake and MPI applications, and `sebobv1_jax` generates a Python/JAX project instead of a C executable.

Claim evidence:
- Claim: `dendro_bssn.py` and `dendro_fccz4.py` generate standalone `Dendro_NRPy_BSSN` and `Dendro_NRPy_fCCZ4` CMake and MPI applications that are built against Dendrolib.
- Role: descriptive behavior
- Deciding authority: [dendro_bssn.py](../../nrpy/examples/dendro_bssn.py) and [dendro_fccz4.py](../../nrpy/examples/dendro_fccz4.py), `main`
- Corroboration: [main.yml](../../.github/workflows/main.yml), `dendro-validation`, which generates, builds, and runs both applications

CI deliberately creates and builds generated outputs as validation evidence. The Ubuntu and macOS codegen jobs install NRPy, generate many projects, and run `make` where applicable. Separate Einstein Toolkit and Charm++ jobs generate backend products inside CI containers, the Dendro job generates them in the runner's temporary directory, and each job runs its configured build, runtime, and required checks.

Non-text files such as binaries, images, archives, and compiled outputs require
maintainer approval before they can be committed. Registering a file as evidence
does not grant that approval. The KB permits selected generated text files to be
registered as frozen evidence, while its exclusions still prohibit binaries,
images, archives, and compiled outputs. Cite the generator and stable symbols
unless a selected generated text file supplies necessary evidence.

Claim evidence:
- Claim: Registering generated text as frozen evidence does not replace maintainer approval for non-text contributions or override the KB's excluded file types.
- Role: normative rule
- Deciding authority: [coding_style.md](../../coding_style.md), `## Python Coding Style`, `### Formatting`; [KB Page Format](../SCHEMA.md), `## Not In The KB`
- Corroboration: [original-agents.md](../../raw/source-docs/original-agents.md), `## Scope`, also requires maintainer approval for non-text additions.

### Identifying Generated Files

Generated-file comments depend on the writing function. `BHaH_defines.h` names
its generating module in a do-not-edit banner. ETLegacy and CarpetX `.ccl`
writers name NRPy. The superB black-hole collision generator also writes
`BHaH_defines.h` through the shared BHaH writer, so that file carries the same
banner in the superB project. BHaH parameter declaration files instead include comments
identifying each parameter's defining module. Dendro's generated main program
includes the upstream source and function name, adaptation notes, copyright,
and license text for its imported puncture initial-data routine, as well as a
do-not-edit banner. These examples do not establish one required banner format
for every generated file.

Claim evidence:
- Claim: Generated-file identification is writer-specific: BHaH definition headers and CCL files carry generation banners, BHaH parameter declarations identify defining modules, and Dendro's main program also preserves upstream attribution and license text.
- Role: descriptive behavior
- Deciding authority: [BHaH_defines_h.py](../../nrpy/infrastructures/BHaH/BHaH_defines_h.py), `output_BHaH_defines_h`; [CodeParameters.py](../../nrpy/infrastructures/BHaH/CodeParameters.py), `write_CodeParameters_h_files`; [main_cpp.py](../../nrpy/infrastructures/Dendro/main_cpp.py), `DENDRO_LICENSE` and `output_main_cpp`
- Corroboration: [superB_two_blackholes_collide.py](../../nrpy/examples/superB_two_blackholes_collide.py) calls the shared BHaH header writer; [ETLegacy interface_ccl.py](../../nrpy/infrastructures/ETLegacy/interface_ccl.py) and [CarpetX interface_ccl.py](../../nrpy/infrastructures/CarpetX/interface_ccl.py), `construct_interface_ccl`, emit the NRPy generation banner.

## Sources

- [superB_two_blackholes_collide.py](../../nrpy/examples/superB_two_blackholes_collide.py) - `BHaH.BHaH_defines_h.output_BHaH_defines_h` call

- [coding_style.md](../../coding_style.md) - `## Python Coding Style`, `### Formatting`
- [BHaH_defines_h.py](../../nrpy/infrastructures/BHaH/BHaH_defines_h.py) - `output_BHaH_defines_h`
- [CodeParameters.py](../../nrpy/infrastructures/BHaH/CodeParameters.py) - `write_CodeParameters_h_files`
- [main_cpp.py](../../nrpy/infrastructures/Dendro/main_cpp.py) - `DENDRO_LICENSE`, `output_main_cpp`
- [ETLegacy interface_ccl.py](../../nrpy/infrastructures/ETLegacy/interface_ccl.py) and [CarpetX interface_ccl.py](../../nrpy/infrastructures/CarpetX/interface_ccl.py) - `construct_interface_ccl`

- [README.md](../../README.md) - `## What Gets Generated?`, `## Project Families and Example Generators`
- [.github/workflows/main.yml](../../.github/workflows/main.yml) - `einsteintoolkit-validation`, `dendro-validation`, `charmpp-validation`, `codegen-ubuntu`, `codegen-mac`
- [dendro_bssn.py](../../nrpy/examples/dendro_bssn.py) and [dendro_fccz4.py](../../nrpy/examples/dendro_fccz4.py) - the Dendro application output directories.
- [raw/source-docs/original-agents.md](../../raw/source-docs/original-agents.md) - `## Scope`

## See Also

- Parent: [Architecture](index.md)
- Depends on: [Build And Run](build-and-run.md)
- Depends on: [Contribution Style And Static Analysis](contribution-style-and-static-analysis.md)
- See also: [Overview](overview.md)
