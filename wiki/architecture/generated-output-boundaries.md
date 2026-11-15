# Generated Output Boundaries

> Boundary between checked-in NRPy generator source files and disposable generated `project/` trees, thorns, binaries, and banners. · Status: confirmed
> Up: [Architecture](index.md)

## Summary

Handwritten Python generators, symbolic equation modules, infrastructure emitters, helpers, packaging files, and CI YAML constitute NRPy's source tree, whereas generated C/CUDA/C++ applications, Einstein Toolkit/CarpetX thorns, Charm++ projects, and JAX packages written under `project/<name>/` are disposable outputs. Binaries, images, archives, and compiled objects are never committed without explicit maintainer approval and are excluded from KB content.

## Detail

### Project Output Trees And Non-Text Exclusions

Example generators write disposable trees under `project/<name>/`: standalone BHaH C/CUDA projects (`Makefile`, `.c`/`.cu`, `.h`, `.par`), Einstein Toolkit and CarpetX thorn directories linked into Cactus trees, `superB` Charm++ chare projects, `Dendro_NRPy_BSSN` / `Dendro_NRPy_fCCZ4` CMake+MPI applications (`dendro_bssn.py`, `dendro_fccz4.py`), and `sebobv1_jax` Python/JAX packages. GitHub Actions (`codegen-ubuntu`, `codegen-mac`, `einsteintoolkit-validation`, `dendro-validation`, `charmpp-validation`) generates and builds these outputs in temporary directories.

Committing binaries, images, archives, or compiled objects requires explicit maintainer approval (`coding_style.md`). Within the KB, cite the handwritten generator module and symbol rather than generated `project/` files unless a specific generated text file is registered as `frozen` in `raw/SOURCES.md` ([KB Rules](../SCHEMA.md#not-in-the-kb)).

Evidence (normative rule): `coding_style.md` (`## Python Coding Style`); `dendro_bssn.py`; `.github/workflows/main.yml`.

### Identifying Generated File Banners

Generated headers and sources carry emitter-specific banners rather than one universal template: `BHaH/BHaH_defines_h.py::output_BHaH_defines_h` writes a do-not-edit banner naming its Python module (shared by `superB_two_blackholes_collide.py`), `BHaH/CodeParameters.py::write_CodeParameters_h_files` annotates each parameter with its registering module, ETLegacy and CarpetX `interface_ccl.py::construct_interface_ccl` name NRPy in `.ccl` headers, and `Dendro/main_cpp.py::output_main_cpp` emits both a do-not-edit banner and upstream `DENDRO_LICENSE` provenance for `punctureDataPhysicalCoord`.

Evidence (descriptive behavior): `BHaH/BHaH_defines_h.py`; `BHaH/CodeParameters.py`; `Dendro/main_cpp.py`; `ETLegacy/interface_ccl.py`; `CarpetX/interface_ccl.py`.

## Sources

- [coding_style.md](../../coding_style.md) - `## Python Coding Style`, `### Formatting`
- [README.md](../../README.md) - `## What Gets Generated?`, `## Project Families and Example Generators`
- [main.yml](../../.github/workflows/main.yml) - `codegen-ubuntu`, `codegen-mac`, `einsteintoolkit-validation`, `dendro-validation`, `charmpp-validation`
- [BHaH_defines_h.py](../../nrpy/infrastructures/BHaH/BHaH_defines_h.py) - `output_BHaH_defines_h`
- [CodeParameters.py](../../nrpy/infrastructures/BHaH/CodeParameters.py) - `write_CodeParameters_h_files`
- [main_cpp.py](../../nrpy/infrastructures/Dendro/main_cpp.py) - `DENDRO_LICENSE`, `output_main_cpp`
- [ETLegacy interface_ccl.py](../../nrpy/infrastructures/ETLegacy/interface_ccl.py) - `construct_interface_ccl`
- [CarpetX interface_ccl.py](../../nrpy/infrastructures/CarpetX/interface_ccl.py) - `construct_interface_ccl`
- [superB_two_blackholes_collide.py](../../nrpy/examples/superB_two_blackholes_collide.py) - shared `output_BHaH_defines_h` call
- [dendro_bssn.py](../../nrpy/examples/dendro_bssn.py) - `Dendro_NRPy_BSSN` output tree
- [dendro_fccz4.py](../../nrpy/examples/dendro_fccz4.py) - `Dendro_NRPy_fCCZ4` output tree

## See Also

- Depends on: [Build And Run](build-and-run.md)
- Depends on: [Contribution Style And Static Analysis](contribution-style-and-static-analysis.md)
- See also: [Overview](overview.md)
