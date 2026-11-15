# Infrastructure Code Style

> `register_CFunction_*` structure, `CodeParameter` registration scope, parallel-codegen guards, `validate_strings` size limits, and BHaH loop/codegen conventions. · Status: provisional
> Up: [Infrastructures](index.md)

## Summary

Infrastructure modules under `nrpy/infrastructures/` center on `register_CFunction_<name>()` functions written in generated execution order: they register `CodeParameter` entries inside one setup function (never at import time), use `pcg_registration_phase()` guards when parallel codegen applies, limit `validate_strings()` baselines to small structural C/C++/CUDA files rather than large symbolic kernels, and prefer `simple_loop()` + `c_codegen(..., automatically_read_gf_data_from_memory=True)` over handwritten memory-load boilerplate.

## Detail

### Existing-Backend Conformance, Module Layout, And `CodeParameter` Scope

New infrastructure code conforms one-way to patterns already established in `>= 3` sites across `nrpy/infrastructures/{BHaH,ETLegacy,CarpetX,superB}` ([New Infrastructure Conformance](new-infrastructure-conformance.md)), following [Python Coding Style](../architecture/python-coding-style.md) and [Code Test Policy](../validation/code-test-policy.md). Each module normally defines `register_CFunction_<name>()` declaring `subdirectory`, `includes`, optional `prefunc` (only for C helpers called `>= 2` times or exceeding `~10–15` logic lines), `desc`, `cfunc_type`, `name`, `params`, `include_CodeParameters_h`, `body`, and optional `postfunc` linearly before `cfc.register_CFunction()`, using `# fmt: off`/`# fmt: on` sparingly (e.g., aligned `CodeParameter` tables).

Because `CodeParameter` registration mutates global header/parfile registries, never register `CodeParameter` objects at import time or redundantly across modules: register each parameter in one central setup/registration function and document downstream dependencies. Functions participating in `nrpy.helpers.parallel_codegen` return `Union[None, pcg.NRPyEnv_type]` and guard entry with `if pcg.pcg_registration_phase(): pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals()); return None` (non-PCG registrars return `None`).

### `Doctests:`, `validate_strings()` Size Rule, And BHaH Codegen Idioms

Place `Doctests:` immediately before `>>>` prompts. Reserve `clang_format()` + `validate_strings()` baselines for small (`~1–80` KB), predominantly structural C/C++/CUDA functions, headers, or parfiles; never golden-test symbolically lowered RHS, Ricci, constraint, initial-data conversion (`Dendro/general_relativity/ADM_to_BSSN.py`), connection, or algebraic-projection kernels, which are validated via `nrpy/equations/` `trusted_dict` suites and compiled CI runs ([Test Oracles And Safe Updates](../validation/test-oracles-and-safe-updates.md)).

In new BHaH per-grid/per-point kernels:
- Use `BHaH.simple_loop.simple_loop()` and `ccg.c_codegen(..., automatically_read_gf_data_from_memory=True)` (`in_gfs`, etc.) rather than handwritten C loops, string replacements on generated C, post-registration `CFunction_dict` mutation, or body-internal `#include "set_CodeParameters.h"` (use `get_params_commondata_symbols_from_expr_list()` + `generate_definition_header()` when needed), assembling `body += ...` top-to-bottom.
- Standard C structs and macros: `commondata_struct *restrict commondata`, `griddata_struct *restrict griddata`, `params_struct *restrict params`, `bc_struct *restrict bcstruct`; groups `EVOL`, `AUXEVOL`, `AUX`, `DIAG` with uppercase `*GF` tokens (`TRKGF`, `CFGF`, `HDD00GF`, `ADD00GF`); `IDX3`, `IDX4`, `IDX4pt`, `IDX2` (plus local `SRC_IDX4`/`DST_IDX4`); and BHaHAHA-style enum return codes (`BHAHAHA_SUCCESS`, `bah_error_message()`, `#pragma omp critical` shared `error_flag`).

Evidence (normative rule): `coding_style.md` (`## Infrastructure Code Patterns`, `### CodeParameter Registration Scope`, `#### validate_strings pattern`); `Dendro/general_relativity/ADM_to_BSSN.py`.

## Sources

- [coding_style.md](../../coding_style.md) - `## Infrastructure Code Patterns`, `### CodeParameter Registration Scope`, `#### validate_strings pattern`
- [ADM_to_BSSN.py](../../nrpy/infrastructures/Dendro/general_relativity/ADM_to_BSSN.py) - symbolic pointwise kernel without generated-source baseline

## See Also

- Depends on: [New Infrastructure Conformance](new-infrastructure-conformance.md)
- Depends on: [C Function Registry](../core/c-function-registry.md)
- Depends on: [Gridfunctions And Parameters](../core/gridfunctions-and-parameters.md)
