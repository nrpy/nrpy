# Project Assembly And Generating Functions

> Explain how NRPy emits complete sibling Dendro-GR applications. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`nrpy.examples.dendro_bssn` and `nrpy.examples.dendro_fccz4` generate
`NRPy_BSSN_GR/` and `NRPy_fCCZ4_GR/` beside `BSSN_GR/`. Each directory is a
complete Dendro application, not a kernel library, mock host, or adapter around
Dendro-GR BSSN sources. With no explicit output directory, each example finds
the nearest ancestor containing `CMakeLists.txt` and `BSSN_GR/`.

## Detail

Each example registers the complete canonical state before parallel work. It
then queues every expensive Ricci, RHS, constraint, conversion, initial-data,
and runtime-service registrar for finite-difference orders 4, 6, and 8. One
`do_parallel_codegen()` call constructs and lowers all queued expressions.
The parent process does not reconstruct those expressions.

After `do_parallel_codegen()`, each example calls
`state_h.validate_registered_state` to confirm that the merged registry holds
exactly the canonical Dendro state; no two parallel tasks register the same name
with different definitions, so the merge order does not matter and repeated
generation produces byte-identical trees. Inexpensive emitters then write
headers, parameter input, context, executable entry point, checkpoint support,
prototypes, and CMake.

One Python module owns each generated numerical operation. Python basename,
registrar suffix, CFunction name, and C++ basename correspond directly. Thus
`Ricci_eval.py` registers `Ricci_eval_order_6` and emits
`Ricci_eval_order_6.cpp`; `ADM_to_BSSN.py` emits order-specific
`ADM_to_BSSN` sources; `initial_data_lambdaU.py` owns only the separate
connection-initialization pass.

`CMakeLists.py` writes an explicit source list. It includes the generated
context, entry point, checkpoint support, local TwoPunctures implementation,
runtime services, conversions, projection, and every order-specific numerical
kernel. It uses no source glob, `bssn_common`, or source from `BSSN_GR/`.
Dendrolib is the shared runtime dependency. Standalone builds default `CPU_ARCH`
to `native` and apply that architecture to the solver and fetched libraries;
`generic_avx2` selects `-mavx2 -mfma`. In-tree builds inherit Dendro-GR's
architecture setting. The examples emit intrinsic-based Ricci and RHS kernels
and package NRPy's `simd_intrinsics.h` under `generated/include`;
`register_CFunction_Ricci_eval` and `register_CFunction_rhs_eval` generate
scalar kernels only when called with `enable_intrinsics=False`. The slow-start
lapse exponential is evaluated once per block kernel call in either mode, while
its runtime `SSL_sigma` parameter remains available to the solver. `Ctx::rhs`
does not clear the unzipped RHS or Ricci buffers: both kernels write every
interior value, the RHS kernel reads Ricci only at points that the Ricci kernel
wrote with the same loop, and `Mesh::zip` reads only interior values.

Run `nrpyBssnSolver --tpid PARFILE` or `nrpyFccz4Solver --tpid PARFILE`
with one MPI task before starting a fresh evolution. This computes the
TwoPunctures spectral solution once and writes
`TPID_FILEPREFIX_nrpy_tpid_sol.bin`. Either generated solver can load that
file for the same puncture inputs. Fresh MPI evolutions read the coefficients
on every rank; they do not solve the puncture equations. A checkpoint restart
does not read the file. The reader rejects a missing, truncated, or
parameter-mismatched file. The NRPy suffix keeps these coefficients separate
from native Dendro-GR's `TPID_FILEPREFIX_tpid_sol.bin` format.

The solver reads every host parameter through a small parameter-file object that
records each key it reads; the generated binding loop assigns registered
CodeParameter keys directly, so they are excluded from the report below. An
absent key takes its default, and a present key of the wrong TOML type, such as
an integer for a real-valued parameter, stops the run with a located type
error. Every parameter is read before the TwoPunctures data are loaded or
solved. At that point rank 0 prints `<solver>: warning: parameter KEY has no
effect` for every remaining key or table member that was never read, on `--tpid`
and evolution runs alike. This covers misspelled keys, native Dendro-GR keys the
generated solver does not implement, and settings that other parameters make
inapplicable, such as the target masses when `TPID_GIVE_BARE_MASS` is 0.
`BSSN_ID_TYPE` must be 0 and `TPID_REPLACE_LAPSE_WITH_SQRT_CHI` must be true, or
the run stops at startup; because the solver always sets the initial lapse to
`sqrt(chi) = W`, native `INITIAL_LAPSE` and `TPID_INITIAL_LAPSE_PSI_EXPONENT`
have no effect and are reported as such. Each step's terminal output reduces the
maximum absolute lapse with a NaN counted as infinity, and the solver stops with
`the lapse became nonfinite` when the reduced value is not finite.

Claim evidence:
- Claim: The generated solver reads every host parameter before the TwoPunctures data are loaded, stops on a present parameter of the wrong TOML type, warns on rank 0 about every parameter-file key or table member it never read, requires `BSSN_ID_TYPE = 0` and `TPID_REPLACE_LAPSE_WITH_SQRT_CHI = true`, and stops when the per-step maximum absolute lapse, with NaN counted as infinity, is not finite.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp` (`ParameterFile`, startup checks, unread-parameter report, evolution loop); `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::terminal_output`.
- Corroboration: `nrpy/infrastructures/Dendro/CodeParameters.py`, `output_toml_bindings`; `nrpy/infrastructures/Dendro/param_toml.py`, `generate_default_parfile`; `nrpy/examples/tests/dendro_application_check.py`, `Leg.universal_checks` (S1) and `Leg.run_negatives` (N1 lapse, type, and blow-up cases; N2).

Claim evidence:
- Claim: Single-rank `--tpid` precomputes reusable spectral coefficients for both generated formulations; fresh evolution loads matching coefficients, while checkpoint restoration needs no TwoPunctures file.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`.
- Corroboration: `nrpy/examples/dendro_bssn.py` and `nrpy/examples/dendro_fccz4.py`, calls to `output_main_cpp`; `nrpy/infrastructures/Dendro/param_toml.py`, `generate_default_parfile`.

Claim evidence:
- Claim: Standalone generated builds apply the selected CPU architecture to the solver and fetched libraries, while in-tree builds inherit the parent build setting.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/CMakeLists.py`, `output_CFunctions_function_prototypes_and_construct_CMakeLists`.
- Corroboration: Dendro-GR `CMakeLists.txt`, `CPU_ARCH` selection and `add_compile_options`.

## Sources

- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - BSSN registration wave and file emission.
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - fCCZ4 registration wave and file emission.
- [parallel_codegen.py](../../../nrpy/helpers/parallel_codegen.py) - worker execution and registry merge.
- [CMakeLists.py](../../../nrpy/infrastructures/Dendro/CMakeLists.py) - prototype and explicit CMake source emission.
- [Dendro-GR CMakeLists.txt](https://github.com/paralab/Dendro-GR/blob/master/CMakeLists.txt) - parent architecture selection for in-tree builds.
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - application entry point, parameter-file reading, and startup checks.
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - `output_toml_bindings`, CodeParameter key binding.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - generated Dendro runtime context.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- Validated by: [Production Validation And Deferred Checks](validation-standalone-host-and-deferral-gates.md)
