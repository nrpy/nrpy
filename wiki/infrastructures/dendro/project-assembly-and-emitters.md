# Project Assembly And Generating Functions

> Explain how NRPy emits complete, standalone Dendrolib applications. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`nrpy.examples.dendro_bssn` and `nrpy.examples.dendro_fccz4` generate
`project/Dendro_NRPy_BSSN/` and `project/Dendro_NRPy_fCCZ4/`; `--project-dir` replaces
`project`. Each directory is a complete, standalone Dendro application, not a
kernel library, mock host, or adapter around Dendro-GR BSSN sources. It needs no
Dendro-GR checkout: its CMake project fetches Dendrolib and toml11, and every
other source is generated or packaged by NRPy. Each example ends by printing
copy-paste commands that build and run the application inside its own directory.

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
kernel. It uses no source glob, `bssn_common`, or source from `BSSN_GR/`. The
CMake project is standalone only: it fetches Dendrolib (`paralab/Dendro-5.01`, branch `master`)
and toml11 with `FetchContent`, and it is
not meant to be added to another CMake tree. `CPU_ARCH` defaults to `native` and
applies to the solver and the fetched libraries; `generic_avx2` selects `-mavx2
-mfma`. The examples emit intrinsic-based Ricci and RHS kernels and package
NRPy's `simd_intrinsics.h` under `generated/include`. The option `--fd-order`
sets only the default run-time order; see [Finite-Difference Profiles And Dendro
Conformance](finite-difference-profiles-and-dendro-conformance.md).

`register_CFunction_Ricci_eval` and `register_CFunction_rhs_eval` take a
required `enable_intrinsics` argument. `True` generates SIMD-intrinsic kernels
and is accepted only with `CoordSystem = "Cartesian"`; `False` generates
scalar kernels. Both examples pass `True`.

Claim evidence:
- Claim: `register_CFunction_Ricci_eval` and `register_CFunction_rhs_eval` take a required keyword argument `enable_intrinsics`; `True` emits SIMD-intrinsic kernels and raises `ValueError` unless `CoordSystem` is `Cartesian`, `False` emits scalar kernels, and both examples pass `True`.
- Role: descriptive behavior
- Deciding authority: [Ricci_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py), `register_CFunction_Ricci_eval`; [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), `register_CFunction_rhs_eval`
- Corroboration: [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) and [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py), `main`, `enable_intrinsics`

The slow-start
lapse exponential is evaluated once per block kernel call in either mode, while
its runtime coefficient `SSL_sigma` (key `BSSN_SSL_SIGMA`) remains available to the solver. `Ctx::rhs`
does not clear the unzipped RHS or Ricci buffers: both kernels write every
interior value, the RHS kernel reads Ricci only at points that the Ricci kernel
wrote with the same loop, and `Mesh::zip` reads only interior values.

Run `nrpyBssnSolver --tpid PARFILE` or `nrpyFccz4Solver --tpid PARFILE`
with one MPI task before starting a fresh evolution. This computes the
TwoPunctures spectral solution once and writes
`TPID_FILEPREFIX_nrpy_tpid_sol.bin`. Either generated solver can load that
file for the same puncture inputs. Fresh MPI evolutions read the coefficients
on every rank; they do not solve the puncture equations. A restore from a
checkpoint does not read the file; a restore request that finds no checkpoint
metadata starts a fresh evolution and does read it. The reader rejects a missing,
truncated, or parameter-mismatched file. The NRPy suffix keeps these
coefficients separate from native Dendro-GR's `TPID_FILEPREFIX_tpid_sol.bin`
format.

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
have no effect and are reported as such. Each step checks the
maximum absolute lapse with a NaN counted as infinity, and the solver stops with
`the lapse became nonfinite` when the reduced value is not finite.
`BSSN_TIME_STEP_OUTPUT_FREQ` controls terminal printing; disabling printing
does not disable this check. [Runtime Parameter Keys](runtime-parameters.md)
lists every key with its fallback and the values the solver rejects.

Claim evidence:
- Claim: The generated solver reads every host parameter before the TwoPunctures data are loaded, stops on a present parameter of the wrong TOML type, warns on rank 0 about every parameter-file key or table member it never read, requires `BSSN_ID_TYPE = 0` and `TPID_REPLACE_LAPSE_WITH_SQRT_CHI = true`, and stops when the per-step maximum absolute lapse, with NaN counted as infinity, is not finite.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp` (`ParameterFile`, startup checks, unread-parameter report, evolution loop); `nrpy/infrastructures/Dendro/solver_context.py`, `Ctx::terminal_output`.
- Corroboration: `nrpy/infrastructures/Dendro/CodeParameters.py`, `output_toml_bindings`; `nrpy/infrastructures/Dendro/param_toml.py`, `generate_default_parfile`; `nrpy/examples/tests/dendro_application_check.py`, `Leg.universal_checks` (S1) and `Leg.run_negatives` (N1 lapse, type, and blow-up cases; N2).

Claim evidence:
- Claim: Single-rank `--tpid` precomputes reusable spectral coefficients for both generated formulations; fresh evolution loads matching coefficients, while checkpoint restoration needs no TwoPunctures file.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`.
- Corroboration: `nrpy/examples/dendro_bssn.py` and `nrpy/examples/dendro_fccz4.py`, calls to `output_main_cpp`; `nrpy/infrastructures/Dendro/param_toml.py`, `generate_default_parfile`.

After generating, each example prints the prerequisites (CMake 3.18 or newer,
because Dendrolib's default BLAS and LAPACK path links the imported targets
`BLAS::BLAS` and `LAPACK::LAPACK`, which CMake provides from 3.18;
GNU compilers, since the solver build passes `-fext-numeric-literals`; MPI with
C, C++, and Fortran bindings; OpenMP; GSL; BLAS and LAPACK; git and network
access for Dendrolib, toml11, and spdlog) and five numbered commands: `cd` into the
application directory, configure with `cmake -S . -B build
-DCMAKE_BUILD_TYPE=Release`, build with `cmake --build build --parallel 4`, solve
the TwoPunctures data with `mpiexec -n 1 build/<executable> --tpid
pars/q1.par.lowres.toml`, and evolve with `mpiexec -n 4 build/<executable>
pars/q1.par.lowres.toml`, noting that the default parameter file is a production binary-black-hole
run whose values differ from the generated defaults (see [BSSN Application
Wiring](bssn-application-wiring.md#q1-versus-generated-defaults)). For a short test, add `BSSN_MAX_ITERATIONS = 100` to the parameter file above its first `[table]` header: a line below a `[table]` header belongs to that table and is not read as `BSSN_MAX_ITERATIONS`, and the packaged q1 file ends inside `[AEH_PARAMS]`, so a line appended at its end has no effect.
The solver's default output prefixes are relative, so a run from the application
directory writes `dat/` diagnostics, the TwoPunctures file, `vtu/`, and `cp/`
there, and `bah/` when `AEH_SOLVER_FREQ` is positive.
Both examples copy `nrpy/examples/q1.par.lowres.toml` next to their generated
`pars/bssn.toml` or `pars/fccz4.toml`. The q1 file is included in the Python
package by `setup.py`, so installed generators have the same input file.

Claim evidence:
- Claim: The generated CMake project is standalone: it declares its own project and CMake minimum 3.18, fetches Dendrolib from `master` and toml11, and applies the selected CPU architecture to the solver and the fetched libraries. The minimum matches the imported targets `BLAS::BLAS` and `LAPACK::LAPACK` that Dendrolib's default `WITH_BLAS_LAPACK` path links; the KB records no configure on an older CMake.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/CMakeLists.py`, `output_CFunctions_function_prototypes_and_construct_CMakeLists`.
- Corroboration: `nrpy/examples/tests/dendro_application_check.py`, `Leg.generate_and_build`, configures and builds each generated tree on its own; Dendrolib `CMakeLists.txt` (`paralab/Dendro-5.01`), the `WITH_BLAS_LAPACK` option and its `target_link_libraries` line.

Claim evidence:
- Claim: Each Dendro example writes its application to `<project-dir>/<SOLVER_NAME>/` (default `project`) and prints the prerequisites and the copy-paste commands that configure, build, solve the TwoPunctures data, and evolve inside that directory, followed by a short-test hint to set `BSSN_MAX_ITERATIONS = 100` above the first `[table]` header of the parameter file.
- Role: descriptive behavior
- Deciding authority: `nrpy/examples/dendro_bssn.py` and `nrpy/examples/dendro_fccz4.py`, `parse_args` and `main`; `nrpy/infrastructures/Dendro/CMakeLists.py`, `build_and_run_instructions`.
- Corroboration: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp`, relative default output prefixes and the `--tpid` requirement; `nrpy/infrastructures/Dendro/param_toml.py`, `generate_default_parfile`.

## Sources

- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - BSSN registration wave and file emission.
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - fCCZ4 registration wave and file emission.
- [parallel_codegen.py](../../../nrpy/helpers/parallel_codegen.py) - worker execution and registry merge.
- [CMakeLists.py](../../../nrpy/infrastructures/Dendro/CMakeLists.py) - prototype and explicit CMake source emission, standalone project, and build/run instructions.
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - application entry point, parameter-file reading, and startup checks.
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - `output_toml_bindings`, CodeParameter key binding.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - generated Dendro runtime context.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- Validated by: [Production Validation And Deferred Checks](production-validation-and-deferred-checks.md)
