# Production Validation And Deferred Checks

> Generation guards, `dendro_application_check.py` MPI profiles and stored-reference checks, native comparison prerequisites, and unverified configurations for NRPy's Dendro solvers. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`nrpy.examples.dendro_bssn` and `nrpy.examples.dendro_fccz4` combine generation-time `ValueError` checks, `static-analysis` linting on generator modules, and the `dendro-validation` / weekly `dendrolib-canary.yml` MPI test script `nrpy/examples/tests/dendro_application_check.py`, supplemented by the manual `dendrolib_capability_test.cpp` host-layout check.

## Detail

### Generation-Time Guards And Native Comparison Prerequisites

Generator modules (`dendro_bssn.py`, `dendro_fccz4.py`) are checked by `static-analysis` (which skips `*/tests/*`, so `dendro_application_check.py` is excluded) and enforce generation-time `ValueError` guards:
- `output_solver_context_cpp` requires all order-specific kernels plus `BSSN_to_ADM`, `floor_the_lapse_and_conformal_factor`, `gravitational_waves`, `adm_quantities`, and `physical_boundary_ghosts` (`twopunctures`, `physical_boundary`, `diagnostics`, `apparent_horizon`, and `enforce_detgbar_equals_detghat_trAzero` have no registration guard).
- `register_CFunction_rhs_eval` enforces canonical RHS key order, `validate_registered_state` enforces canonical `EVOL`/`AUXEVOL` state order, and `output_CFunctions_function_prototypes_and_construct_CMakeLists` rejects unsafe, duplicate, or native `BSSN_GR` paths among application-supplied CMake sources (module-basename-to-registrar matching and raw pointer indexing beyond RHS key order remain review checks).

When comparing generated runs against native `Dendro-GR` (`BSSN_GR`), sharing a TOML file does not isolate the continuum equations: record evolved conformal factor (`W` vs `chi`), full-`psi` lapse replacement, constant `eta` vs radial `RIT_ETA_*`, KO strength and native startup/post-merger `CAKO`, single-step vs Heun puncture tracking (which also alters Nyquist and excision histories), RK4/`min(dx,dy,dz)` vs `BSSN_RK_TYPE`/X-width CFL, omitted-key defaults, lower-index `MD0..MD2` conventions, and equal-weight unexcised unique-node `*_Constraints.dat` (`*_Constraints_volweighted.dat` is distinct) at matched physical times and post-remesh grids.

Evidence (descriptive behavior): `Dendro/solver_context.py::output_solver_context_cpp`; `Dendro/general_relativity/rhs_eval.py`; `Dendro/state_h.py::validate_registered_state`; `Dendro/CMakeLists.py`.

### Automated `dendro-validation` Profiles And CLI Execution

Run `python nrpy/examples/tests/dendro_application_check.py --formulation {bssn,fccz4} --work-dir DIR` (nonexistent `DIR`, deleted on exit; defaults `--launcher "mpiexec --oversubscribe --bind-to none"`, `--ranks 4` `>= 4`, `--build-jobs 4`; running with no arguments executes doctests only). Both `main.yml` (`dendro-validation`, `OMP_NUM_THREADS=1`, 75-minute timeout) and weekly/manual `dendrolib-canary.yml` (`--dendrolib-ref master` shallow-cloning `FETCHCONTENT_SOURCE_DIR_DENDROLIB`) build `W` and `chi` variants with `-DCPU_ARCH=x86-64-v3` after asserting byte-identical repeat generation, solve `--tpid` once per variant, and run:
- **Profile P** (`MAXDEPTH = 13`, FD6, `8` steps, cadence `4`, `BSSN_SCALE_VTU_AND_GW_EXTRACTION = false`): checks finite-radius closed-form `E` and `J_z` plus vanishing `P_i` and `J_{x,y}`, horizon irreducible masses vs puncture ADM masses and reflection symmetry, vanishing odd-`m` modes and `C(l,-m) = (-1)^l conj C(l,m)` across 21 modes (`l = 2..4`, 2 radii), step-4 puncture-center reflection and motion along momenta, byte-identical `dat/`/`bah/`/`vtu/` output across step-4 checkpoint/restore (normalizing `GridInfo` wall time and skipping `dgr__PARAM_DUMP__*.toml`), 3-rank vs 4-rank diagnostic and unexcised node-count agreement (horizon observables through irreducible mass, excluding solver residuals), step-0 `W`-vs-`chi` agreement (`nodes`, `E`, `J_z`, `Psi4`), and Run A step `0, 4, 8` regression against `dendro_application_check_reference.py` (`--update-reference` requires `black`, rejects `--dendrolib-ref`, and updates sequentially per formulation).
- **Profile O** (`MAXDEPTH = 10`, `4` steps): checks initial Hamiltonian-constraint norm halving from FD4 to FD6 to FD8 on an equal-element-count mesh plus a 1-rank repeat, alongside per-run sanity checks (finite rows, expected cadence, node ceiling, zero unread-key warnings).
- **Negative and warning checks**: requires exit `1..123` with expected diagnostics for `W`/`chi` cross-restore, multi-rank `--tpid`, TwoPunctures parameter mismatch or missing file, invalid `BSSN_ELE_ORDER` or `BSSN_REFINEMENT_MODE`, blow-up CFL, nonfinite `BSSN_CFL_FACTOR`/`BSSN_SSL_SIGMA`/`ETA_CONST`, directory-blocked `*_ADM.dat`, `TPID_REPLACE_LAPSE_WITH_SQRT_CHI = false`, integer for float TOML key, and lapse blow-up with constraint output disabled; also verifies unread-key warnings and fresh fallback on `BSSN_RESTORE_SOLVER = 1` without metadata.

Evidence (CI behavior): `.github/workflows/main.yml`; `.github/workflows/dendrolib-canary.yml`; `dendro_application_check.py`; `dendro_application_check_reference.py`.

### Verification Limits And Manual `dendrolib_capability_test.cpp`

CI tolerances in `dendro_application_check.py` (`REFERENCE_RTOL`, `REFERENCE_ATOL`) are literal regression thresholds tied to printed precision (`10` significant digits on constraints/ADM/horizon radii, `16` on area/mass, `3` and `6` fixed decimals on time and centroid where a 1-ULP boundary flip can exceed tolerance; roundoff-zero `P_i`, `J_{x,y}`, and step-0 `Lambda` plus checkpoint write/restore `1e-10` algebraic projection act only as upper bounds). Rank comparisons of wave/GridInfo tables match common steps pairwise, FD-order mesh sharing checks element count rather than octant coordinates, `W`-vs-`chi` omits constraints/horizons, pre-creating `dat/`/`vtu/`/`cp/`/`bah/` skips solver directory creation, `-DCPU_ARCH=x86-64-v3` skips `native`, and `BSSN_SCALE_VTU_AND_GW_EXTRACTION = false` skips native cadence scaling.

Unexercised by CI: `pars/q1.par.lowres.toml` and printed CLI commands, generation with `--fd-order {4,8}`, `enable_KreissOliger_dissipation=False`, `--ybs-gamma`, or `--ybs-momentum`, 3-rank horizon-JSON restore contents, pointwise `alpha = W = sqrt(chi)` field inspection, intra-rank block-boundary isolation, ASan/UBSan, and long-time or merger evolutions.

`nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp` (documented in its sibling `README.md`, not built in CI) manually tests a Dendrolib build's scalar ABI, padded block dimensions/origin, unzip offsets, variable-major `x`-fastest layout, and halo validity (`CAPTEST_ORDERS` selects even orders; `CAPTEST_INJECT` injects zero/shift defects; note `err > 1e-9` evaluates `false` on `NaN`, so `NaN` host writes would pass).

Evidence (descriptive behavior): `dendro_application_check.py`; `Dendro/tests_infra/dendrolib_capability_test.cpp`; `Dendro/tests_infra/README.md`; `Dendro/checkpoint.py`.

## Sources

- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - BSSN generator
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - fCCZ4 generator
- [CMakeLists.py](../../../nrpy/infrastructures/Dendro/CMakeLists.py) - CMake manifest and `BSSN_GR` source screen
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - canonical RHS key-order guard
- [state_h.py](../../../nrpy/infrastructures/Dendro/state_h.py) - `validate_registered_state`
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - `output_solver_context_cpp` kernel guards and diagnostic precision
- [BSSN_constraints.py](../../../nrpy/infrastructures/Dendro/general_relativity/BSSN_constraints.py) - lower-index `MD0..MD2`
- [checkpoint.py](../../../nrpy/infrastructures/Dendro/checkpoint.py) - formulation tags and `1e-10` projection-residual checks
- [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py) - Profiles P and O, reference/rank/restart/negative checks
- [dendro_application_check_reference.py](../../../nrpy/examples/tests/dendro_application_check_reference.py) - `trusted_dict` reference values
- [main.yml](../../../.github/workflows/main.yml) - `static-analysis` and `dendro-validation`
- [dendrolib-canary.yml](../../../.github/workflows/dendrolib-canary.yml) - weekly `dendro-validation-dendrolib-master`
- [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp) - manual Dendrolib ABI/layout checker
- [README.md](../../../nrpy/infrastructures/Dendro/tests_infra/README.md) - `CAPTEST_ORDERS` and `CAPTEST_INJECT` instructions
- [Dendro-GR](https://github.com/paralab/Dendro-GR) - native `parameters.cpp`, `bssngr_main.cpp`, `grUtils.cpp`, `rhs.cpp`

## See Also

- Depends on: [Generated Project CI](../../validation/generated-project-ci.md)
- Depends on: [Code Test Policy](../../validation/code-test-policy.md)
- See also: [Project Assembly And Generating Functions](project-assembly-and-emitters.md)
