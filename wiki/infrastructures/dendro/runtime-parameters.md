# Runtime Parameter Keys

> TOML parameter keys, call-site fallbacks, `pars/<stem>.toml` vs `pars/q1.par.lowres.toml` values, `CodeParameter` bindings, and startup validation in generated Dendro solvers. · Status: confirmed
> Up: [Dendro](index.md)

## Summary

`Dendro_NRPy_BSSN` and `Dendro_NRPy_fCCZ4` read one `[--tpid] PARAM_FILE` TOML configuration. Omitted keys use call-site fallbacks, wrong TOML value types abort with a located error, and unread keys or table members print `warning: parameter KEY has no effect` on rank `0`. Choose intentionally between minimal generated `pars/<stem>.toml` and full BBH `pars/q1.par.lowres.toml` ([q1 versus generated defaults](bssn-application-wiring.md#q1-versus-generated-defaults)).

## Detail

### Key Resolution, Run Control, Domain, And Mesh Refinement

Host keys are read from the TOML root unless qualified by `[table]` (lines placed after `[BSSN_BH1]`, `[BSSN_BH2]`, or `[AEH_PARAMS]` belong to that table and do not set root keys such as `BSSN_MAX_ITERATIONS`). Columns `generated` and `q1` show values in `pars/<stem>.toml` and `pars/q1.par.lowres.toml` (`none` means omitted so the fallback applies).

| Key | Fallback | generated | q1 | Meaning and accepted values |
| --- | --- | --- | --- | --- |
| `BSSN_ID_TYPE` | `0` | `0` | `0` | Must be `0` (TwoPunctures). |
| `BSSN_RESTORE_SOLVER` | `0` | `0` | `0` | Nonzero requests checkpoint restore ([Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md)). |
| `BSSN_ELE_ORDER` | `--fd-order` (`6`) | `--fd-order` | `6` | Element/FD order (`4`, `6`, or `8`); must match checkpoint on restore ([Finite-Difference Profiles](finite-difference-profiles-and-dendro-conformance.md)). |
| `BSSN_RK_TIME_BEGIN`, `BSSN_RK_TIME_END` | `0.0`, `700.0` | `0.0`, `1000000.0` | `0.0`, `1000000.0` | Start and end time (`END > BEGIN`). |
| `BSSN_MAX_ITERATIONS` | unlimited | `none` | `none` | Root-level step cap. |
| `BSSN_CFL_FACTOR` | `0.25` | `none` (`0.25` w/ `--ybs-momentum`) | `0.25` | Positive finite factor multiplying finest `min(dx, dy, dz)` (also `CFL_FACTOR` for `--ybs-momentum`). |
| `BSSN_TIME_STEP_OUTPUT_FREQ` | `80` | `80` | `80` | Terminal print cadence in steps (per-step lapse check always runs). |
| `BSSN_GRID_MIN_{X,Y,Z}`, `BSSN_GRID_MAX_{X,Y,Z}` | `-400.0`, `400.0` | `none` | `-400.0`, `400.0` | Domain bounds (`MAX > MIN`, finite widths). |
| `BSSN_MINDEPTH`, `BSSN_MAXDEPTH` | `4`, `14` | `none` | `4`, `14` | Octree levels (`MIN <= MAX < 31`). |
| `BSSN_REFINEMENT_MODE` | `4` | `none` | `4` | Must be `4` (`BH_WAMR`). |
| `BSSN_USE_WAVELET_TOL_FUNCTION` | `0` | `none` | `6` | `0` constant `BSSN_WAVELET_TOL`; `6` causal radial profile ([Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md)). |
| `BSSN_WAVELET_TOL`, `_MAX`, `BSSN_GW_REFINE_WTOL` | `1e-5`, `WTOL`, `WTOL` | `none` | `1e-5`, `1e-3`, `1e-4` | Base, mode-6 outer/pre-causal max, and extraction-zone wavelet tolerances (positive finite). |
| `BSSN_DENDRO_AMR_FAC`, `_POST_MERGER` | `0.1`, `0.0` | `none` | `0.2`, `0.03125` | Coarsening factor in `(0, 1]` and post-merger-checkpoint override in `[0, 1]` (`0` keeps base). |
| `BSSN_REFINE_VARIABLE_INDICES`, `BSSN_NUM_REFINE_VARS` | all evolved, list len | `none` | `0..23`, `24` | Refinement field indices and active prefix count (`1 .. len`, valid evolved indices). |
| `BSSN_REMESH_TEST_FREQ`, `_AFTER_MERGER` | `50`, `10` | `none` | `50`, `none` | Pre- and post-merger remesh test cadences in steps. |
| `BSSN_INIT_GRID_ITER`, `BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE` | `10`, `true` | `none` | `10`, `true` | Initial refinement cap (`<= N - 1` passes for `N > 1`); boolean must be `true` when `ITER > 0`. |
| `BSSN_DENDRO_GRAIN_SZ`, `BSSN_LOAD_IMB_TOL`, `BSSN_SPLIT_FIX` | `1000`, `0.1`, `2` | `none` | `50`, `0.1`, `none` | Partition grain (`> 0`), load-imbalance tolerance (`>= 0`, finite), and splitter arg. |
| `BSSN_BH{1,2}_AMR_R`, `BSSN_AMR_R_RATIO` | `2.0`, `2.0` | `none` | `1.0`, `1.618033988749895` | Puncture inner floor/tolerance radius (`> 0`) and geometric shell ratio (`> 1`). |
| `BSSN_BH{1,2}_MAX_LEV` | `BSSN_MAXDEPTH` | `none` | `14` | Sets puncture floor `MAX_LEV - 2` and initial depth `min - MAXDEAPTH_LEVEL_DIFF - 2`; requires `max(MAXDEAPTH_LEVEL_DIFF + 2, BSSN_MINDEPTH, 2) <= MAX_LEV <= BSSN_MAXDEPTH`. |
| `BSSN_NYQUIST_M` | `0` | `none` | `7` | Wave-zone Nyquist azimuthal mode (`0` disables). |

Evidence (descriptive behavior): `Dendro/main_cpp.py::output_main_cpp`; `Dendro/solver_context.py::output_solver_context_cpp`; `Dendro/param_toml.py::generate_default_parfile`; `nrpy/examples/q1.par.lowres.toml`.

### Punctures, Initial Data, Outputs, And Apparent-Horizon Finder

See [TwoPunctures inputs](initial-data-and-wave-diagnostics.md#twopunctures-inputs) for which `[BSSN_BH1]`/`[BSSN_BH2]` keys enter the spectral solve versus seed, tracking, and Nyquist estimation.

| Key | Fallback | generated | q1 | Meaning and accepted values |
| --- | --- | --- | --- | --- |
| `BSSN_BH{1,2}.MASS` | `0.5`, `0.5` | `none` | `0.5`, `0.5` | Positive finite puncture masses for `q`, excision, and AMR floors. |
| `BSSN_BH1.{X,Y,Z}`, `BSSN_BH2.{X,Y,Z}` | `(4,0,0)`, `(-4,0,0)` | `none` | `(4,0,0)`, `(-4,0,0)` | Initial excision, seed, and tracker centers. |
| `BSSN_BH{1,2}.V_{X,Y,Z}` | `0.0` (seed/Nyquist); TP `BH1`: `(-0.002284343811437988, 0.11284523509709575)` | `none` | `BH1` TP fallbacks with `V_Z=0`; `BH2 = -BH1` | Initial momenta (`|BH1.V_X|`, `|BH1.V_Y|` in TP without mass division; all 6 in seed/Nyquist; tracker uses evolved `vetU`). |
| `BSSN_BH{1,2}.{SPIN,SPIN_THETA,SPIN_PHI}` | `0.0` | `none` | `0.0` | Seed-only spin inputs (TP solve zeroes spins). |
| `BSSN_BH{1,2}_CONSTRAINT_R` | `1.0` | `none` | `1.0` | Strict (`r < R`) constraint-norm excision radii (`>= 0`, finite). |
| `TPID_PAR_B`, `TPID_NPOINTS_{A,B,PHI}` | `4.0`, `(65, 78, 10)` | `none` | `4.0`, `(65, 78, 10)` | Half-separation and spectral grid counts. |
| `TPID_GIVE_BARE_MASS`, `TPID_TARGET_M_{PLUS,MINUS}` | `1`, `0.48236442246752931` | `none` | `1`, `0.4823644224675293` | Nonzero uses bare masses; `0` iterates ADM masses (`TARGET_M_*` unread). |
| `TPID_NEWTON_TOL`, `TPID_ADM_TOL`, `TPID_FILEPREFIX` | `ID_persist` defaults, `tp` | `none`, `tp` | `1.6e-14`, `3e-16`, `tp_q001` | Solver tolerances and non-empty `<prefix>_nrpy_tpid_sol.bin` stem. |
| `TPID_REPLACE_LAPSE_WITH_SQRT_CHI`, `CHI_FLOOR` | `true`, `0.1` (seed) | `true`, `1e-4` | `true`, `1e-4` | Boolean must be `true`; `CHI_FLOOR` also binds evolved `chi_floor`. |
| `BSSN_PROFILE_FILE_PREFIX`, `BSSN_VTU_FILE_PREFIX`, `BSSN_CHKPT_FILE_PREFIX` | `dat/dgr`, `vtu/<form>_twopunctures`, `cp/<form>_twopunctures` | `dat/dgr`, `none`, `none` | `dat/dgr`, `vtu/bssn_gr`, `cp/bssn_cp` | Output path prefixes. |
| `BSSN_IO_OUTPUT_FREQ`, `BSSN_GW_EXTRACT_FREQ`, `_AFTER_MERGER` | `80`, `80`, `80` | `80`, `80`, `80` | `80`, `80`, `80` | Base cadences for VTU/horizon rows and GW/constraint/ADM/puncture rows. |
| `BSSN_SCALE_VTU_AND_GW_EXTRACTION`, `BSSN_VTU_Z_SLICE_ONLY` | `true`, `true` | `true`, `none` | `true`, `true` | Finest-level cadence scaling ([Constraints And Diagnostic Norms](constraints-and-diagnostic-norms.md#grid-size-and-native-output-cadence)) and `z`-slice VTU toggle. |
| `BSSN_VTU_OUTPUT_{EVOL,CONST}_INDICES`, `BSSN_NUM_{EVOL,CONST}_VARS_VTU_OUTPUT` | `(all, 1)`, `(0..5, 1)` | `none` | `(0..23, 24)`, `(0..5, 6)` | VTU field lists and active prefix counts (`<= len`, valid field indices). |
| `BSSN_CHECKPT_FREQ`, `AEH_SOLVER_FREQ` | `100`, `0` | `none` | `800`, `80` | Checkpoint and apparent-horizon cadences (`0` disables). |
| `BSSN_GW_RADAII`, `BSSN_GW_NUM_RADAII` | `[50.0]`, `len` | `none` | `50..100` step `10`, `6` | Non-empty positive finite extraction radii (`NUM == len`; ADM uses last entry). |
| `BSSN_GW_L_MODES`, `BSSN_GW_NUM_LMODES` | `[2]`, `len` | `none` | `2..8`, `7` | Non-empty `l` list (`NUM == len`, `max(l)` in `[2, 8]`). |
| `[AEH_PARAMS]` `AEH_SAVE_DIR`, `VERBOSITY_LEVEL` | `bah`, `1` | `none` | `bah`, `2` | Horizon output directory and verbosity; 3-element arrays `CFL_FACTOR` (`1.0` / `0.95`), `THETA_L2_M_TOL` (`1e-5`), `THETA_LINF_M_TOL` (`1e-2` / `1e-3`), `MAX_SEARCH_RADIUS` (`1.5` / `0.55,0.55,1.25`), `NR_INTERP_MAX` (`48` / `63`) configure the 3 tracked horizons. |

Evidence (descriptive behavior): `Dendro/main_cpp.py::output_main_cpp`; `nrpy/examples/q1.par.lowres.toml`; `dendro_application_check.py`.

### Registered CodeParameter Coefficients And Startup Rejections

Runtime `CodeParameter` coefficients bind via `Q1_TOML_PARAMETER_NAMES` (writing the unmapped C symbol such as `C_CAHD = 0.15` warns as unread; changed values are rejected on checkpoint restore):

| TOML Key | `CodeParameter` | Default | Scope and notes |
| --- | --- | --- | --- |
| `ETA_CONST` | `eta` | `1.0` | Spatially constant shift-driver damping. |
| `KO_DISS_SIGMA` | `KreissOliger_strength_{gauge,nongauge}` | `0.3` (`0.4` in `pars/<stem>.toml`) | Sets both strengths; exists only when KO is generated. |
| `BSSN_SSL_H`, `BSSN_SSL_SIGMA` | `SSL_h`, `SSL_sigma` | `0.6`, `20.0` | Slow-start lapse amplitude and width. |
| `BSSN_CAHD_C` | `C_CAHD` | `0.06` (BSSN), `0.15` (fCCZ4) | CAHD strength. |
| `CHI_FLOOR` | `chi_floor` | `1.0e-4` | Strictly positive finite floor on `alpha` and `chi` (`sqrt` for `W`); `<= 0` throws `CHI_FLOOR must be finite and positive` when the floor kernel runs. |
| `BSSN_CFL_FACTOR`, `YBS_chi`, `C_YBS_mom` | `CFL_FACTOR`, `YBS_chi`, `C_YBS_mom` | `0.25`, `0.0`, `0.0` | Registered conditionally for `--ybs-momentum` / `--ybs-gamma`. |
| `kappa1`, `kappa2` | `kappa1`, `kappa2` | `0.1`, `0.0` | fCCZ4 Z4 damping coefficients (fCCZ4 only). |

Any key in `pars/q1.par.lowres.toml` omitted from the tables above (including `BSSN_WAVELET_TOL_FUNCTION_R{0,1}` and `BSSN_RK_TYPE`) is unread and triggers a rank-0 warning.

Invalid CLI flags exit with status `2`. Startup aborts all MPI ranks with status `1` before TwoPunctures loading when any range/positivity/finiteness/length constraint in the tables above fails (including mode-6 `r_first <= 8` or `r_last < r_first`, empty `TPID_FILEPREFIX`, or nonfinite `CodeParameter`); in evolution runs (not `--tpid`), `Ctx` construction after mesh build and TwoPunctures loading additionally rejects nonpositive/nonfinite `BSSN_BH{1,2}.MASS` or `BSSN_GW_RADAII` entries, `BSSN_BH{1,2}_MAX_LEV < max(BSSN_MINDEPTH, 2)`, and out-of-range active VTU indices.

Evidence (descriptive behavior): `Dendro/CodeParameters.py::Q1_TOML_PARAMETER_NAMES`; `Dendro/general_relativity/rhs_eval.py`; `Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py`; `dendro_application_check.py::Leg.run_negatives`.

## Sources

- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - `ParameterFile`, host fallbacks, range checks, unread-key warnings
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - `Q1_TOML_PARAMETER_NAMES`, `output_toml_bindings`, `register_CFunctions_parameters`
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - `Ctx` constructor checks
- [param_toml.py](../../../nrpy/infrastructures/Dendro/param_toml.py) - `generate_default_parfile`
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - `CodeParameter` defaults
- [floor_the_lapse_and_conformal_factor.py](../../../nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py) - `chi_floor` positivity check
- [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py) - `kappa1`, `kappa2`
- [dendro_bssn.py](../../../nrpy/examples/dendro_bssn.py) - BSSN `eta` default and CLI
- [dendro_fccz4.py](../../../nrpy/examples/dendro_fccz4.py) - fCCZ4 `eta` default and CLI
- [q1.par.lowres.toml](../../../nrpy/examples/q1.par.lowres.toml) - packaged BBH configuration
- [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py) - `COMMON_OVERRIDES`, `PROFILE_P`, `PROFILE_O`, `Leg.run_negatives`

## See Also

- Depends on: [Project Assembly And Generating Functions](project-assembly-and-emitters.md)
- See also: [BSSN Application Wiring](bssn-application-wiring.md)
- See also: [Initial Data, Wave Extraction, And Output Files](initial-data-and-wave-diagnostics.md)
