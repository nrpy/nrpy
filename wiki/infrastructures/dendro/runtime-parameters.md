# Runtime Parameter Keys

> List the parameter-file keys that the generated Dendro solver reads, their fallback values, and the values it rejects. · Status: confirmed
> Up: [Dendro](index.md)

## Summary

The generated `Dendro_NRPy_BSSN` and `Dendro_NRPy_fCCZ4` executables read one
TOML parameter file, given as `[--tpid] PARAM_FILE`. A key that the file omits
takes the fallback listed below, a key of the wrong TOML type stops the run, and
a key or table member that the solver never reads is reported as having no
effect. Two files are shipped: the generated `pars/<stem>.toml`, which sets a
short list of keys, and the packaged `pars/q1.par.lowres.toml`, which the
printed run commands use and which sets many more. Choose the file knowingly:
the two differ in more than their key lists.

## Detail

### How the keys are resolved

Host keys are read by the generated entry point, and each call site carries its
own fallback. Registered runtime CodeParameters (equation and kernel
coefficients) are bound by name: a CodeParameter uses its own name as the key
unless `Q1_TOML_PARAMETER_NAMES` renames it, as listed under "Registered
coefficients". A host key is read from the root of the file unless the table
below names a table; a line written after a `[table]` header is a member of that
table and is not read as a root key. The columns "generated" and "q1" give the
value the shipped file sets, and "none" means the file omits the key and the
fallback applies. The keys on which the two shipped files give different values
are collected in [BSSN Application
Wiring](bssn-application-wiring.md#q1-versus-generated-defaults).

Claim evidence:
- Claim: The generated solver reads each host key with the call-site fallback listed in this page, accepts an omitted key, stops with a located type error on a wrong TOML type, rejects the values listed under "Values the solver rejects at startup" (all but the solver-context constructor checks before the TwoPunctures data are loaded), and prints `warning: parameter KEY has no effect` on rank 0 for every key or table member that no call site reads. A key under a TOML table header is looked up only through the table-qualified read.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp` (the `ParameterFile` reads in `main`, the startup range checks, and the unread-parameter report); `nrpy/infrastructures/Dendro/solver_context.py`, the `Ctx` constructor within `output_solver_context_cpp` (the context checks).
- Corroboration: `nrpy/examples/tests/dendro_application_check.py`, `Leg.run_negatives` (rejected element order, refinement mode, CFL factor, lapse option, wrong value type, and the unread-key warning) and `Leg.universal_checks` (no unread-parameter warning in the profiles).

### Run control and domain

| Key | Fallback | generated | q1 | Meaning and accepted values |
| --- | --- | --- | --- | --- |
| `BSSN_ID_TYPE` | 0 | 0 | 0 | Must be 0 (TwoPunctures). |
| `BSSN_RESTORE_SOLVER` | 0 | 0 | 0 | A nonzero value requests a checkpoint restore; see [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md). |
| `BSSN_ELE_ORDER` | the `--fd-order` value used at generation (6 by default) | the same value | 6 | Element and finite-difference order; must be 4, 6, or 8. It selects the order-specific kernels; a restore requires it to equal the order stored in the checkpoint and is rejected otherwise; see [Finite-Difference Profiles And Dendro Conformance](finite-difference-profiles-and-dendro-conformance.md). |
| `BSSN_RK_TIME_BEGIN`, `BSSN_RK_TIME_END` | 0.0, 700.0 | 0.0, 1000000.0 | 0.0, 1000000.0 | Evolution runs while time is below the end time; the end time must exceed the begin time. |
| `BSSN_MAX_ITERATIONS` | no limit | none | none | Maximum number of time steps. It is a root key, so a line appended after the last `[table]` header of a file does not set it. |
| `BSSN_CFL_FACTOR` | 0.25 | none | 0.25 | Time step is this factor times the smallest physical axis spacing at the current finest level; positive and finite. This is also the CodeParameter `CFL_FACTOR` when `--ybs-momentum` is generated. |
| `BSSN_TIME_STEP_OUTPUT_FREQ` | 80 | 80 | 80 | Terminal print cadence in steps; the per-step lapse check runs regardless. |
| `BSSN_GRID_MIN_X`, `_Y`, `_Z` | -400.0 | none | -400.0 | Lower domain bounds. |
| `BSSN_GRID_MAX_X`, `_Y`, `_Z` | 400.0 | none | 400.0 | Upper domain bounds; each must exceed the lower bound of its axis, with finite widths. |

### Mesh construction and refinement

| Key | Fallback | generated | q1 | Meaning and accepted values |
| --- | --- | --- | --- | --- |
| `BSSN_MINDEPTH`, `BSSN_MAXDEPTH` | 4, 14 | none | 4, 14 | Octree level range; the minimum must not exceed the maximum, and the maximum must be below 31. |
| `BSSN_REFINEMENT_MODE` | 4 | none | 4 | Only 4 (`BH_WAMR`) is accepted. |
| `BSSN_USE_WAVELET_TOL_FUNCTION` | 0 | none | 6 | Wavelet tolerance mode: 0 is the constant `BSSN_WAVELET_TOL`; 6 is the causal radial profile described in [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md). Other values are rejected. |
| `BSSN_WAVELET_TOL` | 1.0e-5 | none | 1.0e-5 | Wavelet tolerance; positive and finite. |
| `BSSN_WAVELET_TOL_MAX` | `BSSN_WAVELET_TOL` | none | 1.0e-3 | Largest tolerance of the mode-6 profile; positive and finite. |
| `BSSN_GW_REFINE_WTOL` | `BSSN_WAVELET_TOL` | none | 1.0e-4 | Tolerance of the mode-6 profile inside the extraction radii; positive and finite. |
| `BSSN_DENDRO_AMR_FAC` | 0.1 | none | 0.2 | Coarsening factor of the wavelet test; in (0, 1]. |
| `BSSN_DENDRO_AMR_FAC_POST_MERGER` | 0.0 | none | 0.03125 | Coarsening factor used once a post-merger checkpoint has been written, when positive; in [0, 1]. A value of 0 keeps the base factor. |
| `BSSN_REFINE_VARIABLE_INDICES` | every evolved field, in generated order | none | indices 0 to 23 | Evolved fields whose wavelet coefficients enter the refinement test. |
| `BSSN_NUM_REFINE_VARS` | length of the index list | none | 24 | Uses the first entries of the list; must be at least 1 and at most the list length, and every index must be a valid evolved field. |
| `BSSN_REMESH_TEST_FREQ` | 50 | none | 50 | Remesh test cadence in steps. |
| `BSSN_REMESH_TEST_FREQ_AFTER_MERGER` | 10 | none | none | Remesh test cadence once the punctures have merged. |
| `BSSN_INIT_GRID_ITER` | 10 | none | 10 | Initial-grid refinement iterations; at most this many passes (one fewer when the value exceeds 1). |
| `BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE` | true | none | true | Must be true whenever `BSSN_INIT_GRID_ITER` is positive. |
| `BSSN_DENDRO_GRAIN_SZ` | 1000 | none | 50 | Partition grain size passed to the Dendrolib mesh constructor; must be positive. |
| `BSSN_LOAD_IMB_TOL` | 0.1 | none | 0.1 | Load-imbalance tolerance passed to the mesh constructor; must be nonnegative and finite. |
| `BSSN_SPLIT_FIX` | 2 | none | none | Splitter-selection argument passed to the mesh constructor. |
| `BSSN_BH1_AMR_R`, `BSSN_BH2_AMR_R` | 2.0 | none | 1.0 | Radius around each puncture inside which the wavelet tolerance is relaxed and the innermost level floor applies; positive and finite. |
| `BSSN_BH1_MAX_LEV`, `BSSN_BH2_MAX_LEV` | `BSSN_MAXDEPTH` | none | 14 | Level that sets the puncture refinement floors (each puncture's floor starts at this value minus 2); at most `BSSN_MAXDEPTH` and at least `MAXDEAPTH_LEVEL_DIFF + 2`, where `MAXDEAPTH_LEVEL_DIFF` is a Dendrolib constant; the solver context also requires at least `max(BSSN_MINDEPTH, 2)`. The initial octree depth is the smaller of the two values minus `MAXDEAPTH_LEVEL_DIFF` minus 2. |
| `BSSN_AMR_R_RATIO` | 2.0 | none | 1.618033988749895 | Radius ratio between successive level floors around a puncture; must exceed 1. |
| `BSSN_NYQUIST_M` | 0 | none | 7 | Azimuthal mode number of the wave-zone Nyquist refinement; 0 disables it. |

### Punctures and initial data

The tables `[BSSN_BH1]` and `[BSSN_BH2]` hold per-puncture values. Which of them
feed the TwoPunctures solve and which feed only the initial octree seed, the
tracking, and the Nyquist history is stated in [BSSN Application
Wiring](bssn-application-wiring.md#twopunctures-inputs).

| Key | Fallback | generated | q1 | Meaning and accepted values |
| --- | --- | --- | --- | --- |
| `BSSN_BH1.MASS`, `BSSN_BH2.MASS` | 0.5, 0.5 | none | 0.5, 0.5 | Puncture masses used for the mass ratio, the excision and refinement geometry, and the post-merger level floors. |
| `BSSN_BH1.X`, `.Y`, `.Z` | 4.0, 0.0, 0.0 | none | 4.0, 0.0, 0.0 | Initial puncture center for excision, seeding, and tracking. |
| `BSSN_BH2.X`, `.Y`, `.Z` | -4.0, 0.0, 0.0 | none | -4.0, 0.0, 0.0 | Same for the second puncture. |
| `BSSN_BH{1,2}.V_X`, `.V_Y`, `.V_Z` | 0.0 for the initial octree seed and the Nyquist history; for the TwoPunctures momentum `BSSN_BH1.V_X` falls back to -0.002284343811437988 and `BSSN_BH1.V_Y` to 0.11284523509709575 | none | `BSSN_BH1`: -0.002284343811437988, 0.11284523509709575, 0.0; `BSSN_BH2`: the negatives of those values | Initial coordinate velocities. The tracked puncture centers move with the evolved shift, not with these values. |
| `BSSN_BH{1,2}.SPIN`, `.SPIN_THETA`, `.SPIN_PHI` | 0.0 | none | 0.0 | Read, but they feed only the initial octree seed; the TwoPunctures data are nonspinning. |
| `BSSN_BH1_CONSTRAINT_R`, `BSSN_BH2_CONSTRAINT_R` | 1.0 | none | 1.0 | Excision radius of the constraint norms around each tracked puncture; points whose coordinate distance from the center is below the radius are excluded (strict inequality, so 0 excludes nothing). The value must be finite and nonnegative. |
| `TPID_PAR_B` | 4.0 | none | 4.0 | Half of the coordinate separation used by the TwoPunctures solve. |
| `TPID_NPOINTS_A`, `_B`, `_PHI` | 65, 78, 10 | none | 65, 78, 10 | TwoPunctures spectral grid sizes. |
| `TPID_GIVE_BARE_MASS` | 1 | none | 1 | Nonzero: `TPID_TARGET_M_PLUS` and `TPID_TARGET_M_MINUS` are bare masses. Zero: the solve iterates the bare masses to the ADM targets and the two keys are not read. |
| `TPID_TARGET_M_PLUS`, `TPID_TARGET_M_MINUS` | 0.48236442246752931 | none | 0.4823644224675293 | Bare masses when `TPID_GIVE_BARE_MASS` is nonzero. |
| `TPID_NEWTON_TOL`, `TPID_ADM_TOL` | the `ID_persist_struct` defaults | none | 1.6e-14, 3e-16 | Newton and ADM-mass tolerances of the solve. |
| `TPID_FILEPREFIX` | `tp` | `tp` | `tp_q001` | Prefix of the TwoPunctures coefficient file `<prefix>_nrpy_tpid_sol.bin`; must not be empty. |
| `TPID_REPLACE_LAPSE_WITH_SQRT_CHI` | true | true | true | Must be true. |
| `CHI_FLOOR` | 0.1 for the initial-octree seed | 0.0001 | 0.0001 | Also the evolved floor, see "Registered coefficients". |

### Output, checkpoints, and extraction

| Key | Fallback | generated | q1 | Meaning and accepted values |
| --- | --- | --- | --- | --- |
| `BSSN_PROFILE_FILE_PREFIX` | `dat/dgr` | `dat/dgr` | `dat/dgr` | Prefix of the diagnostic files. |
| `BSSN_VTU_FILE_PREFIX` | `vtu/bssn_twopunctures` (BSSN), `vtu/fccz4_twopunctures` (fCCZ4) | none | `vtu/bssn_gr` | Prefix of the VTU files. |
| `BSSN_CHKPT_FILE_PREFIX` | `cp/bssn_twopunctures` (BSSN), `cp/fccz4_twopunctures` (fCCZ4) | none | `cp/bssn_cp` | Prefix of the checkpoint files. |
| `BSSN_IO_OUTPUT_FREQ` | 80 | 80 | 80 | Base cadence of VTU output and of the horizon diagnostics rows. |
| `BSSN_GW_EXTRACT_FREQ`, `BSSN_GW_EXTRACT_FREQ_AFTER_MERGER` | 80, 80 | 80, 80 | 80, 80 | Base cadence of the extraction, constraint, ADM, and puncture-location output before and after the merger. |
| `BSSN_SCALE_VTU_AND_GW_EXTRACTION` | true | true | true | Scales the base cadences with the finest mesh level; see [Constraints And Diagnostic Norms](constraints-and-diagnostic-norms.md#grid-size-and-native-output-cadence). |
| `BSSN_VTU_Z_SLICE_ONLY` | true | none | true | Writes the z-normal slice instead of the full volume. |
| `BSSN_VTU_OUTPUT_EVOL_INDICES`, `BSSN_NUM_EVOL_VARS_VTU_OUTPUT` | every evolved field in generated order; 1 | none | indices 0 to 23; 24 | The first entries of the list are written; the count must not exceed the list length. |
| `BSSN_VTU_OUTPUT_CONST_INDICES`, `BSSN_NUM_CONST_VARS_VTU_OUTPUT` | indices 0 to 5; 1 | none | indices 0 to 5; 6 | Constraint and Psi4 fields in native numbering; the count must not exceed the list length. |
| `BSSN_CHECKPT_FREQ` | 100 | none | 800 | Checkpoint cadence in steps; 0 disables checkpoint writes. |
| `AEH_SOLVER_FREQ` | 0 | none | 80 | Apparent-horizon search cadence in steps; 0 disables the finder. |
| `BSSN_GW_RADAII` | the single radius 50.0 | none | 50.0 to 100.0 in steps of 10.0 | Psi4 extraction radii, nonempty. The ADM integrals use the last listed radius. |
| `BSSN_GW_NUM_RADAII` | length of `BSSN_GW_RADAII` | none | the list length | Must equal the list length. |
| `BSSN_GW_L_MODES` | the single mode 2 | none | 2 to 8 | Psi4 modes l; the largest entry must lie in [2, 8]. |
| `BSSN_GW_NUM_LMODES` | length of `BSSN_GW_L_MODES` | none | the list length | Must equal the list length. |

### Apparent-horizon finder

The table `[AEH_PARAMS]` holds the finder settings. The list-valued keys have
three entries, one for each horizon the finder tracks.

| Key | Fallback | q1 | Meaning |
| --- | --- | --- | --- |
| `AEH_SAVE_DIR` | `bah` | `bah` | Directory of the horizon diagnostics files. |
| `CFL_FACTOR` | 1.0, 1.0, 1.0 | 0.95, 0.95, 0.95 | Hyperbolic-relaxation Courant factor of the finder. |
| `THETA_L2_M_TOL` | 1.0e-5 each | 1.0e-5 each | L2 convergence tolerance of the expansion. |
| `THETA_LINF_M_TOL` | 1.0e-2 each | 1.0e-3 each | Maximum-norm convergence tolerance of the expansion. |
| `MAX_SEARCH_RADIUS` | 1.5 each | 0.55, 0.55, 1.25 | Largest search radius. |
| `NR_INTERP_MAX` | 48 each | 63 each | Largest number of interpolation points. |
| `VERBOSITY_LEVEL` | 1 | 2 | Finder output level. |

Claim evidence:
- Claim: The host keys in the tables above have the stated call-site fallbacks, the stated range checks, and the listed q1 and generated-file values; `BSSN_GW_RADAII`, `BSSN_GW_L_MODES`, and the index lists are read as arrays, and the `[AEH_PARAMS]` members are read only through the table-qualified call.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/main_cpp.py`, `output_main_cpp` (the `parameters.get` calls and range checks in `main`); `nrpy/infrastructures/Dendro/param_toml.py`, `generate_default_parfile` (the generated file's keys); `nrpy/examples/q1.par.lowres.toml` (the q1 values).
- Corroboration: `nrpy/examples/tests/dendro_application_check.py`, `COMMON_OVERRIDES`, `PROFILE_P`, and `PROFILE_O`, which override many of these keys for the short profiles and never set a key the solver does not read.

### Registered coefficients

These CodeParameters keep the default shown when their key is omitted. The
generated file writes the same value, except for `KO_DISS_SIGMA`, which it sets
to 0.4. The solver binds each coefficient to the key shown.

| Key | CodeParameter | Default | Meaning |
| --- | --- | --- | --- |
| `ETA_CONST` | `eta` | 1.0 | Spatially constant shift-driver damping. |
| `KO_DISS_SIGMA` | `KreissOliger_strength_gauge` and `KreissOliger_strength_nongauge` | 0.3 when the key is omitted; the generated file sets 0.4 | One key sets both strengths. The key and both strengths exist only when Kreiss-Oliger dissipation is generated; without it the key is unbound and a line for it is reported as having no effect. |
| `BSSN_SSL_H`, `BSSN_SSL_SIGMA` | `SSL_h`, `SSL_sigma` | 0.6, 20.0 | Amplitude and width of the slow-start lapse relaxation. |
| `BSSN_CAHD_C` | `C_CAHD` | 0.06 (BSSN), 0.15 (fCCZ4) | Coefficient of the constraint-adjusted Hamiltonian damping. |
| `CHI_FLOOR` | `chi_floor` | 1.0e-4 | `alpha` is floored at this value and the evolved conformal factor at the matching chi floor (`sqrt` of it for W). |
| `BSSN_CFL_FACTOR` | `CFL_FACTOR` | 0.25 | Registered only with `--ybs-momentum`; shares the host time-step read. |
| `YBS_chi` | `YBS_chi` | 0.0 | Registered only with `--ybs-gamma`. |
| `C_YBS_mom` | `C_YBS_mom` | 0.0 | Registered only with `--ybs-momentum`. |
| `kappa1`, `kappa2` | `kappa1`, `kappa2` | 0.1, 0.0 | fCCZ4 Z4 damping coefficients; fCCZ4 only. |

A coefficient written under its CodeParameter name when the table lists a
different key is never read: for example a line `C_CAHD = 0.15` is reported as
having no effect, and the key is `BSSN_CAHD_C`. A checkpoint stores these
coefficients and a restore rejects a changed value; see [Octree Grid, AMR, And
Time Stepping](grid-amr-and-time-stepping.md).

Claim evidence:
- Claim: Each registered coefficient is read through the key named in the table (the mapped name when `Q1_TOML_PARAMETER_NAMES` has one, its own name otherwise), keeps the default in the table when its key is omitted (the Kreiss-Oliger strengths 0.3, while the generated file sets `KO_DISS_SIGMA = 0.4`; `eta` 1.0), is bound before the unread-key report, and is rejected at startup with `invalid runtime parameters` when it is nonfinite or fails the generated validation. The Kreiss-Oliger strengths and the key `KO_DISS_SIGMA` exist only when Kreiss-Oliger dissipation is generated.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/Dendro/CodeParameters.py`, `Q1_TOML_PARAMETER_NAMES`, `output_toml_bindings`, and `register_CFunctions_parameters`; `nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py`, `register_CFunction_rhs_eval` (the Kreiss-Oliger call with its 0.3 strengths, and the SSL, CAHD, and YBS registrations and defaults); `nrpy/examples/dendro_bssn.py` and `nrpy/examples/dendro_fccz4.py`, `main` (the `eta` default); `nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py` (`chi_floor`); `nrpy/equations/general_relativity/fCCZ4_RHSs.py` (`kappa1`, `kappa2`); `nrpy/infrastructures/Dendro/param_toml.py`, `generate_default_parfile` (the `ETA_CONST` and `KO_DISS_SIGMA` values written to the generated file).
- Corroboration: `nrpy/examples/tests/dendro_application_check.py`, `Leg.run_negatives`, the nonfinite `BSSN_SSL_SIGMA` and `ETA_CONST` cases.

### Keys of the q1 file that the solver never reads

The packaged q1 file also carries native Dendro-GR keys that the generated
solver has no counterpart for. Each is reported as having no effect: `BSSN_ASYNC_COMM_K`,
`BSSN_DIM`, `BSSN_WAVELET_TOL_FUNCTION_R0`, `BSSN_WAVELET_TOL_FUNCTION_R1`,
`BSSN_RK_TYPE`, `BSSN_RK45_TIME_STEP_SIZE`, `BSSN_RK45_DESIRED_TOL`,
`BSSN_ENABLE_BLOCK_ADAPTIVITY`, `BSSN_BLK_MIN_X`, `BSSN_BLK_MIN_Y`,
`BSSN_BLK_MIN_Z`, `BSSN_BLK_MAX_X`, `BSSN_BLK_MAX_Y`, `BSSN_BLK_MAX_Z`,
`ETA_R0`, `ETA_DAMPING`, `ETA_DAMPING_EXP`, `BSSN_LAMBDA`, `BSSN_LAMBDA_F`,
`BSSN_XI`, `BSSN_TRK0`, `BSSN_ETA_R0`, `BSSN_ETA_POWER`, the table
`TPID_CENTER_OFFSET`, `TPID_VERBOSE`, `TPID_GRID_SETUP_METHOD`, `INITIAL_LAPSE`,
`TPID_INITIAL_LAPSE_PSI_EXPONENT`, `TPID_SOLVE_MOMENTUM_CONSTRAINT`,
`EXTRACTION_VAR_ID`, `EXTRACTION_TOL`, `BSSN_KO_SIGMA_SCALE_BY_CONFORMAL`,
`DENDRO_LOG_FILE`, `DENDRO_LOG_FILE_LEVEL`, `DENDRO_LOG_CONSOLE_LEVEL`, and
`DENDRO_LOG_FORCE_FILE_FLUSH`. In particular the mode-6 wavelet radii of the
native tolerance function are not read: the generated mode-6 profile uses the
fixed radii described in [Octree Grid, AMR, And Time
Stepping](grid-amr-and-time-stepping.md), and `BSSN_RK_TYPE` does not change the
fixed RK4 stepper.

### Values the solver rejects at startup

The solver prints a message and aborts every MPI task with status 1 for each of
these; invalid command-line arguments return status 2.

- `BSSN_ID_TYPE` other than 0, or `TPID_REPLACE_LAPSE_WITH_SQRT_CHI` false.
- `BSSN_ELE_ORDER` other than 4, 6, or 8, and, at a restore, a value that differs from the element order stored in the checkpoint (see [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md)).
- `BSSN_MINDEPTH` above `BSSN_MAXDEPTH`, or `BSSN_MAXDEPTH` of 31 or more.
- `BSSN_NUM_REFINE_VARS` of 0 or above the list length, or a refinement index that is not an evolved field.
- `BSSN_REFINEMENT_MODE` other than 4, or `BSSN_USE_WAVELET_TOL_FUNCTION` other than 0 or 6.
- A nonpositive or nonfinite wavelet tolerance, `BSSN_DENDRO_AMR_FAC` outside (0, 1], `BSSN_DENDRO_AMR_FAC_POST_MERGER` outside [0, 1], a nonpositive or nonfinite `BSSN_CFL_FACTOR`, an end time not above the begin time, or a domain axis with a nonpositive or nonfinite width.
- `BSSN_INIT_GRID_ITER` positive while `BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE` is false.
- A zero `BSSN_DENDRO_GRAIN_SZ`, a negative or nonfinite `BSSN_LOAD_IMB_TOL`, a nonpositive or nonfinite `BSSN_BH{1,2}_AMR_R`, `BSSN_AMR_R_RATIO` not above 1, or `BSSN_BH{1,2}_MAX_LEV` above `BSSN_MAXDEPTH` or below `MAXDEAPTH_LEVEL_DIFF + 2`.
- `BSSN_NUM_EVOL_VARS_VTU_OUTPUT` or `BSSN_NUM_CONST_VARS_VTU_OUTPUT` above the length of its index list.
- A negative or nonfinite `BSSN_BH1_CONSTRAINT_R` or `BSSN_BH2_CONSTRAINT_R`.
- `BSSN_GW_NUM_RADAII` or `BSSN_GW_NUM_LMODES` not equal to its list length, an empty `BSSN_GW_RADAII` or `BSSN_GW_L_MODES`, a largest l outside [2, 8], or, in wavelet mode 6, a first radius not above 8 or a last radius below the first.
- An empty `TPID_FILEPREFIX`, and any registered coefficient that is nonfinite or fails the generated validation.
- The solver-context constructor also rejects a nonpositive or nonfinite `BSSN_BH{1,2}.MASS`, a `BSSN_BH{1,2}_MAX_LEV` below `max(BSSN_MINDEPTH, 2)`, a nonpositive or nonfinite entry of `BSSN_GW_RADAII`, and a VTU field index that is not a valid evolved field or constraint field. These checks run in evolution runs only, after the mesh is built and after the point where a fresh run loads the TwoPunctures file; a `--tpid` run does not reach them.

## Sources

- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - host key reads, fallbacks, range checks, and the unread-key report.
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - key mapping and bindings of registered coefficients.
- [param_toml.py](../../../nrpy/infrastructures/Dendro/param_toml.py) - keys of the generated parameter file.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - coefficient registrations and defaults.
- [floor_the_lapse_and_conformal_factor.py](../../../nrpy/infrastructures/Dendro/general_relativity/floor_the_lapse_and_conformal_factor.py) - `chi_floor`.
- [fCCZ4_RHSs.py](../../../nrpy/equations/general_relativity/fCCZ4_RHSs.py) - `kappa1` and `kappa2`.
- [q1.par.lowres.toml](../../../nrpy/examples/q1.par.lowres.toml) - the packaged parameter file.
- [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py) - profile overrides and negative cases.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Project Assembly And Generating Functions](project-assembly-and-emitters.md)
- See also: [BSSN Application Wiring](bssn-application-wiring.md)
- See also: [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md)
