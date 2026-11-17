# Initial Data, Wave Extraction, And Output Files

> Shared TwoPunctures parameter consumption, finite-radius ADM sphere integrals, Lebedev Psi4 mode and L2 files, puncture tracking, TOML parameter dumps, and VTU output in `Dendro_NRPy_BSSN` and `Dendro_NRPy_fCCZ4`. · Status: provisional
> Up: [Dendro](index.md)

## Summary

`Dendro_NRPy_BSSN` and `Dendro_NRPy_fCCZ4` share the same TwoPunctures input-key partition, BHaHAHA-styled `# column N = ...` ASCII headers, finite-radius `*_ADM.dat` sphere integrals, `LEBEDEV_025` `Psi4` mode (`*_GW_l<l>_m<m>.dat`) and unweighted sample-norm (`*_GW_L2.dat`) files, `*_BHLocations.dat` puncture coordinates, startup `__PARAM_DUMP__*.toml` files, and `VTU` slice/volume output.

## Detail

### TwoPunctures inputs

TwoPunctures consumes only a subset of black-hole TOML keys; the remaining keys are still read by the initial octree seed, puncture tracker, or Nyquist estimator so no unread-key warning is emitted:
- **Positions**: the spectral solve uses separation `2 * TPID_PAR_B` along `x` and mass ratio with an internal center-of-mass shift, ignoring `BSSN_BH{1,2}.{X,Y,Z}` (which set excision centers, octree seed, and initial tracked positions—keep them consistent with `TPID_PAR_B`) and never reading `TPID_CENTER_OFFSET`.
- **Momenta**: with `p_t = |BSSN_BH1.V_Y|` and `p_r = |BSSN_BH1.V_X|` (omitted fallbacks `0.11284523509709575` and `-0.002284343811437988`), the solve sets `P_plus = (-p_r, +p_t, 0)` and `P_minus = -P_plus`; signs, `V_Z`, and all `BSSN_BH2.V_*` keys feed only the octree seed and initial Nyquist wavelength.
- **Spins**: the solve zeroes both spins (nonspinning data); `BSSN_BH{1,2}.{SPIN,SPIN_THETA,SPIN_PHI}` feed only the octree seed.
- **Masses**: `BSSN_BH{1,2}.MASS` define `q = max(M1, M2) / min(M1, M2)` plus excision/refinement geometry. When `TPID_GIVE_BARE_MASS != 0` (default), `TPID_TARGET_M_{PLUS,MINUS}` are bare masses; when `0`, bare masses are iterated to ADM targets `q / (1 + q)` and `1 / (1 + q)` without reading `TPID_TARGET_M_*`.

Evidence (descriptive behavior): `Dendro/main_cpp.py::output_main_cpp`; `TwoPunctures/ID_persist_struct.py::register_CFunction_initialize_ID_persist_struct`; `dendro_application_check.py::Leg.check_run_a`.

### Finite-Radius ADM Sphere Integrals And Labeled ASCII Headers

Apparent-horizon searches interpolate evolved fields and convert search points to ADM variables, while `BSSN_to_ADM` reconstructs grid fields for `*_ADM.dat`. At nonzero effective gravitational-wave cadence (`BSSN_GW_EXTRACT_FREQ`), rank `0` writes step, time, radius `r` (the last entry of non-empty `BSSN_GW_RADAII` about `(0, 0, 0)`; wavelet mode `6` requires `r_first > 8` and `r_last >= r_first` without enforcing monotonicity), and finite-radius Lebedev quadrature (`4 pi r^2 w`, coordinate normal `n^i`, up to 10 significant digits):

```text
E   = (1 / 16 pi) oint (d_j gamma_ij - d_i gamma_jj) n^i dS,
P_i = (1 / 8 pi)  oint (K_ij - gamma_ij K) n^j dS,
J_i = (1 / 8 pi)  oint eps_ijk x^j (K_kl - gamma_kl K) n^l dS.
```

For `gamma_ij = psi^4 delta_ij` with `psi = 1 + E / (2r)`, the finite-radius integral gives `E (1 + E / (2r))^3` (first-order fractional excess `3E / (2r)`), and `dendro_application_check.py` scales solved TwoPunctures `E` and `J_z` by `(1 + E / (2r))^3` and `1 / (1 + E / (2r))^2`.

All ASCII diagnostic tables (`*_ADM.dat`, constraint norms, `*_GW_l<l>_m<m>.dat`, `*_GW_L2.dat`, `*_BHLocations.dat`) use BHaHAHA-style comment headers (`open_labeled_output`): when creating an empty file, rank `0` writes `# <title>`, `#`, and `# column N = <name>: <meaning>` (`TimeStep` and `time`, or `t` in GW files), whereas native `BSSN_GR` writes uncommented first lines (`pandas.read_csv` callers such as `getstrain.py` must skip `#` lines). Restarts append below existing headers.

Evidence (descriptive behavior): `Dendro/solver_context.py::output_solver_context_cpp`; `Dendro/general_relativity/adm_quantities.py`; `Dendro/general_relativity/adm_quantities_surface_data.py`; `BHaHAHA/diagnostics_file_output.py`.

### Lebedev Psi4 Modes, GW L2 Norms, And Puncture Locations

Waveform extraction evaluates `Psi4` directly from evolved BSSN/fCCZ4 fields using the Baker-Campanelli-Lousto tetrad ([Psi4 And Tetrads](../../equations/general-relativity/psi4-and-tetrads.md)) seeded with `v1 = (-y, x, 0)` (`(1, 0, 0)` on the `z` axis) and `v2 = (x, y, z)`. On each origin-centered radius in `BSSN_GW_RADAII`, 25th-order Lebedev quadrature (`LEBEDEV_025`) projects `4 pi sum_i w_i Psi4(theta_i, phi_i) conj(sY_lm(theta_i, phi_i))` with **no** factor of `r` (matching native `Dendro-GR`), writing scientific-notation `(Re, Im)` pairs (`10` decimal places) to `*_GW_l<l>_m<m>.dat` for every `l` in `2 .. max(BSSN_GW_L_MODES)`.

At the same cadence, rank `0` appends `<BSSN_PROFILE_FILE_PREFIX>_GW_L2.dat` containing MPI-reduced unweighted sample norms `sqrt(sum(Re(Psi4)^2))` and `sqrt(sum(Im(Psi4)^2))` per radius (skipping the already-written checkpoint step on restart) and `<BSSN_PROFILE_FILE_PREFIX>_BHLocations.dat` containing tracked puncture/excision `(x, y, z)` centers (not apparent-horizon centers).

Evidence (descriptive behavior): `Dendro/general_relativity/gravitational_waves.py`; `Dendro/general_relativity/psi4_eval.py`; `Dendro/solver_context.py`; [Dendro-GR gwExtract.h](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/include/gwExtract.h).

### Startup Parameter Dump And VTU Output

On normal evolution startup (not `--tpid`), rank `0` writes `<BSSN_PROFILE_FILE_PREFIX>__PARAM_DUMP__YYYY-MM-DD-HH-MM-SS.toml` in full round-trip `double` precision. Missing keys receive defaults only when every host fallback and mapped `CodeParameter` default agree; keys with distinct call-site fallbacks—such as `BSSN_BH1.V_{X,Y}` (octree seed/Nyquist vs TwoPunctures momentum) and `CHI_FLOOR` (`0.1` puncture-seed floor vs `1e-4` evolved-field floor)—remain omitted so rerunning the dump preserves identical behavior and TwoPunctures cache reuse.

`Ctx::write_vtu` writes `z`-plane elements when `BSSN_VTU_Z_SLICE_ONLY=true` (default) or the full octree otherwise. `BSSN_VTU_OUTPUT_EVOL_INDICES` (first `BSSN_NUM_EVOL_VARS_VTU_OUTPUT` entries, default `1` -> `alpha`) indexes NRPy's evolved fields in `BSSN_GR` slot order (`cf`, `lambdaU`, `hDD`, `aDD`, rescaled `betU`, plus `Theta_fCCZ4` at index `24` in fCCZ4); `BSSN_VTU_OUTPUT_CONST_INDICES` (first `BSSN_NUM_CONST_VARS_VTU_OUTPUT` entries, default `1` -> `H`) selects `H`, `MD0..MD2`, `psi4_real`, `psi4_imag` (`H_Z4` and `Z4constraintU` are not VTU-selectable), evaluating constraints or `Psi4` on demand.

Evidence (descriptive behavior): `Dendro/main_cpp.py::ParameterFile`; `Dendro/CodeParameters.py::output_toml_default_assignments`; `Dendro/solver_context.py::output_solver_context_cpp`; [Dendro-GR bssnCtx.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/bssnCtx.cpp).

## Sources

- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - `Ctx::adm_output`, `Ctx::gravitational_wave_output`, `Ctx::black_hole_locations_output`, `Ctx::write_vtu`, `open_labeled_output`
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - `ParameterFile`, `output_main_cpp`
- [CodeParameters.py](../../../nrpy/infrastructures/Dendro/CodeParameters.py) - `output_toml_default_assignments`
- [adm_quantities.py](../../../nrpy/infrastructures/Dendro/general_relativity/adm_quantities.py) - `register_CFunction_adm_quantities`
- [adm_quantities_surface_data.py](../../../nrpy/infrastructures/Dendro/general_relativity/adm_quantities_surface_data.py) - `register_CFunction_adm_quantities_surface_data`
- [gravitational_waves.py](../../../nrpy/infrastructures/Dendro/general_relativity/gravitational_waves.py) - `register_CFunction_gravitational_waves`
- [psi4_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/psi4_eval.py) - `register_CFunction_psi4_eval`
- [ID_persist_struct.py](../../../nrpy/infrastructures/BHaH/general_relativity/TwoPunctures/ID_persist_struct.py) - `register_CFunction_initialize_ID_persist_struct`
- [diagnostics_file_output.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_file_output.py) - column-label header convention
- [dendro_application_check.py](../../../nrpy/examples/tests/dendro_application_check.py) - `Leg.check_run_a`, `parse_table`, `parse_modes`
- [Dendro-GR gwExtract.h](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/include/gwExtract.h) - `GW::extractFarFieldPsi4`
- [Dendro-GR dataUtils.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/dataUtils.cpp) - `writeBHCoordinates`
- [Dendro-GR bssnCtx.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/bssnCtx.cpp) - `BSSNCtx::write_vtu`
- [Dendro-GR bssngr_main.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/bssngr_main.cpp) - `writeParamTOMLFile`
- [Dendro-GR parameters.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/parameters.cpp) - VTU and parameter writer defaults

## See Also

- Depends on: [Psi4 And Tetrads](../../equations/general-relativity/psi4-and-tetrads.md)
- See also: [BSSN Application Wiring](bssn-application-wiring.md)
- Validated by: [Production Validation And Deferred Checks](production-validation-and-deferred-checks.md)
