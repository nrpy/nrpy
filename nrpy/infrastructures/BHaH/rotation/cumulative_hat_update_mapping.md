# Cumulative Hat Update Mapping (NRPy SO(3) Rotation)

Date: 2026-03-03  
Scope: active `nrpy/` tree only

## Owner Struct And Fields
- Struct: `commondata_struct`
- Fields:
  - `cumulative_regrid_xhatU[3]`
  - `cumulative_regrid_yhatU[3]`
  - `cumulative_regrid_zhatU[3]`

## Authoritative Write-Path Audit Result
- Active authoritative writes in `nrpy/`: **none found**.
- Active read-side consumers in `nrpy/infrastructures/BHaH/rotation/`:
  - `unrotate_xCart_to_fixed_frame`
  - `unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame`

## Reproducible Audit Commands And Outputs
1. Command:
   - `rg -n "cumulative_regrid_[xyz]hatU" nrpy/`
   Output (abridged to active files):
   - `nrpy/infrastructures/BHaH/rotation/unrotate_xCart_to_fixed_frame.py:69-71`
   - `nrpy/infrastructures/BHaH/rotation/unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame.py:78-80`
   - `nrpy/infrastructures/BHaH/rotation/so3_matrix_ops.py:26,35,44`
2. Command:
   - `rg -n "commondata\.cumulative_regrid_[xyz]hatU\[[0-2]\]\s*=\s*" nrpy/`
   Output:
   - no matches (exit code `1`)
3. Command:
   - `rg -n "cumulative_regrid_[xyz]hatU\[[0-2]\]\s*=\s*" nrpy/`
   Output:
   - no matches (exit code `1`)
4. Command:
   - `rg -n "commondata\.cumulative_regrid_[xyz]hatU" nrpy/infrastructures/BHaH/rotation/`
   Output:
   - read-only call sites in the two consumer functions above

## Active Read-Path Mapping (File / Function / Call Order / Owner)
| File | Function | Call Order | Owner |
|---|---|---|---|
| `nrpy/infrastructures/BHaH/rotation/unrotate_xCart_to_fixed_frame.py` | `unrotate_xCart_to_fixed_frame` | read hats -> `so3_validate_and_optionally_fix_hats(..., do_fix=0)` -> `build_R_from_cumulative_hats` -> `so3_apply_R_to_vector` | `commondata_struct` cumulative hats |
| `nrpy/infrastructures/BHaH/rotation/unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame.py` | `unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame` | read hats -> `so3_validate_and_optionally_fix_hats(..., do_fix=0)` -> `build_R_from_cumulative_hats` -> `so3_matrix_to_axis_angle` -> set part2 identity | `commondata_struct` cumulative hats |

## Archive Note (Out Of Active Path)
- Rotation-related material exists under `ancient/` in the repository root.
- That tree is not part of active `nrpy/` execution or code generation paths used by this refactor.

## Integration Hook For Future Writer Sites
- If an authoritative cumulative-hat update site is introduced in active `nrpy/` codegen, invoke:
  - `so3_validate_and_optionally_fix_hats(xhatU, yhatU, zhatU, 1)`
  immediately after each update.
- Active consumer policy remains validate-only (`do_fix=0`) with hard-fail on invalid hats.
