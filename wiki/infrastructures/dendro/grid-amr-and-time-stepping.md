# Octree Grid, AMR, And Time Stepping

> Explain Dendro's balanced-octree mesh, octant-to-block decomposition, zipped and unzipped data, remeshing, and time-stepper choices. · Status: confirmed
> Up: [Dendro](index.md)

## Summary

Dendro represents the spatial mesh as a unique, sorted, 2:1 balanced octree.
An octant is the unit that receives a refine, coarsen, or no-change decision.
Dendrolib then groups octants into regular `ot::Block` objects so finite-
difference kernels can run on padded Cartesian arrays. A block is therefore a
stencil-execution and storage unit, not the atomic AMR unit.

Use Dendro's own time-stepping terms. Dendrolib calls the uniform schedule UTS
and the spatially adaptive schedule NUTS. The `ts::ETS` class provides the
global explicit evolution used by the generated driver, while
`ts::ExplicitNUTS` implements the nonuniform block schedule. Dendrolib defines
`ot::Block` as a regular block derived from the balanced octree, not as the
octant refinement unit. Use UTS/NUTS terminology rather than describing the
spatial mesh as an "atomic Berger-Oliger grid." The current NRPy-generated
driver instantiates `ts::ETS` with RK4. Its block callbacks match interfaces
needed by local time stepping, but current qualification does not run
`ts::ExplicitNUTS` or remeshing.

## Detail

### Octants And Regular Blocks

`ot::Mesh` accepts an octree that is 2:1 balanced, unique, and sorted. Dendro's
space-filling-curve ordering and partition controls distribute the octants.
Each local octant is an `ot::TreeNode`. Mesh refinement criteria mark it with
`OCT_NO_CHANGE`, `OCT_SPLIT`, or `OCT_COARSE`. Refinement replaces one octant
with its children. Coarsening requires a complete sibling group whose members
agree to coarsen.

Dendrolib decomposes this adaptive octree into a finite sequence of regular
blocks. Each `ot::Block` records its enclosing tree node, regular-grid level,
contiguous local element interval, component offset, allocation dimensions,
element order, and padding width. Generated finite-difference code traverses
these regular blocks. Refinement decisions still belong to octants, not to the
regular block abstraction.

Claim evidence:
- Claim: Dendro's spatial AMR unit is the octant, while `ot::Block` is a regular padded-grid unit formed from the balanced octree for stencil evaluation.
- Role: descriptive behavior
- Deciding authority: [mesh.h](https://github.com/paralab/Dendro-5.01/blob/master/include/mesh.h), `OCT_NO_CHANGE`, `OCT_SPLIT`, `OCT_COARSE`, `ot::Mesh`, and `setMeshRefinementFlags`; [block.h](https://github.com/paralab/Dendro-5.01/blob/master/include/block.h), `ot::Block`
- Corroboration: [mesh.cpp](https://github.com/paralab/Dendro-5.01/blob/master/src/mesh.cpp), `ot::Mesh::octree2BlockDecomposition`

### Zipped And Unzipped Data

Dendrolib exposes two data layouts used by the generated context:

- `OCT_SHARED_NODES` stores the compact octree nodal vector used for evolution
  state and communication.
- `OCT_LOCAL_WITH_PADDING` stores component-major regular-block arrays with
  halo points for stencil evaluation.

`Mesh::unzip` expands octree data into padded block arrays. Its same-level and
coarse/fine paths populate block interiors and halos using the mesh transfer
rules. `Mesh::zip` compresses computed block data back into the octree nodal
representation. Each generated kernel receives only the selected block's
dimensions, spacing, padded physical origin, padding width, and component
offset. It does not inspect the octree or decide coarse/fine transfer.

Claim evidence:
- Claim: Dendrolib owns octree-to-block zip/unzip and coarse/fine data movement; NRPy-generated kernels consume the resulting padded block geometry without owning mesh topology.
- Role: descriptive behavior
- Deciding authority: [mesh.h](https://github.com/paralab/Dendro-5.01/blob/master/include/mesh.h), `ot::Mesh::unzip`, `ot::Mesh::unzip_scatter`, and `ot::Mesh::zip`; [dvec.h](https://github.com/paralab/Dendro-5.01/blob/master/include/dvec.h), `DVEC_TYPE`
- Corroboration: [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), `Ctx::rhs`, `Ctx::rhs_blkwise`, and `block_geometry`

### Remeshing And State Transfer

Refinement criteria mark octants. Dendro-GR provides several application-level
criteria, including wavelet AMR, event-horizon/lapse-threshold refinement,
black-hole location refinement, and combinations of these choices. Dendrolib's
`Mesh::ReMesh` constructs the replacement balanced mesh and repartitions it.
`Mesh::interGridTransfer` moves variables from the old mesh to the new mesh.
Dendro-GR then replaces allocations associated with the old mesh.

`BSSN_ENABLE_BLOCK_ADAPTIVITY` in Dendro-GR is a separate fixed-construction
mode that disables runtime AMR. It must not be confused with octant remeshing
or with Dendrolib's regular block decomposition.

NRPy-generated equation code supplies pointwise and blockwise operations. It
does not select refinement indicators, call `ReMesh`, or own intergrid state
transfer. Those operations remain responsibilities of Dendrolib and the
Dendro-GR application.

Claim evidence:
- Claim: Dendro-GR chooses refinement criteria, Dendrolib constructs the replacement mesh and transfers variables, and current NRPy-generated solver code does not own this remeshing sequence.
- Role: descriptive behavior
- Deciding authority: [mesh.h](https://github.com/paralab/Dendro-5.01/blob/master/include/mesh.h), `ot::Mesh::ReMesh` and `ot::Mesh::interGridTransfer`; [rkBSSN.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/rkBSSN.cpp), the remesh and `intergridTransferVars` sequence
- Corroboration: [grDef.h](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/include/grDef.h), `bssn::RefinementMode`; [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), generated `Ctx` operations

### Uniform And Nonuniform Time Stepping

Dendrolib names four time-stepper modes: `UTS`, `UTS_ADAP`, `NUTS`, and
`NUTS_ADAP`. UTS advances all spatial refinement levels with one uniform step
size. `UTS_ADAP` remains uniform across the mesh while changing that step size
over time. NUTS means spatially adaptive time stepping, and `NUTS_ADAP` permits
the smallest NUTS step to change over time. `ts::ExplicitNUTS` uses block-level
state, synchronization, correction, block unzip/zip, and block RHS callbacks.

This terminology separates two independent choices:

- Spatial AMR decides which octants exist and maintains 2:1 balance.
- The time stepper decides whether all refinement levels advance uniformly or
  through Dendro's NUTS schedule.

The generated NRPy real-host driver constructs `ts::ETS`, selects RK4, and
calls `evolve()` once per requested step. Generated `rhs_blk` and related block
hooks make the context structurally compatible with Dendrolib blockwise calls,
but that interface alone does not prove `ExplicitNUTS` behavior. The real-host
tests call block callbacks directly and exercise a fixed mesh; they do not run
the NUTS scheduler, remesh, or intergrid transfer.

Claim evidence:
- Claim: Dendrolib distinguishes UTS from NUTS; the current NRPy-generated driver instantiates `ts::ETS`, while Dendrolib provides `ts::ExplicitNUTS` for the NUTS schedule, which current qualification does not exercise.
- Role: descriptive behavior
- Deciding authority: [ets.h](https://github.com/paralab/Dendro-5.01/blob/master/ODE/include/ets.h), `ts::TimeStepperType` and `ts::ETS`; [enuts.h](https://github.com/paralab/Dendro-5.01/blob/master/ODE/include/enuts.h), `ts::ExplicitNUTS`; [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py), `_REAL_MAIN`
- Corroboration: [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp), direct whole-vector and block-callback checks

### Responsibility Boundary

| Operation | Owner |
| --- | --- |
| Octree construction, balance, partition, and regular-block decomposition | Dendrolib |
| Refinement criterion and remesh policy | Dendro-GR application |
| Zip, unzip, halo population, and intergrid transfer | Dendrolib |
| Symbolic equations and generated point/block kernels | NRPy Dendro infrastructure |
| Selected runtime time stepper | Generated or surrounding Dendro-GR driver |

## Sources

- [block.h](https://github.com/paralab/Dendro-5.01/blob/master/include/block.h) - `ot::Block` and octree-to-regular-block description
- [mesh.h](https://github.com/paralab/Dendro-5.01/blob/master/include/mesh.h) - balanced sorted octree input, mesh refinement flags, `ot::Mesh::unzip`, `zip`, `ReMesh`, and `interGridTransfer`
- [mesh.cpp](https://github.com/paralab/Dendro-5.01/blob/master/src/mesh.cpp) - `ot::Mesh::octree2BlockDecomposition`, `setMeshRefinementFlags`, and `ReMesh`
- [dvec.h](https://github.com/paralab/Dendro-5.01/blob/master/include/dvec.h) - `DVEC_TYPE` and `ot::DVector`
- [ets.h](https://github.com/paralab/Dendro-5.01/blob/master/ODE/include/ets.h) - `ts::TimeStepperType` and `ts::ETS`
- [enuts.h](https://github.com/paralab/Dendro-5.01/blob/master/ODE/include/enuts.h) - `ts::ExplicitNUTS`
- [rkBSSN.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/rkBSSN.cpp) - Dendro-GR refinement, remesh, and intergrid-transfer sequence
- [grDef.h](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/include/grDef.h) - `bssn::RefinementMode`
- [main_cpp.py](../../../nrpy/infrastructures/Dendro/main_cpp.py) - generated `ts::ETS` RK4 driver
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - DVector layouts, zip/unzip calls, block callbacks, and geometry extraction
- [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp) - fixed-mesh whole-vector and block-callback qualification

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Project Assembly And Generating Functions](project-assembly-and-emitters.md)
- See also: [Finite-Difference Profiles And Dendro Conformance](finite-difference-profiles-and-dendro-conformance.md)
- See also: [Gridfunctions, Naming, And Loops](gridfunctions-naming-and-loops.md)
- Validated by: [Validation, Standalone Host, And Deferred Tests](validation-standalone-host-and-deferral-gates.md)
