# Finite-Difference Profiles And Dendro Conformance

> Explain why NRPy's Dendro profiles use centered advection and explicit Kreiss-Oliger base orders, how their stencil reach fits Dendro blocks, and where they differ from Dendro-GR. · Status: confirmed
> Up: [Dendro](index.md)

## Summary

NRPy's Dendro profiles use centered derivatives for advection terms and pair
regular finite-difference orders 4, 6, and 8 with Kreiss-Oliger (KO) base
orders 2, 4, and 6. The KO base order is not the order of the resulting
difference: NRPy's `dKOD` construction adds two, so the effective KO
differences have orders 4, 6, and 8. Centered regular derivatives and KO terms
therefore fit the same 2-, 3-, or 4-point block padding.

These are Dendro geometry-compatible profile choices, not demonstrated performance
optimizations. They keep each generated kernel within the padded regular-block
geometry supplied by Dendrolib and agree with important public Dendro-GR
interior-stencil conventions. They do not reproduce every Dendro-GR stencil. In
particular, NRPy's order-8 profile uses an eighth-difference KO term, while
Dendro-GR's order-8, padding-4 selection uses a sixth-difference interior KO
stencil with special physical-boundary closures.

## Detail

### NRPy profile meaning

| Regular order | KO base order | Effective KO difference | Required padding per side |
| --- | --- | --- | --- |
| 4 | 2 | 4 | 2 |
| 6 | 4 | 6 | 3 |
| 8 | 6 | 8 | 4 |

`DENDRO_FD_PROFILES` records the KO base order and required padding for each
regular order. `compute_fdcoeffs_fdstencl` adds two to the explicit KO base
order before it constructs the stencil. Thus, “lower-order KO” describes the
input parameter only. It must not be read as a KO difference two orders below
the regular derivative.

Without the explicit lower base order, NRPy's default `dKOD` construction
would add two to the regular finite-difference order. Its symmetric stencil
would then reach one point farther than a centered derivative of that regular
order. The selected base orders instead give both operator families the same
maximum reach. This keeps the order-4, order-6, and order-8 generated kernels
within Dendro padding widths 2, 3, and 4.

Claim evidence:
- Claim: NRPy's Dendro profiles use regular/KO-base pairs 4/2, 6/4, and 8/6; `dKOD` adds two to the KO base order, so the effective KO differences have orders 4, 6, and 8 and require the same 2-, 3-, and 4-point reach as the centered regular derivatives.
- Role: descriptive behavior
- Deciding authority: [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), `DENDRO_FD_PROFILES` and `build_rhs_eval`; [finite_difference.py](../../../nrpy/finite_difference.py), `compute_fdcoeffs_fdstencl`
- Corroboration: [block_kernel_helpers.py](../../../nrpy/infrastructures/Dendro/block_kernel_helpers.py), `padding_from_derivative_operators`; [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp), padding checks for element orders 4, 6, and 8

### Centered advection

The BSSN and fCCZ4 equation factories express shift advection with directional
derivative symbols. Before C code generation, the Dendro RHS builder replaces
the supported one-point upwind and downwind symbols with centered `dD`
symbols. It then uses the canonical derivative-operator extraction and rejects
any directional operator that remains. Generated Dendro RHS kernels therefore
contain centered advection derivatives for every supported profile.

For interior points, this choice matches formulas reached by the public
Dendro-GR order-6 and order-8 advection functions: `deriv644adv_*` calls
`deriv644_*`, and `deriv8642adv_*` calls `deriv8642_*`. NRPy-generated kernels
do not implement Dendro-GR's `bflag`-dependent physical-boundary derivative
closures.

Dendro-GR retains directional order-4 `deriv42adv_*` implementations and
conditionally generated advection call sites. The reviewed public build does
not wire `BSSN_USE_ADVECTIVE_DERIVS` to an `adv_deriv_*` selection, so these
functions establish a source-level alternative, not a qualified selectable
path. The supported comparison is narrow: NRPy's centered interior treatment
matches Dendro-GR's public high-order advection functions, not that Dendro-GR
universally forbids upwinding.

Claim evidence:
- Claim: supported NRPy Dendro RHS profiles emit centered advection derivatives; on interior points this matches formulas reached by Dendro-GR's order-6 and order-8 advection functions, while NRPy does not implement Dendro-GR's `bflag`-dependent physical-boundary closures and the reviewed public build does not wire its retained directional order-4 functions into a selectable path.
- Role: descriptive behavior
- Deciding authority: [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), centered directional-symbol replacement and remaining-operator rejection in `build_rhs_eval`
- Corroboration: [Dendro-GR derivs.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/derivs.cpp), `deriv42adv_*`, `deriv644adv_*`, and `deriv8642adv_*`; [Dendro-GR derivs.h](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/include/derivs.h), commented `adv_deriv_*` mappings; [Dendro-GR bssnrhs_derivs_adv.h](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/scripts/bssnrhs_derivs_adv.h), `BSSN_USE_ADVECTIVE_DERIVS`; [Dendro-GR CMakeLists.txt](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/CMakeLists.txt), derivative-order options

### KO comparison with Dendro-GR

Dendro-GR's `set_appropriate_derivs` selects `ko_deriv21_*` for order 4 with
padding 2 and `ko_deriv42_*` for order 6 with padding 3. These selections have
the same effective difference order and reach as NRPy's 4/2 and 6/4 profiles.
For order 8 with padding 4, Dendro-GR instead selects
`ko_pw4_deriv42_*`. Its interior expression is a symmetric sixth difference
over offsets -3 through 3, while the routine reserves padding 4 and supplies
separate physical-boundary formulas. NRPy's 8/6 profile uses a symmetric
eighth difference over offsets -4 through 4.

Therefore all three NRPy profiles fit the Dendro block padding selected for the
same regular derivative order. Only the order-4 and order-6 comparisons also
match Dendro-GR's selected KO difference order and interior reach. Order 8
matches host padding, not Dendro-GR's exact KO interior stencil or boundary
closures.

Claim evidence:
- Claim: NRPy's order-4 and order-6 KO profiles match Dendro-GR's selected effective KO difference order and reach; NRPy's order-8 profile matches Dendro padding 4 but uses an eighth-difference KO stencil instead of Dendro-GR's sixth-difference interior `ko_pw4_deriv42_*` stencil and its boundary closures.
- Role: descriptive behavior
- Deciding authority: [Dendro-GR derivs.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/derivs.cpp), `set_appropriate_derivs`, `ko_deriv21_*`, `ko_deriv42_*`, and `ko_pw4_deriv42_*`
- Corroboration: [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), `DENDRO_FD_PROFILES`; [finite_difference.py](../../../nrpy/finite_difference.py), `compute_fdcoeffs_fdstencl`

### Relation to Dendro block geometry

Dendrolib unzips octree data into padded, component-major regular-block
arrays. NRPy-generated kernels receive one block's dimensions, spacing,
padded origin, padding width, and component offset. They do not request another
halo exchange or own coarse/fine interpolation. A generated stencil must
therefore fit the padding already attached to that block.

Centered advection and the explicit KO base orders make this stencil reach
equal to half the regular finite-difference order. Real and standalone host
checks reject any block whose padding does not exactly equal the generated
profile's required padding. This is the documented numerical and storage
reason for these choices. No cited source establishes that they are faster,
more stable, or more accurate than directional advection or another KO profile
in a production evolution.

Claim evidence:
- Claim: current NRPy Dendro stencil profiles require exact agreement with the padding of Dendrolib regular blocks without requiring generated kernels to own halo exchange or coarse/fine transfer; this is a geometry-compatibility statement, not a performance or stability result.
- Role: descriptive behavior
- Deciding authority: [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py), `Ctx::rhs`, `Ctx::rhs_blkwise`, and `block_geometry`; [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py), required-padding check
- Corroboration: [Dendrolib block.h](https://github.com/paralab/Dendro-5.01/blob/master/include/block.h), `ot::Block`; [Dendrolib mesh.h](https://github.com/paralab/Dendro-5.01/blob/master/include/mesh.h), `ot::Mesh::unzip` and `ot::Mesh::zip`

### Qualification boundary

Generated self-tests are parameterized by formulation, regular order, and KO
mode. The `dendro-validation` CI job generates and runs all twelve BSSN/fCCZ4,
order-4/6/8, and KO-off/on combinations. Their independent multiprecision
stencil evaluation checks the NRPy-generated operators. Real-host capability
tests check padding 2, 3, and 4, padded geometry, component layout, unzip/zip
data movement, and the FD4/6/8 KO-on plus FD6 KO-off runtime paths for both
formulations.

The independent KO oracle compares exact coefficients and signs against
Dendro-GR's `ko_deriv21`, `ko_deriv42`, and radius-four `ko_deriv64` interior
functions. The configured Dendro-GR order-8 selector still uses another KO
family; this comparison establishes the required radius-four operator, not
selector equivalence or physical-boundary closure equivalence.

These checks do not compare centered and directional advection, benchmark the
profiles, or reproduce Dendro-GR's `bflag`-dependent physical-boundary closures
for regular, advection, or KO derivatives. The checks do not establish long-time stability and
convergence on a remeshing AMR evolution or exercise Dendrolib's NUTS schedule.
Such claims require separate numerical experiments against the intended
Dendro application.

Claim evidence:
- Claim: current qualification checks the emitted NRPy operators across all twelve formulation/order/KO combinations, compares KO interior coefficients with the corresponding Dendro-GR radius-two, radius-three, and radius-four functions, and checks the selected real-host profiles; it does not qualify comparative performance, directional advection, Dendro-GR's `bflag`-dependent derivative closures or order-8 selector, remeshing, NUTS, or long-time evolution.
- Role: descriptive behavior
- Deciding authority: [main.yml](../../../.github/workflows/main.yml), `dendro-validation`; [general_relativity/self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py), generated numerical checks; [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp), host geometry checks; [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp), direct runtime checks
- Corroboration: [main_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/main_cpp.py), CTest profile registration

## Sources

- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - `DENDRO_FD_PROFILES`, centered advection normalization, operator extraction, and padding checks
- [block_kernel_helpers.py](../../../nrpy/infrastructures/Dendro/block_kernel_helpers.py) - canonical operator extraction and stencil-reach calculation
- [finite_difference.py](../../../nrpy/finite_difference.py) - `compute_fdcoeffs_fdstencl` and explicit `ko_fd_order` semantics
- [general_relativity/self_tests_cpp.py](../../../nrpy/infrastructures/Dendro/general_relativity/self_tests_cpp.py) - independent generated numerical checks
- [dendrolib_capability_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/dendrolib_capability_test.cpp) - real-host padding and layout checks
- [runtime_integration_test.cpp](../../../nrpy/infrastructures/Dendro/tests_infra/runtime_integration_test.cpp) - real-host runtime checks
- [main.yml](../../../.github/workflows/main.yml) - complete standalone formulation/order/KO matrix
- [Dendro-GR derivs.cpp](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/src/derivs.cpp) - regular, advection, and KO derivative implementations and selection
- [Dendro-GR derivs.h](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/include/derivs.h) - derivative declarations and retained commented advection mappings
- [Dendro-GR CMakeLists.txt](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/CMakeLists.txt) - public derivative-order build selections
- [Dendro-GR bssnrhs_derivs_adv.h](https://github.com/paralab/Dendro-GR/blob/master/BSSN_GR/scripts/bssnrhs_derivs_adv.h) - optional advection-derivative work
- [Dendrolib block.h](https://github.com/paralab/Dendro-5.01/blob/master/include/block.h) - regular-block geometry
- [Dendrolib mesh.h](https://github.com/paralab/Dendro-5.01/blob/master/include/mesh.h) - zip/unzip ownership

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Octree Grid, AMR, And Time Stepping](grid-amr-and-time-stepping.md)
- Depends on: [Finite Differences](../../core/finite-difference.md)
- See also: [BSSN Application Wiring](bssn-application-wiring.md)
- See also: [fCCZ4 Application Wiring](fccz4-application-wiring.md)
- Validated by: [Validation, Standalone Host, And Deferred Tests](validation-standalone-host-and-deferral-gates.md)
