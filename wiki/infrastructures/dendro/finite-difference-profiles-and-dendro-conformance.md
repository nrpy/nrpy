# Finite-Difference Profiles And Dendro Conformance

> Define supported centered finite-difference and Kreiss-Oliger profiles. · Status: confirmed
> Up: [Dendro](index.md)

## Summary

Dendro generation emits centered regular orders 4, 6, and 8 with KO orders 2,
4, and 6. Their stencil radii are 2, 3, and 4 points. No upwind derivative is
generated.

## Detail

| Regular order | KO finite-difference order | KO derivative order | Required padding |
| --- | --- | --- | --- |
| 4 | 2 | 4 | 2 |
| 6 | 4 | 6 | 3 |
| 8 | 6 | 8 | 4 |

One generated application contains every profile. Runtime dispatch accepts
Dendro element order 4, 6, or 8, selects the matching Ricci, RHS, constraint,
and derivative-dependent initial-data kernels, and requires padding at least
`order/2`. A mismatched KO order or unsupported element order is a configuration
error.

Kreiss-Oliger dissipation is enabled in the production profile. The generated
default uses order 6 with KO order 4. Validation must exercise each profile in
the complete application; source presence alone does not establish correct
dispatch or block indexing.

## Sources

- [finite_difference.py](../../../nrpy/finite_difference.py) - finite-difference and KO stencil construction.
- [rhs_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/rhs_eval.py) - order-specific RHS registration.
- [Ricci_eval.py](../../../nrpy/infrastructures/Dendro/general_relativity/Ricci_eval.py) - order-specific Ricci registration.
- [solver_context.py](../../../nrpy/infrastructures/Dendro/solver_context.py) - runtime order dispatch and padding checks.
- [Dendrolib block.h](https://github.com/paralab/Dendro-5.01/blob/master/include/block.h) - regular block geometry.

## See Also

- Parent: [Dendro](index.md)
- Depends on: [Finite Difference](../../core/finite-difference.md)
- Validated by: [Production Validation And Deferred Checks](validation-standalone-host-and-deferral-gates.md)
