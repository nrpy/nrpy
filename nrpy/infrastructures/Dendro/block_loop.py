# nrpy/infrastructures/Dendro/block_loop.py
"""
Dendro numerical block-loop emission built on the generic NRPy loop helper.

The block loop iterates the local block list supplied by Dendro and invokes
the registered per-block CFunction for each block.  The loop is emitted by
NRPy; the Dendro runtime only provides the block list and count.
"""

import nrpy.helpers.loop as lp


def block_loop(loop_body: str, count: str = "num_blocks") -> str:
    """
    Emit a Dendro numerical block loop around a loop body.

    :param loop_body: C code executed once per block; has access to ``blk_id``.
    :param count: C expression for the number of local blocks.
    :return: The generated loop C code string.

    Doctests:
    >>> print(block_loop("fccz4_rhs_block(blk_id);", count="NBLK"))
    for (int blk_id = 0; blk_id < static_cast<std::ptrdiff_t>(NBLK); blk_id++) {
    fccz4_rhs_block(blk_id);
    } // END LOOP: for blk_id over [0, static_cast<std::ptrdiff_t>(NBLK))
    <BLANKLINE>
    """
    return str(
        lp.loop(
            idx_var="blk_id",
            lower_bound="0",
            upper_bound=f"static_cast<std::ptrdiff_t>({count})",
            increment="1",
            pragma="",
            loop_body=loop_body,
        )
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
