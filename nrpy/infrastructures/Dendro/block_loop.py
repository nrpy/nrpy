# nrpy/infrastructures/Dendro/block_loop.py
"""
Dendro numerical block-loop emission built on the generic NRPy loop helper.

The block loop iterates the local block list supplied by Dendro and invokes
the registered per-block CFunction for each block.  The loop is emitted by
NRPy; the Dendro runtime only provides the block list and num_blocks.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.helpers.loop as lp


def block_loop(loop_body: str, num_blocks: str = "numBlocks") -> str:
    """
    Emit a Dendro numerical block loop around a loop body.

    :param loop_body: C code executed once per block; has access to ``blk``.
    :param num_blocks: C expression for the number of local blocks.
    :return: The generated loop C code string.

    Doctests:
    >>> print(block_loop("rhs_block(blk);", num_blocks="NBLK"))
    for (int blk = 0; blk < static_cast<std::ptrdiff_t>(NBLK); blk++) {
    rhs_block(blk);
    } // END LOOP: for blk over [0, static_cast<std::ptrdiff_t>(NBLK))
    <BLANKLINE>
    """
    return str(
        lp.loop(
            idx_var="blk",
            lower_bound="0",
            upper_bound=f"static_cast<std::ptrdiff_t>({num_blocks})",
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
