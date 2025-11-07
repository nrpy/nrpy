"""
Generates the C function for unloading a single k-d tree snapshot file.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def unload_kdtree_snapshot() -> None:
    """
    Generate and register the C function for unmapping a k-d tree file.

    This builder generates the C function `unload_kdtree_snapshot()`, which
    uses the `munmap` system call to release the memory mapping of a
    .kdtree.bin file.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h", "<sys/mman.h>"]
    desc = r"""@brief Unloads a memory-mapped k-d tree snapshot.
    @param tree Pointer to the CustomKDTree struct to be unmapped.
    """
    name = "unload_kdtree_snapshot"
    params = "CustomKDTree *tree"
    body = r"""
    if (tree->original_mmap_ptr != NULL) {
        munmap(tree->original_mmap_ptr, tree->file_size);
        tree->original_mmap_ptr = NULL;
    }
    """
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )