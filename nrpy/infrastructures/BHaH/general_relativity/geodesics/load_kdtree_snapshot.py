"""
Generates the C function for loading a single k-d tree snapshot file.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def load_kdtree_snapshot() -> None:
    """
    Generate and register the C function for memory-mapping a k-d tree file.

    This builder generates the C function `load_kdtree_snapshot()`, which uses
    the `mmap` system call to efficiently map a .kdtree.bin file into memory.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = [
        "BHaH_defines.h",
        "<stdio.h>",
        "<stdlib.h>",
        "<sys/mman.h>",
        "<sys/stat.h>",
        "<fcntl.h>",
        "<unistd.h>",
    ]
    desc = r"""@brief Loads a .kdtree.bin snapshot file into memory using mmap.
    @details Memory-mapping is used for performance, allowing the OS to handle
             on-demand paging of the file data from disk.
    @param[in]  filename The path to the .kdtree.bin file.
    @param[out] tree     Pointer to the CustomKDTree struct to be populated.
    @return 0 on success, -1 on failure.
    """
    name = "load_kdtree_snapshot"
    params = "const char *filename, CustomKDTree *tree"
    cfunc_type = "int"
    body = r"""
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
        perror("Error opening k-d tree file");
        return -1; // Failure
    }

    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        perror("Error getting file size");
        close(fd);
        return -1;
    }
    tree->file_size = sb.st_size;

    void *mapped_mem = mmap(NULL, tree->file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped_mem == MAP_FAILED) {
        perror("Error memory-mapping the file");
        close(fd);
        return -1;
    }
    close(fd); // File descriptor no longer needed after mmap

    tree->original_mmap_ptr = mapped_mem;
    char *current_ptr = (char *)mapped_mem;

    // Read header information from the start of the memory-mapped region.
    tree->num_particles = *(uint64_t *)current_ptr;
    current_ptr += sizeof(uint64_t);
    tree->dimensions = *(uint64_t *)current_ptr;
    current_ptr += sizeof(uint64_t);

    // Set pointers to the data payloads within the mapped region.
    tree->node_metadata = (int32_t *)current_ptr;
    current_ptr += sizeof(int32_t) * tree->num_particles;
    tree->particle_data = (MassiveParticle *)current_ptr;

    return 0; // Success
    """
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        cfunc_type=cfunc_type,
    )