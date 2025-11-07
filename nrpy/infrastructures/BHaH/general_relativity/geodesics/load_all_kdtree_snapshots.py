"""
Generates the C orchestrator for loading all k-d tree snapshots.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def load_all_kdtree_snapshots() -> None:
    """
    Generate and register the C orchestrator for finding, sorting, and loading
    all k-d tree snapshot files from a directory into memory.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "stdio.h",
        "stdlib.h",
        "<dirent.h>",
        "<string.h>",
    ]
    desc = r"""@brief Finds, sorts, and loads all .kdtree.bin files from the snapshot directory.
    @details This function encapsulates the logic for memory-mapping the k-d tree
             data needed for disk intersection checks.
    @param[in]  commondata        Pointer to the commondata struct with runtime parameters.
    @param[out] kdtree_snapshots  Pointer to be allocated and filled with snapshot data.
    @param[out] snapshot_times    Pointer to be allocated and filled with snapshot times.
    @return The total number of snapshots successfully loaded.
    """
    name = "load_all_kdtree_snapshots"
    cfunc_type = "int"
    params = """
    const commondata_struct *restrict commondata,
    CustomKDTree **kdtree_snapshots,
    double **snapshot_times
    """

    # The body is algorithmic, not symbolic, so it is defined as a raw C string.
    body = r"""
    const char* snapshot_dir_path = "../processed_snapshots"; // Assumes a relative path
    printf("Loading k-d tree snapshots from '%s'...\n", snapshot_dir_path);
    DIR *dir;
    struct dirent *ent;
    int num_snapshots = 0;

    // First pass: count the number of snapshot files to allocate memory.
    if ((dir = opendir(snapshot_dir_path)) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (strstr(ent->d_name, ".kdtree.bin") != NULL) {
                num_snapshots++;
            }
        }
        closedir(dir);
    } else {
        perror("Could not open snapshot directory");
        exit(1);
    }

    if (num_snapshots == 0) {
        fprintf(stderr, "Warning: No .kdtree.bin snapshot files found in '%s'. Disk intersection will be disabled.\n", snapshot_dir_path);
        *kdtree_snapshots = NULL;
        *snapshot_times = NULL;
        return 0;
    }

    // Second pass: store filenames and sort them.
    char **filenames = (char **)malloc(sizeof(char *) * num_snapshots);
    if (filenames == NULL) { exit(1); }
    dir = opendir(snapshot_dir_path);
    int count = 0;
    while ((ent = readdir(dir)) != NULL) {
        if (strstr(ent->d_name, ".kdtree.bin") != NULL) {
            filenames[count] = strdup(ent->d_name);
            count++;
        }
    }
    closedir(dir);
    qsort(filenames, num_snapshots, sizeof(char *), compare_filenames);

    // Allocate memory for the snapshot data and timestamps.
    *kdtree_snapshots = (CustomKDTree *)malloc(sizeof(CustomKDTree) * num_snapshots);
    *snapshot_times = (double *)malloc(sizeof(double) * num_snapshots);
    if (*kdtree_snapshots == NULL || *snapshot_times == NULL) { exit(1); }

    // Third pass: load each sorted file into memory.
    for (int i = 0; i < num_snapshots; ++i) {
        char filepath[512];
        snprintf(filepath, sizeof(filepath), "%s/%s", snapshot_dir_path, filenames[i]);
        if (load_kdtree_snapshot(filepath, &(*kdtree_snapshots)[i]) != 0) {
            fprintf(stderr, "Error: Failed to load snapshot %s\n", filepath);
            exit(1);
        }
        int snapshot_index;
        sscanf(filenames[i], "mass_blueprint_t_%d.kdtree.bin", &snapshot_index);
        (*snapshot_times)[i] = (double)snapshot_index * commondata->mass_snapshot_every_t;
        free(filenames[i]);
    }
    free(filenames);
    printf("Successfully loaded and sorted %d snapshots.\n", num_snapshots);

    return num_snapshots;
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )