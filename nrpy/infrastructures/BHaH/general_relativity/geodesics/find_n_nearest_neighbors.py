"""
Generates the C engine for k-d tree nearest neighbor search.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def find_n_nearest_neighbors() -> None:
    """
    Generate and register the C functions for recursive k-NN search.

    This builder generates the C function `find_n_nearest_neighbors`, which
    acts as an entry point to a high-performance, recursive search algorithm.
    The recursive worker, `search_recursive`, is defined as a `static` helper
    function in a `prefunc` block. This version is performance-optimized with
    `__builtin_prefetch` to hide memory latency.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h", "<math.h>", "<stdio.h>"]

    # The prefunc contains all the static helper functions for the search.
    prefunc = r"""
// Helper to initialize the WinnersCircle struct, which tracks the best N candidates.
static void wc_init(WinnersCircle *wc, int n_wanted) {
    wc->n_wanted = n_wanted;
    wc->count = 0;
    for (int i = 0; i < n_wanted; ++i) {
        wc->indices[i] = -1;
        wc->sq_distances[i] = 1e300; // Initialize with a very large number
    }
}

// Helper to add a candidate to the WinnersCircle, maintaining sorted order.
static void wc_add(WinnersCircle *wc, int index, double sq_dist) {
    if (wc->count < wc->n_wanted) {
        wc->indices[wc->count] = index;
        wc->sq_distances[wc->count] = sq_dist;
        wc->count++;
    } else if (sq_dist < wc->sq_distances[wc->n_wanted - 1]) {
        wc->indices[wc->n_wanted - 1] = index;
        wc->sq_distances[wc->n_wanted - 1] = sq_dist;
    } else {
        return; // Not a winner
    }

    // Simple insertion sort to keep the winners list sorted.
    for (int i = wc->count - 1; i > 0; --i) {
        if (wc->sq_distances[i] < wc->sq_distances[i - 1]) {
            double temp_d = wc->sq_distances[i];
            int temp_i = wc->indices[i];
            wc->sq_distances[i] = wc->sq_distances[i - 1];
            wc->indices[i] = wc->indices[i - 1];
            wc->sq_distances[i - 1] = temp_d;
            wc->indices[i - 1] = temp_i;
        }
    }
}

// The core recursive search function, optimized with prefetching.
static void search_recursive(const CustomKDTree *tree, const double query_pos[3], int current_idx, WinnersCircle *wc) {
    if (current_idx < 0 || current_idx >= (int)tree->num_particles) {
        return;
    }

    const MassiveParticle *pivot = &tree->particle_data[current_idx];
    const int split_axis = tree->node_metadata[current_idx];

    const double dx = query_pos[0] - pivot->pos[0];
    const double dy = query_pos[1] - pivot->pos[1];
    const double dz = query_pos[2] - pivot->pos[2];
    const double dist_sq = dx*dx + dy*dy + dz*dz;
    wc_add(wc, current_idx, dist_sq);

    if (split_axis == -1) { // Leaf node
        return;
    }

    const double axis_dist = query_pos[split_axis] - pivot->pos[split_axis];
    const int good_side_idx = (axis_dist < 0) ? (2 * current_idx + 1) : (2 * current_idx + 2);
    const int bad_side_idx = (axis_dist < 0) ? (2 * current_idx + 2) : (2 * current_idx + 1);

    // *** PERFORMANCE OPTIMIZATION ***
    // Issue prefetch instructions for the data of the child nodes. This hints to the
    // CPU to start loading this memory into the cache while we process the current node.
    if (good_side_idx < (int)tree->num_particles) {
        __builtin_prefetch(&tree->particle_data[good_side_idx], 0, 1);
    }
    if (bad_side_idx < (int)tree->num_particles) {
        __builtin_prefetch(&tree->particle_data[bad_side_idx], 0, 1);
    }

    // Recurse down the "good" side of the tree.
    search_recursive(tree, query_pos, good_side_idx, wc);

    // If the search sphere overlaps with the splitting plane, we must also search the "bad" side.
    const double search_radius_sq = wc->sq_distances[wc->n_wanted - 1];
    if (axis_dist * axis_dist < search_radius_sq) {
        search_recursive(tree, query_pos, bad_side_idx, wc);
    }
}
"""
    desc = r"""@brief Finds the N nearest neighbors to a query point in a k-d tree.
    @details This function is the public entry point for the k-NN search. It
             initializes a WinnersCircle and starts the recursive search from
             the root of the tree.
    @param[in]  tree             Pointer to the loaded k-d tree data.
    @param[in]  query_pos        The 3D Cartesian point to search around.
    @param[in]  n_neighbors      The number of nearest neighbors to find.
    @param[out] neighbor_results An array to be filled with the found neighbor data.
    """
    name = "find_n_nearest_neighbors"
    params = "const CustomKDTree *tree, const double query_pos[3], int n_neighbors, MassiveParticle *neighbor_results"

    body = r"""
    if (n_neighbors > MAX_NEIGHBORS) {
        fprintf(stderr, "Error: Requested more neighbors than MAX_NEIGHBORS.\\n");
        return;
    }

    WinnersCircle wc;
    wc_init(&wc, n_neighbors);

    // Start the recursive search from the root (index 0).
    search_recursive(tree, query_pos, 0, &wc);

    // Copy the results from the WinnersCircle into the output array.
    for (int i = 0; i < wc.count; ++i) {
        neighbor_results[i] = tree->particle_data[wc.indices[i]];
    }
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )
