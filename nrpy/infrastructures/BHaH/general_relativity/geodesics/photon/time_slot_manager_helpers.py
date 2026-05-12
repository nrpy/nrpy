# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/time_slot_manager_helpers.py
"""
Defines the C data structures and helpers for the Time Slot Manager.

This module sets up a lock-free arena allocator for the orchestration of batched
photon trajectories based on their physical coordinate time. Binning photons
temporally enforces the split-pipeline architecture required for relativistic ray
tracing. Atomic operations and speculative memory writes guarantee thread safety and
mathematical consistency during concurrent orchestration.

Author: Dalton J. Moone
"""

from nrpy.infrastructures.BHaH import BHaH_defines_h as Bdefines_h


def time_slot_manager_helpers() -> None:
    """
    Register the TimeSlotManager struct and associated inline helper functions within the BHaH macro definitions file.

    >>> time_slot_manager_helpers()
    """
    portable_tsm = r"""
    // Structure representing a discrete temporal arena for binning active photon trajectories.
    typedef struct {
        double t_min; // The absolute minimum physical coordinate time boundary $t_{min}$ for the temporal arena.
        double t_max; // The absolute maximum physical coordinate time boundary $t_{max}$.
        double delta_t_slot; // The uniform temporal width $\Delta t_{slot}$ of each individual bin.
        int num_slots; // The total number of discrete time slots available in the arena.
        long int max_capacity; // The global maximum number of photons the manager can track.
        long int *photon_next_ptrs; // Linked-list pointer array tracking subsequent photons in the same slot.
        long int *slot_heads; // Array of pointers indicating the first photon index residing in each time slot.
        long int *slot_counts; // Array tracking the active number of photons currently residing in each time slot.
    } TimeSlotManager; // END STRUCT: TimeSlotManager

    //==========================================
    // ARENA ALLOCATION & INITIALIZATION
    //==========================================
    //  Initializes the temporal arena structures for photon binning.
    static inline void slot_manager_init(TimeSlotManager *tsm, double t_min, double t_max, double delta_t_slot, long int num_rays) {
        tsm->t_min = t_min;
        tsm->t_max = t_max;
        tsm->delta_t_slot = delta_t_slot;
        tsm->num_slots = (int)ceil((t_max - t_min) / delta_t_slot);
        tsm->max_capacity = num_rays;

        // Allocate Host memory for the pointer structures mapping temporal slots to photon chains.
        tsm->photon_next_ptrs = (long int *)malloc(sizeof(long int) * num_rays);
        tsm->slot_heads = (long int *)malloc(sizeof(long int) * tsm->num_slots);
        tsm->slot_counts = (long int *)malloc(sizeof(long int) * tsm->num_slots);

        // Iterator iterating over the available discrete time slots.
        int i;
        for (i = 0; i < tsm->num_slots; i++) {
            tsm->slot_heads[i] = -1;
            tsm->slot_counts[i] = 0;
        }

        // Iterator iterating over the maximum allocated capacity of photon rays.
        long int j;
        for (j = 0; j < num_rays; j++) {
            tsm->photon_next_ptrs[j] = -1;
        }
    } // END FUNCTION: slot_manager_init

    //==========================================
    // ARENA DEALLOCATION
    //==========================================
    //  Frees all dynamically allocated host-side memory tied to the temporal arena.
    static inline void slot_manager_free(TimeSlotManager *tsm) {
        if (!tsm) return;
        free(tsm->photon_next_ptrs);
        free(tsm->slot_heads);
        free(tsm->slot_counts);
    } // END FUNCTION: slot_manager_free

    //==========================================
    // TEMPORAL INDEXING
    //==========================================
    //  Computes the designated slot index for a given physical coordinate time $t$.
    static inline int slot_get_index(const TimeSlotManager *tsm, double t) {
        if (isnan(t) || t < tsm->t_min || t >= tsm->t_max) return -1;
        return (int)floor((t - tsm->t_min) / tsm->delta_t_slot);
    }

    //==========================================
    // PHOTON REGISTRATION
    //==========================================
    //  Atomically pushes a photon's global index into the appropriate temporal bin.
    // Implements a lock-free linked list insertion using a Compare-and-Swap (CAS) loop.
    static inline void slot_add_photon(TimeSlotManager *tsm, int slot_idx, long int photon_idx) {
        // Prevent index propagation beyond the pre-allocated manager capacity.
        if (photon_idx >= tsm->max_capacity) return;

        // Variable representing the expected head of the linked list during the atomic evaluation.
        long int expected_head;

        //==========================================
        // LOCK-FREE LINKED LIST INSERTION
        //==========================================
        do {
            // Capture the current head to use as the speculative $photon_next_ptrs$ link.
            expected_head = tsm->slot_heads[slot_idx];

            // Enforcing the Link-Before-Publish invariant ensures atomicity across the Host CPU threads.
            tsm->photon_next_ptrs[photon_idx] = expected_head;

            // Atomic Operation: Update the slot head only if it remains equal to the captured $expected\_head$.
        } while (!__sync_bool_compare_and_swap(&tsm->slot_heads[slot_idx], expected_head, photon_idx)); // END DO-WHILE: lock-free linked list insertion

        // Atomic Operation: Increment the active photon count for the designated temporal bin.
        __sync_fetch_and_add(&tsm->slot_counts[slot_idx], 1L);
    } // END FUNCTION: slot_add_photon

    //==========================================
    // PHOTON EXTRACTION
    //==========================================
    //  Atomically pops a specific number of photons from a designated time slot into a contiguous buffer.
    // Utilizes a CAS-based pop mechanism to safely extract the head of the linked list.
    static inline void slot_remove_chunk(TimeSlotManager *tsm, int slot_idx, long int *chunk_buffer, long int chunk_size) {
        // Iterator tracking the current number of successfully extracted photons.
        long int i;
        for(i = 0; i < chunk_size; i++) {
            // Variable representing the current head of the linked list.
            long int current_head;
            // Variable representing the subsequent node to promote to the new head.
            long int next_node;

            //==========================================
            // LOCK-FREE LINKED LIST EXTRACTION
            //==========================================
                        do {
                // Capture the current head of the linked list.
                current_head = tsm->slot_heads[slot_idx];

                // If the bin is empty, terminate the extraction loop.
                if (current_head == -1) return;

                // Identify the subsequent node to promote to the new head.
                next_node = tsm->photon_next_ptrs[current_head];

                // Atomic Operation: Attempt to swap the head with the next node in the sequence.
            } while (!__sync_bool_compare_and_swap(&tsm->slot_heads[slot_idx], current_head, next_node)); // END DO-WHILE: lock-free linked list extraction

            // Populate the contiguous buffer with the successfully extracted index.
            chunk_buffer[i] = current_head;

            // Atomic Operation: Decrement the active count for the temporal bin.
            __sync_fetch_and_sub(&tsm->slot_counts[slot_idx], 1L);
        } // END LOOP: for i over chunk_size
    } // END FUNCTION: slot_remove_chunk
    """
    Bdefines_h.register_BHaH_defines("time_slot_manager", portable_tsm)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
