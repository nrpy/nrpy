"""
Module generating the C data structures and helpers for the Time Slot Manager.

Provides a lock-free Arena Allocator designed for Host-side orchestration
of batched photon trajectories based on their physical coordinate time.
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
        double t_min; // The absolute minimum physical coordinate time boundary for the temporal arena.
        double t_max; // The absolute maximum physical coordinate time boundary.
        double delta_t_slot; // The uniform temporal width of each individual bin.
        int num_slots; // The total number of discrete time slots available in the arena.
        long int max_capacity; // The global maximum number of photons the manager can track.
        long int *photon_next_ptrs; // Linked-list pointer array tracking subsequent photons in the same slot.
        long int *slot_heads; // Array of pointers indicating the first photon index residing in each time slot.
        long int *slot_counts; // Array tracking the active number of photons currently residing in each time slot.
    } TimeSlotManager;

    // @brief Initializes the temporal arena structures for photon binning.
    static inline void slot_manager_init(TimeSlotManager *tsm, double t_min, double t_max, double delta_t_slot, long int num_rays) {
        tsm->t_min = t_min;
        tsm->t_max = t_max;
        tsm->delta_t_slot = delta_t_slot;
        tsm->num_slots = (int)ceil((t_max - t_min) / delta_t_slot);
        tsm->max_capacity = num_rays;

        tsm->photon_next_ptrs = (long int *)malloc(sizeof(long int) * num_rays);
        tsm->slot_heads = (long int *)malloc(sizeof(long int) * tsm->num_slots);
        tsm->slot_counts = (long int *)malloc(sizeof(long int) * tsm->num_slots);

        for (int i = 0; i < tsm->num_slots; i++) {
            tsm->slot_heads[i] = -1;
            tsm->slot_counts[i] = 0;
        }

        for (long int i = 0; i < num_rays; i++) {
            tsm->photon_next_ptrs[i] = -1;
        }
    }

    // @brief Frees all dynamically allocated host-side memory tied to the temporal arena.
    static inline void slot_manager_free(TimeSlotManager *tsm) {
        if (!tsm) return;
        free(tsm->photon_next_ptrs);
        free(tsm->slot_heads);
        free(tsm->slot_counts);
    }

    // @brief Computes the designated slot index for a given physical coordinate time.
    static inline int slot_get_index(const TimeSlotManager *tsm, double t) {
        if (isnan(t) || t < tsm->t_min || t >= tsm->t_max) return -1;
        return (int)floor((t - tsm->t_min) / tsm->delta_t_slot);
    }
    
    // @brief Atomically pushes a photon's global index into the appropriate temporal bin.
    /* Algorithmic Step: Implements a lock-free linked list insertion using a Compare-and-Swap (CAS) loop.
       Hardware Justification: Speculatively writing the next pointer before publishing the head ensures 
       that any concurrent Host thread reading the slot always traverses a mathematically consistent state, 
       preventing list severance. */
    static inline void slot_add_photon(TimeSlotManager *tsm, int slot_idx, long int photon_idx) {
        // Architectural Guard: Prevent index propagation beyond the pre-allocated manager capacity.
        if (photon_idx >= tsm->max_capacity) return;

        // Variable representing the expected head of the linked list during the atomic evaluation.
        long int expected_head;

        // --- LOCK-FREE LINKED LIST INSERTION ---
        do {
            // Capture the current head to use as the speculative $photon\_next\_ptrs$ link.
            expected_head = tsm->slot_heads[slot_idx];

            // Functional Justification: Link the current photon to the captured head BEFORE attempting to publish.
            // This enforces the "Link-Before-Publish" invariant to ensure atomicity across the Host CPU threads.
            tsm->photon_next_ptrs[photon_idx] = expected_head;

            // Atomic Operation: Update the slot head only if it remains equal to the captured $expected\_head$.
        } while (!__sync_bool_compare_and_swap(&tsm->slot_heads[slot_idx], expected_head, photon_idx));

        // Atomic Operation: Increment the active photon count for the designated temporal bin.
        __sync_fetch_and_add(&tsm->slot_counts[slot_idx], 1L);
    }
    
    // @brief Atomically pops a specific number of photons from a designated time slot into a contiguous buffer.
    /* Algorithmic Step: Utilizes a CAS-based "pop" mechanism to safely extract the head of the linked list.
       Hardware Justification: Ensuring the extraction is atomic prevents multiple threads from 
       claiming the same $photon\_idx$ during high-throughput Host orchestration. */
    static inline void slot_remove_chunk(TimeSlotManager *tsm, int slot_idx, long int *chunk_buffer, long int chunk_size) {
        for(long int i = 0; i < chunk_size; i++) {
            long int current_head;
            long int next_node;

            // --- LOCK-FREE LINKED LIST EXTRACTION ---
            do {
                // Capture the current head of the linked list.
                current_head = tsm->slot_heads[slot_idx];

                // If the bin is empty, terminate the extraction loop.
                if (current_head == -1) return;

                // Identify the subsequent node to promote to the new head.
                next_node = tsm->photon_next_ptrs[current_head];

                // Atomic Operation: Attempt to swap the head with the next node in the sequence.
            } while (!__sync_bool_compare_and_swap(&tsm->slot_heads[slot_idx], current_head, next_node));

            // Populate the contiguous buffer with the successfully extracted index.
            chunk_buffer[i] = current_head;

            // Atomic Operation: Decrement the active count for the temporal bin.
            __sync_fetch_and_sub(&tsm->slot_counts[slot_idx], 1L);
        }
    }
    """
    Bdefines_h.register_BHaH_defines("time_slot_manager", portable_tsm)
