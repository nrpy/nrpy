"""
Module generating the C data structures and helpers for the Time Slot Manager.

Provides a lock-free Arena Allocator designed for Host-side orchestration 
of batched photon trajectories based on their physical coordinate time.
"""
from nrpy.infrastructures.BHaH import BHaH_defines_h as Bdefines_h

def time_slot_manager_helpers() -> None:
    """
    Registers the TimeSlotManager struct and associated inline helper functions 
    within the BHaH macro definitions file.
    
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
        if (t < tsm->t_min || t >= tsm->t_max) return -1;
        return (int)floor((t - tsm->t_min) / tsm->delta_t_slot);
    }

    // @brief Atomically pushes a photon's global index into the appropriate temporal bin.
    static inline void slot_add_photon(TimeSlotManager *tsm, int slot_idx, long int photon_idx) {
        if (photon_idx >= tsm->max_capacity) return; 

        // Variable tracking the previous head of the linked list during the atomic swap.
        long int old_head;
        #pragma omp atomic capture
        {
            old_head = tsm->slot_heads[slot_idx];
            tsm->slot_heads[slot_idx] = photon_idx;
        }
        tsm->photon_next_ptrs[photon_idx] = old_head;
        #pragma omp atomic
        tsm->slot_counts[slot_idx]++;
    }

    // @brief Pops a specific number of photons from a designated time slot into a contiguous buffer.
    static inline void slot_remove_chunk(TimeSlotManager *tsm, int slot_idx, long int *chunk_buffer, long int chunk_size) {
        for(long int i = 0; i < chunk_size; i++) {
            // The current lead photon index residing in the requested temporal slot.
            long int current_head = tsm->slot_heads[slot_idx];
            if (current_head == -1) break; 
            chunk_buffer[i] = current_head;
            tsm->slot_heads[slot_idx] = tsm->photon_next_ptrs[current_head];
            tsm->slot_counts[slot_idx]--;
        }
    }
    """
    Bdefines_h.register_BHaH_defines("time_slot_manager", portable_tsm)