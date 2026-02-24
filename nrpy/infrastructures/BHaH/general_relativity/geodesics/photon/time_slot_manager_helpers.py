"""
Generates the C data structures and helper functions for the Time Slot Manager.
Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
Includes lock-free Arena Allocator with device-resident memory offloading.
"""
from nrpy.infrastructures.BHaH import BHaH_defines_h as Bdefines_h

def time_slot_manager_helpers() -> None:
    # Project Singularity-Axiom: Portable Body Wrapper
    portable_tsm = r"""
    #ifdef USE_GPU
    #pragma omp declare target
    #endif

    typedef struct {
        double t_min;
        double t_max;
        double delta_t_slot;
        int num_slots;
        long int max_capacity;
        long int *photon_next_ptrs;
        long int *slot_heads;
        long int *slot_counts;
    } TimeSlotManager;

    static inline void slot_manager_init(TimeSlotManager *tsm, double t_min, double t_max, double delta_t_slot, long int num_rays) {
        tsm->t_min = t_min;
        tsm->t_max = t_max;
        tsm->delta_t_slot = delta_t_slot;
        tsm->num_slots = (int)ceil((t_max - t_min) / delta_t_slot);
        tsm->max_capacity = num_rays;

        // Dual-Architecture Memory Allocation
        #ifdef USE_GPU
            tsm->photon_next_ptrs = (long int *)omp_target_alloc(sizeof(long int) * num_rays, omp_get_default_device());
            tsm->slot_heads = (long int *)omp_target_alloc(sizeof(long int) * tsm->num_slots, omp_get_default_device());
            tsm->slot_counts = (long int *)omp_target_alloc(sizeof(long int) * tsm->num_slots, omp_get_default_device());
        #else
            tsm->photon_next_ptrs = (long int *)malloc(sizeof(long int) * num_rays);
            tsm->slot_heads = (long int *)malloc(sizeof(long int) * tsm->num_slots);
            tsm->slot_counts = (long int *)malloc(sizeof(long int) * tsm->num_slots);
        #endif

        // Note: Slot initialization must follow the architecture's memory rules.
        // On GPU, this loop must be offloaded to populate device-resident arrays.
        #ifdef USE_GPU
            #pragma omp target teams distribute parallel for
        #else
            #pragma omp parallel for
        #endif
        for (int i = 0; i < tsm->num_slots; i++) {
            tsm->slot_heads[i] = -1;
            tsm->slot_counts[i] = 0;
        }

        #ifdef USE_GPU
            #pragma omp target teams distribute parallel for
        #else
            #pragma omp parallel for
        #endif
        for (long int i = 0; i < num_rays; i++) {
            tsm->photon_next_ptrs[i] = -1;
        }
    }

    static inline void slot_manager_free(TimeSlotManager *tsm) {
        if (!tsm) return;
        #ifdef USE_GPU
            omp_target_free(tsm->photon_next_ptrs, omp_get_default_device());
            omp_target_free(tsm->slot_heads, omp_get_default_device());
            omp_target_free(tsm->slot_counts, omp_get_default_device());
        #else
            free(tsm->photon_next_ptrs);
            free(tsm->slot_heads);
            free(tsm->slot_counts);
        #endif
    }

    static inline int slot_get_index(const TimeSlotManager *tsm, double t) {
        if (t < tsm->t_min || t >= tsm->t_max) return -1;
        return (int)floor((t - tsm->t_min) / tsm->delta_t_slot);
    }

    static inline void slot_add_photon(TimeSlotManager *tsm, int slot_idx, long int photon_idx) {
        if (photon_idx >= tsm->max_capacity) return; 

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

    static inline void slot_remove_chunk(TimeSlotManager *tsm, int slot_idx, long int *chunk_buffer, long int chunk_size) {
        for(long int i = 0; i < chunk_size; i++) {
            long int current_head = tsm->slot_heads[slot_idx];
            if (current_head == -1) break; 
            chunk_buffer[i] = current_head;
            tsm->slot_heads[slot_idx] = tsm->photon_next_ptrs[current_head];
            tsm->slot_counts[slot_idx]--;
        }
    }

    #ifdef USE_GPU
    #pragma omp end declare target
    #endif
    """
    Bdefines_h.register_BHaH_defines("time_slot_manager", portable_tsm)