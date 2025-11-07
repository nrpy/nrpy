"""
Generates the C data structures and helper functions for the Time Slot Manager.

Author: Dalton J. Moone
"""

from nrpy.infrastructures.BHaH import BHaH_defines_h as Bdefines_h


def time_slot_manager_helpers() -> None:
    """
    Generate and register the TimeSlotManager C code for injection into BHaH_defines.h.

    This function generates the C code for the data structures (`PhotonList`,
    `TimeSlotManager`) and a suite of `static inline` helper functions that
    manage the time slotting algorithm. This includes initialization, memory
    management, and the core logic for adding and removing photons from slots.
    """
    # The entire C code block to be injected into the header.
    # This is algorithmic, not symbolic.
    c_code_for_header = r"""
// =============================================
// NRPy-Generated Time Slot Manager
// =============================================

// --- Data Structures ---
// A dynamically-sized list to hold photon indices for a single time slot.
typedef struct {
    long int *photons;
    long int count;
    long int capacity;
} PhotonList;

// The main manager struct, containing an array of all time slots.
typedef struct {
    double t_min, t_max, delta_t_slot;
    int num_slots;
    PhotonList *slots;
} TimeSlotManager;

// --- Function Definitions (static inline to live in a header) ---

// Initializes the time slot manager, allocating memory for all slots.
static inline void slot_manager_init(TimeSlotManager *tsm, double t_min, double t_max, double delta_t_slot) {
    tsm->t_min = t_min;
    tsm->t_max = t_max;
    tsm->delta_t_slot = delta_t_slot;
    tsm->num_slots = (int)ceil((t_max - t_min) / delta_t_slot);
    if (tsm->num_slots <= 0) { fprintf(stderr, "Error: Invalid TimeSlotManager dimensions.\n"); exit(1); }
    tsm->slots = (PhotonList *)malloc(sizeof(PhotonList) * tsm->num_slots);
    if (tsm->slots == NULL) { fprintf(stderr, "Error: Failed to allocate memory for time slots.\n"); exit(1); }
    for (int i = 0; i < tsm->num_slots; i++) {
        tsm->slots[i].photons = (long int *)malloc(sizeof(long int) * 16); // Initial capacity
        if (tsm->slots[i].photons == NULL) { fprintf(stderr, "Error: Failed to allocate memory for a slot's photon list.\n"); exit(1); }
        tsm->slots[i].count = 0;
        tsm->slots[i].capacity = 16;
    }
}

// Frees all memory associated with the time slot manager.
static inline void slot_manager_free(TimeSlotManager *tsm) {
    for (int i = 0; i < tsm->num_slots; i++) { free(tsm->slots[i].photons); }
    free(tsm->slots);
}

// A fast hash function that maps a coordinate time 't' to a slot index.
static inline int slot_get_index(const TimeSlotManager *tsm, double t) {
    if (t < tsm->t_min || t >= tsm->t_max) return -1; // Out of bounds
    return (int)floor((t - tsm->t_min) / tsm->delta_t_slot);
}

// Adds a photon's index to a slot, resizing the list if necessary.
static inline void slot_add_photon(PhotonList *slot, long int photon_idx) {
    if (slot->count >= slot->capacity) {
        slot->capacity *= 2;
        slot->photons = (long int *)realloc(slot->photons, sizeof(long int) * slot->capacity);
        if (slot->photons == NULL) { fprintf(stderr, "Error: Failed to reallocate memory for a slot's photon list.\n"); exit(1); }
    }
    slot->photons[slot->count++] = photon_idx;
}

// Efficiently removes a batch of photons from the front of a slot's list.
static inline void slot_remove_chunk(PhotonList *slot, long int *chunk_buffer, long int chunk_size) {
    for (long int i = 0; i < chunk_size; ++i) {
        chunk_buffer[i] = slot->photons[i];
    }
    // Use memmove for safe overlapping memory copy to shift remaining elements.
    memmove(slot->photons, slot->photons + chunk_size, (slot->count - chunk_size) * sizeof(long int));
    slot->count -= chunk_size;
}
"""

    # Register this C code block to be injected into the BHaH_defines.h header file.
    Bdefines_h.register_BHaH_defines("time_slot_manager", c_code_for_header)