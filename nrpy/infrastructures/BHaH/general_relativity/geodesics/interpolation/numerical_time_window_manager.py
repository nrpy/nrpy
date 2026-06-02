# nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/numerical_time_window_manager.py
"""
Define C helpers for numerical-spacetime time-window mmap management.

The NumericalTimeWindowManager owns the active mapped group of adjacent 3D-grid
payloads in one trusted combined numerical-spacetime binary container. It
depends on shared TimeSlotManager definitions for slot time bounds, but keeps
mmap and numerical file state out of analytic builds.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from nrpy.infrastructures.BHaH import BHaH_defines_h as Bdefines_h


def numerical_time_window_manager() -> None:
    """
    Register NumericalTimeWindowManager helpers in BHaH_defines.h.

    The time_slot_manager() registration must run before this function so the
    generated C has TimeSlotManager and slot bound helpers available first.
    """
    numerical_time_window_manager_c_code = r"""
    //==========================================
    // NUMERICAL TIME-WINDOW MANAGER
    //==========================================
    // Lightweight trusted-input mmap window over one combined numerical-spacetime container.
#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

#define NUMERICAL_TIME_WINDOW_MANAGER_FIXED_HEADER_BYTES 4096ULL
#define NUMERICAL_TIME_WINDOW_MANAGER_SLICE_TABLE_ENTRY_BYTES 128ULL

#define NUMERICAL_TIME_WINDOW_MANAGER_HEADER_SLICE_TABLE_OFFSET 64U
#define NUMERICAL_TIME_WINDOW_MANAGER_HEADER_NUM_TIME_SLICES 80U
#define NUMERICAL_TIME_WINDOW_MANAGER_HEADER_FIRST_PAYLOAD_OFFSET 104U
#define NUMERICAL_TIME_WINDOW_MANAGER_HEADER_PAYLOAD_BYTES_PER_SLICE 112U
#define NUMERICAL_TIME_WINDOW_MANAGER_HEADER_PAYLOAD_STRIDE_BYTES 120U
#define NUMERICAL_TIME_WINDOW_MANAGER_HEADER_POINT_RECORD_COUNT 144U
#define NUMERICAL_TIME_WINDOW_MANAGER_HEADER_POINT_RECORD_BYTES 160U
#define NUMERICAL_TIME_WINDOW_MANAGER_SLICE_ENTRY_TIME 8U

#define NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS 0
#define NUMERICAL_TIME_WINDOW_MANAGER_ERROR 1

    // Owns the currently mapped group of adjacent 3D-grid payloads.
    typedef struct {
        int fd; // Open read-only file descriptor for the combined numerical .bin file.
        long int page_size; // Runtime page size used for mmap offset alignment.

        uint64_t num_time_slices; // Number of stored 3D grids.
        uint64_t slice_table_offset; // Absolute byte offset of the slice table.
        uint64_t first_payload_offset; // Absolute byte offset of slice 0 payload.
        uint64_t payload_bytes_per_slice; // Valid bytes in one 3D-grid payload.
        uint64_t payload_stride_bytes; // Byte stride between adjacent slice payloads.
        uint64_t point_record_count; // Number of point records in one payload.
        uint32_t point_record_bytes; // Bytes in one point record.

        double time_start; // Coordinate time of slice 0.
        double grid_time_spacing; // Uniform coordinate-time spacing between 3D grids.
        double inv_grid_time_spacing; // Reciprocal grid spacing.

        int temporal_interp_half_width; // Centered temporal interpolation half-width n.
        int temporal_interp_num_points; // Number of temporal interpolation points, 2n + 1.
        double max_backward_dt_lookahead; // Maximum backward coordinate-time RKF45 lookahead from commondata.

        uint64_t mapped_first_slice; // First mapped slice, inclusive.
        uint64_t mapped_last_slice_exclusive; // One past the final mapped slice.
        uint64_t mapped_file_offset; // Page-aligned file offset passed to mmap().
        size_t mapped_length_bytes; // Length passed to mmap()/munmap().
        const unsigned char *mapped_base; // Pointer returned by mmap(), or NULL.
        const unsigned char *active_payload_base; // Pointer to mapped_first_slice payload.
    } NumericalTimeWindowManager; // END STRUCT: NumericalTimeWindowManager

    //==========================================
    // NUMERICAL WINDOW LOW-LEVEL READ HELPERS
    //==========================================
    // These helpers assume the trusted combined container uses native little-endian f64/u64/u32 fields.
    static inline uint64_t numerical_time_window_manager_load_u64(const unsigned char *bytes) {
        uint64_t value;
        memcpy(&value, bytes, sizeof(uint64_t));
        return value;
    } // END FUNCTION: numerical_time_window_manager_load_u64

    static inline uint32_t numerical_time_window_manager_load_u32(const unsigned char *bytes) {
        uint32_t value;
        memcpy(&value, bytes, sizeof(uint32_t));
        return value;
    } // END FUNCTION: numerical_time_window_manager_load_u32

    static inline double numerical_time_window_manager_load_f64(const unsigned char *bytes) {
        double value;
        memcpy(&value, bytes, sizeof(double));
        return value;
    } // END FUNCTION: numerical_time_window_manager_load_f64

    /**
     * Read exactly nbytes from one absolute file offset.
     *
     * @param fd Open file descriptor.
     * @param[out] buffer Destination buffer.
     * @param nbytes Number of bytes to read.
     * @param offset Absolute file offset.
     * @return NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS or NUMERICAL_TIME_WINDOW_MANAGER_ERROR.
     */
    static inline int numerical_time_window_manager_read_exact_at(int fd, void *buffer, size_t nbytes, uint64_t offset) {
        unsigned char *dst = (unsigned char *)buffer;
        size_t bytes_read = 0U;

        while (bytes_read < nbytes) {
            const ssize_t nread = pread(fd, dst + bytes_read, nbytes - bytes_read, (off_t)(offset + bytes_read));
            if (nread <= 0) {
                return NUMERICAL_TIME_WINDOW_MANAGER_ERROR;
            } // END IF: metadata read failed before requested byte count was reached
            bytes_read += (size_t)nread;
        } // END WHILE: read requested metadata bytes
        return NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS;
    } // END FUNCTION: numerical_time_window_manager_read_exact_at

    //==========================================
    // NUMERICAL WINDOW LIFETIME MANAGEMENT
    //==========================================
    /**
     * Reset a numerical time-window manager to an inert state.
     *
     * @param[in,out] ntwm Numerical time-window manager to reset.
     */
    static inline void numerical_time_window_manager_set_inert(NumericalTimeWindowManager *ntwm) {
        ntwm->fd = -1;
        ntwm->page_size = 0;
        ntwm->num_time_slices = 0ULL;
        ntwm->slice_table_offset = 0ULL;
        ntwm->first_payload_offset = 0ULL;
        ntwm->payload_bytes_per_slice = 0ULL;
        ntwm->payload_stride_bytes = 0ULL;
        ntwm->point_record_count = 0ULL;
        ntwm->point_record_bytes = 0U;
        ntwm->time_start = 0.0;
        ntwm->grid_time_spacing = 0.0;
        ntwm->inv_grid_time_spacing = 0.0;
        ntwm->temporal_interp_half_width = 0;
        ntwm->temporal_interp_num_points = 0;
        ntwm->max_backward_dt_lookahead = 0.0;
        ntwm->mapped_first_slice = 0ULL;
        ntwm->mapped_last_slice_exclusive = 0ULL;
        ntwm->mapped_file_offset = 0ULL;
        ntwm->mapped_length_bytes = 0U;
        ntwm->mapped_base = NULL;
        ntwm->active_payload_base = NULL;
    } // END FUNCTION: numerical_time_window_manager_set_inert

    /**
     * Close the active mmap window.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     */
    static inline void numerical_time_window_manager_close_mmap(NumericalTimeWindowManager *ntwm) {
        if (ntwm->mapped_base != NULL) {
            munmap((void *)ntwm->mapped_base, ntwm->mapped_length_bytes);
        } // END IF: active mmap window exists
        ntwm->mapped_first_slice = 0ULL;
        ntwm->mapped_last_slice_exclusive = 0ULL;
        ntwm->mapped_file_offset = 0ULL;
        ntwm->mapped_length_bytes = 0U;
        ntwm->mapped_base = NULL;
        ntwm->active_payload_base = NULL;
    } // END FUNCTION: numerical_time_window_manager_close_mmap

    /**
     * Release all resources owned by the numerical time-window manager.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     */
    static inline void numerical_time_window_manager_free(NumericalTimeWindowManager *ntwm) {
        numerical_time_window_manager_close_mmap(ntwm);
        if (ntwm->fd >= 0) {
            close(ntwm->fd);
        } // END IF: combined numerical container was open
        numerical_time_window_manager_set_inert(ntwm);
    } // END FUNCTION: numerical_time_window_manager_free

    /**
     * Initialize a trusted-input numerical time-window manager.
     *
     * The caller must call numerical_time_window_manager_set_inert() before
     * the first initialization. Reinitialization releases the old mapping and
     * file before opening the replacement container.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     * @param[in] combined_path Path to the combined numerical .bin container.
     * @param[in] commondata Common runtime parameters.
     * @param temporal_interp_half_width Centered temporal interpolation half-width n.
     * @return NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS or NUMERICAL_TIME_WINDOW_MANAGER_ERROR.
     *
     * @note The forward-lookahead distance is read from
     * commondata->rkf45_max_delta_t.
     */
    static inline int numerical_time_window_manager_init(
        NumericalTimeWindowManager *ntwm,
        const char *combined_path,
        const commondata_struct *restrict commondata,
        const int temporal_interp_half_width) {
        unsigned char header_bytes[NUMERICAL_TIME_WINDOW_MANAGER_FIXED_HEADER_BYTES];
        unsigned char time_bytes[sizeof(double)];

        numerical_time_window_manager_free(ntwm);
        ntwm->temporal_interp_half_width = temporal_interp_half_width;
        ntwm->temporal_interp_num_points = 2 * temporal_interp_half_width + 1;
        ntwm->max_backward_dt_lookahead = commondata->rkf45_max_delta_t;
        ntwm->page_size = sysconf(_SC_PAGESIZE);

        ntwm->fd = open(combined_path, O_RDONLY);
        if (ntwm->fd < 0) {
            return NUMERICAL_TIME_WINDOW_MANAGER_ERROR;
        } // END IF: combined numerical container could not be opened
        if (numerical_time_window_manager_read_exact_at(
                ntwm->fd, header_bytes, (size_t)NUMERICAL_TIME_WINDOW_MANAGER_FIXED_HEADER_BYTES, 0ULL) !=
            NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS) {
            numerical_time_window_manager_free(ntwm);
            return NUMERICAL_TIME_WINDOW_MANAGER_ERROR;
        } // END IF: fixed header could not be read

        ntwm->slice_table_offset =
            numerical_time_window_manager_load_u64(header_bytes + NUMERICAL_TIME_WINDOW_MANAGER_HEADER_SLICE_TABLE_OFFSET);
        ntwm->num_time_slices =
            numerical_time_window_manager_load_u64(header_bytes + NUMERICAL_TIME_WINDOW_MANAGER_HEADER_NUM_TIME_SLICES);
        ntwm->first_payload_offset =
            numerical_time_window_manager_load_u64(header_bytes + NUMERICAL_TIME_WINDOW_MANAGER_HEADER_FIRST_PAYLOAD_OFFSET);
        ntwm->payload_bytes_per_slice =
            numerical_time_window_manager_load_u64(header_bytes + NUMERICAL_TIME_WINDOW_MANAGER_HEADER_PAYLOAD_BYTES_PER_SLICE);
        ntwm->payload_stride_bytes =
            numerical_time_window_manager_load_u64(header_bytes + NUMERICAL_TIME_WINDOW_MANAGER_HEADER_PAYLOAD_STRIDE_BYTES);
        ntwm->point_record_count =
            numerical_time_window_manager_load_u64(header_bytes + NUMERICAL_TIME_WINDOW_MANAGER_HEADER_POINT_RECORD_COUNT);
        ntwm->point_record_bytes =
            numerical_time_window_manager_load_u32(header_bytes + NUMERICAL_TIME_WINDOW_MANAGER_HEADER_POINT_RECORD_BYTES);

        if (numerical_time_window_manager_read_exact_at(
                ntwm->fd, time_bytes, sizeof(double),
                ntwm->slice_table_offset + NUMERICAL_TIME_WINDOW_MANAGER_SLICE_ENTRY_TIME) !=
            NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS) {
            numerical_time_window_manager_free(ntwm);
            return NUMERICAL_TIME_WINDOW_MANAGER_ERROR;
        } // END IF: first slice time could not be read
        ntwm->time_start = numerical_time_window_manager_load_f64(time_bytes);

        if (numerical_time_window_manager_read_exact_at(
                ntwm->fd, time_bytes, sizeof(double),
                ntwm->slice_table_offset + NUMERICAL_TIME_WINDOW_MANAGER_SLICE_TABLE_ENTRY_BYTES +
                    NUMERICAL_TIME_WINDOW_MANAGER_SLICE_ENTRY_TIME) != NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS) {
            numerical_time_window_manager_free(ntwm);
            return NUMERICAL_TIME_WINDOW_MANAGER_ERROR;
        } // END IF: second slice time could not be read
        ntwm->grid_time_spacing = numerical_time_window_manager_load_f64(time_bytes) - ntwm->time_start;
        ntwm->inv_grid_time_spacing = 1.0 / ntwm->grid_time_spacing;
        return NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS;
    } // END FUNCTION: numerical_time_window_manager_init

    //==========================================
    // NUMERICAL WINDOW SLOT RANGE SELECTION
    //==========================================
    /**
     * Compute the required time-slice range for one TimeSlotManager slot.
     *
     * @param[in] ntwm Initialized numerical time-window manager.
     * @param[in] tsm Base time-slot manager.
     * @param slot_idx Slot index.
     * @param[out] first_slice First required slice, inclusive.
     * @param[out] last_slice_exclusive One past the final required slice.
     * @return NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS or NUMERICAL_TIME_WINDOW_MANAGER_ERROR.
     */
    static inline int numerical_time_window_manager_required_grid_range(
        const NumericalTimeWindowManager *ntwm,
        const TimeSlotManager *tsm,
        const int slot_idx,
        uint64_t *first_slice,
        uint64_t *last_slice_exclusive) {
        const double t_lower = slot_lower_time(tsm, slot_idx);
        const double t_upper = slot_upper_time(tsm, slot_idx);
        const double query_min = t_lower - ntwm->max_backward_dt_lookahead;
        const double query_max = t_upper;
        const double scaled_min = (query_min - ntwm->time_start) * ntwm->inv_grid_time_spacing;
        const double scaled_max = (query_max - ntwm->time_start) * ntwm->inv_grid_time_spacing;
        const long int center_first = (long int)floor(scaled_min);
        const long int center_last = (long int)ceil(scaled_max);
        const long int halo = (long int)ntwm->temporal_interp_half_width + 1L;
        const long int required_first = center_first - halo;
        const long int required_last_exclusive = center_last + halo + 1L;

        if (required_first < 0L || required_last_exclusive > (long int)ntwm->num_time_slices ||
            required_last_exclusive <= required_first) {
            *first_slice = 0ULL;
            *last_slice_exclusive = 0ULL;
            return NUMERICAL_TIME_WINDOW_MANAGER_ERROR;
        } // END IF: requested temporal interpolation window is outside the trusted file

        *first_slice = (uint64_t)required_first;
        *last_slice_exclusive = (uint64_t)required_last_exclusive;
        return NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS;
    } // END FUNCTION: numerical_time_window_manager_required_grid_range

    /**
     * Map exactly one slot's required contiguous time-slice payload window.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     * @param[in] tsm Base time-slot manager.
     * @param slot_idx Slot index.
     * @return NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS or NUMERICAL_TIME_WINDOW_MANAGER_ERROR.
     *
     * @pre No interpolation reader may hold pointers returned by
     * numerical_time_window_manager_grid_ptr during this call.
     */
    static inline int numerical_time_window_manager_mmap_for_slot(
        NumericalTimeWindowManager *ntwm,
        const TimeSlotManager *tsm,
        const int slot_idx) {
        uint64_t first_slice = 0ULL;
        uint64_t last_slice_exclusive = 0ULL;
        if (numerical_time_window_manager_required_grid_range(ntwm, tsm, slot_idx, &first_slice, &last_slice_exclusive) !=
            NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS) {
            return NUMERICAL_TIME_WINDOW_MANAGER_ERROR;
        } // END IF: requested slot window is outside the trusted file

        if (ntwm->mapped_base != NULL && ntwm->mapped_first_slice == first_slice &&
            ntwm->mapped_last_slice_exclusive == last_slice_exclusive) {
            return NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS;
        } // END IF: requested slot window is already active

        const uint64_t requested_offset = ntwm->first_payload_offset + first_slice * ntwm->payload_stride_bytes;
        const uint64_t last_payload_offset =
            ntwm->first_payload_offset + (last_slice_exclusive - 1ULL) * ntwm->payload_stride_bytes;
        const uint64_t requested_end = last_payload_offset + ntwm->payload_bytes_per_slice;
        const uint64_t page_size = (uint64_t)ntwm->page_size;
        const uint64_t map_offset = requested_offset - (requested_offset % page_size);
        const size_t map_length = (size_t)(requested_end - map_offset);

        numerical_time_window_manager_close_mmap(ntwm);
        void *new_mapping = mmap(NULL, map_length, PROT_READ, MAP_PRIVATE, ntwm->fd, (off_t)map_offset);
        if (new_mapping == MAP_FAILED) {
            return NUMERICAL_TIME_WINDOW_MANAGER_ERROR;
        } // END IF: requested slot window could not be mapped

        ntwm->mapped_first_slice = first_slice;
        ntwm->mapped_last_slice_exclusive = last_slice_exclusive;
        ntwm->mapped_file_offset = map_offset;
        ntwm->mapped_length_bytes = map_length;
        ntwm->mapped_base = (const unsigned char *)new_mapping;
        ntwm->active_payload_base = ntwm->mapped_base + (requested_offset - map_offset);
        return NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS;
    } // END FUNCTION: numerical_time_window_manager_mmap_for_slot

    /**
     * Return a pointer to one mapped 3D-grid payload.
     *
     * @param[in] ntwm Numerical time-window manager.
     * @param slice_index Global slice index.
     * @return Pointer to the slice payload.
     *
     * @pre slice_index is inside the currently mapped window.
     */
    static inline const double *numerical_time_window_manager_grid_ptr(
        const NumericalTimeWindowManager *ntwm,
        const uint64_t slice_index) {
        const uint64_t local_slice = slice_index - ntwm->mapped_first_slice;
        return (const double *)(ntwm->active_payload_base + local_slice * ntwm->payload_stride_bytes);
    } // END FUNCTION: numerical_time_window_manager_grid_ptr

    /**
     * Return the mapped temporal interpolation stencil for one photon time.
     *
     * @param[in] ntwm Numerical time-window manager with an active mapped slot window.
     * @param photon_time Photon coordinate time.
     * @param temporal_interp_half_width Centered temporal interpolation half-width n.
     * @param[out] slice_indices Global slice indices for the temporal stencil.
     * @param[out] slice_times Physical coordinate times for the temporal stencil.
     * @param[out] slice_payloads Pointers to mapped 3D-grid payloads for the stencil.
     * @return NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS or NUMERICAL_TIME_WINDOW_MANAGER_ERROR.
     *
     * @pre The caller provides arrays of length at least `2*temporal_interp_half_width + 1`.
     */
    static inline int numerical_time_window_manager_stencil_for_time(
        const NumericalTimeWindowManager *ntwm,
        const double photon_time,
        const int temporal_interp_half_width,
        uint64_t *restrict slice_indices,
        REAL *restrict slice_times,
        const double **restrict slice_payloads) {
        const int temporal_num_points = 2 * temporal_interp_half_width + 1;
        const double scaled_time = (photon_time - ntwm->time_start) * ntwm->inv_grid_time_spacing;
        const long int center_slice = (long int)floor(scaled_time + 0.5);
        const long int first_slice = center_slice - (long int)temporal_interp_half_width;
        const long int last_slice_exclusive = first_slice + (long int)temporal_num_points;

        if (first_slice < (long int)ntwm->mapped_first_slice ||
            last_slice_exclusive > (long int)ntwm->mapped_last_slice_exclusive ||
            first_slice < 0L || last_slice_exclusive > (long int)ntwm->num_time_slices) {
            return NUMERICAL_TIME_WINDOW_MANAGER_ERROR;
        } // END IF: requested photon stencil is outside the active mapped window

        for (int s = 0; s < temporal_num_points; s++) {
            const uint64_t slice_index = (uint64_t)(first_slice + (long int)s);
            slice_indices[s] = slice_index;
            slice_times[s] =
                (REAL)(ntwm->time_start + ((double)slice_index) * ntwm->grid_time_spacing);
            slice_payloads[s] = numerical_time_window_manager_grid_ptr(ntwm, slice_index);
        } // END LOOP: for s over temporal interpolation stencil slices
        return NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS;
    } // END FUNCTION: numerical_time_window_manager_stencil_for_time

    """
    Bdefines_h.register_BHaH_defines(
        "numerical_time_window_manager",
        numerical_time_window_manager_c_code,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
