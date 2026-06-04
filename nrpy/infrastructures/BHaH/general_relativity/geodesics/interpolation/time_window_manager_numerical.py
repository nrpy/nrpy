# nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/time_window_manager_numerical.py
"""
Define C helpers for numerical-spacetime time-window mmap management.

The NumericalTimeWindowManager owns the active mapped group of adjacent 3D-grid
payloads in one trusted combined numerical-spacetime binary container. It
depends on shared TimeSlotManager definitions for slot time bounds, but keeps
mmap and numerical file state out of analytic builds.

The key design constraint is that numerical-spacetime data is expensive to
load, while photons are already grouped by coordinate time through the
TimeSlotManager. This manager therefore chooses one contiguous range of
numerical 3D-grid time slices for an entire photon slot instead of loading
data ray-by-ray. The mapped range for one slot is intentionally conservative:
reverse ray tracing receives additional lower-time lookahead below the slot
boundary, and then both ends are widened by the centered temporal
interpolation halo so RKF45 steps can remain inside the active mapping while
photons cross slot boundaries.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.params as par
from nrpy.infrastructures.BHaH import BHaH_defines_h as Bdefines_h


def time_window_manager_numerical() -> None:
    """
    Register NumericalTimeWindowManager helpers in BHaH_defines.h.

    This registration owns the `rkf45_max_delta_t` CodeParameter consumed by
    the generated time-window C helpers. It also ensures the shared
    TimeSlotManager helper is registered first so the emitted code sees
    `TimeSlotManager`, `slot_lower_time()`, and `slot_upper_time()`, and so
    the commondata slot-lattice parameters already exist.

    The generated slot-window selector intentionally maps a conservative slice
    interval rather than the exact minimal interpolation stencil. The upper
    side only needs the centered temporal-interpolation halo, while the lower
    side also needs the additional backward-time RKF45 lookahead promised by
    `rkf45_max_delta_t`. The matching RKF45 step-size cap is enforced in the
    companion finalization kernel.

    Doctests:
    >>> time_window_manager_numerical()
    >>> generated = par.glb_extras_dict["BHaH_defines"]["time_window_manager_numerical"]
    >>> "time_window_manager_numerical_validate_uniform_cadence" in generated
    True
    >>> "time_window_manager_numerical_checked_floor_to_long" in generated
    True
    """
    from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.time_slot_manager_helpers import (  # pylint: disable=import-outside-toplevel
        time_slot_manager_helpers,
    )

    if "time_slot_manager" not in par.glb_extras_dict.get("BHaH_defines", {}):
        time_slot_manager_helpers()

    _ = par.register_CodeParameter(
        "REAL",
        __name__,
        "rkf45_max_delta_t",
        10.0,
        commondata=True,
        add_to_parfile=True,
        description=(
            "Maximum backward coordinate-time lookahead used by the numerical "
            "time-window manager to keep reverse ray tracing inside the mapped "
            "time-slice window."
        ),
    )

    time_window_manager_numerical_c_code = r"""
    //==========================================
    // NUMERICAL TIME-WINDOW MANAGER
    //==========================================
    // Lightweight trusted-input mmap window over one combined numerical-spacetime container.
    // Time slices are assumed ordered by increasing coordinate time, while
    // reverse ray tracing moves from larger t toward smaller t. The manager
    // maps one contiguous slot-derived time window so many photons can reuse
    // the same numerical 3D-grid payloads. The mapped window is wider than
    // the slot itself because it includes temporal-interpolation halo slices
    // and lower-time RKF45 lookahead.
#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

#define TIME_WINDOW_MANAGER_NUMERICAL_FIXED_HEADER_BYTES 4096ULL
#define TIME_WINDOW_MANAGER_NUMERICAL_SLICE_TABLE_ENTRY_BYTES 128ULL

#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_SLICE_TABLE_OFFSET 64U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NUM_TIME_SLICES 80U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_FIRST_PAYLOAD_OFFSET 104U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_PAYLOAD_BYTES_PER_SLICE 112U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_PAYLOAD_STRIDE_BYTES 120U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_POINT_RECORD_COUNT 144U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_POINT_RECORD_BYTES 160U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NUM_GRIDS 188U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NGHOSTS 192U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NXX 220U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NXX_PLUS_2NGHOSTS 244U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_DXX 340U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_INVDXX 364U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_XXMIN 388U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_XXMAX 412U
#define TIME_WINDOW_MANAGER_NUMERICAL_HEADER_CART_ORIGIN 436U
#define TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_TIME 8U

#define TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS 0
#define TIME_WINDOW_MANAGER_NUMERICAL_ERROR 1
#define TIME_WINDOW_MANAGER_NUMERICAL_MAX_TEMPORAL_INTERP_HALF_WIDTH 32

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
      uint32_t num_grids; // Writer-side grid count stored in the combined header.
      uint32_t nghosts; // Writer-side NGHOSTS stored in the combined header.
      uint64_t Nxx[3]; // Interior point counts stored in the combined header.
      uint64_t Nxx_plus_2NGHOSTS[3]; // Full logical-grid counts stored in the combined header.
      double dxx[3]; // Native coordinate spacings stored in the combined header.
      double invdxx[3]; // Reciprocal native coordinate spacings stored in the combined header.
      double xxmin[3]; // Native lower bounds stored in the combined header.
      double xxmax[3]; // Native upper bounds stored in the combined header.
      double cart_origin[3]; // Cartesian grid origin stored in the combined header.

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
    static inline uint64_t time_window_manager_numerical_load_u64(const unsigned char *bytes) {
      uint64_t value;
      memcpy(&value, bytes, sizeof(uint64_t));
      return value;
    } // END FUNCTION: time_window_manager_numerical_load_u64

    static inline uint32_t time_window_manager_numerical_load_u32(const unsigned char *bytes) {
      uint32_t value;
      memcpy(&value, bytes, sizeof(uint32_t));
      return value;
    } // END FUNCTION: time_window_manager_numerical_load_u32

    static inline double time_window_manager_numerical_load_f64(const unsigned char *bytes) {
      double value;
      memcpy(&value, bytes, sizeof(double));
      return value;
    } // END FUNCTION: time_window_manager_numerical_load_f64

    /**
     * Read exactly nbytes from one absolute file offset.
     *
     * @param fd Open file descriptor.
     * @param[out] buffer Destination buffer.
     * @param nbytes Number of bytes to read.
     * @param offset Absolute file offset.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_read_exact_at(
        int fd, void *buffer, size_t nbytes, uint64_t offset) {
      unsigned char *dst = (unsigned char *)buffer;
      size_t bytes_read = 0U;

      while (bytes_read < nbytes) {
        const ssize_t nread = pread(
            fd, dst + bytes_read, nbytes - bytes_read,
            (off_t)(offset + bytes_read));
        if (nread <= 0)
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        bytes_read += (size_t)nread;
      } // END WHILE: read requested metadata bytes
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_read_exact_at

    /**
     * Validate that stored slice times follow one uniform cadence.
     *
     * @param fd Open file descriptor.
     * @param slice_table_offset Absolute byte offset of the slice table.
     * @param num_time_slices Number of stored slice-table entries.
     * @param first_time Physical time stored for slice 0.
     * @param grid_time_spacing Expected uniform spacing between adjacent slices.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_validate_uniform_cadence(
        int fd,
        uint64_t slice_table_offset,
        uint64_t num_time_slices,
        double first_time,
        double grid_time_spacing) {
      unsigned char time_bytes[sizeof(double)];
      double previous_time = first_time;

      for (uint64_t slice_index = 2ULL; slice_index < num_time_slices; slice_index++) {
        const uint64_t time_offset =
            slice_table_offset +
            slice_index * TIME_WINDOW_MANAGER_NUMERICAL_SLICE_TABLE_ENTRY_BYTES +
            TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_TIME;
        if (time_window_manager_numerical_read_exact_at(
                fd, time_bytes, sizeof(double), time_offset) !=
            TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        } // END IF: slice time could not be read during uniform-cadence validation

        const double this_time = time_window_manager_numerical_load_f64(time_bytes);
        const double expected_time =
            first_time + (double)slice_index * grid_time_spacing;
        const double time_tolerance =
            1.0e-12 * fmax(1.0, fabs(expected_time));
        if (!isfinite(this_time) || !isfinite(expected_time) ||
            this_time <= previous_time ||
            fabs(this_time - expected_time) > time_tolerance) {
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        } // END IF: stored slice times were not uniformly spaced
        previous_time = this_time;
      } // END LOOP: for slice_index over slice-table times during uniform-cadence validation
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_validate_uniform_cadence

    /**
     * Convert one finite double to long int using floor().
     *
     * @param value Finite source value.
     * @param[out] out Converted long int destination.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_checked_floor_to_long(
        double value,
        long int *out) {
      if (out == NULL || !isfinite(value) ||
          value < (double)LONG_MIN || value > (double)LONG_MAX) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: floor() source value was not representable as long int
      *out = (long int)floor(value);
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_checked_floor_to_long

    /**
     * Convert one finite double to long int using ceil().
     *
     * @param value Finite source value.
     * @param[out] out Converted long int destination.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_checked_ceil_to_long(
        double value,
        long int *out) {
      if (out == NULL || !isfinite(value) ||
          value < (double)LONG_MIN || value > (double)LONG_MAX) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: ceil() source value was not representable as long int
      *out = (long int)ceil(value);
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_checked_ceil_to_long

    //==========================================
    // NUMERICAL WINDOW LIFETIME MANAGEMENT
    //==========================================
    /**
     * Reset a numerical time-window manager to an inert state.
     *
     * @param[in,out] ntwm Numerical time-window manager to reset.
     */
    static inline void time_window_manager_numerical_set_inert(NumericalTimeWindowManager *ntwm) {
      ntwm->fd = -1;
      ntwm->page_size = 0;
      ntwm->num_time_slices = 0ULL;
      ntwm->slice_table_offset = 0ULL;
      ntwm->first_payload_offset = 0ULL;
      ntwm->payload_bytes_per_slice = 0ULL;
      ntwm->payload_stride_bytes = 0ULL;
      ntwm->point_record_count = 0ULL;
      ntwm->point_record_bytes = 0U;
      ntwm->num_grids = 0U;
      ntwm->nghosts = 0U;
      for (int dirn = 0; dirn < 3; dirn++) {
        ntwm->Nxx[dirn] = 0ULL;
        ntwm->Nxx_plus_2NGHOSTS[dirn] = 0ULL;
        ntwm->dxx[dirn] = 0.0;
        ntwm->invdxx[dirn] = 0.0;
        ntwm->xxmin[dirn] = 0.0;
        ntwm->xxmax[dirn] = 0.0;
        ntwm->cart_origin[dirn] = 0.0;
      } // END LOOP: for dirn over stored grid metadata arrays during inert reset
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
    } // END FUNCTION: time_window_manager_numerical_set_inert

    /**
     * Close the active mmap window.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     */
    static inline void time_window_manager_numerical_close_mmap(NumericalTimeWindowManager *ntwm) {
      if (ntwm->mapped_base != NULL)
        munmap((void *)ntwm->mapped_base, ntwm->mapped_length_bytes);
      ntwm->mapped_first_slice = 0ULL;
      ntwm->mapped_last_slice_exclusive = 0ULL;
      ntwm->mapped_file_offset = 0ULL;
      ntwm->mapped_length_bytes = 0U;
      ntwm->mapped_base = NULL;
      ntwm->active_payload_base = NULL;
    } // END FUNCTION: time_window_manager_numerical_close_mmap

    /**
     * Release all resources owned by the numerical time-window manager.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     */
    static inline void time_window_manager_numerical_free(NumericalTimeWindowManager *ntwm) {
      time_window_manager_numerical_close_mmap(ntwm);
      if (ntwm->fd >= 0)
        close(ntwm->fd);
      time_window_manager_numerical_set_inert(ntwm);
    } // END FUNCTION: time_window_manager_numerical_free

    /**
     * Copy combined-file grid metadata into a runtime params_struct.
     *
     * Only the grid fields serialized by the numerical-data writer are
     * overwritten here. Higher-level caller-owned params metadata unrelated to
     * the spatial lookup contract remains untouched.
     *
     * @param[in] ntwm Initialized numerical time-window manager.
     * @param[out] params Destination params_struct to populate.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     */
    static inline int time_window_manager_numerical_apply_metadata_to_params(
        const NumericalTimeWindowManager *ntwm,
        params_struct *restrict params) {
      if (ntwm == NULL || params == NULL)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      if (ntwm->num_grids != 1U || ntwm->nghosts != (uint32_t)NGHOSTS)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      for (int dirn = 0; dirn < 3; dirn++) {
        if (ntwm->Nxx[dirn] > (uint64_t)INT_MAX ||
            ntwm->Nxx_plus_2NGHOSTS[dirn] > (uint64_t)INT_MAX)
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END LOOP: for dirn over integer-valued grid metadata before params assignment

      params->Nxx0 = (int)ntwm->Nxx[0];
      params->Nxx1 = (int)ntwm->Nxx[1];
      params->Nxx2 = (int)ntwm->Nxx[2];
      params->Nxx_plus_2NGHOSTS0 = (int)ntwm->Nxx_plus_2NGHOSTS[0];
      params->Nxx_plus_2NGHOSTS1 = (int)ntwm->Nxx_plus_2NGHOSTS[1];
      params->Nxx_plus_2NGHOSTS2 = (int)ntwm->Nxx_plus_2NGHOSTS[2];
      params->xxmin0 = (REAL)ntwm->xxmin[0];
      params->xxmin1 = (REAL)ntwm->xxmin[1];
      params->xxmin2 = (REAL)ntwm->xxmin[2];
      params->xxmax0 = (REAL)ntwm->xxmax[0];
      params->xxmax1 = (REAL)ntwm->xxmax[1];
      params->xxmax2 = (REAL)ntwm->xxmax[2];
      params->dxx0 = (REAL)ntwm->dxx[0];
      params->dxx1 = (REAL)ntwm->dxx[1];
      params->dxx2 = (REAL)ntwm->dxx[2];
      params->invdxx0 = (REAL)ntwm->invdxx[0];
      params->invdxx1 = (REAL)ntwm->invdxx[1];
      params->invdxx2 = (REAL)ntwm->invdxx[2];
      params->Cart_originx = (REAL)ntwm->cart_origin[0];
      params->Cart_originy = (REAL)ntwm->cart_origin[1];
      params->Cart_originz = (REAL)ntwm->cart_origin[2];
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_apply_metadata_to_params

    /**
     * Initialize a trusted-input numerical time-window manager.
     *
     * The caller must call time_window_manager_numerical_set_inert() before
     * the first initialization. Reinitialization releases the old mapping and
     * file before opening the replacement container.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     * @param[in] combined_path Path to the combined numerical .bin container.
     * @param[in] commondata Common runtime parameters.
     * @param temporal_interp_half_width Centered temporal interpolation half-width n.
     * @param[out] params Optional params_struct updated from the combined-file grid metadata.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     *
     * @note The backward-lookahead distance is read from
     * commondata->rkf45_max_delta_t.
     */
    static inline int time_window_manager_numerical_init(
        NumericalTimeWindowManager *ntwm,
        const char *combined_path,
        const commondata_struct *restrict commondata,
        const int temporal_interp_half_width,
        params_struct *restrict params) {
      unsigned char header_bytes[TIME_WINDOW_MANAGER_NUMERICAL_FIXED_HEADER_BYTES];
      unsigned char time_bytes[sizeof(double)];
      const uint64_t expected_point_record_bytes = 53ULL * (uint64_t)sizeof(double);

      if (ntwm == NULL)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      time_window_manager_numerical_free(ntwm);
      if (combined_path == NULL || commondata == NULL)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      if (temporal_interp_half_width < 0 ||
          temporal_interp_half_width >
              TIME_WINDOW_MANAGER_NUMERICAL_MAX_TEMPORAL_INTERP_HALF_WIDTH)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      const double max_backward_dt_lookahead = commondata->rkf45_max_delta_t;
      if (!isfinite(max_backward_dt_lookahead) || max_backward_dt_lookahead < 0.0)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      const long int page_size = sysconf(_SC_PAGESIZE);
      if (page_size <= 0)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      ntwm->fd = open(combined_path, O_RDONLY);
      if (ntwm->fd < 0)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      ntwm->temporal_interp_half_width = temporal_interp_half_width;
      ntwm->temporal_interp_num_points = 2 * temporal_interp_half_width + 1;
      ntwm->max_backward_dt_lookahead = max_backward_dt_lookahead;
      ntwm->page_size = page_size;

      if (time_window_manager_numerical_read_exact_at(
              ntwm->fd, header_bytes,
              (size_t)TIME_WINDOW_MANAGER_NUMERICAL_FIXED_HEADER_BYTES, 0ULL) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: fixed header could not be read

      ntwm->slice_table_offset = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_SLICE_TABLE_OFFSET);
      ntwm->num_time_slices = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NUM_TIME_SLICES);
      ntwm->first_payload_offset = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_FIRST_PAYLOAD_OFFSET);
      ntwm->payload_bytes_per_slice = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_PAYLOAD_BYTES_PER_SLICE);
      ntwm->payload_stride_bytes = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_PAYLOAD_STRIDE_BYTES);
      ntwm->point_record_count = time_window_manager_numerical_load_u64(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_POINT_RECORD_COUNT);
      ntwm->point_record_bytes = time_window_manager_numerical_load_u32(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_POINT_RECORD_BYTES);
      ntwm->num_grids = time_window_manager_numerical_load_u32(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NUM_GRIDS);
      ntwm->nghosts = time_window_manager_numerical_load_u32(
          header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NGHOSTS);
      for (int dirn = 0; dirn < 3; dirn++) {
        const uint64_t offset_u64 = (uint64_t)dirn * (uint64_t)sizeof(uint64_t);
        const uint64_t offset_f64 = (uint64_t)dirn * (uint64_t)sizeof(double);
        ntwm->Nxx[dirn] = time_window_manager_numerical_load_u64(
            header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NXX + offset_u64);
        ntwm->Nxx_plus_2NGHOSTS[dirn] = time_window_manager_numerical_load_u64(
            header_bytes +
            TIME_WINDOW_MANAGER_NUMERICAL_HEADER_NXX_PLUS_2NGHOSTS + offset_u64);
        ntwm->dxx[dirn] = time_window_manager_numerical_load_f64(
            header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_DXX + offset_f64);
        ntwm->invdxx[dirn] = time_window_manager_numerical_load_f64(
            header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_INVDXX + offset_f64);
        ntwm->xxmin[dirn] = time_window_manager_numerical_load_f64(
            header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_XXMIN + offset_f64);
        ntwm->xxmax[dirn] = time_window_manager_numerical_load_f64(
            header_bytes + TIME_WINDOW_MANAGER_NUMERICAL_HEADER_XXMAX + offset_f64);
        ntwm->cart_origin[dirn] = time_window_manager_numerical_load_f64(
            header_bytes +
            TIME_WINDOW_MANAGER_NUMERICAL_HEADER_CART_ORIGIN + offset_f64);
      } // END LOOP: for dirn over fixed-header grid metadata arrays during initialization

      if (ntwm->num_time_slices < 2ULL || ntwm->payload_bytes_per_slice == 0ULL ||
          ntwm->payload_stride_bytes < ntwm->payload_bytes_per_slice ||
          ntwm->num_grids != 1U || ntwm->nghosts != (uint32_t)NGHOSTS ||
          ntwm->point_record_count == 0ULL ||
          (uint64_t)ntwm->point_record_bytes != expected_point_record_bytes ||
          ntwm->point_record_count > UINT64_MAX / expected_point_record_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted file payload metadata could not describe interpolation records

      uint64_t expected_point_record_count = 1ULL;
      for (int dirn = 0; dirn < 3; dirn++) {
        if (ntwm->Nxx[dirn] == 0ULL ||
            ntwm->Nxx[dirn] > UINT64_MAX / expected_point_record_count ||
            ntwm->Nxx_plus_2NGHOSTS[dirn] !=
                ntwm->Nxx[dirn] + 2ULL * (uint64_t)ntwm->nghosts ||
            !isfinite(ntwm->dxx[dirn]) || ntwm->dxx[dirn] <= 0.0 ||
            !isfinite(ntwm->invdxx[dirn]) || ntwm->invdxx[dirn] <= 0.0 ||
            !isfinite(ntwm->xxmin[dirn]) || !isfinite(ntwm->xxmax[dirn]) ||
            ntwm->xxmax[dirn] <= ntwm->xxmin[dirn] ||
            !isfinite(ntwm->cart_origin[dirn])) {
          time_window_manager_numerical_free(ntwm);
          return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
        } // END IF: combined-file grid metadata was not numerically usable
        expected_point_record_count *= ntwm->Nxx[dirn];
      } // END LOOP: for dirn over combined-file grid metadata during validation
      if (ntwm->point_record_count != expected_point_record_count) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: point-record count disagreed with combined-file grid dimensions

      const uint64_t expected_payload_bytes =
          ntwm->point_record_count * expected_point_record_bytes;
      if (ntwm->payload_bytes_per_slice != expected_payload_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted payload length disagreed with point-record metadata

      const off_t file_size = lseek(ntwm->fd, 0, SEEK_END);
      if (file_size < (off_t)TIME_WINDOW_MANAGER_NUMERICAL_FIXED_HEADER_BYTES) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: combined numerical container size was not usable
      const uint64_t file_size_bytes = (uint64_t)file_size;

      if (ntwm->num_time_slices >
          UINT64_MAX / TIME_WINDOW_MANAGER_NUMERICAL_SLICE_TABLE_ENTRY_BYTES) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted slice-table byte count would overflow
      const uint64_t slice_table_bytes =
          ntwm->num_time_slices *
          TIME_WINDOW_MANAGER_NUMERICAL_SLICE_TABLE_ENTRY_BYTES;
      if (ntwm->slice_table_offset > UINT64_MAX - slice_table_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted slice-table end offset would overflow
      const uint64_t slice_table_end =
          ntwm->slice_table_offset + slice_table_bytes;
      if (slice_table_end > file_size_bytes ||
          slice_table_end > ntwm->first_payload_offset) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted slice table was outside the usable file layout

      if (ntwm->num_time_slices > UINT64_MAX / ntwm->payload_stride_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted payload-region byte count would overflow
      const uint64_t payload_region_bytes =
          ntwm->num_time_slices * ntwm->payload_stride_bytes;
      if (ntwm->first_payload_offset > UINT64_MAX - payload_region_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted payload-region end offset would overflow
      const uint64_t payload_region_end =
          ntwm->first_payload_offset + payload_region_bytes;
      if (payload_region_end > file_size_bytes) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted payload region exceeded the combined container

      const uint64_t first_time_offset =
          ntwm->slice_table_offset +
          TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_TIME;
      const uint64_t second_time_offset =
          ntwm->slice_table_offset +
          TIME_WINDOW_MANAGER_NUMERICAL_SLICE_TABLE_ENTRY_BYTES +
          TIME_WINDOW_MANAGER_NUMERICAL_SLICE_ENTRY_TIME;
      if (time_window_manager_numerical_read_exact_at(
              ntwm->fd, time_bytes, sizeof(double), first_time_offset) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: first slice time could not be read
      const double first_time = time_window_manager_numerical_load_f64(time_bytes);

      if (time_window_manager_numerical_read_exact_at(
              ntwm->fd, time_bytes, sizeof(double), second_time_offset) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: second slice time could not be read
      const double second_time = time_window_manager_numerical_load_f64(time_bytes);
      const double grid_time_spacing = second_time - first_time;
      if (!isfinite(first_time) || !isfinite(second_time) ||
          !isfinite(grid_time_spacing) || grid_time_spacing <= 0.0) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted file time samples were not finite and strictly increasing

      const double inv_grid_time_spacing = 1.0 / grid_time_spacing;
      if (!isfinite(inv_grid_time_spacing)) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: trusted file time spacing produced an unusable reciprocal
      if (time_window_manager_numerical_validate_uniform_cadence(
              ntwm->fd, ntwm->slice_table_offset, ntwm->num_time_slices,
              first_time, grid_time_spacing) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: stored slice times were not uniformly spaced

      ntwm->time_start = first_time;
      ntwm->grid_time_spacing = grid_time_spacing;
      ntwm->inv_grid_time_spacing = inv_grid_time_spacing;
      if (params != NULL &&
          time_window_manager_numerical_apply_metadata_to_params(ntwm, params) !=
              TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        time_window_manager_numerical_free(ntwm);
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: caller requested params_struct population from combined metadata
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_init

    //==========================================
    // NUMERICAL WINDOW SLOT RANGE SELECTION
    //==========================================
    /**
     * Compute the conservative mapped time-slice range for one TimeSlotManager slot.
     *
     * @param[in] ntwm Initialized numerical time-window manager.
     * @param[in] tsm Base time-slot manager.
     * @param slot_idx Slot index.
     * @param[out] first_slice First mapped slice, inclusive.
     * @param[out] last_slice_exclusive One past the final mapped slice.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     *
     * @note This mapped interval is intentionally wider than the exact slice
     * set requested later by time_window_manager_numerical_stencil_for_time().
     * The upper side covers the centered temporal-interpolation halo above the
     * slot, while the lower side covers both the centered halo and the extra
     * backward-time RKF45 lookahead stored in
     * ntwm->max_backward_dt_lookahead. This asymmetry matches backward
     * raytracing, where photons move toward lower coordinate time.
     */
    static inline int time_window_manager_numerical_required_grid_range(
        const NumericalTimeWindowManager *ntwm,
        const TimeSlotManager *tsm,
        const int slot_idx,
        uint64_t *first_slice,
        uint64_t *last_slice_exclusive) {
      const double t_lower = slot_lower_time(tsm, slot_idx);
      const double t_upper = slot_upper_time(tsm, slot_idx);
      // The mapped query is biased toward lower coordinate time because the
      // photons are evolved backward in time.
      const double query_min = t_lower - ntwm->max_backward_dt_lookahead;
      const double query_max = t_upper;
      const double scaled_min =
          (query_min - ntwm->time_start) * ntwm->inv_grid_time_spacing;
      const double scaled_max =
          (query_max - ntwm->time_start) * ntwm->inv_grid_time_spacing;
      long int center_first = 0L;
      long int center_last = 0L;
      if (time_window_manager_numerical_checked_floor_to_long(
              scaled_min, &center_first) !=
              TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS ||
          time_window_manager_numerical_checked_ceil_to_long(
              scaled_max, &center_last) !=
              TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        *first_slice = 0ULL;
        *last_slice_exclusive = 0ULL;
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: slot query bounds could not be represented as long int
      // Conceptually, for temporal half-width n, centered interpolation uses
      // 2n + 1 slices. A photon anywhere in the slot may therefore need
      // about n + 1 halo slices outside the queried coordinate-time
      // interval. Since photons evolve backward in time, the queried
      // interval is extended below the slot before the halo is applied:
      //   query_min = slot_lower - max_backward_dt_lookahead
      //   query_max = slot_upper
      //   mapped range = [floor(query_min / dt_grid) - (n + 1),
      //                   ceil(query_max / dt_grid) + (n + 1))
      // The exclusive upper bound in the half-open slice indexing adds the
      // final +1 seen below.
      const long int halo = (long int)ntwm->temporal_interp_half_width + 1L;
      const long int required_first = center_first - halo;
      const long int required_last_exclusive = center_last + halo + 1L;

      if (required_first < 0L ||
          required_last_exclusive > (long int)ntwm->num_time_slices ||
          required_last_exclusive <= required_first) {
        *first_slice = 0ULL;
        *last_slice_exclusive = 0ULL;
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: requested conservative slot mapping window is outside the trusted file

      *first_slice = (uint64_t)required_first;
      *last_slice_exclusive = (uint64_t)required_last_exclusive;
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_required_grid_range

    /**
     * Map exactly one slot's required contiguous time-slice payload window.
     *
     * @param[in,out] ntwm Numerical time-window manager.
     * @param[in] tsm Base time-slot manager.
     * @param slot_idx Slot index.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     *
     * @pre No interpolation reader may hold pointers returned by
     * time_window_manager_numerical_grid_ptr during this call.
     */
    static inline int time_window_manager_numerical_mmap_for_slot(
        NumericalTimeWindowManager *ntwm,
        const TimeSlotManager *tsm,
        const int slot_idx) {
      if (ntwm == NULL || tsm == NULL || ntwm->fd < 0 || ntwm->page_size <= 0 ||
          ntwm->num_time_slices < 2ULL || ntwm->payload_bytes_per_slice == 0ULL ||
          ntwm->payload_stride_bytes < ntwm->payload_bytes_per_slice ||
          !isfinite(ntwm->time_start) || !isfinite(ntwm->grid_time_spacing) ||
          ntwm->grid_time_spacing <= 0.0 || !isfinite(ntwm->inv_grid_time_spacing))
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      uint64_t first_slice = 0ULL;
      uint64_t last_slice_exclusive = 0ULL;
      if (time_window_manager_numerical_required_grid_range(
              ntwm, tsm, slot_idx, &first_slice, &last_slice_exclusive) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      if (ntwm->mapped_base != NULL && ntwm->mapped_first_slice == first_slice &&
          ntwm->mapped_last_slice_exclusive == last_slice_exclusive)
        return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;

      if (last_slice_exclusive <= first_slice ||
          last_slice_exclusive > ntwm->num_time_slices)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      const uint64_t last_slice = last_slice_exclusive - 1ULL;
      if (first_slice > (UINT64_MAX - ntwm->first_payload_offset) /
                            ntwm->payload_stride_bytes ||
          last_slice > (UINT64_MAX - ntwm->first_payload_offset) /
                           ntwm->payload_stride_bytes)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      const uint64_t requested_offset =
          ntwm->first_payload_offset + first_slice * ntwm->payload_stride_bytes;
      const uint64_t last_payload_offset =
          ntwm->first_payload_offset + last_slice * ntwm->payload_stride_bytes;
      if (last_payload_offset > UINT64_MAX - ntwm->payload_bytes_per_slice)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      const uint64_t requested_end =
          last_payload_offset + ntwm->payload_bytes_per_slice;
      const uint64_t page_size = (uint64_t)ntwm->page_size;
      const uint64_t map_offset =
          requested_offset - (requested_offset % page_size);
      if (requested_end <= map_offset)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      const uint64_t map_length_u64 = requested_end - map_offset;
      if (map_length_u64 > (uint64_t)SIZE_MAX)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      const size_t map_length = (size_t)map_length_u64;

      time_window_manager_numerical_close_mmap(ntwm);
      void *new_mapping = mmap(NULL, map_length, PROT_READ, MAP_PRIVATE, ntwm->fd,
                               (off_t)map_offset);
      if (new_mapping == MAP_FAILED)
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;

      ntwm->mapped_first_slice = first_slice;
      ntwm->mapped_last_slice_exclusive = last_slice_exclusive;
      ntwm->mapped_file_offset = map_offset;
      ntwm->mapped_length_bytes = map_length;
      ntwm->mapped_base = (const unsigned char *)new_mapping;
      ntwm->active_payload_base = ntwm->mapped_base + (requested_offset - map_offset);
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_mmap_for_slot

    /**
     * Return a pointer to one mapped 3D-grid payload.
     *
     * @param[in] ntwm Numerical time-window manager.
     * @param slice_index Global slice index.
     * @return Pointer to the slice payload.
     *
     * @pre slice_index is inside the currently mapped window.
     */
    static inline const double *time_window_manager_numerical_grid_ptr(
        const NumericalTimeWindowManager *ntwm,
        const uint64_t slice_index) {
      const uint64_t local_slice = slice_index - ntwm->mapped_first_slice;
      return (const double *)(ntwm->active_payload_base +
                              local_slice * ntwm->payload_stride_bytes);
    } // END FUNCTION: time_window_manager_numerical_grid_ptr

    /**
     * Return the mapped temporal interpolation stencil for one photon time.
     *
     * @param[in] ntwm Numerical time-window manager with an active mapped slot window.
     * @param photon_time Photon coordinate time.
     * @param temporal_interp_half_width Centered temporal interpolation half-width n.
     * @param[out] slice_indices Global slice indices for the temporal stencil.
     * @param[out] slice_times Physical coordinate times for the temporal stencil.
     * @param[out] slice_payloads Pointers to mapped 3D-grid payloads for the stencil.
     * @return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS or TIME_WINDOW_MANAGER_NUMERICAL_ERROR.
     *
     * @pre The caller provides arrays of length at least `2*temporal_interp_half_width + 1`.
     *
     * @note This routine returns the exact centered stencil for one photon,
     * whereas time_window_manager_numerical_required_grid_range() maps the
     * larger slot-level window needed so every photon in the slot can request
     * its own stencil without forcing a remap.
     */
    static inline int time_window_manager_numerical_stencil_for_time(
        const NumericalTimeWindowManager *ntwm,
        const double photon_time,
        const int temporal_interp_half_width,
        uint64_t *restrict slice_indices,
        REAL *restrict slice_times,
        const double **restrict slice_payloads) {
      if (temporal_interp_half_width < 0 ||
          temporal_interp_half_width >
              TIME_WINDOW_MANAGER_NUMERICAL_MAX_TEMPORAL_INTERP_HALF_WIDTH) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: requested temporal interpolation half-width was outside the supported range
      const int temporal_num_points = 2 * temporal_interp_half_width + 1;
      if (temporal_interp_half_width != ntwm->temporal_interp_half_width ||
          temporal_num_points != ntwm->temporal_interp_num_points) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: requested temporal interpolation stencil disagreed with the active mapped window
      const double scaled_time =
          (photon_time - ntwm->time_start) * ntwm->inv_grid_time_spacing;
      long int center_slice = 0L;
      if (time_window_manager_numerical_checked_floor_to_long(
              scaled_time + 0.5, &center_slice) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: photon-centered slice index could not be represented as long int
      const long int first_slice =
          center_slice - (long int)temporal_interp_half_width;
      const long int last_slice_exclusive =
          first_slice + (long int)temporal_num_points;

      if (first_slice < (long int)ntwm->mapped_first_slice ||
          last_slice_exclusive >
              (long int)ntwm->mapped_last_slice_exclusive ||
          first_slice < 0L ||
          last_slice_exclusive > (long int)ntwm->num_time_slices) {
        return TIME_WINDOW_MANAGER_NUMERICAL_ERROR;
      } // END IF: requested photon stencil is outside the active mapped window

      for (int s = 0; s < temporal_num_points; s++) {
        const uint64_t slice_index = (uint64_t)(first_slice + (long int)s);
        slice_indices[s] = slice_index;
        slice_times[s] =
            (REAL)(ntwm->time_start + ((double)slice_index) * ntwm->grid_time_spacing);
        slice_payloads[s] =
            time_window_manager_numerical_grid_ptr(ntwm, slice_index);
      } // END LOOP: for s over temporal interpolation stencil slices
      return TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS;
    } // END FUNCTION: time_window_manager_numerical_stencil_for_time

    """
    Bdefines_h.register_BHaH_defines(
        "time_window_manager_numerical",
        time_window_manager_numerical_c_code,
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
