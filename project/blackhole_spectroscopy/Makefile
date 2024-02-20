CC ?= gcc  # assigns the value CC to gcc only if environment variable CC is not already set
CFLAGS = -O2 -march=native -g -Wall -Wno-unused-variable -std=gnu99 $(shell gsl-config --cflags)
#CFLAGS = -O3 -funroll-loops -march=native -g -Wall -Wno-unused-variable -std=gnu99 -std=gnu99 $(shell gsl-config --cflags)
#CFLAGS = -O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas -std=gnu99 $(shell gsl-config --cflags)

INCLUDEDIRS =
LDFLAGS = $(shell gsl-config --libs) -lm

# Check for OpenMP support
OPENMP_FLAG = -fopenmp
COMPILER_SUPPORTS_OPENMP := $(shell echo | $(CC) $(OPENMP_FLAG) -E - >/dev/null 2>&1 && echo YES || echo NO)

ifeq ($(COMPILER_SUPPORTS_OPENMP), YES)
    CFLAGS += $(OPENMP_FLAG)
    LDFLAGS += $(OPENMP_FLAG)  # -lgomp does not work with clang in termux
endif

OBJ_FILES = Cart_to_xx_and_nearest_i0i1i2.o CoordSystem_hash_setup.o MoL_free_memory_non_y_n_gfs.o MoL_free_memory_y_n_gfs.o MoL_malloc_non_y_n_gfs.o MoL_malloc_y_n_gfs.o MoL_step_forward_in_time.o NRPyPN_quasicircular_momenta.o Ricci_eval.o SinhSpherical/Cart_to_xx_and_nearest_i0i1i2__rfm__SinhSpherical.o SinhSpherical/Ricci_eval__rfm__SinhSpherical.o SinhSpherical/apply_bcs_outerradiation_and_inner__rfm__SinhSpherical.o SinhSpherical/bcstruct_set_up__rfm__SinhSpherical.o SinhSpherical/cfl_limited_timestep__rfm__SinhSpherical.o SinhSpherical/constraints_eval__rfm__SinhSpherical.o SinhSpherical/diagnostics_nearest_1d_y_axis__rfm__SinhSpherical.o SinhSpherical/diagnostics_nearest_1d_z_axis__rfm__SinhSpherical.o SinhSpherical/diagnostics_nearest_2d_xy_plane__rfm__SinhSpherical.o SinhSpherical/diagnostics_nearest_2d_yz_plane__rfm__SinhSpherical.o SinhSpherical/diagnostics_nearest_grid_center__rfm__SinhSpherical.o SinhSpherical/enforce_detgammabar_equals_detgammahat__rfm__SinhSpherical.o SinhSpherical/initial_data_reader__convert_ADM_Cartesian_to_BSSN__rfm__SinhSpherical.o SinhSpherical/numerical_grid_params_Nxx_dxx_xx__rfm__SinhSpherical.o SinhSpherical/psi4_part0__rfm__SinhSpherical.o SinhSpherical/psi4_part1__rfm__SinhSpherical.o SinhSpherical/psi4_part2__rfm__SinhSpherical.o SinhSpherical/psi4_tetrad__rfm__SinhSpherical.o SinhSpherical/rfm_precompute_defines__rfm__SinhSpherical.o SinhSpherical/rfm_precompute_free__rfm__SinhSpherical.o SinhSpherical/rfm_precompute_malloc__rfm__SinhSpherical.o SinhSpherical/rhs_eval__rfm__SinhSpherical.o SinhSpherical/xx_to_Cart__rfm__SinhSpherical.o TwoPunctures/TP_CoordTransf.o TwoPunctures/TP_Equations.o TwoPunctures/TP_FuncAndJacobian.o TwoPunctures/TP_Interp.o TwoPunctures/TP_Newton.o TwoPunctures/TP_solve.o TwoPunctures/TP_utilities.o TwoPunctures/initialize_ID_persist_struct.o apply_bcs_inner_only.o apply_bcs_outerextrap_and_inner.o apply_bcs_outerradiation_and_inner.o bcstruct_set_up.o cfl_limited_timestep.o cmdline_input_and_parfile_parser.o commondata_struct_set_to_default.o constraints_eval.o diagnostics.o diagnostics_nearest_1d_y_axis.o diagnostics_nearest_1d_z_axis.o diagnostics_nearest_2d_xy_plane.o diagnostics_nearest_2d_yz_plane.o diagnostics_nearest_grid_center.o enforce_detgammabar_equals_detgammahat.o initial_data.o initial_data_reader__convert_ADM_Cartesian_to_BSSN.o main.o numerical_grid_params_Nxx_dxx_xx.o numerical_grids_and_timestep.o params_struct_set_to_default.o progress_indicator.o psi4_part0.o psi4_part1.o psi4_part2.o psi4_spinweightm2_decomposition_on_sphlike_grids.o psi4_tetrad.o read_checkpoint.o rfm_precompute_defines.o rfm_precompute_free.o rfm_precompute_malloc.o rhs_eval.o spin_weight_minus2_sph_harmonics.o write_checkpoint.o xx_to_Cart.o

all: blackhole_spectroscopy

%.o: %.c $(COMMON_HEADERS)
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

blackhole_spectroscopy: $(OBJ_FILES)
	$(CC) $^ -o $@ $(LDFLAGS)

# Use $(RM) to be cross-platform compatible.
clean:
	$(RM) *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png blackhole_spectroscopy
