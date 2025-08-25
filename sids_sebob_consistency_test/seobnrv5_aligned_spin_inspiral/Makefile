CC ?= gcc  # assigns the value CC to gcc only if environment variable CC is not already set

CFLAGS = -std=gnu99 -O2 -march=native -g -Wall $(shell gsl-config --cflags)
VALGRIND_CFLAGS = -std=gnu99 -O2 -g -Wall -Wno-unknown-pragmas $(shell gsl-config --cflags)
INCLUDEDIRS =
LDFLAGS = $(shell gsl-config --libs) -lm

# Check for OpenMP support
OPENMP_FLAG = -fopenmp
COMPILER_SUPPORTS_OPENMP := $(shell echo | $(CC) $(OPENMP_FLAG) -E - >/dev/null 2>&1 && echo YES || echo NO)

ifeq ($(COMPILER_SUPPORTS_OPENMP), YES)
    CFLAGS += $(OPENMP_FLAG)
    LDFLAGS += $(OPENMP_FLAG)
endif

OBJ_FILES = BOB_aligned_spin_NQC_rhs.o BOB_aligned_spin_waveform.o BOB_aligned_spin_waveform_from_times.o cmdline_input_and_parfile_parser.o commondata_io.o commondata_struct_set_to_default.o eval_abs_deriv.o find_local_minimum_index.o handle_gsl_return_status.o main.o params_struct_set_to_default.o SEOBNRv5_aligned_multidimensional_root_wrapper.o SEOBNRv5_aligned_spin_argrelmin.o SEOBNRv5_aligned_spin_augments.o SEOBNRv5_aligned_spin_coefficients.o SEOBNRv5_aligned_spin_flux.o SEOBNRv5_aligned_spin_gamma_wrapper.o SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit.o SEOBNRv5_aligned_spin_IMR_waveform.o SEOBNRv5_aligned_spin_initial_conditions_conservative.o SEOBNRv5_aligned_spin_initial_conditions_dissipative.o SEOBNRv5_aligned_spin_interpolate_dynamics.o SEOBNRv5_aligned_spin_interpolate_modes.o SEOBNRv5_aligned_spin_iterative_refinement.o SEOBNRv5_aligned_spin_NQC_corrections.o SEOBNRv5_aligned_spin_ode_integration.o SEOBNRv5_aligned_spin_radial_momentum_conditions.o SEOBNRv5_aligned_spin_right_hand_sides.o SEOBNRv5_aligned_spin_unwrap.o SEOBNRv5_aligned_spin_waveform.o SEOBNRv5_aligned_spin_waveform_from_dynamics.o

all: seobnrv5_aligned_spin_inspiral

%.o: %.c $(COMMON_HEADERS)
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

seobnrv5_aligned_spin_inspiral: $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)


valgrind: clean
	$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" all
	valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all -s ./seobnrv5_aligned_spin_inspiral

# Use $(RM) to be cross-platform compatible.
clean:
	$(RM) *.o */*.o *~ */*~ ./#* *.txt *.gp *.dat *.avi *.png seobnrv5_aligned_spin_inspiral
