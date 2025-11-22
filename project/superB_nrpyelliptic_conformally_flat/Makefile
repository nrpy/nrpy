CC = charmc
CFLAGS = -O2 -march=native -g -fpermissive
#CFLAGS = -O3 -funroll-loops -march=native -g -Wall -fpermissive
#CFLAGS = -O2 -g -Wall -Wno-unknown-pragmas -fpermissive

VALGRIND_CFLAGS = -fpermissive -O0 -g -fno-tree-vectorize -mno-avx512f -mno-avx512vl -mno-avx512dq -mno-avx512bw
INCLUDEDIRS =
LDFLAGS = -module CkIO -lm

OBJ_FILES = MoL/MoL_free_intermediate_stage_gfs.o MoL/MoL_malloc_intermediate_stage_gfs.o MoL_step_forward_in_time.o MoL_sync_data_defines.o SinhSpherical/apply_bcs_outerradiation_and_inner__rfm__SinhSpherical.o SinhSpherical/auxevol_gfs_set_to_constant__rfm__SinhSpherical.o SinhSpherical/bcstruct_chare_set_up__rfm__SinhSpherical.o SinhSpherical/bcstruct_set_up__rfm__SinhSpherical.o SinhSpherical/charecommstruct_set_up__rfm__SinhSpherical.o SinhSpherical/ds_min_single_pt__rfm__SinhSpherical.o SinhSpherical/numerical_grid_params_Nxx_dxx_xx__rfm__SinhSpherical.o SinhSpherical/numerical_grid_params_Nxx_dxx_xx_chare__rfm__SinhSpherical.o SinhSpherical/rfm_precompute_defines__rfm__SinhSpherical.o SinhSpherical/rfm_precompute_free__rfm__SinhSpherical.o SinhSpherical/rfm_precompute_malloc__rfm__SinhSpherical.o SinhSpherical/rhs_eval__rfm__SinhSpherical.o SinhSpherical/xx_to_Cart__rfm__SinhSpherical.o apply_bcs_inner_only.o apply_bcs_inner_only_nonlocal.o apply_bcs_inner_only_specific_auxgfs.o apply_bcs_outerextrap_and_inner.o apply_bcs_outerradiation_and_inner.o auxevol_gfs_set_to_constant.o bcstruct_chare_set_up.o bcstruct_set_up.o cfl_limited_timestep.o charecommstruct_set_up.o cmdline_input_and_parfile_parser.o commondata_struct_set_to_default.o diagnostics/SinhSpherical/diagnostics_nearest_1d_y_and_z_axes__rfm__SinhSpherical.o diagnostics/SinhSpherical/diagnostics_nearest_2d_xy_and_yz_planes__rfm__SinhSpherical.o diagnostics/SinhSpherical/diagnostics_nearest_grid_center__rfm__SinhSpherical.o diagnostics/SinhSpherical/residual_H_compute_all_points__rfm__SinhSpherical.o diagnostics/SinhSpherical/sqrt_detgammahat_d3xx_volume_element__rfm__SinhSpherical.o diagnostics/diagnostic_gfs_set.o diagnostics/diagnostics.o diagnostics/diagnostics_nearest.o diagnostics/diagnostics_nearest_1d_y_and_z_axes.o diagnostics/diagnostics_nearest_2d_xy_and_yz_planes.o diagnostics/diagnostics_nearest_grid_center.o diagnostics/diagnostics_volume_integration.o diagnostics/residual_H_compute_all_points.o diagnostics/sqrt_detgammahat_d3xx_volume_element.o ds_min_single_pt.o griddata_free.o initial_data.o initial_guess_single_point.o initialize_yn_and_non_yn_gfs_to_nan.o log10_L2norm_gf.o numerical_grid_params_Nxx_dxx_xx.o numerical_grid_params_Nxx_dxx_xx_chare.o numerical_grids_and_timestep.o numerical_grids_chare.o params_struct_set_to_default.o progress_indicator.o rfm_precompute_defines.o rfm_precompute_free.o rfm_precompute_malloc.o rhs_eval.o stop_conditions_check.o superB/superB_pup_routines.o timestepping_free_memory_tmpBuffer.o timestepping_malloc_tmpBuffer.o xx_to_Cart.o

all: superB_nrpyelliptic_conformally_flat

%.o: %.cpp $(COMMON_HEADERS)
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

timestepping.o: timestepping.cpp timestepping.h main.h timestepping.def.h
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

timestepping.h: timestepping.decl.h
main.h: timestepping.decl.h main.decl.h

timestepping.decl.h timestepping.def.h: timestepping.ci
	$(CC) $(CFLAGS) $(INCLUDEDIRS) timestepping.ci

main.o: main.cpp main.h main.decl.h main.def.h timestepping.decl.h
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

main.decl.h main.def.h: main.ci
	$(CC) $(CFLAGS) $(INCLUDEDIRS) main.ci

superB_nrpyelliptic_conformally_flat: $(OBJ_FILES) timestepping.o main.o
	$(CC) -language charm++ $^ -o $@ $(LDFLAGS)

# Use $(RM) to be cross-platform compatible.
#~ clean:
#~ 	$(RM) *.o */*.o *~ */*~ ./#* *.txt *.gp *.dat *.avi *.png superB_nrpyelliptic_conformally_flat *.decl.h *.def.h charmrun
#~ 	$(RM) -r log
clean:
	$(RM) *.o */*.o */*/*.o *~ */*~ ./#* *.txt *.gp *.dat *.avi *.png superB_nrpyelliptic_conformally_flat *.decl.h *.def.h charmrun
	$(RM) -r log


#~ valgrind: clean
#~ 	$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" all
#~ 	valgrind --track-origins=yes --leak-check=full -s --exit-on-first-error=yes ./superB_nrpyelliptic_conformally_flat
#~ valgrind: clean
#~ 	$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" all
#~ 	valgrind --track-origins=yes --leak-check=full -s \
#~ 	         --exit-on-first-error=yes --error-exitcode=99 \
#~ 	         ./superB_nrpyelliptic_conformally_flat
valgrind: clean
	$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" all
	valgrind --track-origins=yes --leak-check=full -s \
	         ./superB_nrpyelliptic_conformally_flat

