"""
Register CFunctions read_checkpoint and write_checkpoint using OpenMP parallelization.

Provides checkpointing capabilities to BHaH simulations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

from typing import Tuple

import nrpy.infrastructures.gpu.checkpoints.base_checkpointing as base_chkpt


class register_CFunction_read_checkpoint(
    base_chkpt.base_register_CFunction_read_checkpoint
):
    """
    Register read_checkpoint CFunction for reading checkpoints.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.

    DOCTEST:
    >>> import nrpy.c_function as cfc, json
    >>> from nrpy.helpers.generic import decompress_base64_to_string, diff_strings
    >>> with open("nrpy/infrastructures/gpu/checkpoints/tests/DOCTEST-cuda__register_CFunction_read_checkpoint.json",'r') as f:
    ...     expected_str_dict = json.load(f)
    >>> key = "default"
    >>> _ = register_CFunction_read_checkpoint()
    >>> registered_str = cfc.CFunction_dict["read_checkpoint"].full_function
    >>> expected_str = decompress_base64_to_string(expected_str_dict[key])
    >>> if registered_str != expected_str:
    ...     raise ValueError(f"{key}: {diff_strings(expected_str, registered_str)}")
    """

    def __init__(
        self,
        filename_tuple: Tuple[str, str] = (
            r"checkpoint-conv_factor%.2f.dat",
            "commondata->convergence_factor",
        ),
    ) -> None:
        super().__init__(filename_tuple=filename_tuple)

        self.params += ", griddata_struct *restrict griddata_GPU"
        self.loop_body = r"""
    FREAD(&griddata[grid].params, sizeof(params_struct), 1, cp_file);

    // Copy params to griddata that is used with the device
    memcpy(&griddata_GPU[grid].params, &griddata[grid].params, sizeof(params_struct));

    int count;
    FREAD(&count, sizeof(int), 1, cp_file);

    int * out_data_indices = (int *)malloc(sizeof(int) * count);
    REAL * compact_out_data = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);

    const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
    const int ntot = griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2;
    fprintf(stderr, "Reading checkpoint: grid = %d | pts = %d / %d | %d\n", grid, count, ntot, Nxx_plus_2NGHOSTS2);
    FREAD(out_data_indices, sizeof(int), count, cp_file);
    FREAD(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file);

    // Malloc for both GFs
    CUDA__free_host_gfs(&griddata[grid].gridfuncs);
    CUDA__malloc_host_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);

    MoL_free_memory_y_n_gfs(&griddata_GPU[grid].gridfuncs);
    MoL_malloc_y_n_gfs(commondata, &griddata_GPU[grid].params, &griddata_GPU[grid].gridfuncs);
#pragma omp parallel for
    for (int i = 0; i < count; i++) {
      for (int gf = 0; gf < NUM_EVOL_GFS; gf++) {
        griddata[grid].gridfuncs.y_n_gfs[IDX4pt(gf, out_data_indices[i])] = compact_out_data[i * NUM_EVOL_GFS + gf];
      }
    }
    free(out_data_indices);
    free(compact_out_data);
"""
        self.loop_body += """

    // Set gridfunctions aliases for HOST data
    REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;

    // Set gridfunctions aliases for GPU data
    REAL *restrict y_n_gfs_GPU = griddata_GPU[grid].gridfuncs.y_n_gfs;
    for(int gf = 0; gf < NUM_EVOL_GFS; ++gf) {
      cpyHosttoDevice__gf(commondata, &griddata[grid].params, y_n_gfs, y_n_gfs_GPU, gf, gf);
    }
"""
        self.post_loop_body += """
// local stream syncs?
cudaDeviceSynchronize();
"""
        self.register()


class register_CFunction_write_checkpoint(
    base_chkpt.base_register_CFunction_write_checkpoint
):
    r"""
    Register write_checkpoint CFunction for writing checkpoints.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
    :param default_checkpoint_every: The default checkpoint interval in physical time units.

    DOCTEST:
    >>> import nrpy.c_function as cfc, json
    >>> from nrpy.helpers.generic import decompress_base64_to_string, diff_strings
    >>> with open("nrpy/infrastructures/gpu/checkpoints/tests/DOCTEST-cuda__register_CFunction_write_checkpoint.json",'r') as f:
    ...     expected_str_dict = json.load(f)
    >>> key = "default"
    >>> _ = register_CFunction_write_checkpoint()
    >>> registered_str = cfc.CFunction_dict["write_checkpoint"].full_function
    >>> expected_str = decompress_base64_to_string(expected_str_dict[key])
    >>> if registered_str != expected_str:
    ...     raise ValueError(f"\\n{key}: {diff_strings(expected_str, registered_str)}")
    """

    def __init__(
        self,
        default_checkpoint_every: float = 2.0,
        filename_tuple: Tuple[str, str] = (
            "checkpoint-conv_factor%.2f.dat",
            "commondata->convergence_factor",
        ),
    ) -> None:
        super().__init__(
            default_checkpoint_every=default_checkpoint_every,
            filename_tuple=filename_tuple,
        )
        self.params += ", griddata_struct *restrict griddata_GPU"

        self.loop_body = r"""
      // Set gridfunctions aliases for HOST data
      REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.y_n_gfs;

      // Set gridfunctions aliases for GPU data
      REAL *restrict y_n_gfs = griddata_GPU[grid].gridfuncs.y_n_gfs;

    // Make sure host griddata has correct params
    memcpy(&griddata[grid].params, &griddata_GPU[grid].params, sizeof(params_struct));
    for(int gf = 0; gf < NUM_EVOL_GFS; ++gf) {
      cpyDevicetoHost__gf(commondata, &griddata[grid].params, diagnostic_output_gfs, y_n_gfs, gf, gf);
    }

    fwrite(&griddata[grid].params, sizeof(params_struct), 1, cp_file);
    const int ntot =
        (griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2);

    // Does this need to be on GPU?  Where is MASK?
    int count = 0;
    const int maskval = 1; // to be replaced with griddata[grid].mask[i].
#pragma omp parallel for reduction(+:count)
    for(int i=0;i<ntot;i++) {
      if(maskval >= +0) count++;
    }
    fwrite(&count, sizeof(int), 1, cp_file);

    int  * out_data_indices = (int  *)malloc(sizeof(int)                 * count);
    REAL * compact_out_data = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);
    int which_el = 0;
    // Should be a local sync?
    cudaDeviceSynchronize();
    for (int i = 0; i < ntot; i++) {
      if (maskval >= +0) {
        out_data_indices[which_el] = i;
        for (int gf = 0; gf < NUM_EVOL_GFS; gf++)
          compact_out_data[which_el * NUM_EVOL_GFS + gf] = diagnostic_output_gfs[ntot * gf + i];
        which_el++;
      }
    }
    // printf("HEY which_el = %d | count = %d\n", which_el, count);
    fwrite(out_data_indices, sizeof(int), count, cp_file);
    fwrite(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file);
    free(out_data_indices);
    free(compact_out_data);
"""
        self.register()


def register_CFunctions(
    filename_tuple: Tuple[str, str] = (
        "checkpoint-conv_factor%.2f.dat",
        "commondata->convergence_factor",
    ),
    default_checkpoint_every: float = 2.0,
) -> None:
    """
    Register CFunctions for checkpointing.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
    :param default_checkpoint_every: The default checkpoint interval in physical time units.
    """
    register_CFunction_read_checkpoint(filename_tuple=filename_tuple)
    register_CFunction_write_checkpoint(
        filename_tuple=filename_tuple, default_checkpoint_every=default_checkpoint_every
    )
