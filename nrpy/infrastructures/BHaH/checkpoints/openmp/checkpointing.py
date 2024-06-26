"""
Register CFunctions read_checkpoint and write_checkpoint using OpenMP parallelization.

Provides checkpointing capabilities to BHaH simulations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

from typing import Tuple

import nrpy.infrastructures.BHaH.checkpoints.base_checkpointing as base_chkpt


class register_CFunction_read_checkpoint(
    base_chkpt.base_register_CFunction_read_checkpoint
):
    r"""
    Register read_checkpoint CFunction for reading checkpoints.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.

    DOCTEST:
    >>> import nrpy.c_function as cfc, json
    >>> from nrpy.helpers.generic import decompress_base64_to_string, diff_strings
    >>> with open("nrpy/infrastructures/BHaH/checkpoints/tests/DOCTEST-openmp__register_CFunction_read_checkpoint.json",'r') as f:
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
    >>> with open("nrpy/infrastructures/BHaH/checkpoints/tests/DOCTEST-openmp__register_CFunction_write_checkpoint.json",'r') as f:
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

        self.loop_body = r"""fwrite(&griddata[grid].params, sizeof(params_struct), 1, cp_file);
    const int ntot = ( griddata[grid].params.Nxx_plus_2NGHOSTS0*
                       griddata[grid].params.Nxx_plus_2NGHOSTS1*
                       griddata[grid].params.Nxx_plus_2NGHOSTS2 );
    // First we free up memory so we can malloc more: copy y_n_gfs to diagnostic_output_gfs & then free y_n_gfs.
#pragma omp parallel for
    for(int i=0;i<ntot*NUM_EVOL_GFS;i++) griddata[grid].gridfuncs.diagnostic_output_gfs[i] = griddata[grid].gridfuncs.y_n_gfs[i];
    MoL_free_memory_y_n_gfs(&griddata[grid].gridfuncs);

    int count = 0;
    const int maskval = 1; // to be replaced with griddata[grid].mask[i].
#pragma omp parallel for reduction(+:count)
    for(int i=0;i<ntot;i++) {
      if(maskval >= +0) count++;
    }
    fwrite(&count, sizeof(int), 1, cp_file);

    int  *restrict out_data_indices = (int  *restrict)malloc(sizeof(int)                 * count);
    REAL *restrict compact_out_data = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);
    int which_el = 0;
    for(int i=0;i<ntot;i++) {
      if(maskval >= +0) {
        out_data_indices[which_el] = i;
        for(int gf=0; gf<NUM_EVOL_GFS; gf++)
          compact_out_data[which_el*NUM_EVOL_GFS + gf] = griddata[grid].gridfuncs.diagnostic_output_gfs[ntot*gf + i];
        which_el++;
      }
    }
    //printf("HEY which_el = %d | count = %d\n", which_el, count);
    fwrite(out_data_indices, sizeof(int) , count               , cp_file);
    fwrite(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file);
    free(out_data_indices); free(compact_out_data);
    MoL_malloc_y_n_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
#pragma omp parallel for
    for(int i=0;i<ntot*NUM_EVOL_GFS;i++) griddata[grid].gridfuncs.y_n_gfs[i] = griddata[grid].gridfuncs.diagnostic_output_gfs[i];
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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
