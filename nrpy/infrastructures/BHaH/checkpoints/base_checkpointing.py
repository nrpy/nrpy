"""
Base classes for generating CFunctions read_checkpoint and write_checkpoint.

Provides checkpointing capabilities to BHaH simulations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

from typing import Tuple

import nrpy.c_function as cfc
import nrpy.params as par


class base_register_CFunction_read_checkpoint:
    """
    Base class to generate read_checkpoint CFunction for reading checkpoints.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
    """

    def __init__(
        self,
        filename_tuple: Tuple[str, str] = (
            r"checkpoint-conv_factor%.2f.dat",
            "commondata->convergence_factor",
        ),
    ) -> None:
        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "unistd.h"]
        self.prefunc = r"""
#define FREAD(ptr, size, nmemb, stream) { const int numitems=fread((ptr), (size), (nmemb), (stream)); }
"""
        self.desc = "Read a checkpoint file"
        self.cfunc_type = "int"
        self.name = "read_checkpoint"
        self.params = (
            "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
        )
        self.body = rf"""
  char filename[256];
  snprintf(filename, 256, "{filename_tuple[0]}", {filename_tuple[1]});
  // If the checkpoint doesn't exist then return 0.
  if (access(filename, F_OK) != 0)
    return 0;

  FILE *cp_file = fopen(filename, "r");
  FREAD(commondata, sizeof(commondata_struct), 1, cp_file);
  fprintf(stderr, "cd struct size = %ld time=%e\n", sizeof(commondata_struct), commondata->time);
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {{
"""
        self.loop_body = r"""
    FREAD(&griddata[grid].params, sizeof(params_struct), 1, cp_file);

    int count;
    FREAD(&count, sizeof(int), 1, cp_file);

    int *restrict out_data_indices = (int *restrict)malloc(sizeof(int) * count);
    REAL *restrict compact_out_data = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * count);

    const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
    const int ntot = griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2;
    fprintf(stderr, "Reading checkpoint: grid = %d | pts = %d / %d | %d\n", grid, count, ntot, Nxx_plus_2NGHOSTS2);
    FREAD(out_data_indices, sizeof(int), count, cp_file);
    FREAD(compact_out_data, sizeof(REAL), count * NUM_EVOL_GFS, cp_file);

    MoL_malloc_y_n_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
    int which_el = 0;
#pragma omp parallel for
    for (int i = 0; i < count; i++) {
      for (int gf = 0; gf < NUM_EVOL_GFS; gf++) {
        griddata[grid].gridfuncs.y_n_gfs[IDX4pt(gf, out_data_indices[i])] = compact_out_data[i * NUM_EVOL_GFS + gf];
      }
    }
    free(out_data_indices);
    free(compact_out_data);
"""
        self.post_loop_body = r"""
  fclose(cp_file);
  fprintf(stderr, "FINISHED WITH READING\n");

  // Next set t_0 and n_0
  commondata->t_0 = commondata->time;
  commondata->nn_0 = commondata->nn;
"""

    def register(self) -> None:
        """Register function."""
        self.body += self.loop_body
        self.body += rf"""  }}
  {self.post_loop_body}
  return 1;
"""
        cfc.register_CFunction(
            prefunc=self.prefunc,
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=False,
            body=self.body,
        )


class base_register_CFunction_write_checkpoint:
    """
    Base class to generate write_checkpoint CFunction for writing checkpoints.

    :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
    :param default_checkpoint_every: The default checkpoint interval in physical time units.
    """

    def __init__(
        self,
        default_checkpoint_every: float = 2.0,
        filename_tuple: Tuple[str, str] = (
            "checkpoint-conv_factor%.2f.dat",
            "commondata->convergence_factor",
        ),
    ) -> None:
        par.register_CodeParameter(
            "REAL",
            __name__,
            "checkpoint_every",
            default_checkpoint_every,
            commondata=True,
        )
        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
        self.desc = "Write a checkpoint file"
        self.cfunc_type = "void"
        self.name = "write_checkpoint"
        self.params = "const commondata_struct *restrict commondata, griddata_struct *restrict griddata"
        self.body = rf"""
  char filename[256];
  snprintf(filename, 256, "{filename_tuple[0]}", {filename_tuple[1]});
  const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->checkpoint_every;
// Explanation of the if() below:
// Step 1: round(currtime / outevery) rounds to the nearest integer multiple of currtime.
// Step 2: Multiplying by outevery yields the exact time we should output again, t_out.
// Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible!
if (fabs(round(currtime / outevery) * outevery - currtime) < 0.5 * currdt) {{
  FILE *cp_file = fopen(filename, "w+");
  fwrite(commondata, sizeof(commondata_struct), 1, cp_file);
  fprintf(stderr, "WRITING CHECKPOINT: cd struct size = %ld time=%e\n", sizeof(commondata_struct), commondata->time);
  for(int grid=0; grid<commondata->NUMGRIDS; grid++) {{"""
        self.loop_body = ""

    def register(self) -> None:
        """Register function."""
        self.body += self.loop_body
        self.body += r"""  }
  fclose(cp_file);
  fprintf(stderr, "FINISHED WRITING CHECKPOINT\n");
  }
"""
        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=False,
            body=self.body,
        )
