"""
Base classes for generating CFunctions read_checkpoint and write_checkpoint.

Provides checkpointing capabilities to BHaH simulations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com        
"""

from typing import Tuple

import nrpy.params as par


class base_register_CFunction_read_checkpoint:
    def __init__(
        self,
        filename_tuple: Tuple[str, str] = (
            r"checkpoint-conv_factor%.2f.dat",
            "commondata->convergence_factor",
        ),
    ) -> None:
        """
        Base class to generate read_checkpoint CFunction for reading checkpoints.

        :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
        """
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
"""


class base_register_CFunction_write_checkpoint:
    def __init__(
        self,
        default_checkpoint_every: float = 2.0,
        filename_tuple: Tuple[str, str] = (
            "checkpoint-conv_factor%.2f.dat",
            "commondata->convergence_factor",
        ),
    ) -> None:
        """
        Base class to generate write_checkpoint CFunction for writing checkpoints.

        :param filename_tuple: A tuple containing the filename format and the variables to be inserted into the filename.
        :param default_checkpoint_every: The default checkpoint interval in physical time units.
        """
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
"""
