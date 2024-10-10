"""
Base class for generating griddata_free method.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot** com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

import nrpy.c_function as cfc


class base_register_CFunction_griddata_free:
    """Base class for registering the C function griddata_free() to free all memory within the griddata struct."""

    def __init__(
        self,
    ) -> None:

        self.desc = """Free all memory within the griddata struct,
    except perhaps non_y_n_gfs (e.g., after a regrid, in which non_y_n_gfs are freed first)."""
        self.cfunc_type = "void"
        self.name = "griddata_free"
        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

        # Default params
        self.params = "const commondata_struct *restrict commondata, griddata_struct *restrict griddata, const bool enable_free_non_y_n_gfs"
        self.body: str = ""

    def register(self) -> None:
        """Register CFunction based on class specifications."""
        cfc.register_CFunction(
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            body=self.body,
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
