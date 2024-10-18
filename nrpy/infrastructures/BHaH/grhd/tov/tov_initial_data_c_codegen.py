"""
Reads in TOV initial data from file and store in memory.

The reader expects the following ID_persist struct to be defined:

typedef struct ID_persist {
  const char *filename;
  int N_r, star_radius_idx;
  double star_radius;
  double *r_Schw_arr, *rho_arr, *rho_baryon_arr, *P_arr;
  double *M_arr, *expnu_arr, *exp4phi_arr, *rbar_arr;
} ID_persist;

Author: Leonardo Rosa Werneck
        wernecklr **at** gmail **dot* com
"""

import nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library as BCl
from nrpy.c_codegen import c_codegen
from nrpy.c_function import CFunction_dict, register_CFunction
from nrpy.equations.grhd.tov.tov_quantities import tov_ADM_quantities, tov_T4UU
from nrpy.params import set_parval_from_str


def register_CFunction_tov_reader() -> None:
    """Register tov_reader function, to be used by the initial_data function."""
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = "TOV initial data reader"
    cfunc_type = "void"
    name = "tov_reader"
    params = "ID_persist_struct *restrict ID_persist"
    include_CodeParameters_h = False
    prefunc = r"""
static int count_num_lines_in_file(const char *filename) {

  // Open the input file
  FILE *fp = fopen(filename, "r");

  // Check the file exists
  if(!fp) {
    fprintf(stderr, "Could not open file %s\n", filename);
    exit(1);
  }

  // Count the number of lines in the file
  size_t len = 0;
  int numlines_in_file = 0;
  char *line = NULL;
  while(getline(&line, &len, fp) != -1) {
    numlines_in_file++;
  }

  // Clean up
  free(line);
  fclose(fp);

  return numlines_in_file;
}

    """
    body = r"""
// Get the number of points in the file
ID_persist->N_r = count_num_lines_in_file("tov.asc");

// Allocate memory for the file data
ID_persist->r_Schw_arr     = malloc(sizeof(REAL) * ID_persist->N_r);
ID_persist->rho_arr        = malloc(sizeof(REAL) * ID_persist->N_r);
ID_persist->rho_baryon_arr = malloc(sizeof(REAL) * ID_persist->N_r);
ID_persist->P_arr          = malloc(sizeof(REAL) * ID_persist->N_r);
ID_persist->M_arr          = malloc(sizeof(REAL) * ID_persist->N_r);
ID_persist->expnu_arr      = malloc(sizeof(REAL) * ID_persist->N_r);
ID_persist->exp4phi_arr    = malloc(sizeof(REAL) * ID_persist->N_r);
ID_persist->rbar_arr       = malloc(sizeof(REAL) * ID_persist->N_r);

// Read in the file
const char *delimiters = " \t";
char *line = NULL;
size_t len = 0;
FILE *fp = fopen("tov.asc", "r");

for(int which_line = 0; which_line < ID_persist->N_r; which_line++) {
  if(getline(&line, &len, fp) == -1) {
    fprintf(stderr, "Error: failed to read line %d from file\n", which_line);
  }

  char *token=strtok(line, delimiters); if(!token) { fprintf(stderr, "BADDDD\n"); exit(1); }
  ID_persist->r_Schw_arr[which_line]     = strtod(token, NULL); token = strtok(NULL, delimiters);
  ID_persist->rho_arr[which_line]        = strtod(token, NULL); token = strtok(NULL, delimiters);
  ID_persist->rho_baryon_arr[which_line] = strtod(token, NULL); token = strtok(NULL, delimiters);
  ID_persist->P_arr[which_line]          = strtod(token, NULL); token = strtok(NULL, delimiters);
  ID_persist->M_arr[which_line]          = strtod(token, NULL); token = strtok(NULL, delimiters);
  ID_persist->expnu_arr[which_line]      = strtod(token, NULL); token = strtok(NULL, delimiters);
  ID_persist->exp4phi_arr[which_line]    = strtod(token, NULL); token = strtok(NULL, delimiters);
  ID_persist->rbar_arr[which_line]       = strtod(token, NULL);
}

// Clean up
free(line);
fclose(fp);

// Finally set Rbar and Rbar_idx
ID_persist->star_radius     = -100.0;
ID_persist->star_radius_idx = -100;
for(int i=1;i<ID_persist->N_r;i++) {
  if(ID_persist->rho_arr[i-1] > 0  &&  ID_persist->rho_arr[i] == 0) {
    ID_persist->star_radius = ID_persist->rbar_arr[i-1];
    ID_persist->star_radius_idx = i-1;
  }
}
if(ID_persist->star_radius < 0) {
  fprintf(stderr,"Error: could not find the star radius from data file.\n");
  exit(1);
}
"""

    register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        prefunc=prefunc,
        body=body,
    )


def register_CFunction_tov_free() -> None:
    """Register a function to free memory for the ID_persist struct."""
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = "TOV initial data interpolator"
    cfunc_type = "void"
    name = "tov_free"
    params = "ID_persist_struct *restrict ID_persist"
    include_CodeParameters_h = False
    body = """
free(ID_persist->r_Schw_arr);
free(ID_persist->rho_arr);
free(ID_persist->rho_baryon_arr);
free(ID_persist->P_arr);
free(ID_persist->M_arr);
free(ID_persist->expnu_arr);
free(ID_persist->exp4phi_arr);
free(ID_persist->rbar_arr);
"""

    register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )


def register_CFunction_tov_interp() -> None:
    """Register a function to perform radial interpolation of TOV initial data."""
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = "TOV initial data interpolator"
    cfunc_type = "void"
    name = "tov_interp"
    prefunc = r"""
static inline int
bisection_idx_finder(const REAL rrbar, const int numlines_in_file, const REAL *restrict rbar_arr) {
  int x1 = 0;
  int x2 = numlines_in_file-1;
  REAL y1 = rrbar-rbar_arr[x1];
  REAL y2 = rrbar-rbar_arr[x2];
  if(y1*y2 >= 0) {
    fprintf(stderr,"INTERPOLATION BRACKETING ERROR %e | %e %e\n",rrbar,y1,y2);
    exit(1);
  }
  for(int i=0;i<numlines_in_file;i++) {
    int x_midpoint = (x1+x2)/2;
    REAL y_midpoint = rrbar-rbar_arr[x_midpoint];
    if(y_midpoint*y1 < 0) {
      x2 = x_midpoint;
      y2 = y_midpoint;
    } else {
      x1 = x_midpoint;
      y1 = y_midpoint;
    }
    if( abs(x2-x1) == 1 ) {
      // If rbar_arr[x1] is closer to rrbar than rbar_arr[x2] then return x1:
      if(fabs(rrbar-rbar_arr[x1]) < fabs(rrbar-rbar_arr[x2])) return x1;
      // Otherwiser return x2:
      return x2;
    }
  }
  fprintf(stderr,"INTERPOLATION BRACKETING ERROR: DID NOT CONVERGE.\n");
  exit(1);
}

void TOV_interpolate_1D(REAL rrbar, const ID_persist_struct *restrict ID_persist,
                        REAL *restrict rho,REAL *restrict rho_baryon,REAL *restrict P,
                        REAL *restrict M,REAL *restrict expnu,REAL *restrict exp4phi) {
                        
  const int interp_stencil_size = 12;
  const int numlines_in_file = ID_persist->N_r;
  const int Rbar_idx = ID_persist->star_radius_idx;
  const REAL Rbar = ID_persist->star_radius;
  const REAL *rbar_arr       = ID_persist->rbar_arr;
  const REAL *r_Schw_arr     = ID_persist->r_Schw_arr;
  const REAL *rho_arr        = ID_persist->rho_arr;
  const REAL *rho_baryon_arr = ID_persist->rho_baryon_arr;
  const REAL *P_arr          = ID_persist->P_arr;
  const REAL *M_arr          = ID_persist->M_arr;
  const REAL *expnu_arr      = ID_persist->expnu_arr;
  const REAL *exp4phi_arr    = ID_persist->exp4phi_arr;

  // For this case, we know that for all functions, f(r) = f(-r)
  if(rrbar < 0) rrbar = -rrbar;

  // First find the central interpolation stencil index:
  int idx = bisection_idx_finder(rrbar, ID_persist->N_r, rbar_arr);

  int idxmin = MAX(0,idx-interp_stencil_size/2-1);

  // -= Do not allow the interpolation stencil to cross the star's surface =-
  // max index is when idxmin + (interp_stencil_size-1) = Rbar_idx
  //  -> idxmin at most can be Rbar_idx - interp_stencil_size + 1
  if(rrbar < Rbar) {
    idxmin = MIN(idxmin,Rbar_idx - interp_stencil_size + 1);
  } else {
    idxmin = MAX(idxmin,Rbar_idx+1);
    idxmin = MIN(idxmin,numlines_in_file - interp_stencil_size + 1);
  }
  // Now perform the Lagrange polynomial interpolation:

  // First set the interpolation coefficients:
  REAL rbar_sample[interp_stencil_size];
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    rbar_sample[i-idxmin] = rbar_arr[i];
  }
  REAL l_i_of_r[interp_stencil_size];
  for(int i=0;i<interp_stencil_size;i++) {
    REAL numer = 1.0;
    REAL denom = 1.0;
    for(int j=0;j<i;j++) {
      numer *= rrbar - rbar_sample[j];
      denom *= rbar_sample[i] - rbar_sample[j];
    }
    for(int j=i+1;j<interp_stencil_size;j++) {
      numer *= rrbar - rbar_sample[j];
      denom *= rbar_sample[i] - rbar_sample[j];
    }
    l_i_of_r[i] = numer/denom;
  }

  // Then perform the interpolation:
  *rho = 0.0;
  *rho_baryon = 0.0;
  *P = 0.0;
  *M = 0.0;
  *expnu = 0.0;
  *exp4phi = 0.0;

  REAL r_Schw = 0.0;
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    r_Schw      += l_i_of_r[i-idxmin] * r_Schw_arr[i];
    *rho        += l_i_of_r[i-idxmin] * rho_arr[i];
    *rho_baryon += l_i_of_r[i-idxmin] * rho_baryon_arr[i];
    *P          += l_i_of_r[i-idxmin] * P_arr[i];
    *M          += l_i_of_r[i-idxmin] * M_arr[i];
    *expnu      += l_i_of_r[i-idxmin] * expnu_arr[i];
    *exp4phi    += l_i_of_r[i-idxmin] * exp4phi_arr[i];
  }

  if(rrbar > Rbar) {
    *rho        = 0;
    *rho_baryon = 0;
    *P          = 0;
    *M          = M_arr[Rbar_idx+1];
    *expnu      = 1. - 2.*(*M) / r_Schw;
    *exp4phi    = pow(r_Schw / rrbar,2.0);
  }
}
"""
    params = """
const commondata_struct *restrict commondata,
const params_struct *restrict params,
const REAL xCart[3],
const ID_persist_struct *restrict ID_persist,
initial_data_struct *restrict initial_data
"""
    include_CodeParameters_h = False
    body = r"""
  int dummy[3];
  REAL xx[3];
  Cart_to_xx_and_nearest_i0i1i2(commondata, params, xCart, xx, dummy);
  const REAL rbar = xx[0]; // Note that we use coordinates such that r = rbar
  const REAL theta = xx[1];
  const REAL phi = xx[2];

  // Next set gamma_{ij} in spherical basis
  REAL rho,rho_baryon,P,M,expnu,exp4phi;
  TOV_interpolate_1D(rbar, ID_persist, &rho, &rho_baryon, &P, &M, &expnu, &exp4phi);
"""
    alpha, betaU, BU, gammaDD, KDD = tov_ADM_quantities()
    T4UU = tov_T4UU()

    list_of_output_exprs = [alpha]
    list_of_output_varnames = ["initial_data->alpha"]
    for i in range(3):
        list_of_output_exprs += [betaU[i]]
        list_of_output_varnames += ["initial_data->betaSphorCartU" + str(i)]
        list_of_output_exprs += [BU[i]]
        list_of_output_varnames += ["initial_data->BSphorCartU" + str(i)]
        for j in range(i, 3):
            list_of_output_exprs += [gammaDD[i][j]]
            list_of_output_varnames += [
                "initial_data->gammaSphorCartDD" + str(i) + str(j)
            ]
            list_of_output_exprs += [KDD[i][j]]
            list_of_output_varnames += ["initial_data->KSphorCartDD" + str(i) + str(j)]
    for mu in range(4):
        for nu in range(mu, 4):
            list_of_output_exprs += [T4UU[mu][nu]]
            list_of_output_varnames += [
                "initial_data->T4SphorCartUU" + str(mu) + str(nu)
            ]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (
        list(t)
        for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs)))
    )

    body += c_codegen(
        list_of_output_exprs,
        list_of_output_varnames,
        include_braces=False,
        verbose=False,
    )

    register_CFunction(
        includes=includes,
        desc=desc,
        prefunc=prefunc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )


def register_CFunction_initial_data(CoordSystem: str) -> None:
    """
    Register all functions needed to set TOV initial data.

    :param CoordSystem: Coordinate system used in time evolution.
    """
    if CFunction_dict.get("tov_reader", None) is None:
        register_CFunction_tov_reader()
        register_CFunction_tov_free()
        register_CFunction_tov_interp()

    set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)
    BCl.register_CFunction_initial_data(
        CoordSystem=CoordSystem,
        IDtype="tov_interp",
        IDCoordSystem="Spherical",
        ID_persist_struct_str="""
        const char *filename;
        int N_r, star_radius_idx;
        REAL star_radius;
        REAL *r_Schw_arr, *rho_arr, *rho_baryon_arr, *P_arr;
        REAL *M_arr, *expnu_arr, *exp4phi_arr, *rbar_arr;
        """,
        populate_ID_persist_struct_str="tov_reader(&ID_persist);",
        free_ID_persist_struct_str="tov_free(&ID_persist);",
        enable_T4munu=True,
    )
