"""
Python module for parsing command-line arguments for TwoPunctures initial data

License: Lesser GNU Public License, version 2.0+

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
import nrpy.c_function as cfc


def add_to_Cfunction_dict__TwoPunctures_parse_cmdline_arguments():
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Parse command-line arguments for TwoPunctures initial data"""
    c_type = "void"
    name = "TwoPunctures_parse_cmdline_arguments"
    params = "const int argc, const char *argv[], ID_persist_struct *restrict ID_persist, commondata_struct *restrict commondata"
    body = r"""
  // Step 0.d: Set ID_type and TP_ID_type
  char TP_ID_type[100];
  sprintf(TP_ID_type, "NRPyPN");


  // Default BBH parameters: q=1 nonspinning, separation=2.0
  commondata->bbh_physical_params.initial_orbital_separation = 2.0;
  commondata->bbh_physical_params.mass_ratio = 1.0;
  for(int i=0;i<3;i++) commondata->bbh_physical_params.chi_BH_m[i] =
                         commondata->bbh_physical_params.chi_BH_M[i] = 0.0;

  // Step 0.e: Read command-line input if not reading from checkpoint; error out if nonconformant
  if(argc == 1) {
    // do nothing; using default parameters :)
  } else if(argc == 2) {
    snprintf(TP_ID_type, 100, "%s", argv[1]);
  } else if(argc == 9 || argc == 10) {
    // PUNCTURES ORBIT IN THE z-x PLANE; y-AXIS IS PERPENDICULAR TO INITIAL ORBITAL PLANE
    int ww=1;
    commondata->bbh_physical_params.initial_orbital_separation = strtod(argv[ww], NULL); ww++;
    commondata->bbh_physical_params.mass_ratio                 = strtod(argv[ww], NULL); ww++;
    if(commondata->bbh_physical_params.mass_ratio < 1.0) {
      fprintf(stderr, "ERROR: we use the mass ratio convention such that mass_ratio = q = M/m >= 1.0. q=%.5f is invalid\n", commondata->bbh_physical_params.mass_ratio);
      exit(1);
    }

    // The inputs assume the BHs are initially (instantaneously) orbiting on the xy plane.
    //   This is by convention. Our actual simulations are performed on the zx plane,
    //   and conversions are performed in TwoPunctures_BBH_params_library__set_BBH_params.c
    commondata->bbh_physical_params.chi_BH_M[0]                = strtod(argv[ww], NULL); ww++;
    commondata->bbh_physical_params.chi_BH_M[1]                = strtod(argv[ww], NULL); ww++;
    commondata->bbh_physical_params.chi_BH_M[2]                = strtod(argv[ww], NULL); ww++;
    commondata->bbh_physical_params.chi_BH_m[0]                = strtod(argv[ww], NULL); ww++;
    commondata->bbh_physical_params.chi_BH_m[1]                = strtod(argv[ww], NULL); ww++;
    commondata->bbh_physical_params.chi_BH_m[2]                = strtod(argv[ww], NULL); ww++;
    if(argc == 10) sprintf(TP_ID_type, "TP_BL_custom");
    if(commondata->bbh_physical_params.initial_orbital_separation <= 0 ||
       commondata->bbh_physical_params.initial_orbital_separation >= 30) {
      fprintf(stderr, "Error: commondata->bbh_physical_params.initial_orbital_separation = %e is not supported.\n",
              commondata->bbh_physical_params.initial_orbital_separation);
      exit(1);
    }
  } else {
    fprintf(stderr, "Error: Expected 0, 1, or 8 command-line arguments: ./BlackHolesatHome_Playground (default parameters),\n");
    fprintf(stderr, " 0 arguments: ./BlackHolesatHome_Playground <- (default parameters in executable)\n");
    fprintf(stderr, " 1 argument:  ./BlackHolesatHome_Playground [description] <- description can be QC0, GW150914ET, etc.\n");
    fprintf(stderr, " 8 (9) arguments: ./BlackHolesatHome_Playground initial_diameter_of_separation q m_chi_x m_chi_y m_chi_z M_chi_x M_chi_y M_chi_z (ctp; optional if want cheap TP)\n");
    exit(1);
  }
  TwoPunctures_BBH_params_library__set_BBH_params(commondata, TP_ID_type, ID_persist);

  fprintf(stderr, "#################################\n");
  fprintf(stderr, "-={ INITIAL BINARY PARAMETERS }=-\n");
  fprintf(stderr, "M=1 (sum of individual ADM masses as defined in TwoPunctures)\n");
  fprintf(stderr, "d_initial/M = %.15f, q = %.15f\n", commondata->bbh_physical_params.initial_orbital_separation,commondata->bbh_physical_params.mass_ratio);
  fprintf(stderr, "chi_BH_m = %.15f %.15f %.15f\n", commondata->bbh_physical_params.chi_BH_m[0],commondata->bbh_physical_params.chi_BH_m[1],commondata->bbh_physical_params.chi_BH_m[2]);
  fprintf(stderr, "chi_BH_M = %.15f %.15f %.15f\n", commondata->bbh_physical_params.chi_BH_M[0],commondata->bbh_physical_params.chi_BH_M[1],commondata->bbh_physical_params.chi_BH_M[2]);
  fprintf(stderr, "p_t = %.15f, p_r = %.15f\n", commondata->bbh_physical_params.initial_p_t,commondata->bbh_physical_params.initial_p_r);
  fprintf(stderr, "TP resolution: %d  %d  %d\n", ID_persist->npoints_A,ID_persist->npoints_B,ID_persist->npoints_phi);
  fprintf(stderr, "#################################\n");
"""
    cfc.register_CFunction(
        includes=includes, desc=desc, c_type=c_type, name=name, params=params, body=body
    )
