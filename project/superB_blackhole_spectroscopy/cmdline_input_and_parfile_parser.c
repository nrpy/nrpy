#include "BHaH_defines.h"
#include <ctype.h>
#include <errno.h>
#include <string.h>
#define NUM_PARAMETERS 27 // Define the number of parameters
#define LINE_SIZE 1024    // Define the max size of a line
#define PARAM_SIZE 128    // Define the max param string size

static char *trim_space(char *str) {
  char *end;

  // Trim leading spaces
  while (isspace((unsigned char)*str))
    str++;

  // Trim trailing spaces
  end = str + strlen(str) - 1;
  while (end > str && isspace((unsigned char)*end))
    end--;

  *(end + 1) = '\0';

  return str;
}

static void safe_copy(char *dest, const char *src, size_t size) {
  if (src == NULL) {
    fprintf(stderr, "Error: Source string is NULL.\n");
    exit(1);
  }
  if (dest == NULL) {
    fprintf(stderr, "Error: Destination string is NULL.\n");
    exit(1);
  }
  if (size == 0) {
    fprintf(stderr, "Error: Size is zero.\n");
    exit(1);
  }
  size_t src_len = strlen(src);
  if (src_len >= size) {
    fprintf(stderr, "Error: Buffer overflow detected.\n");
    exit(1);
  }
  strncpy(dest, src, size - 1);
  dest[size - 1] = '\0';
}
// Function to print usage instructions
static void print_usage() {
  fprintf(stderr, "Usage option 0: ./superB_blackhole_spectroscopy [--help or -h] <- Outputs this usage command\n");
  fprintf(stderr, "Usage option 1: ./superB_blackhole_spectroscopy <- reads in parameter file superB_blackhole_spectroscopy.par\n");
  fprintf(stderr, "Usage option 2: ./superB_blackhole_spectroscopy [parfile] <- reads in parameter file [parfile]\n");
  fprintf(stderr, "Usage option 3: ./superB_blackhole_spectroscopy [convergence_factor] <- overwrites parameters in list after reading in "
                  "superB_blackhole_spectroscopy.par\n");
  fprintf(stderr, "Usage option 4: ./superB_blackhole_spectroscopy [parfile] [convergence_factor] <- overwrites list of steerable parameters after "
                  "reading in [parfile]\n");
}
static void read_integer(const char *value, int *result, const char *param_name) {
  char *endptr;
  errno = 0; // To detect overflow
  long int_val = strtol(value, &endptr, 10);

  if (endptr == value || *endptr != '\0' || errno == ERANGE) {
    fprintf(stderr, "Error: Invalid integer value for %s: %s.\n", param_name, value);
    exit(1);
  }

  *result = (int)int_val;
}

static void read_double(const char *value, double *result, const char *param_name) {
  char *endptr;
  errno = 0; // To detect overflow
  double double_val = strtod(value, &endptr);

  if (endptr == value || *endptr != '\0' || errno == ERANGE) {
    fprintf(stderr, "Error: Invalid double value for %s: %s.\n", param_name, value);
    exit(1);
  }

  *result = double_val;
}

static void read_chararray(const char *value, char *result, const char *param_name, size_t size) {
  if (strlen(value) >= size) {
    fprintf(stderr, "Error: Buffer overflow detected for %s.\n", param_name);
    exit(1);
  }
  safe_copy(result, value, size);
}

/*
 * AUTOMATICALLY GENERATED BY parameter_file_read_and_parse.py
 * parameter_file_read_and_parse() function:
 * Reads and parses a parameter file to populate commondata_struct commondata.
 *
 * This function takes in the command-line arguments and a pointer to a commondata_struct.
 * It reads the provided file and extracts the parameters defined in the file, populating
 * the commondata_struct with the values. The file is expected to contain key-value pairs
 * separated by an equals sign (=), and it may include comments starting with a hash (#).
 * The function handles errors such as file opening failure, duplicate parameters, and
 * invalid parameter names.
 *
 * @param griddata_params: Pointer to the commondata struct to be populated.
 * @param argc: The argument count from the command-line input.
 * @param argv: The argument vector containing command-line arguments.
 */
void cmdline_input_and_parfile_parser(commondata_struct *restrict commondata, int argc, const char *argv[]) {
  const int number_of_steerable_parameters = 1;

  int option;

  // Check for "-h" or "--help" options
  if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {
    print_usage();
    exit(0);
  }

  // Determine the usage option based on argc
  if (argc == 1) {
    option = 1; // Usage option 1: Process default parameter file "superB_blackhole_spectroscopy.par"
  } else if (argc == 2) {
    // Check if the argument is a file
    FILE *file_check = fopen(argv[1], "r");
    if (file_check != NULL) {
      fclose(file_check);
      option = 2; // Usage option 2: Process parameter file provided in argv[1]
    } else if (argc == 1 + number_of_steerable_parameters) {
      option = 3;
    } else {
      fprintf(stderr, "Error: Invalid number of arguments or file cannot be opened.\n");
      print_usage();
      exit(1);
    }
  } else if (argc == 1 + number_of_steerable_parameters) {
    option = 3; // Usage option 3: Overwrite steerable parameters after processing "superB_blackhole_spectroscopy.par"
  } else if (argc == 2 + number_of_steerable_parameters) {
    // Check if the first argument is a file
    FILE *file_check = fopen(argv[1], "r");
    if (file_check != NULL) {
      fclose(file_check);
      option = 4; // Usage option 4: Overwrite steerable parameters after processing parameter file provided in argv[1]
    } else {
      fprintf(stderr, "Error: File cannot be opened for option 4.\n");
      print_usage();
      exit(1);
    }
  } else {
    fprintf(stderr, "Error: Invalid number of arguments\n");
    print_usage();
    exit(1);
  }

  // fprintf(stderr, "Using option %d\n", option);

  const char *filename = (option == 1 || option == 3) ? "superB_blackhole_spectroscopy.par" : argv[1];
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    print_usage();
    exit(1);
  }

  char line[LINE_SIZE];
  char param[PARAM_SIZE];
  char value[PARAM_SIZE];
  int params_set[NUM_PARAMETERS] = {0}; // Record of parameters set (one for each parameter in the struct)

  while (fgets(line, sizeof(line), file)) {
    // Removing comments from the line
    char *comment_start = strchr(line, '#');
    if (comment_start != NULL) {
      *comment_start = '\0';
    }

    char *p = strtok(line, "=");
    if (p) {
      safe_copy(param, trim_space(p), sizeof(param));
      p = strtok(NULL, "=");
      if (p) {
        safe_copy(value, trim_space(p), sizeof(value));

        // Check for naming convention violations
        for (int i = 0; param[i]; i++) {
          if (!isalnum(param[i]) && param[i] != '_') {
            fprintf(stderr, "Error: Invalid parameter name %s.\n", param);
            exit(1);
          }
        }

        int param_index = -1;
        if (1 == 0)
          param_index = -2;
        else if (strcmp(param, "CFL_FACTOR") == 0)
          param_index = 0;
        else if (strcmp(param, "KreissOliger_strength_gauge") == 0)
          param_index = 1;
        else if (strcmp(param, "KreissOliger_strength_nongauge") == 0)
          param_index = 2;
        else if (strcmp(param, "NUMGRIDS") == 0)
          param_index = 3;
        else if (strcmp(param, "TP_BBH_description") == 0)
          param_index = 4;
        else if (strcmp(param, "TP_bare_mass_M") == 0)
          param_index = 5;
        else if (strcmp(param, "TP_bare_mass_m") == 0)
          param_index = 6;
        else if (strcmp(param, "TP_npoints_A") == 0)
          param_index = 7;
        else if (strcmp(param, "TP_npoints_B") == 0)
          param_index = 8;
        else if (strcmp(param, "TP_npoints_phi") == 0)
          param_index = 9;
        else if (strcmp(param, "bbhxy_BH_M_chix") == 0)
          param_index = 10;
        else if (strcmp(param, "bbhxy_BH_M_chiy") == 0)
          param_index = 11;
        else if (strcmp(param, "bbhxy_BH_M_chiz") == 0)
          param_index = 12;
        else if (strcmp(param, "bbhxy_BH_m_chix") == 0)
          param_index = 13;
        else if (strcmp(param, "bbhxy_BH_m_chiy") == 0)
          param_index = 14;
        else if (strcmp(param, "bbhxy_BH_m_chiz") == 0)
          param_index = 15;
        else if (strcmp(param, "checkpoint_every") == 0)
          param_index = 16;
        else if (strcmp(param, "convergence_factor") == 0)
          param_index = 17;
        else if (strcmp(param, "diagnostics_output_every") == 0)
          param_index = 18;
        else if (strcmp(param, "eta") == 0)
          param_index = 19;
        else if (strcmp(param, "initial_p_r") == 0)
          param_index = 20;
        else if (strcmp(param, "initial_p_t") == 0)
          param_index = 21;
        else if (strcmp(param, "initial_sep") == 0)
          param_index = 22;
        else if (strcmp(param, "mass_ratio") == 0)
          param_index = 23;
        else if (strcmp(param, "outer_bc_type") == 0)
          param_index = 24;
        else if (strcmp(param, "swm2sh_maximum_l_mode_to_compute") == 0)
          param_index = 25;
        else if (strcmp(param, "t_final") == 0)
          param_index = 26;
        else
          fprintf(stderr, "Warning: Unrecognized parameter %s.\n", param);

        // Check for duplicates
        if (param_index != -1 && params_set[param_index] == 1) {
          fprintf(stderr, "Error: Duplicate parameter %s.\n", param);
          exit(1);
        }
        if (param_index != -1)
          params_set[param_index] = 1;

        // Assign values
        if (param_index == -2)
          exit(1); // impossible.
        else if (param_index == 0)
          read_double(value, &commondata->CFL_FACTOR, "CFL_FACTOR");
        else if (param_index == 1)
          read_double(value, &commondata->KreissOliger_strength_gauge, "KreissOliger_strength_gauge");
        else if (param_index == 2)
          read_double(value, &commondata->KreissOliger_strength_nongauge, "KreissOliger_strength_nongauge");
        else if (param_index == 3)
          read_integer(value, &commondata->NUMGRIDS, "NUMGRIDS");
        else if (param_index == 4)
          read_chararray(value, commondata->TP_BBH_description, "TP_BBH_description", 100);
        else if (param_index == 5)
          read_double(value, &commondata->TP_bare_mass_M, "TP_bare_mass_M");
        else if (param_index == 6)
          read_double(value, &commondata->TP_bare_mass_m, "TP_bare_mass_m");
        else if (param_index == 7)
          read_integer(value, &commondata->TP_npoints_A, "TP_npoints_A");
        else if (param_index == 8)
          read_integer(value, &commondata->TP_npoints_B, "TP_npoints_B");
        else if (param_index == 9)
          read_integer(value, &commondata->TP_npoints_phi, "TP_npoints_phi");
        else if (param_index == 10)
          read_double(value, &commondata->bbhxy_BH_M_chix, "bbhxy_BH_M_chix");
        else if (param_index == 11)
          read_double(value, &commondata->bbhxy_BH_M_chiy, "bbhxy_BH_M_chiy");
        else if (param_index == 12)
          read_double(value, &commondata->bbhxy_BH_M_chiz, "bbhxy_BH_M_chiz");
        else if (param_index == 13)
          read_double(value, &commondata->bbhxy_BH_m_chix, "bbhxy_BH_m_chix");
        else if (param_index == 14)
          read_double(value, &commondata->bbhxy_BH_m_chiy, "bbhxy_BH_m_chiy");
        else if (param_index == 15)
          read_double(value, &commondata->bbhxy_BH_m_chiz, "bbhxy_BH_m_chiz");
        else if (param_index == 16)
          read_double(value, &commondata->checkpoint_every, "checkpoint_every");
        else if (param_index == 17)
          read_double(value, &commondata->convergence_factor, "convergence_factor");
        else if (param_index == 18)
          read_double(value, &commondata->diagnostics_output_every, "diagnostics_output_every");
        else if (param_index == 19)
          read_double(value, &commondata->eta, "eta");
        else if (param_index == 20)
          read_double(value, &commondata->initial_p_r, "initial_p_r");
        else if (param_index == 21)
          read_double(value, &commondata->initial_p_t, "initial_p_t");
        else if (param_index == 22)
          read_double(value, &commondata->initial_sep, "initial_sep");
        else if (param_index == 23)
          read_double(value, &commondata->mass_ratio, "mass_ratio");
        else if (param_index == 24)
          read_chararray(value, commondata->outer_bc_type, "outer_bc_type", 50);
        else if (param_index == 25)
          read_integer(value, &commondata->swm2sh_maximum_l_mode_to_compute, "swm2sh_maximum_l_mode_to_compute");
        else if (param_index == 26)
          read_double(value, &commondata->t_final, "t_final");
        else {
          fprintf(stderr, "Error: Unrecognized parameter %s.\n", param);
          exit(1); // Exit on unrecognized parameter
        }
      }
    }
  }

  fclose(file);
  // Handling options 3 and 4: Overwriting steerable parameters
  if (option == 3 || option == 4) {
    // For options 3 and 4, we extract the last three arguments as steerable parameters
    read_double(argv[argc - number_of_steerable_parameters + 0], &commondata->convergence_factor, "convergence_factor");
  }
}
