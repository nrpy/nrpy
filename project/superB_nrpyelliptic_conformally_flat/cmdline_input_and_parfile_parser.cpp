#include "BHaH_defines.h"
#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <string.h>

#define NUM_PARAMETERS 30 // Define the number of parameters

#define LINE_SIZE 1024 // Define the max size of a line
#define PARAM_SIZE 128 // Define the max param string size

/**
 * Trims leading and trailing whitespace from a string.
 *
 * @param str - The input string to be trimmed.
 * @return - A pointer to the trimmed string.
 *
 * This function removes any leading and trailing whitespace characters from the input string.
 * It modifies the original string by inserting a null terminator after the last non-space character.
 */
static char *trim_space(char *str) {
  char *end;

  // Trim leading spaces
  while (isspace((unsigned char)*str))
    str++;

  if (*str == 0) // All spaces?
    return str;

  // Trim trailing spaces
  end = str + strlen(str) - 1;
  while (end > str && isspace((unsigned char)*end))
    end--;

  // Write new null terminator
  *(end + 1) = '\0';

  return str;
}

/**
 * Safely copies a source string to a destination buffer, preventing buffer overflows.
 *
 * @param dest - The destination buffer where the string will be copied.
 * @param src - The source string to be copied.
 * @param size - The size of the destination buffer.
 *
 * This function ensures that the source string fits within the destination buffer.
 * It checks for NULL pointers and verifies that the source string does not exceed the buffer size.
 * If any of these conditions fail, the function prints an error message and terminates the program.
 */
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

/**
 * Prints usage instructions for the program.
 *
 * This function outputs the different usage options available for running the
 * superB_nrpyelliptic_conformally_flat executable. It guides users on how to provide parameter
 * files and overwrite steerable parameters through command-line arguments.
 */
static void print_usage() {
  fprintf(stderr, "Usage option 0: ./superB_nrpyelliptic_conformally_flat [--help or -h] - Outputs this usage command\n");
  fprintf(stderr, "Usage option 1: ./superB_nrpyelliptic_conformally_flat - Reads in parameter file superB_nrpyelliptic_conformally_flat.par\n");
  fprintf(stderr, "Usage option 2: ./superB_nrpyelliptic_conformally_flat [parfile] - Reads in parameter file [parfile]\n");
  fprintf(stderr, "Usage option 3: ./superB_nrpyelliptic_conformally_flat [convergence_factor] - Overwrites parameters in list after reading in "
                  "superB_nrpyelliptic_conformally_flat.par\n");
  fprintf(stderr, "Usage option 4: ./superB_nrpyelliptic_conformally_flat [parfile] [convergence_factor] - Overwrites list of steerable parameters "
                  "after reading in [parfile]\n");
}

#define MAX_ARRAY_SIZE 100 // Adjust as needed

// Parameter types
typedef enum { PARAM_REAL, PARAM_INT, PARAM_BOOL, PARAM_CHARARRAY, PARAM_REAL_ARRAY, PARAM_INT_ARRAY } param_type;

// Parameter descriptor struct
typedef struct {
  const char *name;
  int index;
  param_type type;
  int array_size;     // For array parameters
  size_t buffer_size; // For char arrays
} param_descriptor;

// param_table[] is a list of parameter descriptors, each containing the parameter's name, unique index, type, array size, and buffer size.
param_descriptor param_table[] = {{"CFL_FACTOR", 0, PARAM_REAL, 0, 0},
                                  {"MINIMUM_GLOBAL_WAVESPEED", 1, PARAM_REAL, 0, 0},
                                  {"Nchare0", 2, PARAM_INT, 0, 0},
                                  {"Nchare1", 3, PARAM_INT, 0, 0},
                                  {"Nchare2", 4, PARAM_INT, 0, 0},
                                  {"P0_x", 5, PARAM_REAL, 0, 0},
                                  {"P0_y", 6, PARAM_REAL, 0, 0},
                                  {"P0_z", 7, PARAM_REAL, 0, 0},
                                  {"P1_x", 8, PARAM_REAL, 0, 0},
                                  {"P1_y", 9, PARAM_REAL, 0, 0},
                                  {"P1_z", 10, PARAM_REAL, 0, 0},
                                  {"S0_x", 11, PARAM_REAL, 0, 0},
                                  {"S0_y", 12, PARAM_REAL, 0, 0},
                                  {"S0_z", 13, PARAM_REAL, 0, 0},
                                  {"S1_x", 14, PARAM_REAL, 0, 0},
                                  {"S1_y", 15, PARAM_REAL, 0, 0},
                                  {"S1_z", 16, PARAM_REAL, 0, 0},
                                  {"bare_mass_0", 17, PARAM_REAL, 0, 0},
                                  {"bare_mass_1", 18, PARAM_REAL, 0, 0},
                                  {"convergence_factor", 19, PARAM_REAL, 0, 0},
                                  {"diagnostics_output_every", 20, PARAM_REAL, 0, 0},
                                  {"eta_damping", 21, PARAM_REAL, 0, 0},
                                  {"log10_current_residual", 22, PARAM_REAL, 0, 0},
                                  {"log10_residual_tolerance", 23, PARAM_REAL, 0, 0},
                                  {"nn_max", 24, PARAM_INT, 0, 0},
                                  {"outer_bc_type", 25, PARAM_CHARARRAY, 0, 50},
                                  {"output_progress_every", 26, PARAM_INT, 0, 0},
                                  {"stop_relaxation", 27, PARAM_BOOL, 0, 0},
                                  {"t_final", 28, PARAM_REAL, 0, 0},
                                  {"zPunc", 29, PARAM_REAL, 0, 0}};
#define NUM_PARAMS (int)(sizeof(param_table) / sizeof(param_descriptor))

/**
 * Searches for a parameter descriptor by its name.
 *
 * @param param_name - The name of the parameter to search for.
 * @return - A pointer to the corresponding param_descriptor if found; otherwise, NULL.
 *
 * This function iterates through the param_table to find a descriptor that matches the given
 * parameter name. It facilitates parameter validation and assignment during parsing.
 */
param_descriptor *find_param_descriptor(const char *param_name) {
  for (int i = 0; i < NUM_PARAMS; i++) {
    if (strcmp(param_table[i].name, param_name) == 0) {
      return &param_table[i];
    }
  }
  return NULL;
}

/**
 * Parses a parameter string to extract the parameter name and array size if applicable.
 *
 * @param param_str - The input parameter string, potentially containing array notation.
 * @param param_name - Buffer to store the extracted parameter name.
 * @param array_size - Pointer to store the extracted array size; set to 0 for scalar parameters.
 *
 * This function analyzes the parameter string to determine if it represents an array. If array
 * notation (e.g., param[10]) is detected, it extracts the base parameter name and the specified
 * array size. For scalar parameters, it copies the parameter name directly and sets the array
 * size to zero.
 */
static void parse_param(const char *param_str, char *param_name, int *array_size) {
  const char *bracket_start = strchr(param_str, '[');
  if (bracket_start != NULL) {
    // It's an array parameter
    size_t name_len = bracket_start - param_str;
    strncpy(param_name, param_str, name_len);
    param_name[name_len] = '\0';
    const char *bracket_end = strchr(bracket_start + 1, ']');
    if (bracket_end == NULL) {
      fprintf(stderr, "Error: Missing closing bracket in parameter %s.\n", param_str);
      exit(1);
    } // END IF (bracket_end == NULL): Check for missing closing bracket in parameter
    char size_str[16];
    size_t size_len = bracket_end - bracket_start - 1;
    strncpy(size_str, bracket_start + 1, size_len);
    size_str[size_len] = '\0';
    *array_size = atoi(size_str);
    if (*array_size <= 0) {
      fprintf(stderr, "Error: Invalid array size in parameter %s.\n", param_name);
      exit(1);
    } // END IF (*array_size <= 0): Validate that array size is positive
  } else {
    // Scalar parameter
    safe_copy(param_name, param_str, PARAM_SIZE);
    *array_size = 0;
  } // END IF (bracket_start != NULL): Distinguish between array and scalar parameter
} // END FUNCTION parse_param

// Function to parse value string into array of values
static void parse_value(const char *value_str, char values[][PARAM_SIZE], int *value_count) {
  if (value_str[0] == '{') {
    // Array value
    size_t len = strlen(value_str);
    if (value_str[len - 1] != '}') {
      fprintf(stderr, "Error: Missing closing brace in value %s.\n", value_str);
      exit(1);
    } // END IF (value_str[len - 1] != '}'): Check for missing closing brace

    // Extract the values inside the braces
    char value_copy[LINE_SIZE];
    safe_copy(value_copy, value_str + 1, LINE_SIZE); // Skip the opening brace
    value_copy[len - 2] = '\0';                      // Remove the closing brace
    // Now split value_copy by ','
    char *val_token;
    int count = 0;
    char *saveptr;
    val_token = strtok_r(value_copy, ",", &saveptr);
    while (val_token != NULL) {
      safe_copy(values[count], trim_space(val_token), PARAM_SIZE);
      count++;
      if (count > MAX_ARRAY_SIZE) {
        fprintf(stderr, "Error: Array size exceeds maximum allowed.\n");
        exit(1);
      } // END IF (count > MAX_ARRAY_SIZE): Check for exceeding maximum array size
      val_token = strtok_r(NULL, ",", &saveptr);
    } // END while (val_token != NULL): Split value_copy into tokens using comma as delimiter
    *value_count = count;
  } else {
    // Scalar value
    char mutable_value[PARAM_SIZE];
    safe_copy(mutable_value, value_str, PARAM_SIZE);
    char *trimmed = trim_space(mutable_value);
    safe_copy(values[0], trimmed, PARAM_SIZE);
    *value_count = 1;
  } // END IF (value_str[0] == '{'): Handle array and scalar values
} // END FUNCTION parse_value

// Function to read an integer value
static void read_integer(const char *value, int *result, const char *param_name) {
  char *endptr;
  errno = 0; // To detect overflow
  long int_val = strtol(value, &endptr, 10);

  if (endptr == value || *endptr != '\0' || errno == ERANGE) {
    fprintf(stderr, "Error: Invalid integer value for %s: %s.\n", param_name, value);
    exit(1);
  } // END IF (endptr == value || *endptr != '\0' || errno == ERANGE): Validate integer conversion and check for errors

  *result = (int)int_val;
} // END FUNCTION read_integer

// Function to read a REAL (double) value
static void read_REAL(const char *value, REAL *result, const char *param_name) {
  char *endptr;
  errno = 0; // To detect overflow
  double double_val = strtod(value, &endptr);

  if (endptr == value || *endptr != '\0' || errno == ERANGE) {
    fprintf(stderr, "Error: Invalid double value for %s: %s.\n", param_name, value);
    exit(1);
  } // END IF (endptr == value || *endptr != '\0' || errno == ERANGE): Validate double conversion and check for errors

  *result = (REAL)double_val;
} // END FUNCTION read_REAL

// Function to read a character array
static void read_chararray(const char *value, char *result, const char *param_name, size_t size) {
  // Remove surrounding quotes if present
  size_t len = strlen(value);
  char trimmed_value[PARAM_SIZE];
  if (value[0] == '"' && value[len - 1] == '"') {
    if (len - 2 >= size) {
      fprintf(stderr, "Error: Buffer overflow detected for %s.\n", param_name);
      exit(1);
    } // END IF (len - 2 >= size): Check for buffer overflow before trimming quotes
    strncpy(trimmed_value, value + 1, len - 2);
    trimmed_value[len - 2] = '\0';
  } else {
    safe_copy(trimmed_value, value, PARAM_SIZE);
  } // END IF: Determine if surrounding quotes exist and handle accordingly

  if (strlen(trimmed_value) >= size) {
    fprintf(stderr, "Error: Buffer overflow detected for %s.\n", param_name);
    exit(1);
  } // END IF (strlen(trimmed_value) >= size): Check for buffer overflow after trimming
  safe_copy(result, trimmed_value, size);
} // END FUNCTION read_chararray

// Function to read a boolean value
static void read_boolean(const char *value, bool *result, const char *param_name) {
  // To allow case-insensitive comparison
  char *lower_value = strdup(value);
  if (lower_value == NULL) {
    fprintf(stderr, "Error: Memory allocation failed for boolean value of %s.\n", param_name);
    exit(1);
  } // END memory allocation error check

  // Convert the string to lowercase
  for (char *p = lower_value; *p != '\0'; p++) {
    *p = tolower((unsigned char)*p);
  } // END lowercase conversion loop

  // Strip out any single or double quotation marks
  char *src = lower_value, *dst = lower_value;
  while (*src != '\0') {
    if (*src != '\"' && *src != '\'') {
      *dst++ = *src;
    } // END IF check for quotation marks
    src++;
  } // END while loop for stripping quotation marks
  *dst = '\0';

  // Check if the input is (case-insensitive) "true" or "false"; or "0" or "1"
  if (strcmp(lower_value, "true") == 0 || strcmp(lower_value, "yes") == 0 || strcmp(lower_value, "1") == 0) {
    *result = true;
  } else if (strcmp(lower_value, "false") == 0 || strcmp(lower_value, "no") == 0 || strcmp(lower_value, "0") == 0) {
    *result = false;
  } else {
    fprintf(stderr, "Error: Invalid boolean value for %s: %s.\n", param_name, value);
    free(lower_value);
    exit(1);
  } // END boolean value check (if-else chain)

  // Free the allocated memory for the lowercase copy of the value
  free(lower_value);
} // END FUNCTION read_boolean

/**
 * AUTOMATICALLY GENERATED BY cmdline_input_and_parfiles.py
 *
 * Reads and parses a parameter file to populate the commondata_struct.
 *
 * This function processes command-line arguments and reads parameters from a specified
 * parameter file or a default file. It supports various usage options, including displaying
 * help information, reading different parameter files, and overwriting steerable parameters
 * with provided convergence factors. The function ensures that all parameters are valid,
 * correctly formatted, and not duplicated.
 *
 * @param commondata - Pointer to the commondata_struct to be populated.
 * @param argc - The argument count from the command-line input.
 * @param argv - The argument vector containing command-line arguments.
 */
void cmdline_input_and_parfile_parser(commondata_struct *restrict commondata, int argc, const char *argv[]) {
  const int number_of_steerable_parameters = 1;

  int option;

  // Check for "-h" or "--help" options
  if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {
    print_usage();
    exit(0);
  } // END IF: Checking for help option.

  // Determine the usage option based on argc
  if (argc == 1) {
    option = 1; // Usage option 1: Process default parameter file "superB_nrpyelliptic_conformally_flat.par"
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
    option = 3; // Usage option 3: Overwrite steerable parameters after processing "superB_nrpyelliptic_conformally_flat.par"
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
  } // END IF-ELSE: Determining usage option.

  // fprintf(stderr, "Using option %d\n", option);

  // Determine the filename based on the selected option.
  const char *filename = (option == 1 || option == 3) ? "superB_nrpyelliptic_conformally_flat.par" : argv[1];
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    print_usage();
    exit(1);
  } // END IF: Checking if the parameter file opened successfully.

  char line[LINE_SIZE];
  char param[PARAM_SIZE];
  char value[PARAM_SIZE];
  int params_set[NUM_PARAMETERS] = {0}; // Track which parameters have been set.

  // Iterate through each line of the parameter file.
  while (fgets(line, sizeof(line), file)) {
    // Remove comments from the line.
    char *comment_start = strchr(line, '#');
    if (comment_start != NULL) {
      *comment_start = '\0';
    }

    // Split the line into parameter and value based on '=' delimiter.
    char *p = strtok(line, "=");
    if (p) {
      safe_copy(param, trim_space(p), sizeof(param));
      p = strtok(NULL, "=");
      if (p) {
        safe_copy(value, trim_space(p), sizeof(value));

        char param_name[PARAM_SIZE];
        int array_size = 0;
        parse_param(param, param_name, &array_size);

        // Validate characters in the parameter name.
        for (int i = 0; param_name[i]; i++) {
          if (!isalnum(param_name[i]) && param_name[i] != '_') {
            fprintf(stderr, "Error: Invalid parameter name %s.\n", param_name);
            exit(1);
          }
        }

        // Parse the value string into individual values.
        char values_array[MAX_ARRAY_SIZE][PARAM_SIZE];
        int value_count = 0;
        parse_value(value, values_array, &value_count);

        // Lookup the parameter descriptor.
        param_descriptor *param_desc = find_param_descriptor(param_name);
        if (param_desc == NULL) {
          fprintf(stderr, "Warning: Unrecognized parameter %s.\n", param_name);
          continue; // Decide whether to exit or ignore
        }

        // Check for duplicate parameter definitions.
        if (params_set[param_desc->index] == 1) {
          fprintf(stderr, "Error: Duplicate parameter %s.\n", param_name);
          exit(1);
        }
        params_set[param_desc->index] = 1;

        // Validate array size if applicable.
        if (param_desc->type != PARAM_CHARARRAY && param_desc->array_size > 0) {
          // It's an array parameter
          if (array_size != param_desc->array_size) {
            fprintf(stderr, "Error: Array size mismatch for parameter %s.\n", param_name);
            exit(1);
          }
          if (value_count != param_desc->array_size) {
            fprintf(stderr, "Error: Number of values does not match array size for parameter %s.\n", param_name);
            exit(1);
          }
        } else {
          // It's a scalar parameter, including PARAM_CHARARRAY
          if (array_size > 0) {
            fprintf(stderr, "Error: Unexpected array size for scalar parameter %s.\n", param_name);
            exit(1);
          }
          if (value_count != 1) {
            fprintf(stderr, "Error: Expected a single value for parameter %s.\n", param_name);
            exit(1);
          }
        }

        // Assign the parsed values to the corresponding fields in commondata.

        if (param_desc->index == 0) {
          read_REAL(values_array[0], &commondata->CFL_FACTOR, "CFL_FACTOR");
        } else if (param_desc->index == 1) {
          read_REAL(values_array[0], &commondata->MINIMUM_GLOBAL_WAVESPEED, "MINIMUM_GLOBAL_WAVESPEED");
        } else if (param_desc->index == 2) {
          read_integer(values_array[0], &commondata->Nchare0, "Nchare0");
        } else if (param_desc->index == 3) {
          read_integer(values_array[0], &commondata->Nchare1, "Nchare1");
        } else if (param_desc->index == 4) {
          read_integer(values_array[0], &commondata->Nchare2, "Nchare2");
        } else if (param_desc->index == 5) {
          read_REAL(values_array[0], &commondata->P0_x, "P0_x");
        } else if (param_desc->index == 6) {
          read_REAL(values_array[0], &commondata->P0_y, "P0_y");
        } else if (param_desc->index == 7) {
          read_REAL(values_array[0], &commondata->P0_z, "P0_z");
        } else if (param_desc->index == 8) {
          read_REAL(values_array[0], &commondata->P1_x, "P1_x");
        } else if (param_desc->index == 9) {
          read_REAL(values_array[0], &commondata->P1_y, "P1_y");
        } else if (param_desc->index == 10) {
          read_REAL(values_array[0], &commondata->P1_z, "P1_z");
        } else if (param_desc->index == 11) {
          read_REAL(values_array[0], &commondata->S0_x, "S0_x");
        } else if (param_desc->index == 12) {
          read_REAL(values_array[0], &commondata->S0_y, "S0_y");
        } else if (param_desc->index == 13) {
          read_REAL(values_array[0], &commondata->S0_z, "S0_z");
        } else if (param_desc->index == 14) {
          read_REAL(values_array[0], &commondata->S1_x, "S1_x");
        } else if (param_desc->index == 15) {
          read_REAL(values_array[0], &commondata->S1_y, "S1_y");
        } else if (param_desc->index == 16) {
          read_REAL(values_array[0], &commondata->S1_z, "S1_z");
        } else if (param_desc->index == 17) {
          read_REAL(values_array[0], &commondata->bare_mass_0, "bare_mass_0");
        } else if (param_desc->index == 18) {
          read_REAL(values_array[0], &commondata->bare_mass_1, "bare_mass_1");
        } else if (param_desc->index == 19) {
          read_REAL(values_array[0], &commondata->convergence_factor, "convergence_factor");
        } else if (param_desc->index == 20) {
          read_REAL(values_array[0], &commondata->diagnostics_output_every, "diagnostics_output_every");
        } else if (param_desc->index == 21) {
          read_REAL(values_array[0], &commondata->eta_damping, "eta_damping");
        } else if (param_desc->index == 22) {
          read_REAL(values_array[0], &commondata->log10_current_residual, "log10_current_residual");
        } else if (param_desc->index == 23) {
          read_REAL(values_array[0], &commondata->log10_residual_tolerance, "log10_residual_tolerance");
        } else if (param_desc->index == 24) {
          read_integer(values_array[0], &commondata->nn_max, "nn_max");
        } else if (param_desc->index == 25) {
          read_chararray(values_array[0], commondata->outer_bc_type, "outer_bc_type", 50);
        } else if (param_desc->index == 26) {
          read_integer(values_array[0], &commondata->output_progress_every, "output_progress_every");
        } else if (param_desc->index == 27) {
          read_boolean(values_array[0], &commondata->stop_relaxation, "stop_relaxation");
        } else if (param_desc->index == 28) {
          read_REAL(values_array[0], &commondata->t_final, "t_final");
        } else if (param_desc->index == 29) {
          read_REAL(values_array[0], &commondata->zPunc, "zPunc");
        } else {
          fprintf(stderr, "Error: Unknown parameter index for %s.\n", param_name);
          exit(1);
        }
      }
    }

    // Handling options 3 and 4: Overwriting steerable parameters
    if (option == 3 || option == 4) {
      // For options 3 and 4, we extract the last arguments as steerable parameters
      read_REAL(argv[argc - number_of_steerable_parameters + 0], &commondata->convergence_factor, "convergence_factor");
    } // END IF (option == 3 || option == 4)
  } // END WHILE LOOP over all lines in the file
  fclose(file); // Close the parameter file.
} // END FUNCTION cmdline_input_and_parfile_parser
