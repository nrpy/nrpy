"""
Functions for parsing command-line input and parameter files.

This module provides functions to auto-generate C code for parsing parameter
files and to create a default parameter file based on registered CodeParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Optional, Set, Tuple

import nrpy.c_function as cfc
import nrpy.params as par

# region C Code Snippets for cmdline_input_and_parfile_parser
# Each C helper function and template is isolated in its own constant for maximum clarity and modularity.

_C_DEFINES = r"""#define LINE_SIZE 1024 // Define the max size of a line
#define PARAM_SIZE 128 // Define the max param string size
"""

_C_TRIM_SPACE_FUNC = r"""
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
"""

_C_SAFE_COPY_FUNC = r"""
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
"""

_C_PARAM_DEFS = r"""
#define MAX_ARRAY_SIZE 100 // Adjust as needed

// Parameter types
typedef enum { PARAM_REAL, PARAM_INT, PARAM_BOOL, PARAM_CHARARRAY, PARAM_REAL_ARRAY, PARAM_INT_ARRAY } param_type;

// Parameter descriptor struct
typedef struct {
    const char *name;
    int index;
    param_type type;
    int array_size;    // For array parameters
    size_t buffer_size; // For char arrays
} param_descriptor;
"""

_C_FIND_PARAM_DESCRIPTOR_FUNC = r"""
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
"""

_C_PARSE_PARAM_FUNC = r"""
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
"""

_C_PARSE_VALUE_FUNC = r"""
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
"""

_C_READ_INTEGER_FUNC = r"""
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
"""

_C_READ_REAL_FUNC = r"""
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
"""

_C_READ_CHARARRAY_FUNC = r"""
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
"""

_C_READ_BOOLEAN_FUNC = r"""
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
"""

# endregion

# region C Code Generation for cmdline_input_and_parfile_parser


class CParameter(NamedTuple):
    """
    A structured representation of a parameter for C code generation.

    This provides a clear data contract between analysis and generation functions.
    """

    param_index: int
    name: str
    enum_type: str
    array_size: int
    buffer_size: int
    primitive_type: str
    original_cparam_type: str


def _process_parameters_for_c_parser() -> Tuple[List[CParameter], Set[str]]:
    """
    Analyze global CodeParameters and translate them into a structured list for C code generation.

    :return: A tuple containing:
             - A list of CParameter objects.
             - A set of strings indicating which C 'read_*' functions are needed.
    """
    parameters_list: List[CParameter] = []
    needed_readers: Set[str] = set()
    param_index = 0

    sorted_params = sorted(par.glb_code_params_dict.items())
    for name, codeparam in sorted_params:
        if not (codeparam.add_to_parfile and codeparam.commondata):
            continue

        cparam_type = codeparam.cparam_type.strip()
        info: Dict[str, Any] = {"array_size": 0, "buffer_size": 0}

        if "int[" in cparam_type or "int [" in cparam_type:
            info.update(
                enum_type="PARAM_INT_ARRAY",
                array_size=int(cparam_type.split("[")[1].split("]")[0]),
                primitive_type="integer",
            )
        elif "REAL[" in cparam_type or "REAL [" in cparam_type:
            info.update(
                enum_type="PARAM_REAL_ARRAY",
                array_size=int(cparam_type.split("[")[1].split("]")[0]),
                primitive_type="REAL",
            )
        elif cparam_type == "int":
            info.update(enum_type="PARAM_INT", primitive_type="integer")
        elif cparam_type == "REAL":
            info.update(enum_type="PARAM_REAL", primitive_type="REAL")
        elif cparam_type == "bool":
            info.update(enum_type="PARAM_BOOL", primitive_type="boolean")
        elif "char" in cparam_type and "[" in cparam_type and "]" in cparam_type:
            info.update(
                enum_type="PARAM_CHARARRAY",
                buffer_size=int(cparam_type.split("[")[1].split("]")[0]),
                primitive_type="chararray",
            )
        else:
            continue  # Skip unsupported types

        needed_readers.add(info["primitive_type"])
        parameters_list.append(
            CParameter(
                param_index=param_index,
                name=name,
                enum_type=info["enum_type"],
                array_size=info["array_size"],
                buffer_size=info["buffer_size"],
                primitive_type=info["primitive_type"],
                original_cparam_type=cparam_type,
            )
        )
        param_index += 1

    return parameters_list, needed_readers


def _generate_c_usage_function(project_name: str, cmdline_inputs: List[str]) -> str:
    """
    Generate the C code for the print_usage() function.

    :param project_name: The name of the project.
    :param cmdline_inputs: A list of steerable command-line inputs.
    :return: A string containing the C print_usage() function.
    """
    list_of_steerable_params_str = " ".join(cmdline_inputs)
    return rf"""
/**
 * Prints usage instructions for the program.
 *
 * This function outputs the different usage options available for running the
 * {project_name} executable. It guides users on how to provide parameter
 * files and overwrite steerable parameters through command-line arguments.
 */
static void print_usage() {{
  fprintf(stderr, "Usage option 0: ./{project_name} [--help or -h] - Outputs this usage command\n");
  fprintf(stderr, "Usage option 1: ./{project_name} - Reads in parameter file {project_name}.par\n");
  fprintf(stderr, "Usage option 2: ./{project_name} [parfile] - Reads in parameter file [parfile]\n");
  fprintf(stderr, "Usage option 3: ./{project_name} [{list_of_steerable_params_str}] - Overwrites parameters in list after reading in {project_name}.par\n");
  fprintf(stderr, "Usage option 4: ./{project_name} [parfile] [{list_of_steerable_params_str}] - Overwrites list of steerable parameters after reading in [parfile]\n");
}}"""


def _generate_c_param_table(parameters_list: List[CParameter]) -> str:
    """
    Generate the C code for the param_descriptor table.

    :param parameters_list: A list of CParameter objects.
    :return: A string containing the C param_descriptor table definition.
    """
    param_table_entries = [
        f'    {{"{p.name}", {p.param_index}, {p.enum_type}, {p.array_size}, {p.buffer_size}}}'
        for p in parameters_list
    ]
    return (
        "param_descriptor param_table[] = {\n"
        + ",\n".join(param_table_entries)
        + "\n};"
    )


def _generate_c_preamble(
    project_name: str,
    cmdline_inputs: List[str],
    parameters_list: List[CParameter],
    needed_readers: Set[str],
) -> str:
    """
    Assemble the entire C preamble, including all helper functions and definitions.

    :param project_name: The name of the project.
    :param cmdline_inputs: A list of steerable command-line inputs.
    :param parameters_list: A list of CParameter objects.
    :param needed_readers: A set of strings indicating required 'read_*' functions.
    :return: A string containing the complete C preamble.
    """
    num_commondata_params = len(parameters_list)
    prefunc_parts = [
        f"#define NUM_PARAMETERS {num_commondata_params} // Define the number of parameters\n"
    ]
    prefunc_parts.append(_C_DEFINES)
    prefunc_parts.append(_C_TRIM_SPACE_FUNC)
    prefunc_parts.append(_C_SAFE_COPY_FUNC)
    prefunc_parts.append(_generate_c_usage_function(project_name, cmdline_inputs))
    prefunc_parts.append(_C_PARAM_DEFS)
    prefunc_parts.append(
        "// param_table[] is a list of parameter descriptors, each containing the parameter's name, unique index, type, array size, and buffer size."
    )
    prefunc_parts.append(_generate_c_param_table(parameters_list))
    prefunc_parts.append(
        "#define NUM_PARAMS (int)(sizeof(param_table) / sizeof(param_descriptor))\n"
    )
    prefunc_parts.append(_C_FIND_PARAM_DESCRIPTOR_FUNC)
    prefunc_parts.append(_C_PARSE_PARAM_FUNC)
    prefunc_parts.append(_C_PARSE_VALUE_FUNC)

    # Conditionally add the required C 'read_*' functions in the original order.
    reader_map = {
        "integer": _C_READ_INTEGER_FUNC,
        "REAL": _C_READ_REAL_FUNC,
        "chararray": _C_READ_CHARARRAY_FUNC,
        "boolean": _C_READ_BOOLEAN_FUNC,
    }
    canonical_reader_order = ["integer", "REAL", "chararray", "boolean"]
    for reader_key in canonical_reader_order:
        if reader_key in needed_readers:
            prefunc_parts.append(reader_map[reader_key])

    return "\n".join(prefunc_parts)


def _generate_c_body(
    project_name: str, cmdline_inputs: List[str], parameters_list: List[CParameter]
) -> str:
    """
    Generate the main body of the C parser function by filling in a template.

    :param project_name: The name of the project.
    :param cmdline_inputs: A list of steerable command-line inputs.
    :param parameters_list: A list of CParameter objects.
    :return: A string containing the body of the main C parser function.
    """
    assignment_code = ""
    first_param = True
    for p in parameters_list:
        prefix = "if" if first_param else "else if"
        first_param = False
        action = ""
        if p.enum_type == "PARAM_REAL":
            action = f'read_REAL(values_array[0], &commondata->{p.name}, "{p.name}");'
        elif p.enum_type == "PARAM_INT":
            action = (
                f'read_integer(values_array[0], &commondata->{p.name}, "{p.name}");'
            )
        elif p.enum_type == "PARAM_BOOL":
            action = (
                f'read_boolean(values_array[0], &commondata->{p.name}, "{p.name}");'
            )
        elif p.enum_type == "PARAM_CHARARRAY":
            action = f'read_chararray(values_array[0], commondata->{p.name}, "{p.name}", {p.buffer_size});'
        elif p.enum_type == "PARAM_REAL_ARRAY":
            action = f'for (int i = 0; i < {p.array_size}; i++) {{\n            read_REAL(values_array[i], &commondata->{p.name}[i], "{p.name}");\n        }}'
        elif p.enum_type == "PARAM_INT_ARRAY":
            action = f'for (int i = 0; i < {p.array_size}; i++) {{\n            read_integer(values_array[i], &commondata->{p.name}[i], "{p.name}");\n        }}'

        assignment_code += f"""
    {prefix} (param_desc->index == {p.param_index}) {{
        {action}
    }}"""

    steerable_c_code = ""
    for i, key in enumerate(cmdline_inputs):
        c_param = next((p for p in parameters_list if p.name == key), None)
        if not c_param:
            continue

        arg = f"argv[argc - number_of_steerable_parameters + {i}]"
        line = ""
        if c_param.primitive_type == "chararray":
            line = f'read_chararray({arg}, commondata->{key}, "{key}", {c_param.buffer_size});'
        elif c_param.primitive_type == "integer":
            line = f'read_integer({arg}, &commondata->{key}, "{key}");'
        elif c_param.primitive_type == "REAL":
            line = f'read_REAL({arg}, &commondata->{key}, "{key}");'
        elif c_param.primitive_type == "boolean":
            line = f'read_boolean({arg}, &commondata->{key}, "{key}");'
        steerable_c_code += f"        {line}\n"

    body = rf"""
  const int number_of_steerable_parameters = {len(cmdline_inputs)};

  int option;

  // Check for "-h" or "--help" options
  if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {{
    print_usage();
    exit(0);
  }}  // END IF: Checking for help option.

  // Determine the usage option based on argc
  if (argc == 1) {{
    option = 1; // Usage option 1: Process default parameter file "{project_name}.par"
  }} else if (argc == 2) {{
    // Check if the argument is a file
    FILE *file_check = fopen(argv[1], "r");
    if (file_check != NULL) {{
      fclose(file_check);
      option = 2; // Usage option 2: Process parameter file provided in argv[1]
    }} else if (argc == 1 + number_of_steerable_parameters) {{
      option = 3;
    }} else {{
      fprintf(stderr, "Error: Invalid number of arguments or file cannot be opened.\n");
      print_usage();
      exit(1);
    }}
  }} else if (argc == 1 + number_of_steerable_parameters) {{
    option = 3; // Usage option 3: Overwrite steerable parameters after processing "{project_name}.par"
  }} else if (argc == 2 + number_of_steerable_parameters) {{
    // Check if the first argument is a file
    FILE *file_check = fopen(argv[1], "r");
    if (file_check != NULL) {{
      fclose(file_check);
      option = 4; // Usage option 4: Overwrite steerable parameters after processing parameter file provided in argv[1]
    }} else {{
      fprintf(stderr, "Error: File cannot be opened for option 4.\n");
      print_usage();
      exit(1);
    }}
  }} else {{
    fprintf(stderr, "Error: Invalid number of arguments\n");
    print_usage();
    exit(1);
  }}  // END IF-ELSE: Determining usage option.

  // fprintf(stderr, "Using option %d\n", option);

  // Determine the filename based on the selected option.
  const char *filename = (option == 1 || option == 3) ? "{project_name}.par" : argv[1];
  FILE *file = fopen(filename, "r");
  if (file == NULL) {{
    print_usage();
    exit(1);
  }}  // END IF: Checking if the parameter file opened successfully.

  char line[LINE_SIZE];
  char param[PARAM_SIZE];
  char value[PARAM_SIZE];
  int params_set[NUM_PARAMETERS] = {{0}}; // Track which parameters have been set.

  // Iterate through each line of the parameter file.
  while (fgets(line, sizeof(line), file)) {{
    // Remove comments from the line.
    char *comment_start = strchr(line, '#');
    if (comment_start != NULL) {{
      *comment_start = '\0';
    }}

    // Split the line into parameter and value based on '=' delimiter.
    char *p = strtok(line, "=");
    if (p) {{
      safe_copy(param, trim_space(p), sizeof(param));
      p = strtok(NULL, "=");
      if (p) {{
        safe_copy(value, trim_space(p), sizeof(value));

        char param_name[PARAM_SIZE];
        int array_size = 0;
        parse_param(param, param_name, &array_size);

        // Validate characters in the parameter name.
        for (int i = 0; param_name[i]; i++) {{
          if (!isalnum(param_name[i]) && param_name[i] != '_') {{
            fprintf(stderr, "Error: Invalid parameter name %s.\n", param_name);
            exit(1);
          }}
        }}

        // Parse the value string into individual values.
        char values_array[MAX_ARRAY_SIZE][PARAM_SIZE];
        int value_count = 0;
        parse_value(value, values_array, &value_count);

        // Lookup the parameter descriptor.
        param_descriptor *param_desc = find_param_descriptor(param_name);
        if (param_desc == NULL) {{
          fprintf(stderr, "Warning: Unrecognized parameter %s.\n", param_name);
          continue; // Decide whether to exit or ignore
        }}

        // Check for duplicate parameter definitions.
        if (params_set[param_desc->index] == 1) {{
          fprintf(stderr, "Error: Duplicate parameter %s.\n", param_name);
          exit(1);
        }}
        params_set[param_desc->index] = 1;

        // Validate array size if applicable.
        if (param_desc->type != PARAM_CHARARRAY && param_desc->array_size > 0) {{
          // It's an array parameter
          if (array_size != param_desc->array_size) {{
            fprintf(stderr, "Error: Array size mismatch for parameter %s.\n", param_name);
            exit(1);
          }}
          if (value_count != param_desc->array_size) {{
            fprintf(stderr, "Error: Number of values does not match array size for parameter %s.\n", param_name);
            exit(1);
          }}
        }} else {{
          // It's a scalar parameter, including PARAM_CHARARRAY
          if (array_size > 0) {{
            fprintf(stderr, "Error: Unexpected array size for scalar parameter %s.\n", param_name);
            exit(1);
          }}
          if (value_count != 1) {{
            fprintf(stderr, "Error: Expected a single value for parameter %s.\n", param_name);
            exit(1);
          }}
        }}

        // Assign the parsed values to the corresponding fields in commondata.
{assignment_code}
        else {{
            fprintf(stderr, "Error: Unknown parameter index for %s.\n", param_name);
            exit(1);
        }}
    }}
}}

    // Handling options 3 and 4: Overwriting steerable parameters
    if (option == 3 || option == 4) {{
        // For options 3 and 4, we extract the last arguments as steerable parameters
{steerable_c_code.rstrip()}
    }} // END IF (option == 3 || option == 4)
  }} // END WHILE LOOP over all lines in the file
  fclose(file); // Close the parameter file.
"""
    return body


def register_CFunction_cmdline_input_and_parfile_parser(
    project_name: str,
    cmdline_inputs: Optional[List[str]] = None,
) -> None:
    """
    Register a C function to handle command-line input and parse parameter files for a given project.

    This function defines and registers a series of C functions and constants that facilitate reading
    and parsing parameter files containing key-value pairs. It supports handling whitespace and comments,
    and has specific error handling for buffer overflows, invalid integers, and invalid doubles.

    It also incorporates various usage options for handling command-line arguments and integrates
    steerable parameters that can be overwritten from the command line.

    :param project_name: Name of the project. Used for file naming and error messaging.
    :param cmdline_inputs: Optional list of command-line inputs that can be used to overwrite specific
                           parameters defined in the parameter file.
    """
    if cmdline_inputs is None:
        cmdline_inputs = []

    # Step 1: Analyze parameters and determine C generation requirements.
    parameters_list, needed_readers = _process_parameters_for_c_parser()

    # Step 2: Generate the C preamble with all helper functions.
    prefunc = _generate_c_preamble(
        project_name, cmdline_inputs, parameters_list, needed_readers
    )

    # Step 3: Generate the main C function body.
    body = _generate_c_body(project_name, cmdline_inputs, parameters_list)

    # Step 4: Register the complete C function with NRPy+.
    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "<string.h>",
            "<ctype.h>",
            "<errno.h>",
            "<stdbool.h>",
        ],
        prefunc=prefunc,
        desc=r"""AUTOMATICALLY GENERATED BY cmdline_input_and_parfiles.py

Reads and parses a parameter file to populate the commondata_struct.

This function processes command-line arguments and reads parameters from a specified
parameter file or a default file. It supports various usage options, including displaying
help information, reading different parameter files, and overwriting steerable parameters
with provided convergence factors. The function ensures that all parameters are valid,
correctly formatted, and not duplicated.

@param commondata - Pointer to the commondata_struct to be populated.
@param argc - The argument count from the command-line input.
@param argv - The argument vector containing command-line arguments.""",
        cfunc_type="void",
        name="cmdline_input_and_parfile_parser",
        params="commondata_struct *restrict commondata, int argc, const char *argv[]",
        body=body,
    )


# endregion

# region Default Parameter File Generation


def _parse_c_array_type_string(cparam_type: str) -> Optional[Dict[str, Any]]:
    """
    Parse a C type string (e.g., "REAL[3]") to extract array info.

    :param cparam_type: The C type string to parse.
    :return: A dictionary with 'base_type' and 'size' if it's an array, otherwise None.
    """
    if "[" in cparam_type and "]" in cparam_type:
        base_type, rest = cparam_type.split("[", 1)
        size_str = rest.split("]", 1)[0]
        if size_str.isdigit():
            return {"base_type": base_type.strip(), "size": int(size_str)}
    return None


def _format_parfile_line(parname: str, code_param: par.CodeParameter) -> str:
    """
    Format a single parameter entry for the parameter file.

    :param parname: The name of the parameter.
    :param code_param: The CodeParameter object containing its details.
    :raises ValueError: For unsupported or misconfigured parameter types.
    :return: A formatted string for one line of the parameter file.
    """
    cp_type = code_param.cparam_type
    description = code_param.description.strip()
    desc_suffix = f" {description}" if description else ""
    default_val = code_param.defaultvalue

    array_info = _parse_c_array_type_string(cp_type)

    if array_info:
        base_type = array_info["base_type"]
        size = array_info["size"]
        if base_type.lower() not in ("real", "int", "char"):
            raise ValueError(
                f"Unsupported array base type '{base_type}' for parameter '{parname}'. Only 'char', 'int', and 'REAL' are supported."
            )

        if base_type.lower() in ("real", "int"):
            display_type = "REAL" if base_type.lower() == "real" else "int"
            if isinstance(default_val, list):
                defaults = ", ".join(map(str, default_val))
            else:
                cast_val = (
                    float(default_val)
                    if base_type.lower() == "real"
                    else int(default_val)
                )
                defaults = ", ".join(str(cast_val) for _ in range(size))
            return (
                f"{parname}[{size}] = {{ {defaults} }}  # ({display_type}){desc_suffix}"
            )

        # char array
        if not isinstance(default_val, str):
            raise ValueError(
                f"Default value for char array parameter '{parname}' must be a string."
            )
        escaped_val = default_val.replace('"', '\\"')
        return f'{parname} = "{escaped_val}"  # (char[{size}]){desc_suffix}'

    # Scalar parameter
    if cp_type.lower().startswith("char"):
        if not isinstance(default_val, str):
            raise ValueError(
                f"Default value for char parameter '{parname}' must be a string."
            )
        escaped_val = default_val.replace('"', '\\"')
        return f'{parname} = "{escaped_val}"  # ({cp_type}){desc_suffix}'

    return f"{parname} = {default_val}  # ({cp_type}){desc_suffix}"


def _align_parfile_comments(text_block: str) -> str:
    """
    Align the '#' comment symbols in a block of text for neatness.

    :param text_block: A multi-line string.
    :return: The same block of text with comments aligned.
    """
    lines = text_block.split("\n")
    max_hash_pos = 0
    for line in lines:
        if "#" in line and not line.strip().startswith("#"):
            max_hash_pos = max(max_hash_pos, line.find("#"))

    adjusted_lines = []
    for line in lines:
        if "#" in line and not line.strip().startswith("#"):
            hash_pos = line.find("#")
            spaces_to_add = max_hash_pos - hash_pos
            adjusted_lines.append(
                line[:hash_pos] + " " * spaces_to_add + line[hash_pos:]
            )
        else:
            adjusted_lines.append(line)
    return "\n".join(adjusted_lines)


def generate_default_parfile(project_dir: str, project_name: str) -> None:
    """
    Generate a default parameter file with sorted modules and parameters.

    :param project_dir: The directory where the parameter file will be stored.
    :param project_name: The name of the project.

    Doctest:
    >>> import shutil
    >>> from pathlib import Path
    >>> # Clear existing parameters
    >>> par.glb_code_params_dict.clear()
    >>> # Register scalar REAL parameters
    >>> _, __ = par.register_CodeParameters(
    ...     "REAL", "CodeParameters_c_files",
    ...     ["a", "pi_three_sigfigs"],
    ...     [1.0, 3.14],
    ...     commondata=True,
    ...     descriptions=["The value of a", "Pi to three significant figures"]
    ... )
    >>> # Register a #define parameter
    >>> ___ = par.register_CodeParameter(
    ...     "#define", "CodeParameters_c_files", "b", 0, description="A define parameter"
    ... )
    >>> # Register parameters that should be excluded from the parfile
    >>> _leaveitbe = par.register_CodeParameter(
    ...     "REAL", "CodeParameters_c_files", "leaveitbe",
    ...     add_to_parfile=False, add_to_set_CodeParameters_h=False
    ... )
    >>> _leaveitoutofparfile = par.register_CodeParameter(
    ...     "REAL", "CodeParameters_c_files", "leaveitoutofparfile",
    ...     add_to_parfile=False
    ... )
    >>> # Register a char array parameter
    >>> _str = par.register_CodeParameter(
    ...     "char[100]", "CodeParameters_c_files", "string",
    ...     "cheese", commondata=True, description="A string parameter"
    ... )
    >>> # Register another char array parameter
    >>> _str2 = par.register_CodeParameter(
    ...     "char[50]", "CodeParameters_c_files", "outer_bc_type",
    ...     "radiation", commondata=True, description="A bc string parameter"
    ... )
    >>> # Register an int parameter
    >>> _int = par.register_CodeParameter(
    ...     "int", "CodeParameters_c_files", "blahint", -1,
    ...     commondata=True, add_to_parfile=True, add_to_set_CodeParameters_h=False,
    ...     description="An integer parameter"
    ... )
    >>> # Register a bool parameter
    >>> _bool = par.register_CodeParameter(
    ...     "bool", "CodeParameters_c_files", "BHaH_is_amazing",
    ...     "true", description="A boolean parameter"
    ... )
    >>> # Register a REAL array parameter
    >>> _real_array = par.register_CodeParameter(
    ...     cparam_type="REAL[3]", module="CodeParameters_c_files",
    ...     name="bah_initial_x_center", defaultvalue=0.0,
    ...     commondata=True, add_to_parfile=True, add_to_set_CodeParameters_h=False,
    ...     description="Initial X centers"
    ... )
    >>> # Register an int array parameter
    >>> _int_array = par.register_CodeParameter(
    ...     cparam_type="int[2]", module="CodeParameters_c_files",
    ...     name="initial_levels", defaultvalue=[4, 2],
    ...     commondata=True, add_to_parfile=True, add_to_set_CodeParameters_h=False,
    ...     description=""
    ... )

    >>> # Clear any existing CFunction_dict if necessary
    >>> cfc.CFunction_dict.clear()
    >>> # Setup project directory
    >>> project_dir = Path("/tmp/tmp_BHaH_parfile")
    >>> if project_dir.exists():
    ...     shutil.rmtree(project_dir)
    >>> project_dir.mkdir(parents=True, exist_ok=True)
    >>> # Register an unsupported array type parameter
    >>> _unsupported_array = par.register_CodeParameter(
    ...     cparam_type="double[5]", module="CodeParameters_c_files",
    ...     name="unsupported_param", defaultvalue=1.0,
    ...     commondata=True, add_to_parfile=True, add_to_set_CodeParameters_h=False,
    ...     description="Unsupported array type"
    ... )
    >>> # Generate the parfile
    >>> generate_default_parfile(str(project_dir), "example_project")
    Traceback (most recent call last):
    ...
    ValueError: Unsupported array base type 'double' for parameter 'unsupported_param'. Only 'char', 'int', and 'REAL' are supported.
    >>> # Now remove the unsupported parameter and generate the parfile successfully
    >>> del par.glb_code_params_dict["unsupported_param"]
    >>> # Generate the parfile
    >>> generate_default_parfile(str(project_dir), "example_project")
    >>> # Read and print the generated parfile
    >>> print((project_dir / 'example_project.par').read_text())
    #### example_project BH@H parameter file. NOTE: only commondata CodeParameters appear here ###
    ###########################
    ###########################
    ### Module: CodeParameters_c_files
    a = 1.0                                      # (REAL) The value of a
    bah_initial_x_center[3] = { 0.0, 0.0, 0.0 }  # (REAL) Initial X centers
    blahint = -1                                 # (int) An integer parameter
    initial_levels[2] = { 4, 2 }                 # (int)
    outer_bc_type = "radiation"                  # (char[50]) A bc string parameter
    pi_three_sigfigs = 3.14                      # (REAL) Pi to three significant figures
    string = "cheese"                            # (char[100]) A string parameter
    <BLANKLINE>
    """
    parfile_output_by_module: Dict[str, List[str]] = defaultdict(list)
    sorted_params = sorted(
        par.glb_code_params_dict.items(), key=lambda item: (item[1].module, item[0])
    )

    for parname, code_param in sorted_params:
        if code_param.commondata and code_param.add_to_parfile:
            line = _format_parfile_line(parname, code_param)
            parfile_output_by_module[code_param.module].append(line)

    output_parts = [
        f"#### {project_name} BH@H parameter file. NOTE: only commondata CodeParameters appear here ###"
    ]
    for module, params in sorted(parfile_output_by_module.items()):
        output_parts.append("###########################")
        output_parts.append("###########################")
        output_parts.append(f"### Module: {module}")
        output_parts.extend(params)

    # Add a blank line at the end to match doctest
    output_parts.append("")

    final_output = _align_parfile_comments("\n".join(output_parts))
    with (Path(project_dir) / f"{project_name}.par").open(
        "w", encoding="utf-8"
    ) as file:
        file.write(final_output)


# endregion

if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
