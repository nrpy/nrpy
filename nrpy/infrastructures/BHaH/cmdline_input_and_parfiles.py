"""
Functions for parsing command-line input and parameter files.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import nrpy.c_function as cfc
import nrpy.params as par


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

    # Count the number of parameters to be included
    num_commondata_params = sum(
        1
        for CodeParam in par.glb_code_params_dict.values()
        if CodeParam.commondata and CodeParam.add_to_parfile
    )

    # Initialize the preamble for the C code
    prefunc = f"#define NUM_PARAMETERS {num_commondata_params} // Define the number of parameters\n\n"
    prefunc += r"""
#define LINE_SIZE 1024 // Define the max size of a line
#define PARAM_SIZE 128 // Define the max param string size

// Trim leading and trailing spaces
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

// Safe copy function to prevent buffer overflows
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

    # Generate the usage instructions string
    list_of_steerable_params_str = " ".join(cmdline_inputs)
    prefunc += rf"""
// Function to print usage instructions
static void print_usage() {{
  fprintf(stderr, "Usage option 0: ./{project_name} [--help or -h] - Outputs this usage command\n");
  fprintf(stderr, "Usage option 1: ./{project_name} - Reads in parameter file {project_name}.par\n");
  fprintf(stderr, "Usage option 2: ./{project_name} [parfile] - Reads in parameter file [parfile]\n");
  fprintf(stderr, "Usage option 3: ./{project_name} [{list_of_steerable_params_str}] - Overwrites parameters in list after reading in {project_name}.par\n");
  fprintf(stderr, "Usage option 4: ./{project_name} [parfile] [{list_of_steerable_params_str}] - Overwrites list of steerable parameters after reading in [parfile]\n");
}}"""

    # Define parameter types and descriptor struct
    prefunc += """
#define MAX_ARRAY_SIZE 100 // Adjust as needed

// Parameter types
typedef enum { PARAM_REAL, PARAM_INT, PARAM_CHARARRAY, PARAM_REAL_ARRAY, PARAM_INT_ARRAY } param_type;

// Parameter descriptor struct
typedef struct {
    const char *name;
    int index;
    param_type type;
    int array_size; // For array parameters
} param_descriptor;

// param_table[] is a list of parameter descriptors, each containing the parameter's name, unique index, type, and array size.
// This table is used to lookup and manage parameters during parsing, ensuring that each parameter is correctly identified
// and assigned the appropriate value in the common data structure.
"""

    # Define the parameter table entries
    parameters_list: List[Dict[str, Any]] = []
    param_table_entries = []
    param_index = 0
    found_integer = False
    found_REAL = False
    found_chararray = False
    found_boolean = False

    for key in sorted(par.glb_code_params_dict.keys()):
        CodeParam = par.glb_code_params_dict[key]
        if CodeParam.add_to_parfile and CodeParam.commondata:
            param_name = key
            cparam_type = CodeParam.cparam_type.strip()
            array_size = 0
            if "int[" in cparam_type or "int [" in cparam_type:
                param_type = "PARAM_INT_ARRAY"
                array_size = int(cparam_type.split("[")[1].split("]")[0])
                found_integer = True
            elif "REAL[" in cparam_type or "REAL [" in cparam_type:
                param_type = "PARAM_REAL_ARRAY"
                array_size = int(cparam_type.split("[")[1].split("]")[0])
                found_REAL = True
            elif cparam_type == "int":
                param_type = "PARAM_INT"
                found_integer = True
            elif cparam_type == "REAL":
                param_type = "PARAM_REAL"
                found_REAL = True
            elif "char" in cparam_type and "[" in cparam_type and "]" in cparam_type:
                param_type = "PARAM_CHARARRAY"
                array_size = int(cparam_type.split("[")[1].split("]")[0])
                found_chararray = True
            else:
                continue  # Skip unsupported types

            # Append to parameters_list
            parameters_list.append(
                {
                    "index": param_index,
                    "name": param_name,
                    "type": param_type,
                    "array_size": array_size,
                }
            )
            # Append to param_table_entries
            param_table_entries.append(
                f'    {{"{param_name}", {param_index}, {param_type}, {array_size}}}'
            )
            param_index += 1

    # Generate the param_table and NUM_PARAMS
    param_table_str = (
        "param_descriptor param_table[] = {\n"
        + ",\n".join(param_table_entries)
        + "\n};\n"
    )
    param_table_str += (
        "#define NUM_PARAMS (sizeof(param_table) / sizeof(param_descriptor))\n\n"
    )

    prefunc += param_table_str

    # Define the functions to find parameter descriptors and parse parameters
    prefunc += r"""
// Function to find parameter descriptor by name
param_descriptor *find_param_descriptor(const char *param_name) {
  for (int i = 0; i < NUM_PARAMS; i++) {
    if (strcmp(param_table[i].name, param_name) == 0) {
      return &param_table[i];
    }
  }
  return NULL;
}

// Function to parse parameter name and array size
static void parse_param(const char *param_str, char *param_name, int *array_size) {
  char *bracket_start = strchr(param_str, '[');
  if (bracket_start != NULL) {
    // It's an array parameter
    size_t name_len = bracket_start - param_str;
    strncpy(param_name, param_str, name_len);
    param_name[name_len] = '\0';
    char *bracket_end = strchr(bracket_start + 1, ']');
    if (bracket_end == NULL) {
      fprintf(stderr, "Error: Missing closing bracket in parameter %s.\n", param_str);
      exit(1);
    }
    char size_str[16];
    size_t size_len = bracket_end - bracket_start - 1;
    strncpy(size_str, bracket_start + 1, size_len);
    size_str[size_len] = '\0';
    *array_size = atoi(size_str);
    if (*array_size <= 0) {
      fprintf(stderr, "Error: Invalid array size in parameter %s.\n", param_name);
      exit(1);
    }
  } else {
    // Scalar parameter
    safe_copy(param_name, param_str, PARAM_SIZE);
    *array_size = 0;
  }
}

// Function to parse value string into array of values
static void parse_value(const char *value_str, char values[][PARAM_SIZE], int *value_count) {
  if (value_str[0] == '{') {
    // Array value
    size_t len = strlen(value_str);
    if (value_str[len - 1] != '}') {
      fprintf(stderr, "Error: Missing closing brace in value %s.\n", value_str);
      exit(1);
    }
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
      }
      val_token = strtok_r(NULL, ",", &saveptr);
    }
    *value_count = count;
  } else {
    // Scalar value
    char mutable_value[PARAM_SIZE];
    safe_copy(mutable_value, value_str, PARAM_SIZE);
    char *trimmed = trim_space(mutable_value);
    safe_copy(values[0], trimmed, PARAM_SIZE);
    *value_count = 1;
  }
}
"""

    # Add reading functions if needed
    if found_integer:
        prefunc += r"""
// Function to read an integer value
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
"""

    if found_REAL:
        prefunc += r"""
// Function to read a REAL (double) value
static void read_REAL(const char *value, REAL *result, const char *param_name) {
  char *endptr;
  errno = 0; // To detect overflow
  double double_val = strtod(value, &endptr);

  if (endptr == value || *endptr != '\0' || errno == ERANGE) {
    fprintf(stderr, "Error: Invalid double value for %s: %s.\n", param_name, value);
    exit(1);
  }

  *result = (REAL)double_val;
}
"""

    if found_chararray:
        prefunc += r"""
// Function to read a character array
static void read_chararray(const char *value, char *result, const char *param_name, size_t size) {
  if (strlen(value) >= size) {
    fprintf(stderr, "Error: Buffer overflow detected for %s.\n", param_name);
    exit(1);
  }
  safe_copy(result, value, size);
}
"""

    if found_boolean:
        prefunc += r"""
// Function to read a boolean value
static void read_boolean(const char *value, bool *result, const char *param_name) {
  // To allow case-insensitive comparison
  char *lower_value = strdup(value);
  if (lower_value == NULL) {
    fprintf(stderr, "Error: Memory allocation failed for boolean value of %s.\n", param_name);
    exit(1);
  }
  for (char *p = lower_value; *p != '\0'; p++) {
    *p = tolower((unsigned char)*p);
  }

  // Check if the input is "true", "false", "0", or "1"
  if (strcmp(lower_value, "true") == 0 || strcmp(lower_value, "1") == 0) {
    *result = true;
  } else if (strcmp(lower_value, "false") == 0 || strcmp(lower_value, "0") == 0) {
    *result = false;
  } else {
    fprintf(stderr, "Error: Invalid boolean value for %s: %s.\n", param_name, value);
    free(lower_value);
    exit(1);
  }

  // Free the allocated memory for the lowercase copy of the value
  free(lower_value);
}
"""

    # Start building the function body
    cfunc_type = "void"
    name = "cmdline_input_and_parfile_parser"
    params_str = "commondata_struct *restrict commondata, int argc, const char *argv[]"

    # Generate the main function body
    body = rf"""
  const int number_of_steerable_parameters = {len(cmdline_inputs)};

  int option;

  // Check for "-h" or "--help" options
  if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {{
    print_usage();
    exit(0);
  }}

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
  }}

  // fprintf(stderr, "Using option %d\n", option);

  const char *filename = (option == 1 || option == 3) ? "{project_name}.par" : argv[1];
  FILE *file = fopen(filename, "r");
  if (file == NULL) {{
    print_usage();
    exit(1);
  }}

  char line[LINE_SIZE];
  char param[PARAM_SIZE];
  char value[PARAM_SIZE];
  int params_set[NUM_PARAMETERS] = {{0}}; // Record of parameters set

  while (fgets(line, sizeof(line), file)) {{
    // Removing comments from the line
    char *comment_start = strchr(line, '#');
    if (comment_start != NULL) {{
      *comment_start = '\0';
    }}

    char *p = strtok(line, "=");
    if (p) {{
      safe_copy(param, trim_space(p), sizeof(param));
      p = strtok(NULL, "=");
      if (p) {{
        safe_copy(value, trim_space(p), sizeof(value));

        char param_name[PARAM_SIZE];
        int array_size = 0;
        parse_param(param, param_name, &array_size);

        // Check for invalid characters in parameter name
        for (int i = 0; param_name[i]; i++) {{
          if (!isalnum(param_name[i]) && param_name[i] != '_') {{
            fprintf(stderr, "Error: Invalid parameter name %s.\n", param_name);
            exit(1);
          }}
        }}

        // Parse value
        char values_array[MAX_ARRAY_SIZE][PARAM_SIZE];
        int value_count = 0;
        parse_value(value, values_array, &value_count);

        // Lookup parameter descriptor
        param_descriptor *param_desc = find_param_descriptor(param_name);
        if (param_desc == NULL) {{
          fprintf(stderr, "Warning: Unrecognized parameter %s.\n", param_name);
          continue; // Decide whether to exit or ignore
        }}

        // Check for duplicate parameters
        if (params_set[param_desc->index] == 1) {{
          fprintf(stderr, "Error: Duplicate parameter %s.\n", param_name);
          exit(1);
        }}
        params_set[param_desc->index] = 1;

        // Check array size
        if (param_desc->array_size > 0) {{
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
          // It's a scalar parameter
          if (array_size > 0) {{
            fprintf(stderr, "Error: Unexpected array size for scalar parameter %s.\n", param_name);
            exit(1);
          }}
          if (value_count != 1) {{
            fprintf(stderr, "Error: Expected a single value for parameter %s.\n", param_name);
            exit(1);
          }}
        }}

        // Now, assign values
"""

    # Build assignment code
    assignment_code = ""
    first_param = True
    for param in parameters_list:
        index = param["index"]
        param_name = param["name"]
        param_type = param["type"]
        array_size = param["array_size"]
        if first_param:
            prefix = "if"
            first_param = False
        else:
            prefix = "else if"
        if param_type == "PARAM_REAL":
            assignment_code += (
                f"                {prefix} (param_desc->index == {index}) {{\n"
                f'                    read_REAL(values_array[0], &commondata->{param_name}, "{param_name}");\n'
                f"                }}\n"
            )
        elif param_type == "PARAM_INT":
            assignment_code += (
                f"                {prefix} (param_desc->index == {index}) {{\n"
                f'                    read_integer(values_array[0], &commondata->{param_name}, "{param_name}");\n'
                f"                }}\n"
            )
        elif param_type == "PARAM_CHARARRAY":
            # Assuming size from array_size
            size = array_size
            assignment_code += (
                f"                {prefix} (param_desc->index == {index}) {{\n"
                f'                    read_chararray(values_array[0], commondata->{param_name}, "{param_name}", {size});\n'
                f"                }}\n"
            )
        elif param_type == "PARAM_REAL_ARRAY":
            assignment_code += (
                f"                {prefix} (param_desc->index == {index}) {{\n"
                f"                    for (int i = 0; i < {array_size}; i++) {{\n"
                f'                        read_REAL(values_array[i], &commondata->{param_name}[i], "{param_name}");\n'
                f"                    }}\n"
                f"                }}\n"
            )
        elif param_type == "PARAM_INT_ARRAY":
            assignment_code += (
                f"                {prefix} (param_desc->index == {index}) {{\n"
                f"                    for (int i = 0; i < {array_size}; i++) {{\n"
                f'                        read_integer(values_array[i], &commondata->{param_name}[i], "{param_name}");\n'
                f"                    }}\n"
                f"                }}\n"
            )

    # Add the assignment code and the else block
    body += assignment_code
    body += r"""                else {
                    fprintf(stderr, "Error: Unknown parameter index for %s.\n", param_name);
                    exit(1);
                }
            }
        }

        fclose(file);
        // Handling options 3 and 4: Overwriting steerable parameters
        if (option == 3 || option == 4) {
            // For options 3 and 4, we extract the last arguments as steerable parameters
"""

    # Handle steerable parameters overwriting
    steerable_body = ""
    for i, key in enumerate(cmdline_inputs):
        CodeParam = par.glb_code_params_dict[key]
        if CodeParam.add_to_parfile and CodeParam.commondata:
            if CodeParam.cparam_type == "int":
                steerable_body += f'            read_integer(argv[argc - number_of_steerable_parameters + {i}], &commondata->{key}, "{key}");\n'
            elif CodeParam.cparam_type == "REAL":
                steerable_body += f'            read_REAL(argv[argc - number_of_steerable_parameters + {i}], &commondata->{key}, "{key}");\n'
    body += steerable_body
    body += "        }\n"  # Close the if (option == 3 || option == 4) block
    body += "}\n"  # Close the main function

    # Register the C function with the constructed components
    cfc.register_CFunction(
        includes=[
            "BHaH_defines.h",
            "<string.h>",
            "<ctype.h>",
            "<errno.h>",
            "<stdbool.h>",
        ],
        prefunc=prefunc,
        desc=r"""AUTOMATICALLY GENERATED BY parameter_file_read_and_parse.py
parameter_file_read_and_parse() function:
Reads and parses a parameter file to populate commondata_struct commondata.

This function takes in the command-line arguments and a pointer to a commondata_struct.
It reads the provided file and extracts the parameters defined in the file, populating
the commondata_struct with the values. The file is expected to contain key-value pairs
separated by an equals sign (=), and it may include comments starting with a hash (#).
The function handles errors such as file opening failure, duplicate parameters, and
invalid parameter names.

@param commondata: Pointer to the commondata struct to be populated.
@param argc: The argument count from the command-line input.
@param argv: The argument vector containing command-line arguments.""",
        cfunc_type=cfunc_type,
        name=name,
        params=params_str,
        body=body,
    )


def generate_default_parfile(project_dir: str, project_name: str) -> None:
    """
    Generate a default parameter file with sorted modules and parameters.

    :param project_dir: The parameter file will be stored in project_dir.
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
    ...     "char", "CodeParameters_c_files", "string[100]",
    ...     "cheese", commondata=True, description="A string parameter"
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
    ...     name="initial_levels", defaultvalue=4,
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
    >>> # Generate the parfile
    >>> generate_default_parfile(str(project_dir), "example_project")
    >>> # Read and print the generated parfile
    >>> print((project_dir / 'example_project.par').read_text())
    #### example_project BH@H parameter file. NOTE: only commondata CodeParameters appear here ###
    ###########################
    ###########################
    ### Module: CodeParameters_c_files
    a = 1.0                                      # (REAL) The value of a
    bah_initial_x_center[3] = { 0.0, 0.0, 0.0 }  # (REAL[3]) Initial X centers
    blahint = -1                                 # (int) An integer parameter
    initial_levels[2] = { 4, 4 }                 # (int[2])
    pi_three_sigfigs = 3.14                      # (REAL) Pi to three significant figures
    string[100] = cheese                         # (char) A string parameter
    <BLANKLINE>
    """
    parfile_output_dict: Dict[str, List[str]] = defaultdict(list)

    # Function to parse array types without using regex
    def parse_array_type(cparam_type: str) -> Union[None, Dict[str, Any]]:
        if "[" in cparam_type and "]" in cparam_type:
            base_type = cparam_type[: cparam_type.find("[")]
            size_str = cparam_type[cparam_type.find("[") + 1 : cparam_type.find("]")]
            if size_str.isdigit():
                size = int(size_str)
                return {"base_type": base_type, "size": size}
        return None

    # Sorting by module name
    for parname, CodeParam in sorted(
        par.glb_code_params_dict.items(), key=lambda x: x[1].module
    ):
        if CodeParam.commondata and CodeParam.add_to_parfile:
            CPtype = CodeParam.cparam_type
            array_info = parse_array_type(CPtype)
            description = CodeParam.description.strip()
            description_suffix = f" {description}" if description else ""
            if array_info:
                base_type = array_info["base_type"]
                size = array_info["size"]
                default_val = CodeParam.defaultvalue
                # Generate a list with the default value repeated 'size' times
                if base_type.lower() in ["real", "int"]:
                    # Format based on base type
                    if base_type.lower() == "real":
                        default_vals = ", ".join([f"{float(default_val)}"] * size)
                    else:  # int
                        default_vals = ", ".join([str(int(default_val))] * size)
                    # **Append the size to the variable name**
                    parfile_output_dict[CodeParam.module].append(
                        f"{parname}[{size}] = {{ {default_vals} }}  # ({CPtype}){description_suffix}\n"
                    )
                else:
                    # Handle other array types if necessary
                    parfile_output_dict[CodeParam.module].append(
                        f"{parname} = {{ {default_val} }}  # ({CPtype}){description_suffix}\n"
                    )
            else:
                if CPtype.startswith("char") and "[" in CPtype and "]" in CPtype:
                    # Handle char arrays
                    chararray_name = parname.split("[")[0]
                    parfile_output_dict[CodeParam.module].append(
                        f'{chararray_name} = "{CodeParam.defaultvalue}"  # (char array){description_suffix}\n'
                    )
                else:
                    # Handle scalar parameters
                    parfile_output_dict[CodeParam.module].append(
                        f"{parname} = {CodeParam.defaultvalue}  # ({CPtype}){description_suffix}\n"
                    )

    # Sorting the parameters within each module
    for module in parfile_output_dict:
        parfile_output_dict[module] = sorted(parfile_output_dict[module])

    output_str = f"#### {project_name} BH@H parameter file. NOTE: only commondata CodeParameters appear here ###\n"
    for module, params in sorted(parfile_output_dict.items()):
        output_str += "###########################\n"
        output_str += "###########################\n"
        output_str += f"### Module: {module}\n"
        output_str += "".join(params)

    def align_by_hash(original_string: str) -> str:
        lines = original_string.split("\n")
        max_length = max(
            (
                line.find("#")
                for line in lines
                if "#" in line and not line.strip().startswith("#")
            ),
            default=0,
        )

        adjusted_lines = []
        for line in lines:
            if "#" in line and not line.strip().startswith("#"):
                index = line.find("#")
                spaces_needed = max_length - index
                if spaces_needed > 0:
                    adjusted_line = line[:index] + " " * spaces_needed + line[index:]
                else:
                    adjusted_line = line
                adjusted_lines.append(adjusted_line)
            else:
                adjusted_lines.append(line)

        return "\n".join(adjusted_lines)

    with (Path(project_dir) / f"{project_name}.par").open(
        "w", encoding="utf-8"
    ) as file:
        file.write(align_by_hash(output_str))


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
