"""
Constructs a Python project structure for JAX-based functions from registered PyFunctions.

This module generates a complete Python package structure from functions registered in PyFunction_dict,
including proper package layout, dependencies, and documentation. It's specifically designed for
JAX-accelerated numerical computing projects.

Example:
    >>> from nrpy.infrastructures.JAX.jax_project_generator import output_PyFunction_files_and_construct_project
    >>> output_PyFunction_files_and_construct_project("project/blank_project", "blank_project")

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import logging
from pathlib import Path
from typing import Dict, List, Set, Union

from nrpy.infrastructures.JAX.commondata import (
    commondata_params_dict,
    generate_commondata_dataclass,
    register_commondata_param,
)
from nrpy.py_function import PyFunction_dict, register_PyFunction

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)

# Define constants for filenames to avoid magic strings
_JAX_PROTOTYPES_H = "JAX_function_prototypes.h"
_PYTHON_INIT = "__init__.py"
_REQUIREMENTS = ["jax", "jaxlib>=0.4.0", "numpy>=1.21.0"]
PYTHON_MIN_VERSION = "3.8"


def _validate_inputs(
    project_dir: Union[str, Path],
    project_name: str,
) -> None:
    """
    Validate the inputs for project generation.

    :param project_dir: The target directory for the project.
    :param project_name: The name of the project.
    :raises TypeError: If input types are incorrect.
    :raises ValueError: If project_name is invalid.
    :raises OSError: If project_dir cannot be created or written to.
    """
    logger.info("Validating inputs...")

    # Type validation
    if not isinstance(project_name, str):
        raise TypeError(
            f"project_name must be a string, got {type(project_name).__name__}"
        )

    # Project name validation
    if not project_name.replace("_", "").isalnum():
        raise ValueError(
            "Project name must only contain alphanumeric characters or underscores"
        )
    if project_name[0].isnumeric():
        raise ValueError("Project name cannot start with a digit")

    # Directory validation
    project_path = Path(project_dir).resolve()
    if project_path.exists() and not project_path.is_dir():
        raise OSError(f"{project_path} exists and is not a directory")

    # Check write permissions
    test_file = project_path / ".write_test"
    try:
        project_path.mkdir(parents=True, exist_ok=True)
        test_file.touch()
        test_file.unlink()
    except (OSError, PermissionError) as e:
        raise OSError(f"Cannot write to {project_path}: {e}") from e


def _generate_python_files_and_init(
    project_path: Path, lib_function_prefix: str
) -> List[Path]:
    """
    Generate Python source files and the main __init__.py file.

    :param project_path: The path to the project directory.
    :param lib_function_prefix: Prefix to add to library function names.
    :return: A list of Path objects to the generated Python source files.
    :raises Exception: If there are issues creating directories or writing files.
    """
    if not PyFunction_dict:
        logger.warning(
            "No functions found in PyFunction_dict. Only Commondata module will be generated."
        )

    src_path = project_path / "src" / project_path.name
    src_path.mkdir(parents=True, exist_ok=True)
    py_files: List[Path] = []
    created_dirs: Set[Path] = set()

    logger.info("Generating Python files in %s...", src_path)

    try:
        # Create Python code files and directory structure
        for name, pyfunc in list(
            PyFunction_dict.items()
        ):  # Create a list to avoid modifying dict during iteration
            try:
                # Apply prefix if specified
                if lib_function_prefix:
                    name = f"{lib_function_prefix}{name}"
                    pyfunc.name = name
                    pyfunc.full_function = pyfunc.generate_full_function()

                # Determine file path and create necessary directories
                subdir = Path(pyfunc.subdirectory or ".")
                module_dir = src_path / subdir

                # Skip if module_dir is None or empty
                if not module_dir:
                    logger.warning("Skipping %s: Invalid module directory", name)
                    continue

                if module_dir not in created_dirs:
                    module_dir.mkdir(parents=True, exist_ok=True)
                    created_dirs.add(module_dir)

                    # Create __init__.py in new directories
                    init_path = module_dir / _PYTHON_INIT
                    with init_path.open("w", encoding="utf-8") as f:
                        f.write("# Package initialization\n")
                    logger.debug("Created package: %s", init_path)

                # Write the Python module
                py_file_path = module_dir / f"{name}.py"
                module_docstring = f'"""{name} module."""'
                try:
                    with py_file_path.open("w", encoding="utf-8") as f:
                        f.write(f"{module_docstring}\n\n{pyfunc.full_function}\n")
                except Exception as e:
                    logger.error("Failed to write %s: %s", py_file_path, e)
                    raise
                py_files.append(py_file_path.relative_to(project_path))
                logger.debug("Generated module: %s", py_file_path)

            except Exception as e:
                logger.error(
                    "Failed to generate Python file for %s: %s", name, e, exc_info=True
                )
                raise

        # Generate Commondata.py
        commondata_content = generate_commondata_dataclass()
        commondata_file_path = src_path / "Commondata.py"
        with commondata_file_path.open("w", encoding="utf-8") as f:
            f.write(commondata_content)
        py_files.append(commondata_file_path.relative_to(project_path))
        logger.debug("Generated module: %s", commondata_file_path)

        # Generate top-level __init__.py with all exports
        init_content = [
            f'"""{project_path.name} - JAX-accelerated numerical computations."""\n\n__version__ = "0.1.0"\n\n',
            "# Import all functions to make them available at package level\n",
        ]

        # Group imports by subdirectory
        imports_by_dir: Dict[str, List[str]] = {}
        for name in sorted(PyFunction_dict.keys()):
            pyfunc = PyFunction_dict[name]
            if pyfunc.subdirectory and pyfunc.subdirectory != ".":
                dir_name = pyfunc.subdirectory
                func_name = (
                    f"{lib_function_prefix}{name}" if lib_function_prefix else name
                )
                if dir_name not in imports_by_dir:
                    imports_by_dir[dir_name] = []
                imports_by_dir[dir_name].append(func_name)

        # Add imports from submodules
        for dir_name, funcs in imports_by_dir.items():
            for func in funcs:
                init_content.append(
                    f"from .{dir_name.replace('/', '.')}.{func} import {func}"
                )
            init_content.append("")

        # Add direct imports (no subdirectory)
        direct_imports = []
        for name in sorted(PyFunction_dict.keys()):
            pyfunc = PyFunction_dict[name]
            if not pyfunc.subdirectory or pyfunc.subdirectory == ".":
                func_name = (
                    f"{lib_function_prefix}{name}" if lib_function_prefix else name
                )
                direct_imports.append(func_name)
                init_content.append(f"from .{func_name} import {func_name}")

        init_content.append("from .Commondata import Commondata")

        # Add __all__ for better IDE support
        if direct_imports:
            all_exports = direct_imports + ["Commondata"]
            init_content.append(
                "\n__all__ = [" + ", ".join(f'"{f}"' for f in all_exports) + "]"
            )

        # Write the main __init__.py
        init_path = src_path / "__init__.py"
        with init_path.open("w", encoding="utf-8") as f:
            f.write("\n".join(init_content))

        logger.info("Generated %d Python modules in %s", len(py_files), src_path)
        return py_files

    except Exception as e:
        logger.error("Failed to generate Python files: %s", e, exc_info=True)
        raise


def _generate_project_metadata(project_path: Path, project_name: str) -> None:
    """
    Generate project metadata files like README, setup.cfg, etc.

    :param project_path: Path to the project directory.
    :param project_name: Name of the project.
    :raises Exception: If there are other errors during project generation.
    """
    logger.info("Generating project metadata...")

    # Create pyproject.toml
    pyproject_content = """[build-system]
requires = ["setuptools>=42", "wheel"]
build-backend = "setuptools.build_meta"
"""
    (project_path / "pyproject.toml").write_text(pyproject_content, encoding="utf-8")

    # Create setup.cfg
    setup_cfg_content = """[metadata]
name = {name}
version = 0.1.0
author = Your Name
author_email = your.email@example.com
description = JAX-accelerated numerical computations
long_description = file: README.md
long_description_content_type = text/markdown
url = https://github.com/username/{name}
project_urls =
    Bug Tracker = https://github.com/username/{name}/issues
classifiers =
    Programming Language :: Python :: 3
    License :: OSI Approved :: MIT License
    Operating System :: OS Independent

[options]
package_dir =
    = src
packages = find:
python_requires = >=3.8
install_requires =
    jax
    jaxlib>=0.4.0
    numpy>=1.21.0

[options.packages.find]
where = src

[options.extras_require]
test =
    pytest

[egg_info]
tag_build = 
tag_date = 0
""".format(name=project_name)
    (project_path / "setup.cfg").write_text(setup_cfg_content, encoding="utf-8")

    # Create requirements.txt
    requirements_content = "\n".join(_REQUIREMENTS) + "\n"
    (project_path / "requirements.txt").write_text(
        requirements_content, encoding="utf-8"
    )

    # Create README.md
    readme_content = f"""# {project_name}

JAX-accelerated numerical computations.

## Installation

### From source
```bash
git clone https://github.com/username/{project_name}.git
cd {project_name}
pip install -e .  # For development
# or
pip install .     # For regular installation
```

### From PyPI
```bash
pip install {project_name.lower().replace("-", "")}
```

## Usage

```python
import {project_name}

# Your code here
```

## Development

### Running tests
```bash
pip install -e ".[test]"
pytest
```

### Code formatting
```bash
pip install black isort
black .
isort .
```
"""
    (project_path / "README.md").write_text(readme_content, encoding="utf-8")

    # Create .gitignore
    gitignore_content = """# Byte-compiled / optimized / DLL files
__pycache__/
*.py[cod]
*$py.class

# C extensions
*.so

# Distribution / packaging
.Python
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
share/python-wheels/
*.egg-info/
.installed.cfg
*.egg
MANIFEST

# PyInstaller
*.manifest
*.spec

# Installer logs
pip-log.txt
pip-delete-this-directory.txt

# Unit test / coverage reports
htmlcov/
.tox/
.nox/
.coverage
.coverage.*
.cache
nosetests.xml
coverage.xml
*.cover
*.py,cover
.hypothesis/
.pytest_cache/
cover/

# Jupyter Notebook
.ipynb_checkpoints

# Environments
.env
.venv
env/
venv/
ENV/
env.bak/
venv.bak/

# IDE
.vscode/
.idea/
*.swp
*.swo

# Project specific
*.h5
*.dat
*.npy
*.npz
"""
    (project_path / ".gitignore").write_text(gitignore_content, encoding="utf-8")

    # Create tests directory
    tests_dir = project_path / "tests"
    tests_dir.mkdir(exist_ok=True)
    (tests_dir / "__init__.py").touch()

    # Create a basic test file
    test_file = tests_dir / "test_basic.py"
    try:
        test_content = [
            f'"""Basic tests for {project_name} package."""',
            "\n\nimport pytest",
            f"\nimport {project_name}\n",
            "\ndef test_import():",
            '    """Test that the package can be imported."""',
            f'    assert hasattr({project_name},"__version__")\n',
        ]
        test_file.write_text("\n".join(test_content), encoding="utf-8")
        logger.debug("Created test file: %s", test_file)
    except Exception as e:
        logger.error("Failed to create test file: %s", e)
        raise

    logger.info("Generated project metadata files")


def output_PyFunction_files_and_construct_project(
    project_dir: Union[str, Path],
    project_name: str,
    lib_function_prefix: str = "",
) -> None:
    """
    Generate a complete Python project from registered PyFunctions.

    This function creates a well-structured Python package with proper packaging,
    documentation, and development setup. It's specifically designed for
    JAX-accelerated numerical computing projects.

    :param project_dir: Directory where the project will be created.
    :param project_name: Name of the Python package (used for imports).
    :param lib_function_prefix: Optional prefix for all generated function names.
    :raises RuntimeError: If project generation fails for other reasons.
    """
    logger.info("Starting project generation for %s in %s", project_name, project_dir)

    try:
        # Convert to Path object and resolve to absolute path
        project_path = Path(project_dir).resolve()

        # Validate inputs and check permissions
        _validate_inputs(project_path, project_name)

        # Create project directory structure
        project_path.mkdir(parents=True, exist_ok=True)

        # Generate Python source files
        py_files = _generate_python_files_and_init(project_path, lib_function_prefix)

        # Generate project metadata and configuration files
        _generate_project_metadata(project_path, project_name)

        logger_str = f"""
Successfully generated project at {project_path}
Generated {len(py_files)} Python modules
"""
        logger.info(logger_str)

    except Exception as e:
        logger.error("Failed to generate project: %s", e, exc_info=True)
        raise RuntimeError(f"Project generation failed: {e}") from e


if __name__ == "__main__":
    # Configure logging
    log_level = logging.DEBUG
    logging.basicConfig(level=log_level)

    logger.info("Running in test mode with sample functions")

    # Clear dictionaries to ensure a clean slate for test run
    PyFunction_dict.clear()
    commondata_params_dict.clear()

    def create_test_function() -> None:
        """Create a test function in PyFunction_dict for testing purposes."""
        # Define function components
        name = "test_function"
        desc = "A test function that doubles the input.\n\n    Args:\n        x: Input array.\n        \n    Returns:\n        The input array multiplied by 2."

        # Create and register the function
        register_PyFunction(
            name=name,
            desc=desc,
            imports=["import jax", "import jax.numpy as jnp"],
            params="x: jnp.ndarray",
            body="return 2.0 * x",
            subdirectory="utils",
            pyfunc_decorators="@jax.jit",
        )

    def create_test_commondata() -> None:
        """Create a test commondata parameter for testing purposes."""
        register_commondata_param(
            name="ETA",
            dtype="float",
            default=2.0,
            description="Symmetric mass ratio",
        )

    create_test_function()
    create_test_commondata()

    output_PyFunction_files_and_construct_project(
        project_dir="project/tmp_project",
        project_name="tmp_project",
        lib_function_prefix="",
    )
