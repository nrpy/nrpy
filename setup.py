"""
Install the nrpy package.

To install the nrpy package, navigate to this directory and execute:
    pip install .
This will install nrpy and its required dependencies.

Instructions for uploading latest release to PyPI:
    rm -rf build dist && python setup.py sdist bdist_wheel && twine check dist/*
    twine upload dist/*
"""

import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

from setuptools import find_packages, setup

# pylint: disable=consider-using-f-string


def check_python_version() -> None:
    """
    Check for the minimum Python version (3.6 or newer).

    :raises SystemExit: If the Python version is less than 3.6.
    """
    if sys.version_info < (3, 6):
        raise SystemExit(
            "This project requires Python 3.6 and newer. Python {0}.{1} detected.".format(
                sys.version_info[0], sys.version_info[1]
            )
        )


def clang_format_is_installed() -> bool:
    """
    Check if clang-format is installed.

    :return: True if clang-format is installed, False otherwise.
    """
    try:
        subprocess.run(
            ["clang-format", "--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def read_requirements_file() -> List[str]:
    """
    Read the contents of the requirements.txt file.

    :return: List of strings containing the required packages.
    """
    with open("requirements.txt", "r", encoding="utf-8") as file:
        return file.read().splitlines()


def get_nrpy_version(pkg_root_directory: str) -> str:
    """
    Fetch the version from the nrpy package.

    :param pkg_root_directory: Root directory where 'release.txt' is located.
    :return: Version as a string.
    :raises ValueError: When version information could not be found.
    """
    with open(
        Path(pkg_root_directory) / "nrpy" / "release.txt", encoding="utf-8"
    ) as file:
        for line in file:
            if line.startswith("version ="):
                return line.split("=")[1].strip().strip("\"'")
    raise ValueError("Version information could not be found in 'release.txt'")


def discover_header_package_data(
    pkg_root_directory: str, top_package: str = "nrpy"
) -> Dict[str, List[str]]:
    """
    Discover C header files under the package directory and map them into package_data.

    Recursively searches for ``*.h`` under ``top_package`` within ``pkg_root_directory``,
    skipping any paths that contain a ``tests`` component. Each header is attached to the
    nearest ancestor directory that is a Python package (i.e., contains ``__init__.py``).

    :param pkg_root_directory: Repository root containing the top-level package directory.
    :param top_package: The top-level package name to scan (defaults to "nrpy").
    :return: Mapping from dotted package name to a sorted list of header paths relative to that package directory.
    """
    root_path = Path(pkg_root_directory).resolve()
    top_path = (root_path / top_package).resolve()

    # Build a map from absolute package directories to dotted package names.
    packages = find_packages(where=pkg_root_directory)
    pkg_dir_to_name = {
        (root_path / pkg.replace(".", "/")).resolve(): pkg for pkg in packages
    }

    result: Dict[str, List[str]] = defaultdict(list)

    # Find all headers under nrpy/, skip anything in a tests/ directory.
    for header in top_path.rglob("*.h"):
        if any(part == "tests" for part in header.parts):
            continue

        # Find the nearest ancestor that is a Python package dir.
        cur = header.parent.resolve()
        nearest_pkg_dir = None
        while True:
            if cur in pkg_dir_to_name:
                nearest_pkg_dir = cur
                break
            if cur == root_path:
                break
            cur = cur.parent

        if nearest_pkg_dir is None:
            # Fallback: attach relative to top_package if no ancestor has __init__.py
            pkg_name = top_package
            rel = header.relative_to(top_path).as_posix()
        else:
            pkg_name = pkg_dir_to_name[nearest_pkg_dir]
            rel = header.relative_to(nearest_pkg_dir).as_posix()

        result[pkg_name].append(rel)

    # Sort for determinism
    return {pkg: sorted(files) for pkg, files in result.items()}


if __name__ == "__main__":
    # Don't install NRPy if this is run from a doctest.
    if "DOCTEST_MODE" in os.environ:
        sys.exit(0)

    check_python_version()

    dir_setup = os.path.dirname(os.path.realpath(__file__))

    requirements = read_requirements_file()
    if not clang_format_is_installed():
        requirements.append("clang-format")

    # Auto-discover header files
    auto_pkg_data = discover_header_package_data(dir_setup, top_package="nrpy")

    # Ensure py.typed is included for PEP 561 typing
    auto_pkg_data.setdefault("nrpy", [])
    if "py.typed" not in auto_pkg_data["nrpy"]:
        auto_pkg_data["nrpy"].append("py.typed")

    setup(
        name="nrpy",
        version=get_nrpy_version(dir_setup),
        license="BSD-2-Clause",
        data_files=[("license", ["LICENSE"])],
        description="Python/SymPy-based code generation for numerical relativity... and beyond!",
        long_description=(Path(__file__).parent / "README.md").read_text(
            encoding="utf-8"
        ),
        long_description_content_type="text/markdown",
        python_requires=">=3.6",
        packages=find_packages(),
        package_data=auto_pkg_data,  # ← use the discovered headers
        classifiers=[
            "License :: OSI Approved :: BSD License",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Astronomy",
            "Topic :: Scientific/Engineering :: Mathematics",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering :: Visualization",
            "Topic :: Software Development :: Code Generators",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python :: 3.12",
            "Programming Language :: Python :: 3 :: Only",
            "Programming Language :: Python :: Implementation :: CPython",
            "Programming Language :: Python :: Implementation :: PyPy",
        ],
        url="https://github.com/nrpy/nrpy",
        author="Zachariah Etienne",
        install_requires=requirements,
        entry_points={
            "console_scripts": [
                "nrpyinline=bin.nrpyinline:main",
            ],
        },
        include_package_data=True,
        zip_safe=False,
    )
