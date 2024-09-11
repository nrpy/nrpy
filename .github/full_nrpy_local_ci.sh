#!/bin/bash

# Clear cache
rm -rf ~/.cache/nrpy

# Display Python version
python_version=$(python --version)
echo "Python version: $python_version"

# Install dependencies
python -m pip install --upgrade pip
pip install -U -r requirements.txt

if [[ "$python_version" != *"3.6.7"* && "$python_version" != *"3.7.13"* ]]; then
  pip install -U -r requirements-dev.txt
else
  # Install packages from requirements-dev.txt except those with a fixed version.
  grep -v '==' requirements-dev.txt | xargs pip install -U
  pip install -U clang-format
fi

# Install NRPyLaTeX and additional dependencies
pip install -U git+https://github.com/nrpy/nrpylatex.git
pip install -U ipython setuptools

# Install SymPy or development SymPy based on Python version
# if [[ "$python_version" != *"3.6.7"* && "$python_version" != *"3.7.13"* ]]; then
#   pip install git+https://github.com/sympy/sympy.git
# else
#   pip install sympy
# fi

# Display SymPy version
echo "Running CI tests with SymPy version = $(isympy --version)"

# Initialize array to store failed tests
failed_tests=()

# Use find to locate python files
python_files=$(find . -name '*.py' -not -name '__init__.py' -not -path './project/*' -not -path './build/*' -not -path '*/tests/*' -not -path './nrpy/examples/visualization_scripts/*')

# Loop through each python file
for python_file in $python_files; do
  # Initialize step counter
  step_counter=1

  echo ""
  echo "-={ $python_file }=-"

  if [[ ! $python_file =~ nrpy/examples/.* ]]; then
    echo "-={ Step $step_counter: Doctests/run Python module }=-"
    DOCTEST_MODE=1 PYTHONPATH=.:$PYTHONPATH python $python_file || { failed_tests+=("doctest in $python_file"); break; }
    ((step_counter++))
  fi

  if [[ "$python_version" != *"3.6.7"* && "$python_version" != *"3.7.13"* && "$python_version" != *"3.8.18"* ]]; then
    echo "-={ Step $step_counter: black }=-"
    black --check $python_file || { failed_tests+=("black in $python_file"); break; }
    ((step_counter++))

    echo "-={ Step $step_counter: isort }=-"
    isort --check-only $python_file || { failed_tests+=("isort in $python_file"); break; }
    ((step_counter++))

    echo "-={ Step $step_counter: mypy }=-"
    PYTHONPATH=.:$PYTHONPATH mypy --strict --pretty --allow-untyped-calls $python_file || { failed_tests+=("mypy in $python_file"); break; }
    ((step_counter++))
  fi

  echo "-={ Step $step_counter: pylint }=-"
  pylint_score="0"
  if [[ "$python_version" == *"3.6.7"* || "$python_version" == *"3.7.13"* ]]; then
    pylint_score=$(PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc_python36 $python_file | tail -2 | grep -Eo '[0-9\.]+' | head -1 || echo "0")
  else
    pylint_score=$(PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc $python_file | tail -2 | grep -Eo '[0-9\.]+' | head -1 || echo "0")
  fi
  echo "Pylint score is $pylint_score"
  if (( $(echo "$pylint_score < 9.91" | bc -l) )); then
    PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc $python_file || true
    echo "Pylint score is below 9.91, failing..."
    failed_tests+=("pylint in $python_file")
    break
  fi
  ((step_counter++))

  echo "-={ Step $step_counter: pydocstyle }=-"
  pydocstyle $python_file || { failed_tests+=("pydocstyle in $python_file"); break; }
  ((step_counter++))

  echo "-={ Step $step_counter: darglint }=-"
  darglint -v 2 $python_file || { failed_tests+=("darglint in $python_file"); break; }
  ((step_counter++))
done

# Exit with failure if any tests failed
if [ ${#failed_tests[@]} -ne 0 ]; then
  echo "The following tests failed: ${failed_tests[*]}"
  exit 1
fi

# Run example scripts and compile related projects
example_scripts=(
  "nrpy/examples/wave_equation_cartesian.py project/wavetoy"
  "nrpy/examples/wave_equation_curvilinear.py project/curviwavetoy"
  "nrpy/examples/nrpyelliptic_conformally_flat.py project/nrpyelliptic_conformally_flat"
  "nrpy/examples/two_blackholes_collide.py project/two_blackholes_collide"
  "nrpy/examples/blackhole_spectroscopy.py project/blackhole_spectroscopy"
  "nrpy/examples/spinning_blackhole.py project/spinning_blackhole"
  "nrpy/examples/nrpypn_quasicircular_momenta.py project/nrpypn_quasicircular_momenta"
  "nrpy/examples/wave_equation_multicoord_wavetoy.py project/multicoords_curviwavetoy"
  "nrpy/examples/carpet_baikal_thorns.py project/carpet_baikal_thorns"
  "nrpy/examples/carpet_wavetoy_thorns.py project/carpet_wavetoy_thorns"
  "nrpy/examples/carpetx_baikal_thorns.py project/carpetx_baikal_thorns"
  "nrpy/examples/carpetx_wavetoy_thorns.py project/carpetx_wavetoy_thorns"
  "nrpy/examples/seobnrv5_aligned_spin_inspiral.py project/seobnrv5_aligned_spin_inspiral"
  "nrpy/examples/superB_two_blackholes_collide.py project/superB_two_blackholes_collide"
  "nrpy/examples/superB_blackhole_spectroscopy.py project/superB_blackhole_spectroscopy"
  "nrpy/examples/tovola_neutron_star.py project/tovola_neutron_star"
  "nrpy/examples/hydro_without_hydro project/hydro_without_hydro"
  "nrpy/examples/manga_bhah_lib project/bhah_lib"
)

for script in "${example_scripts[@]}"; do
  IFS=' ' read -r script_path project_path <<< "$script"
  PYTHONPATH=.:$PYTHONPATH python $script_path
  if [[ $? -ne 0 ]]; then
    echo "Error: Python script $script_path failed."
    exit 1
  fi
  if [[ $script_path != *"superB"* && $script_path != *"carpet"* ]]; then
    (cd $project_path && make && make clean)
    if [[ $? -ne 0 ]]; then
      echo "Error: Compilation in $project_path failed."
      exit 1
    fi
  fi
done
