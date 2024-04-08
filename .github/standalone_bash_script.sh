#!/bin/bash

# Clear cache
rm -rf ~/.cache/nrpy

# Display Python version
python_version=$(python --version)
echo "Python version: $python_version"

# Install dependencies
python -m pip install --upgrade pip
pip install -U -r requirements.txt
pip install -U -r requirements-dev.txt

# Display sympy version
echo "Running CI tests with SymPy version = $(isympy --version)"

# Initialize array to store failed tests
failed_tests=()

# Use find to locate python files
python_files=$(find . -name '*.py' -not -name '__init__.py' -not -path './project/*' -not -path '*/tests/*' -not -path './nrpy/examples/visualization_scripts/*')

# Loop through each python file
for python_file in $python_files; do
  echo ""
  echo "-={ $python_file }=-"
  echo "-={ Step 1: Doctests/run Python module }=-"
  PYTHONPATH=.:$PYTHONPATH DOCTEST_MODE=1 python $python_file || { failed_tests+=("doctest in $python_file"); break; }

  if [[ "$python_version" != *"3.6.7"* ]]; then
    echo "-={ Step 2: black }=-"
    black --check $python_file || { failed_tests+=("black in $python_file"); break; }
    echo "-={ Step 3: mypy }=-"
    PYTHONPATH=.:$PYTHONPATH mypy --strict --pretty --allow-untyped-calls $python_file || { failed_tests+=("mypy in $python_file"); break; }
  fi

  echo "-={ Step 4: pylint }=-"
  pylint_score=$(PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc $python_file | tail -2 | grep -Eo '[0-9\.]+' | head -1 || echo "0")
  echo "Pylint score is $pylint_score"
  if (( $(echo "$pylint_score < 9.91" | bc -l) )); then
    PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc $python_file || true
    echo "Pylint score is below 9.91, failing..."
    failed_tests+=("pylint in $python_file")
    break
  fi

  echo "-={ Step 5: pydocstyle }=-"
  pydocstyle $python_file || { failed_tests+=("pydocstyle in $python_file"); break; }

  echo "-={ Step 6: darglint }=-"
  darglint -v 2 $python_file || { failed_tests+=("darglint in $python_file"); break; }
done

# Exit with failure if any tests failed
if [ ${#failed_tests[@]} -ne 0 ]; then
  echo "The following tests failed: ${failed_tests[*]}"
  exit 1
fi

PYTHONPATH=.:$PYTHONPATH python nrpy/examples/wave_equation_cartesian.py   && (cd project/wavetoy && make && make clean) &&
PYTHONPATH=.:$PYTHONPATH python nrpy/examples/wave_equation_curvilinear.py && (cd project/curviwavetoy && make && make clean) &&
PYTHONPATH=.:$PYTHONPATH python nrpy/examples/nrpyelliptic_conformally_flat.py && (cd project/nrpyelliptic_conformally_flat && make && make clean) &&
PYTHONPATH=.:$PYTHONPATH python nrpy/examples/two_blackholes_collide.py    && (cd project/two_blackholes_collide && make && make clean) &&
PYTHONPATH=.:$PYTHONPATH python nrpy/examples/blackhole_spectroscopy.py    && (cd project/blackhole_spectroscopy && make && make clean) &&
PYTHONPATH=.:$PYTHONPATH python nrpy/examples/spinning_blackhole.py        && (cd project/spinning_blackhole && make && make clean) &&
PYTHONPATH=.:$PYTHONPATH python nrpy/examples/nrpypn_quasicircular_momenta.py && (cd project/nrpypn_quasicircular_momenta && make && make clean) &&
PYTHONPATH=.:$PYTHONPATH python nrpy/examples/wave_equation_multicoord_wavetoy.py && (cd project/multicoords_curviwavetoy && make && make clean)
