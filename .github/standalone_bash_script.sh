#!/bin/bash

# Clear cache
rm -rf ~/.cache/nrpy

# Display Python version
python --version

# Install dependencies
python -m pip install --upgrade pip
if [ -f requirements.txt ]; then
  pip install -U -r requirements.txt
fi
pip install -U mypy==1.8.0 pylint black==24.1.1 clang-format ipython setuptools

# Display sympy version
echo "Running CI tests with SymPy version = $(isympy --version)"

# Initialize array to store failed tests
failed_tests=()

# Use find to locate python files based on pattern or directory structure.
#   Don't analyze Python scripts in tests/ (though they should pass!)
#   Ignore CarpetX infrastructure for now as well...
python_files=$(find . -name '*.py' -not -name '__init__.py' -not -path './build/*' -not -path './nrpy/infrastructures/CarpetX/*' -not -path '*/tests/*')

# Loop through each python file
for python_file in $python_files; do
  echo ""
  echo "-={ $python_file }=-"
  echo "-={ Step 1: Doctests/run Python module }=-"
  PYTHONPATH=.:$PYTHONPATH DOCTEST_MODE=1 python $python_file || { failed_tests+=("doctest in $python_file"); break; }
  if [ "$(python --version 2>&1 | cut -d ' ' -f2)" != "3.6.7" ]; then
    # Turns out that black in Python 3.6.7 has a heart attack when parsing equations/general_relativity/BSSN_quantities.py:
    # INTERNAL ERROR: Black produced code that is not equivalent to the source. Please report a bug on .... [HOW ABOUT NO. BEGGING FOR WONTFIX]
    echo "-={ Step 2: black }=-"
    black --check $python_file || { failed_tests+=("black in $python_file"); break; }
    echo "-={ Step 3: mypy }=-"
    PYTHONPATH=.:$PYTHONPATH mypy --strict --pretty --allow-untyped-calls $python_file || { failed_tests+=("mypy in $python_file"); break; }
  fi
  echo "-={ Step 4: pylint }=-"
  PYTHONPATH=.:$PYTHONPATH pylint_score=$(pylint --rcfile=.pylintrc $python_file | tail -2 | grep -Eo '[0-9\.]+' | head -1 || echo "0")
  echo "Pylint score is $pylint_score"
  if (( $(echo "$pylint_score < 9.91" | bc -l) )); then
    PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc $python_file || true
    echo "Pylint score is below 9.91, failing..."
    failed_tests+=("pylint in $python_file")
    break
  fi
  if (( $(echo "$pylint_score < 10.00" | bc -l) )); then
    PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc $python_file || true
  fi
  echo "-={ Step 5: pydocstyle }=-"  # Added pydocstyle step
  pydocstyle $python_file || { failed_tests+=("pydocstyle in $python_file"); break; }
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
