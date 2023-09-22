#!/bin/bash

# Display Python version
python --version

# Install dependencies
python -m pip install --upgrade pip
if [ -f requirements.txt ]; then
  pip install -r requirements.txt
fi
pip install mypy pylint black clang-format

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
  if (( $(echo "$pylint_score < 9.5" | bc -l) )); then
    PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc $python_file || true
    echo "Pylint score is below 9.5, failing..."
    failed_tests+=("pylint in $python_file")
    break
  fi
  if (( $(echo "$pylint_score < 10.00" | bc -l) )); then
    PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc $python_file || true
  fi
done

# Exit with failure if any tests failed
if [ ${#failed_tests[@]} -ne 0 ]; then
  echo "The following tests failed: ${failed_tests[*]}"
  exit 1
fi
