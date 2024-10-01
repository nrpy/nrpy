#!/bin/bash

# Ensure the script is called with a single file argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <python_file>"
  exit 1
fi

python_file=$1

# Check if the file exists
if [ ! -f "$python_file" ]; then
  echo "File $python_file does not exist."
  exit 1
fi

# Display Python version
python_version=$(python --version)
echo "Python version: $python_version"

# Initialize array to store failed tests
failed_tests=()

# Initialize step counter
step_counter=1

echo ""
echo "-={ $python_file }=-"

# Static analysis steps (doctests, black, isort, mypy, pylint, pydocstyle, darglint)
echo "-={ Step $step_counter: black }=-"
black --check $python_file || { failed_tests+=("black in $python_file"); }
((step_counter++))

echo "-={ Step $step_counter: isort }=-"
isort --check-only $python_file || { failed_tests+=("isort in $python_file"); }
((step_counter++))

echo "-={ Step $step_counter: mypy }=-"
PYTHONPATH=.:$PYTHONPATH mypy --strict --pretty --allow-untyped-calls $python_file || { failed_tests+=("mypy in $python_file"); }
((step_counter++))

echo "-={ Step $step_counter: pylint }=-"
pylint_score=$(PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc $python_file | tail -2 | grep -Eo '[0-9\.]+' | head -1 || echo "0")
echo "Pylint score is $pylint_score"
if (( $(echo "$pylint_score < 9.91" | bc -l) )); then
  PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc $python_file || true
  echo "Pylint score is below 9.91, failing..."
  failed_tests+=("pylint in $python_file")
fi
((step_counter++))

echo "-={ Step $step_counter: pydocstyle }=-"
pydocstyle $python_file || { failed_tests+=("pydocstyle in $python_file"); }
((step_counter++))

echo "-={ Step $step_counter: darglint }=-"
darglint -v 2 $python_file || { failed_tests+=("darglint in $python_file"); }
((step_counter++))

# Exit with failure if any tests failed
if [ ${#failed_tests[@]} -ne 0 ]; then
  echo "The following tests failed: ${failed_tests[*]}"
  exit 1
fi

echo "All static analysis tests passed successfully."
