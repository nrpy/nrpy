#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Ensure the script is called with a single file argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <python_file>"
  exit 1
fi

python_file=$1

# Check if the file exists
if [ ! -f "$python_file" ]; then
  echo "File '$python_file' does not exist."
  exit 1
fi

# Display Python version
python_version=$(python --version 2>&1)
echo "Python version: $python_version"

# Initialize array to store failed tests
failed_tests=()

# Initialize step counter
step_counter=1

echo ""
echo "-={ $python_file }=-"

# Function to handle each test step
run_test_step() {
  local step_name=$1
  local command=$2
  local failure_message=$3

  echo "-={ Step $step_counter: $step_name }=-"
  eval "$command"
  if [ $? -ne 0 ]; then
    failed_tests+=("$failure_message in $python_file")
  fi
  ((step_counter++))
}

# Static analysis steps: black, isort, mypy, pylint, pydocstyle, darglint

# Step 1: black
run_test_step "black" "black --check \"$python_file\"" "black"

# Step 2: isort
run_test_step "isort" "isort --check-only \"$python_file\"" "isort"

# Step 3: mypy
run_test_step "mypy" "PYTHONPATH=.:$PYTHONPATH mypy --strict --pretty --allow-untyped-calls \"$python_file\"" "mypy"

# Step 4: pylint
echo "-={ Step $step_counter: pylint }=-"
pylint_score=$(PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc "$python_file" | tail -2 | grep -Eo '[0-9\.]+' | head -1 || echo "0")
echo "Pylint score is $pylint_score"
if (( $(echo "$pylint_score < 9.91" | bc -l) )); then
  PYTHONPATH=.:$PYTHONPATH pylint --rcfile=.pylintrc "$python_file" || true
  echo "Pylint score is below 9.91, failing..."
  failed_tests+=("pylint in $python_file")
fi
((step_counter++))

# Step 5: pydocstyle
run_test_step "pydocstyle" "pydocstyle \"$python_file\"" "pydocstyle"

# Step 6: darglint
run_test_step "darglint" "darglint -v 2 \"$python_file\"" "darglint"

# Step 7: doctests
echo "-={ Step $step_counter: doctests }=-"
# Capture the output of the doctest run
doctest_output=$(python "$python_file" 2>&1)
doctest_exit_code=$?
echo "$doctest_output"

if [ $doctest_exit_code -ne 0 ]; then
  failed_tests+=("doctests in $python_file")
fi
((step_counter++))

# Exit with failure if any tests failed
if [ ${#failed_tests[@]} -ne 0 ]; then
  echo ""
  echo "The following tests failed: ${failed_tests[*]}"
  exit 1
fi

echo ""
echo "All static analysis tests passed successfully."
