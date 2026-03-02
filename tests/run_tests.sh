#!/bin/bash

# Check that all headers compile separately. This avoids problems like headers
# that need to be included in a specific order for the code to compile
rm -rf tests/compile/build
cmake -S tests/compile -B tests/compile/build || exit 1
cmake --build tests/compile/build || exit 1

exit 0
