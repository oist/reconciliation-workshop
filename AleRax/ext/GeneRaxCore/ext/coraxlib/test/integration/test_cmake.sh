#!/bin/bash

EXIT_SUCCESS=0
EXIT_FAIL=1

set -e
set -o pipefail

if [ ! -d "src/corax/" ]; then
    echo "Please call this script from the project's root directory."
    exit 1
fi

# shellcheck disable=SC1091
source test/integration/assert.sh

fail() {
    echo "$1"
    exit $EXIT_FAIL
}

assert_file_exists() {
    file="$1"
    msg="$2"
    if [ ! -f "$file" ]; then
        echo "$msg"
        exit $EXIT_FAIL
    fi
}

BUILD_DIR="$(mktemp --directory "$(pwd)/build.XXXXX")"
LIB_OUT_DIR="$(mktemp --directory "$(pwd)/lib.XXXXX")"

# First, try the default configuration
echo -n "Configuring and building default configuration ... "
configure_output="$(cmake -B "$BUILD_DIR" -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY="$LIB_OUT_DIR" || fail "Configuration failed")"
assert_contain "$(echo "$configure_output" | grep "Build type:")" "Release" "The default build type should be Release"
assert_contain "$(echo "$configure_output" | grep "Build coraxlib as a ")" "static library" "We should default to building a static library"

cmake --build "$BUILD_DIR" -j >/dev/null 2>&1 || fail "Building failed"
assert_file_exists "$LIB_OUT_DIR/libcorax.a" "Static library not found"

rm -r "$BUILD_DIR"
echo "done"

# Build as shared library, in Debug mode
echo -n "Configuring and building shared library configuration ... "
configure_output="$(cmake -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Debug -DCORAX_BUILD_SHARED=On -DCMAKE_LIBRARY_OUTPUT_DIRECTORY="$LIB_OUT_DIR" || fail "Configuration failed")"
assert_contain "$(echo "$configure_output" | grep "Build type:")" "Debug" "Did not build in Debug mode"
assert_contain "$(echo "$configure_output" | grep "Build coraxlib as a ")" "shared library" "Did not build as shared library"

cmake --build "$BUILD_DIR" -j >/dev/null 2>&1 || fail "Building failed"
assert_file_exists "$LIB_OUT_DIR/libcorax.so" "Shared library not found"

rm -r "$BUILD_DIR"
echo "done"

# Build as benchmarks, tests, documentation, and difficulty prediction
echo -n "Configuring and building with everything enabled ... "
configure_output="$(cmake \
    -B "$BUILD_DIR" \
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY="$LIB_OUT_DIR"\
    -DCORAX_BUILD_BENCHMARKS=On\
    -DCORAX_BUILD_DIFFICULTY_PREDICTION=On\
    -DCORAX_BUILD_DOCS=On\
    -DCORAX_BUILD_TESTS=On || \
    -DCMAKE_BUILD_TYPE=Debug \
    fail "Configuration failed")"
assert_contain "$(echo "$configure_output" | grep "Building benchmarks ")" "On" "Did not build benchmarks"
assert_contain "$(echo "$configure_output" | grep "Building documentation ")" "On" "Did not build documentation"
assert_contain "$(echo "$configure_output" | grep "Building tests ")" "On" "Did not build tests"
assert_contain "$(echo "$configure_output" | grep "Building difficulty prediction ")" "On" "Did not build difficulty prediction"

cmake --build "$BUILD_DIR" -j >/dev/null 2>&1 || fail "Building failed"
assert_file_exists "$LIB_OUT_DIR/libcoraxlib_difficulty_prediction_lib.a" "Difficulty prediction library not found"
assert_file_exists "$BUILD_DIR/test/unit/test_dummy" "test_dummy not found -> test not built?"

rm -r "$BUILD_DIR"
echo "done"

# Build the example cmake-project
echo -n "Configuring and building the example cmake-project ... "
cd examples/cmake-project || exit $EXIT_FAIL
cmake -B build >/dev/null 2>&1 || exit $EXIT_FAIL
cmake --build build >/dev/null 2>&1 || exit $EXIT_FAIL
assert_contain "$(build/newton | grep -e "-\*- .* -\*-")" "Optimizing branch length" "Running newton in the example project failed" 
rm -r build
cd ../..
echo "done"

exit $EXIT_SUCCESS
