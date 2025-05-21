list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(CoraxlibTestHelper)

# Load the GoogleTest included in the git submodule.
add_subdirectory("${PROJECT_SOURCE_DIR}/lib/googletest" "lib/googletest")
set(INSTALL_GTEST Off)
include(GoogleTest)

# Keep cache nad IDEs which support folders clean. See:
# https://cliutils.gitlab.io/modern-cmake/chapters/testing/googletest.html
mark_as_advanced(
    BUILD_GMOCK
    BUILD_GTEST
    BUILD_SHARED_LIBS
    gmock_build_tests
    gtest_build_samples
    gtest_build_tests
    gtest_disable_pthreads
    gtest_force_shared_crt
    gtest_hide_internal_symbols
    benchmark_DIR
    INSTALL_GTEST
)

set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)
set_target_properties(gmock PROPERTIES FOLDER extern)
set_target_properties(gmock_main PROPERTIES FOLDER extern)

# Add the unit tests.
add_subdirectory(unit)
add_subdirectory(integration EXCLUDE_FROM_ALL)
