add_custom_target(
    test_cmake
    COMMAND test/integration/test_cmake.sh
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    COMMENT "Testing CMake build process"
)
