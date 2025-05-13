# Some code taken from https://vicrucann.github.io/tutorials/quick-cmake-doxygen/

find_package(Doxygen)

if (DOXYGEN_FOUND)
    # Set input and output files.
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # Request to configure the Doxygen file.
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    # The option ALL allows to build the docs together with the application.
    add_custom_target(
        docs ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
else (DOXYGEN_FOUND)
    message(WARNING "Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
