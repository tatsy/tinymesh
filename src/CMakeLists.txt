add_subdirectory(tinymesh)

# ----------
# Python module
# ----------
if (BUILD_PYTHON_MODULE)
    message(STATUS "Build Python module")
    include_directories(${CMAKE_CURRENT_LIST_DIR}/tinymesh)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/tinymesh/ext/pybind11)
    add_subdirectory(tinymesh/ext/pybind11)

    add_library(${TINYMESH_PYTHON_MODULE} SHARED pybind11.cpp)
    target_link_libraries(${TINYMESH_PYTHON_MODULE} PRIVATE pybind11::module ${TINYMESH_LIBRARY})
    set_target_properties(${TINYMESH_PYTHON_MODULE} PROPERTIES
                          PREFIX "${PYTHON_MODULE_PREFIX}"
                          SUFFIX "${PYTHON_MODULE_EXTENSION}")

    if (NOT PYTHON_MODULE_OUTPUT_DIRECTORY)
        set(PYTHON_MODULE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python)
    endif()

    set_target_properties(${TINYMESH_PYTHON_MODULE} PROPERTIES
                          OUTPUT_NAME "tinymesh"
                          RUNTIME_OUTPUT_DIRECTORY "${PYTHON_MODULE_OUTPUT_DIRECTORY}"
                          RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PYTHON_MODULE_OUTPUT_DIRECTORY}"
                          RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PYTHON_MODULE_OUTPUT_DIRECTORY}"
                          ARCHIVE_OUTPUT_DIRECTORY "${PYTHON_MODULE_OUTPUT_DIRECTORY}"
                          ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${PYTHON_MODULE_OUTPUT_DIRECTORY}"
                          ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${PYTHON_MODULE_OUTPUT_DIRECTORY}"
                          LIBRARY_OUTPUT_DIRECTORY "${PYTHON_MODULE_OUTPUT_DIRECTORY}"
                          LIBRARY_OUTPUT_DIRECTORY_DEBUG "${PYTHON_MODULE_OUTPUT_DIRECTORY}"
                          LIBRARY_OUTPUT_DIRECTORY_RELEASE "${PYTHON_MODULE_OUTPUT_DIRECTORY}")
endif()
