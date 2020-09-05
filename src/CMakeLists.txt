add_subdirectory(tinymesh)

# ----------
# Python module
# ----------
if (BUILD_PYTHON_MODULE)
    message(STATUS "Build Python module")
    include_directories(${CMAKE_CURRENT_LIST_DIR}/tinymesh)
    include_directories(${CMAKE_CURRENT_LIST_DIR}/tinymesh/ext/pybind11)
    add_subdirectory(tinymesh/ext/pybind11)
    pybind11_add_module(${TINYMESH_PYTHON_MODULE} pybind11.cpp)
    set_target_properties(${TINYMESH_PYTHON_MODULE} PROPERTIES OUTPUT_NAME "tinymesh")
    target_link_libraries(${TINYMESH_PYTHON_MODULE} PRIVATE ${TINYMESH_LIBRARY})
endif()
