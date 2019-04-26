cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(TinyMesh)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(Macros)

# ----------
# Options
# ----------
option(BUILD_STATIC_LIB "Build static link library" OFF)
option(BUILD_EXAMPLES "Build examples." OFF)
option(BUILD_TESTS "Build unit tests." OFF)

if (NOT MSVC)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS " -Xlinker -w")
    endif()
endif()

# ----------
# Common target names
# ----------
set(TINYMESH_LIBRARY "tinymesh")
set(TINYMESH_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/src")

# ----------
# Output paths
# ----------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_DEBUG_POSTFIX "-debug")

# ----------
# OS specific settings
# ----------
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if (APPLE)
    set(CMAKE_MACOSX_RPATH 1)
endif()

# ----------
# Traverse subdirectories
# ----------
add_subdirectory(src)

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if (BUILD_TESTS)
    find_package(GoogleTest REQUIRED)
    add_subdirectory(tests)
endif()
