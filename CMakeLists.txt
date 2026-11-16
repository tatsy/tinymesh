cmake_minimum_required(VERSION 3.24.0 FATAL_ERROR)
project(TinyMesh LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(Macros)

# ----------
# Options
# ----------
option(TINYMESH_WITH_OPENMP "Use OpenMP or not" ON)
option(TINYMESH_BUILD_TESTS "Build unit tests." OFF)
option(TINYMESH_BUILD_EXAMPLES "Build examples." OFF)
option(TINYMESH_BUILD_SHARED "Build shared library." OFF)
option(TINYMESH_BUILD_PYTHON "Build the Python module." OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ----------
# Common target names
# ----------
set(TINYMESH_LIBRARY "tinymesh")
set(TINYMESH_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/src")

# ----------
# OpenMP
# ----------
if(TINYMESH_WITH_OPENMP)
  find_package(OpenMP COMPONENTS CXX)
endif()

# ----------
# Output paths
# ----------
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin)
endif()

if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin)
endif()

if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/lib)
endif()

set(CMAKE_DEBUG_POSTFIX "-debug")

# ----------
# OS specific settings
# ----------
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(APPLE)
  set(CMAKE_MACOSX_RPATH 1)
endif()

# ----------
# Traverse subdirectories
# ----------
if (TINYMESH_BUILD_PYTHON)
  if (TINYMESH_BUILD_TESTS OR TINYMESH_BUILD_EXAMPLES)
    message(FATAL_ERROR "Build Python bindings separately from C++ tests and examples.")
  endif()

  find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
endif()

include(Dependencies)

add_subdirectory(src)

if(TINYMESH_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(TINYMESH_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
