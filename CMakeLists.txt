cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(tinymesh)

# ----------
# Options
# ----------
option(BUILD_STATIC_LIB "Build static link library" OFF)
option(BUILD_UI "Build user interface." OFF)
option(BUILD_EXAMPLES "Build examples." OFF)

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
