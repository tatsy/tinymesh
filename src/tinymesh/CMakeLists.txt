# ----------
# Components
# ----------
set(SOURCE_FILES ${SOURCE_FILES} "tinymesh.h")
add_folder(SOURCE_FILES "core")
add_folder(SOURCE_FILES "filters")
add_folder(SOURCE_FILES "ops")
add_folder(SOURCE_FILES "remesh")
add_folder(SOURCE_FILES "restore")

if (TINYMESH_BUILD_SHARED)
  add_library(${TINYMESH_LIBRARY} SHARED ${SOURCE_FILES})
  target_compile_definitions(${TINYMESH_LIBRARY} PUBLIC TINYMESH_SHARED PRIVATE TINYMESH_API_EXPORT)
else()
  add_library(${TINYMESH_LIBRARY} STATIC ${SOURCE_FILES})
endif()

target_compile_features(${TINYMESH_LIBRARY} PUBLIC cxx_std_20)
set_target_properties(${TINYMESH_LIBRARY} PROPERTIES POSITION_INDEPENDENT_CODE ON)

if (TINYMESH_BUILD_PYTHON)
  target_compile_definitions(${TINYMESH_LIBRARY}
    PUBLIC TINYMESH_PYTHON_MODULE=1)

  target_link_libraries(${TINYMESH_LIBRARY} PUBLIC pybind11::module)

  if (MSVC)
    target_compile_options(${TINYMESH_LIBRARY} PUBLIC /utf-8)
    target_compile_definitions(${TINYMESH_LIBRARY} PUBLIC _CRT_SECURE_NO_WARNINGS)
  endif()
endif()

if(MSVC)
  target_compile_options(${TINYMESH_LIBRARY} PRIVATE /MP /wd4251 /wd4819)
else()
  target_compile_options(${TINYMESH_LIBRARY} PRIVATE -Wall)
endif()

source_group("Source Files" FILES ${SOURCE_FILES})

target_include_directories(${TINYMESH_LIBRARY}
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/.."
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(${TINYMESH_LIBRARY}
  PUBLIC
    tinymesh::eigen
  PRIVATE
    tinymesh::spectra
    tinymesh::tinyobjloader
    tinymesh::tinyply
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(${TINYMESH_LIBRARY} PUBLIC OpenMP::OpenMP_CXX)
endif()

# Disable warning for external files
set(WARNING_DISABLED_CXX_FLAGS)
if(MSVC)
  set(WARNING_DISABLED_CXX_FLAGS "${CMAKE_CXX_FLAGS} /w")
else()
  set(WARNING_DISABLED_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif()

set_source_files_properties(
  ${CMAKE_CURRENT_LIST_DIR}/core/external.cpp
  PROPERTIES
  COMPILE_FLAGS ${WARNING_DISABLED_CXX_FLAGS}
)
