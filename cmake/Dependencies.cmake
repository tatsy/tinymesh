include_guard(GLOBAL)
include(FetchContent)

# Preserve the commits previously recorded by the Git submodules. These four
# dependencies are used as headers; do not run their upstream CMake projects.
# tinyobjloader and tinyply implementations remain in core/external.cpp.
FetchContent_Declare(eigen
  URL https://gitlab.com/libeigen/eigen/-/archive/3147391d946bb4b6c68edd901f2add6ac1f31f8c/eigen-3147391d946bb4b6c68edd901f2add6ac1f31f8c.tar.gz
  URL_HASH SHA256=0c8c490764f9c2a793133491adca0cd073b73e0bde965c68cbe58d91b5ed4261
  DOWNLOAD_EXTRACT_TIMESTAMP FALSE
  SOURCE_SUBDIR _tinymesh_headers_only
)
FetchContent_Declare(spectra
  URL https://github.com/yixuan/spectra/archive/8c7242e08f0fb7f6a0022cfe232e3dc5b5bd4eb4.tar.gz
  URL_HASH SHA256=2328f90afc2594ccc82de06b09ff927b9caf5b569e85db52be4bccffa6c919d9
  DOWNLOAD_EXTRACT_TIMESTAMP FALSE
  SOURCE_SUBDIR _tinymesh_headers_only
)
FetchContent_Declare(tinyobjloader
  URL https://github.com/syoyo/tinyobjloader/archive/f37fed32f3eb0912cc10a970f78774cd98598ef6.tar.gz
  URL_HASH SHA256=040b850d033909dc150eadb7b6c9a3167bdaf852e2b7498e3ed024dcafbb6091
  DOWNLOAD_EXTRACT_TIMESTAMP FALSE
  SOURCE_SUBDIR _tinymesh_headers_only
)
FetchContent_Declare(tinyply
  URL https://github.com/ddiakopoulos/tinyply/archive/a6ca3c89480b89d271b6909e16d8423ef0025336.tar.gz
  URL_HASH SHA256=d6f5d03d4593721df57a0712816c878faa5436f0d5b579bf921acdaeb42c6f29
  DOWNLOAD_EXTRACT_TIMESTAMP FALSE
  SOURCE_SUBDIR _tinymesh_headers_only
)
FetchContent_MakeAvailable(eigen spectra tinyobjloader tinyply)

add_library(tinymesh_eigen INTERFACE)
add_library(tinymesh::eigen ALIAS tinymesh_eigen)
target_include_directories(tinymesh_eigen SYSTEM INTERFACE "${eigen_SOURCE_DIR}")

add_library(tinymesh_spectra INTERFACE)
add_library(tinymesh::spectra ALIAS tinymesh_spectra)
target_include_directories(tinymesh_spectra SYSTEM INTERFACE "${spectra_SOURCE_DIR}/include")
target_link_libraries(tinymesh_spectra INTERFACE tinymesh::eigen)

add_library(tinymesh_tinyobjloader INTERFACE)
add_library(tinymesh::tinyobjloader ALIAS tinymesh_tinyobjloader)
target_include_directories(tinymesh_tinyobjloader SYSTEM INTERFACE "${tinyobjloader_SOURCE_DIR}")

add_library(tinymesh_tinyply INTERFACE)
add_library(tinymesh::tinyply ALIAS tinymesh_tinyply)
target_include_directories(tinymesh_tinyply SYSTEM INTERFACE "${tinyply_SOURCE_DIR}/source")

if(TINYMESH_BUILD_TESTS)
  if(MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  endif()
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    URL_HASH SHA256=edd885a1ab32b6999515a880f669efadb80b3f880215f315985fa3f6eca7c4d3
    DOWNLOAD_EXTRACT_TIMESTAMP FALSE
  )
  FetchContent_MakeAvailable(googletest)
endif()

if(TINYMESH_BUILD_PYTHON)
  set(PYBIND11_INSTALL OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(pybind11
    URL https://github.com/pybind/pybind11/archive/refs/tags/v3.1.0.tar.gz
    URL_HASH SHA256=ef712655692a2e9bf7bb7874c022564a45f91d847ddee987e720cd9e28849665
    DOWNLOAD_EXTRACT_TIMESTAMP FALSE
  )
  FetchContent_MakeAvailable(pybind11)
endif()
