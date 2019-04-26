set(BUILD_TARGET ${TINYMESH_LIBRARY})

# ----------
# Common target names
# ----------
include_directories(${CMAKE_CURRENT_LIST_DIR}/ext/tinyobjloader)
set(SOURCE_FILES ${SOURCE_FILES};${CMAKE_CURRENT_LIST_DIR}/ext/tinyobjloader/tiny_obj_loader.h)
set(SOURCE_FILES ${SOURCE_FILES};${CMAKE_CURRENT_LIST_DIR}/ext/tinyobjloader/tiny_obj_loader.cc)
set(SOURCE_FILES ${SOURCE_FILES};${CMAKE_CURRENT_LIST_DIR}/ext/tinyply/source/tinyply.h)
set(SOURCE_FILES ${SOURCE_FILES};${CMAKE_CURRENT_LIST_DIR}/ext/tinyply/source/tinyply.cpp)

# ----------
# Components
# ----------
set(SOURCE_FILES ${SOURCE_FILES} "tinymesh.h")
add_folder(SOURCE_FILES "core")
add_folder(SOURCE_FILES "math")
add_folder(SOURCE_FILES "trimesh")
add_folder(SOURCE_FILES "geoproc")

if (BUILD_STATIC_LIB)
    add_library(${BUILD_TARGET} STATIC ${SOURCE_FILES})
else()
    add_library(${BUILD_TARGET} SHARED ${SOURCE_FILES})
endif()

target_include_directories(${BUILD_TARGET} PRIVATE ${CMAKE_CURRENT_LIST_DIR})

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi")
    set_property(TARGET ${BUILD_TARGET} APPEND PROPERTY LINK_FLAGS "/DEBUG /PROFILE")
endif()
