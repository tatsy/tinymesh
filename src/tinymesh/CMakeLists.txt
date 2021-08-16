# ----------
# Common target names
# ----------
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_LIST_DIR}/ext/tinyobjloader)
set(SOURCE_FILES ${SOURCE_FILES};${CMAKE_CURRENT_LIST_DIR}/ext/tinyobjloader/tiny_obj_loader.h)
include_directories(${CMAKE_CURRENT_LIST_DIR}/ext/tinyply/source)
set(SOURCE_FILES ${SOURCE_FILES};${CMAKE_CURRENT_LIST_DIR}/ext/tinyply/source/tinyply.h)

# ----------
# Components
# ----------
set(SOURCE_FILES ${SOURCE_FILES} "tinymesh.h")
add_folder(SOURCE_FILES "core")
add_folder(SOURCE_FILES "filters")
add_folder(SOURCE_FILES "remesh")

add_library(${TINYMESH_LIBRARY} STATIC ${SOURCE_FILES})

source_group("Source Files" FILES ${SOURCE_FILES})

target_include_directories(${TINYMESH_LIBRARY} PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(${TINYMESH_LIBRARY} ${CXX_FS_LIBRARY})

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi")
    set_property(TARGET ${TINYMESH_LIBRARY} APPEND PROPERTY LINK_FLAGS "/DEBUG /PROFILE")
endif()
