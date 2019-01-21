set(BUILD_TARGET "tinymesh")

file(GLOB SOURCE_FILES "*.cpp" "*.h")
include_directories(${CMAKE_CURRENT_LIST_DIR}/ext/tinyobjloader)
set(${SOURCE_FILES} ${SOURCE_FILES};${CMAKE_CURRENT_LIST_DIR}/ext/tinyobjloader/tiny_obj_loader.h)
set(${SOURCE_FILES} ${SOURCE_FILES};${CMAKE_CURRENT_LIST_DIR}/ext/tinyobjloader/tiny_obj_loader.cc)

add_library(${BUILD_TARGET} ${SOURCE_FILES})

source_group("Source Files" FILES ${SOURCE_FILES})

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi")
    set_property(TARGET ${BUILD_TARGET} APPEND PROPERTY LINK_FLAGS "/DEBUG /PROFILE")
endif()
