# ----------
# Add source codes in the folder
# ----------
function(add_folder SOURCES FOLDER_NAME)
    file(GLOB LOCAL_FILES
         "${FOLDER_NAME}/*.cpp" "${FOLDER_NAME}/*.hpp"
         "${FOLDER_NAME}/*.c" "${FOLDER_NAME}/*.h")
    set(${SOURCES} ${${SOURCES}} ${LOCAL_FILES} PARENT_SCOPE)

    source_group(${FOLDER_NAME} FILES ${LOCAL_FILES})
endfunction()

# ----------
# Add example module
# ----------
function(ADD_EXAMPLE EXPNAME)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs)
    cmake_parse_arguments(ADD_EXAMPLE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(TARGET_EXAMPLE_DIR "${CMAKE_CURRENT_LIST_DIR}/${EXPNAME}")

    file(GLOB SOURCE_FILES "${EXPNAME}/*.cpp" "${EXPNAME}/*.h")

    include_directories(${TINYMESH_INCLUDE_DIR})
    link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

    message(STATUS "Example: ${EXPNAME}")
    add_executable(${EXPNAME} ${SOURCE_FILES} ${SHADER_FILES})
    add_dependencies(${EXPNAME} ${TINYMESH_LIBRARY})

    target_link_libraries(${EXPNAME} ${TINYMESH_LIBRARY})

    set_target_properties(${EXPNAME} PROPERTIES FOLDER "Examples")
    set_target_properties(${EXPNAME} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
    source_group("Source Files" FILES ${SOURCE_FILES})

    if (MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi")
        set_property(TARGET ${EXPNAME} APPEND PROPERTY LINK_FLAGS "/DEBUG /PROFILE")
    endif()
endfunction(ADD_EXAMPLE)

# ----------
# Add test program folder
# ----------
function(add_test_folder FOLDER_NAME)
    file(GLOB SOURCE_FILES
         "${FOLDER_NAME}/*.cpp" "${FOLDER_NAME}/*.hpp"
         "${FOLDER_NAME}/*.c" "${FOLDER_NAME}/*.h")
    include_directories(${TINYMESH_INCLUDE_DIR} ${GTEST_INCLUDE_DIRS})
    add_executable(${FOLDER_NAME} ${SOURCE_FILES})
    target_link_libraries(${FOLDER_NAME} ${TINYMESH_LIBRARY} ${GTEST_LIBRARIES})

    set_target_properties(${FOLDER_NAME} PROPERTIES FOLDER "Tests")
    set_target_properties(${FOLDER_NAME} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
    source_group("Source Files" FILES ${SOURCE_FILES})

    add_test(NAME "testing" COMMAND ${FOLDER_NAME})
    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose --gtest_shuffle DEPENDS ${FOLDER_NAME})
endfunction(add_test_folder)
