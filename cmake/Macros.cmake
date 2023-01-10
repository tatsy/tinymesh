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
# Add example program
# ----------
function(ADD_EXAMPLE)
    set(options)
    set(oneValueArgs NAME)
    set(multiValueArgs SOURCES)
    cmake_parse_arguments(ADD_EXAMPLE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    message(STATUS "Example: ${ADD_EXAMPLE_NAME}")
    set(EXP_EXE_NAME "example_${ADD_EXAMPLE_NAME}")
    add_executable(${EXP_EXE_NAME} ${ADD_EXAMPLE_SOURCES})
    add_dependencies(${EXP_EXE_NAME} ${TINYMESH_LIBRARY})

    target_include_directories(${EXP_EXE_NAME} PRIVATE ${TINYMESH_INCLUDE_DIR})
    target_link_libraries(${EXP_EXE_NAME} PRIVATE ${TINYMESH_LIBRARY})

    set_target_properties(${EXP_EXE_NAME} PROPERTIES FOLDER "Examples")
    set_target_properties(${EXP_EXE_NAME} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
    source_group("Source Files" FILES ${SOURCE_FILES})

    if (MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi")
        set_property(TARGET ${EXP_EXE_NAME} APPEND PROPERTY LINK_FLAGS "/DEBUG /PROFILE")
    endif()
endfunction(ADD_EXAMPLE)

# ----------
# Add test program
# ----------
function(ADD_UNIT_TEST)
    set(options)
    set(oneValueArgs NAME DEPS)
    set(multiValueArgs SOURCES)
    cmake_parse_arguments(ADD_UNIT_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    message(STATUS "Test: ${ADD_UNIT_TEST_NAME}")
    add_executable(${ADD_UNIT_TEST_NAME})
    add_dependencies(${ADD_UNIT_TEST_NAME} ${TINYMESH_LIBRARY})

    target_sources(${ADD_UNIT_TEST_NAME} PRIVATE ${ADD_UNIT_TEST_SOURCES})
    target_include_directories(${ADD_UNIT_TEST_NAME} PRIVATE ${TINYMESH_INCLUDE_DIR} ${GTEST_INCLUDE_DIRS})
    target_link_libraries(${ADD_UNIT_TEST_NAME} PRIVATE ${TINYMESH_LIBRARY} ${GTEST_LIBRARIES})

    set_target_properties(${ADD_UNIT_TEST_NAME} PROPERTIES FOLDER "Tests")
    set_target_properties(${ADD_UNIT_TEST_NAME} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
    source_group("Source Files" FILES ${SOURCE_FILES})

    set(${ADD_UNIT_TEST_DEPS} ${${ADD_UNIT_TEST_DEPS}} ${ADD_UNIT_TEST_NAME} PARENT_SCOPE)

    add_test(NAME ${ADD_UNIT_TEST_NAME} COMMAND ${ADD_UNIT_TEST_NAME})
    set_tests_properties(${ADD_UNIT_TEST_NAME} PROPERTIES TIMEOUT 120)
endfunction(ADD_UNIT_TEST)
