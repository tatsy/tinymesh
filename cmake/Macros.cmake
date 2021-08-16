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
    target_link_libraries(${EXP_EXE_NAME} PRIVATE ${TINYMESH_LIBRARY} ${CXX_FS_LIBRARY})

    set_target_properties(${EXP_EXE_NAME} PROPERTIES FOLDER "Examples")
    set_target_properties(${EXP_EXE_NAME} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
    source_group("Source Files" FILES ${SOURCE_FILES})

    if (MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi")
        set_property(TARGET ${EXP_EXE_NAME} APPEND PROPERTY LINK_FLAGS "/DEBUG /PROFILE")
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
    add_definitions(-DGTEST_LANG_CXX11)
    add_executable(${FOLDER_NAME} ${SOURCE_FILES})
    target_link_libraries(${FOLDER_NAME} ${TINYMESH_LIBRARY} ${GTEST_LIBRARIES})
    target_include_directories(${FOLDER_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR})

    set_target_properties(${FOLDER_NAME} PROPERTIES FOLDER "Tests")
    set_target_properties(${FOLDER_NAME} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
    source_group("Source Files" FILES ${SOURCE_FILES})

    add_test(NAME ${FOLDER_NAME} COMMAND ${FOLDER_NAME})
endfunction(add_test_folder)
