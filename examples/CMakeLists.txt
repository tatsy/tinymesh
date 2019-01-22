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

add_example(hello_tinymesh)
add_example(simplify)