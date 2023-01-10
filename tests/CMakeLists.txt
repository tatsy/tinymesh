set(BUILD_TARGET "run_all_tests")

configure_file("${CMAKE_CURRENT_LIST_DIR}/test_config.h.in"
               "${CMAKE_CURRENT_LIST_DIR}/test_config.h" @ONLY)

set(TEST_HEADERS test_config.h test_utils.h)

set(SOURCE_FILES
  test_config.h
  test_utils.h
  gtest.cpp
  test_filters.cpp
  test_mesh.cpp
  test_ops.cpp
  test_remesh.cpp
  test_vec.cpp)

set(TEST_DEPS)
add_unit_test(DEPS TEST_DEPS NAME "test_filters" SOURCES test_filters.cpp ${TEST_HEADERS})
add_unit_test(DEPS TEST_DEPS NAME "test_mesh" SOURCES test_mesh.cpp ${TEST_HEADERS})
add_unit_test(DEPS TEST_DEPS NAME "test_ops" SOURCES test_ops.cpp ${TEST_HEADERS})
add_unit_test(DEPS TEST_DEPS NAME "test_remesh" SOURCES test_remesh.cpp ${TEST_HEADERS})
add_unit_test(DEPS TEST_DEPS NAME "test_vec" SOURCES test_vec.cpp ${TEST_HEADERS})
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose --gtest_shuffle DEPENDS ${TEST_DEPS})


# include_directories(${TINYMESH_INCLUDE_DIR} ${GTEST_INCLUDE_DIRS})
# add_executable(${BUILD_TARGET})

# target_sources(${BUILD_TARGET} PRIVATE ${SOURCE_FILES})
# target_include_directories(${BUILD_TARGET} PRIVATE ${CMAKE_CURRENT_LIST_DIR})
# target_link_libraries(${BUILD_TARGET} PRIVATE ${TINYMESH_LIBRARY} ${GTEST_LIBRARIES})

# set_target_properties(${BUILD_TARGET} PROPERTIES FOLDER "Tests")
# set_target_properties(${BUILD_TARGET} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
# source_group("Source Files" FILES ${SOURCE_FILES})

# add_test(NAME ${BUILD_TARGET} COMMAND ${BUILD_TARGET})
# add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose --gtest_shuffle DEPENDS ${BUILD_TARGET})
