configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test_config.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/test_config.h" @ONLY)

set(TEST_HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/test_config.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/test_utils.h"
)

set(TEST_DEPS)
add_unit_test(DEPS TEST_DEPS NAME "test_filters" SOURCES test_filters.cpp ${TEST_HEADERS})
add_unit_test(DEPS TEST_DEPS NAME "test_mesh" SOURCES test_mesh.cpp ${TEST_HEADERS})
add_unit_test(DEPS TEST_DEPS NAME "test_ops" SOURCES test_ops.cpp ${TEST_HEADERS})
add_unit_test(DEPS TEST_DEPS NAME "test_remesh" SOURCES test_remesh.cpp ${TEST_HEADERS})
add_unit_test(DEPS TEST_DEPS NAME "test_vec" SOURCES test_vec.cpp ${TEST_HEADERS})

add_custom_target(check
  COMMAND "${CMAKE_CTEST_COMMAND}"
    -C "$<CONFIG>"
    --output-on-failure
  DEPENDS ${TEST_DEPS}
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  VERBATIM
)

