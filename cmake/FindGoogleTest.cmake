include(FindPackageHandleStandardArgs)

set(GoogleTest_DIR "GoogleTest_DIR" CACHE PATH "")

if (WIN32)
    find_path(GTEST_INCLUDE_DIR
              NAMES gtest/gtest.h
              PATHS
              ${GoogleTest_DIR}
              ${GoogleTest_DIR}/include)

    find_library(GTEST_LIBRARY
                 NAMES gtest
                 PATHS
                 ${GoogleTest_DIR}
                 ${GoogleTest_DIR}/lib
                 ${GoogleTest_DIR}/lib/Release)

    find_library(GTEST_MAIN_LIBRARY
                 NAMES gtest_main
                 PATHS
                 ${GoogleTest_DIR}
                 ${GoogleTest_DIR}/lib
                 ${GoogleTest_DIR}/lib/Release)
else()
    find_path(GTEST_INCLUDE_DIR
              NAMES gtest/gtest.h
              PATHS
              /usr/include
              /usr/local/include
              ${GoogleTest_DIR}
              ${GoogleTest_DIR}/include)

    find_library(GTEST_LIBRARY
                 NAMES gtest
                 PATHS
                 /usr/lib
                 /usr/local/lib
                 /usr/lib/x86_64-linux-gnu
                 ${GoogleTest_DIR}
                 ${GoogleTest_DIR}/lib)

    find_library(GTEST_MAIN_LIBRARY
                 NAMES gtest_main
                 PATHS
                 /usr/lib
                 /usr/local/lib
                 /usr/lib/x86_64-linux-gnu
                 ${GoogleTest_DIR}
                 ${GoogleTest_DIR}/lib
                 ${GoogleTest_DIR}/lib/Release)
endif()

find_package_handle_standard_args(
    GoogleTest DEFAULT_MSG
    GTEST_INCLUDE_DIR
    GTEST_LIBRARY
    GTEST_MAIN_LIBRARY
)

if (GoogleTest_FOUND)
    message(STATUS "GoogleTest include: ${GTEST_INCLUDE_DIR}")
    message(STATUS "GoogleTest library: ${GTEST_LIBRARY}")
    message(STATUS "GoogleTest main library: ${GTEST_MAIN_LIBRARY}")
    set(GTEST_INCLUDE_DIRS ${GTEST_INCLUDE_DIR})
    set(GTEST_LIBRARIES ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
endif()

mark_as_advanced(GTEST_INCLUDE_DIR GTEST_LIBRARY GTEST_MAIN_LIBRARY)
