# ----------
# GoogleTest
# ----------
if(TINYMESH_BUILD_TESTS)
  enable_testing()

  find_package(GTest QUIET)

  if(NOT GTest_FOUND)
    include(FetchContent)

    if(MSVC)
      set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    endif()

    FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    )
    FetchContent_MakeAvailable(googletest)
  endif()
endif()
