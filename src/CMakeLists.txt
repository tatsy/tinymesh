add_subdirectory(tinymesh)

if(TINYMESH_BUILD_PYTHON)
  pybind11_add_module(tinymesh_python MODULE pybind11.cpp)

  target_link_libraries(tinymesh_python
    PRIVATE ${TINYMESH_LIBRARY}
  )

  set_target_properties(tinymesh_python PROPERTIES
    OUTPUT_NAME _core
    ARCHIVE_OUTPUT_NAME tinymesh_python
    DEBUG_POSTFIX ""
  )

  install(TARGETS tinymesh_python
    LIBRARY DESTINATION tinymesh
    RUNTIME DESTINATION tinymesh
  )
endif()
