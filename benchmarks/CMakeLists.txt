add_executable(benchmark_construct benchmark_construct.cpp)
target_link_libraries(benchmark_construct PRIVATE ${TINYMESH_LIBRARY})
target_compile_definitions(benchmark_construct PRIVATE TINYMESH_BENCH_BUILD_CONFIG="$<CONFIG>")
