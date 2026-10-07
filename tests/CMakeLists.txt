# Unit tests and solver benchmarks (configure with -DCARDILLO_BUILD_TESTS=ON). The MOSEK tests and
# benchmark additionally require -DCARDILLO_WITH_MOSEK=ON and a valid MOSEK license at run time.

if(CARDILLO_WITH_MOSEK)
    add_executable(mosek_tests mosek_tests.cpp)
    target_link_libraries(mosek_tests PRIVATE cardillo)
    add_test(NAME mosek_tests COMMAND mosek_tests)

    add_executable(mosek_benchmark mosek_benchmark.cpp)
    target_link_libraries(mosek_benchmark PRIVATE cardillo)
endif()
