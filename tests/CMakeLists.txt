# Unit tests and solver benchmarks (configure with -DCARDILLO_BUILD_TESTS=ON). The MOSEK tests and
# benchmark additionally require -DCARDILLO_WITH_MOSEK=ON and a valid MOSEK license at run time; the
# SCS ones require -DCARDILLO_WITH_SCS=ON.

if(CARDILLO_WITH_MOSEK)
    add_executable(mosek_tests mosek_tests.cpp)
    target_link_libraries(mosek_tests PRIVATE cardillo)
    add_test(NAME mosek_tests COMMAND mosek_tests)

    add_executable(mosek_benchmark mosek_benchmark.cpp)
    target_link_libraries(mosek_benchmark PRIVATE cardillo)
endif()

if(CARDILLO_WITH_SCS)
    add_executable(scs_tests scs_tests.cpp)
    target_link_libraries(scs_tests PRIVATE cardillo)
    add_test(NAME scs_tests COMMAND scs_tests)

    add_executable(scs_benchmark scs_benchmark.cpp)
    target_link_libraries(scs_benchmark PRIVATE cardillo)
endif()

# Replays cone programs dumped by ScsSolver (scs.dump_dir) through all conic backends.
if(CARDILLO_WITH_SCS)
    add_executable(conic_replay conic_replay.cpp)
    target_link_libraries(conic_replay PRIVATE cardillo)
endif()
