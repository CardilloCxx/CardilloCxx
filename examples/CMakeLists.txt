add_executable(example main.cpp)
target_link_libraries(example PRIVATE cardillo)

add_executable(free_flight free_flight.cpp)
target_link_libraries(free_flight PRIVATE cardillo)

add_executable(pm_midpoint pm_midpoint.cpp)
target_link_libraries(pm_midpoint PRIVATE cardillo)