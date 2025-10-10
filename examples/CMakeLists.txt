add_executable(example main.cpp)
target_link_libraries(example PRIVATE cardillo)

add_executable(free_flight free_flight.cpp)
target_link_libraries(free_flight PRIVATE cardillo)
