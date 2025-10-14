# add_executable(example main.cpp)
# target_link_libraries(example PRIVATE cardillo)
# 
# add_executable(free_flight free_flight.cpp)
# target_link_libraries(free_flight PRIVATE cardillo)

# add_executable(pm_midpoint pm_midpoint.cpp)
# target_link_libraries(pm_midpoint PRIVATE cardillo)

# Collisions demo writing scene and contacts VTK
add_executable(collisions_demo collisions_demo.cpp)
target_link_libraries(collisions_demo PRIVATE cardillo)

# ECS physics engine example (requires EnTT)
# add_executable(physics_engine physics_engine.cpp)
# target_link_libraries(physics_engine PRIVATE cardillo)
# target_include_directories(physics_engine PRIVATE ${entt_SOURCE_DIR}/src)