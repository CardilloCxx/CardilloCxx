function(cardillo_add_test target source)
	add_executable(${target} ${source})
	target_link_libraries(${target} PRIVATE cardillo)
	add_test(NAME ${target} COMMAND ${target})
endfunction()

cardillo_add_test(test_joint_jacobian test_joint_jacobian.cpp)
