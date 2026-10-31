# Install rules and the CMake package: after `cmake --install`, other projects use
#   find_package(cardillo REQUIRED)
#   target_link_libraries(app PRIVATE cardillo::cardillo)
# Fetched dependencies (Eigen, EnTT, Coal, conic backends) install themselves into the same prefix.
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CARDILLO_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/cardillo")

install(TARGETS cardillo
    EXPORT cardilloTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT cardilloTargets NAMESPACE cardillo:: DESTINATION ${CARDILLO_CMAKE_DIR})

# Clarabel's Rust core is a plain shared library built by cargo, without install rules of its own.
if(CARDILLO_WITH_CLARABEL)
    install(FILES $<TARGET_PROPERTY:libclarabel_c_shared,INTERFACE_LINK_LIBRARIES> DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cardilloConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/cardilloConfig.cmake
    INSTALL_DESTINATION ${CARDILLO_CMAKE_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/cardilloConfigVersion.cmake COMPATIBILITY SameMinorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cardilloConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/cardilloConfigVersion.cmake
    DESTINATION ${CARDILLO_CMAKE_DIR})
