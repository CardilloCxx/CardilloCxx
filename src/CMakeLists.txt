# The engine library. Users include its headers as <cardillo/...>. Sources of the optional conic
# backends (cardillo/physics/solver/conic/) are added only when their CARDILLO_WITH_<NAME> option is on;
# their headers are private to the library.
set(CARDILLO_CONIC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cardillo/physics/solver/conic")

file(GLOB_RECURSE CARDILLO_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/cardillo/*.cpp")
file(GLOB_RECURSE CARDILLO_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/cardillo/*.hpp")
list(FILTER CARDILLO_SOURCES EXCLUDE REGEX "/physics/solver/conic/")
list(FILTER CARDILLO_HEADERS EXCLUDE REGEX "/physics/solver/conic/")

if(CARDILLO_WITH_CLARABEL OR CARDILLO_WITH_SCS)
    # Cone program shared by Clarabel and SCS.
    list(APPEND CARDILLO_SOURCES "${CARDILLO_CONIC_DIR}/cone_program_assembler.cpp")
endif()
foreach(_backend CLARABEL QOCO CONICXX SCS MOSEK)
    if(CARDILLO_WITH_${_backend})
        string(TOLOWER ${_backend} _dir)
        file(GLOB _backend_sources CONFIGURE_DEPENDS "${CARDILLO_CONIC_DIR}/${_dir}/*.cpp")
        list(APPEND CARDILLO_SOURCES ${_backend_sources})
    endif()
endforeach()

add_library(cardillo SHARED ${CARDILLO_SOURCES})
add_library(cardillo::cardillo ALIAS cardillo)
target_sources(cardillo PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${CARDILLO_HEADERS})
target_compile_features(cardillo PUBLIC cxx_std_20)
if(CARDILLO_HAS_MARCH_NATIVE AND CARDILLO_ENABLE_NATIVE_ARCH)
    # Code using cardillo must be compiled for the same ISA: the alignment of fixed-size Eigen members
    # (e.g. 32 bytes with AVX, 16 with SSE2) is part of the ABI of the headers.
    target_compile_options(cardillo PUBLIC -march=native)
endif()
set_target_properties(cardillo PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    # The installed library finds the (fetched and installed) shared dependencies next to it.
    INSTALL_RPATH "$ORIGIN")

# Eigen, EnTT and Coal types appear in the public headers; VTK is used only in the implementation.
target_link_libraries(cardillo
    PUBLIC Eigen3::Eigen EnTT::EnTT coal::coal
    PRIVATE VTK::CommonCore VTK::CommonDataModel VTK::IOXML)

if(OpenMP_CXX_FOUND)
    target_link_libraries(cardillo PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(cardillo PUBLIC CARDILLO_HAVE_OPENMP=1)
else()
    message(STATUS "OpenMP not found; the condensed solver's parallel sweep modes run single-threaded.")
endif()

# The backends are linked as build-only dependencies: they are not part of the exported interface,
# the installed libcardillo finds their shared libraries next to itself (INSTALL_RPATH).
if(CARDILLO_WITH_CLARABEL)
    target_link_libraries(cardillo PRIVATE $<BUILD_INTERFACE:libclarabel_c_shared>)
    target_compile_definitions(cardillo PUBLIC CARDILLO_HAVE_CLARABEL=1)
endif()
if(CARDILLO_WITH_QOCO)
    # QOCO installs its library without a CMake package, so it is linked by file (its target would have
    # to be exported with cardillo). qoco_assembler.hpp fills QOCO's structs from its internal headers.
    add_dependencies(cardillo qoco)
    target_link_libraries(cardillo PRIVATE $<BUILD_INTERFACE:$<TARGET_FILE:qoco>>)
    get_property(_qoco_include GLOBAL PROPERTY QOCO_INCLUDE)
    target_include_directories(cardillo SYSTEM PRIVATE ${_qoco_include})
    target_compile_definitions(cardillo PUBLIC CARDILLO_HAVE_QOCO=1)
endif()
if(CARDILLO_WITH_CONICXX)
    target_link_libraries(cardillo PRIVATE $<BUILD_INTERFACE:conicxx::conicxx>)
    target_compile_definitions(cardillo PUBLIC CARDILLO_HAVE_CONICXX=1)
endif()
if(CARDILLO_WITH_SCS)
    target_link_libraries(cardillo PRIVATE $<BUILD_INTERFACE:scs::scsdir>)
    target_compile_definitions(cardillo PUBLIC CARDILLO_HAVE_SCS=1)
endif()
if(CARDILLO_WITH_MOSEK)
    # Linked by file like QOCO: MOSEK is a user installation, not part of cardillo's package. The
    # installed library finds it through its RUNPATH.
    target_link_libraries(cardillo PRIVATE $<BUILD_INTERFACE:${MOSEK_LIBRARY}>)
    target_include_directories(cardillo SYSTEM PRIVATE ${MOSEK_INCLUDE_DIR})
    get_filename_component(_mosek_lib_dir "${MOSEK_LIBRARY}" DIRECTORY)
    set_property(TARGET cardillo APPEND PROPERTY INSTALL_RPATH "${_mosek_lib_dir}")
    target_compile_definitions(cardillo PUBLIC CARDILLO_HAVE_MOSEK=1)
endif()
