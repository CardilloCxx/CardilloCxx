# ffclib export (src/cardillo/io/ffclib/): the one-step problems of a simulation written as ffclib files (fclib's
# global problem plus the compliance of the bilateral rows; docs/chapters/using_the_engine/ffclib_export.rst) with
# the ffclib library (Eigen and the HDF5 C library). An installed ffclib is used if found (find_package(ffclib),
# e.g. -DCMAKE_PREFIX_PATH=<prefix>), otherwise it is fetched from GitHub without its tests, examples and conic
# solvers. A local checkout can be used with -DFETCHCONTENT_SOURCE_DIR_FFCLIB=<path>.
include(FetchContent)

if(CARDILLO_WITH_FFCLIB)
    set(FFCLIB_BUILD_TESTS OFF)
    set(FFCLIB_BUILD_EXAMPLES OFF)
    set(FFCLIB_INSTALL OFF)
    foreach(_solver CLARABEL SCS QOCO CONICXX CUOPT)
        set(FFCLIB_WITH_${_solver} OFF)
    endforeach()
    FetchContent_Declare(ffclib GIT_REPOSITORY https://github.com/FlexibleFrictionalContactLibrary/ffclib.git GIT_TAG main
        FIND_PACKAGE_ARGS CONFIG)
    FetchContent_MakeAvailable(ffclib)
    if(ffclib_SOURCE_DIR)
        message(STATUS "Using the ffclib library from ${ffclib_SOURCE_DIR}")
    else()
        message(STATUS "Using the installed ffclib library (${ffclib_DIR})")
    endif()

    # Recorded in the description of every exported file; taken at configure time.
    find_package(Git QUIET)
    set(CARDILLO_GIT_DESCRIBE "unknown")
    if(GIT_FOUND)
        execute_process(COMMAND "${GIT_EXECUTABLE}" describe --always --dirty --abbrev=12
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" OUTPUT_VARIABLE _describe OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE _result ERROR_QUIET)
        if(_result EQUAL 0)
            set(CARDILLO_GIT_DESCRIBE "${_describe}")
        endif()
    endif()
endif()
