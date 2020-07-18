cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -O0 -D_DEBUG")

project(alm)

# Version numbers
file(READ ${PROJECT_SOURCE_DIR}/../include/version.h version_file)
string(REGEX MATCH "ALAMODE_VERSION = \"([0-9]+\\.[0-9]+\\.[0-9]+)\"" alm_version ${version_file})
set(alm_version ${CMAKE_MATCH_1})
MESSAGE("Setting up for ALM code (ALAMODE version: ${CMAKE_MATCH_1})")

if (SPGLIB_PATH)
    include_directories("${SPGLIB_PATH}/include")
    set(spglib "-L${SPGLIB_PATH}/lib -lsymspg")
else()
    message(FATAL_ERROR "SPGLIB_PATH must be set via -DSPGLIB_PATH option.")
endif()
#
# Add openmp flag if available
#
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Boost)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    if (BOOST_INCLUDE)
        include_directories(${BOOST_INCLUDE})
    else()
        message(FATAL_ERROR "Boost was not found.
        Please specify the location of boost include directories
        via -DBOOST_INCLUDE option.")
    endif()
endif()

find_package(Eigen3)
if (Eigen3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    if (EIGEN3_INCLUDE)
        include_directories(${EIGEN3_INCLUDE})
    else()
        message(FATAL_ERROR "Eigen3 was not found.
        Please specify the location of boost include directories
        via -DEIGEN3_INCLUDE option.")
    endif()
endif()

find_package(LAPACK REQUIRED)
include_directories(${Lapack_INCLUDE_DIRS})

# Source code
include_directories("${PROJECT_SOURCE_DIR}/../include")
set(SOURCES ${PROJECT_SOURCE_DIR}/alm.cpp
            ${PROJECT_SOURCE_DIR}/cluster.cpp
            ${PROJECT_SOURCE_DIR}/constraint.cpp
            ${PROJECT_SOURCE_DIR}/fcs.cpp
            ${PROJECT_SOURCE_DIR}/files.cpp
            ${PROJECT_SOURCE_DIR}/input_parser.cpp
            ${PROJECT_SOURCE_DIR}/input_setter.cpp
            ${PROJECT_SOURCE_DIR}/optimize.cpp
            ${PROJECT_SOURCE_DIR}/patterndisp.cpp
            ${PROJECT_SOURCE_DIR}/rref.cpp
            ${PROJECT_SOURCE_DIR}/symmetry.cpp
            ${PROJECT_SOURCE_DIR}/system.cpp
            ${PROJECT_SOURCE_DIR}/timer.cpp
            ${PROJECT_SOURCE_DIR}/writer.cpp)

set(HEADERS ${PROJECT_SOURCE_DIR}/alm.h
            ${PROJECT_SOURCE_DIR}/cluster.h
            ${PROJECT_SOURCE_DIR}/constraint.h
            ${PROJECT_SOURCE_DIR}/fcs.h
            ${PROJECT_SOURCE_DIR}/files.h
            ${PROJECT_SOURCE_DIR}/input_parser.h
            ${PROJECT_SOURCE_DIR}/input_setter.h
            ${PROJECT_SOURCE_DIR}/optimize.h
            ${PROJECT_SOURCE_DIR}/patterndisp.h
            ${PROJECT_SOURCE_DIR}/rref.h
            ${PROJECT_SOURCE_DIR}/symmetry.h
            ${PROJECT_SOURCE_DIR}/system.h
            ${PROJECT_SOURCE_DIR}/timer.h
            ${PROJECT_SOURCE_DIR}/writer.h)

# # Trick to use gcc to compile *.cpp: This doesn't work because of boost
# SET_SOURCE_FILES_PROPERTIES(${SOURCES} PROPERTIES LANGUAGE C)
# set(CMAKE_C_COMPILER "gcc")

# Executable
add_executable(alm ${PROJECT_SOURCE_DIR}/main.cpp
                   ${PROJECT_SOURCE_DIR}/alm_cui.cpp
                   ${SOURCES} ${HEADERS})
target_link_libraries(alm ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${spglib})
target_link_libraries(alm ${LAPACK_LIBRARIES})

if (BUILD_LIBRARIES)
    set(serial "${alm_version}")
    set(soserial "1")
    # Shared library
    add_library(almcxx SHARED ${SOURCES})
    target_link_libraries(almcxx ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${spglib})
    set_property(TARGET almcxx PROPERTY VERSION ${serial})
    set_property(TARGET almcxx PROPERTY SOVERSION ${soserial})
    install(TARGETS almcxx DESTINATION ${PROJECT_SOURCE_DIR}/lib)

    # Static link library
    add_library(almcxx_static STATIC ${SOURCES})
    set_property(TARGET almcxx_static PROPERTY VERSION ${serial})
    set_property(TARGET almcxx_static PROPERTY SOVERSION ${soserial})
    set_property(TARGET almcxx_static PROPERTY OUTPUT_NAME almcxx)
    install(TARGETS almcxx_static ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib)

    # Header file
    install(FILES ${PROJECT_SOURCE_DIR}/alm.h DESTINATION ${PROJECT_SOURCE_DIR}/include)
endif()
