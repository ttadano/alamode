cmake_minimum_required(VERSION 3.1)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -O2")
#set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_BUILD_TYPE RELEASE)

project(alm)

include_directories("/Users/tadano/src/spglib/include")
set(spglib "-L/Users/tadano/src/spglib/lib -lsymspg")

#
# Add openmp flag if available
#
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
find_package(LAPACK REQUIRED)
include_directories(${Lapack_INCLUDE_DIRS})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Version numbers
file(READ ${PROJECT_SOURCE_DIR}/../include/version.h version_file)
string(REGEX MATCH "ALAMODE_VERSION = \"([0-9]+\\.[0-9]+\\.[0-9]+)\"" alm_version ${version_file})
set(alm_version ${CMAKE_MATCH_1})
MESSAGE("ALM version: ${CMAKE_MATCH_1}")
set(serial "${alm_version}")
set(soserial "1")

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
set_property(TARGET alm PROPERTY CXX_STANDARD 11)
set_property(TARGET alm PROPERTY CXX_STANDARD_REQUIRED ON)

# Shared library
add_library(almcxx SHARED ${SOURCES})
target_link_libraries(almcxx ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${spglib})
set_property(TARGET almcxx PROPERTY VERSION ${serial})
set_property(TARGET almcxx PROPERTY SOVERSION ${soserial})
set_property(TARGET almcxx PROPERTY CXX_STANDARD 11)
set_property(TARGET almcxx PROPERTY CXX_STANDARD_REQUIRED ON)
install(TARGETS almcxx DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# Static link library
add_library(almcxx_static STATIC ${SOURCES})
set_property(TARGET almcxx_static PROPERTY VERSION ${serial})
set_property(TARGET almcxx_static PROPERTY SOVERSION ${soserial})
set_property(TARGET almcxx_static PROPERTY OUTPUT_NAME almcxx)
set_property(TARGET almcxx_static PROPERTY CXX_STANDARD 11)
set_property(TARGET almcxx_static PROPERTY CXX_STANDARD_REQUIRED ON)
install(TARGETS almcxx_static ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# Header file
install(FILES ${PROJECT_SOURCE_DIR}/alm.h DESTINATION ${PROJECT_SOURCE_DIR}/include)
