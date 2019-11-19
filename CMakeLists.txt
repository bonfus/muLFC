cmake_minimum_required(VERSION 3.2)

project (LFC)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)

include(cxx)

OPTION (USE_OpenMP "Use OpenMP" OFF)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()

# find Eigen library
message(STATUS "Trying to include Eigen library")
include(eigen)
include_directories(${EIGEN3_INCLUDE_DIR})

# find pybind11
add_subdirectory(third-party/pybind11)
#find_package (pybind11 2.0 REQUIRED NO_MODULE)

add_subdirectory (src)

# need eigen before we start compiling the module
if (NOT EIGEN3_FOUND)
add_dependencies(lfclib eigen3)
endif()
