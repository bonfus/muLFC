# -*- mode: cmake -*-
include(ExternalProject)
include(CMakePushCheckState)

# Check for Eigen
find_package(Eigen3 3.0)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if(EIGEN3_FOUND)

  message(STATUS "Using system version of Eigen library")

  cmake_push_check_state()

  # Perform a compile check with Eigen
  list(APPEND CMAKE_REQUIRED_INCLUDES ${EIGEN3_INCLUDE_DIR})
  # compile and see...
  CHECK_CXX_SOURCE_COMPILES("
    #include <Eigen/Core>
    #include <Eigen/Dense>
    #include <iostream>
    int main(int argc, char* argv[]){
      Eigen::MatrixXd m = Eigen::MatrixXd::Random(5, 5);
      m = m.transpose() + m;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(m);
      Eigen::MatrixXd m_invsqrt = eig.operatorInverseSqrt();
      std::cout << m_invsqrt << std::endl;
    }"
    EIGEN_COMPILES)
      
  cmake_pop_check_state()

  if (NOT EIGEN_COMPILES)
    message(FATAL_ERROR "Eigen found at ${EIGEN3_INCLUDE_DIR}, but could not compile test program")
  endif()

else()

  # Create a cache entry for Eigen build variables.
  set(EIGEN_URL "https://github.com/eigenteam/eigen-git-mirror.git" CACHE STRING 
      "Path to the Eigen repository")
  set(EIGEN_TAG "3.2.4" CACHE STRING "The Eigen revision tag")

  message("** Will build Eigen from ${EIGEN_URL}")

  ExternalProject_Add(eigen3
        GIT_REPOSITORY ${EIGEN_URL}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/eigen
        GIT_SHALLOW 1)

  set(EIGEN3_INCLUDE_DIR ${PROJECT_BINARY_DIR}/eigen/include/eigen3)
  
endif()
