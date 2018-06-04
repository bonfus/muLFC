# find Eigen library
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable (Test test.cc)

target_include_directories ( Test PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/eigen3)
target_link_libraries (Test LFClib)
