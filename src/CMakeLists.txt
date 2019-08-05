# find Eigen library
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories( ${EIGEN3_INCLUDE_DIRS} )
find_package (pybind11 2.0 REQUIRED NO_MODULE)

# set c++11 rnd features if available
if ( HAVE_CXX11RANDOM )
  add_definitions( -DCXX11RANDOM )
endif()

# library version
set(LFClib_VERSION_MAJOR 0)
set(LFClib_VERSION_MINOR 2)
set(LFClib_VERSION_BUGFIX 0)

pybind11_add_module(lfclib pybinds.cc lattice.cc dipolesum.cc pile.cc)
