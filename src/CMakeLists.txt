# set source files
set (sources simplesum.cc lattice.cc pile.cc dipolesum.cc fastincommsum.cc rotatesum.cc dipolartensor.cc)
set (devel-headers simplesum.h fastincommsum.h rotatesum.h dipolartensor.h config.h)

#Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# library version
set(LFClib_VERSION_MAJOR 0)
set(LFClib_VERSION_MINOR 2)
set(LFClib_VERSION_BUGFIX 0)

# dynamic library
add_library (LFClib SHARED ${sources})
target_link_libraries (LFClib Eigen3::Eigen)

set_target_properties(LFClib PROPERTIES
	VERSION ${LFClib_VERSION_MAJOR}.${LFClib_VERSION_MINOR}.${LFClib_VERSION_BUGFIX}
	SOVERSION ${LFClib_VERSION_MAJOR} )

# static library
add_library (LFClib-static STATIC  ${sources})
target_link_libraries (LFClib-static Eigen3::Eigen)
set_target_properties(LFClib-static PROPERTIES OUTPUT_NAME "lfc")

# install libraries
install(TARGETS LFClib DESTINATION lib)

# install headers
install(FILES ${devel-headers} DESTINATION include/liblfc)

