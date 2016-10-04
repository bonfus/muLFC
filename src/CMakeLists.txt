# set source files
set (sources simplesum.c fastincommsum.c pile.c rotatesum.c utility.c dipolartensor.c mat3.c vec3.c)
set (devel-headers simplesum.h fastincommsum.h rotatesum.h dipolartensor.h config.h)


# library version
set(LFClib_VERSION_MAJOR 0)
set(LFClib_VERSION_MINOR 1)
set(LFClib_VERSION_BUGFIX 0)

# dynamic library
add_library (LFClib SHARED ${sources})
set_target_properties(LFClib PROPERTIES
	VERSION ${LFClib_VERSION_MAJOR}.${LFClib_VERSION_MINOR}.${LFClib_VERSION_BUGFIX}
	SOVERSION ${LFClib_VERSION_MAJOR} )

# static library
add_library (LFClib-static STATIC  ${sources})
set_target_properties(LFClib-static PROPERTIES OUTPUT_NAME "lfc")

# install libraries
install(TARGETS LFClib DESTINATION lib)

# install headers
install(FILES ${devel-headers} DESTINATION include/liblfc)

