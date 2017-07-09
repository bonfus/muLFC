# set source files
set (sources simplesum.c pile.c mat3.c vec3.c)
set (devel-headers simplesum.h fastincommsum.h rotatesum.h dipolartensor.h config.h)

FIND_PACKAGE(GSL)
if (GSL_FOUND)
  add_definitions(-D__GSL)
  SET(LIBS ${LIBS} ${GSL_LIBRARIES})
  INCLUDE_DIRECTORIES(${INCLUDE_DIRS})
else()
  message("GSL not found, using suboptimal portable version.")
endif()


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

