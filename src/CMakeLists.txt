# set source files
set (sources simplesum.c rotatesum.c fastincommsum.c dipolartensor.c pile.c mat3.c vec3.c)
set (devel-headers simplesum.h fastincommsum.h rotatesum.h dipolartensor.h config.h)

FIND_PACKAGE(GSL)
if (GSL_FOUND)
  add_definitions(-D__GSL)
  SET(LIBS ${LIBS} ${GSL_LIBRARIES})
  INCLUDE_DIRECTORIES(${INCLUDE_DIRS})
else()
  message("GSL not found, using suboptimal portable version.")
endif()

#show problems
if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
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

