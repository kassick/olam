cmake_minimum_required (VERSION 2.6)
project (trace2paje)

find_library(RASTRO_LIBRARY NAMES Rastro PATHS ENV RASTRO_PREFIX PATH_SUFFIXES lib/)
find_path(RASTRO_INCLUDES NAMES rastro.h PATHS ENV RASTRO_PREFIX PATH_SUFFIXES include/)
get_filename_component(rastro_libdir ${RASTRO_LIBRARY} PATH)




set(RASTRO_LIBDIR ${rastro_libdir} CACHE PATH "Rastro libdir" )
set(RASTRO_LIBRARIES ${RASTRO_LIBRARY} CACHE FILEPATH "Rastro library1")
link_directories(${rastro_libdir})

find_package(FLEX)
find_package(BISON)
find_package(Boost)

add_subdirectory(src/)
set_property(DIRECTORY src/ PROPERTY COMPILE_FLAGS -std=c++0x)

