cmake_minimum_required (VERSION 2.6)
include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${RASTRO_INCLUDES}" )

add_library(rastro SHARED IMPORTED)
set_property(TARGET rastro PROPERTY IMPORTED_LOCATION "${RASTRO_LIBRARY}")

add_executable(olam_trace2paje olam_trace2paje.c aky_paje.c aky_names.c aky_keys.c)
target_link_libraries(olam_trace2paje rastro)

add_executable(aky_names aky_names.c)
set_property(TARGET aky_names PROPERTY COMPILE_DEFINITIONS UNIT_TEST=1)


install (TARGETS olam_trace2paje DESTINATION bin)
