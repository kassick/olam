cmake_minimum_required (VERSION 2.6)
project (olamfs-trace-converter)

find_library(RASTRO_LIBRARY NAMES Rastro PATHS ENV RASTRO_PREFIX PATH_SUFFIXES lib/)
find_path(RASTRO_INCLUDES NAMES rastro.h PATHS ENV RASTRO_PREFIX
  PATH_SUFFIXES include/)
add_subdirectory(src/)
