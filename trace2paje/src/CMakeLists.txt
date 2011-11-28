cmake_minimum_required (VERSION 2.6)
include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${RASTRO_INCLUDES}" )




#add_library(Rastro SHARED IMPORTED)






BISON_TARGET(DescriptionFileParser descriptionfile.y ${CMAKE_CURRENT_BINARY_DIR}/descriptionfileparser.cpp)
FLEX_TARGET(DescriptionFileScanner descriptionfile.l ${CMAKE_CURRENT_BINARY_DIR}/descriptionfilescanner.cpp)
ADD_FLEX_BISON_DEPENDENCY(DescriptionFileScanner DescriptionFileParser)

include_directories(tree)
include_directories(${CMAKE_CURRENT_BINARY_DIR})



############################
## Main exe

add_executable(rastro2paje
  baseevent.cc
  baseeventtype.cc
  statetype.cc
  linktype.cc
  link.cc
  dummyevent.cc
  state.cc
  event.cc
  containertrigger.cc
  container.cc
  semantics.cc
  paje.cc
  paje_functions.cc
  pajeelement.cc
  symbols.cc
  rastro_loop.cc
  rastro_generate.c
  usermap.cc
  tree/tree.hpp # my tree implementation
  rastro2paje.cc
  logging.cc
  ${BISON_DescriptionFileParser_OUTPUTS}
  ${FLEX_DescriptionFileScanner_OUTPUTS}
)
set_property(TARGET rastro2paje PROPERTY COMPILE_FLAGS "-std=c++0x")
set_property(TARGET rastro2paje PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
#set_property(TARGET rastro2paje PROPERTY BUILD_WITH_INSTALL_RPATH true)
#set_property(TARGET rastro2paje PROPERTY IMPORTED_LOCATION "${RASTRO_PREFIX}/lib")

set_property(SOURCE 
  container.cc semantics.cc paje.cc rastro2paje.cc symbols.cc
  PROPERTY COMPILE_FLAGS "-std=c++0x")

target_link_libraries(rastro2paje fl boost_regex boost_filesystem boost_system ${RASTRO_LIBRARIES})




install (TARGETS rastro2paje DESTINATION bin)



#####################################
## Symbol test
add_executable(symbols_test
  symbols.cc
  symbols.hh
  )

target_link_libraries(symbols_test boost_regex)
set_property(TARGET symbols_test PROPERTY COMPILE_DEFINITIONS SYMBOL_UNIT_TEST=1)



###############################
## Scanner Test
#add_executable(scanner_test 
#  baseevent.cc
#  baseeventtype.cc
#  linktype.cc
#  statetype.cc
#  link.cc
#  dummyevent.cc
#  state.cc
#  containertrigger.cc
#  symbols.cc
#  event.cc
#  container.cc
#  semantics.cc
#  paje.cc
#  paje_functions.cc
#  pajeelement.cc
#  usermap.cc
#  logging.cc
#  tree/tree.hpp # my tree implementation
#  ${FLEX_DescriptionFileScanner_OUTPUTS} 
#  ${BISON_DescriptionFileParser_OUTPUTS}
#)

#set_property(TARGET scanner_test PROPERTY COMPILE_FLAGS "-std=c++0x")
#set_property(TARGET scanner_test PROPERTY COMPILE_DEFINITIONS SCANNER_UNIT_TEST=1)
#target_link_libraries(scanner_test fl)



############################
### Tree test
add_executable(tree/tree_test
  tree/tree.cpp
  tree/tree.hpp
)
set_property(TARGET tree/tree_test PROPERTY COMPILE_FLAGS "-std=c++0x")

#FLEX_TARGET(test_tokens descriptionfile.l descriptionfile.yy.c) #[COMPILE_FLAGS <string>])

#add_executable(test_tokens descriptionfile.l)
#set_property(TARGET test_tokens PROPERTY COMPILE_DEFINITIONS UNIT_TEST=1)

#add_executable(olam_trace2paje olam_trace2paje.c aky_paje.c aky_names.c aky_keys.c)
#target_link_libraries(olam_trace2paje rastro)

#add_executable(aky_names aky_names.c)
#set_property(TARGET aky_names PROPERTY COMPILE_DEFINITIONS UNIT_TEST=1)

#add_executable(gen_fortran_names gen_fortran_names.c aky_names.c)


#install (TARGETS olam_trace2paje DESTINATION bin)
