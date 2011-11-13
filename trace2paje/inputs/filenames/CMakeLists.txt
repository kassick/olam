
# Import rastro2paje location from the current project
GET_TARGET_PROPERTY(TRACE_CONVERTER rastro2paje LOCATION)

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/ids.txt
   COMMAND ${TRACE_CONVERTER} -i ${CMAKE_CURRENT_SOURCE_DIR}/filenames.txt --gen-ids ${CMAKE_CURRENT_SOURCE_DIR}/ids.txt
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/filenames.txt
   )

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/evt_defs.h ${CMAKE_CURRENT_SOURCE_DIR}/rst_functions.c ${CMAKE_CURRENT_SOURCE_DIR}/rst_functions.h
  COMMAND  ${TRACE_CONVERTER} -i ${CMAKE_CURRENT_SOURCE_DIR}/main.txt --c-header ${CMAKE_CURRENT_SOURCE_DIR}/evt_defs.h --rst-c-functions ${CMAKE_CURRENT_SOURCE_DIR}/rst_functions.c --rst-c-header ${CMAKE_CURRENT_SOURCE_DIR}/rst_functions.h
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/ids.txt
)






add_executable(filenames filenames.c rst_functions.c evt_defs.h)
include_directories ("${RASTRO_INCLUDES}" )
target_link_libraries(filenames ${RASTRO_LIBRARIES} )

add_custom_target(filenames-run
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/filenames
  DEPENDS filenames
)

add_custom_target(filenames-paje 
  COMMAND ${TRACE_CONVERTER} -i main.txt -o main.paje *.rst
  DEPENDS filenames-run
)


