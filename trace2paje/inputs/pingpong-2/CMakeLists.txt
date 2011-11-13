
# Import rastro2paje location from the current project
GET_TARGET_PROPERTY(TRACE_CONVERTER rastro2paje LOCATION)

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/ids.txt
   COMMAND ${TRACE_CONVERTER} -i ${CMAKE_CURRENT_SOURCE_DIR}/pingpong-link.txt --gen-ids ${CMAKE_CURRENT_SOURCE_DIR}/ids.txt
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/pingpong-link.txt
   )

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/evt_defs.h ${CMAKE_CURRENT_SOURCE_DIR}/rst_functions.c ${CMAKE_CURRENT_SOURCE_DIR}/rst_functions.h
  COMMAND  ${TRACE_CONVERTER} -i ${CMAKE_CURRENT_SOURCE_DIR}/pingpong.txt --c-header ${CMAKE_CURRENT_SOURCE_DIR}/evt_defs.h --rst-c-functions ${CMAKE_CURRENT_SOURCE_DIR}/rst_functions.c --rst-c-header ${CMAKE_CURRENT_SOURCE_DIR}/rst_functions.h
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/ids.txt
)






add_executable(mpi-pingpong mpi-pingpong.c rst_functions.c evt_defs.h)
include_directories(${MPI_C_INCLUDE_PATH})
include_directories ("${RASTRO_INCLUDES}" )
set_property(TARGET mpi-pingpong PROPERTY COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
set_property(TARGET mpi-pingpong PROPERTY LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_C_LINK_FLAGS})
target_link_libraries(mpi-pingpong ${MPI_C_LIBRARIES} ${RASTRO_LIBRARIES} )

add_custom_target(mpi-run
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 8 ${CMAKE_CURRENT_BINARY_DIR}/mpi-pingpong
  DEPENDS mpi-pingpong
)

add_custom_target(mpi-pingpong-paje
  COMMAND ${TRACE_CONVERTER} -i pingpong.txt -o pingpong.paje *.rst
  DEPENDS mpi-run
)


