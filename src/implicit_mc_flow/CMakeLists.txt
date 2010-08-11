INCLUDE(CheckLibraryExists)

SET( MCFLOW_SRCS apply_mc_flow.cpp ../arg_parser.cpp )

SET( MCFLOW_CHOLESKY_SRCS apply_mc_flow_cholesky.cpp ../arg_parser.cpp )

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -DTAUCS_MEMORY_TEST_yes" )

INCLUDE_DIRECTORIES(
  ${BOOST_INCLUDEDIR}
  /usr/local/include
  ${PROJECT_ROOT}/include/
  ${CMAKE_SOURCE_DIR}
  ./	
)


LINK_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/include
  ${BOOST_LIBRARYDIR}
  /usr/local/lib/OpenMesh
  /usr/lib
  /usr/local/lib
  ${PROJECT_ROOT}/lib
)

LINK_LIBRARIES( m pthread ${Boost_LIBRARIES}
  OpenMeshCore OpenMeshTools 
  taucs
  gfortran
  cblas
  lapack
  metis
  )
              
ADD_EXECUTABLE( apply_mc_flow ${MCFLOW_SRCS} )

ADD_EXECUTABLE( apply_mc_flow_cholesky ${MCFLOW_CHOLESKY_SRCS} )