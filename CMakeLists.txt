CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
SET( CMAKE_CXX_COMPILER g++ )
SET( CMAKE_C_COMPILER gcc )

MESSAGE( STATUS "Building in ${BUILDENV} Environment" )

Project( mcflow )

SET(CMAKE_VERBOSE_MAKEFILE ON)

SET( CMAKE_BUILD_TYPE Debug )
SET( CMAKE_CXX_FLAGS "-O3 -DBOOST_NO_HASH -g -funroll-loops -frounding-math ${CMAKE_CXX_FLAGS}")


FIND_PACKAGE( Boost 1.42 COMPONENTS regex program_options )

# FIND_PACKAGE( Qt4 REQUIRED )

# SET(QT_USE_QTMAIN 1)
# SET(QT_USE_QTOPENGL 1)
# SET(QT_USE_QTXMLPATTERNS 1)

SET( PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR} )
SET( CMAKE_INSTALL_PREFIX ${PROJECT_ROOT} )
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_ROOT}/build/bin )


#execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/shaders
#                                                           ${CMAKE_BINARY_DIR}/shaders )

#execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/data
#                                                           ${CMAKE_BINARY_DIR}/data )

SUBDIRS( src )


