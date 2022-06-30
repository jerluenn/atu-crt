cmake_minimum_required(VERSION 3.10)
project(atu-crt)

set (CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC -fopenmp  -DBOOST_LOG_DYN_LINK")
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
SET( LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/shared )
SET( SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src )
SET( TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests )
SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib/static)
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_PATH} ) 
SET( acados_include "~/acados/include" )
SET( acados_lib "~/acados/lib")
SET( c_functions "${SOURCE_DIR}/tether_units/c_generated_code_tether_units")

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

find_package(Eigen3 3.4 REQUIRED)
if(NOT EIGEN3_FOUND)
  # Fallback to cmake_modules
  find_package(cmake_modules REQUIRED)
  find_package(Eigen REQUIRED)
  set(EIGEN3_INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS})
  set(EIGEN3_LIBRARIES ${EIGEN_LIBRARIES})
endif()
find_package(PythonLibs REQUIRED)

## System dependencies are found with CMake's conventions
find_package(Boost REQUIRED COMPONENTS system log_setup log REQUIRED)
find_package(catkin REQUIRED COMPONENTS message_generation std_msgs sensor_msgs geometry_msgs roscpp)

INCLUDE_DIRECTORIES(

	${SOURCE_DIR}/tether_units
	${acados_include}
	${acados_include}/blasfeo/include/
	${acados_include}/hpipm/include/
	${acados_include}/acados/
	${acados_include}/qpOASES_e/
	${PROJECT_SOURCE_DIR}/utils
	${LIB_PATH}
	${PYTHON_INCLUDE_DIRS}

)


FILE( GLOB SOURCES ${SOURCE_DIR}/tether_units/*.cpp )
FILE( GLOB UTIL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cpp )
FILE( GLOB TEST_FILES ${TEST_DIR}/*.cpp)
SUBDIRLIST(C_FUNC_TETHER_LIST ${c_functions})

FOREACH(SUBDIR ${C_FUNC_TETHER_LIST})

	SUBDIRLIST(SUBSUBDIRS ${c_functions}/${SUBDIR})
	
	FOREACH(SUBSUBDIR ${SUBSUBDIRS})

		INCLUDE_DIRECTORIES(${c_functions}/${SUBDIR}/${SUBSUBDIR})
		INCLUDE_DIRECTORIES(${c_functions}/${SUBDIR})

	ENDFOREACH()

ENDFOREACH(SUBDIR ${c_functions})

ADD_LIBRARY(MathUtils ${CMAKE_CURRENT_SOURCE_DIR}/utils/MathUtils.cpp)
TARGET_LINK_LIBRARIES(MathUtils 
PRIVATE Eigen3::Eigen 
${PYTHON_LIBRARIES} 
${acados_lib}/libacados.so)

ADD_LIBRARY(IntegratorInterface ${CMAKE_CURRENT_SOURCE_DIR}/utils/IntegratorInterface.cpp)
TARGET_LINK_LIBRARIES(IntegratorInterface 
PUBLIC MathUtils
PRIVATE Eigen3::Eigen 
${PYTHON_LIBRARIES} 
${acados_lib}/libacados.so)

ADD_LIBRARY(TetherUnit_Solver ${CMAKE_CURRENT_SOURCE_DIR}/src/tether_units/TetherUnit_Solver.cpp)
TARGET_LINK_LIBRARIES(TetherUnit_Solver
PUBLIC IntegratorInterface 
PUBLIC MathUtils
PRIVATE Eigen3::Eigen 
${PYTHON_LIBRARIES} 
${acados_lib}/libacados.so)

FOREACH(TEST_FILE ${TEST_FILES})

  FILE( GLOB LIBFILES ${LIB_PATH}/*.so)
	GET_FILENAME_COMPONENT( EXEC_NAME ${TEST_FILE} NAME_WE)
	ADD_EXECUTABLE(${EXEC_NAME} ${TEST_FILE})
	TARGET_LINK_LIBRARIES(${EXEC_NAME} 
  PUBLIC TetherUnit_Solver
  PUBLIC MathUtils
  PUBLIC IntegratorInterface 
  PRIVATE Eigen3::Eigen 
  ${PYTHON_LIBRARIES} 
  ${acados_lib}/libacados.so
  ${LIBFILES}
)

ENDFOREACH(TEST_FILE ${TEST_FILES})

