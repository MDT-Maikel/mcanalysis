#########################################
### CMake file for mcanalysis project ###
#########################################
cmake_minimum_required(VERSION 2.6)
project(mcanalysis)

## CMake search paths
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}" ${CMAKE_MODULE_PATH})

## Compiler settings
set(CMAKE_CXX_FLAGS "-std=c++0x")

## Library: Zlib
find_package(ZLIB REQUIRED)

## Library: GZStream
set(GZSTREAM_INCLUDE_DIRS "deps/gzstream")
add_library(gzstream STATIC IMPORTED)
set_property(TARGET gzstream PROPERTY IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/deps/gzstream/libgzstream.a)
set(GZSTRREAM_LIBRARIES gzstream)

## Library: Boost
find_package(Boost 1.40.0 COMPONENTS system filesystem python REQUIRED)

## Library: Python
find_package(PythonLibs REQUIRED)

## Library: Root
find_package(ROOT REQUIRED)
set(ROOT_LIBRARIES ${ROOT_LIBRARIES} -lMathMore) ##TODO: make this work better
set(ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIR})

## Library: Pythia8
find_package(Pythia8 REQUIRED)

## Library: FastJet
set(fastjet_USER_INCLUDE_PATH "$ENV{FASTJET_DIR}/include")
set(fastjet_USER_LIB_PATH "$ENV{FASTJET_DIR}/lib")
find_package(FastJet REQUIRED)
set(FASTJET_INCLUDE_DIRS ${fastjet_INCLUDE_DIRS})
set(FASTJET_LIBRARIES ${fastjet_LIBRARIES})

## Library: mcanalysis
set(MCANALYSIS_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/source")
add_subdirectory(source)
set(MCANALYSIS_LIBRARIES analysis)

## Include directories
message("====================================")
message("Using following include dirs:")
message("zlib: " "${ZLIB_INCLUDE_DIRS}")
message("gzstream: " "${GZSTREAM_INCLUDE_DIRS}")
message("python: " "${PYTHON_INCLUDE_DIRS}")
message("boost: " "${Boost_INCLUDE_DIRS}")
message("root: " "${ROOT_INCLUDE_DIRS}")
message("pythia8: " "${PYTHIA8_INCLUDE_DIRS}")
message("fastjet: " "${fastjet_INCLUDE_DIRS}")
message("mcanalysis: " "${MCANALYSIS_INCLUDE_DIRS}")
message("====================================")
#include_directories(
	#"${ZLIB_INCLUDE_DIRS}"
	#"${GZSTREAM_INCLUDE_DIRS}"
	#"${PYTHON_INCLUDE_DIRS}"
	#"${Boost_INCLUDE_DIRS}"
	#"${ROOT_INCLUDE_DIRS}"
	#"${PYTHIA8_INCLUDE_DIRS}"
	#"${fastjet_INCLUDE_DIRS}"
	#"${MCANALYSIS_INCLUDE_DIRS}"
#)

## Print libraries
message("====================================")
message("zlib: " "${ZLIB_LIBRARIES}")
message("gzstream: " "${GZSTRREAM_LIBRARIES}")
message("python: " "${PYTHON_LIBRARIES}")
message("boost: " "${Boost_LIBRARIES}")
message("root: " "${ROOT_LIBRARIES}")
message("pythia8: " "${PYTHIA8_LIBRARIES}")
message("fastjet: " "${FASTJET_LIBRARIES}")
message("mcanalysis: " "${MCANALYSIS_LIBRARIES}")
message("====================================")

## Executables: tests
add_subdirectory(tests)

## Executables: tools
add_subdirectory(tools)

