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

## Print the status of each of the libraries
message("=======================================================")
message("==================== Library Status ===================")
message("=======================================================")
message("zlib: " "${ZLIB_FOUND}")
message("inc: " "${ZLIB_INCLUDE_DIRS}")
message("lib: " "${ZLIB_LIBRARIES}")
message("=======================================================")
message("gzstream: " "${GZSTRREAM_FOUND}")
message("inc: " "${GZSTREAM_INCLUDE_DIRS}")
message("lib: " "${GZSTRREAM_LIBRARIES}")
message("=======================================================")
message("python: " "${PYTHON_FOUND}")
message("inc: " "${PYTHON_INCLUDE_DIRS}")
message("lib: " "${PYTHON_LIBRARIES}")
message("=======================================================")
message("boost: " "${Boost_FOUND}")
message("inc: " "${Boost_INCLUDE_DIRS}")
message("lib: " "${Boost_LIBRARIES}")
message("=======================================================")
message("root: " "${ROOT_FOUND}")
message("inc: " "${ROOT_INCLUDE_DIRS}")
message("lib: " "${ROOT_LIBRARIES}")
message("=======================================================")
message("pythia8: " "${ROOT_FOUND}")
message("inc: " "${PYTHIA8_INCLUDE_DIRS}")
message("lib: " "${PYTHIA8_LIBRARIES}")
message("=======================================================")
message("fastjet: " "${ROOT_FOUND}")
message("inc: " "${FASTJET_INCLUDE_DIRS}")
message("lib: " "${FASTJET_LIBRARIES}")
message("=======================================================")
message("mcanalysis: " "${MCANALYSIS_FOUND}")
message("inc: " "${MCANALYSIS_INCLUDE_DIRS}")
message("lib: " "${MCANALYSIS_LIBRARIES}")
message("=======================================================")

## Executables: tests
add_subdirectory(tests)

## Executables: tools
add_subdirectory(tools)

