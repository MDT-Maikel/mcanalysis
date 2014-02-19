#########################################
### CMake file for mcanalysis project ###
#########################################
cmake_minimum_required(VERSION 2.6)
project(mcanalysis)

## Set CMake search paths for user defined modules
## These are in the project's source dir
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

## Compiler settings
if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS "-O2 -std=c++0x")
endif()

## Library: Zlib
find_package(ZLIB REQUIRED)

## Library: GZStream
set(GZSTREAM_INCLUDE_DIRS "deps/gzstream")
add_library(gzstream STATIC IMPORTED)
set_property(TARGET gzstream PROPERTY IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/deps/gzstream/libgzstream.a)
set(GZSTREAM_LIBRARIES gzstream)
if (GZSTREAM_INCLUDE_DIRS AND GZSTREAM_LIBRARIES)
	set(GZSTREAM_FOUND true)
endif()

## Library: Boost
find_package(Boost 1.40.0 COMPONENTS system filesystem python REQUIRED)

## Library: Python
find_package(Python)
if (NOT PYTHON_FOUND)
	find_package(PythonLibs REQUIRED)
	if (PYTHON_INCLUDE_DIRS AND PYTHON_LIBRARIES)
		set(PYTHON_FOUND true)
	endif()
endif()

## Library: Root
find_package(ROOT REQUIRED)
set(ROOT_LIBRARIES ${ROOT_LIBRARIES} -lMathMore) ##TODO: make this work better
set(ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIR})

## Library: HepMC
find_package(HepMC)

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
if (MCANALYSIS_INCLUDE_DIRS AND MCANALYSIS_LIBRARIES)
	set(MCANALYSIS_FOUND true)
endif()

## Print the status of each of the libraries
message("=======================================================")
message("==================== Library Status ===================")
message("=======================================================")
if(ZLIB_FOUND)
	message("zlib: FOUND")
else() 
	message("zlib: NOT FOUND")
endif()
message("inc: " "${ZLIB_INCLUDE_DIRS}")
message("lib: " "${ZLIB_LIBRARIES}")
message("=======================================================")
if(GZSTREAM_FOUND)
	message("gzstream: FOUND")
else() 
	message("gzstream: NOT FOUND")
endif()
message("inc: " "${GZSTREAM_INCLUDE_DIRS}")
message("lib: " "${GZSTREAM_LIBRARIES}")
message("=======================================================")
if(PYTHON_FOUND)
	message("python: FOUND")
else() 
	message("python: NOT FOUND")
endif()
message("inc: " "${PYTHON_INCLUDE_DIRS}")
message("lib: " "${PYTHON_LIBRARIES}")
message("=======================================================")
if(Boost_FOUND)
	message("boost: FOUND")
else() 
	message("boost: NOT FOUND")
endif()
message("inc: " "${Boost_INCLUDE_DIRS}")
message("lib: " "${Boost_LIBRARIES}")
message("=======================================================")
if(ROOT_FOUND)
	message("root: FOUND")
else() 
	message("root: NOT FOUND")
endif()
message("inc: " "${ROOT_INCLUDE_DIRS}")
message("lib: " "${ROOT_LIBRARIES}")
message("=======================================================")
if(HEPMC_FOUND)
	message("hepmc: FOUND")
else() 
	message("hepmc: NOT FOUND")
endif()
message("inc: " "${HEPMC_INCLUDE_DIRS}")
message("lib: " "${HEPMC_LIBRARIES}")
message("=======================================================")
if(PYTHIA8_FOUND)
	message("pythia8: FOUND")
else() 
	message("pythia8: NOT FOUND")
endif()
message("inc: " "${PYTHIA8_INCLUDE_DIRS}")
message("lib: " "${PYTHIA8_LIBRARIES}")
message("=======================================================")
if(FASTJET_FOUND)
	message("fastjet: FOUND")
else() 
	message("fastjet: NOT FOUND")
endif()
message("inc: " "${FASTJET_INCLUDE_DIRS}")
message("lib: " "${FASTJET_LIBRARIES}")
message("=======================================================")
if(MCANALYSIS_FOUND)
	message("mcanalysis: FOUND")
else() 
	message("mcanalysis: NOT FOUND")
endif()
message("inc: " "${MCANALYSIS_INCLUDE_DIRS}")
message("lib: " "${MCANALYSIS_LIBRARIES}")
message("=======================================================")

## Executables: tests
add_subdirectory(tests)

## Executables: tools
add_subdirectory(tools)

## Executables: private
## this directory is for the users private programs
## a CMakeLists.txt must be included by the user
## only add to compilation if both are present
if(EXISTS "${PROJECT_SOURCE_DIR}/private/CMakeLists.txt")
	add_subdirectory(private)
endif()

## Executables: gres
if(INCLUDE_GRES)
	add_subdirectory(gres)
endif()
