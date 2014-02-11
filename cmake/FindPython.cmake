#################################
### CMake find Python library ###
#################################
# This module defines
# PYTHON_INCLUDE_DIRS  where to locate Pythia.h file
# PYTHON_LIBRARIES     the libraries 
# PYTHON_FOUND         whether Pythia8 was found

## Python not found as default
set(PYTHON_FOUND false)

## Find the HepMC include directory
find_path(PYTHON_INCLUDE_DIRS pyconfig.h
	$ENV{PYTHON_DIR}/include/python2.7
	$ENV{PYTHON_DIR}/include
	DOC "Specify the Python directory containing pyconfig.h."
)

## Find the HepMC library
find_library(PYTHON_LIBRARIES NAMES python2.7 PATHS
	$ENV{PYTHON_DIR}/lib
	DOC "Specify the Python library here."
)

## Log whether the libraries where found
if(PYTHON_INCLUDE_DIRS AND PYTHON_LIBRARIES)
	set(PYTHON_FOUND true)
	message(STATUS "Found Python library at ${PYTHON_LIBRARIES}")
endif()

## Mark as advanced settings 
MARK_AS_ADVANCED(PYTHON_FOUND PYTHON_INCLUDE_DIRS PYTHON_LIBRARIES)
