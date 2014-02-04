###################################
### CMake find Pythia 8 library ###
###################################
# This module defines
# PYTHIA8_INCLUDE_DIR  where to locate Pythia.h file
# PYTHIA8_LIBRARIES    the libraries 
# PYTHIA8_FOUND        whether Pythia8 was found

## Pythia 8 not found as default
set(PYTHIA8_FOUND false)

## Find the Pythia 8 include directory
find_path(PYTHIA8_INCLUDE_DIRS Pythia.h
	$ENV{PYTHIA8_DIR}/include/Pythia8
	$ENV{PYTHIA8_DIR}/include
	DOC "Specify the directory containing Pythia.h."
)
string(REPLACE "/Pythia8" "" PYTHIA8_INCLUDE_DIRS ${PYTHIA8_INCLUDE_DIRS})

## Find the Pythia 8 library
find_library(PYTHIA8_LIBRARY NAMES Pythia8 pythia8 PATHS
	$ENV{PYTHIA8_DIR}/lib/archive
	$ENV{PYTHIA8_DIR}/lib
	DOC "Specify the Pythia8 library here."
)

## Find the Pythia8 LHAPDF library
find_library(PYTHIA8_LIBRARY_LHAPDF NAMES lhapdfdummy PATHS
	$ENV{PYTHIA8_DIR}/lib/archive
	$ENV{PYTHIA8_DIR}/lib
	DOC "Specify the Pythia8 library here."
)

## Combine all the Pythia8 libraries into one
set(PYTHIA8_LIBRARIES ${PYTHIA8_LIBRARY} ${PYTHIA8_LIBRARY_LHAPDF})

## Log whether the libraries where found
if(PYTHIA8_INCLUDE_DIRS AND PYTHIA8_LIBRARIES)
	set(PYTHIA8_FOUND true)
	message(STATUS "Found Pythia8 library at ${PYTHIA8_LIBRARIES}")
endif()

## Mark as advanced settings 
MARK_AS_ADVANCED(PYTHIA8_FOUND PYTHIA8_LIBRARIES PYTHIA8_INCLUDE_DIRS)
