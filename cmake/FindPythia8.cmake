###################################
### CMake find Pythia 8 library ###
###################################
# This module defines
# PYTHIA8_INCLUDE_DIRS where to locate Pythia.h file
# PYTHIA8_LIBRARIES    the libraries 
# PYTHIA8_FOUND        whether Pythia8 was found

## Pythia8 not found as default
set(PYTHIA8_FOUND false)
set(PYTHIA8_HEPMC_FOUND false)

## Find the Pythia8 include directory
find_path(PYTHIA8_INCLUDE_DIRS Pythia8/Pythia.h
	$ENV{PYTHIA8_DIR}/include
	DOC "Specify the directory containing Pythia.h."
)

## Find the Pythia8 library
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

## Find the Pythia8 HepMC library
## but only if HepMC library has been found
if (HEPMC_FOUND)
	find_library(PYTHIA8_LIBRARY_HEPMC NAMES pythia8tohepmc PATHS
		$ENV{PYTHIA8_DIR}/lib/archive
		$ENV{PYTHIA8_DIR}/lib
		DOC "Specify the Pythia8 library here."
	)
endif()

## Combine all the Pythia8 libraries into one
set(PYTHIA8_LIBRARIES ${PYTHIA8_LIBRARY} ${PYTHIA8_LIBRARY_LHAPDF} ${PYTHIA8_LIBRARY_HEPMC})

## Log whether the libraries where found
if(PYTHIA8_INCLUDE_DIRS AND PYTHIA8_LIBRARY AND PYTHIA8_LIBRARY_LHAPDF)
	set(PYTHIA8_FOUND true)
	message(STATUS "Found Pythia8 library at ${PYTHIA8_LIBRARIES}")
	if (PYTHIA8_LIBRARY_HEPMC)
		set(PYTHIA8_HEPMC_FOUND true)
		message(STATUS "Found Pythia8 - HepMC library as well")
	endif()
endif()

## Mark as advanced settings 
MARK_AS_ADVANCED(PYTHIA8_FOUND PYTHIA8_HEPMC_FOUND PYTHIA8_INCLUDE_DIRS PYTHIA8_LIBRARIES)
