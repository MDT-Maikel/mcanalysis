################################
### CMake find HepMC library ###
################################
# This module defines
# HEPMC_INCLUDE_DIR  where to locate Pythia.h file
# HEPMC_LIBRARIES    the libraries 
# HEPMC_FOUND        whether Pythia8 was found

## HepMC not found as default
set(HEPMC_FOUND false)

## Find the HepMC include directory
find_path(HEPMC_INCLUDE_DIRS HepMC/GenEvent.h
	$ENV{HEPMC_DIR}/include
	DOC "Specify the directory containing GenEvent.h."
)

## Find the HepMC library
find_library(HEPMC_LIBRARIES NAMES HepMC PATHS
	$ENV{HEPMC_DIR}/lib
	DOC "Specify the Pythia8 library here."
)

## Log whether the libraries where found
if(HEPMC_INCLUDE_DIRS AND HEPMC_LIBRARIES)
	set(HEPMC_FOUND true)
	message(STATUS "Found HepMC library at ${HEPMC_LIBRARIES}")
endif()

## Mark as advanced settings 
MARK_AS_ADVANCED(HEPMC_FOUND HEPMC_INCLUDE_DIRS HEPMC_LIBRARIES)
