############################
### CMake file for tools ###
############################

## Mandatory include directories
include_directories(
	"${ZLIB_INCLUDE_DIRS}"
	"${GZSTREAM_INCLUDE_DIRS}"
	"${PYTHON_INCLUDE_DIRS}"
	"${Boost_INCLUDE_DIRS}"
	"${ROOT_INCLUDE_DIRS}"
	"${PYTHIA8_INCLUDE_DIRS}"
	"${FASTJET_INCLUDE_DIRS}"
	"${MCANALYSIS_INCLUDE_DIRS}"
)

## Optional include directories
if(HEPMC_FOUND)
	include_directories("${HEPMC_INCLUDE_DIRS}")
endif()	

## Executable: mergemc
add_executable(mergemc source/mergemc.cpp)
target_link_libraries(
	mergemc
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${PYTHON_LIBRARIES}
	${Boost_LIBRARIES}
)

## Executable: gen_hepmc
## only build if HepMC is available
if(HEPMC_FOUND AND PYTHIA8_HEPMC_FOUND)
	add_executable(gen_hepmc source/gen_hepmc.cpp)
	target_link_libraries(
		gen_hepmc
		${HEPMC_LIBRARIES}
		${PYTHIA8_LIBRARIES}
	)
endif()

## Executable: ttag
add_executable(ttag source/ttag.cpp)
target_link_libraries(
	ttag
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${PYTHON_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	${PYTHIA8_LIBRARIES}
	${FASTJET_LIBRARIES}
)
