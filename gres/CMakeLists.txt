###########################
### CMake file for gres ###
###########################

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

## Executable: gres_comb
add_executable(gres_comb gres_comb.cpp)

## Executable: gres_basic
add_executable(gres_basic gres_basic.cpp)
target_link_libraries(
	gres_basic
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${PYTHON_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
)
