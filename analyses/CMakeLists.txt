###############################
### CMake file for analyses ###
###############################

## Mandatory include directories
include_directories(
	"${ZLIB_INCLUDE_DIRS}"
	"${GZSTREAM_INCLUDE_DIRS}"
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

## Executable: dijet_atl_conf-12-148
add_executable(dijet_atl_conf-12-148 dijet_atl_conf-12-148.cpp)
target_link_libraries(
	dijet_atl_conf-12-148
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
)

## Executable: cms_dijet
add_executable(dijet_cms_exo-12-059 dijet_cms_exo-12-059.cpp)
target_link_libraries(
	dijet_cms_exo-12-059
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
)
