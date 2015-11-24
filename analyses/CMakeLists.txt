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

## Executable: dijet_cms_exo-12-059
add_executable(dijet_cms_exo-12-059 dijet_cms_exo-12-059.cpp)
target_link_libraries(
	dijet_cms_exo-12-059
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
)

## Executable: paired_dijet_cms_exo-11-016
add_executable(paired_dijet_cms_exo-11-016 paired_dijet_cms_exo-11-016.cpp)
target_link_libraries(
	paired_dijet_cms_exo-11-016
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
)

## Executable: jetmass_atl_susy-2013-07
add_executable(jetmass_atl_susy-2013-07 jetmass_atl_susy-2013-07.cpp)
target_link_libraries(
	jetmass_atl_susy-2013-07
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
)

## Executable: fchi_atl
add_executable(fchi_atl fchi_atl.cpp)
target_link_libraries(
	fchi_atl
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	${PYTHIA8_LIBRARIES}
	${FASTJET_LIBRARIES}
)
