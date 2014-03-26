############################
### CMake file for tools ###
############################

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

## Executable: mergemc
add_executable(mergemc mergemc.cpp)
target_link_libraries(
	mergemc
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
)

## Executable: plotmc
add_executable(plotmc plotmc.cpp)
target_link_libraries(
	plotmc
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
)

## Executable: cutmc
add_executable(cutmc cutmc.cpp)
target_link_libraries(
	cutmc
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
)

## Executable: gen_hepmc
## only build if HepMC is available
if(HEPMC_FOUND AND PYTHIA8_HEPMC_FOUND)
	add_executable(gen_hepmc gen_hepmc.cpp)
	target_link_libraries(
		gen_hepmc
		${HEPMC_LIBRARIES}
		${PYTHIA8_LIBRARIES}
		${Boost_LIBRARIES}
	)
	## Copy the gen_hepmc.cmnd file as well
	add_custom_target(copy_gen_hepmc ALL COMMENT "Copying gen_hepmc.cmnd")
	add_custom_command(TARGET copy_gen_hepmc COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tools/gen_hepmc.cmnd ${PROJECT_BINARY_DIR}/tools/gen_hepmc.cmnd)
endif()

## Executable: tztag
add_executable(tztag tztag.cpp)
target_link_libraries(
	tztag
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	${PYTHIA8_LIBRARIES}
	${FASTJET_LIBRARIES}
)

## Executable: tztag_basic
add_executable(tztag_basic tztag_basic.cpp)
target_link_libraries(
	tztag_basic
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	${PYTHIA8_LIBRARIES}
	${FASTJET_LIBRARIES}
)
## Copy the tztag_basic.cmnd file as well
add_custom_target(copy_tztag ALL COMMENT "Copying tztag_basic.cmnd")
add_custom_command(TARGET copy_tztag COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tools/tztag_basic.cmnd ${PROJECT_BINARY_DIR}/tools/tztag_basic.cmnd)

