###########################
### CMake file for ttag ###
###########################

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

## Executable: tztag_1top
add_executable(tztag_1top tztag_1top.cpp)
target_link_libraries(
	tztag_1top
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	${PYTHIA8_LIBRARIES}
	${FASTJET_LIBRARIES}
)
## Copy the tztag.cmnd file as well
add_custom_target(copy_tztag ALL COMMENT "Copying tztag.cmnd")
add_custom_command(TARGET copy_tztag COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/ttag/tztag.cmnd ${PROJECT_BINARY_DIR}/ttag/tztag.cmnd)

## Executable: tztag_hunt
add_executable(tztag_hunt tztag_hunt.cpp)
target_link_libraries(
	tztag_hunt
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	${PYTHIA8_LIBRARIES}
	${FASTJET_LIBRARIES}
)
## Copy the tztag_hunt_1top.cmnd file as well
add_custom_target(copy_tztag_hunt_1top ALL COMMENT "Copying tztag_hunt_1top.cmnd")
add_custom_command(TARGET copy_tztag_hunt_1top COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/ttag/tztag_hunt_1top.cmnd ${PROJECT_BINARY_DIR}/ttag/tztag_hunt_1top.cmnd)
## Copy the tztag_hunt_kin_dist.cmnd file as well
add_custom_target(copy_tztag_hunt_kin_dist ALL COMMENT "Copying tztag_hunt_kin_dist.cmnd")
add_custom_command(TARGET copy_tztag_hunt_kin_dist COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/ttag/tztag_hunt_kin_dist.cmnd ${PROJECT_BINARY_DIR}/ttag/tztag_hunt_kin_dist.cmnd)
