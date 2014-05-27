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

## Executable: calc_val
add_executable(calc_val calc_val.cpp)
target_link_libraries(
	calc_val
	${Boost_LIBRARIES}
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
## Copy the tztag_1top.cmnd file as well
add_custom_target(copy_tztag_1top ALL COMMENT "Copying tztag_1top.cmnd")
add_custom_command(TARGET copy_tztag_1top COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/ttag/tztag_1top.cmnd ${PROJECT_BINARY_DIR}/ttag/tztag_1top.cmnd)

## Executable: tztag_dist
add_executable(tztag_dist tztag_dist.cpp)
target_link_libraries(
	tztag_dist
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	${PYTHIA8_LIBRARIES}
	${FASTJET_LIBRARIES}
)
## Copy the tztag_dist_pre_cuts.cmnd file as well
add_custom_target(copy_tztag_dist_pre_cuts ALL COMMENT "Copying tztag_dist_pre_cuts.cmnd")
add_custom_command(TARGET copy_tztag_dist_pre_cuts COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/ttag/tztag_dist_pre_cuts.cmnd ${PROJECT_BINARY_DIR}/ttag/tztag_dist_pre_cuts.cmnd)
## Copy the tztag_dist_after_cuts.cmnd file as well
add_custom_target(copy_tztag_dist_after_cuts ALL COMMENT "Copying tztag_dist_after_cuts.cmnd")
add_custom_command(TARGET copy_tztag_dist_after_cuts COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/ttag/tztag_dist_after_cuts.cmnd ${PROJECT_BINARY_DIR}/ttag/tztag_dist_after_cuts.cmnd)

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
## Copy the tztag_hunt.cmnd file as well
add_custom_target(copy_tztag_hunt ALL COMMENT "Copying tztag_hunt.cmnd")
add_custom_command(TARGET copy_tztag_hunt COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/ttag/tztag_hunt.cmnd ${PROJECT_BINARY_DIR}/ttag/tztag_hunt.cmnd)
