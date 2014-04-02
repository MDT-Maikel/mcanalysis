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
add_custom_target(copy_tztag_basic ALL COMMENT "Copying tztag_basic.cmnd")
add_custom_command(TARGET copy_tztag_basic COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/ttag/tztag_basic.cmnd ${PROJECT_BINARY_DIR}/ttag/tztag_basic.cmnd)

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
