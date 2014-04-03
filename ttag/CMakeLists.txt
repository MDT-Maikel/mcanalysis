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

## Executable: tztag_2top
add_executable(tztag_2top tztag_2top.cpp)
target_link_libraries(
	tztag_2top
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	${PYTHIA8_LIBRARIES}
	${FASTJET_LIBRARIES}
)
## Copy the tztag_2top.cmnd file as well
add_custom_target(copy_tztag_2top ALL COMMENT "Copying tztag_2top.cmnd")
add_custom_command(TARGET copy_tztag_2top COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/ttag/tztag_2top.cmnd ${PROJECT_BINARY_DIR}/ttag/tztag_2top.cmnd)

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
