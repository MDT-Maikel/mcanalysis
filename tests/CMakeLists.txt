############################
### CMake file for tests ###
############################

## Include directories
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

## Executable: test_all
add_executable(test_all source/test_all.cpp)

## Executable: test_lhco
add_executable(test_lhco source/test_lhco.cpp)
target_link_libraries(
	test_lhco
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
)

## Executable: test_lhe
add_executable(test_lhe source/test_lhe.cpp)
target_link_libraries(
	test_lhe
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
)

## Executable: test_event
add_executable(test_event source/test_event.cpp)
target_link_libraries(
	test_event
	${MCANALYSIS_LIBRARIES}
)

## Executable: test_cuts
add_executable(test_cuts source/test_cuts.cpp)
target_link_libraries(
	test_cuts
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${Boost_LIBRARIES}
)

## Executable: test_histogram
add_executable(test_histogram source/test_histogram.cpp)
target_link_libraries(
	test_histogram
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${PYTHON_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
)

## Executable: test_plot
add_executable(test_plot source/test_plot.cpp)
target_link_libraries(
	test_plot
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${PYTHON_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
)

## Executable: test_bumphunter
add_executable(test_bumphunter source/test_bumphunter.cpp)
target_link_libraries(
	test_bumphunter
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${PYTHON_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
)

## Executable: test_jets
add_executable(test_jets source/test_jets.cpp)
target_link_libraries(
	test_jets
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${PYTHON_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	${PYTHIA8_LIBRARIES}
	${FASTJET_LIBRARIES}
)

## Executable: test_heptop
add_executable(test_heptop source/test_heptop.cpp)
target_link_libraries(
	test_heptop
	${MCANALYSIS_LIBRARIES}
	${GZSTREAM_LIBRARIES}
	${ZLIB_LIBRARIES}
	${PYTHON_LIBRARIES}
	${Boost_LIBRARIES}
	${ROOT_LIBRARIES}
	${PYTHIA8_LIBRARIES}
	${FASTJET_LIBRARIES}
)


