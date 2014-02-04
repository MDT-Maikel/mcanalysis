#########################################
### CMake file for mcanalysis library ###
#########################################

## library source
set(MCANALYSIS_SOURCES
	particle/particle.h
	particle/particle.cpp
	particle/lhe.h
	particle/lhe.cpp
	particle/lhco.h
	particle/lhco.cpp
	event/event.h
	event/event.cpp
	event/mt2_bisect.h
	event/mt2_bisect.cpp
	utility/utility.h
	utility/utility.cpp
	utility/utility.tpp
	cuts/cuts.h
	cuts/cuts.cpp
	cuts/cuts_imp.h
	histogram/histogram.h
	histogram/histogram.cpp
	histogram/histogram2D.h
	histogram/histogram2D.cpp
	histogram/interpolate.py
	plot/plot.h
	plot/plot.cpp
	plot/plot2d.h
	plot/plot2d.cpp
	bumphunter/bumphunter.h
	bumphunter/bumphunter.cpp
	jet_analysis/jet_analysis.h
	jet_analysis/jet_analysis.cpp
	jet_analysis/HEPTopTagger.h
	jet_analysis/HEPTopTagger.cpp
)

## Include directories
include_directories(
	"${Boost_INCLUDE_DIRS}"
	"${PYTHON_INCLUDE_DIRS}"
	"${ROOT_INCLUDE_DIRS}"
	"${PYTHIA8_INCLUDE_DIRS}"
	"${FASTJET_INCLUDE_DIRS}"
)

## Make the library
add_library(analysis ${MCANALYSIS_SOURCES})