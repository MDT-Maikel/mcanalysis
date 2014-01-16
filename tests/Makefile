### Compiler & Options ###
CC = g++
COPTS = -std=c++11
CFLAGS = -Wall -g 
LDFLAGS =
COMPILER = $(CC) $(COPTS)

### Locations ###
PYTHIA8LOCATION=/home/maikel/Programs/pythia8180
FASTJETLOCATION=/home/maikel/Programs/fastjet-install
PYTHONDIR=/usr

### Libraries ###
LBOOST = -lboost_filesystem -lboost_system -lboost_python
LPYTHON = -L$(PYTHONDIR)/lib -lpython2.7
LGZSTREAM = -L../deps/gzstream/ -lz -lgzstream
LROOT = -L$(ROOTSYS)/lib -lCore -lRIO -lNet -lHist -lGpad -lGraf3d -lGraf -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl
LPYTHIA = -L$(PYTHIA8LOCATION)/lib/archive -lpythia8 -llhapdfdummy
LFASTJET = -L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins`
LINKS = $(LBOOST) $(LPYTHON) $(LGZSTREAM) $(LROOT) $(LPYTHIA) $(LFASTJET)

### Dependencies ###
IPYTHON = -I$(PYTHONDIR)/include/python2.7
IGZSTREAM = ../deps/gzstream/gzstream.h ../deps/gzstream/gzstream.C 
IROOT = -I$(ROOTSYS)/include
IPYTHIA = -I$(PYTHIA8LOCATION)/include
IFASTJET = `$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`
DEPS = $(IPYTHON) $(IGZSTREAM) $(IROOT) $(IPYTHIA) $(IFASTJET)

### Object Files ###
OBJPART = ../source/particle/lhco.o ../source/particle/particle.o ../source/particle/lhe.o
OBJEVNT = ../source/event/event.o 
OBJUTIL = ../source/utility/utility.o
OBJHIST = ../source/histogram/histogram.o ../source/histogram/histogram2D.o
OBJPLOT = ../source/plot/plot.o ../source/plot/plot_standard.o
OBJCUTS = ../source/cuts/cuts.o ../source/cuts/cuts_standard.o ../source/cuts/cuts_advanced.o
OBJHUNT = ../source/bumphunter/bumphunter.o
OBJJETS = ../source/jet_analysis/jet_analysis.o
OBJ = $(OBJPART) $(OBJEVNT) $(OBJUTIL) $(OBJHIST) $(OBJPLOT) $(OBJCUTS) $(OBJHUNT) $(OBJJETS)

### Make EXES ###
EXE = test_histogram test_jetanalysis test_bumphunter
OBJEXES = source/test_histogram.o source/test_jetanalysis.o source/test_bumphunter.o

all: $(EXE)

test_histogram: source/test_histogram.o $(OBJHIST)
	$(COMPILER) $(LDFLAGS) \
	$(IROOT) \
	-o $@ $+ \
	$(LBOOST) \
	$(LROOT) \
	$(LPYTHON) \

test_jetanalysis: source/test_jetanalysis.o $(OBJJETS)
	$(COMPILER) $(LDFLAGS) \
	$(IPYTHIA) \
	$(IFASTJET) \
	-o $@ $+ \
	$(LPYTHIA) \
	$(LFASTJET) \

test_bumphunter: source/test_bumphunter.o $(OBJHUNT)
	$(COMPILER) $(LDFLAGS) \
	$(IROOT) \
	-o $@ $+ \
	$(LBOOST) \
	$(LROOT) \

### Make Clean ###
clean:
	rm -f $(EXE) \
	rm -f $(OBJEXES) \
	rm -f $(OBJ)

### Make Targets ###

source/test_histogram.o: source/test_histogram.cpp
	$(COMPILER) $(CFLAGS) \
	$(IROOT) \
	$(IPYTHON) \
	-c -o $@ $< \

source/test_jetanalysis.o: source/test_jetanalysis.cpp
	$(COMPILER) $(CFLAGS) \
	$(IPYTHIA) \
	$(IFASTJET) \
	-c -o $@ $< \

source/test_bumphunter.o: source/test_bumphunter.cpp 
	$(COMPILER) $(CFLAGS) \
	$(IROOT) \
	-c -o $@ $< \
	
### Make Sources ###	

../source/particle/particle.o: ../source/particle/particle.cpp ../source/particle/particle.h 
	$(COMPILER) $(CFLAGS) \ 
	-c -o $@ $< \

../source/particle/lhco.o: ../source/particle/lhco.cpp ../source/particle/lhco.h 
	$(COMPILER) $(CFLAGS) \
	-c -o $@ $< \

../source/particle/lhe.o: ../source/particle/lhe.cpp ../source/particle/lhe.h 
	$(COMPILER) $(CFLAGS) \
	-c -o $@ $< \

../source/event/event.o: ../source/event/event.cpp ../source/event/event.h 
	$(COMPILER) $(CFLAGS) \
	-c -o $@ $< \

../source/utility/utility.o: ../source/utility/utility.cpp ../source/utility/utility.h ../source/utility/utility.tpp
	$(COMPILER) $(CFLAGS) \
	-c -o $@ $< \

../source/histogram/histogram.o: ../source/histogram/histogram.cpp ../source/histogram/histogram.h
	$(COMPILER) $(CFLAGS) \
	$(IROOT) \
	-c -o $@ $< \
	
../source/histogram/histogram2D.o: ../source/histogram/histogram2D.cpp ../source/histogram/histogram2D.h
	$(COMPILER) $(CFLAGS) \
	$(IROOT) \
	$(IPYTHON) \
	-c -o $@ $< \

../source/plot/plot.o: ../source/plot/plot.cpp ../source/plot/plot.h
	$(COMPILER) $(CFLAGS) \
	$(IROOT) \
	-c -o $@ $< \

../source/plot/plot_standard.o: ../source/plot/plot_standard.cpp ../source/plot/plot_standard.h
	$(COMPILER) $(CFLAGS) \
	$(IROOT) \
	-c -o $@ $< \

../source/cuts/cuts.o: ../source/cuts/cuts.cpp ../source/cuts/cuts.h 
	$(COMPILER) $(CFLAGS) \
	-c -o $@ $< \

../source/cuts/cuts_standard.o: ../source/cuts/cuts_standard.cpp ../source/cuts/cuts_standard.h 
	$(COMPILER) $(CFLAGS) \
	-c -o $@ $< \

../source/cuts/cuts_advanced.o: ../source/cuts/cuts_advanced.cpp ../source/cuts/cuts_advanced.h 
	$(COMPILER) $(CFLAGS) \
	-c -o $@ $< \
	
../source/bumphunter/bumphunter.o: ../source/bumphunter/bumphunter.cpp ../source/bumphunter/bumphunter.h
	$(COMPILER) $(CFLAGS) \
	$(IROOT) \
	-c -o $@ $< \
	
../source/jet_analysis/jet_analysis.o: ../source/jet_analysis/jet_analysis.cpp ../source/jet_analysis/jet_analysis.h
	$(COMPILER) $(CFLAGS) \
	$(IPYTHIA) \
	$(IFASTJET) \
	-c -o $@ $< \


#$@ refers to the target file
#$< refers to first dependency
#$+ refers to all dependencies

