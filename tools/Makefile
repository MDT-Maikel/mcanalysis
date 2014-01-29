### Compiler & Options ###
CC = g++
COPTS = -std=c++11
CFLAGS = -Wall -g 
LDFLAGS =
COMPILER = $(CC) $(COPTS)

### Libraries ###
LBOOST = -lboost_filesystem -lboost_system
LGZSTREAM = -L../deps/gzstream/ -lz -lgzstream
LROOT = -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGpad -lGraf3d -lGraf -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl
LINKS = $(LGZSTREAM) $(LBOOST) $(LROOT)
# -lTree -lRint 

### Dependencies ###
DEPSLOCAL = ../deps/gzstream/gzstream.h ../deps/gzstream/gzstream.C 
DEPSGLOBAL = -I$(ROOTSYS)/include
DEPS = $(DEPSLOCAL) $(DEPSGLOBAL)

### Object Files ###
OBJPART = ../source/particle/lhco.o ../source/particle/particle.o ../source/particle/lhe.o
OBJEVNT = ../source/event/event.o 
OBJUTIL = ../source/utility/utility.o
OBJHIST = ../source/histogram/histogram.o
OBJPLOT = ../source/plot/plot.o ../source/plot/plot_standard.o
OBJCUTS = ../source/cuts/cuts.o ../source/cuts/cuts_standard.o
OBJ = $(OBJPART) $(OBJEVNT) $(OBJUTIL) $(OBJHIST) $(OBJPLOT) $(OBJCUTS)

### Make EXES ###
all: mergemc cutmc plotmc

mergemc: merge/mergemc.o $(OBJ)
	$(COMPILER) $(LDFLAGS) $(DEPS) -o $@ $+ $(LINKS) 

cutmc: cut/cutmc.o $(OBJ)
	$(COMPILER) $(LDFLAGS) $(DEPS) -o $@ $+ $(LINKS) 

plotmc: plot/plotmc.o $(OBJ)
	$(COMPILER) $(LDFLAGS) $(DEPS) -o $@ $+ $(LINKS) 

### Make Clean ###
clean:
	rm -f mergemc merge/mergemc.o cutmc cut/cutmc.o plotmc plot/plotmc.o $(OBJ)

### Make Targets ###

merge/mergemc.o: merge/mergemc.cpp 
	$(COMPILER) $(CFLAGS) $(DEPSGLOBAL) -c -o $@ $<

cut/cutmc.o: cut/cutmc.cpp 
	$(COMPILER) $(CFLAGS) $(DEPSGLOBAL) -c -o $@ $<

plot/plotmc.o: plot/plotmc.cpp 
	$(COMPILER) $(CFLAGS) $(DEPSGLOBAL) -c -o $@ $<

../source/particle/particle.o: ../source/particle/particle.cpp ../source/particle/particle.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/particle/lhco.o: ../source/particle/lhco.cpp ../source/particle/lhco.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/particle/lhe.o: ../source/particle/lhe.cpp ../source/particle/lhe.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/event/event.o: ../source/event/event.cpp ../source/event/event.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/utility/utility.o: ../source/utility/utility.cpp ../source/utility/utility.h ../source/utility/utility.tpp
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/histogram/histogram.o: ../source/histogram/histogram.cpp ../source/histogram/histogram.h
	$(COMPILER) $(CFLAGS) $(DEPSGLOBAL) -c -o $@ $<

../source/plot/plot.o: ../source/plot/plot.cpp ../source/plot/plot.h
	$(COMPILER) $(CFLAGS) $(DEPSGLOBAL) -c -o $@ $<

../source/plot/plot_standard.o: ../source/plot/plot_standard.cpp ../source/plot/plot_standard.h
	$(COMPILER) $(CFLAGS) $(DEPSGLOBAL) -c -o $@ $<

../source/cuts/cuts.o: ../source/cuts/cuts.cpp ../source/cuts/cuts.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/cuts/cuts_standard.o: ../source/cuts/cuts_standard.cpp ../source/cuts/cuts_standard.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

#$@ refers to the target file
#$< refers to first dependency
#$+ refers to all dependencies

