### Compiler & Options ###
CC = g++
COPTS = -std=c++11
CFLAGS = -Wall -g 
LDFLAGS =
COMPILER = $(CC) $(COPTS)

### Libraries ###
LBOOST = -lboost_filesystem -lboost_system
LGZSTREAM = -L../deps/gzstream/ -lz -lgzstream
LINKS = $(LGZSTREAM) $(LBOOST) 

### Dependencies ###
DEPS = ../deps/gzstream/gzstream.h ../deps/gzstream/gzstream.C

### Object Files ###
OBJ = \
	../source/particle/lhco.o ../source/particle/particle.o  ../source/particle/lhe.o ../source/event/event.o \
	../source/utility/utility.o ../source/utility/settings.o \
	../source/cuts/cuts.o ../source/cuts/cuts_standard.o ../source/cuts/cuts_advanced.o

### Make EXES ###

all: mergemc cutmc

mergemc: mergemc.o $(OBJ)
	$(COMPILER) $(LDFLAGS) $(DEPS) -o $@ $+ $(LINKS) 

cutmc: cutmc.o $(OBJ)
	$(COMPILER) $(LDFLAGS) $(DEPS) -o $@ $+ $(LINKS) 

### Make Clean ###
clean:
	rm -f mergemc mergemc.o cutmc cutmc.o $(OBJ)

### Make Targets ###

mergemc.o: mergemc.cpp 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

cutmc.o: cutmc.cpp 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/particle/particle.o: ../source/particle/particle.cpp ../source/particle/particle.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/particle/lhco.o: ../source/particle/lhco.cpp ../source/particle/lhco.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/particle/lhe.o: ../source/particle/lhe.cpp ../source/particle/lhe.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/event/event.o: ../source/event/event.cpp ../source/event/event.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/utility/utility.o: ../source/utility/utility.cpp ../source/utility/utility.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/utility/settings.o: ../source/utility/settings.cpp ../source/utility/settings.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/cuts/cuts.o: ../source/cuts/cuts.cpp ../source/cuts/cuts.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/cuts/cuts_standard.o: ../source/cuts/cuts_standard.cpp ../source/cuts/cuts_standard.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/cuts/cuts_advanced.o: ../source/cuts/cuts_advanced.cpp ../source/cuts/cuts_advanced.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

#$@ refers to the target file
#$< refers to first dependency
#$+ refers to all dependencies

