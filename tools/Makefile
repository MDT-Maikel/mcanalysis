### Compiler & Options ###
CC = g++
COPTS = -std=c++11
CFLAGS = -Wall
LDFLAGS =
COMPILER = $(CC) $(COPTS)

### Libraries ###
LBOOST = -lboost_filesystem -lboost_system
LGZSTREAM = -L../deps/gzstream/ -lz -lgzstream
LINKS = $(LGZSTREAM) $(LBOOST) 

### Dependencies ###
DEPS = ../deps/gzstream/gzstream.h ../deps/gzstream/gzstream.C

### Object Files ###
OBJ = ../source/particle/lhco.o ../source/particle/particle.o  ../source/particle/lhe.o ../source/event/event.o ../source/utility/utility.o

### Make EXE ###

all: mergemc

mergemc: mergemc.o $(OBJ)
	$(COMPILER) $(LDFLAGS) $(DEPS) -o $@ $+ $(LINKS) 

### Make Clean ###
clean:
	rm -f $(EXE) $(OBJ)

### Make Targets ###

mergemc.o: mergemc.cpp 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/particle/particle.o: ../source/particle/particle.cpp ../source/particle/particle.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/particle/lhco.o: ../source/particle/lhco.cpp ../source/particle/lhco.h 
	$(COMPILER) $(CFLAGS) -c $< -o $@

../source/particle/lhe.o: ../source/particle/lhe.cpp ../source/particle/lhe.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/event/event.o: ../source/event/event.cpp ../source/event/event.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<

../source/utility/utility.o: ../source/utility/utility.cpp ../source/utility/utility.h 
	$(COMPILER) $(CFLAGS) -c -o $@ $<


#$@ refers to the target file
#$< refers to first dependency
#$+ refers to all dependencies

