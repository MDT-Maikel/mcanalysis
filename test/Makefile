### Compiler & Options ###
CC = g++
COPTS = -std=c++0x
CFLAGS = -Wall
LDFLAGS =
COMPILER = $(CC) $(COPTS)

### Libraries ###
LBOOST = -lboost_filesystem -lboost_system
LGZSTREAM = -L../deps/gzstream/ -lz -lgzstream
LROOT = -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl
LINKS = $(LROOT)

### Dependencies ###
DEPS = -I$(ROOTSYS)/include

### Object Files ###
OBJ = ../source/histogram/histogram.o

### Make EXE ###

all: test

test: test.o $(OBJ)
	$(COMPILER) $(LDFLAGS) $(DEPS) -o $@ $+ $(LINKS) 

### Make Clean ###
clean:
	rm -f $(EXE) $(OBJ)

### Make Targets ###

test.o: test.cpp 
	$(COMPILER) $(CFLAGS) $(DEPS) -c -o $@ $<

../source/histogram/histogram.o: ../source/histogram/histogram.cpp ../source/histogram/histogram.h 
	$(COMPILER) $(CFLAGS) $(DEPS) -c -o $@ $<


#$@ refers to the target file
#$< refers to first dependency
#$+ refers to all dependencies

