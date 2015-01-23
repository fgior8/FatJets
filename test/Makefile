CC=g++
#wildcard, this finds all *.c file and make *.o files
objects = $(patsubst %.c,%.o,$(wildcard *.c))

#for linking fftw3 routines
FFTWDIR := $(PWD)/fftw
FFTWIND := -I$(FFTWDIR)/include
FFTWLIB := -L$(FFTWDIR)/lib -lfftw3

ROOTCFLAGS    := $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      := $(shell $(ROOTSYS)/bin/root-config --libs)
#ROOTGLIBS     := $(shell $(ROOTSYS)/bin/root-config --glibs) -L TMVA/lib
ROOTGLIBS     := $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = g++ -m64
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2 -Xlinker -zmuldefs 
LD            = g++ -m64
LDFLAGS       = -g
SOFLAGS       = -shared


all: analyzer
covariance.o : matrix_inverse.h

#rule for basic objects
%.o: %.c %.h
	$(CXX) -Wall -g -c $< -o $@

dataset.o: dataset.c dataset.h
	$(CXX) -Wall -g -c $< \
	$(FFTWIND) $(FFTWLIB) \
	-o $@

example: $(objects) example.C
	$(CXX) -O3 $(objects) tensor.h example.C \
	$(ROOTFLAGS) \
	$(ROOTLIBS) \
	$(FFTWIND) $(FFTWLIB) \
	-o example.exe

analyzer: $(objects) runAnalyzer.cc 
	$(CXX) $(objects) tensor.h Data.cc ChainMaker.cc SignalPlots.cc BTagSFUtil.C Analyzer.cc runAnalyzer.cc \
	$(ROOTCFLAGS) \
	$(ROOTLIBS) \
	$(FFTWIND) $(FFTWLIB) \
	-o runAnalyzer.exe 

clean:
	rm -f *.o
