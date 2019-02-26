#------------------------------------------------------------------------------------#
#                                  INSTRUCTIONS
#------------------------------------------------------------------------------------#
#  1. Type `make crab_install`
#  2. Type `make all`
#  3. Type `./domain` or `./analysis`
#------------------------------------------------------------------------------------#

PWD=$(shell pwd)
CRAB_INSTALL=${PWD}/crab-install

UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
    DYLIB=dylib
else
    DYLIB=so
endif
CRAB_LIBS = $(CRAB_INSTALL)/crab/lib/libCrab.${DYLIB}

LDLIBS := \
     $(CRAB_LIBS) -lgmpxx -lgmp -lm -lstdc++ 

CXXFLAGS := \
     -Wall -Werror -Wfatal-errors \
     -Wno-unused-local-typedefs -Wno-inconsistent-missing-override -Wno-unused-const-variable \
     -Wno-uninitialized -Wno-deprecated \
     -DBSD -DHAVE_IEEE_754 -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8 \
     -I $(CRAB_INSTALL)/crab/include \
     -O2 -g -std=c++11

# Lookup paths for crab shared libraries
LDFLAGS := -Wl,-rpath,$(CRAB_INSTALL)/crab/lib/

all: domain analysis

%.o:  $(CXX) $(CXXFLAGS) $< -c -o $@

domain: domain.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) domain.o $(LDLIBS) -o domain

analysis: analysis.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) analysis.o $(LDLIBS) -o analysis

clean:
	rm -f domain analysis  *.a *.o

crab_install:
	git clone https://github.com/seahorn/crab.git crab
	mkdir -p crab/build
	cd crab/build \
	    && cmake -DCMAKE_INSTALL_PREFIX=$(CRAB_INSTALL)  ../ \
	    && cmake --build . --target install
