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
CRAB_LIBS = $(CRAB_INSTALL)/crab/lib/libCrab.a
LDLIBS := $(CRAB_LIBS) -lgmp -lm -lstdc++ 

CXXFLAGS := \
     -Wall -Werror -Wfatal-errors \
     -Wno-unused-local-typedefs \
     -Wno-uninitialized \
     -DBSD -DHAVE_IEEE_754 -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8 \
     -I $(CRAB_INSTALL)/crab/include \
     -O2 -g -std=c++11

ifeq ($(UNAME),Darwin)
## For Mojave
CXXFLAGS += -Wno-nullability-completeness
CXXFLAGS += -isystem /usr/local/include 
CXXFLAGS += -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk
endif

# Lookup paths for crab shared libraries
# LDFLAGS := -Wl,-rpath,$(CRAB_INSTALL)/crab/lib/
# For gmp library (you might need to adapt this)
LDFLAGS := -L/usr/local/lib

all: domain analysis

%.o:  $(CXX) $(CXXFLAGS) $< -c -o $@

domain: domain.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) domain.o $(LDLIBS) -o domain

analysis: analysis.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) analysis.o $(LDLIBS) -o analysis

clean:
	rm -f domain analysis *.a *.o

crab_install:
	git clone https://github.com/seahorn/crab.git crab
	mkdir -p crab/build
	cd crab/build \
	    && cmake -DCMAKE_INSTALL_PREFIX=$(CRAB_INSTALL)  \
	       -DCRAB_BUILD_LIBS_SHARED=OFF ../ \
	    && cmake --build . --target install
