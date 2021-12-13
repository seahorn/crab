#------------------------------------------------------------------------------------#
#                                  INSTRUCTIONS
#------------------------------------------------------------------------------------#
#  0. Define variables BOOST_ROOT and GMP_DIR
#  1. Type `make crab_install`
#  2. Type `make all`
#  3. Type `./domain` or `./analysis`
#------------------------------------------------------------------------------------#

BOOST_ROOT=
GMP_DIR=

PWD=$(shell pwd)
UNAME := $(shell uname)

CRAB_INSTALL=${PWD}/crab-install
CRAB_LIBS = $(CRAB_INSTALL)/crab/lib/libCrab.a

CXXFLAGS := \
     -Wall -Werror -Wfatal-errors \
     -Wno-unused-local-typedefs \
     -Wno-uninitialized \
     -DBSD -DHAVE_IEEE_754 -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8 \
     -I $(CRAB_INSTALL)/crab/include \
     -I $(BOOST_ROOT)/include \
     -I $(GMP_DIR)/include \
     -O2 -g -std=c++11

ifeq ($(UNAME),Darwin)
CXXFLAGS += -Wno-nullability-completeness
CXXFLAGS += -isysroot $(shell xcrun --show-sdk-path)
endif

LDFLAGS := $(CRAB_LIBS)
LDFLAGS += -L$(GMP_DIR)/lib -lgmp

all: domain analysis

%.o:  $(CXX) $(CXXFLAGS) $< -c -o $@

domain: domain.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) domain.o -o domain

analysis: analysis.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) analysis.o -o analysis

clean:
	rm -f domain analysis *.a *.o

crab_install:
	git clone https://github.com/seahorn/crab.git crab
	mkdir -p crab/build
	cd crab/build \
	    && cmake -DCMAKE_INSTALL_PREFIX=$(CRAB_INSTALL)  \
	       -DCRAB_BUILD_LIBS_SHARED=OFF ../ \
	    && cmake --build . --target install
