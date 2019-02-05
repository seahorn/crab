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

LDD_LIBS := \
    $(CRAB_INSTALL)/ldd/lib/libldd.a \
    $(CRAB_INSTALL)/ldd/lib/libtvpi.a \
    $(CRAB_INSTALL)/ldd/lib/libcudd.a \
    $(CRAB_INSTALL)/ldd/lib/libst.a \
    $(CRAB_INSTALL)/ldd/lib/libutil.a \
    $(CRAB_INSTALL)/ldd/lib/libmtr.a \
    $(CRAB_INSTALL)/ldd/lib/libepd.a 

# Elina shared libraries have extension .so even on mac
ELINA_LIBS := \
    $(CRAB_INSTALL)/elina/lib/libelinalinearize.so \
    $(CRAB_INSTALL)/elina/lib/libelinaux.so \
    $(CRAB_INSTALL)/elina/lib/liboptoct.so \
    $(CRAB_INSTALL)/elina/lib/liboptpoly.so \
    $(CRAB_INSTALL)/elina/lib/liboptzones.so \
    $(CRAB_INSTALL)/elina/lib/libpartitions.so \
    $(CRAB_INSTALL)/elina/lib/libzonotope.so

# Lookup paths for crab and elina shared libraries
LDFLAGS := -Wl,-rpath,$(CRAB_INSTALL)/crab/lib/,-rpath,$(CRAB_INSTALL)/elina/lib/

LDLIBS := \
    $(CRAB_LIBS) $(LDD_LIBS) $(ELINA_LIBS) \
    -lmpfr -lgmpxx -lgmp -lm -lstdc++ 

CXXFLAGS := \
    -Wall -Werror -Wfatal-errors \
    -Wno-unused-function \
    -Wno-unused-local-typedefs -Wno-inconsistent-missing-override -Wno-unused-const-variable \
    -Wno-uninitialized -Wno-deprecated \
    -DBSD -DHAVE_IEEE_754 -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8 \
    -I $(CRAB_INSTALL)/crab/include \
    -I $(CRAB_INSTALL)/ldd/include/ldd \
    -I $(CRAB_INSTALL)/elina/include \
    -O2 -g -std=c++11

all: domain analysis
# Ugly hack: add elina lib directory to rpath is not enough
	@echo "Compilation done. \nType \"export LD_LIBRARY_PATH=$(CRAB_INSTALL)/elina/lib\""

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
	    && cmake -DCMAKE_INSTALL_PREFIX=$(CRAB_INSTALL) -DUSE_LDD=ON -DUSE_ELINA=ON ../ \
	    && cmake --build . --target ldd && cmake ../ \
	    && cmake --build . --target elina && cmake ../ \
	    && cmake --build . --target install
