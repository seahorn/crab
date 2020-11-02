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
ifeq ($(UNAME),Linux)
    LDFLAGS += -Wl,--disable-new-dtags 
endif

LDLIBS := \
    $(CRAB_LIBS) $(ELINA_LIBS) \
    -lmpfr -lgmp -lm -lstdc++ 

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

ifeq ($(UNAME),Darwin)
## For Mojave
CXXFLAGS += -Wno-nullability-completeness
CXXFLAGS += -isystem /usr/local/include 
CXXFLAGS += -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk
endif
# For gmp library (you might need to adapt this)
LDFLAGS := -L/usr/local/lib

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
	git clone https://github.com/seahorn/crab.git crab -b dev
	mkdir -p crab/build
	cd crab/build \
	    && cmake -DCMAKE_INSTALL_PREFIX=$(CRAB_INSTALL) -DCRAB_BUILD_LIBS_SHARED=OFF \
	             -DCRAB_USE_ELINA=ON ../ \
	    && cmake --build . --target elina && cmake ../ \
	    && cmake --build . --target install
