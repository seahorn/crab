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

LDD_LIBS := \
    $(CRAB_INSTALL)/ldd/lib/libldd.a \
    $(CRAB_INSTALL)/ldd/lib/libtvpi.a \
    $(CRAB_INSTALL)/ldd/lib/libcudd.a \
    $(CRAB_INSTALL)/ldd/lib/libst.a \
    $(CRAB_INSTALL)/ldd/lib/libutil.a \
    $(CRAB_INSTALL)/ldd/lib/libmtr.a \
    $(CRAB_INSTALL)/ldd/lib/libepd.a 

APRON_LIBS := \
    $(CRAB_INSTALL)/apron/lib/libpolkaMPQ.a \
    $(CRAB_INSTALL)/apron/lib/liboctD.a \
    $(CRAB_INSTALL)/apron/lib/liboptoct.a \
    $(CRAB_INSTALL)/apron/lib/liblinkedlistapi.a \
    $(CRAB_INSTALL)/apron/lib/libapron.a \
    $(CRAB_INSTALL)/apron/lib/libboxMPQ.a \
    $(CRAB_INSTALL)/apron/lib/libitvMPQ.a

LDLIBS := \
     $(CRAB_LIBS) $(LDD_LIBS) $(APRON_LIBS) \
     -lmpfr -lgmp -lm -lstdc++ 

CXXFLAGS := \
     -Wall -Werror -Wfatal-errors \
     -Wno-unused-local-typedefs -Wno-inconsistent-missing-override -Wno-unused-const-variable \
     -Wno-uninitialized -Wno-deprecated \
     -DBSD -DHAVE_IEEE_754 -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8 \
     -I $(CRAB_INSTALL)/crab/include \
     -I $(CRAB_INSTALL)/ldd/include/ldd \
     -I $(CRAB_INSTALL)/apron/include \
     -O2 -g -std=c++11

ifeq ($(UNAME),Darwin)
## For Mojave
CXXFLAGS += -Wno-nullability-completeness
CXXFLAGS += -isystem /usr/local/include 
CXXFLAGS += -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk
endif

# Lookup paths for crab shared libraries
#LDFLAGS := -Wl,-rpath,$(CRAB_INSTALL)/crab/lib/
# For gmp library (you might need to adapt this)
LDFLAGS := -L/usr/local/lib

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
	@echo "Installing crab in ${CRAB_INSTALL} ... "
	cd crab/build \
	     && cmake -DCMAKE_INSTALL_PREFIX=$(CRAB_INSTALL) -DCRAB_BUILD_LIBS_SHARED=OFF \
	              -DCRAB_USE_LDD=ON -DCRAB_USE_APRON=ON ../ \
	     && cmake --build . --target ldd && cmake ../ \
	     && cmake --build . --target apron && cmake ../ \
	     && cmake --build . --target install
