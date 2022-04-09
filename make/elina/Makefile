#------------------------------------------------------------------------------------#
#                                  INSTRUCTIONS
#------------------------------------------------------------------------------------#
#  0. Define variables BOOST_ROOT, GMP_DIR in Makefile.config and MPFR_DIR
#  1. Type `make crab_elina_install`
#  2. Type `make all`
#  3. Type `./domain` or `./analysis`
#------------------------------------------------------------------------------------#

include ../Makefile.config

MPFR_DIR=

# Elina shared libraries have extension .so even on mac
ELINA_LIBS := \
    $(CRAB_INSTALL)/elina/lib/libelinalinearize.so \
    $(CRAB_INSTALL)/elina/lib/libelinaux.so \
    $(CRAB_INSTALL)/elina/lib/liboptoct.so \
    $(CRAB_INSTALL)/elina/lib/liboptpoly.so \
    $(CRAB_INSTALL)/elina/lib/liboptzones.so \
    $(CRAB_INSTALL)/elina/lib/libpartitions.so \
    $(CRAB_INSTALL)/elina/lib/libzonotope.so

CXXFLAGS += \
     -I $(MPFR_DIR)/include \
     -I $(CRAB_INSTALL)/elina/include 

LDFLAGS += $(ELINA_LIBS)
LDFLAGS += -L$(MPFR_DIR)/lib -lmpfr

# Lookup paths for crab and elina shared libraries
LDFLAGS += -Wl,-rpath,$(CRAB_INSTALL)/crab/lib/,-rpath,$(CRAB_INSTALL)/elina/lib/
ifeq ($(UNAME),Linux)
    LDFLAGS += -Wl,--disable-new-dtags 
endif

all: domain analysis
# Ugly hack: add elina lib directory to rpath is not enough
	@echo "Compilation done. \nType \"export LD_LIBRARY_PATH=$(CRAB_INSTALL)/elina/lib\""

%.o:  $(CXX) $(CXXFLAGS) $< -c -o $@

domain: domain.o
	$(CXX) domain.o $(LDFLAGS) -o domain

analysis: analysis.o
	$(CXX) analysis.o $(LDFLAGS) -o analysis

clean:
	rm -f domain analysis  *.a *.o

.PHONY: crab_elina_install
crab_elina_install:
	git clone https://github.com/seahorn/crab.git crab -b dev
	mkdir -p crab/build
	cd crab/build \
	    && cmake -DCMAKE_INSTALL_PREFIX=$(CRAB_INSTALL) -DCRAB_BUILD_LIBS_SHARED=OFF \
	             -DCRAB_USE_ELINA=ON ../ \
	    && cmake --build . --target elina && cmake ../ \
	    && cmake --build . --target install
