#------------------------------------------------------------------------------------#
#                                  INSTRUCTIONS
#------------------------------------------------------------------------------------#
#  0. Define variables BOOST_ROOT, GMP_DIR in Makefile.config and MPFR_DIR
#  1. Type `make crab_apron_ldd_install`
#  2. Type `make all`
#  3. Type `./domain` or `./analysis`
#------------------------------------------------------------------------------------#

include ../Makefile.config

MPFR_DIR=

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

CXXFLAGS += \
     -I $(MPFR_DIR)/include \
     -I $(CRAB_INSTALL)/ldd/include/ldd \
     -I $(CRAB_INSTALL)/apron/include 

LDFLAGS += $(LDD_LIBS) $(APRON_LIBS)
LDFLAGS += -L$(MPFR_DIR)/lib -lmpfr


all: domain analysis

%.o:  $(CXX) $(CXXFLAGS) $< -c -o $@

domain: domain.o
	$(CXX) domain.o $(LDFLAGS) -o domain

analysis: analysis.o
	$(CXX) analysis.o $(LDFLAGS) -o analysis

clean:
	rm -f domain analysis  *.a *.o

.PHONY: crab_apron_ldd_install
crab_apron_ldd_install:
	git clone https://github.com/seahorn/crab.git crab
	mkdir -p crab/build
	@echo "Installing crab+apron+ldd in ${CRAB_INSTALL} ... "
	cd crab/build \
	     && cmake -DCMAKE_INSTALL_PREFIX=$(CRAB_INSTALL) -DCRAB_BUILD_LIBS_SHARED=OFF \
	              -DCRAB_USE_LDD=ON -DCRAB_USE_APRON=ON ../ \
	     && cmake --build . --target ldd && cmake ../ \
	     && cmake --build . --target apron && cmake ../ \
	     && cmake --build . --target install

