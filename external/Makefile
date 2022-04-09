#------------------------------------------------------------------------------------#
#                                  INSTRUCTIONS
#------------------------------------------------------------------------------------#
#  0. Define variables BOOST_ROOT and GMP_DIR in Makefile.config
#  1. Type `make crab_install`
#  2. Type `make all`
#  3. Type `./domain` or `./analysis`
#------------------------------------------------------------------------------------#

include Makefile.config

all: domain analysis

%.o:  $(CXX) $(CXXFLAGS) $< -c -o $@

domain: domain.o
	$(CXX) domain.o $(LDFLAGS)  -o domain

analysis: analysis.o
	$(CXX) analysis.o $(LDFLAGS) -o analysis

clean:
	rm -f domain analysis *.a *.o
