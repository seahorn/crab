#
# Dockerfile for Crab image with elina library.
#
# produces package in /crab/build
# Arguments:
#  - UBUNTU:     xenial, bionic (default)
#  - BUILD_TYPE: Debug, Release (default)
#  - BRANCH:     master (default)
#

ARG UBUNTU

# Pull base image.
FROM seahorn/buildpack-deps-crab:$UBUNTU

ARG BRANCH=master
RUN cd / && rm -rf /crab && \
    git clone -b $BRANCH https://github.com/seahorn/crab crab --depth=10 ; \
    mkdir -p /crab/build
WORKDIR /crab/build

ARG BUILD_TYPE
# Build configuration.
RUN cmake -GNinja \
          -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
          -DCMAKE_INSTALL_PREFIX=run \
          -DCMAKE_CXX_COMPILER=g++-6 \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
          -DCRAB_USE_ELINA=ON \
	  -DCRAB_ENABLE_TESTS=ON \
          ../ && \
    cmake --build . --target elina  && cmake .. && \
    cmake --build . --target install

# To find elina dynamic libraries
ENV LD_LIBRARY_PATH "/crab/build/run/elina/lib:$LD_LIBRARY_PATH"

# Run tests
RUN /crab/tests/run_tests.sh /crab/tests/expected_results.elina.out /crab/build

WORKDIR /crab

