#
# Dockerfile for Crab binary with apron and boxes libraries.
#
# produces package in /crab/build
# Arguments:
#  - UBUNTU:     xenial, bionic
#  - BUILD_TYPE: debug, release
#  - BRANCH
#

ARG UBUNTU

# Pull base image.
FROM buildpack-deps:$UBUNTU

ARG BUILD_TYPE
RUN echo "Build type set to: $BUILD_TYPE" && \
     # Install deps.
    apt-get update && \
    apt-get install -yqq software-properties-common && \
    add-apt-repository -y ppa:mhier/libboost-latest && \         
    apt-get update && \
    apt-get install -yqq cmake cmake-data g++-5 \
                         ninja-build libstdc++5 \
			 libboost1.68-dev \
                         libgmp-dev libmpfr-dev \
			 lcov ggcov 

WORKDIR /tmp/dockerutils

# Create a helper script that works as switch (VAL) { Key0 : Val0, ...}.
# This is to work around docker limitations and pass right correct flag to the
# python configuration script.
RUN echo '#!/bin/sh' > switch.sh && \ 
    echo 'VAL=$1;shift;while test $# -gt 0;do if [ "$1" = "$VAL" ];then echo $2;exit 0;fi;shift;shift;done' >> switch.sh && \
    chmod +x switch.sh && \
    /tmp/dockerutils/switch.sh $BUILD_TYPE Debug "debug" Release "rel" Coverage "rel" \
    > /tmp/dockerutils/dt_out.txt && \
    export BT=$(cat /tmp/dockerutils/dt_out.txt) && \
    export UB=$(lsb_release --a 2>&1 | cut -f2 | tail -n 1) && \
    echo "$UB"_"$BT" > /tmp/dockerutils/prefix.txt && \
    cat /tmp/dockerutils/prefix.txt 

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
          -DCMAKE_CXX_COMPILER=g++-5 \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
          -DCRAB_USE_LDD=ON \
          -DCRAB_USE_APRON=ON \
	  -DCRAB_ENABLE_TESTS=ON \
          ../ && \
    cmake --build . --target ldd  && cmake .. && \
    cmake --build . --target apron  && cmake .. && \
    cmake --build . --target install

# Run tests
RUN /crab/tests/run_tests.sh /crab/tests/expected_results.apron.out /crab/build

WORKDIR /crab

