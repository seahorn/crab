#
# Dockerfile for Crab binary
# produces package in /crab/build
# Arguments:
#  - UBUNTU:     trusty, xenial, bionic
#  - BUILD_TYPE: debug, release
#

ARG UBUNTU

# Pull base image.
FROM buildpack-deps:$UBUNTU

ARG BUILD_TYPE
RUN echo "Build type set to: $BUILD_TYPE" && \
     # Install deps.
    apt-get update && \
    apt-get install -yqq software-properties-common && \
    apt-get update && \
    apt-get install -yqq cmake cmake-data g++-5 \
                         ninja-build libstdc++5 \
                         libgmp-dev libmpfr-dev


WORKDIR /tmp/dockerutils

# Create a helper script that works as switch (VAL) { Key0 : Val0, ...}.
# This is to work around docker limitations and pass right correct flag to the
# python configuration script.
RUN echo '#!/bin/sh' > switch.sh && \ 
    echo 'VAL=$1;shift;while test $# -gt 0;do if [ "$1" = "$VAL" ];then echo $2;exit 0;fi;shift;shift;done' >> switch.sh && \
    chmod +x switch.sh && \
    /tmp/dockerutils/switch.sh $BUILD_TYPE Debug "debug" Release "rel" \
    > /tmp/dockerutils/dt_out.txt && \
    export BT=$(cat /tmp/dockerutils/dt_out.txt) && \
    export UB=$(lsb_release --a 2>&1 | cut -f2 | tail -n 1) && \
    echo "$UB"_"$BT" > /tmp/dockerutils/prefix.txt && \
    cat /tmp/dockerutils/prefix.txt && \
    mkdir -p /deps

WORKDIR /deps
RUN export PREFIX=$(cat /tmp/dockerutils/prefix.txt) && \
    export DEPS_BASE=$(echo https://github.com/seahorn/seahorn-ext-deps/releases/download/crab/"$PREFIX") && \
    curl -sSOL "$DEPS_BASE"_boost_1_68.tar.gz && \
    tar -xf "$PREFIX"_boost_1_68.tar.gz 

RUN cd / && rm -rf /crab && \
    git clone https://github.com/seahorn/crab crab --depth=10 ; \
    mkdir -p /crab/build
WORKDIR /crab/build

ARG BUILD_TYPE
# Build configuration.
RUN cmake -GNinja \
          -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
          -DBOOST_ROOT=/deps/boost \
          -DCMAKE_INSTALL_PREFIX=run \
          -DCMAKE_CXX_COMPILER=g++-5 \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
          -DUSE_LDD=ON \
          -DUSE_APRON=ON \
	  -DENABLE_TESTS=ON \
          ../ && \
    cmake --build . --target ldd  && cmake .. && \
    cmake --build . --target apron  && cmake .. && \
    cmake --build . --target install

ENV PATH "/crab/build/run/bin:$PATH"

WORKDIR /crab

