#
# Dockerfile for building Crab dependencies for Ubuntu xenial and bionic.
#
# Arguments:
#  - UBUNTU:     xenial, bionic
#

ARG UBUNTU

# Pull base image.
FROM buildpack-deps:$UBUNTU

RUN apt-get update && \
    apt-get install -yqq software-properties-common && \
    add-apt-repository -y ppa:mhier/libboost-latest && \     
    apt-get update && \
    apt-get install -yqq cmake cmake-data g++-6 \
                         ninja-build libstdc++6 \
			 libboost1.74-dev \
                         libgmp-dev libmpfr-dev libflint-dev \
			 lcov ggcov 

