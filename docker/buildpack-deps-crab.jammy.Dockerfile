#
# Dockerfile for building Crab dependencies for Ubuntu 22.04 jammy
#

# Pull base image.
FROM buildpack-deps:jammy

RUN apt-get update && \
    apt-get install -yqq software-properties-common && \
    apt-get install -yqq build-essential && \
    add-apt-repository -y ppa:mhier/libboost-latest && \     
    apt-get update && \
    apt-get install -yqq cmake cmake-data \
                         ninja-build \
			 g++-12 \
			 libboost1.74-dev libboost-program-options1.74-dev \
                         libgmp-dev libmpfr-dev libflint-dev 

