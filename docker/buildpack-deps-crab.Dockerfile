#
# Dockerfile for building Crab dependencies.
#
# Arguments:
#  - UBUNTU:      jammy (default), bionic, xenial
#  - GCC_VERSION: version of g++ to install (12 by default)
#
# Usage:
#   docker build -t seahorn/buildpack-deps-crab:jammy \
#                -f docker/buildpack-deps-crab.Dockerfile .
#   docker build --build-arg UBUNTU=bionic --build-arg GCC_VERSION=6 \
#                -t seahorn/buildpack-deps-crab:bionic \
#                -f docker/buildpack-deps-crab.Dockerfile .
#

ARG UBUNTU=jammy

# Pull base image.
FROM buildpack-deps:${UBUNTU}

ARG GCC_VERSION=12

RUN apt-get update && \
    apt-get install -yqq software-properties-common && \
    apt-get install -yqq build-essential && \
    add-apt-repository -y ppa:mhier/libboost-latest && \
    apt-get update && \
    apt-get install -yqq cmake cmake-data \
                         ninja-build libstdc++6 \
                         g++-${GCC_VERSION} \
                         libboost1.74-dev libboost-program-options1.74-dev \
                         libgmp-dev libmpfr-dev libflint-dev \
                         lcov
