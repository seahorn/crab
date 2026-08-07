#
# Parameterized Dockerfile for all Crab configurations.
#
# It builds Crab, installs it under /crab/build/run and runs the test suites.
# The configuration (which external libraries are enabled, which external
# targets must be built, which tests are run) is entirely given by build
# arguments. Do not add per-configuration Dockerfiles: add a new entry to
# docker/configs.sh instead and build it with docker/build.sh.
#
# Arguments:
#  - BASE_TAG:    tag of seahorn/buildpack-deps-crab (jammy by default)
#  - CXX:         C++ compiler available in the base image (g++-12 by default)
#  - BUILD_TYPE:  Release (default), Debug, Coverage
#  - CRAB_FLAGS:  extra cmake options, e.g. "-DCRAB_USE_APRON=ON"
#  - EXT_TARGETS: external dependencies to build, in configure order
#  - TESTS:       expected-results files under /crab/tests
#  - BUILD_JOBS:  parallel compilation jobs (empty means as many as cores).
#                 Useful on memory-constrained machines: the default can run the
#                 build out of memory and get the compiler killed.
#
# Usage (prefer docker/build.sh, which fills these in from docker/configs.sh):
#
#   docker build --build-arg CRAB_FLAGS="-DCRAB_USE_ELINA=ON" \
#                --build-arg EXT_TARGETS="elina" \
#                --build-arg TESTS="expected_results.elina.out" \
#                -t seahorn/crab_elina -f docker/crab.Dockerfile .
#

ARG BASE_TAG=jammy

# Pull base image.
FROM seahorn/buildpack-deps-crab:${BASE_TAG}

ARG CXX=g++-12
ARG BUILD_TYPE=Release
ARG CRAB_FLAGS=
ARG EXT_TARGETS=
ARG TESTS=expected_results.out
ARG BUILD_JOBS=

# Assume that docker-build is ran in the top-level Crab directory
COPY . /crab
# Re-create the build directory that might have been present in the source tree
RUN rm -rf /crab/build && mkdir /crab/build
WORKDIR /crab/build

# Build configuration.
#
# Each external dependency needs its own build+reconfigure round: the top-level
# CMakeLists.txt returns early from the configuration as soon as it finds a
# missing dependency, so one `cmake ..` per target in EXT_TARGETS is required.
# EXT_TARGETS must therefore list the targets in the order in which CMake looks
# for them (ldd, pplite, apron, elina).
RUN set -ex; \
    cmake -GNinja \
          -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
          -DCMAKE_INSTALL_PREFIX=run \
          -DCMAKE_CXX_COMPILER=$CXX \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
          -DCRAB_ENABLE_TESTS=ON \
          $CRAB_FLAGS \
          ../; \
    for target in $EXT_TARGETS; do \
        cmake --build . --target $target ${BUILD_JOBS:+-j $BUILD_JOBS}; \
        cmake ..; \
    done; \
    cmake --build . --target install ${BUILD_JOBS:+-j $BUILD_JOBS}

# To find elina dynamic libraries. ENV cannot be made conditional, but the entry
# is harmless for the configurations that do not build elina.
ENV LD_LIBRARY_PATH="/crab/build/run/elina/lib:$LD_LIBRARY_PATH"

# Run tests
RUN set -ex; \
    for expected in $TESTS; do \
        /crab/tests/run_tests.sh /crab/tests/$expected /crab/build; \
    done

WORKDIR /crab
