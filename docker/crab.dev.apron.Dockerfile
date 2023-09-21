#
# Dockerfile for Crab image with apron and boxes libraries.
#

# Pull base image.
FROM seahorn/buildpack-deps-crab:jammy

# Assume that docker-build is ran in the top-level Crab directory
COPY . /crab
# Re-create the build directory that might have been present in the source tree
RUN rm -rf /crab/build /crab/debug /crab/release && mkdir /crab/build
WORKDIR /crab/build

ARG BUILD_TYPE
# Build configuration.
RUN cmake -GNinja \
          -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
          -DCMAKE_INSTALL_PREFIX=run \
          -DCMAKE_CXX_COMPILER=g++-12 \	  
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
RUN /crab/tests/run_tests.sh /crab/tests/expected_results.boxes.out /crab/build

WORKDIR /crab

