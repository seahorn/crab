# Building Crab's dependencies with Docker #

The only external dependency that needs to be prebuilt is Boost.  This
is because old Linux distributions (trusty) installs Boost 1.55 which
is too old for Crab.

Boost 1.68 can be prebuilt using docker:


```shell
cd deps
docker build --build-arg UBUNTU=xenial --build-arg BUILD_TYPE=Release -t crab_deps_xenial_rel .
docker run -v $(pwd):/host -it crab_deps_xenial_rel
```

This will automatically create a `boost.tar` file in the current working directory.

For all the dependencies, the possible build arguments are:
- UBUNTU: trusty, xenial, bionic
- BUILD_TYPE: Release, Debug

Note that both `UBUNTU` and `BUILD_TYPE` are required arguments.
