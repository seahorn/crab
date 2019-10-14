# Building Crab with Docker and running tests #

```shell
docker build --build-arg UBUNTU=xenial --build-arg BUILD_TYPE=Release -t seahorn/crab_apron:xenial -f docker/crab.apron.Dockerfile .
docker run -v `pwd`:/host -it seahorn/crab_apron:xenial
```

This will automatically download all dependencies from a base image
and build Crab under `/crab/build`.

Crab's install directory is added to `PATH`.

Build arguments (required):
- UBUNTU: xenial, bionic
- BUILD_TYPE: Release, Debug

