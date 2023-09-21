# Building Docker image for Crab dependencies #


```shell
$ docker build -t seahorn/buildpack-deps-crab:bionic -f docker/buildpack-deps-crab.Dockerfile .
$ docker push seahorn/buildpack-deps-crab:bionic
```
or 

```shell
$ docker build -t seahorn/buildpack-deps-crab:bionic -f docker/buildpack-deps-crab.jammy.Dockerfile .
$ docker push seahorn/buildpack-deps-crab:jammy
```

# Building Docker image for Crab and running tests #

```shell
docker build --build-arg BUILD_TYPE=Release -t seahorn/crab -f docker/crab.Dockerfile .
docker run -v `pwd`:/host -it seahorn/crab:latest
```

This will automatically download all dependencies from a base image
and build Crab under `/crab/build`.

Crab's install directory is added to `PATH`.

Build arguments (required):
- BUILD_TYPE: Release, Debug

