# Building Docker image for Crab dependencies #


```shell
$ docker build --build-arg UBUNTU=bionic -t seahorn/buildpack-deps-crab:bionic -f docker/buildpack-deps-crab.Dockerfile .
$ docker push seahorn/buildpack-deps-crab:bionic
```


# Building Docker image for Crab and running tests #

```shell
docker build --build-arg UBUNTU=bionic --build-arg BUILD_TYPE=Release -t seahorn/crab:bionic -f docker/crab.Dockerfile .
docker run -v `pwd`:/host -it seahorn/crab:bionic
```

This will automatically download all dependencies from a base image
and build Crab under `/crab/build`.

Crab's install directory is added to `PATH`.

Build arguments (required):
- UBUNTU: xenial, bionic
- BUILD_TYPE: Release, Debug

