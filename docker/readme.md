# Building Docker image for Crab dependencies #

```shell
$ docker build -t seahorn/buildpack-deps-crab:jammy -f docker/buildpack-deps-crab.Dockerfile .
$ docker push seahorn/buildpack-deps-crab:jammy
```

The base image can be built for other Ubuntu releases with the `UBUNTU` and
`GCC_VERSION` build arguments:

```shell
$ docker build --build-arg UBUNTU=bionic --build-arg GCC_VERSION=6 \
               -t seahorn/buildpack-deps-crab:bionic -f docker/buildpack-deps-crab.Dockerfile .
```

# Building Docker image for Crab and running tests #

All Crab images are built from the single parameterized
`docker/crab.Dockerfile`. The configurations (which external libraries are
enabled, which external targets are built and which tests are run) live in
`docker/configs.sh` and are selected by name:

| Configuration   | Libraries                          | Image                          |
|-----------------|------------------------------------|--------------------------------|
| `default`       | none                               | `seahorn/crab`                 |
| `apron`         | Apron + Ldd (boxes)                | `seahorn/crab_apron_boxes`     |
| `elina`         | Elina                              | `seahorn/crab_elina`           |
| `pplite`        | PPLite through the Apron interface | `seahorn/crab_pplite`          |
| `pplite-native` | PPLite through the native interface| `seahorn/crab_pplite_native`   |

Use `docker/build.sh`, which fills in the build arguments from
`docker/configs.sh`:

```shell
$ docker/build.sh default
$ docker run -v `pwd`:/host -it seahorn/crab:latest
```

This will automatically download all dependencies from a base image, build Crab
under `/crab/build` and run the test suites of the configuration.

Options (see `docker/build.sh --help`):

- `--build-type TYPE`: `Release` (default), `Debug`, `Coverage`
- `--jobs N`: parallel compilation jobs (default: as many as cores). Lower it if
  the compiler gets killed on a memory-constrained machine
- `--base-tag TAG`: tag of `seahorn/buildpack-deps-crab` (default: `jammy`)
- `--cxx COMPILER`: C++ compiler of the base image (default: `g++-12`)
- `--image NAME`, `--image-tag TAG`: override the image name and tag
- `--push`: push the image after a successful build

```shell
$ docker/build.sh apron --build-type Debug
$ docker/build.sh elina --base-tag bionic --cxx g++-6
```

# Adding a new configuration #

Add a case to `crab_config` in `docker/configs.sh` and its name to
`CRAB_CONFIGS`. No new Dockerfile is needed. The GitHub workflows build every
name listed in their matrix, so add it there too.
