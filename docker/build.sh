#!/usr/bin/env bash
#
# Build (and optionally push) a Crab docker image.
#
# The configuration is looked up in docker/configs.sh and passed to
# docker/crab.Dockerfile as build arguments.
#
# Examples:
#   docker/build.sh default
#   docker/build.sh apron --build-type Debug
#   docker/build.sh elina --base-tag bionic --cxx g++-6
#   docker/build.sh apron --image-tag nightly --push
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# shellcheck source=configs.sh
. "$SCRIPT_DIR/configs.sh"

usage() {
    cat <<EOF
Usage: $0 CONFIG [options]

  CONFIG is one of: $CRAB_CONFIGS

Options:
  --build-type TYPE   Release (default), Debug, Coverage
  --jobs N            parallel compilation jobs (default: as many as cores)
  --base-tag TAG      tag of seahorn/buildpack-deps-crab (default: jammy)
  --cxx COMPILER      C++ compiler in the base image (default: g++-12)
  --image NAME        image name, overrides the one from docker/configs.sh
  --image-tag TAG     image tag (default: latest)
  --push              docker push the image after a successful build
  --print-image       print the image name:tag and exit without building
  -h, --help          print this message
EOF
}

if [ $# -eq 0 ]; then
    usage >&2
    exit 1
fi

CONFIG=""
BUILD_TYPE="Release"
BUILD_JOBS=""
BASE_TAG="jammy"
CXX="g++-12"
IMAGE_OVERRIDE=""
IMAGE_TAG="latest"
PUSH=false
PRINT_IMAGE=false

while [ $# -gt 0 ]; do
    case "$1" in
        --build-type)  BUILD_TYPE="$2"; shift 2 ;;
        --jobs)        BUILD_JOBS="$2"; shift 2 ;;
        --base-tag)    BASE_TAG="$2"; shift 2 ;;
        --cxx)         CXX="$2"; shift 2 ;;
        --image)       IMAGE_OVERRIDE="$2"; shift 2 ;;
        --image-tag)   IMAGE_TAG="$2"; shift 2 ;;
        --push)        PUSH=true; shift ;;
        --print-image) PRINT_IMAGE=true; shift ;;
        -h|--help)     usage; exit 0 ;;
        -*)            echo "error: unknown option '$1'" >&2; usage >&2; exit 1 ;;
        *)
            if [ -n "$CONFIG" ]; then
                echo "error: unexpected argument '$1'" >&2
                usage >&2
                exit 1
            fi
            CONFIG="$1"; shift ;;
    esac
done

if [ -z "$CONFIG" ]; then
    echo "error: no configuration given" >&2
    usage >&2
    exit 1
fi

# Sets FLAGS, TARGETS, TESTS and IMAGE, or fails if CONFIG is unknown.
crab_config "$CONFIG"

if [ -n "$IMAGE_OVERRIDE" ]; then
    IMAGE="$IMAGE_OVERRIDE"
fi
IMAGE_REF="$IMAGE:$IMAGE_TAG"

if $PRINT_IMAGE; then
    echo "$IMAGE_REF"
    exit 0
fi

# The build context is the top-level Crab directory (see .dockerignore).
docker build \
    --build-arg BASE_TAG="$BASE_TAG" \
    --build-arg CXX="$CXX" \
    --build-arg BUILD_TYPE="$BUILD_TYPE" \
    --build-arg BUILD_JOBS="$BUILD_JOBS" \
    --build-arg CRAB_FLAGS="$FLAGS" \
    --build-arg EXT_TARGETS="$TARGETS" \
    --build-arg TESTS="$TESTS" \
    -t "$IMAGE_REF" \
    -f "$SCRIPT_DIR/crab.Dockerfile" \
    "$REPO_ROOT"

if $PUSH; then
    docker push "$IMAGE_REF"
fi
