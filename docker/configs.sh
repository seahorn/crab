#
# Single source of truth for Crab's Docker build configurations.
#
# Sourced by docker/build.sh and used (through it) by .github/workflows/*.
# Defines:
#   CRAB_CONFIGS          space-separated list of valid configuration names
#   crab_config <name>    sets FLAGS, TARGETS, TESTS and IMAGE for <name>
#
# Per-configuration variables:
#   FLAGS    extra -D options appended to the cmake configure line
#   TARGETS  external dependencies to build, IN CONFIGURE ORDER. The top-level
#            CMakeLists.txt aborts the configuration (return()) once per missing
#            dependency, so each target needs its own
#            `cmake --build . --target X && cmake ..` round. The order in which
#            CMake looks for them is ldd, pplite, apron, elina (see the comment
#            "PPLite should be built before Apron" in CMakeLists.txt).
#   TESTS    expected-results files under tests/, run by tests/run_tests.sh
#   IMAGE    default image name (docker/build.sh --image overrides it)
#
# To add a new configuration: add a case below and its name to CRAB_CONFIGS.
#

CRAB_CONFIGS="default apron elina pplite pplite-native"

crab_config() {
    case "$1" in
        default)
            # No external libraries.
            FLAGS=""
            TARGETS=""
            TESTS="expected_results.out"
            IMAGE="seahorn/crab"
            ;;
        apron)
            # Apron and ldd (boxes) libraries.
            FLAGS="-DCRAB_USE_LDD=ON -DCRAB_USE_APRON=ON"
            TARGETS="ldd apron"
            TESTS="expected_results.apron.out expected_results.boxes.out"
            IMAGE="seahorn/crab_apron_boxes"
            ;;
        elina)
            # Elina library.
            FLAGS="-DCRAB_USE_ELINA=ON"
            TARGETS="elina"
            TESTS="expected_results.elina.out"
            IMAGE="seahorn/crab_elina"
            ;;
        pplite)
            # PPLite library through the Apron interface.
            FLAGS="-DCRAB_USE_APRON=ON -DCRAB_USE_PPLITE=ON"
            TARGETS="pplite apron"
            TESTS="expected_results.pplite.out"
            IMAGE="seahorn/crab_pplite"
            ;;
        pplite-native)
            # PPLite library through the native interface.
            # CRAB_USE_PPLITE_NATIVE implies CRAB_USE_PPLITE (CMakeLists.txt).
            # The pplite-native target clones the wrapper headers into
            # include/crab/domains/pplite and is only defined while that
            # directory is absent, hence its entry in .dockerignore.
            FLAGS="-DCRAB_USE_PPLITE_NATIVE=ON"
            TARGETS="pplite pplite-native"
            TESTS="expected_results.pplite_native.out"
            IMAGE="seahorn/crab_pplite_native"
            ;;
        *)
            echo "error: unknown configuration '$1' (valid: $CRAB_CONFIGS)" >&2
            return 1
            ;;
    esac
}
