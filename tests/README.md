# Add a test for a new abstract domain (for developers) #

1. Include the header file in `crab_dom.hpp`:

        #include <crab/domains/wrapped_interval_domain.hpp>

2. Create an instantiation of the `run` function for the new domain:

        Z_RUNNER(crab::domain_impl::z_wrapped_interval_domain_t)

The `run` function should be called from your test (see step 3).

If the domain is defined over rational instead of integers then use

	    Q_RUNNER(...)
		
3. Add your test `wrapped_interval_domain_test.cc`. You can use
   `domains/test1.cc` as a template.

4. Include the directory where your test is located in
   `CMakeLists.txt`. For instance, if your test is under a
   subdirectory called `domains/wrapint`:

        AddTestDir(domains)        
		...
        AddTestDir (domains/wrapint) # directory where your test is located

Of course, if you put your test in the rest of the above directories
you can skip step 4.

**IMPORTANT:** Tests are only compiled if option `-DCRAB_ENABLE_TESTS=ON`
is enabled in the `cmake` command.

# Run the tests (for developers) #

There is no `make test` target. Tests are checked by running every
compiled binary and diffing its output against an expected-output file.

1. Build all the test binaries, e.g.:

        cmake --build build -j4

2. Run the suite with `run_tests.sh EXPECTED_OUTPUT BUILD_DIR`:

        ./tests/run_tests.sh tests/expected_results.out build

   It runs every `BUILD_DIR/test-bin/*`, concatenates their stdout, and
   `diff`s the result against `EXPECTED_OUTPUT` (ignoring `=== ...`
   header lines and `CRAB WARNING:*` lines). A test is considered
   correct iff its output matches the expected file: there are no
   asserts.

The expected-output files are:

  - `expected_results.out` — the default suite. Binaries whose name
    contains `apron`, `pplite`, `elina`, or `boxes` are **excluded**.
  - `expected_results.apron.out`, `expected_results.pplite.out`,
    `expected_results.elina.out`, `expected_results.boxes.out` — one per
    external library. For these, only binaries whose name contains the
    corresponding token are run.

**IMPORTANT:** When you add a test, or change the output of an existing
one, you must regenerate the matching `expected_results.*` file(s) or the
diff will fail. A test binary with no library token in its name (e.g.
`unittests-zones`) only affects `expected_results.out`.
		
