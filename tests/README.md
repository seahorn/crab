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
		
