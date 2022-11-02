# Crab: A C++ Library for Building Program Static Analyses #

Crab is a C++ library for building program static analyses based on
[Abstract Interpretation](https://en.wikipedia.org/wiki/Abstract_interpretation). Crab provides
a rich set of abstract domains, Kleene-based fixpoint solvers, as well
as different analyses such as dataflow, inter-procedural and
backward. The design of Crab is quite modular so that it is easy to
plugin new abstract domains and solvers or build new analyses.

Crab abstract domains can reason about memory contents, C-like arrays
and numerical properties. Crab uses efficient implementations of
popular numerical domains such as [Zones and
Octagons](https://dl.acm.org/doi/abs/10.1145/3457885) and novel
domains to reason, for instance, about [symbolic
terms](https://dl.acm.org/doi/10.1007/978-3-662-49122-5_4) (aka
uninterpreted functions). Crab also implements popular non-relational
domains such as interval or congruences using [efficient environment
maps](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.37.5452),
and allows the combination of arbitrary domains via standard reduced
product constructions. Crab also provides non-convex domains such as
specialized disjunctive intervals called
[Boxes](https://link.springer.com/chapter/10.1007/978-3-642-15769-1_18)
based on linear decision diagrams and a more general value
partitioning strategy that lifts an arbitrary domain to an
over-approximation of its disjunction completion. In addition to these
domains, all developed by Crab authors, the Crab library integrates
popular abstract domain libraries such as
[Apron](https://github.com/antoinemine/apron),
[Elina](http://elina.ethz.ch/), and [PPLite](https://github.com/ezaffanella/PPLite).

Crab provides the state-of-the-art [interleaved fixpoint
solver](https://link.springer.com/chapter/10.1007/978-3-642-38856-9_4)
that uses Bourdoncle's [Weak Topological
Ordering](https://link.springer.com/chapter/10.1007/BFb0039704) to
select the set of widening points. To mitigate precision losses during
widening, Crab implements some popular techniques such as widening
with thresholds and [lookahead
widening](https://link.springer.com/chapter/10.1007/11817963_41).

Crab provides two different implementations of inter-procedural
analyses: a top-down with memoization inter-procedural analysis with
support for recursive calls, and a hybrid of bottom-up + top down
analysis. Last but not least, Crab also implements a more experimental
backward analysis that can be used to compute necessary preconditions
and/or reduce the number of false alarms.

Crab does not analyze directly a mainstream programming language but
instead it analyzes its own CFG-based intermediate representation
called
[CrabIR](https://link.springer.com/chapter/10.1007/978-3-030-95561-8_8).
CrabIR is three-address code and it is strongly typed. In CrabIR,
control flow is defined via non-deterministic goto instructions. Apart
from standard boolean and arithmetic operations, CrabIR provides
special assume and assert statements. The former can be used to refine
the control flow and the latter provides a simple mechanishm to check
for user-defined properties. In spite of its simple design, CrabIR is
rich enough to represent languages such as
[LLVM](https://github.com/seahorn/clam).

**Crab is actively under development. If you find a bug please open an
Github issue. Pull requests with new features are very welcome.  The
available documentation can be found in our
[wiki](https://github.com/seahorn/crab/wiki/Home). If you use this library please cite this [paper](https://dblp.uni-trier.de/rec/conf/vstte/GurfinkelN21.html).**

<br/>

<table>
  <tr>
    <th>Windows</th><th>Ubuntu</th><th>OS X</th><th>Coverage</th>
  </tr>
    <td>TBD</td>
    <td> <a href="https://github.com/seahorn/crab/actions"><img src="https://github.com/seahorn/crab/workflows/CI/badge.svg?branch=master" title="Ubuntu 18.04 LTS 64bit, g++-6.0"/></a> </td>
    <td>TBD</td>
    <td><a href="https://codecov.io/gh/seahorn/crab"><img src="https://codecov.io/gh/seahorn/crab/branch/master/graph/badge.svg" /></a></td>
  </tr>
</table>

# Docker # 

A (nightly) pre-built version of Crab that runs all tests can be
obtained using Docker:


``` shell
docker pull seahorn/crab:bionic
docker run -v `pwd`:/host -it seahorn/crab:bionic
```

# Requirements #

Crab is written in C++ and relies on the Boost library. The main
requirements are:

- C++11 compiler 
- Boost >= 1.65
- GMP 
- MPFR (if `-DCRAB_USE_APRON=ON` or `-DCRAB_USE_ELINA=ON`)
- FLINT (only if `-DCRAB_USE_PPLITE=ON`) 

In linux, you can install requirements typing the commands:

	sudo apt-get install libboost-all-dev libboost-program-options-dev
    sudo apt-get install libgmp-dev
    sudo apt-get install libmpfr-dev	
	sudo apt-get install libflint-dev

# Compilation and Installation #

To install Crab, type:

     1. mkdir build && cd build
     2. cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR ../
     3. cmake --build . --target install 

The `tests` directory contains many examples of how to build programs
written in CrabIR and compute invariants using different analyses and
abstract domains. To compile these tests add option `-DCRAB_ENABLE_TESTS=ON` to line 2.

and then, for instance, to run `test1`:

    build/test-bin/test1
    
## Include third-party abstract domain libraries ##

The [Boxes](https://github.com/seahorn/ldd)/[Apron](https://github.com/antoinemine/apron)/[Elina](https://github.com/eth-sri/ELINA)/[PPLite](https://github.com/ezaffanella/PPLite) domains require third-party libraries. To avoid
the burden to users who are not interested in those domains, the
installation of the libraries is optional.

- If you want to use the Boxes domain then add `-DCRAB_USE_LDD=ON` option.

- If you want to use the Apron library domains then add
  `-DCRAB_USE_APRON=ON` option.

- If you want to use the Elina library domains then add
  `-DCRAB_USE_ELINA=ON` option.

- If you want to use the PPLite library domains then add
  `-DCRAB_USE_PPLITE=ON` option.


**Important:** Apron and Elina are currently not compatible so you
cannot enable `-DCRAB_USE_APRON=ON` and `-DCRAB_USE_ELINA=ON` at the same time. 
	
For instance, to install Crab with Boxes and Apron, type:

     1. mkdir build && cd build
     2. cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DCRAB_USE_LDD=ON -DCRAB_USE_APRON=ON ../
     3. cmake --build . --target ldd && cmake ..
     4. cmake --build . --target apron && cmake ..	
     5. cmake --build . --target install 	

Lines 3 and 4 will download, compile and install the Boxes and Apron domains, respectively. Replace `apron` at line 4 with `elina` or `pplite` if you want to use Elina or PPLite instead. If you have already compiled and installed these libraries in your machine then skip commands at line 3 and 4 and add the following options at line 2.

- For Apron: `-DAPRON_ROOT=$APRON_INSTALL_DIR`
- For Elina: `-DELINA_ROOT=$ELINA_INSTALL_DIR`
- For Boxes: `-DCUDD_ROOT=$CUDD_INSTALL_DIR -DLDD_ROOT=$LDD_INSTALL_DIR`
- For PPLite: `-DPPLITE_ROOT=$PPLITE_INSTALL_DIR -DFLINT_ROOT=$FLINT_INSTALL_DIR`

# Using Crab library in other C++ projects #

To include Crab in your C++ application you need to:

- Include the C++ header files located at the
`$INSTALL_DIR/crab/include`, and
 
- Link your application with the Crab libraries installed in
`$INSTALL_DIR/crab/lib`.

If you compile with Boxes/Apron/Elina/PPLite you need also to include
`$INSTALL_DIR/EXT/include` and link with `$INSTALL_DIR/EXT/lib`
where `EXT=apron|elina|ldd|pplite`.

## CMake ## 

If your project uses `cmake` then you just need to add in your project's `CMakeLists.txt`:

```
add_subdirectory(crab)
include_directories(${CRAB_INCLUDE_DIRS})
```

And then link your executable with `${CRAB_LIBS}`

## Make ## 

If your project uses `make`, read this
sample [Makefile](https://github.com/seahorn/crab/blob/master/make/Makefile).
