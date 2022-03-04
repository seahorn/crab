# Crab: A C++ Library for Building Program Static Analyses #

<img src="http://i.imgur.com/IDKhq5h.png" alt="crab logo" width=280 height=200 />

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


**If you use this library please cite this [paper](https://dblp.uni-trier.de/rec/conf/vstte/GurfinkelN21.html).**

---

# Description #

Crab (CoRnucopia of ABstractions) is a C++ library for building
program static analyses based on
[Abstract Interpretation](https://en.wikipedia.org/wiki/Abstract_interpretation).

Crab does not analyze directly a mainstream programming language such
as C/C++ or Java but instead it analyzes its own CFG-based
intermediate representation (CrabIR). This allows building analyses
for different programming languages assuming a translator to CrabIR is
available. In spite of its simple design, CrabIR is rich enough to
represent languages such as [LLVM](https://llvm.org/) bitcode.

Crab has been designed to have two kind of users:

1.  Program analysis/verification tools that want to use invariants
    computed by abstract interpretation.

2.  Researchers on abstract interpretation who would like to
    experiment with new abstract domains or fixpoint algorithms.

The available documentation can be found in
our [wiki](https://github.com/seahorn/crab/wiki/Home).

# Docker # 

A (nightly) pre-built version of Crab that runs all tests can be
obtained using Docker:


``` shell
docker pull seahorn/crab:bionic
docker run -v `pwd`:/host -it seahorn/crab:bionic
```

# Requirements for compiling from sources #

Crab is written in C++ and relies on the Boost library. The main
requirements are:

- C++11 compiler 
- Boost >= 1.65
- GMP 
- MPFR (if `-DCRAB_USE_APRON=ON` or `-DCRAB_USE_ELINA=ON`)

In linux, you can install requirements typing the commands:

	sudo apt-get install libboost-all-dev libboost-program-options-dev
    sudo apt-get install libgmp-dev
    sudo apt-get install libmpfr-dev	

# Building from sources and installation #

To install Crab, type:

	mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=_INSTALL_DIR_ ../
    cmake --build . --target install 

The `tests` directory contains many examples of how to build programs
written in CrabIR and compute invariants using different analyses and
abstract domains. To compile these tests type:

	mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=_INSTALL_DIR_ -DCRAB_ENABLE_TESTS=ON ../	
    cmake --build . --target install 	

and then, for instance, to run `test1`:

    build/test-bin/test1
    
## Building from sources with third-party abstract domain libraries ##

The [Boxes](https://github.com/seahorn/ldd)/[Apron](https://github.com/antoinemine/apron)/[Elina](https://github.com/eth-sri/ELINA) domains require third-party libraries. To avoid
the burden to users who are not interested in those domains, the
installation of the libraries is optional.

- If you want to use the Boxes domain then add `-DCRAB_USE_LDD=ON` option.

- If you want to use the Apron library domains then add
  `-DCRAB_USE_APRON=ON` option.

- If you want to use the Elina library domains then add
  `-DCRAB_USE_ELINA=ON` option.

**Important:** Apron and Elina are currently not compatible so you
cannot enable `-DCRAB_USE_APRON=ON` and `-DCRAB_USE_ELINA=ON` at the same time. 
	
For instance, to install Crab with Boxes and Apron, type:

	mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=_INSTALL_DIR_ -DCRAB_USE_LDD=ON -DCRAB_USE_APRON=ON ../
	cmake --build . --target ldd && cmake ..
	cmake --build . --target apron && cmake ..	
    cmake --build . --target install 	


# Using Crab in other C++ projects #

To include Crab in your C++ application you need to:

- Include the C++ header files located at the
`_INSTALL_DIR_/crab/include`, and
 
- Link your application with the Crab libraries installed in
`_INSTALL_DIR_/crab/lib`.

If you compile with Boxes/Apron/Elina you need also to include
`_INSTALL_DIR_/EXT/include` and link with `_INSTALL_DIR_/EXT/lib`
where `EXT=apron|elina|ldd`.

## CMake ## 

If your project uses `cmake` then you just need to add in your project's `CMakeLists.txt`:

```
add_subdirectory(crab)
include_directories(${CRAB_INCLUDE_DIRS})
```

And then link your executable with `${CRAB_LIBS}`

## Make ## 

If your project uses `make`, read this
sample [Makefile](https://github.com/seahorn/crab/blob/master/external/Makefile).
