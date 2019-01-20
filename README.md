# Crab: A Language-Agnostic Library for Static Analysis #

<a href="https://travis-ci.org/seahorn/crab"><img src="https://travis-ci.org/seahorn/crab.svg?branch=master" title="Ubuntu 12.04 LTS 64bit, g++-5.0"/></a>

<img src="http://i.imgur.com/IDKhq5h.png" alt="crab logo" width=280 height=200 />

# Description #


Crab (CoRnucopia of ABstractions) performs static analysis of programs
based on [Abstract Interpretation](https://en.wikipedia.org/wiki/Abstract_interpretation).

Crab does not analyze directly a mainstream programming language such
as C, C++, or Java but instead it analyzes a simplified
Control-Flow-Graph (CFG) based language which is
language-independent. This allows Crab analyzing different programming
languages assuming a translator to the CFG-based language is
available.

Crab has been designed to have two kind of users:

1.  Analysis/verification tools that want to compute invariants using
    abstract interpretation.

2.  Researchers on abstract interpretation who would like to
    experiment with new abstract domains and fixpoint iterators.

In spite of its simple design, Crab can scale with large real programs
and its CFG-based language is rich enough to represent programs with
loops, functions, pointers, etc.

The foundations of Crab is based on a collection of abstract domains
and fixpoint iterators built on the top of
[Ikos](http://ti.arc.nasa.gov/opensource/ikos/) (Inference Kernel for
Open Static Analyzers) developed by NASA Ames Research Center.

# Crab architecture #

![Crab Architecture](https://github.com/seahorn/crab/blob/master/Crab_arch.jpg?raw=true "Crab Architecture")

# Requirements #

Crab is written in C++ and relies on the Boost library. The main
requirements are:

- Modern C++ compiler supporting c++11
- Boost
- GMP 
- MPFR (if `-DUSE_APRON=ON` or `-DUSE_ELINA=ON`)

In linux, you can install requirements typing the commands:

	sudo apt-get install libboost-all-dev libboost-program-options-dev
    sudo apt-get install libgmp-dev
    sudo apt-get install libmpfr-dev	

# Installation #

To install Crab, type:

	mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=_INSTALL_DIR_ ../
    cmake --build . --target install 

The Boxes/Apron/Elina domains require third-party libraries. To avoid
the burden to users who are not interested in those domains, the
installation of the libraries is optional.

- If you want to use the Boxes domain then add `-DUSE_LDD=ON` option.

- If you want to use the Apron library domains then add
  `-DUSE_APRON=ON` option.

- If you want to use the Elina library domains then add
  `-DUSE_ELINA=ON` option.

**Important:** Apron and Elina are currently not compatible so you
cannot enable `-DUSE_APRON=ON` and `-DUSE_ELINA=ON` at the same time. 

To use Elina on Linux, you will need to add `_INSTALL_DIR_/lib` in the
environment variable `LD_LIBRARY_PATH` if Elina is installed in a
non-standard directory:

    export LD_LIBRARY_PATH=_INSTALL_DIR_/lib
	
For instance, to install Crab with Boxes and Apron, type:

	mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=_INSTALL_DIR_ -DUSE_LDD=ON -DUSE_APRON=ON ../
	cmake --build . --target ldd && cmake ..
	cmake --build . --target apron && cmake ..	
    cmake --build . --target install 	

The `tests` directory contains many examples of how to build CFGs and
compute invariants using different abstract domains. To compile these tests
type:

	mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=_INSTALL_DIR_ -DUSE_LDD=ON -DUSE_APRON=ON -DENABLE_TESTS=ON ../
	cmake --build . --target ldd && cmake ..
	cmake --build . --target apron && cmake ..	
    cmake --build . --target install 	

and then, for instance, to run `test1`:

    build/test-bin/test1

# Usage #

To include Crab in your C++ application you need to:

- include the C++ header files located at the
`_INSTALL_DIR_/crab/include`, and
 
- link your application with the Crab libraries installed in
`_INSTALL_DIR_/crab/lib`.


If you compile with Boxes/Apron/Elina you need also to include
`_INSTALL_DIR_/xxx/include` and link with `_INSTALL_DIR_/xxx/lib`
where `xxx=apron|elina|ldd`.

Read [this](https://github.com/seahorn/crab/blob/master/external/Makefile) for a real Makefile.

## Example ## 

Assume we want to perform static analysis on the following C-like
program:

```c
	int i,x,y;
	i=0;
	x=1;
	y=0;
	while (i < 100) {
	   x=x+y;
	   y=y+1;
	   i=i+1;
	}	 
``` 

Next, we show a simplified version of the C++ code to build the
corresponding Crab CFG and run the analysis using the Zones domain.

**Note**: this code has been simplified for presentation purposes and
it might not compile like it is. Read [this](https://github.com/seahorn/crab/blob/master/external/analysis.cpp) for real example.

```c++
    // CFG-based language
    #include <crab/cfg/cfg.hpp>
    // Variable factory	
    #include <crab/cfg/var_factory.hpp>
    // Intra forward analyzer	
    #include <crab/analysis/fwd_analyzer.hpp>
    // Zones domain
    #include <crab/domains/split_dbm.hpp
    /* 
      To define a Control-Flow Graph (CFG) users need to define :
      (1) Type for a variable 
      (2) Type for a basic block label
      (3) Choose between integers or rationals (Crab cannot mix them)
    */
    // (1) A variable factory based on strings
    typedef cfg::var_factory_impl::str_variable_factory variable_factory_t;
    typedef typename variable_factory_t::varname_t varname_t;
    // (2) CFG basic block labels
    typedef std::string basic_block_label_t;
    // (3) CFG over integers
    typedef cfg::Cfg<basic_block_label_t, varname_t, z_number> z_cfg_t;
    // Convenient wrapper for a CFG
    typedef cfg:cfg_ref<z_cfg_t> z_cfg_ref_t;

    // Abstract domain: zones or difference-constraints domain
    typedef SplitDBM<z_number, varname_t> zones_domain_t;

    /* 
      Crab provides both intra- and inter-procedural analyses which 
      are parametric on the abstract domain: we choose an
      intra-procedural forward analysis with zones domain
    */
    typedef intra_fwd_analyzer<z_cfg_ref_t, zones_domain_t> intra_zones_analyzer_t;	

    int main (int argc, char**argv) {
       // Create variable factory. 
       // Important: only one variable factory should be used to build a CFG. 
       // Moreover, the variable factory should be alive while the CFG is in use.
       variable_factory_t vfac;	
       // Declare variables i,x, and y
       z_var i (vfac ["i"], INT_TYPE, 32);
       z_var x (vfac ["x"], INT_TYPE, 32);
       z_var y (vfac ["y"], INT_TYPE, 32);
       // Create an empty CFG marking "entry" and "exit" are the labels
       // for the entry and exit blocks.
       cfg_t cfg ("entry","ret");
       // Add blocks
       basic_block_t& entry = cfg.insert ("entry");
       basic_block_t& bb1   = cfg.insert ("bb1");
       basic_block_t& bb1_t = cfg.insert ("bb1_t");
       basic_block_t& bb1_f = cfg.insert ("bb1_f");
       basic_block_t& bb2   = cfg.insert ("bb2");
       basic_block_t& ret   = cfg.insert ("ret");
       // Add control flow 
       entry >> bb1; bb1 >> bb1_t; bb1 >> bb1_f;
       bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
       // Add statements
       entry.assign(i, 0);
       entry.assign(x, 1);
       entry.assign(y, 0);
       bb1_t.assume(i <= 99);
       bb1_f.assume(i >= 100);
       bb2.add(x,x,y);
       bb2.add(y,y,1);
       bb2.add(i,i,1);

       // Build an analyzer and run the zones domain
       zones_domain_t inv;  // initially top
       intra_zones_analyzer_t a (cfg, inv, ...);
       a.run();
       cout << "Invariants using " << zones_domain_t::getDomainName() << "\n";
	
       // Scan all CFG basic blocks and print the invariants that hold
       // at their entries
       for (auto &b : cfg) {
         auto inv = a[b.label()];
         cout << get_label_str(b.label()) << "=" << inv << "\n";
       }
       return 0;
    }
```

The Crab output of this program, showing the invariants that hold at
the entry of each basic block, should be something like this:

    Invariants using SplitDBM
	
	entry={}
	bb1={i -> [0, 100], x -> [1, +oo], y -> [0, 100], y-i<=0, y-x<=0, i-x<=0, i-y<=0}
    bb1_t={i -> [0, 100], x -> [1, +oo], y -> [0, 100], y-i<=0, y-x<=0, i-x<=0, i-y<=0}
    bb1_f={i -> [0, 100], x -> [1, +oo], y -> [0, 100], y-i<=0, y-x<=0, i-x<=0, i-y<=0}
    bb2={i -> [0, 99], x -> [1, +oo], y -> [0, 99], y-i<=0, y-x<=0, i-x<=0, i-y<=0}
	ret={i -> [100, 100], x -> [100, +oo], y -> [100, 100], y-i<=0, y-x<=0, i-x<=0, i-y<=0}

# Integration of Crab in other analysis tools #

Crab has been integrated in these static analysis tools:

- [Crab-Llvm](https://github.com/seahorn/crab-llvm) is a static
analyzer that infers invariants from LLVM-based languages using Crab.

- [SeaHorn](https://github.com/seahorn) is a verification framework
that uses Crab-Llvm to supply invariants to the back-end solvers.

# References #

- "Exploiting Sparsity in Difference-Bounds Matrices" [(PDF)](https://jorgenavas.github.io/papers/zones-SAS16.pdf) by G. Gange, Jorge A. Navas, P. Schachte, H. Sondergaard, and P. Stuckey. SAS'16.
- "An Abstract Domain of Uninterpreted Functions" [(PDF)](https://jorgenavas.github.io/papers/terms-vmcai16.pdf) by G. Gange, Jorge A. Navas, P. Schachte, H. Sondergaard, and P. Stuckey. VMCAI'16.
- "Signedness-Agnostic Program Analysis: Precise Integer Bounds for Low-Level Code" [(PDF)](https://jorgenavas.github.io/papers/wrapped-intervals-aplas12.pdf) by Jorge A. Navas, P. Schachte, H. Sondergaard, and P. Stuckey. APLAS'12.
- "Boxes: A Symbolic Abstract Domain of Boxes" [(PDF)](https://pdfs.semanticscholar.org/93da/8102c5ca512126d1a45ee81da1ab0b0fd47c.pdf) by A. Gurfinkel and S. Chaki. SAS'10.

