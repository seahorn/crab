# (Private) Ikos-core #

This project contains the core of Ikos, a static analyzer based on
abstract interpretation developed at NASA Ames Research Center. This
core consists of an implementation of a fixpoint algorithm and several
abstract numerical domains.

Ikos-core also provides a CFG that interfaces with IKOS as well as
adaptors to execute backward and forward analyses on it. This
repository has been extended with more abstract domains as well as
capabilities to perform abstract interpretation on concurrent systems.

**NOTES**: This repository is intended to be PRIVATE. A public version
  is available [here](https://github.com/seahorn/ikos-core). Currently
  this private repository has more functionality that the public
  one. Eventually, all the functionality of this repository should be
  moved to the public one.


# License #

Some of the software is released under the terms and conditions of the
NASA Open Source Agreement (NOSA) Version 1.3 or later:

- `include\algorithms`
- `include\common`
- `include\domains` except:
    - `dbm`
    - `term`
    - `array_graph`
	- `array_smashing` 
- `iterators`

# Prerequisites #

- The C++ compiler must support c++11
- Boost and gmp 

# Installation #

    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=my_install_dir DEVMODE=ON ../
	cmake --build . --target install 

# Usage #

The tests directory contains some examples of how to build CFGs and
how to compute invariants using different abstract domains.

*Important* : the option `DEVMODE` must be enabled to compile all the
 tests.

`
build/tests/prog-1
`

# How to add a non-standard abstract operations #

Some abstract operations may need some non-standard abstract operations (e.g., array read/write, expand, etc). The signature of these operations are described in `domain_traits.hpp`. Default implementations are in `domain_traits_impl.hpp`. If an abstract domain needs a different implementation from the default one then it should implement the operation under the namespace `ikos::domain_traits`. The file `domain_traits.hpp`  must also include the header file of the abstract domain (if not already there).

# Notes for implementing new abstract domains #

The main task is to implement the following API required by the
fixpoint algorithm:
  
    static AbsDomain top();
    
    static AbsDomain bottom();
    
    bool is_bottom ();

    bool is_top ();

    // Less or equal
    bool operator<=(AbsDomain o);

    // Join
    AbsDomain operator|(AbsDomain o);

    // Meet
    AbsDomain operator&(AbsDomain o);

    // Widening
    AbsDomain operator||(AbsDomain o);

    // Narrowing 
    AbsDomain operator&&(AbsDomain o);
    
# Notes for implementing new numerical abstract domains #

In addition to the previous API, for numerical domains it is also required to implement the API described in
numerical_domains_api.hpp:

    typedef linear_expression< Number, VariableName > linear_expression_t;
    typedef linear_constraint< Number, VariableName > linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
  
    // x = y op z
    virtual void apply(operation_t op, VariableName x, VariableName y, VariableName z) = 0; 

    // x = y op k
    virtual void apply(operation_t op, VariableName x, VariableName y, Number k) = 0; 

    // x = e
    virtual void assign(VariableName x, linear_expression_t e) = 0; 

    // assume (c);
    virtual void operator+=(linear_constraint_system_t csts) = 0;

    // Forget
    virtual void operator-=(VariableName v) = 0;

    virtual ~numerical_domain() { }
      
This API assumes the manipulation of linear expressions and linear
constraints both defined in linear_constraints.hpp so it is good to be
familiar with.

For non-relational domains it is highly recommend to build on the top
of separate_domains which provides an efficient implementation of a
fast mergeable integer map based on patricia trees. This map can be
used to map variable names to abstract values. The implementation of
intervals and congruences use it.

#People#

* [Arie Gurfinkel](arieg.bitbucket.org)
* [Jorge Navas](http://ti.arc.nasa.gov/profile/jorge/)
* [Temesghen Kahsai](http://www.lememta.info/)
* Graeme Gange (The University of Melbourne)