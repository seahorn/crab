# Crab #

<img src="http://i.imgur.com/IDKhq5h.png" alt="crab logo" width=280 height=200 />

#About#

Crab is a language agnostic engine to perform static analysis of
programs based on
[abstract interpretation](https://en.wikipedia.org/wiki/Abstract_interpretation).

At its core, Crab is a collection of abstract domains and fixpoint
iterators built on the top of
[Ikos](http://ti.arc.nasa.gov/opensource/ikos/) (Inference Kernel for
Open Static Analyzers) developed by NASA Ames Research Center.

Crab has been designed to have two kind of users:

1.  Analysis/verification tools that want to compute invariants using
    abstract interpretation.

2.  Researchers on abstract interpretation who would like to
    experiment with new abstract domains and fixpoint iterators.

## Licenses ##

Ikos is distributed under NASA Open Source Agreement (NOSA)
Version 1.3 or later. Crab is distributed under MIT license.

See [Crab_LICENSE.txt](Crab_LICENSE.txt) and
[Ikos_LICENSE.pdf](Ikos_LICENSE.pdf) for details.

## Compilation ##

Crab is written in C++ and uses heavily the Boost library. You will
need:

- C++ compiler supporting c++11
- Boost and Gmp 

Then, just type:

    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=my_install_dir ../
    cmake --build . --target install 

If you want to run the tests add `-DDEVMODE=ON` option.

If you want to use the BOXES domain then add `-DUSE_LDD=ON` option.

If you want to use the Apron library domains then add `-DUSE_APRON=ON` option.

If you want to enable statistics about the analysis then add `-DUSE_STATS=ON` option.

If you want to enable logging for debugging purposes then add `-DUSE_LOG=ON` option.

By default, all these options are set to `OFF` (i.e., disabled).

## Crab input language and output ##

The main input of Crab is a Control Flow Graph (CFG) language that
interfaces with the abstract domains and iterators for the purpose of
generating invariants.

In its simplest form, the CFG consists only of:

- assume,
- havoc, 
- arithmetic and bitwise operations, and
- goto instructions

but it also supports other instructions such as

- load, store, pointer arithmetic, and function pointers
- function calls and returns

To support inter-procedural analysis, Crab can also take as input a
whole program represented as a bunch of functions. For Crab, a
function is a CFG (as described above) with a signature that consists
of the function name, its formal parameters, and return value.

The output of Crab is a map from CFG basic blocks to invariants
expressed in the underlying abstract domain.

## Examples ##

The tests directory contains some examples of how to build CFGs and
how to compute invariants using different abstract domains.

Important: the option `DEVMODE` must be enabled to compile all the
tests.

## How to integrate Crab in other verification tools ##

Check these projects:

- [Crab-Llvm](https://github.com/seahorn/crab-llvm) is a static
analyzer that infers invariants from LLVM-based languages using Crab.

- [SeaHorn](https://github.com/seahorn) is a verification framework
that uses Crab-Llvm to supply invariants to the back-end solvers.

## How to implement new fixpoint iterators ##

The new fixpoint iterator must follow this API:

    template< class NodeName, class AbstractValue >
    class forward_fixpoint_iterator {
     public:
	 
     virtual AbstractValue analyze(NodeName, AbstractValue) = 0;
    
     virtual void process_pre(NodeName, AbstractValue) = 0;
    
     virtual void process_post(NodeName, AbstractValue) = 0;
    
    }; 

## How to implement new abstract domains ##

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
    
### How to implement new numerical abstract domains ###

In addition to the previous API, for numerical domains it is also
required to implement the API described in `numerical_domains_api.hpp`:

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

    // forget
    virtual void operator-=(VariableName v) = 0;

      
This API assumes the manipulation of linear expressions and linear
constraints both defined in `linear_constraints.hpp` so it is good to be
familiar with.

For non-relational domains it is highly recommend to build on the top
of separate_domains which provides an efficient implementation of a
fast mergeable integer map based on patricia trees. This map can be
used to map variable names to abstract values. The implementation of
intervals and congruences use it.

## Check properties with Crab ##

The check of properties is done in two phases. During the first phase,
the invariants are computed for the whole program and stored for now
in memory (in the future it can be done in some external database). In
a second phase, the invariants (computed during the first phase) are
used to prove or disprove the property checks.

First, we need some terminology.

An _analyzer_ object inherits from a fixpoint iterator and keeps in
memory the cfg or call graph as well as the invariants computed by the
fixpoint iterator.

A _checker_ object is a wrapper for an analyzer object that contains
in addition a bunch of property checkers (see below). The checker
scans all the cfg blocks and checks the analyzer's invariants wrt the
properties of interest. A checker object is associated to a unique
analyzer. Since Crab keeps only in memory the invariants that hold
only at the entry of the blocks the checker also needs to propagate
the invariants to each program point. This is done actually via the
analyzer. A checker object can run multiple property checkers as long
as all are associated to the same analyzer.

A property of interest is defined via _property checkers_. The
property checker object scans each block statement wrt to the
analyzer's invariants. Thus, it is the one that actually does the
job. A property checker object is also associated with a unique
analyzer.

At some point, the analyzer(s) and checker should be stable and
fixed. Therefore, new properties should imply only to implement new
property checkers.

As a proof of concept, we have implemented so far some simple property
checkers for nullity, division by zero and user-definable assertions
(similar to C asserts). All the new property checkers should inherit
from a base property checker class called `BaseProperty.hpp`.

## TODO List ##

- No domains for supporting floating point computations.
- No domains for supporting machine arithmetic.
- ...

