# Crab #

<img src="http://i.imgur.com/IDKhq5h.png" alt="crab logo" width=280 height=200 />

#About#

Crab is a language agnostic framework to perform static analysis using
[abstract interpretation](https://en.wikipedia.org/wiki/Abstract_interpretation)
techniques.

At its core, Crab is a collection of abstract domains and fixpoint
iterators built on the top of
[Ikos](http://ti.arc.nasa.gov/opensource/ikos/) (Inference Kernel for
Open Static Analyzers) developed by NASA Ames Research Center.

Crab has been designed to have two kind of users:

1.  Analysis/verification tools that want to compute invariants using
    abstract interpretation.

2.  Researchers on abstract interpretation who would like to
    experiment with new abstract domains and fixpoint iterators.

## License ##

Ikos is distributed under NASA Open Source Agreement (NOSA)
Version 1.3 or later. Crab is distributed under MIT license. See
[Crab_LICENSE.txt](Crab_LICENSE.txt) and
[Ikos_LICENSE.pdf](Ikos_LICENSE.pdf) for details.

## Compilation ##

Crab is written in C++ and uses heavily the Boost library. You will
need:

- C++ compiler supporting c++11
- Boost and Gmp 

Then, just type:

    cmake -DDEVMODE=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=my_install_dir ../
    cmake --build . --target install 

If you want to use the BOXES domain then add `-DUSE_LDD=ON` option.

If you want to use the Apron library domains then add `-DUSE_APRON=ON` option.


## Input/Output ##

The input of Crab is a Control Flow Graph (CFG) language that
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

The output of Crab is a map from CFG basic blocks to invariants. The
format of the output invariants can be chosen by the client.

## Examples ##

The tests directory contains some examples of how to build CFGs and
how to compute invariants using different abstract domains.

Important: the option `DEVMODE` must be enabled to compile all the
tests.

## How to integrate Crab in other verification tools ##

Check these projects:

- [crab-llvm](https://github.com/seahorn/crab-llvm) is a standalone
application that computes invariants from LLVM-based languages using
crab.

- [SeaHorn](https://github.com/seahorn) is a verification framework
that uses crab-llvm to supply invariants to the back-end solvers.

## How to implement new fixpoint iterators ##

The new fixpoint iterator must follow this API:

    template< class NodeName, class AbstractValue >
    class forward_fixpoint_iterator {
     public:
     virtual AbstractValue analyze(NodeName, AbstractValue) = 0;
    
     virtual void process_pre(NodeName, AbstractValue) = 0;
    
     virtual void process_post(NodeName, AbstractValue) = 0;
    
     virtual ~forward_fixpoint_iterator() { }
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

    // Forget
    virtual void operator-=(VariableName v) = 0;

    virtual ~numerical_domain() { }
      
This API assumes the manipulation of linear expressions and linear
constraints both defined in `linear_constraints.hpp` so it is good to be
familiar with.

For non-relational domains it is highly recommend to build on the top
of separate_domains which provides an efficient implementation of a
fast mergeable integer map based on patricia trees. This map can be
used to map variable names to abstract values. The implementation of
intervals and congruences use it.

### How to add non-standard abstract operations ###

Some abstract domains may need some non-standard abstract operations
(e.g., array read/write, expand, etc). The signature of these
operations are described in `domain_traits.hpp`. Default
implementations are in `domain_traits_impl.hpp`. If an abstract domain
needs a different implementation from the default one then it should
implement the operation under the namespace `crab::domain_traits`. The
file `domain_traits.hpp` must also include the header file of the
abstract domain (if not already there).

##People##

* [Arie Gurfinkel](arieg.bitbucket.org)
* [Jorge Navas](http://ti.arc.nasa.gov/profile/jorge/)
* [Temesghen Kahsai](http://www.lememta.info/)
* Graeme Gange 
