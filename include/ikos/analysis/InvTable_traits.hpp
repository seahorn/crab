#ifndef INVARIANT_TABLE_VALUES_TRAITS_HPP
#define INVARIANT_TABLE_VALUES_TRAITS_HPP

#include <type_traits>
#include <ikos/common/types.hpp>

/* 
   This class manipulates the values of the invariant table entries.
   Invariant tables are how ikos-core communicates with other
   clients. Clients can choose the type of the entry values (i.e.,
   invariants).
 */

namespace analyzer {

   using namespace ikos;

   template<typename AbsDomain, typename InvTableTy>
   class inv_tbl_traits {
    public:
     typedef AbsDomain  type_from;
     typedef InvTableTy type_to;

     //! Convert an abstract domain to the format stored in the
     //  invariant tables.
     //  FIXME: should be `const AbsDomain&` but we need to fix first
     //  constness of the abstract domain methods.
     static InvTableTy marshall (AbsDomain& abs) {
       if (!std::is_same<AbsDomain,InvTableTy>::value)
         IKOS_ERROR ("inv_tbl_traits expect two parameters of same type");
       return abs;
     }
     //! Convert from the format stored in the invariant table to an
     //! abstract domain.
     static AbsDomain unmarshall (InvTableTy inv) {
       if (!std::is_same<AbsDomain,InvTableTy>::value)
         IKOS_ERROR ("inv_tbl_traits expect two parameters of same type");
       return inv;
     }
     //! Return true 
     static InvTableTy top () {
       if (!std::is_same<AbsDomain,InvTableTy>::value)
         IKOS_ERROR ("inv_tbl_traits expect two parameters of same type");
       return AbsDomain::top ();
     }
     //! Return false
     static InvTableTy bot () {
       if (!std::is_same<AbsDomain,InvTableTy>::value)
         IKOS_ERROR ("inv_tbl_traits expect two parameters of same type");
       return AbsDomain::bottom ();
     }
   };

   // conjunctive linear constraints
   template<typename AbsDomain>
   class inv_tbl_traits <AbsDomain,
                         typename AbsDomain::linear_constraint_system_t>  {

     typedef typename AbsDomain::linear_expression_t lin_exp_t;

    public:

     typedef typename AbsDomain::linear_constraint_t lin_cst_t;
     typedef typename AbsDomain::linear_constraint_system_t lin_cst_sys_t;

     typedef AbsDomain  type_from;
     typedef lin_cst_sys_t type_to;

    private:

     static lin_cst_t mk_true() {
       return lin_cst_t (lin_exp_t (1) == lin_exp_t (1)); 
     }

     static lin_cst_t mk_false() {
       return lin_cst_t (lin_exp_t (1) == lin_exp_t (0)); 
     }

    public:
     
     static lin_cst_sys_t marshall (AbsDomain& abs) {
       return abs.to_linear_constraint_system ();
     }

     static AbsDomain unmarshall (lin_cst_sys_t csts) {
       AbsDomain inv = AbsDomain::top ();
       for (auto cst : csts)
         inv += cst;
       return inv;
     }

     static lin_cst_sys_t top () {
       return mk_true ();
     }

     static lin_cst_sys_t bot () {
       return mk_false ();
     }
   };

}
#endif 
