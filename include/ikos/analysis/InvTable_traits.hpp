#ifndef INVARIANT_TABLE_ENTRIES_TRAITS_HPP
#define INVARIANT_TABLE_ENTRIES_TRAITS_HPP

/* 
   This class manipulates the values stored in the invariant
   tables. Invariant tables are how ikos-core communicates with other
   clients. This class allows clients to determine the type of the
   table values.
 */

namespace analyzer {

   template<typename AbsDomain, typename InvTableTy>
   class inv_tbl_traits {
    public:
     //! Convert an abstract domain to the format stored in the
     //  invariant tables.
     // FIXME: should be `const AbsDomain&` but we need to fix first
     // constness of the abstract domain methods.
     static InvTableTy marshall (AbsDomain& abs);
     //! Convert from the format stored in the invariant table to an
     //! abstract domain.
     // FIXME: should be `const InvTableTy&` but we need to fix first
     // constness of linear_constraint
     static AbsDomain unmarshall (InvTableTy& inv);
     //! Return true 
     static InvTableTy top ();
     //! Return false
     static InvTableTy bot ();
   };


   // Default implementation: conjunctive linear constraints
   template<typename AbsDomain>
   class inv_tbl_traits <AbsDomain,
                         typename AbsDomain::linear_constraint_system_t>  {

     typedef typename AbsDomain::linear_expression_t lin_exp_t;

    public:

     typedef typename AbsDomain::linear_constraint_t lin_cst_t;
     typedef typename AbsDomain::linear_constraint_system_t lin_cst_sys_t;

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

     static AbsDomain unmarshall (lin_cst_sys_t& csts) {
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
