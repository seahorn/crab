#ifndef ARRAY_SEGMENTATION_HPP
#define ARRAY_SEGMENTATION_HPP

/* 
   Quick analysis to collect the set of variables and constants that might
   define array segment boundaries.
*/

#include <crab/common/types.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>

#include <crab/cfg/Cfg.hpp> 
#include <crab/iterators/killgen_fixpoint_iterator.hpp> 

#include <boost/unordered_map.hpp>
#include <boost/range/iterator_range.hpp>

namespace crab {

   namespace analyzer{

     template<class CFG, class V>
     class array_segment_analysis: public crab::iterators::killgen_analysis<CFG, V> {
       
       typedef crab::iterators::killgen_analysis<CFG, V> killgen_analysis_t;
       
      public:
       
       typedef typename CFG::basic_block_label_t basic_block_label_t;
       typedef typename killgen_analysis_t::killgen_domain_t array_segment_domain_t;       

      private:
       
       class array_segment_visitor: public cfg::StatementVisitor<V> {
         typedef typename cfg::StatementVisitor<V>::z_bin_op_t z_bin_op_t;
         typedef typename cfg::StatementVisitor<V>::z_assign_t z_assign_t;
         typedef typename cfg::StatementVisitor<V>::z_assume_t z_assume_t;
         typedef typename cfg::StatementVisitor<V>::havoc_t havoc_t;
         typedef typename cfg::StatementVisitor<V>::unreach_t unreach_t;
         typedef typename cfg::StatementVisitor<V>::z_select_t z_select_t;
         typedef typename cfg::StatementVisitor<V>::callsite_t callsite_t;
         typedef typename cfg::StatementVisitor<V>::z_assert_t z_assert_t;
         typedef typename cfg::StatementVisitor<V>::z_arr_load_t z_arr_load_t;
         typedef typename cfg::StatementVisitor<V>::z_arr_store_t z_arr_store_t;

         // assume all statements have the same type expression_t;
         typedef typename z_bin_op_t::linear_expression_t linear_expression_t;
         typedef typename z_assume_t::linear_constraint_t linear_constraint_t;
         
         array_segment_domain_t get_variables (linear_expression_t e) {
           array_segment_domain_t res;
           for (auto v: e.variables()) res += v.name();
           return res;
         }

         array_segment_domain_t get_variables (linear_constraint_t c) {
           array_segment_domain_t res;
           for (auto v: c.variables()) res += v.name();
           return res;
         }

        public:

         array_segment_domain_t _indexes;

         array_segment_visitor (array_segment_domain_t inv)
             : _indexes(inv) { }
         
         void visit(z_bin_op_t &s){ 
           if (!(_indexes & s.lhs().name()).is_bottom()) {
             _indexes += get_variables(s.left());
             _indexes += get_variables(s.right());
           }
         }  
         
         void visit(z_assign_t &s) {
           if (!(_indexes & s.lhs().name()).is_bottom()) {
             _indexes += get_variables(s.rhs());
           }
         }
         
         void visit(z_assume_t &s) { 
           auto vars = get_variables(s.constraint());
           if (!(_indexes & vars).is_bottom())
             _indexes += vars;
         }
         
         void visit(z_arr_load_t &s) {
           _indexes += get_variables(s.index());
         }
         
         void visit(z_arr_store_t &s) {
           _indexes += get_variables(s.index());
         }

         void visit(unreach_t&) { }
         // FIXME: implement these 
         void visit(z_select_t&) { }
         void visit(callsite_t&) { }
         void visit(havoc_t&) { }
         void visit(z_assert_t&) { }
       };
       
      public:
       
       array_segment_analysis (CFG cfg): killgen_analysis_t(cfg) { }

       virtual bool is_forward() override { return false;}
       
       virtual std::string name () override { return "array_segment";}

       virtual void init_fixpoint() override { }

       virtual array_segment_domain_t entry() override {
         return array_segment_domain_t::bottom();
       }
       
       virtual array_segment_domain_t merge(array_segment_domain_t d1,  
                                            array_segment_domain_t d2)  override {
         return d1 | d2; 
       }
       
       virtual array_segment_domain_t analyze (basic_block_label_t bb_id, 
                                               array_segment_domain_t in) override {
         auto &bb = this->_cfg.get_node(bb_id);
         array_segment_visitor vis(in);
         for (auto &s: bb) { s.accept (&vis); }
         return vis._indexes;
       }
     };
   
     template<class CFG, class V>
     class ArraySegmentation: public boost::noncopyable, 
                              public crab::iterators::killgen_fixpoint_iterator
         <CFG, array_segment_analysis<CFG,V> >{
      public:

       typedef array_segment_analysis<CFG,V> array_segment_analysis_t;
       typedef typename CFG::basic_block_label_t basic_block_label_t;
       typedef typename CFG::statement_t statement_t;
       typedef typename CFG::varname_t varname_t;       
       typedef typename array_segment_analysis_t::array_segment_domain_t array_segment_domain_t;

      private:

       typedef crab::iterators::killgen_fixpoint_iterator<CFG,array_segment_analysis_t> 
                killgen_fixpoint_iterator_t;
       typedef boost::unordered_map<basic_block_label_t, array_segment_domain_t> segment_map_t;
       
       segment_map_t _segment_map;
       
      public:
       
       ArraySegmentation (CFG cfg)
           : killgen_fixpoint_iterator_t(cfg) { }
       
       void exec() { 
         this->run();
         
         CRAB_LOG("array-segment",
                  crab::outs() << "Array segment variables alive at the block entries\n";);
         
         for (auto p : boost::make_iterator_range(this->out_begin(), this->out_end())) {
           CRAB_LOG("array-segment",
                    crab::outs() << p.first << ":" << p.second << "\n";);
           // XXX: we want the variables from the IN set but because
           // we are working on the reversed view of the cfg we get
           // them instead from the OUT set.
           _segment_map.insert(std::make_pair(p.first, p.second)); 
         } 
         this->release_memory();
       }
       
       array_segment_domain_t get_variables (basic_block_label_t bb) const {
         auto it = _segment_map.find(bb);
         if (it == _segment_map.end()) 
           return array_segment_domain_t();
         else 
           return it->second; 
       }
       
       void write (ostream &o) const { }
     }; 

     // Visitor for finding constants that might appear as array
     // segment boundaries.
     template<class ArraySegmentDom>
     class array_constant_segment_visitor: 
         public cfg::StatementVisitor<typename ArraySegmentDom::element_t> {
       typedef typename ArraySegmentDom::element_t E;
       typedef typename cfg::StatementVisitor<E>::z_bin_op_t z_bin_op_t;
       typedef typename cfg::StatementVisitor<E>::z_assign_t z_assign_t;
       typedef typename cfg::StatementVisitor<E>::z_assume_t z_assume_t;
       typedef typename cfg::StatementVisitor<E>::havoc_t havoc_t;
       typedef typename cfg::StatementVisitor<E>::unreach_t unreach_t;
       typedef typename cfg::StatementVisitor<E>::z_select_t z_select_t;
       typedef typename cfg::StatementVisitor<E>::callsite_t callsite_t;
       typedef typename cfg::StatementVisitor<E>::z_assert_t z_assert_t;
       typedef typename cfg::StatementVisitor<E>::z_arr_load_t z_arr_load_t;
       typedef typename cfg::StatementVisitor<E>::z_arr_store_t z_arr_store_t;
       
       typedef typename z_bin_op_t::linear_expression_t linear_expression_t;
       typedef typename linear_expression_t::number_t number_t;

      public:

       typedef std::vector<number_t> constant_set_t;

      private:

       ArraySegmentDom _dom;
       constant_set_t _csts;
       
      public:

       array_constant_segment_visitor (ArraySegmentDom dom): _dom(dom) { }

       constant_set_t get_constants () { return _csts;}

       /// XXX: we focus for now only on assignments
       void visit(z_assign_t &s) {
         if (!(_dom & s.lhs().name()).is_bottom()) {
           if (s.rhs().is_constant() && s.rhs().constant() >= 0 &&
               std::find(_csts.begin(), _csts.end(), s.rhs().constant()) == _csts.end())
             _csts.push_back(s.rhs().constant());
         }
       }

       void visit(z_bin_op_t &s){}         
       void visit(z_assume_t &s) {}
       void visit(z_arr_load_t &s) {}
       void visit(z_arr_store_t &s) {}
       void visit(z_select_t&) {}
       void visit(callsite_t&) {}
       void visit(havoc_t&) {}
       void visit(z_assert_t&) {}
       void visit(unreach_t&) {}
     };

  } // end namespace analyzer
} // end namespace crab
#endif 
