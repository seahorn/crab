#pragma once

#include <crab/domains/region_domain.hpp>

namespace crab {
namespace domains {

////////////////////////////////////////////////////////////////////////
// Smashing-based abstract domain for regions and references
//
// It does not generate ghost variables so it should be faster. Ghost
// variables are expensive to maintain because
// join/widening/meet/narrowing need to rename them.
////////////////////////////////////////////////////////////////////////
namespace region_without_ghost_domain_impl {
class Params {
public:
  // Keep track of allocation sites 
  enum { allocation_sites = 1 };  
  // Reason about deallocation
  enum { deallocation = 1 };
  // Improve precision but it might not be sound
  enum { refine_uninitialized_regions = 1 };
  // Enable tag analysis
  enum { tag_analysis = 1};
};
} // end namespace region_without_ghost_domain_impl
  
template <typename BaseAbsDom, typename Params = region_without_ghost_domain_impl::Params>
class region_without_ghost_domain final :
    public abstract_domain_api<region_without_ghost_domain<BaseAbsDom, Params>> {
  using region_domain_t = region_without_ghost_domain<BaseAbsDom, Params>;
  using abstract_domain_t = abstract_domain_api<region_domain_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;  
  using number_t = typename BaseAbsDom::number_t;
  using varname_t = typename BaseAbsDom::varname_t;

private:
  using base_abstract_domain_t = BaseAbsDom;
  // Environment domains: map regions to finite domain
  using rgn_counting_env_t = ikos::separate_domain<variable_t, small_range>;
  using rgn_bool_env_t = ikos::separate_domain<variable_t, boolean_value>;
  // Union-find where equivalence classes are attached to boolean values
  using rgn_dealloc_t = union_find_domain<variable_t, boolean_value>;
  // Map each reference or region variable to a set of allocation sites
  using alloc_site_env_t = separate_discrete_domain<variable_t, allocation_site>;
  using allocation_sites = typename alloc_site_env_t::value_type;
  // Tags
  // JN: we don't use crab::tag because we want to have the
  // flexibility of creating tags without a tag manager.
  class tag_t : public indexable {
    ikos::index_t m_id;
  public:
    tag_t(number_t n): m_id(0) {
      if (n < 0)
      { CRAB_ERROR("Cannot use negative numbers for tags");}
      if (!n.fits_int64())
      { CRAB_ERROR("Too large value for a tag"); }
      m_id = (int64_t) n;
    }
    bool operator<(const tag_t &as) const { return m_id < as.m_id; }
    bool operator==(const tag_t &as) const { return m_id == as.m_id; }
    virtual ikos::index_t index() const override { return m_id; }
    void write(crab_os &o) const override { o << "TAG_" << m_id; }
    friend crab_os &operator<<(crab_os &o, const tag_t &as) {
      as.write(o);
      return o;
    }
  };
  // Map variables to sets of tags
  using tag_env_t = separate_discrete_domain<variable_t, tag_t>;
  using tag_set = typename tag_env_t::value_type;
  
  bool m_is_bottom; // special symbol for bottom

  // The base abstract domain: all the heavy lifting is done here.
  // m_base_dom can have variables of any type, included region, array
  // and reference types.
  base_abstract_domain_t m_base_dom;
  
  // Abstract domain to count how many addresses are in a region.
  // This allows us to decide when strong update is sound: only if
  // one address per region (i.e., singleton).
  rgn_counting_env_t m_rgn_counting_dom;
  // Whether some data might have been written to any address within
  // the region.
  rgn_bool_env_t m_rgn_init_dom;
  // Map each reference to its set of possible allocation sites.
  // 
  // It also maps each region variable to the union of all
  // possible allocation sites from all possible references stored in
  // that region. Note that we don't keep track of which allocation
  // sites allocated the memory where the region lives. This is a
  // different question. We are only interested in references. Regions
  // need to be tracked because references are stored in regions.
  alloc_site_env_t m_alloc_site_dom;
  // Keep track of whether some memory within a region has been
  // deallocated.
  // 
  // To reason about deallocation, we need to know which regions might
  // belong to the same allocated memory object.  Each call to
  // ref_make models a new allocation returning a reference to the
  // base address of the allocated block. Then, ref_gep is used to
  // perform pointer arithmetic within an allocated block. Very
  // importantly, ref_gep can switch between regions although we
  // assume that those regions always belong to the same allocated
  // block.
  //
  // We partition region variables into equivalence classes attached
  // to a boolean value. Two region variables are in the same
  // equivalence class if they might belong to the same allocated
  // block. The boolean flag indicates whether the allocated block
  // might have been deallocated.
  //
  // The partitioning is done as follows:
  //
  //   region_init(rgn): create a singleton equivalence class with rgn.
  //   ref_gep(reg1, rgn1, ref2, rgn2, o): join together the
  //                                       equivalence classes of rgn1
  //                                       and rgn2.
  rgn_dealloc_t m_rgn_dealloc_dom;
  // Tag analysis: map each (any type) variable to a set of tags
  tag_env_t m_tag_env;
  
  region_without_ghost_domain(base_abstract_domain_t &&base_dom,
			      rgn_counting_env_t &&rgn_counting_dom,
			      rgn_bool_env_t &&rgn_init_dom, 
			      alloc_site_env_t &&alloc_site_dom,
			      rgn_dealloc_t &&rgn_dealloc_dom,
			      tag_env_t &&tag_env)
      : m_is_bottom(base_dom.is_bottom()), 
	m_base_dom(std::move(base_dom)),	
        m_rgn_counting_dom(std::move(rgn_counting_dom)),
        m_rgn_init_dom(std::move(rgn_init_dom)),
	m_alloc_site_dom(std::move(alloc_site_dom)),
        m_rgn_dealloc_dom(std::move(rgn_dealloc_dom)),
	m_tag_env(std::move(tag_env)) {}

  using base_dom_binop_t = std::function<base_abstract_domain_t(
      base_abstract_domain_t, base_abstract_domain_t)>;

  bool can_propagate_initialized_regions(
      const rgn_bool_env_t &left_init_dom,
      const rgn_bool_env_t &right_init_dom,
      variable_vector_t &regions) const {
    bool propagate = false;
    for (auto kv : left_init_dom) {
      if (kv.second.is_false()) {
	// note that this is maybe-initialized
	if (!right_init_dom[kv.first].is_false()) {
	  regions.push_back(kv.first);
          propagate = true;
	}
      }
    }
    return propagate;
  }

  // Special step: if one region is uninitialized on left operand but
  // initialized on the right then we extract the information from the
  // right and add it to the left. This can happen if the first store
  // to a region happens inside a loop. At the entry of the loop
  // nothing is known about the region so the join there would lose
  // all the information about the region.
  //
  // This step is not perfect because it keeps only non-relational
  // information about regions but no relational information between
  // regions and other variables.
  void refine_regions(const variable_vector_t &regions,
                      /* left operand */
                      base_abstract_domain_t &left_dom,
                      /* right operand */
                      const base_abstract_domain_t &right_dom) const {
    
    if (regions.empty()) {
      return;
    } 

    base_abstract_domain_t dom(right_dom);
    dom.project(regions);
    left_dom = left_dom & dom;
  }

  // Perform *this = join(*this, right)
  void do_join(const region_domain_t &right) {
    rgn_counting_env_t out_rgn_counting_dom(m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(m_rgn_init_dom | right.m_rgn_init_dom);
    alloc_site_env_t out_alloc_site_dom(m_alloc_site_dom |
					right.m_alloc_site_dom);
    rgn_dealloc_t out_rgn_dealloc_dom(m_rgn_dealloc_dom |
                                      right.m_rgn_dealloc_dom);
    tag_env_t out_tag_env(m_tag_env | right.m_tag_env);

    m_base_dom |= right.m_base_dom;    
    m_is_bottom = m_base_dom.is_bottom();
    std::swap(m_rgn_counting_dom, out_rgn_counting_dom);
    std::swap(m_rgn_init_dom, out_rgn_init_dom);
    std::swap(m_alloc_site_dom, out_alloc_site_dom);    
    std::swap(m_rgn_dealloc_dom, out_rgn_dealloc_dom);
    std::swap(m_tag_env, out_tag_env);

    CRAB_LOG("region", crab::outs() << *this << "\n");
  }

  region_domain_t do_join_or_widening(const region_domain_t &left,
                                      const region_domain_t &right,
                                      base_dom_binop_t base_dom_op) const {

    // rgn_counting_dom does not require common renaming (i.e., no
    // ghost variables).  The domain is finite
    rgn_counting_env_t out_rgn_counting_dom(left.m_rgn_counting_dom |
                                            right.m_rgn_counting_dom);
    // rgn_init_dom does not require common renaming (i.e., no ghost
    // variables).  The domain is finite
    rgn_bool_env_t out_rgn_init_dom(left.m_rgn_init_dom | right.m_rgn_init_dom);

    // alloc_site_dom does not require common renaming (i.e., no ghost
    // variables).  The domain is finite
    alloc_site_env_t out_alloc_site_dom(left.m_alloc_site_dom | right.m_alloc_site_dom);
    
    // rgn_dealloc_dom does not require common renaming (i.e., no
    // ghost variables).  The domain is finite
    rgn_dealloc_t out_rgn_dealloc_dom(left.m_rgn_dealloc_dom |
                                      right.m_rgn_dealloc_dom);
    // tag_env does not require common renaming (i.e., no ghost
    // variables).  The domain is finite
    tag_env_t out_tag_env(left.m_tag_env | right.m_tag_env);
    
    // join or widening in the base domain
    base_abstract_domain_t out_base_dom(base_dom_op(left.m_base_dom, right.m_base_dom));

    region_domain_t res(
        std::move(out_base_dom),			
	std::move(out_rgn_counting_dom),
        std::move(out_rgn_init_dom), 
	std::move(out_alloc_site_dom), std::move(out_rgn_dealloc_dom),
	std::move(out_tag_env));
    return res;
  }

  region_domain_t do_meet_or_narrowing(const region_domain_t &left,
                                       const region_domain_t &right,
                                       base_dom_binop_t base_dom_op) const {

    // these domains are finite
    rgn_counting_env_t out_rgn_counting_dom(left.m_rgn_counting_dom &
                                            right.m_rgn_counting_dom);
    rgn_bool_env_t out_rgn_init_dom(left.m_rgn_init_dom & right.m_rgn_init_dom);
    alloc_site_env_t out_alloc_site_dom(left.m_alloc_site_dom &
					right.m_alloc_site_dom);    
    rgn_dealloc_t out_rgn_dealloc_dom(left.m_rgn_dealloc_dom &
                                      right.m_rgn_dealloc_dom);
    tag_env_t out_tag_env(left.m_tag_env & right.m_tag_env);
    
    // This shouldn't happen but just in case ...
    if (out_rgn_counting_dom.is_bottom() || out_rgn_init_dom.is_bottom() ||
        out_rgn_dealloc_dom.is_bottom() || out_alloc_site_dom.is_bottom()) {
      return make_bottom();
    }

    base_abstract_domain_t out_base_dom(base_dom_op(left.m_base_dom,
						    right.m_base_dom));

    region_domain_t res(
	std::move(out_base_dom),			
        std::move(out_rgn_counting_dom),
        std::move(out_rgn_init_dom), 
	std::move(out_alloc_site_dom), std::move(out_rgn_dealloc_dom),
	std::move(out_tag_env));

    return res;
  }
  template<class RangeVars>
  void merge_tags(const variable_t &x, RangeVars vars) {
    tag_set tags = tag_set::bottom();
    for (auto const &v: vars) {
      tags = tags | m_tag_env[v];
    }
    m_tag_env.set(x, tags);
  }

  linear_constraint_t
  convert_ref_cst_to_linear_cst(const reference_constraint_t &ref_cst) {
    if (ref_cst.is_tautology()) {
      return linear_constraint_t::get_true();
    } else if (ref_cst.is_contradiction()) {
      return linear_constraint_t::get_false();
    } else {
      if (ref_cst.is_unary()) {
        assert(ref_cst.lhs().get_type().is_reference());
        variable_t x = ref_cst.lhs();
        if (ref_cst.is_equality()) {
          return linear_constraint_t(x == number_t(0));
        } else if (ref_cst.is_disequality()) {
          return linear_constraint_t(x != number_t(0));
        } else if (ref_cst.is_less_or_equal_than()) {
          return linear_constraint_t(x <= number_t(0));
        } else if (ref_cst.is_less_than()) {
          return linear_constraint_t(x < number_t(0));
        } else if (ref_cst.is_greater_or_equal_than()) {
          return linear_constraint_t(x >= number_t(0));
        } else if (ref_cst.is_greater_than()) {
          return linear_constraint_t(x > number_t(0));
        }
      } else {
        assert(ref_cst.lhs().get_type().is_reference());
        assert(ref_cst.rhs().get_type().is_reference());
        variable_t x = ref_cst.lhs();
        variable_t y = ref_cst.rhs();
        number_t offset = ref_cst.offset();
        if (ref_cst.is_equality()) {
          return linear_constraint_t(x == y + offset);
        } else if (ref_cst.is_disequality()) {
          return linear_constraint_t(x != y + offset);
        } else if (ref_cst.is_less_or_equal_than()) {
          return linear_constraint_t(x <= y + offset);
        } else if (ref_cst.is_less_than()) {
          return linear_constraint_t(x < y + offset);
        } else if (ref_cst.is_greater_or_equal_than()) {
          return linear_constraint_t(x >= y + offset);
        } else if (ref_cst.is_greater_than()) {
          return linear_constraint_t(x > y + offset);
        }
      }
    }
    CRAB_ERROR("unexpected reference constraint");
  }

  // m_base_domain contains variables of type region, array, and
  // reference since we didn't want to add ghost variables.  We don't
  // expose a linear constraint if any of its variables is of type
  // region or array.
  linear_constraint_system_t
  filter_nonscalar_vars(linear_constraint_system_t &&csts) const {
    linear_constraint_system_t res;
    for (auto const &cst : csts) {
      if (std::all_of(
              cst.expression().variables_begin(),
              cst.expression().variables_end(), [](const variable_t &v) {
						  return v.get_type().is_scalar();
              })) {
        res += cst;
      }
    }
    return res;
  }
  
public:
  region_domain_t make_top() const override { return region_domain_t(true); }

  region_domain_t make_bottom() const override {
    return region_domain_t(false);
  }

  void set_to_top() override {
    region_domain_t abs(true);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    region_domain_t abs(false);
    std::swap(*this, abs);
  }

  region_without_ghost_domain(bool is_top = true) : m_is_bottom(!is_top) {}

  region_without_ghost_domain(const region_domain_t &o)
      : m_is_bottom(o.m_is_bottom), 
	m_base_dom(o.m_base_dom),	
        m_rgn_counting_dom(o.m_rgn_counting_dom),
        m_rgn_init_dom(o.m_rgn_init_dom),
	m_alloc_site_dom(o.m_alloc_site_dom),
        m_rgn_dealloc_dom(o.m_rgn_dealloc_dom),
	m_tag_env(o.m_tag_env) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }
  region_without_ghost_domain(region_domain_t &&o)
      : m_is_bottom(o.m_is_bottom), 
	m_base_dom(std::move(o.m_base_dom)),	
        m_rgn_counting_dom(std::move(o.m_rgn_counting_dom)),
        m_rgn_init_dom(std::move(o.m_rgn_init_dom)),
	m_alloc_site_dom(std::move(o.m_alloc_site_dom)),
        m_rgn_dealloc_dom(std::move(o.m_rgn_dealloc_dom)),
	m_tag_env(std::move(o.m_tag_env)) {}

  region_domain_t &operator=(const region_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_base_dom = o.m_base_dom;      
      m_rgn_counting_dom = o.m_rgn_counting_dom;
      m_rgn_init_dom = o.m_rgn_init_dom;
      m_alloc_site_dom = o.m_alloc_site_dom;
      m_rgn_dealloc_dom = o.m_rgn_dealloc_dom;
      m_tag_env = o.m_tag_env;
    }
    return *this;
  }

  region_domain_t &operator=(region_domain_t &&o) {
    if (this != &o) {
      m_is_bottom = std::move(o.m_is_bottom);
      m_base_dom = std::move(o.m_base_dom);      
      m_rgn_counting_dom = std::move(o.m_rgn_counting_dom);
      m_rgn_init_dom = std::move(o.m_rgn_init_dom);
      m_alloc_site_dom = std::move(o.m_alloc_site_dom);
      m_rgn_dealloc_dom = std::move(o.m_rgn_dealloc_dom);
      m_tag_env = std::move(o.m_tag_env);
    }
    return *this;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override {
    bool res = (!is_bottom() && m_base_dom.is_top() && m_rgn_counting_dom.is_top());
    if (Params::allocation_sites) {
      res = res && m_alloc_site_dom.is_top();
    }
    if (Params::deallocation) {
      res = res && m_rgn_dealloc_dom.is_top();
    }
    if (Params::tag_analysis) {
      res = res && m_tag_env.is_top();
    }
    return res;
  }

  bool operator<=(const region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    CRAB_LOG("region-leq", crab::outs() << "Inclusion test:\n\t" << *this
                                        << "\n\t" << o << "\n";);
    if (is_bottom() || o.is_top()) {
      CRAB_LOG("region-leq", crab::outs() << "Result1=1\n";);
      return true;
    } else if (is_top() || o.is_bottom()) {
      CRAB_LOG("region-leq", crab::outs() << "Result2=0\n";);
      return false;
    }

    if (!(m_rgn_counting_dom <= o.m_rgn_counting_dom)) {
      CRAB_LOG("region-leq", crab::outs() << "Result3=0\n";);
      return false;
    }

    if (!(m_rgn_init_dom <= o.m_rgn_init_dom)) {
      CRAB_LOG("region-leq", crab::outs() << "Result4=0\n";);
      return false;
    }

    if (Params::allocation_sites) {
      if (!(m_alloc_site_dom <= o.m_alloc_site_dom)) {
        CRAB_LOG("region-leq", crab::outs() << "Result5=0\n";);
        return false;
      }
    }
    
    if (Params::deallocation) {
      if (!(m_rgn_dealloc_dom <= o.m_rgn_dealloc_dom)) {
        CRAB_LOG("region-leq", crab::outs() << "Result6=0\n";);
        return false;
      }
    }

    if (Params::tag_analysis) {
      if (!(m_tag_env <= o.m_tag_env)) {
        CRAB_LOG("region-leq", crab::outs() << "Result7=0\n";);
        return false;
      }
    }

    bool res = m_base_dom <= o.m_base_dom;
    CRAB_LOG("region-leq", crab::outs() << "Result8=" << res << "\n";);
    return res;
  }

  void operator|=(const region_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    // Trivial cases first
    if (is_bottom()) {
      *this = o;
      return;
    } else if (o.is_bottom()) {
      return;
    } else if (is_top() | o.is_top()) {
      set_to_top();
      return;
    }

    CRAB_LOG("region", crab::outs() << "Join " << *this << " and " << o << "=");

    variable_vector_t left_regions, right_regions;
    bool refine_left = false;
    bool refine_right = false;
    if (Params::refine_uninitialized_regions) {
      refine_left = can_propagate_initialized_regions(
          m_rgn_init_dom, o.m_rgn_init_dom, left_regions);
      refine_right = can_propagate_initialized_regions(
          o.m_rgn_init_dom, m_rgn_init_dom, right_regions);
    }

    // The code is complicated to achieve a zero-cost abstraction. If
    // we cannot improve invariants on the right operand we avoid
    // making a copy of it.
    if (refine_left && !refine_right) {
      refine_regions(left_regions, m_base_dom, o.m_base_dom);
      do_join(o);
    } else if (!refine_left && refine_right) {
      region_domain_t right(o);
      refine_regions(right_regions, right.m_base_dom, m_base_dom);
      do_join(right);
    } else if (refine_left && refine_right) {
      region_domain_t right(o);
      refine_regions(left_regions, m_base_dom, right.m_base_dom);
      refine_regions(right_regions, right.m_base_dom, m_base_dom);
      do_join(right);
    } else {
      do_join(o);
    }
  }

  region_domain_t operator|(const region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    // Trivial cases first
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    } else if (is_top() | o.is_top()) {
      region_domain_t abs;
      abs.set_to_top();
      return abs;
    }

    CRAB_LOG("region", crab::outs() << "Join " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) { return v1 | v2; };

    variable_vector_t left_regions, right_regions;

    bool refine_left = false;
    bool refine_right = false;
    if (Params::refine_uninitialized_regions) {
      refine_left = can_propagate_initialized_regions(
          m_rgn_init_dom, o.m_rgn_init_dom, left_regions);
      refine_right = can_propagate_initialized_regions(
          o.m_rgn_init_dom, m_rgn_init_dom, right_regions);
    }

    // The code is complicated to achieve a zero-cost abstraction. If
    // we cannot improve invariants then we try to avoid making copies
    // of left and/or right operands.

    if (refine_left && !refine_right) {
      // Refine left by propagating information from right's regions
      region_domain_t left(*this);
      refine_regions(left_regions,left.m_base_dom, o.m_base_dom);
      region_domain_t res(std::move(do_join_or_widening(left, o, base_dom_op)));
      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;

    } else if (!refine_left && refine_right) {
      // Refine right by propagating information from left's regions
      region_domain_t right(o);
      refine_regions(right_regions, right.m_base_dom, m_base_dom);
      region_domain_t res(
          std::move(do_join_or_widening(*this, right, base_dom_op)));
      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;

    } else if (refine_left && refine_right) {
      // Refine both left and right
      region_domain_t left(*this);
      region_domain_t right(o);
      refine_regions(left_regions, left.m_base_dom, right.m_base_dom);
      refine_regions(right_regions, right.m_base_dom, left.m_base_dom);
      region_domain_t res(
          std::move(do_join_or_widening(left, right, base_dom_op)));
      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;

    } else {
      region_domain_t res(
          std::move(do_join_or_widening(*this, o, base_dom_op)));

      CRAB_LOG("region", crab::outs() << res << "\n");
      return res;
    }

  }

  region_domain_t operator&(const region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("region", crab::outs() << "Meet " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) {
      return (v1 & v2);
    };
    region_domain_t res(std::move(do_meet_or_narrowing(*this, o, base_dom_op)));

    CRAB_LOG("region", crab::outs() << res << "\n");
    return res;
  }

  region_domain_t operator||(const region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    // Trivial cases first: we don't cover cases where one operand is
    // top because is_top() calls the base domain which we don't know
    // whether it will perform some normalization or not.
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    }

    CRAB_LOG("region", crab::outs()
                           << "Widening " << *this << " and " << o << "=");

    auto base_dom_op = [](const base_abstract_domain_t &v1,
                          const base_abstract_domain_t &v2) {
      return v1 || v2;
    };
    region_domain_t res(std::move(do_join_or_widening(*this, o, base_dom_op)));

    CRAB_LOG("region", crab::outs() << res << "\n");
    return res;
  }

  region_domain_t widening_thresholds(
      const region_domain_t &o,
      const iterators::thresholds<number_t> &thresholds) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    // Trivial cases first: we don't cover cases where one operand is
    // top because is_top() calls the base domain which we don't know
    // whether it will perform some normalization or not.
    if (is_bottom()) {
      return o;
    } else if (o.is_bottom()) {
      return *this;
    }

    CRAB_LOG("region", crab::outs()
                           << "Widening " << *this << " and " << o << "=");

    auto base_dom_op = [&thresholds](const base_abstract_domain_t &v1,
                                     const base_abstract_domain_t &v2) {
      return v1.widening_thresholds(v2, thresholds);
    };

    region_domain_t res(std::move(do_join_or_widening(*this, o, base_dom_op)));

    CRAB_LOG("region", crab::outs() << res << "\n");
    return res;
  }

  region_domain_t operator&&(const region_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (o.is_bottom() || is_top()) {
      return o;
    }

    CRAB_LOG("region", crab::outs() << "Meet " << *this << " and " << o << "=");

    auto base_dom_op = [](base_abstract_domain_t v1,
                          base_abstract_domain_t v2) { return (v1 && v2); };
    region_domain_t res(std::move(do_meet_or_narrowing(*this, o, base_dom_op)));

    CRAB_LOG("region", crab::outs() << res << "\n");
    return res;
  }

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (!is_bottom()) {
      if (v.get_type().is_region()) {
        m_rgn_counting_dom -= v;
        m_rgn_init_dom -= v;	
        if (Params::deallocation) {
          m_rgn_dealloc_dom.forget(v);
        }
	
      }
      if (Params::allocation_sites) {
	if (v.get_type().is_reference() || v.get_type().is_region()) {
	  m_alloc_site_dom -= v;
	}
      }
      if (Params::tag_analysis) {
	m_tag_env -= v;
      }

      m_base_dom -= v;
    }
  }

  // Initialize region rgn.
  void region_init(const variable_t &rgn) override {
    crab::CrabStats::count(domain_name() + ".count.region_init");
    crab::ScopedCrabStats __st__(domain_name() + ".region_init");

    if (is_bottom()) {
      return;
    }

    small_range count_num = m_rgn_counting_dom[rgn];
    if (count_num <= small_range::oneOrMore()) {
      CRAB_ERROR("region_domain::init_region: ", rgn,
                 " cannot be initialized twice");
    }

    // Set to zero the number of references
    m_rgn_counting_dom.set(rgn, small_range::zero());

    // No stores in the region
    m_rgn_init_dom.set(rgn, boolean_value::get_false());

    if (Params::deallocation) {
      // No deallocated objects in the region
      m_rgn_dealloc_dom.set(rgn, boolean_value::get_false());
    }
    if (Params::allocation_sites) {
      // Region does not contain any allocation site
      m_alloc_site_dom.set(rgn, alloc_site_env_t::value_type::bottom());
    }
    if (Params::tag_analysis) {
      // Region does not contain any tag
      m_tag_env.set(rgn, tag_env_t::value_type::bottom());
    }
    CRAB_LOG("region", crab::outs() << "After region_init(" << rgn
                                    << ")=" << *this << "\n";);
  }

  void region_copy(const variable_t &lhs_rgn,
                   const variable_t &rhs_rgn) override {
    crab::CrabStats::count(domain_name() + ".count.region_copy");
    crab::ScopedCrabStats __st__(domain_name() + ".region_copy");

    if (is_bottom()) {
      return;
    }

    if (!lhs_rgn.get_type().is_region()) {
      CRAB_ERROR(domain_name() + "::region_copy ", lhs_rgn,
                 " must have a region type");
    }
    if (!rhs_rgn.get_type().is_region()) {
      CRAB_ERROR(domain_name() + "::region_copy ", rhs_rgn,
                 " must have a region type");
    }

    m_rgn_counting_dom.set(lhs_rgn, m_rgn_counting_dom[rhs_rgn]);
    m_rgn_init_dom.set(lhs_rgn, m_rgn_init_dom[rhs_rgn]);
    if (Params::allocation_sites) {
      m_alloc_site_dom.set(lhs_rgn, m_alloc_site_dom[rhs_rgn]);
    }
    if (Params::deallocation) {
      m_rgn_dealloc_dom.add(rhs_rgn, lhs_rgn);      
    }
    if (Params::tag_analysis) {
      m_tag_env.set(lhs_rgn, m_tag_env[rhs_rgn]);
    }
    
    if (lhs_rgn.get_type().is_unknown_region() ||
        rhs_rgn.get_type().is_unknown_region()) {
      return;
    }

    if (lhs_rgn.get_type() != rhs_rgn.get_type()) {
      CRAB_ERROR(domain_name() + "::region_copy ", lhs_rgn, ":=", rhs_rgn,
                 " with different types");
    }

    m_base_dom.forget({lhs_rgn});
    m_base_dom.expand(rhs_rgn, lhs_rgn);
  }


  void region_cast(const variable_t &src_rgn,
		   const variable_t &dst_rgn) override {
                   
    crab::CrabStats::count(domain_name() + ".count.region_cast");
    crab::ScopedCrabStats __st__(domain_name() + ".region_cast");

    if (is_bottom()) {
      return;
    }

    if (!dst_rgn.get_type().is_region()) {
      CRAB_ERROR(domain_name() + "::region_cast ", dst_rgn,
                 " must have a region type");
    }
    if (!src_rgn.get_type().is_region()) {
      CRAB_ERROR(domain_name() + "::region_cast ", src_rgn,
                 " must have a region type");
    }

    m_rgn_counting_dom.set(dst_rgn, m_rgn_counting_dom[src_rgn]);
    m_rgn_init_dom.set(dst_rgn, m_rgn_init_dom[src_rgn]);
    if (Params::allocation_sites) {
      m_alloc_site_dom.set(dst_rgn, m_alloc_site_dom[src_rgn]);
    }
    if (Params::deallocation) {
      m_rgn_dealloc_dom.add(src_rgn, dst_rgn);
    }
    if (Params::tag_analysis) {
      m_tag_env.set(dst_rgn, m_tag_env[src_rgn]);
    }
    
    if (dst_rgn.get_type() == src_rgn.get_type()) {
      // nothing to do
      return;
    }

    CRAB_LOG("region",
	     CRAB_WARN("TODO region_cast ", src_rgn, ":", src_rgn.get_type(), " to ",
		       dst_rgn, ":", dst_rgn.get_type(), " in base domain"););    
  }
  
  // Create a new reference ref to region rgn.
  void ref_make(const variable_t &ref, const variable_t &rgn,
		const variable_or_constant_t &size /*unused*/,
		const allocation_site &as) override {
    crab::CrabStats::count(domain_name() + ".count.ref_make");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_make");

    if (is_bottom()) {
      return;
    }

    if (!ref.get_type().is_reference()) {
      CRAB_LOG("region", CRAB_WARN("region_domain::ref_make: ", ref,
                                   " must be a reference"););
      return;
    }

    // Update region counting
    auto num_refs = m_rgn_counting_dom[rgn];
    m_rgn_counting_dom.set(rgn, num_refs.increment());

    if (Params::allocation_sites) {
      // Associate allocation site as to ref
      m_alloc_site_dom.set(ref, as);
    }
    
    CRAB_LOG("region", crab::outs() << "After ref_make(" << ref << "," << rgn << ","
	                            << as << ")=" << *this << "\n";);
  }

  void ref_free(const variable_t &rgn, const variable_t &ref) override {
    crab::CrabStats::count(domain_name() + ".count.ref_free");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_free");

    if (is_bottom()) {
      return;
    }

    assert(rgn.get_type().is_region());      
    assert(ref.get_type().is_reference());

    if (Params::allocation_sites) {
      m_alloc_site_dom -= ref;
    }
    
    // We conservatively mark the region's equivalence class as
    // possibly deallocated
    if (Params::deallocation) {
      if (!m_rgn_dealloc_dom.contains(rgn)) {
        CRAB_LOG("region", CRAB_WARN("lost track of dealloc status of ", rgn,
                                     " in ", m_rgn_dealloc_dom));
      } else {
        m_rgn_dealloc_dom.set(rgn, boolean_value::top());
      }
    }
    
    CRAB_LOG("region", crab::outs() << "After ref_free(" << rgn << "," << ref
                                    << ")=" << *this << "\n";);
  }

  // Read the content of reference ref within rgn. The content is
  // stored in res.
  void ref_load(const variable_t &ref, const variable_t &rgn,
                const variable_t &res) override {
    crab::CrabStats::count(domain_name() + ".count.ref_load");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_load");

    if (is_bottom()) {
      return;
    }

    if (Params::allocation_sites) {
      if (res.get_type().is_reference()) {
	m_alloc_site_dom.set(res, m_alloc_site_dom[rgn]);
      }
    }

    if (Params::tag_analysis) {
      m_tag_env.set(res, m_tag_env[rgn]);
    }
    
    // Non-scalar regions should be handled by ref_load_from_array
    if (!rgn.get_type().is_scalar_region()) {
      CRAB_LOG("region-load", CRAB_WARN("region_domain::ref_load: ", rgn,
                                        " must be a scalar region"););
      m_base_dom -= res;
      return;
    } else {
      if ((rgn.get_type().is_bool_region() && !res.get_type().is_bool()) ||
          (rgn.get_type().is_integer_region() &&
           !res.get_type().is_integer()) ||
          (rgn.get_type().is_real_region() && !res.get_type().is_real()) ||
          (rgn.get_type().is_reference_region() &&
           !res.get_type().is_reference())) {
        CRAB_ERROR(domain_name(), "::ref_load: type of lhs ", res, " (",
                   res.get_type(), ") is not compatible with region ", rgn,
                   " (", rgn.get_type(), ")");
      }
    }

    if (!ref.get_type().is_reference()) {
      // This shouldn't happen
      CRAB_LOG("region-load", CRAB_WARN("region_domain::ref_load: ", ref,
                                        " must be a reference"););
      m_base_dom -= res;
      return;
    }

    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("region-load",
               CRAB_WARN("region_domain::ref_load: reference ", ref, " is null."
                         //, " Set to bottom ..."
               ););
      // set_to_bottom();
      m_base_dom -= res;
      return;
    }

    //======================================//
    // Decide whether weak or strong read
    //======================================//
    
    if (rgn.get_type().is_unknown_region()) {
      // TODO: handle unknown regions by casting always to an integer. 
      CRAB_LOG(
          "region-load",
          CRAB_WARN("region_domain::ref_load: skip unknown region ", rgn););
      m_base_dom -= res;
      return;
    }
    
    auto ref_load_aux = [this](const variable_t &res, const variable_t &rgn_var) {
      // At this point references and region of references have been
      // translated to integers and region of integers, respectively.
      // UPDATE: this is not true anymore since we don't use ghost variables.
      // if (rgn_var.get_type() != base_res.get_type()) {
      //   CRAB_ERROR("region_domain::ref_load: ", "Type of region ", rgn,
      //              " does not match with ", base_res);
      // }
      if (res.get_type().is_bool()) {
        this->m_base_dom.assign_bool_var(res, rgn_var, false);
      } else if (res.get_type().is_integer() ||
                 res.get_type().is_real() ||
		 res.get_type().is_reference()) {
        this->m_base_dom.assign(res, rgn_var);
      } else {
        // variables of type base_variable_t cannot be REF_TYPE or ARR_TYPE
	// UPDATE: not sure if this true anymore
        // CRAB_ERROR("region_domain::ref_load: unsupported type in ", base_res);
      }
    };
    
    auto num_refs = m_rgn_counting_dom[rgn];
    if (num_refs.is_zero() || num_refs.is_one()) {
      CRAB_LOG("region-load", crab::outs() << "Reading from singleton\n";);
      ref_load_aux(res, rgn);
    } else {
      CRAB_LOG("region-load", crab::outs() << "Reading from non-singleton\n";);
      auto &vfac = const_cast<varname_t *>(&(rgn.name()))->get_var_factory();
      variable_t fresh_region_var(vfac.get(), rgn.get_type());
      m_base_dom.expand(rgn, fresh_region_var);
      ref_load_aux(res, fresh_region_var);
      m_base_dom -= fresh_region_var;
    }

    CRAB_LOG("region-load", crab::outs() << "After " << res << ":="
                                         << "ref_load(" << ref << "," << rgn
                                         << ")=" << *this << "\n";);
  }

  // Write the content of val to the address pointed by ref in region
  // rgn.
  void ref_store(const variable_t &ref, const variable_t &rgn,
                 const variable_or_constant_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.ref_store");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_store");

    if (is_bottom()) {
      return;
    }

    auto forget_region_if_tracked = [this](const variable_t &v) {
      if (v.get_type().is_region() && !v.get_type().is_unknown_region()) {
	m_base_dom -= v;
      }
    };
				    
    bool is_uninitialized_rgn = m_rgn_init_dom[rgn].is_false();
    // We conservatively mark the region as may-initialized
    m_rgn_init_dom.set(rgn, boolean_value::top());

    // Non-scalar regions should be handled by ref_store_to_array
    if (!rgn.get_type().is_scalar_region()) {
      CRAB_LOG("region-store", CRAB_WARN("region_domain::ref_store: ", rgn,
                                         " must be a scalar region"););
      forget_region_if_tracked(rgn);
      return;
    } else {
      if ((rgn.get_type().is_bool_region() && !val.get_type().is_bool()) ||
          (rgn.get_type().is_integer_region() &&
           !val.get_type().is_integer()) ||
          (rgn.get_type().is_real_region() && !val.get_type().is_real()) ||
          (rgn.get_type().is_reference_region() &&
           !val.get_type().is_reference())) {
        CRAB_ERROR(domain_name(), "::ref_store: type of value ", val, " (",
                   val.get_type(), ") is not compatible with region ", rgn,
                   " (", rgn.get_type(), ")");
      }
    }

    if (!ref.get_type().is_reference()) {
      // This shouldn't happen
      CRAB_LOG("region-store", CRAB_WARN("region_domain::ref_store: ", ref,
                                         " must be a reference"););
      forget_region_if_tracked(rgn);      
      return;
    }

    if (is_null_ref(ref).is_true()) {
      CRAB_LOG("region-store",
               CRAB_WARN("region_domain::ref_store: reference ", ref,
                         " is null. " //, " Set to bottom ..."
               ););
      // set_to_bottom();
      forget_region_if_tracked(rgn);      
      return;
    }

    //======================================//
    // Decide whether weak or strong write
    //======================================//
        
    auto ref_store = [&val, this](base_abstract_domain_t &base_dom,
				  const variable_t &base_rgn_var) {
      // At this point region of references has been translated to region of
      // integers. val can be any scalar including a reference.
      // UPDATE: this is not true anymore
      if (val.get_type().is_bool()) {
        if (val.is_constant()) {
          base_dom.assign_bool_cst(
              base_rgn_var,
              (val.is_bool_true() ? linear_constraint_t::get_true()
                                  : linear_constraint_t::get_false()));
        } else {
          base_dom.assign_bool_var(base_rgn_var, val.get_variable(), false);
        }
      } else if (val.get_type().is_integer() || val.get_type().is_real() ||
                 val.get_type().is_reference()) {
        if (val.is_constant()) {
          base_dom.assign(base_rgn_var, val.get_constant());
        } else {
          base_dom.assign(base_rgn_var, val.get_variable());
        }
      } else {
        CRAB_ERROR(domain_name(), "::ref_store: unsupported type ", val.get_type());
                   
      }
    };

    auto num_refs = m_rgn_counting_dom[rgn];
    // A region can have more than one reference but we can still
    // perform a strong update as long as nobody wrote yet in the
    // region. The use of "is_uninitialized_rgn" avoids fooling the
    // abstract domain with code like this:
    //
    //   r1 := make_ref(R);
    //   r2 := make_ref(R);
    //   store_ref(R, r1, v1);
    //   store_ref(R, r2, v2);
    //
    if (is_uninitialized_rgn || num_refs.is_zero() || num_refs.is_one()) {
      /* strong update */
      CRAB_LOG("region-store", crab::outs() << "Performing strong update\n";);

      if (!rgn.get_type().is_unknown_region()) {
	ref_store(m_base_dom, rgn);
      } else {
	// TODO: handle unknown regions by casting to integers
	CRAB_LOG(
          "region-store",
          CRAB_WARN("region_domain::ref_store: skip unknown region ", rgn););
      } 
      
      if (Params::allocation_sites) {
	if (val.get_type().is_reference()) {
	  if (val.is_reference_null()) {
	    // reset to empty set because of strong update
	    m_alloc_site_dom.set(rgn, allocation_sites::bottom());
	  } else {
	    assert(val.is_variable());
	    m_alloc_site_dom.set(rgn, m_alloc_site_dom[val.get_variable()]);
	  }
	}
      }

      if (Params::tag_analysis) {
	if (val.is_variable()) {
	  m_tag_env.set(rgn, m_tag_env[val.get_variable()]);
	}
      }
	
    } else {
      /* weak update */
      CRAB_LOG("region-store", crab::outs() << "Performing weak update\n";);

      if (!rgn.get_type().is_unknown_region()) {
	base_abstract_domain_t tmp(m_base_dom);
	ref_store(tmp, rgn);
	m_base_dom |= tmp;
      } else {
	// TODO: handle unknown regions by casting to integers
	CRAB_LOG( 
          "region-store",
          CRAB_WARN("region_domain::ref_store: skip unknown region ", rgn););
      }
      
      if (Params::allocation_sites) {
	if (val.get_type().is_reference()) {
	  if (val.is_variable()) {
	    m_alloc_site_dom.set(rgn, m_alloc_site_dom[rgn] | m_alloc_site_dom[val.get_variable()]);
	  }
	}
      }
      if (Params::tag_analysis) {
	if (val.is_variable()) {
	  m_tag_env.set(rgn, m_tag_env[rgn] | m_tag_env[val.get_variable()]);
	}
      }
      
    }
    CRAB_LOG("region-store", crab::outs()
                                 << "After ref_store(" << ref << "," << rgn
                                 << "," << val << ")=" << *this << "\n";);
  }

  // Create a new reference ref2 to region rgn2.
  // The reference ref2 is created by adding offset to ref1.
  void ref_gep(const variable_t &ref1, const variable_t &rgn1,
               const variable_t &ref2, const variable_t &rgn2,
               const linear_expression_t &offset) override {
    crab::CrabStats::count(domain_name() + ".count.ref_gep");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_gep");

    auto eval = [this](const linear_expression_t &e) {
      interval_t r = e.constant();
      for (auto p : e) {
        r += p.first * m_base_dom.operator[](p.second);
      }
      return r;
    };

    if (!is_bottom()) {
      m_base_dom.assign(ref2, ref1 + offset);
      if (!(rgn1 == rgn2 && (eval(offset) == (number_t(0))))) {
        // Update region counting
        auto num_refs = m_rgn_counting_dom[rgn2];
        m_rgn_counting_dom.set(rgn2, num_refs.increment());
      }

      if (Params::allocation_sites) {
	m_alloc_site_dom.set(ref2, m_alloc_site_dom[ref1]);
      }

      if (Params::tag_analysis) {
	m_tag_env.set(ref2, m_tag_env[ref1]);
      }
      
      if (Params::deallocation) {
        bool found1 = m_rgn_dealloc_dom.contains(rgn1);
        bool found2 = m_rgn_dealloc_dom.contains(rgn2);
        if (!found1) {
          CRAB_LOG("region", CRAB_WARN("lost track of dealloc status of ", rgn1,
                                       " in ", m_rgn_dealloc_dom));
        }
        if (!found2) {
          CRAB_LOG("region", CRAB_WARN("lost track of dealloc status of ", rgn2,
                                       " in ", m_rgn_dealloc_dom));
        }
	
        if (found1 && found2) {
          m_rgn_dealloc_dom.join(rgn1, rgn2);
        } else if (!found1 && found2) {
	  //m_rgn_dealloc_dom.remove_equiv_class(rgn2);
	} else if (found1 && !found2) {
	  //m_rgn_dealloc_dom.remove_equiv_class(rgn1);	  
	} else {
	  // do nothing since rgn1 and rgn2 are already untracked.
	}
      }
    }

    CRAB_LOG("region", crab::outs()
                           << "After (" << rgn2 << "," << ref2
                           << ") := ref_gep(" << rgn1 << "," << ref1 << " + "
                           << offset << ")=" << *this << "\n";);
  }

  // Treat memory pointed by ref  as an array and perform an array load.
  void ref_load_from_array(const variable_t &lhs,
                           const variable_t &ref /*unused*/,
                           const variable_t &rgn,
                           const linear_expression_t &index,
                           const linear_expression_t &elem_size) override {
    crab::CrabStats::count(domain_name() + ".count.ref_load_from_array");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_load_from_array");

    if (is_bottom()) {
      return;
    }

    /***
     * These syntactic checks are only needed if the abstract domain
     * API is called directly.
     **/
    if (!lhs.get_type().is_integer() && !lhs.get_type().is_bool() &&
        !lhs.get_type().is_real()) {
      CRAB_LOG("region", CRAB_WARN("region_domain::ref_load_from_array: ", lhs,
                                   " must be bool/int/real type"););
      m_base_dom -= lhs;
    }

    if (!rgn.get_type().is_array_region()) {
      CRAB_LOG("region", CRAB_WARN("region_domain::ref_load_from_array: ", rgn,
                                   " must be an array region"););
      m_base_dom -= lhs;
      return;
    }

    auto num_refs = m_rgn_counting_dom[rgn];
    if (num_refs.is_zero() || num_refs.is_one()) {
      CRAB_LOG("region", crab::outs() << "Reading from singleton\n";);
      m_base_dom.array_load(lhs, rgn, elem_size, index);
    } else {
      // TODO: get the bitwidth from index
      auto &vfac = const_cast<varname_t *>(&(rgn.name()))->get_var_factory();
      variable_t unknown_index(vfac.get(), crab::INT_TYPE, 32);
      m_base_dom.array_load(lhs, rgn, elem_size, unknown_index);
    }

    // TODO: tag analysis
  }

  // Treat region as an array and perform an array store.
  void ref_store_to_array(const variable_t &ref /*unused*/,
                          const variable_t &rgn,
                          const linear_expression_t &index,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.ref_store_to_array");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_store_to_array");

    if (is_bottom()) {
      return;
    }

    /***
     * These syntactic checks are only needed if the abstract domain
     * API is called directly.
     **/
    if (!rgn.get_type().is_array_region()) {
      CRAB_LOG("region", CRAB_WARN("region_domain::ref_store_to_array: ", rgn,
                                   " must be an array region"););
      return;
    }

    // We conservatively mark the region as may-initialized
    m_rgn_init_dom.set(rgn, boolean_value::top());

    auto num_refs = m_rgn_counting_dom[rgn];
    if (num_refs.is_zero() || num_refs.is_one()) {
      CRAB_LOG("region", crab::outs() << "Reading from singleton\n";);
      m_base_dom.array_store(rgn, elem_size, index, val, false);

    } else {
      // TODO: get the bitwidth from index
      auto &vfac = const_cast<varname_t *>(&(rgn.name()))->get_var_factory();      
      variable_t unknown_index(vfac.get(), crab::INT_TYPE, 32);      
      m_base_dom.array_store(rgn, elem_size, unknown_index, val, false);
    }

    // TODO: tag analysis    
  }

  void ref_assume(const reference_constraint_t &ref_cst) override {
    crab::CrabStats::count(domain_name() + ".count.ref_assume");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_assume");

    if (!is_bottom()) {
      if (ref_cst.is_tautology()) {
        return;
      }
      if (ref_cst.is_contradiction()) {
        set_to_bottom();
        return;
      }

      if (Params::allocation_sites) {      
	if (ref_cst.is_equality()) {
	  if (ref_cst.is_binary()) {
	    auto lhs_as = m_alloc_site_dom[ref_cst.lhs()];
	    auto rhs_as = m_alloc_site_dom[ref_cst.rhs()];
	    allocation_sites inter = lhs_as & rhs_as;
	    // FIXME: the lhs or rhs operands might have empty
	    // allocation sites if the program environment is not
	    // modeled precisely enough. In those cases, they might have
	    // no allocation sites but we don't want to set the whole
	    // abstract state to bottom.
	    if (!lhs_as.is_bottom() && !rhs_as.is_bottom() &&
		inter.is_bottom()) {
	      // if they do not have any common allocation site then
	      // they cannot be the same address.
	      // 
	      // However, if ref_cst is a disequality we cannot say
	      // anything by looking only at allocation sites. Even if
	      // both references point to singletons they can still be
	      // different.
	      set_to_bottom();
	      return;
	    }
	  }
	}
      }
      
      auto lin_csts = convert_ref_cst_to_linear_cst(ref_cst);
      m_base_dom += lin_csts;
      m_is_bottom = m_base_dom.is_bottom();
    }
    CRAB_LOG("region",
             crab::outs() << "ref_assume(" << ref_cst << ")" << *this << "\n";);
  }

  void ref_to_int(const variable_t &rgn, const variable_t &ref_var,
                  const variable_t &int_var) override {
    crab::CrabStats::count(domain_name() + ".count.ref_to_int");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_to_int");

    assert(ref_var.get_type().is_reference());
    assert(int_var.get_type().is_integer());

    if (!is_bottom()) {
      m_base_dom.assign(int_var, ref_var);
      if (Params::tag_analysis) {
	m_tag_env.set(int_var, m_tag_env[ref_var]);
      }
    }
  }

  void int_to_ref(const variable_t &int_var, const variable_t &rgn,
                  const variable_t &ref_var) override {
    crab::CrabStats::count(domain_name() + ".count.int_to_ref");
    crab::ScopedCrabStats __st__(domain_name() + ".int_to_ref");

    assert(ref_var.get_type().is_reference());
    assert(int_var.get_type().is_integer());

    if (!is_bottom()) {
      m_base_dom.assign(ref_var, int_var);

      if (Params::allocation_sites) {
	m_alloc_site_dom -= ref_var;
      }

      if (Params::tag_analysis) {
	m_tag_env.set(ref_var, m_tag_env[int_var]);
      }
      
      // Update region counting
      auto num_refs = m_rgn_counting_dom[rgn];
      m_rgn_counting_dom.set(rgn, num_refs.increment());
    }
  }

  // This default implementation is expensive because it will call the
  // join.
  DEFAULT_SELECT_REF(region_domain_t)
  
  // arithmetic operations
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {

    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, x, y, z);
      if (Params::tag_analysis) {
	m_tag_env.set(x, m_tag_env[y] | m_tag_env[z]);
      }
    }
  }
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, x, y, k);
      if (Params::tag_analysis) {
	m_tag_env.set(x, m_tag_env[y]);
      }      
    }
  }
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (!is_bottom()) {
      m_base_dom.assign(x, e);
      if (Params::tag_analysis) {
	merge_tags(x, e.variables());
      }
    }
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
	      const linear_expression_t &e1,  const linear_expression_t &e2) override {
    if (!is_bottom()) {
      m_base_dom.select(lhs, cond, e1, e2);
      if (Params::tag_analysis) {
	tag_set tags = tag_set::bottom();
	for (auto const &v: e1.variables()) {
	  tags = tags | m_tag_env[v];
	}
	for (auto const &v: e2.variables()) {
	  tags = tags | m_tag_env[v];
	}
	m_tag_env.set(lhs, tags);
      }      
    }
  }
  
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");

    if (!is_bottom()) {
      m_base_dom.backward_assign(x, e, invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }      
    }
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    if (!is_bottom()) {
      m_base_dom.backward_apply(op, x, y, z, invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    if (!is_bottom()) {
      m_base_dom.backward_apply(op, x, y, z, invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }
  void operator+=(const linear_constraint_system_t &csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    if (!is_bottom()) {
      m_base_dom += csts;
      m_is_bottom = m_base_dom.is_bottom();
    }
  }

  // not part of the numerical_domains api but it should be
  void set(const variable_t &x, interval_t intv) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (!is_bottom()) {
      if (intv.is_bottom()) {
        set_to_bottom();
      } else {
        m_base_dom.set(x, intv);
      }
      
      if (Params::tag_analysis) {
	m_tag_env.set(x, tag_env_t::value_type::bottom()); 
      }      
    }
  }

  interval_t operator[](const variable_t &x) override {
    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    } else {
      return m_base_dom[x];
    }
  }

  // int cast operations and bitwise operations
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, dst, src);
      if (Params::tag_analysis) {
	m_tag_env.set(dst, m_tag_env[src]);
      }            
    }
  }
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, x, y, z);
      if (Params::tag_analysis) {
	m_tag_env.set(x, m_tag_env[y] | m_tag_env[z]);
      }            
    }
  }
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      m_base_dom.apply(op, x, y, z);
      if (Params::tag_analysis) {
	m_tag_env.set(x, m_tag_env[y]);
      }            
    }
  }

  // boolean operations
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (!is_bottom()) {
      m_base_dom.assign_bool_cst(lhs, rhs);
      if (Params::tag_analysis) {
	merge_tags(lhs, rhs.variables());
      }            
    }
  }

  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    if (!is_bottom()) {

      if (Params::allocation_sites) {      
	if (rhs.is_equality()) {
	  if (rhs.is_binary()) {
	    allocation_sites inter =
	      m_alloc_site_dom[rhs.lhs()] & m_alloc_site_dom[rhs.rhs()];
	    if (inter.is_bottom()) {
	      // if they do not have any common allocation site then
	      // they cannot be the same address.
	    auto false_cst =  linear_constraint_t::get_false();
	    m_base_dom.assign_bool_cst(lhs, false_cst);
	    return;
	    }
	  }
	}
      }
      
      auto lin_rhs = convert_ref_cst_to_linear_cst(rhs);
      m_base_dom.assign_bool_cst(lhs, lin_rhs);
      if (Params::tag_analysis) {
	merge_tags(lhs, rhs.variables());
      }            
    }
  }

  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_var");

    if (!is_bottom()) {
      m_base_dom.assign_bool_var(lhs, rhs, is_not_rhs);
      if (Params::tag_analysis) {
	m_tag_env.set(lhs, m_tag_env[rhs]);
      }            
    }
  }
  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".apply_binary_bool");

    if (!is_bottom()) {
      m_base_dom.apply_binary_bool(op, x, y, z);
      if (Params::tag_analysis) {
	m_tag_env.set(x, m_tag_env[y] | m_tag_env[z]);
      }            
    }
  }
  void assume_bool(const variable_t &v, bool is_negated) override {
    crab::CrabStats::count(domain_name() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".assume_bool");

    if (!is_bottom()) {
      m_base_dom.assume_bool(v, is_negated);
      m_is_bottom = m_base_dom.is_bottom();
    }
  }
  void select_bool(const variable_t &lhs, const variable_t &cond,
		   const variable_t &b1, const variable_t &b2) override {
    if (!is_bottom()) {
      m_base_dom.select_bool(lhs, cond, b1, b2);
      if (Params::tag_analysis) {
	// it can be improve if we know if cond is true or false
	m_tag_env.set(lhs, m_tag_env[b1] | m_tag_env[b2]);
      }                  
    }
  }  
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_cst");

    if (!is_bottom()) {
      m_base_dom.backward_assign_bool_cst(lhs, rhs, invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }      
    }
  }

  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_cst");

    if (!is_bottom()) {
      auto lin_rhs = convert_ref_cst_to_linear_cst(rhs);
      m_base_dom.backward_assign_bool_cst(lhs, lin_rhs,
                                          invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }

  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign_bool_var");

    if (!is_bottom()) {
      m_base_dom.backward_assign_bool_var(lhs, rhs, is_not_rhs, invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }
  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply_binary_bool");

    if (!is_bottom()) {
      m_base_dom.backward_apply_binary_bool(op, x, y, z, invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }

  // array operations
  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx,
                  const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {
    crab::CrabStats::count(domain_name() + ".count.array_init");
    crab::ScopedCrabStats __st__(domain_name() + ".array_init");

    if (!is_bottom()) {
      m_base_dom.array_init(a, elem_size, lb_idx, ub_idx, val);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }      
      
    }
  }
  void array_load(const variable_t &lhs, const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {
    crab::CrabStats::count(domain_name() + ".count.array_load");
    crab::ScopedCrabStats __st__(domain_name() + ".array_load");

    if (!is_bottom()) {
      m_base_dom.array_load(lhs, a, elem_size, i);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &v,
                   bool is_strong_update) override {
    crab::CrabStats::count(domain_name() + ".count.array_store");
    crab::ScopedCrabStats __st__(domain_name() + ".array_store");

    if (!is_bottom()) {
      m_base_dom.array_store(a, elem_size, i, v, is_strong_update);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }
  void array_store_range(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &j,
                         const linear_expression_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.array_store_range");
    crab::ScopedCrabStats __st__(domain_name() + ".array_store_range");

    if (!is_bottom()) {
      m_base_dom.array_store_range(a, elem_size, i, j, v);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }
  void array_assign(const variable_t &lhs, const variable_t &rhs) override {
    crab::CrabStats::count(domain_name() + ".count.array_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".array_assign");

    if (!is_bottom()) {
      m_base_dom.array_assign(lhs, rhs);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }
  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_init");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_init");

    if (!is_bottom()) {
      m_base_dom.backward_array_init(a, elem_size, lb_idx,
                                     ub_idx, val, invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_load");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_load");

    if (!is_bottom()) {
      m_base_dom.backward_array_load(lhs, a, elem_size, i, invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_store");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_store");

    if (!is_bottom()) {
      m_base_dom.backward_array_store(a, elem_size, i, v,
                                      is_strong_update, invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }
  void backward_array_store_range(const variable_t &a,
                                  const linear_expression_t &elem_size,
                                  const linear_expression_t &i,
                                  const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_store_range");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_store_range");

    if (!is_bottom()) {
      m_base_dom.backward_array_store_range(a, elem_size, i, j, v, 
                                            invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             const region_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_array_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_array_assign");

    if (!is_bottom()) {
      m_base_dom.backward_array_assign(lhs, rhs, invariant.m_base_dom);
      if (Params::tag_analysis) {
	// TODO tag analysis
      }            
    }
  }

  void forget(const variable_vector_t &variables) override {
    if (is_bottom() || is_top()) {
      return;
    }

    for (auto &v : variables) {
      if (Params::allocation_sites) {
	if (v.get_type().is_reference() || v.get_type().is_region()) {
	  m_alloc_site_dom -= v;
	}
      }
      if (v.get_type().is_region()) {
        m_rgn_counting_dom -= v;
        m_rgn_init_dom -= v;
        if (Params::deallocation) {
          m_rgn_dealloc_dom.forget(v);
        }	
      }
      if (Params::tag_analysis) {
	m_tag_env -= v;
      }
    }
    m_base_dom.forget(variables);
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }

    variable_vector_t sorted_variables(variables.begin(), variables.end());
    std::sort(sorted_variables.begin(), sorted_variables.end());

    // -- project in the base domain
    m_base_dom.project(sorted_variables);

    m_rgn_counting_dom.project(sorted_variables);
    m_rgn_init_dom.project(sorted_variables);
    if (Params::allocation_sites) {
      m_alloc_site_dom.project(sorted_variables);
    }
    if (Params::deallocation) {
      m_rgn_dealloc_dom.project(sorted_variables);
    }
    if (Params::tag_analysis) {
      m_tag_env.project(sorted_variables);
    }    
  }

  void normalize() override {}
  void minimize() override {}

  void expand(const variable_t &x, const variable_t &new_x) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }
    CRAB_ERROR("region_domain::expand not implemented");
  }

  // The code might be simpler and more efficient if m_var_map is a
  // std::map.  However, this would slow down lookups.  The
  // expectation is that this operation should not been called often
  // since the region domain is usually at the root of the hierarchy
  // of domains.
  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    if (is_bottom() || is_top()) {
      return;
    }
    
    if (from.size() != to.size()) {
      CRAB_ERROR(domain_name(), "::rename different lengths");
    }

    m_base_dom.rename(from, to);
    m_rgn_counting_dom.rename(from, to);
    m_rgn_init_dom.rename(from, to);
    if (Params::allocation_sites) {
      m_alloc_site_dom.rename(from, to);
    }
    if (Params::deallocation) {
      m_rgn_dealloc_dom.rename(from, to);
    }
    if (Params::tag_analysis) {
      m_tag_env.rename(from, to);
    }
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    //=================================================================//
    //       Special intrinsics supported by the region domain
    //=================================================================//
    // ---Use-after-free analysis----
    // 
    // b := is_unfreed_or_null(rgn, ref)
    //       b is true if the reference does not point to a freed
    //       memory object or it is null.
    // unfreed_or_null(rgn, ref)
    //       ensures that the reference does not point to a freed memory object
    // 
    // ---Nullity analysis---
    // 
    // nonnull(ref)
    //       ensures that the reference is not null
    // 
    // ---Tag analysis---
    /// 
    // add_tag(rgn, ref, TAG)
    //       Add TAG to the set of tags already associated with the
    //       data pointed by ref within rgn.
    // b := does_not_have_tag(rgn, ref, TAG)
    //       b is true if the data pointed by ref within rgn does
    //       *definitely* not have tag TAG.
    //=================================================================//
    auto error_if_not_arity = [&name, &inputs, &outputs](unsigned num_inputs,
                                                         unsigned num_outputs) {
      if (inputs.size() != num_inputs || outputs.size() != num_outputs) {
        CRAB_ERROR("Intrinsics ", name, " unexpected number of parameters");
      }
    };
    auto error_if_not_variable = [&name](const variable_or_constant_t &vc) {
      if (!vc.is_variable()) {
	CRAB_ERROR("Intrinsics ", name, " expected a variable input");
      }
    };      
    auto error_if_not_constant = [&name](const variable_or_constant_t &vc) {
      if (!vc.is_constant()) {
	CRAB_ERROR("Intrinsics ", name, " expected a constant input");
      }
    };      
    
    auto error_if_not_bool = [&name](const variable_t &var) {
      if (!var.get_type().is_bool()) {
        CRAB_ERROR("Intrinsics ", name, " parameter ", var, " should be Bool");
      }
    };
    auto error_if_not_rgn = [&name](const variable_t &var) {
      if (!var.get_type().is_region()) {
        CRAB_ERROR("Intrinsics ", name, " parameter ", var,
                   " should be a region");
      }
    };
    auto error_if_not_ref = [&name](const variable_t &var) {
      if (!var.get_type().is_reference()) {
        CRAB_ERROR("Intrinsics ", name, " parameter ", var,
                   " should be a reference");
      }
    };

    auto set_bool_var_to_true = [this](const variable_t &bool_var) {
      /// Require that the base domain can reason about booleans
      operator-=(bool_var);
      assume_bool(bool_var, false /*not negated*/);
    };
    
    if (is_bottom()) {
      return;
    }
    
    if (name == "is_unfreed_or_null") {
      if (Params::deallocation) {
	error_if_not_arity(2, 1);
	error_if_not_variable(inputs[0]);
	error_if_not_variable(inputs[1]);      
	variable_t rgn = inputs[0].get_variable();
	variable_t ref = inputs[1].get_variable();
	variable_t bv = outputs[0];
	error_if_not_rgn(rgn);
	error_if_not_ref(ref);
	error_if_not_bool(bv);
	// the reference is definitely null
	if (is_null_ref(ref).is_true()) {
	  set_bool_var_to_true(bv);
	} else {
	  bool is_allocated = (m_rgn_dealloc_dom.contains(rgn) &&
			       m_rgn_dealloc_dom[rgn].is_false());
	  /// the reference belongs to a memory region that doesn't have
	  /// any deallocated memory object.
	  if (is_allocated) {
	    set_bool_var_to_true(bv);
	  }
	}
      }
    } else if (name == "unfreed_or_null") {
      if (Params::deallocation) {
	error_if_not_arity(2, 0);
	error_if_not_variable(inputs[0]);
	error_if_not_variable(inputs[1]);
	variable_t rgn = inputs[0].get_variable();
	variable_t ref = inputs[1].get_variable();
	error_if_not_rgn(rgn);
	error_if_not_ref(ref);
	if (!is_null_ref(ref).is_true()) {
	  // We can only mark the region as "deallocated" if it doesn't
	  // have any reference yet. Even if it has only one we cannot
	  // tell whether it's the same ref than ref.
	  auto num_refs = m_rgn_counting_dom[rgn];
	  if (num_refs.is_zero()) {
	    m_rgn_dealloc_dom.set(rgn, boolean_value::get_false());
	  }
	}
      }
    } else if (name == "nonnull") {
      error_if_not_arity(1, 0);
      error_if_not_variable(inputs[0]);      
      variable_t ref = inputs[0].get_variable();
      error_if_not_ref(ref);
      auto nonnull_cst = reference_constraint_t::mk_gt_null(ref);
      ref_assume(nonnull_cst);
    } else if (name == "add_tag") {
      if (Params::tag_analysis) {
	error_if_not_arity(3, 0);
	error_if_not_variable(inputs[0]);
	error_if_not_variable(inputs[1]);
	error_if_not_constant(inputs[2]);
	variable_t rgn = inputs[0].get_variable();
	error_if_not_rgn(rgn);
	variable_t ref = inputs[1].get_variable();
	error_if_not_ref(ref);
	tag_t tag(inputs[2].get_constant());

	// We ignore ref and merge the tag with the tags already in
	// the region
	m_tag_env.set(rgn, m_tag_env[rgn] | tag_set(tag));
      }
    } else if (name == "does_not_have_tag") {
      if (Params::tag_analysis) {
	error_if_not_arity(3, 1);
	error_if_not_variable(inputs[0]);
	error_if_not_variable(inputs[1]);
	error_if_not_constant(inputs[2]);
	error_if_not_bool(outputs[0]);
	variable_t rgn = inputs[0].get_variable();
	error_if_not_rgn(rgn);	
	variable_t ref = inputs[1].get_variable();
	error_if_not_ref(ref);	
	tag_t tag(inputs[2].get_constant());
	variable_t bv = outputs[0];
	tag_set tags = m_tag_env[rgn];
	if (!(tag_set(tag) <= tags)) {
	  set_bool_var_to_true(bv);
	} else {
	  operator-=(bv);			  
	} 
      }
    } else {
      // pass the intrinsics to the base domain
      m_base_dom.intrinsic(name, inputs, outputs);
    }
  }

  void backward_intrinsic(std::string name,
			  const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const region_domain_t &invariant) override {
    if (!is_bottom()) {
      m_base_dom.backward_intrinsic(name, inputs, outputs, invariant.m_base_dom);
    }
  }
  /* end intrinsics operations */

  /* Begin abstract_domain_results_api */

  // Return a 3-valued boolean. If true then ref is definitely
  // null. If false then ref is definitely non-null. Otherwise, we do
  // not know.
  boolean_value is_null_ref(const variable_t &ref) override {
    if (is_bottom()) {
      return boolean_value::bottom();
    }
    
    if (!ref.get_type().is_reference()) {
      return boolean_value::get_false();
    }
    
    interval_t ival = m_base_dom[ref];
    number_t zero(0);
    
    if (!(interval_t(zero) <= ival)) {
      return boolean_value::get_false();
    }
    
    boost::optional<number_t> x = ival.lb().number();
    boost::optional<number_t> y = ival.ub().number();
    if (x && y && *x == zero && *y == zero) {
      return boolean_value::get_true();
    }
    
    return boolean_value::top();
  }

  virtual bool get_allocation_sites(const variable_t &ref,
				    std::vector<allocation_site> &alloc_sites) override {
    allocation_sites out = m_alloc_site_dom[ref];
    
    if (out.is_top() || out.is_bottom()) {
      alloc_sites.clear();
      return false;
    }
    for (auto it = out.begin(), et = out.end(); it!=et; ++it) {
      alloc_sites.push_back(*it);
    }
    return true;
  }

  virtual bool get_tags(const variable_t &rgn,
			const variable_t &ref /*unused*/,
			std::vector<uint64_t> &tags) override {
    if (!Params::tag_analysis ||
	(!ref.get_type().is_reference() || !rgn.get_type().is_region())) {
      return false;
    }

    tag_set tag_set = m_tag_env[rgn];    
    
    if (tag_set.is_top()) {
      return false;
    }
    
    if (!tag_set.is_bottom()) {
      for (auto it = tag_set.begin(), et = tag_set.end(); it!=et; ++it) {
	tags.push_back((*it).index());
      }
    }
    return true;
  }
  
  /* End abstract_domain_results_api */
  
  // Only integer or boolean variables are exposed.
  //
  // Variables that shadow memory regions and reference variables are
  // ignored.
  linear_constraint_system_t to_linear_constraint_system() const override {
    if (is_bottom()) {
      return linear_constraint_t::get_false();
    } else if (is_top()) {
      return linear_constraint_t::get_true();
    } else {
      return filter_nonscalar_vars(
		 std::move(m_base_dom.to_linear_constraint_system()));
    }
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_ERROR("region_domain::to_disjunctive_linear_constraint_system not "
               "implemented");
  }

  std::string domain_name() const override {
    return "FastRegionDomain(" + m_base_dom.domain_name() + ")";
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "{}";
    } else {
      CRAB_LOG("region-print",
	       o << "(RgnCounter=" << m_rgn_counting_dom << ","
	       << "RgnInit=" << m_rgn_init_dom;
	       if (Params::allocation_sites) {
		 o << "," << "AllocSites=" << m_alloc_site_dom;
	       }
	       if (Params::deallocation) {
		 o << "," << "RgnDealloc=" << m_rgn_dealloc_dom;
	       }
	       if (Params::tag_analysis) {
		 o << "," << "Tags=" << m_tag_env;		 
	       }
	       o << "," << "BaseDom=" << m_base_dom << ")\n";);

      o << "(RgnCounter=" << m_rgn_counting_dom << ",BaseDom=" << m_base_dom << ")";
    }
  }
}; // class region_without_ghost_domain

template <typename Dom, typename Params>
struct abstract_domain_traits<region_without_ghost_domain<Dom, Params>> {
  using number_t = typename Dom::number_t;
  using varname_t = typename Dom::varname_t;
};
  
} // end namespace domains
} // end namespace crab
