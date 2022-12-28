#pragma once

#include <crab/support/os.hpp>
#include <crab/support/debug.hpp>

#include <vector>

namespace crab {
namespace domains {

class crab_domain_params;

class elina_domain_params {
  bool m_use_tree_expressions;

  friend class crab_domain_params;
public:
  elina_domain_params()
    : m_use_tree_expressions(true) {}
  elina_domain_params(bool use_tree_expressions)
    : m_use_tree_expressions(use_tree_expressions) {}

  
  bool elina_use_tree_expressions() const {
    return m_use_tree_expressions;
  }
  void update_params(const elina_domain_params &params);
  void write(crab::crab_os &o) const;  
};

class array_adaptive_domain_params {
  // options for array smashing
  bool m_is_smashable;
  bool m_smash_at_nonzero_offset;
  // -- if the number of the cells is larger than this number then the
  //    array is not smashed.
  unsigned m_max_smashable_cells;
  // options for array expansion
  // -- maximum size of the array: beyond this number the domain will
  //    try to smash the array if the smashing options allow that.
  unsigned m_max_array_size;
  
  friend class crab_domain_params;  
public:
  array_adaptive_domain_params()
    : m_is_smashable(true),
      m_smash_at_nonzero_offset(true),
      m_max_smashable_cells(64),
      m_max_array_size(64) {
    if (m_max_smashable_cells > m_max_array_size) {
      CRAB_ERROR("array_adaptive_domain params: max_smashable_cells must be <= max_array_size");            
    }
  }
  array_adaptive_domain_params(bool is_smashable,
			       bool smash_at_nonzero_offset,
			       unsigned max_smashable_cells,
			       unsigned max_array_size)
    : m_is_smashable(is_smashable),
      m_smash_at_nonzero_offset(smash_at_nonzero_offset),
      m_max_smashable_cells(max_smashable_cells),
      m_max_array_size(max_array_size) {
    if (m_max_smashable_cells > m_max_array_size) {    
      CRAB_ERROR("array_adaptive_domain params: max_smashable_cells must be <= max_array_size");      
    }
  }
  
  bool array_adaptive_is_smashable() const {
    return m_is_smashable;
  }
  bool array_adaptive_smash_at_nonzero_offset() const {
    return m_smash_at_nonzero_offset;
  }
  unsigned array_adaptive_max_smashable_cells() const{
    return m_max_smashable_cells;
  }
  unsigned array_adaptive_max_array_size() const {
    return m_max_array_size;
  }
  void update_nonsmashable_params();
  void update_params(const array_adaptive_domain_params &p);
  void write(crab::crab_os &o) const;
};

class boxes_domain_params {
  unsigned m_ldd_size;
  int m_convexify_threshold;
  bool m_dynamic_reordering;

  friend class crab_domain_params;  
public:
  boxes_domain_params()
    : m_ldd_size(3000),
      m_convexify_threshold(-1),
      m_dynamic_reordering(false) {
  }
  boxes_domain_params(unsigned ldd_size,
		      int convexify_threshold,
		      bool dynamic_reordering)
    : m_ldd_size(ldd_size),
      m_convexify_threshold(convexify_threshold),
      m_dynamic_reordering(dynamic_reordering) {
  }
  
  unsigned boxes_ldd_size() const {
    return m_ldd_size;
  }
  int boxes_convexify_threshold() const {
    return m_convexify_threshold;
  }
  bool boxes_dynamic_reordering() const {
    return m_dynamic_reordering;
  }
  void update_params(const boxes_domain_params &p);
  void write(crab::crab_os &o) const;  
};

class powerset_domain_params {
  bool m_exact_meet;
  // Smash if the number of disjunctions exceeds this threshold.
  unsigned m_max_disjuncts;

  friend class crab_domain_params;  
public:
  powerset_domain_params()
    : m_exact_meet(false),
      m_max_disjuncts(99999) {
  }
  powerset_domain_params(bool exact_meet,
			 unsigned max_disjuncts)
    : m_exact_meet(exact_meet),
      m_max_disjuncts(max_disjuncts) {
  }
  
  bool powerset_exact_meet() const {
    return m_exact_meet;
  }
  unsigned powerset_max_disjuncts() const {
    return m_max_disjuncts;
  }
  void update_params(const powerset_domain_params &p);
  void write(crab::crab_os &o) const;  
};

class region_domain_params {
  // keep track of allocation sites
  bool m_allocation_sites;
  // reason about is_unfreed_or_null and unfreed_or_null intrinsics
  bool m_deallocation;
  // reason about add_tag and does_not_have_tag intrinsics
  bool m_tag_analysis;
  // reason about is_dereferenceable intrinsics
  bool m_is_dereferenceable;
  // reason about unknown regions
  bool m_skip_unknown_regions;

  friend class crab_domain_params;  
public:
  region_domain_params()
    : m_allocation_sites(true),
      m_deallocation(false),
      m_tag_analysis(true),
      m_is_dereferenceable(false),
      m_skip_unknown_regions(true) {
  }
  region_domain_params(bool allocation_sites,
		       bool deallocation,
		       bool tag_analysis,
		       bool is_dereferenceable,
		       bool skip_unknown_regions)
    : m_allocation_sites(allocation_sites),
      m_deallocation(deallocation),
      m_tag_analysis(tag_analysis),
      m_is_dereferenceable(is_dereferenceable),
      m_skip_unknown_regions(skip_unknown_regions) {
  }

  bool region_allocation_sites() const {
    return m_allocation_sites;
  }
  bool region_deallocation() const {
    return m_deallocation;
  }
  bool region_tag_analysis() const {
    return m_tag_analysis;
  }
  bool region_is_dereferenceable() const {
    return m_is_dereferenceable;
  }
  bool region_skip_unknown_regions() const {
    return m_skip_unknown_regions;
  }
  void update_params(const region_domain_params& p);
  void write(crab::crab_os &o) const;
};

class zones_domain_params {
  bool m_chrome_dijkstra;
  bool m_widen_restabilize;
  bool m_special_assign;
  bool m_close_bounds_inline;

  friend class crab_domain_params;  
public:
  zones_domain_params()
    : m_chrome_dijkstra(true),
      m_widen_restabilize(true),
      m_special_assign(true),
      m_close_bounds_inline(false) {
  }
  zones_domain_params(bool chrome_dijkstra,
		      bool widen_restabilize,
		      bool special_assign,
		      bool close_bounds_inline)
    : m_chrome_dijkstra(chrome_dijkstra),
      m_widen_restabilize(widen_restabilize),
      m_special_assign(special_assign),
      m_close_bounds_inline(close_bounds_inline) {
  }
  
  bool zones_chrome_dijkstra() const {
    return m_chrome_dijkstra;
  }
  bool zones_widen_restabilize() const {
    return m_widen_restabilize;
  }
  bool zones_special_assign() const {
    return m_special_assign;
  }
  bool zones_close_bounds_inline() const {
    return m_close_bounds_inline;
  }
  void update_params(const zones_domain_params& p);
  void write(crab::crab_os &o) const;
};

class oct_domain_params {
  bool m_chrome_dijkstra;
  bool m_widen_restabilize;
  bool m_special_assign;
  bool m_close_bounds_inline;

  friend class crab_domain_params;  
public:
  oct_domain_params()
    : m_chrome_dijkstra(true),
      m_widen_restabilize(true),
      m_special_assign(true),
      m_close_bounds_inline(false) {
  }
  oct_domain_params(bool chrome_dijkstra,
		      bool widen_restabilize,
		      bool special_assign,
		      bool close_bounds_inline)
    : m_chrome_dijkstra(chrome_dijkstra),
      m_widen_restabilize(widen_restabilize),
      m_special_assign(special_assign),
      m_close_bounds_inline(close_bounds_inline) {
  }
  
  bool oct_chrome_dijkstra() const {
    return m_chrome_dijkstra;
  }
  bool oct_widen_restabilize() const {
    return m_widen_restabilize;
  }
  bool oct_special_assign() const {
    return m_special_assign;
  }
  bool oct_close_bounds_inline() const {
    return m_close_bounds_inline;
  }
  void update_params(const oct_domain_params& p);
  void write(crab::crab_os &o) const;
};

class fixed_tvpi_domain_params {
  std::vector<unsigned> m_coefficients;

  friend class crab_domain_params;  
public:
  fixed_tvpi_domain_params() {}
  
  fixed_tvpi_domain_params(const std::vector<unsigned> &coefficients)
    : m_coefficients(coefficients) {}

  const std::vector<unsigned>& coefficients() const;
  std::vector<unsigned>& coefficients();  
  void update_params(const fixed_tvpi_domain_params& p);
  void write(crab::crab_os &o) const;
};

class crab_domain_params: public elina_domain_params,
			  public array_adaptive_domain_params,
			  public boxes_domain_params,
			  public powerset_domain_params,
			  public region_domain_params,
			  public zones_domain_params,
			  public oct_domain_params,
			  public fixed_tvpi_domain_params {
public:
  // To resolve ambiguous name lookup
  using elina_domain_params::update_params;
  using array_adaptive_domain_params::update_params;
  using boxes_domain_params::update_params;
  using powerset_domain_params::update_params;
  using region_domain_params::update_params;
  using zones_domain_params::update_params;
  using oct_domain_params::update_params;
  using fixed_tvpi_domain_params::update_params;
  
  crab_domain_params()
    : elina_domain_params(),
      array_adaptive_domain_params(),
      boxes_domain_params(),
      powerset_domain_params(),
      region_domain_params(),
      zones_domain_params(),
      oct_domain_params(),
      fixed_tvpi_domain_params() {
  }

  /* Supported parameter strings:

     - elina.use_tree_expressions: bool
     - array_adaptive.is_smashable: bool
     - array_adaptive.smash_at_nonzer_offset: bool
     - array_adaptive.max_smashable_cells: unsigned
     - array_adaptive.max_array_size: unsigned
     - boxes.ldd_size: unsigned
     - boxes.convexify_threshold: int
     - boxes.dynamic_reordering: bool
     - powerset.exact_meet: bool
     - powerset.max_disjuncts: unsigned
     - region.allocations_sites: bool
     - region.deallocation: bool
     - region.tag_analysis: bool
     - region.is_dereferenceable: bool
     - region.skip_unknown_regions: bool
     - zones.chrome_dijkstra: bool
     - zones.widen_restabilize: bool
     - zones.special_assign: bool
     - zones.close_bounds_inline: bool
     - oct.chrome_dijkstra: bool
     - oct.widen_restabilize: bool
     - oct.special_assign: bool
     - oct.close_bounds_inline: bool
     - fixed_tvpi.coefficients: list(unsigned)
   */
  void set_param(const std::string &param, const std::string &val);
  void update_params(const crab_domain_params& p);
  void write(crab::crab_os &o) const;
};

class crab_domain_params_man {
public:
  static crab_domain_params& get() {
    static crab_domain_params m_params;
    return m_params;
  }
};


} // end namespace domains
} // end namespace crab
