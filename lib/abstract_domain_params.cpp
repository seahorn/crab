#include <crab/domains/abstract_domain_params.hpp>
#include <string>
#include <stdexcept>
#include <sstream>

namespace crab {
namespace domains {

void elina_domain_params::update_params(const elina_domain_params &params) {
  m_use_tree_expressions = params.elina_use_tree_expressions();
}

void elina_domain_params::write(crab::crab_os &o) const {
  o << "Elina domains parameters:\n";
  o << "\tuse_tree_expressions="<< m_use_tree_expressions<<"\n";  
  
}

void array_adaptive_domain_params::update_nonsmashable_params() {
  m_is_smashable = false;
  m_smash_at_nonzero_offset = false;
  m_max_smashable_cells = 512;
  m_max_array_size = 512;
}
void array_adaptive_domain_params::update_params(const array_adaptive_domain_params &params) {
  m_is_smashable = params.array_adaptive_is_smashable();
  m_smash_at_nonzero_offset = params.array_adaptive_smash_at_nonzero_offset();
  m_max_smashable_cells = params.array_adaptive_max_smashable_cells();
  m_max_array_size = params.array_adaptive_max_array_size();
}

void array_adaptive_domain_params::write(crab::crab_os &o) const {
  o << "Array adaptive parameters:\n";
  o << "\tis_smashable=" << m_is_smashable << "\n";
  o << "\tsmash_at_nonzero_offset=" << m_smash_at_nonzero_offset << "\n";
  o << "\tmax_smashable_cells=" << m_max_smashable_cells << "\n";
  o << "\tmax_array_size=" << m_max_array_size << "\n";
}

void boxes_domain_params::update_params(const boxes_domain_params &params) {
  m_ldd_size = params.boxes_ldd_size();
  m_convexify_threshold = params.boxes_convexify_threshold();
  m_dynamic_reordering = params.boxes_dynamic_reordering();
}

void boxes_domain_params::write(crab::crab_os &o) const {
  o << "Boxes parameters:\n";
  o << "\tldd_size=" << m_ldd_size << "\n";
  o << "\tconvexify_threshold=" << m_convexify_threshold << "\n";
  o << "\tdynamic_reordering=" << m_dynamic_reordering << "\n";
}

void powerset_domain_params::update_params(const powerset_domain_params &params) {
  m_exact_meet = params.powerset_exact_meet();
  m_max_disjuncts = params.powerset_max_disjuncts();
}

void powerset_domain_params::write(crab::crab_os &o) const {
  o << "Powerset parameters:\n";
  o << "\texact_meet=" << m_exact_meet << "\n";
  o << "\tmax_disjuncts=" << m_max_disjuncts << "\n";
}

void region_domain_params::update_params(const region_domain_params &params) {
  m_allocation_sites = params.region_allocation_sites();
  m_deallocation = params.region_deallocation();
  m_tag_analysis = params.region_tag_analysis();
  m_is_dereferenceable = params.region_is_dereferenceable();
  m_skip_unknown_regions = params.region_skip_unknown_regions(); 
}

void region_domain_params::write(crab::crab_os &o) const {
  o << "Region parameters:\n";
  o << "\tallocation_sites=" << m_allocation_sites << "\n";
  o << "\tdeallocation=" << m_deallocation << "\n";
  o << "\ttag_analysis=" << m_tag_analysis << "\n";
  o << "\tis_dereferenceable=" << m_is_dereferenceable << "\n";
  o << "\tskip_unknown_regions=" << m_skip_unknown_regions << "\n";
}

void zones_domain_params::update_params(const zones_domain_params &params) {
  m_chrome_dijkstra = params.zones_chrome_dijkstra();
  m_widen_restabilize = params.zones_widen_restabilize();
  m_special_assign = params.zones_special_assign();
  m_close_bounds_inline = params.zones_close_bounds_inline();
}

void zones_domain_params::write(crab::crab_os &o) const {
  o << "Zones parameters:\n";
  o << "\tchrome_dijkstra=" << m_chrome_dijkstra << "\n";
  o << "\twiden_restabilize=" << m_widen_restabilize << "\n";
  o << "\tspecial_assign=" << m_special_assign << "\n";
  o << "\tclose_bounds_inline=" << m_close_bounds_inline << "\n";  
}

void oct_domain_params::update_params(const oct_domain_params &params) {
  m_chrome_dijkstra = params.oct_chrome_dijkstra();
  m_widen_restabilize = params.oct_widen_restabilize();
  m_special_assign = params.oct_special_assign();
  m_close_bounds_inline = params.oct_close_bounds_inline();
  
}

void oct_domain_params::write(crab::crab_os &o) const {
  o << "Octagon parameters:\n";
  o << "\tchrome_dijkstra=" << m_chrome_dijkstra << "\n";
  o << "\twiden_restabilize=" << m_widen_restabilize << "\n";
  o << "\tspecial_assign=" << m_special_assign << "\n";
  o << "\tclose_bounds_inline=" << m_close_bounds_inline << "\n";  
}

const std::vector<unsigned>& fixed_tvpi_domain_params::coefficients() const {
  return m_coefficients;
}

 std::vector<unsigned>& fixed_tvpi_domain_params::coefficients() {
  return m_coefficients;
}
  
void fixed_tvpi_domain_params::update_params(const fixed_tvpi_domain_params& p) {
  m_coefficients.clear();
  m_coefficients.insert(m_coefficients.end(),
			p.m_coefficients.begin(), p.m_coefficients.end());
}
  
void fixed_tvpi_domain_params::write(crab::crab_os &o) const {
  o << "Fixed-TVPI parameters:\n";
  o << "\ttracked coefficients={";
  for (unsigned i=0, sz=m_coefficients.size();i<sz;) {
    o << m_coefficients[i];
    ++i;
    if (i < sz) {
      o << ",";
    }
  }
  o << "}\n";
}

  
void crab_domain_params::update_params(const crab_domain_params &p) {
  elina_domain_params e_p(p.elina_use_tree_expressions());
  array_adaptive_domain_params aa_p(p.array_adaptive_is_smashable(),
				    p.array_adaptive_smash_at_nonzero_offset(),
				    p.array_adaptive_max_smashable_cells(),
				    p.array_adaptive_max_array_size());
  boxes_domain_params b_p(p.boxes_ldd_size(),
			  p.boxes_convexify_threshold(),
			  p.boxes_dynamic_reordering());
  powerset_domain_params p_p(p.powerset_exact_meet(),
			     p.powerset_max_disjuncts());
  region_domain_params r_p(p.region_allocation_sites(),
			   p.region_deallocation(),
			   p.region_tag_analysis(),
			   p.region_is_dereferenceable(),
			   p.region_skip_unknown_regions());
  zones_domain_params z_p(p.zones_chrome_dijkstra(),
			  p.zones_widen_restabilize(),
			  p.zones_special_assign(),
			  p.zones_close_bounds_inline());
  oct_domain_params o_p(p.oct_chrome_dijkstra(),
			p.oct_widen_restabilize(),
			p.oct_special_assign(),
			p.oct_close_bounds_inline());
  fixed_tvpi_domain_params tvpi_p(p.coefficients());
  
  elina_domain_params::update_params(e_p);
  array_adaptive_domain_params::update_params(aa_p);
  boxes_domain_params::update_params(b_p);
  powerset_domain_params::update_params(p_p);
  region_domain_params::update_params(r_p);
  zones_domain_params::update_params(z_p);
  oct_domain_params::update_params(o_p);
  fixed_tvpi_domain_params::update_params(tvpi_p);
}

static bool to_bool(const std::string &val) {
  if (val == "true") {
    return true;
  } else if (val == "false") {
    return false;
  } else {
    CRAB_ERROR("parameter value ", val, "cannot be converted to \"true\" or \"false\"");
  }
}

static int to_int(const std::string &val) {
  try {
    int res = std::stoi(val);
    return res;
  } catch(std::invalid_argument const& e) {
    CRAB_ERROR("parameter value ", val, " cannot be converted to an integer");
  } catch(std::out_of_range const& e) {
    CRAB_ERROR("parameter value ", val, " is out of range for an integer");
  } 
}

static unsigned to_uint(const std::string &val) {
  try {
    unsigned res =  std::stoul(val);
    return res;
  } catch(std::invalid_argument const& e) {
    CRAB_ERROR("parameter value ", val, " cannot be converted to an unsigned integer");
  } catch(std::out_of_range const& e) {
    CRAB_ERROR("parameter value ", val, " is out of range for an unsigned integer");
  } 
}

static std::vector<unsigned> to_list_of_unsigned(const std::string &val) {
  std::vector<unsigned> res;
  std::stringstream ss(val);
  while (ss.good()) {
    std::string substr;
    getline(ss, substr, ',');
    res.push_back(to_uint(substr));
  }
  return res;
}
  
void crab_domain_params::set_param(const std::string &param, const std::string &val) {
  if (param == "elina.use_tree_expressions") {
    elina_domain_params::m_use_tree_expressions = to_bool(val);
  } else if (param == "array_adaptive.is_smashable") {
    array_adaptive_domain_params::m_is_smashable = to_bool(val);
  } else if (param == "array_adaptive.smash_at_nonzero_offset") {
    array_adaptive_domain_params::m_smash_at_nonzero_offset = to_bool(val);
  } else if (param == "array_adaptive.max_smashable_cells") {
    array_adaptive_domain_params::m_max_smashable_cells = to_uint(val);
  } else if (param == "array_adaptive.max_array_size") {
    array_adaptive_domain_params::m_max_array_size = to_uint(val);
  } else if (param == "boxes.ldd_size") {
    boxes_domain_params::m_ldd_size = to_uint(val);
  } else if (param == "boxes.convexify_threshold") {
    boxes_domain_params::m_convexify_threshold = to_int(val);
  } else if (param == "boxes.dynamic_reordering") {
    boxes_domain_params::m_dynamic_reordering = to_bool(val);
  } else if (param == "powerset.exact_meet") {
    powerset_domain_params::m_exact_meet = to_bool(val);
  } else if (param == "powerset.max_disjuncts") {
    powerset_domain_params::m_max_disjuncts = to_uint(val);
  } else if (param == "region.allocation_sites") {
    region_domain_params::m_allocation_sites = to_bool(val);
  } else if (param == "region.deallocation") {
    region_domain_params::m_deallocation = to_bool(val);
  } else if (param == "region.tag_analysis") {
    region_domain_params::m_tag_analysis = to_bool(val);
  } else if (param == "region.is_dereferenceable") {
    region_domain_params::m_is_dereferenceable = to_bool(val);
  } else if (param == "region.skip_unknown_regions") {
    region_domain_params::m_skip_unknown_regions = to_bool(val);
  } else if (param == "zones.chrome_dijkstra") {
    zones_domain_params::m_chrome_dijkstra = to_bool(val);
  } else if (param == "zones.widen_restabilize") {
    zones_domain_params::m_widen_restabilize = to_bool(val);
  } else if (param == "zones.special_assign") {
    zones_domain_params::m_special_assign = to_bool(val);
  } else if (param == "zones.close_bounds_inline") {
    zones_domain_params::m_close_bounds_inline = to_bool(val);
  } else if (param == "oct.chrome_dijkstra") {
    oct_domain_params::m_chrome_dijkstra = to_bool(val);    
  } else if (param == "oct.widen_restabilize") {
    oct_domain_params::m_widen_restabilize = to_bool(val);    
  } else if (param == "oct.special_assign") {
    oct_domain_params::m_special_assign = to_bool(val);    
  } else if (param == "oct.close_bounds_inline") {
    oct_domain_params::m_close_bounds_inline = to_bool(val);    
  } else if (param == "fixed_tvpi.coefficients") {
    fixed_tvpi_domain_params::m_coefficients = to_list_of_unsigned(val);    
  } else {
    CRAB_WARN("Ignored unsupported parameter ",  param);
  }
}

void crab_domain_params::write(crab::crab_os &o) const {
  elina_domain_params::write(o);
  array_adaptive_domain_params::write(o);
  boxes_domain_params::write(o);
  powerset_domain_params::write(o);
  region_domain_params::write(o);
  zones_domain_params::write(o);
  oct_domain_params::write(o);
  fixed_tvpi_domain_params::write(o);
}
} // end namespace domains
} // end namespace crab
