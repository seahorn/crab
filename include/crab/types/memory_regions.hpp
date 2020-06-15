#pragma once

#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>
#include <crab/types/indexable.hpp> 

#include <functional>

namespace crab {

class memory_region: public indexable {
public:
  
  // A region can contain data with one of these types:
  enum class type_t {BOOL,
		     INT,
		     REAL,
		     REF,
		     ANY};

  // Create an untype region
  static memory_region make_memory_region(ikos::index_t id) {
    return memory_region(id, type_t::ANY);
  }
  
  // Create a region that contains only booleans
  static memory_region make_bool_memory_region(ikos::index_t id) {
    return memory_region(id, type_t::BOOL);
  }

  // Create a region that contains only integers
  static memory_region make_int_memory_region(ikos::index_t id, unsigned bitwidth) {
    return memory_region(id, type_t::INT, bitwidth);
  }

  // Create a region that contains only reals
  static memory_region make_real_memory_region(ikos::index_t id) {
    return memory_region(id, type_t::REAL);
  }

  // Create a region that contains only references
  static memory_region make_ref_memory_region(ikos::index_t id) {
    return memory_region(id, type_t::REF);
  }
  
  virtual ikos::index_t index() const override { return m_id;}
  type_t get_type() const { return m_type;}
  unsigned get_bitwidth() const { return m_bitwidth;}
  
  bool operator<(const memory_region &o) const { return m_id < o.m_id;}
  bool operator==(const memory_region &o) const { return m_id == o.m_id;}
  bool operator!=(const memory_region &o) const { return m_id != o.m_id;}    
  size_t hash() const {
    // casting to size_t may overflow but it shouldn't affect
    // correctness
    return std::hash<size_t>{}(static_cast<size_t>(m_id));
  }
  
  void write(crab::crab_os &o) const {
    o << "region_" << m_id;
    CRAB_LOG("crab-print-types",
	     o << ":";
	     switch(m_type) {
	     case type_t::BOOL: o << "bool"; break;
	     case type_t::INT:  o << "int" << ":" << m_bitwidth; break;
	     case type_t::REAL: o << "real"; break;
	     case type_t::REF: o << "ref"; break;      
	     default: o << "untyped";
	     });
  }

  friend crab::crab_os &operator<<(crab::crab_os &o, const memory_region &m) {
    m.write(o);
    return o;
  }
  
private:
  // identifier for the memory region
  ikos::index_t m_id;
  // type for the memory region
  type_t m_type;
  // bitwidth if m_type is INT, otherwise 0
  unsigned m_bitwidth; 

  memory_region(ikos::index_t id, type_t type)
    : m_id(id), m_type(type), m_bitwidth(0) {}
  
  memory_region(ikos::index_t id, type_t type, unsigned bitwidth)
    : m_id(id), m_type(type), m_bitwidth(bitwidth) {}
  
};  
} // end namespace crab

/** specialization for std::hash for memory_region **/
namespace std {
template<>
struct hash<crab::memory_region> {
  size_t operator()(const crab::memory_region &mem) const { return mem.hash(); }
};
} // end namespace std


