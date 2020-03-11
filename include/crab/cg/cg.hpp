#pragma once

/*
   Build a call graph (CG)
*/

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/optional.hpp>

#include <crab/cfg/cfg.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>

#include <functional> // for wrapper_reference and hash
#include <memory>
#include <unordered_map>

namespace crab {
namespace cg {
// Wrapper for call graph nodes
template <typename CFG> class cg_node {
  typedef typename CFG::basic_block_t::callsite_t callsite_t;
  typedef typename CFG::fdecl_t fdecl_t;

  CFG m_cfg;
  int m_id;

public:
  typedef CFG cfg_t;
  typedef typename CFG::varname_t varname_t;

  cg_node() {} // needed for BGL

  cg_node(CFG cfg, int id) : m_cfg(cfg), m_id(id) {}

  CFG get_cfg() const { return m_cfg; }

  int index() const { return m_id; }

  std::string name() const {
    if (!m_cfg.has_func_decl()) {
      CRAB_ERROR("No function name found");
    }
    return m_cfg.get_func_decl().get_func_name();
  }

  bool operator==(const cg_node &o) const { return index() == o.index(); }

  bool operator!=(const cg_node &o) const { return !(*this == o); }

  size_t hash() const {
    std::hash<int> hasher;
    return hasher(m_id);
  }

  bool operator<(const cg_node &o) const { return index() < o.index(); }

  friend crab_os &operator<<(crab_os &o, cg_node n) {
    o << n.name();
    return o;
  }
};
} // end namespace cg
} // end namespace crab

/**  specialization of std::hash for callgraph nodes **/
namespace std {
template <typename CFG> struct hash<crab::cg::cg_node<CFG>> {
  using cg_node_t = crab::cg::cg_node<CFG>;
  size_t operator()(const cg_node_t &n) const { return n.hash(); }
};
} // end namespace std

namespace crab {
namespace cg {
// Class to build a call graph
// Important: this class assumes that all function calls have been
// resolved. This must be ensured by the client.
template <typename CFG> class call_graph {
  // Wrapper for call graph edges
  // BGL complains if we use std::pair<cg_node,cg_node>
  template <typename T> struct cg_edge {
    T m_s;
    T m_d;
    cg_edge() {}
    cg_edge(T s, T d) : m_s(s), m_d(d) {}
    T src() const { return m_s; }
    T dest() const { return m_d; }
    bool operator==(const cg_edge<T> &o) const {
      return (m_s == o.src() && m_d == o.dest());
    }
    bool operator!=(const cg_edge<T> &o) const { return !(*this == o); }
  };

  /// --- begin internal representation of the call graph
  struct vertex_t {
    cg_node<CFG> func;
  };
  typedef boost::adjacency_list<
      boost::setS, // disallow parallel edges
      boost::vecS, boost::bidirectionalS,
      boost::property<boost::vertex_color_t, boost::default_color_type,
                      vertex_t>>
      cg_t;
  typedef
      typename boost::graph_traits<cg_t>::vertex_descriptor vertex_descriptor_t;
  typedef typename boost::graph_traits<cg_t>::edge_descriptor edge_descriptor_t;
  typedef typename boost::graph_traits<cg_t>::vertex_iterator vertex_iterator;
  typedef
      typename boost::graph_traits<cg_t>::out_edge_iterator out_edge_iterator;
  typedef typename boost::graph_traits<cg_t>::in_edge_iterator in_edge_iterator;
  /// --- end internal representation of the call graph

  typedef typename CFG::varname_t varname_t;
  typedef typename CFG::number_t number_t;
  typedef crab::cfg::statement_visitor<number_t, varname_t> stmt_visitor_t;

  typedef std::unordered_map<std::size_t, vertex_descriptor_t> vertex_map_t;
  typedef std::unordered_map<typename stmt_visitor_t::callsite_t *,
                             cg_node<CFG>>
      callee_map_t;
  typedef std::unordered_map<cg_node<CFG>, vertex_descriptor_t>
      node_vertex_id_map_t;

  struct mk_edge_vis : public stmt_visitor_t {
    typedef typename stmt_visitor_t::bin_op_t bin_op_t;
    typedef typename stmt_visitor_t::assign_t assign_t;
    typedef typename stmt_visitor_t::assume_t assume_t;
    typedef typename stmt_visitor_t::havoc_t havoc_t;
    typedef typename stmt_visitor_t::unreach_t unreach_t;
    typedef typename stmt_visitor_t::select_t select_t;
    typedef typename stmt_visitor_t::callsite_t callsite_t;
    typedef typename CFG::fdecl_t fdecl_t;

    cg_t &m_cg;
    vertex_map_t &m_vertex_map;
    callee_map_t &m_callee_map;
    std::size_t m_from;

    mk_edge_vis(cg_t &cg, vertex_map_t &vertex_map, callee_map_t &callee_map,
                fdecl_t &from)
        : m_cg(cg), m_vertex_map(vertex_map), m_callee_map(callee_map),
          m_from(crab::cfg::cfg_hasher<CFG>::hash(from)) {}

    void visit(callsite_t &cs) {
      size_t to = crab::cfg::cfg_hasher<CFG>::hash(cs);
      auto it_from = m_vertex_map.find(m_from);
      auto it_to = m_vertex_map.find(to);

      CRAB_LOG("cg", crab::outs() << "Visiting call site " << cs << "\n";);

      if (it_from == m_vertex_map.end()) {
        CRAB_LOG("cg", crab::outs() << "Not found caller \n";);
        return;
      }

      if (it_to == m_vertex_map.end()) {
        CRAB_LOG("cg", crab::outs() << "Not found callee \n";);
        return;
      }

      // -- add edge in the call graph.
      auto res = add_edge(it_from->second, it_to->second, m_cg);
      if (res.second) {
        CRAB_LOG("cg", crab::outs() << "Added cg edge " << it_from->second
                                    << " --> " << it_to->second << "\n";);
      }

      // -- record the callee's cfg with the callsite
      m_callee_map.insert({&cs, m_cg[it_to->second].func});
    }
  };

  struct mk_node
      : public std::unary_function<vertex_descriptor_t, cg_node<CFG>> {
    cg_t *_cg;
    mk_node() : _cg(nullptr) {}
    mk_node(cg_t *cg) : _cg(cg) {}
    cg_node<CFG> &operator()(const vertex_descriptor_t &v) const {
      assert(_cg);
      return (*_cg)[v].func;
    }
  };

  struct mk_edge
      : public std::unary_function<edge_descriptor_t, cg_edge<cg_node<CFG>>> {
    cg_t *_cg;
    mk_edge() : _cg(nullptr) {}
    mk_edge(cg_t *cg) : _cg(cg) {}
    cg_edge<cg_node<CFG>> operator()(const edge_descriptor_t &e) const {
      assert(_cg);
      cg_node<CFG> &s = (*_cg)[boost::source(e, (*_cg))].func;
      cg_node<CFG> &t = (*_cg)[boost::target(e, (*_cg))].func;
      return cg_edge<cg_node<CFG>>(s, t);
    }
  };

public:
  typedef cg_node<CFG> node_t;
  typedef cg_edge<node_t> edge_t;
  typedef boost::transform_iterator<mk_node, vertex_iterator> node_iterator;
  typedef boost::transform_iterator<mk_edge, in_edge_iterator> pred_iterator;
  typedef boost::transform_iterator<mk_edge, out_edge_iterator> succ_iterator;

  typedef typename node_t::cfg_t cfg_t;
  typedef typename stmt_visitor_t::callsite_t callsite_t;
  typedef typename cfg_t::fdecl_t fdecl_t;

private:
  // call graph
  std::shared_ptr<cg_t> m_cg;
  // map from callsite to callee's CFG
  callee_map_t m_callee_map;

  //// INTERNAL STATE USED DURING CG CONSTRUCTION

  // map hashed values to internal BGL vertex descriptor
  vertex_map_t m_vertex_map;
  // map cg_node to internal BGL vertex descriptor
  node_vertex_id_map_t m_node_vertex_id_map;
  // counter to generate unique ids
  int m_id;

  vertex_descriptor_t get_vertex(const node_t &n) const {
    auto It = m_node_vertex_id_map.find(n);
    if (It != m_node_vertex_id_map.end())
      return It->second;
    CRAB_ERROR("Call graph could not find node");
  }

  template <typename CFGIt> void build_call_graph(CFGIt I, CFGIt E) {
    crab::ScopedCrabStats __st__("call_graph");

    // --- add vertices in the call graph
    for (auto cfg : boost::make_iterator_range(I, E)) {
      if (!cfg.has_func_decl()) {
        CRAB_ERROR("Could not compute call graph: function info is missing.");
      }

      auto decl = cfg.get_func_decl();
      size_t k = crab::cfg::cfg_hasher<CFG>::hash(decl);
      vertex_descriptor_t v = add_vertex(*m_cg);
      m_vertex_map.insert({k, v});
      node_t f(cfg, m_id++);
      m_node_vertex_id_map.insert({f, v});
      (*m_cg)[v].func = f;

      CRAB_LOG("cg", crab::outs() << "Added call graph node " << decl
                                  << "--- id=" << v << "\n";);
    }

    // --- add edges in the call graph
    for (auto cfg : boost::make_iterator_range(I, E)) {
      assert(cfg.has_func_decl());
      auto decl = cfg.get_func_decl();
      for (auto const &bb :
           boost::make_iterator_range(cfg.begin(), cfg.end())) {
        mk_edge_vis vis(*m_cg, m_vertex_map, m_callee_map, decl);
        for (auto it = bb.begin(); it != bb.end(); ++it) {
          it->accept(&vis);
        }
      }
    }
  }

public:
  call_graph(std::vector<CFG> &cfgs) : m_cg(new cg_t()), m_id(0) {
    build_call_graph(cfgs.begin(), cfgs.end());
  }

  template <typename CFGIt>
  call_graph(CFGIt I, CFGIt E) : m_cg(new cg_t()), m_id(0) {
    build_call_graph(I, E);
  }

  call_graph(const call_graph<CFG> &o) = delete;

  call_graph<CFG> &operator=(const call_graph<CFG> &o) = delete;

  // Check type consistency between function declaration and callsite.
  void type_check() const {
    for (auto const &kv : m_callee_map) {
      CFG callee_cfg = kv.second.get_cfg();
      if (!callee_cfg.has_func_decl()) {
        CRAB_ERROR("CFG without function declaration");
      }

      /// Crab only needs a CFG to have an exit block when performing
      /// inter-procedural or backward analysis. Thus, a CFG without
      /// exit block is still considered well formed. We delegate to
      /// the corresponding analysis to deal with a CFG without an
      /// exit block.

      // if (!callee_cfg.has_exit()) {
      //   CRAB_ERROR("CFG has no exit");
      // }

      const callsite_t &cs = *kv.first;
      fdecl_t fdecl = callee_cfg.get_func_decl();

      if (fdecl.get_num_inputs() != cs.get_num_args()) {
        crab::errs() << "Callsite: " << cs << "\n";
        crab::errs() << "Function declaration: " << fdecl << "\n";
        crab::errs() << callee_cfg << "\n";
        CRAB_ERROR(
            "Mismatch between number of callsite and function parameters");
      }
      if (fdecl.get_num_outputs() != cs.get_lhs().size()) {
        crab::errs() << "Callsite: " << cs << "\n";
        crab::errs() << "Function declaration: " << fdecl << "\n";
        CRAB_ERROR(
            "Mismatch between number of callsite and function return values");
      }
      for (unsigned i = 0; i < cs.get_num_args(); i++) {
        if (fdecl.get_input_type(i) != cs.get_arg_type(i)) {
          crab::errs() << "Callsite: " << cs << "\n";
          crab::errs() << "Function declaration: " << fdecl << "\n";
          CRAB_ERROR(
              "Mismatch between type of callsite and function parameter");
        }
      }
      for (unsigned i = 0; i < cs.get_lhs().size(); i++) {
        if (fdecl.get_output_type(i) != cs.get_lhs()[i].get_type()) {
          crab::errs() << "Callsite: " << cs << "\n";
          crab::errs() << "Function declaration: " << fdecl << "\n";
          CRAB_ERROR(
              "Mismatch between type of callsite and function return value");
        }
      }
    }

    // check each callsite has a corresponding function
    for (auto &cg_node : boost::make_iterator_range(nodes())) {
      CFG cfg = cg_node.get_cfg();
      for (auto &bb : boost::make_iterator_range(cfg.begin(), cfg.end())) {
        for (auto &s : boost::make_iterator_range(bb.begin(), bb.end())) {
          if (s.is_callsite()) {
            auto cs = static_cast<callsite_t *>(&s);
            if (!has_callee(*cs)) {
              CRAB_ERROR("Function not found for callsite ", *cs);
            }
          }
        }
      }
    }
  }

  node_t entry() const {
    // Any node without incoming edges should be considered an entry
    // point. In addition, if all nodes have some incoming edges then
    // it can be the case that the analysis should start from some SCC
    // with multiple nodes. In that case, all SCC's components should
    // be considered as entry point.

    // FIXME: for now, we assume that the call graph has exactly one
    // node without incoming edges. For libraries, we can transform
    // the program in such way that we create a node that calls all
    // library's entry points.

    std::vector<node_t> entries;
    for (node_iterator it = nodes().first, et = nodes().second; it != et;
         ++it) {
      if (num_preds(*it) == 0) {
        entries.push_back(*it);
      }
    }
    size_t num_entries = entries.size();
    if (num_entries == 0) {
      CRAB_ERROR("cannot find entry point of the call graph");
    } else if (num_entries > 1) {
      for (unsigned i = 0, e = entries.size(); i < e; i++) {
        if (entries[i].name() == "main") {
          return entries[i];
        }
      }
      CRAB_ERROR("do not support call graphs with multiple entry points");
    } else {
      return entries[0];
    }
  }

  std::vector<node_t> entries() const {
    // TODO: any node without incoming edges should be considered an
    // entry point. In addition, if all nodes have some incoming edges
    // then it can be the case that the analysis should start from
    // some SCC with multiple nodes. In that case, all SCC's
    // components should be considered as entry point.

    std::vector<node_t> out;
    for (node_iterator it = nodes().first, et = nodes().second; it != et;
         ++it) {
      if (num_preds(*it) == 0) {
        out.push_back(*it);
      }
    }
    return out;
  }

  bool has_callee(callsite_t &cs) const {
    return m_callee_map.find(&cs) != m_callee_map.end();
  }

  node_t get_callee(callsite_t &cs) const {
    auto it = m_callee_map.find(&cs);
    assert(it != m_callee_map.end());
    return it->second;
  }

  std::pair<node_iterator, node_iterator> nodes() const {
    auto p = boost::vertices(*m_cg);
    return std::make_pair(make_transform_iterator(p.first, mk_node(&*m_cg)),
                          make_transform_iterator(p.second, mk_node(&*m_cg)));
  }

  std::pair<succ_iterator, succ_iterator> succs(const node_t &n) const {
    vertex_descriptor_t v = get_vertex(n);
    auto p = boost::out_edges(v, *m_cg);
    return std::make_pair(make_transform_iterator(p.first, mk_edge(&*m_cg)),
                          make_transform_iterator(p.second, mk_edge(&*m_cg)));
  }

  std::pair<pred_iterator, pred_iterator> preds(const node_t &n) const {
    vertex_descriptor_t v = get_vertex(n);
    auto p = boost::in_edges(v, *m_cg);
    return std::make_pair(make_transform_iterator(p.first, mk_edge(&*m_cg)),
                          make_transform_iterator(p.second, mk_edge(&*m_cg)));
  }

  std::size_t num_nodes() const { return boost::num_vertices(*m_cg); }

  std::size_t num_succs(const node_t &n) const {
    vertex_descriptor_t v = get_vertex(n);
    return boost::out_degree(v, *m_cg);
  }

  std::size_t num_preds(const node_t &n) const {
    vertex_descriptor_t v = get_vertex(n);
    return boost::in_degree(v, *m_cg);
  }

  void write(crab_os &o) const {
    for (auto f : boost::make_iterator_range(nodes())) {
      if (num_succs(f) > 0) {
        for (auto e : boost::make_iterator_range(succs(f))) {
          o << e.src() << "--> " << e.dest() << "\n";
        }
      }
    }
  }

}; // end class call_graph<CFG>

template <typename CFG>
inline crab_os &operator<<(crab_os &o, const call_graph<CFG> &cg) {
  cg.write(o);
  return o;
}

// A lightweight object that wraps a reference to a call_graph into a
// copyable, assignable object.
template <class CG> class call_graph_ref {
public:
  typedef typename CG::node_t node_t;
  typedef typename node_t::cfg_t cfg_t;
  typedef typename CG::edge_t edge_t;
  typedef typename CG::node_iterator node_iterator;
  typedef typename CG::pred_iterator pred_iterator;
  typedef typename CG::succ_iterator succ_iterator;
  typedef typename CG::callsite_t callsite_t;

private:
  boost::optional<std::reference_wrapper<CG>> _ref;

public:
  call_graph_ref(CG &cg) : _ref(std::reference_wrapper<CG>(cg)) {}

  const CG &get() const {
    assert(_ref);
    return *_ref;
  }

  CG &get() {
    assert(_ref);
    return *_ref;
  }

  void type_check() const {
    assert(_ref);
    return (*_ref).get().type_check();
  }

  node_t entry() const {
    assert(_ref);
    return (*_ref).get().entry();
  }

  std::vector<node_t> entries() const {
    assert(_ref);
    return (*_ref).get().entries();
  }

  bool has_callee(callsite_t &cs) const {
    assert(_ref);
    return (*_ref).get().has_callee(cs);
  }

  node_t get_callee(callsite_t &cs) const {
    assert(_ref);
    return (*_ref).get().get_callee(cs);
  }

  std::pair<node_iterator, node_iterator> nodes() const {
    assert(_ref);
    return (*_ref).get().nodes();
  }

  std::pair<succ_iterator, succ_iterator> succs(const node_t &n) const {
    assert(_ref);
    return (*_ref).get().succs(n);
  }

  std::pair<pred_iterator, pred_iterator> preds(const node_t &n) const {
    assert(_ref);
    return (*_ref).get().preds(n);
  }

  std::size_t num_nodes() const {
    assert(_ref);
    return (*_ref).get().num_nodes();
  }

  std::size_t num_succs(const node_t &n) const {
    assert(_ref);
    return (*_ref).get().num_succs(n);
  }

  std::size_t num_preds(const node_t &n) const {
    assert(_ref);
    return (*_ref).get().num_preds(n);
  }

  void write(crab_os &o) const {
    assert(_ref);
    (*_ref).get().write(o);
  }
};

template <typename CG>
inline crab_os &operator<<(crab_os &o, const call_graph_ref<CG> &cg) {
  cg.write(o);
  return o;
}

} // end namespace cg
} // end namespace crab
