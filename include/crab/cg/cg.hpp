#pragma once

/*
   Build a call graph (CG)
*/

#include <crab/cfg/cfg.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <boost/range/iterator_range.hpp>

#include <cassert>
#include <functional> // for wrapper_reference and hash
#include <string>
#include <unordered_map>
#include <vector>

namespace crab {
namespace cg {
// Wrapper for call graph nodes
template <typename CFG> class cg_node {
  using callsite_t = typename CFG::basic_block_t::callsite_t;
  using fdecl_t = typename CFG::fdecl_t;

  CFG m_cfg;
  int m_id = -1;

public:
  using cfg_t = CFG;
  using varname_t = typename CFG::varname_t;

  cg_node() {} // needed for BGL

  cg_node(CFG cfg, int id) : m_cfg(cfg), m_id(id) {}

  CFG get_cfg() const { return m_cfg; }

  int index() const { return m_id; }

  const std::string &name() const {
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
  // Wrapper for call graph edges (a source/destination node pair)
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

  using varname_t = typename CFG::varname_t;
  using number_t = typename CFG::number_t;
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using stmt_visitor_t =
      crab::cfg::statement_visitor<basic_block_label_t, number_t, varname_t>;
  using callsite_or_fdecl_t = crab::cfg::callsite_or_fdecl<CFG>;
  // map a callsite or function declaration to a dense vertex id
  using vertex_map_t = crab::cfg::callsite_or_fdecl_map<CFG, std::size_t>;
  using callee_map_t =
      std::unordered_map<const typename stmt_visitor_t::callsite_t *,
                         cg_node<CFG>>;

public:
  using node_t = cg_node<CFG>;
  using edge_t = cg_edge<node_t>;
  // Nodes and (per-node) edges are stored contiguously, so iterating them
  // is just iterating the underlying vectors.
  using node_iterator = typename std::vector<node_t>::const_iterator;
  using succ_iterator = typename std::vector<edge_t>::const_iterator;
  using pred_iterator = typename std::vector<edge_t>::const_iterator;

  using cfg_t = typename node_t::cfg_t;
  using callsite_t = typename stmt_visitor_t::callsite_t;
  using fdecl_t = typename cfg_t::fdecl_t;

private:
  // Statement visitor that adds one call graph edge per resolved callsite.
  struct mk_edge_vis : public stmt_visitor_t {
    using bin_op_t = typename stmt_visitor_t::bin_op_t;
    using assign_t = typename stmt_visitor_t::assign_t;
    using assume_t = typename stmt_visitor_t::assume_t;
    using havoc_t = typename stmt_visitor_t::havoc_t;
    using unreach_t = typename stmt_visitor_t::unreach_t;
    using select_t = typename stmt_visitor_t::select_t;
    using callsite_t = typename stmt_visitor_t::callsite_t;
    using fdecl_t = typename CFG::fdecl_t;

    call_graph &m_parent;
    const fdecl_t &m_from;

    mk_edge_vis(call_graph &parent, const fdecl_t &from)
        : m_parent(parent), m_from(from) {}

    virtual void visit(callsite_t &cs) override {
      auto &vertex_map = m_parent.m_vertex_map;
      auto it_from = vertex_map.find(&m_from);
      auto it_to = vertex_map.find(&cs);

      CRAB_LOG("cg", crab::outs() << "Visiting call site " << cs << "\n";);

      if (it_from == vertex_map.end()) {
        CRAB_LOG("cg", crab::outs() << "Not found caller \n";);
        return;
      }

      if (it_to == vertex_map.end()) {
        CRAB_LOG("cg", crab::outs() << "Not found callee \n";);
        return;
      }

      // -- add edge in the call graph.
      m_parent.add_edge(it_from->second, it_to->second);

      // -- record the callee's cfg with the callsite
      m_parent.m_callee_map.insert({&cs, m_parent.m_nodes[it_to->second]});
    }
  };

  // --- internal representation of the call graph
  // Vertices are dense ids 0..N-1; a node's id (cg_node::index()) is its
  // position in these vectors.
  std::vector<node_t> m_nodes;              // all call graph nodes
  std::vector<std::vector<edge_t>> m_succs; // outgoing edges per node
  std::vector<std::vector<edge_t>> m_preds; // incoming edges per node

  // map from callsite to callee's CFG
  callee_map_t m_callee_map;

  // map a function declaration to its vertex id. Only used while building
  // the graph; cleared afterwards.
  vertex_map_t m_vertex_map;
  // counter to generate unique ids (kept in sync with m_nodes.size())
  int m_id;

  std::size_t get_vertex(const node_t &n) const {
    std::size_t id = static_cast<std::size_t>(n.index());
    if (id >= m_nodes.size() || !(m_nodes[id] == n)) {
      CRAB_ERROR("Call graph could not find node");
    }
    return id;
  }

  // Add a caller -> callee edge, disallowing parallel edges.
  void add_edge(std::size_t from, std::size_t to) {
    for (const edge_t &e : m_succs[from]) {
      if (e.dest().index() == static_cast<int>(to)) {
        return;
      }
    }
    edge_t e(m_nodes[from], m_nodes[to]);
    m_succs[from].push_back(e);
    m_preds[to].push_back(e);
    CRAB_LOG("cg", crab::outs()
                       << "Added cg edge " << from << " --> " << to << "\n";);
  }

  template <typename CFGIt> void build_call_graph(CFGIt I, CFGIt E) {
    //crab::ScopedCrabStats __st__("call_graph", false);

    // --- add vertices in the call graph
    for (auto cfg : boost::make_iterator_range(I, E)) {
      if (!cfg.has_func_decl()) {
        CRAB_ERROR("Could not compute call graph: function info is missing.");
      }

      auto const &decl = cfg.get_func_decl();
      std::size_t v = m_nodes.size();
      m_nodes.emplace_back(cfg, m_id++);
      m_succs.emplace_back();
      m_preds.emplace_back();
      m_vertex_map.insert({callsite_or_fdecl_t(&decl), v});

      CRAB_LOG("cg", crab::outs() << "Added call graph node " << decl
                                  << "--- id=" << v << "\n";);
    }

    // --- add edges in the call graph
    for (auto cfg : boost::make_iterator_range(I, E)) {
      assert(cfg.has_func_decl());
      auto const &decl = cfg.get_func_decl();
      for (auto const &bb :
           boost::make_iterator_range(cfg.begin(), cfg.end())) {
        mk_edge_vis vis(*this, decl);
        for (auto it = bb.begin(); it != bb.end(); ++it) {
          it->accept(&vis);
        }
      }
    }

    // m_vertex_map is only needed while building the graph; runtime vertex
    // lookups go through cg_node::index() (see get_vertex).
    m_vertex_map.clear();
  }

public:
  call_graph(const std::vector<CFG> &cfgs) : m_id(0) {
    build_call_graph(cfgs.begin(), cfgs.end());
  }

  template <typename CFGIt>
  call_graph(CFGIt I, CFGIt E) : m_id(0) {
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
      const fdecl_t &fdecl = callee_cfg.get_func_decl();

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
    // FIXME: for now, we assume that the call graph has exactly one
    // node without incoming edges. For libraries, we can transform
    // the program in such way that we create a node that calls all
    // library's entry points.

    std::vector<node_t> es = entries();
    size_t num_entries = es.size();
    if (num_entries == 0) {
      CRAB_ERROR("cannot find entry point of the call graph");
    } else if (num_entries == 1) {
      return es[0];
    } else {
      for (unsigned i = 0, e = es.size(); i < e; i++) {
        if (es[i].name() == "main") {
          return es[i];
        }
      }
      CRAB_ERROR("do not support call graphs with multiple entry points");
    }
  }

  std::vector<node_t> entries() const {
    // Any node without incoming edges is considered an entry point.
    //
    // TODO: if all nodes have some incoming edges then it can be the
    // case that the analysis should start from some SCC with multiple
    // nodes. In that case, all SCC's components should be considered as
    // entry point.

    std::vector<node_t> out;
    for (node_iterator it = nodes().first, et = nodes().second; it != et;
         ++it) {
      if (num_preds(*it) == 0) {
        out.push_back(*it);
      }
    }
    return out;
  }

  bool has_callee(const callsite_t &cs) const {
    return m_callee_map.find(&cs) != m_callee_map.end();
  }

  node_t get_callee(const callsite_t &cs) const {
    auto it = m_callee_map.find(&cs);
    if (it == m_callee_map.end()) {
      CRAB_ERROR("Call graph could not find callee for callsite");
    }
    return it->second;
  }

  std::pair<node_iterator, node_iterator> nodes() const {
    return std::make_pair(m_nodes.begin(), m_nodes.end());
  }

  std::pair<succ_iterator, succ_iterator> succs(const node_t &n) const {
    std::size_t v = get_vertex(n);
    return std::make_pair(m_succs[v].begin(), m_succs[v].end());
  }

  std::pair<pred_iterator, pred_iterator> preds(const node_t &n) const {
    std::size_t v = get_vertex(n);
    return std::make_pair(m_preds[v].begin(), m_preds[v].end());
  }

  std::size_t num_nodes() const { return m_nodes.size(); }

  std::size_t num_succs(const node_t &n) const {
    return m_succs[get_vertex(n)].size();
  }

  std::size_t num_preds(const node_t &n) const {
    return m_preds[get_vertex(n)].size();
  }

  void write(crab_os &o) const {
    for (auto f : boost::make_iterator_range(nodes())) {
      for (auto e : boost::make_iterator_range(succs(f))) {
        o << e.src() << "--> " << e.dest() << "\n";
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
  using node_t = typename CG::node_t;
  using cfg_t = typename node_t::cfg_t;
  using edge_t = typename CG::edge_t;
  using node_iterator = typename CG::node_iterator;
  using pred_iterator = typename CG::pred_iterator;
  using succ_iterator = typename CG::succ_iterator;
  using callsite_t = typename CG::callsite_t;

private:
  // The reference is always valid: the constructor requires a CG&.
  std::reference_wrapper<CG> _ref;

public:
  call_graph_ref(CG &cg) : _ref(cg) {}

  const CG &get() const { return _ref.get(); }

  CG &get() { return _ref.get(); }

  void type_check() const { return _ref.get().type_check(); }

  node_t entry() const { return _ref.get().entry(); }

  std::vector<node_t> entries() const { return _ref.get().entries(); }

  bool has_callee(const callsite_t &cs) const {
    return _ref.get().has_callee(cs);
  }

  node_t get_callee(const callsite_t &cs) const {
    return _ref.get().get_callee(cs);
  }

  std::pair<node_iterator, node_iterator> nodes() const {
    return _ref.get().nodes();
  }

  std::pair<succ_iterator, succ_iterator> succs(const node_t &n) const {
    return _ref.get().succs(n);
  }

  std::pair<pred_iterator, pred_iterator> preds(const node_t &n) const {
    return _ref.get().preds(n);
  }

  std::size_t num_nodes() const { return _ref.get().num_nodes(); }

  std::size_t num_succs(const node_t &n) const {
    return _ref.get().num_succs(n);
  }

  std::size_t num_preds(const node_t &n) const {
    return _ref.get().num_preds(n);
  }

  void write(crab_os &o) const { _ref.get().write(o); }
};

template <typename CG>
inline crab_os &operator<<(crab_os &o, const call_graph_ref<CG> &cg) {
  cg.write(o);
  return o;
}

} // end namespace cg
} // end namespace crab
