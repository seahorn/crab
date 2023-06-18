#pragma once

#include <crab/cfg/cfg.hpp>
#include <crab/checkers/base_property.hpp>
#include <crab/support/os.hpp>

#include <boost/algorithm/string/replace.hpp>
#include <boost/range/iterator_range.hpp>

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

/* Utilities to print a CFG into dot format */

namespace crab {
namespace cfg {

namespace dot_impl {
template <typename CFG>
void cfg_body_to_dot(const CFG &cfg, const checker::checks_db &db,
                     crab_string_os &os) {
  using basic_block_t = typename CFG::basic_block_t;
  using statement_t = typename CFG::statement_t;
  auto rewrite = [](std::string &str,
                    const std::map<std::pair<unsigned, std::string>,
                                   std::string> &rewrite_map) {
    for (auto &kv : rewrite_map) {
      boost::replace_all(str, kv.first.second, kv.second);
    }
  };

  // The order in which the rewriting happens matters
  std::map<std::pair<unsigned, std::string>, std::string> op_rewrite_map = {
      {{1, "<="}, "&le;"}, {{2, ">="}, "&ge;"}, {{3, "!="}, "&ne;"},
      {{4, "<"}, "&lt;"},  {{5, ">"}, "&gt;"},
  };

  auto print_checks = [&db](const statement_t &stmt, crab_string_os &os) {
    const crab::cfg::debug_info &di = stmt.get_debug_info();
    if (db.has_checks(di)) {
      os << "\t//// File:" << di.get_file() << " line:" << di.get_line()
         << " col:" << di.get_column() << " "
         << "Result: ";
      auto const &checks = db.get_checks(di);
      unsigned safe = 0;
      unsigned warning = 0;
      unsigned error = 0;
      for (unsigned i = 0, num_checks = checks.size(); i < num_checks; ++i) {
        switch (checks[i]) {
        case checker::check_kind::CRAB_SAFE:
        case checker::check_kind::CRAB_UNREACH:
          safe++;
          break;
        case checker::check_kind::CRAB_ERR:
          error++;
          break;
        default:
          warning++;
          break;
        }
      }
      if (error == 0 && warning == 0) {
        os << " OK";
      } else {
        os << " FAIL -- ";
        if (safe > 0)
          os << "num of safe=" << safe << " ";
        if (error > 0)
          os << "num of errors=" << error << " ";
        if (warning > 0)
          os << "num of warnings=" << warning << " ";
      }
      os << "\\l";
    }
  };

  auto statement_to_str = [&op_rewrite_map, &rewrite, &print_checks]
    (const statement_t &stmt) {
    crab_string_os os_s;
    print_checks(stmt, os_s);
    os_s << stmt;
    std::string ss = os_s.str();
    rewrite(ss, op_rewrite_map);
    return ss;
  };

  auto print_node = [&statement_to_str](const basic_block_t *node,
                                        crab_string_os &os) {
    os << "\tNode" << node << " ";
    os << "[shape=record,label=\"{" << node->label() << ":";
    for (auto const &stmt : *node) {
      os << "\\l " << statement_to_str(stmt);
      ;
    }
    os << "\\l}\"];\n";
  };

  auto print_edge = [](const basic_block_t *src, const basic_block_t *dest,
                       crab_string_os &os) {
    os << "\tNode" << src << " -> "
       << "Node" << dest << ";\n";
  };

  std::vector<const basic_block_t *> worklist;
  std::set<const basic_block_t *> visited;

  const basic_block_t *entry = &(cfg.get_node(cfg.entry()));
  worklist.push_back(entry);
  visited.insert(entry);
  while (!worklist.empty()) {
    const basic_block_t *curr_node = worklist.back();
    worklist.pop_back();
    print_node(curr_node, os);
    for (auto succ : boost::make_iterator_range(curr_node->next_blocks())) {
      const basic_block_t *succ_node = &(cfg.get_node(succ));
      print_edge(curr_node, succ_node, os);
      if (visited.insert(succ_node).second) {
        worklist.push_back(succ_node);
      }
    }
  }
}
} // end namespace dot_impl

template <typename CFG> void cfg_to_dot(const CFG &cfg) {
  if (cfg.has_func_decl()) {
    crab_string_os os;
    auto const &decl = cfg.get_func_decl();
    os << "digraph \"CFG for \'" << decl.get_func_name() << "\' function\" {\n";
    os << "\tlabel=\"CFG for \'" << decl.get_func_name() << "\' function\";\n";
    checker::checks_db db;
    dot_impl::cfg_body_to_dot(cfg, db, os);
    os << "}\n";
    std::string outfilename(decl.get_func_name() + ".crab.dot");
    std::ofstream outfile;
    outfile.open(outfilename);
    outfile << os.str();
    outfile.close();
    crab::outs() << "Writing \'" << outfilename << "\'...\n";
  }
}

template <typename CFG, typename Dom>
void cfg_to_dot(const CFG &cfg,
                std::function<boost::optional<Dom>(
                    const typename CFG::basic_block_label_t &)>
                    pre_fn,
                std::function<boost::optional<Dom>(
                    const typename CFG::basic_block_label_t &)>
                    post_fn,
                const checker::checks_db &db) {
  using basic_block_t = typename CFG::basic_block_t;
  auto print_invariants = [](const basic_block_t *node,
                             boost::optional<Dom> pre,
                             boost::optional<Dom> post, crab_string_os &os) {
    if (pre) {
      os << "\tNodePreInv" << node << " ";
      os << "[style=\"filled\",fillcolor=cornsilk,label=\"" << *pre << "\"];\n";
      os << "\tNodePreInv" << node << " -> "
         << "Node" << node << " [arrowhead= diamond];\n";
    }
    if (post) {
      os << "\tNodePostInv" << node << " ";
      os << "[style=\"filled\",fillcolor=azure,label=\"" << *post << "\"];\n";
      os << "\tNode" << node << " -> "
         << "NodePostInv" << node << " [arrowhead= diamond];\n";
    }
  };

  if (cfg.has_func_decl()) {
    crab_string_os os;
    auto const &decl = cfg.get_func_decl();
    os << "digraph \"CFG for \'" << decl.get_func_name() << "\' function\" {\n";
    os << "\tlabel=\"CFG for \'" << decl.get_func_name() << "\' function\";\n";
    dot_impl::cfg_body_to_dot(cfg, db, os);
    // Print invariants
    for (auto it = cfg.begin(), et = cfg.end(); it != et; ++it) {
      const basic_block_t &node = *it;
      print_invariants(&node, pre_fn(node.label()), post_fn(node.label()), os);
    }
    os << "}\n";
    std::string outfilename(decl.get_func_name() + ".crab.dot");
    std::ofstream outfile;
    outfile.open(outfilename);
    outfile << os.str();
    outfile.close();
    crab::outs() << "Writing \'" << outfilename << "\'...\n";
  }
}

} // end namespace cfg
} // end namespace crab
