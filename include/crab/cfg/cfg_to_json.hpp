#pragma once

/**
 * JSON serialization for Crab CFGs.
 *
 * The JSON counterpart of cfg_to_dot.hpp: where that produces something to look
 * at, this produces something to compute with. It emits the CFG *as Crab holds
 * it in memory* -- after any front-end lowering and after `simplify()` -- which
 * is the CFG the analyses actually see, and therefore the one whose labels
 * match the invariants exported by an abstract domain.
 *
 * Shape
 * -----
 *   {"name": "foo",
 *    "entry": "start",
 *    "exit": "___exit",                 // null if the CFG has no exit block
 *    "declaration": {"inputs": [var, ...], "outputs": [var, ...]},   // or null
 *    "blocks": [{"label": "start",
 *                "stmts": [stmt, ...],
 *                "succs": ["loop"]},
 *               ...]}
 *
 * Every block of the CFG is emitted, not only those reachable from the entry.
 * Blocks are sorted by label so that the output is stable: Crab stores them in
 * an unordered_map, whose iteration order is neither insertion order nor
 * portable across standard libraries, which would make the JSON undiffable.
 *
 * Statements
 * ----------
 * Every statement is an object with a "stmt" discriminator. All statement kinds
 * are exported, including the region and reference statements; a statement is
 * never silently skipped, because dropping one would produce a JSON document
 * describing a *different program* than the one Crab analyzed, which is the
 * most dangerous failure mode available here.
 *
 * The one deliberate exception is the `value_partition_start`/
 * `value_partition_end` intrinsics, which are directives to the value
 * partitioning domain with no effect on the concrete semantics; they are
 * omitted. Any other intrinsic raises CRAB_ERROR.
 *
 * Naming and order are uniform across every statement kind:
 *
 *   - "lhs" is always the variable the statement writes, "rhs" the one it
 *     reads, "left"/"right" when it reads two, "cond" for a condition. A
 *     qualifier is a suffix: "lhs_region", "rhs_width", "left_region".
 *   - Operands are printed sources first, destinations last, with "loc" (when
 *     present) always closing the object.
 *   - In a region/reference pair the region precedes the reference that lives
 *     in it, and both take the prefix of the side they belong to.
 *
 * So int_to_ref and ref_to_int read as exact mirrors of each other, the
 * reference operand carrying a "_ref" suffix on whichever side is written.
 */

#include <crab/cfg/cfg.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/json.hpp>
#include <crab/types/linear_constraints_to_json.hpp>
#include <crab/types/reference_constraints_to_json.hpp>
#include <crab/types/variable_to_json.hpp>

#include <boost/optional.hpp>

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace crab {
namespace cfg {

namespace json_impl {

template <typename CFG>
class statement_to_json_visitor
    : public const_statement_visitor<typename CFG::basic_block_label_t,
                                     typename CFG::number_t,
                                     typename CFG::varname_t> {

  using basic_block_label_t = typename CFG::basic_block_label_t;
  using number_t = typename CFG::number_t;
  using varname_t = typename CFG::varname_t;
  using const_statement_visitor_t =
      const_statement_visitor<basic_block_label_t, number_t, varname_t>;

  using bin_op_t = typename const_statement_visitor_t::bin_op_t;
  using assign_t = typename const_statement_visitor_t::assign_t;
  using assume_t = typename const_statement_visitor_t::assume_t;
  using select_t = typename const_statement_visitor_t::select_t;
  using assert_t = typename const_statement_visitor_t::assert_t;
  using int_cast_t = typename const_statement_visitor_t::int_cast_t;
  using unreach_t = typename const_statement_visitor_t::unreach_t;
  using havoc_t = typename const_statement_visitor_t::havoc_t;
  using callsite_t = typename const_statement_visitor_t::callsite_t;
  using intrinsic_t = typename const_statement_visitor_t::intrinsic_t;
  using arr_init_t = typename const_statement_visitor_t::arr_init_t;
  using arr_store_t = typename const_statement_visitor_t::arr_store_t;
  using arr_load_t = typename const_statement_visitor_t::arr_load_t;
  using arr_assign_t = typename const_statement_visitor_t::arr_assign_t;
  using bool_bin_op_t = typename const_statement_visitor_t::bool_bin_op_t;
  using bool_assign_cst_t = typename const_statement_visitor_t::bool_assign_cst_t;
  using bool_assign_var_t = typename const_statement_visitor_t::bool_assign_var_t;
  using bool_assume_t = typename const_statement_visitor_t::bool_assume_t;
  using bool_select_t = typename const_statement_visitor_t::bool_select_t;
  using bool_assert_t = typename const_statement_visitor_t::bool_assert_t;
  using region_init_t = typename const_statement_visitor_t::region_init_t;
  using region_copy_t = typename const_statement_visitor_t::region_copy_t;
  using region_cast_t = typename const_statement_visitor_t::region_cast_t;
  using make_ref_t = typename const_statement_visitor_t::make_ref_t;
  using remove_ref_t = typename const_statement_visitor_t::remove_ref_t;
  using load_from_ref_t = typename const_statement_visitor_t::load_from_ref_t;
  using store_to_ref_t = typename const_statement_visitor_t::store_to_ref_t;
  using gep_ref_t = typename const_statement_visitor_t::gep_ref_t;
  using assume_ref_t = typename const_statement_visitor_t::assume_ref_t;
  using assert_ref_t = typename const_statement_visitor_t::assert_ref_t;
  using select_ref_t = typename const_statement_visitor_t::select_ref_t;
  using int_to_ref_t = typename const_statement_visitor_t::int_to_ref_t;
  using ref_to_int_t = typename const_statement_visitor_t::ref_to_int_t;

  crab::json::writer &m_w;

  /** Attach the source location of a statement, when it has one. */
  void write_debug_info(const crab::cfg::debug_info &di) {
    if (!di.has_debug()) {
      m_w.kv_null("loc");
      return;
    }
    m_w.key("loc");
    m_w.begin_object(true /*compact*/);
    m_w.kv_string("file", di.get_file());
    m_w.kv_int("line", di.get_line());
    m_w.kv_int("col", di.get_column());
    m_w.kv_int("id", di.get_id());
    m_w.end_object();
  }

  void var(const char *key, const crab::variable<number_t, varname_t> &v) {
    m_w.key(key);
    crab::json::write(m_w, v);
  }

  void type(const char *key, const crab::variable_type &ty) {
    m_w.key(key);
    crab::json::write(m_w, ty);
  }

  void expr(const char *key,
            const ikos::linear_expression<number_t, varname_t> &e) {
    m_w.key(key);
    crab::json::write(m_w, e);
  }

  void cst(const char *key,
           const ikos::linear_constraint<number_t, varname_t> &c) {
    m_w.key(key);
    crab::json::write(m_w, c);
  }

  void ref_cst(const char *key,
               const crab::reference_constraint<number_t, varname_t> &c) {
    m_w.key(key);
    crab::json::write(m_w, c);
  }

  void var_or_cst(const char *key,
                  const crab::variable_or_constant<number_t, varname_t> &vc) {
    m_w.key(key);
    crab::json::write(m_w, vc);
  }

  /** An operand that a statement may leave unset; null when it is absent. */
  void opt_var(const char *key,
               const boost::optional<crab::variable<number_t, varname_t>> &v) {
    if (!v) {
      m_w.kv_null(key);
      return;
    }
    m_w.key(key);
    crab::json::write(m_w, *v);
  }

  void var_array(const char *key,
                 const std::vector<crab::variable<number_t, varname_t>> &vs) {
    m_w.key(key);
    m_w.begin_array();
    for (auto const &v : vs) {
      crab::json::write(m_w, v);
    }
    m_w.end_array();
  }

public:
  explicit statement_to_json_visitor(crab::json::writer &w) : m_w(w) {}

  // -- integer statements ----------------------------------------------------

  virtual void visit(const bin_op_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "binop");
    m_w.kv_string("op", op_name(s.op()));
    expr("left", s.left());
    expr("right", s.right());
    var("lhs", s.lhs());
    m_w.end_object();
  }

  virtual void visit(const assign_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "assign");
    expr("rhs", s.rhs());
    var("lhs", s.lhs());
    m_w.end_object();
  }

  virtual void visit(const assume_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "assume");
    cst("cond", s.constraint());
    m_w.end_object();
  }

  virtual void visit(const assert_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "assert");
    cst("cond", s.constraint());
    write_debug_info(s.get_debug_info());
    m_w.end_object();
  }

  virtual void visit(const select_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "select");
    cst("cond", s.cond());
    expr("left", s.left());
    expr("right", s.right());
    var("lhs", s.lhs());
    m_w.end_object();
  }

  virtual void visit(const int_cast_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "cast");
    m_w.kv_string("op", op_name(s.op()));
    var("rhs", s.src());
    m_w.kv_unsigned("rhs_width", s.src_width());
    var("lhs", s.dst());
    m_w.kv_unsigned("lhs_width", s.dst_width());
    m_w.end_object();
  }

  virtual void visit(const havoc_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "havoc");
    var("lhs", s.get_variable());
    m_w.end_object();
  }

  virtual void visit(const unreach_t &) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "unreachable");
    m_w.end_object();
  }

  // -- boolean statements ----------------------------------------------------

  virtual void visit(const bool_assign_cst_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "bool_assign_cst");
    // The right-hand side is either a linear or a reference constraint. They
    // are told apart by "cst_kind" rather than by the shape of "rhs", so a
    // consumer does not have to guess from the keys inside it.
    if (s.is_rhs_linear_constraint()) {
      m_w.kv_string("cst_kind", "linear");
      cst("rhs", s.rhs_as_linear_constraint());
    } else {
      m_w.kv_string("cst_kind", "reference");
      ref_cst("rhs", s.rhs_as_reference_constraint());
    }
    var("lhs", s.lhs());
    m_w.end_object();
  }

  virtual void visit(const bool_assign_var_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "bool_assign_var");
    var("rhs", s.rhs());
    m_w.kv_bool("negated", s.is_rhs_negated());
    var("lhs", s.lhs());
    m_w.end_object();
  }

  virtual void visit(const bool_bin_op_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "bool_binop");
    m_w.kv_string("op", op_name(s.op()));
    var("left", s.left());
    var("right", s.right());
    var("lhs", s.lhs());
    m_w.end_object();
  }

  virtual void visit(const bool_assume_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "bool_assume");
    var("cond", s.cond());
    m_w.kv_bool("negated", s.is_negated());
    m_w.end_object();
  }

  virtual void visit(const bool_assert_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "bool_assert");
    var("cond", s.cond());
    write_debug_info(s.get_debug_info());
    m_w.end_object();
  }

  virtual void visit(const bool_select_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "bool_select");
    var("cond", s.cond());
    var("left", s.left());
    var("right", s.right());
    var("lhs", s.lhs());
    m_w.end_object();
  }

  // -- arrays ----------------------------------------------------------------

  virtual void visit(const arr_init_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "array_init");
    expr("lb", s.lb_index());
    expr("ub", s.ub_index());
    expr("val", s.val());
    expr("elem_size", s.elem_size());
    var("array", s.array());
    m_w.end_object();
  }

  virtual void visit(const arr_store_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "array_store");
    expr("lb", s.lb_index());
    expr("ub", s.ub_index());
    expr("value", s.value());
    expr("elem_size", s.elem_size());
    var("array", s.array());
    m_w.end_object();
  }

  virtual void visit(const arr_load_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "array_load");
    var("array", s.array());
    expr("index", s.index());
    expr("elem_size", s.elem_size());
    var("lhs", s.lhs());
    m_w.end_object();
  }

  virtual void visit(const arr_assign_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "array_assign");
    var("rhs", s.rhs());
    var("lhs", s.lhs());
    m_w.end_object();
  }

  // -- calls -----------------------------------------------------------------

  virtual void visit(const callsite_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "callsite");
    m_w.kv_string("callee", s.get_func_name());
    var_array("args", s.get_args());
    var_array("lhs", s.get_lhs());
    m_w.end_object();
  }

  virtual void visit(const intrinsic_t &s) override {
    const std::string &name = s.get_intrinsic_name();
    // Directives for the value partitioning domain. They constrain how the
    // analysis explores the program, not what the program does, so they have no
    // concrete semantics and are omitted.
    if (name == "value_partition_start" || name == "value_partition_end") {
      return;
    }
    CRAB_ERROR("cfg_to_json: intrinsic '", name,
               "' is not supported by the JSON exporter");
  }

  // -- regions and references ------------------------------------------------

  virtual void visit(const region_init_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "region_init");
    var("lhs", s.region());
    m_w.end_object();
  }

  virtual void visit(const region_copy_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "region_copy");
    var("rhs", s.rhs_region());
    var("lhs", s.lhs_region());
    m_w.end_object();
  }

  virtual void visit(const region_cast_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "region_cast");
    var("rhs", s.src());
    type("rhs_type", s.src_type());
    var("lhs", s.dst());
    type("lhs_type", s.dst_type());
    m_w.end_object();
  }

  virtual void visit(const make_ref_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "make_ref");
    var_or_cst("size", s.size());
    // The allocation site identifies which abstract object the new reference
    // points to; only its index is meaningful outside Crab.
    m_w.kv_unsigned("alloc_site", s.alloc_site().index());
    var("lhs_region", s.region());
    var("lhs", s.lhs());
    write_debug_info(s.get_debug_info());
    m_w.end_object();
  }

  virtual void visit(const remove_ref_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "remove_ref");
    var("region", s.region());
    var("ref", s.ref());
    m_w.end_object();
  }

  virtual void visit(const load_from_ref_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "load_from_ref");
    var("region", s.region());
    var("ref", s.ref());
    var("lhs", s.lhs());
    write_debug_info(s.get_debug_info());
    m_w.end_object();
  }

  virtual void visit(const store_to_ref_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "store_to_ref");
    var("region", s.region());
    var("ref", s.ref());
    var_or_cst("value", s.val());
    write_debug_info(s.get_debug_info());
    m_w.end_object();
  }

  virtual void visit(const gep_ref_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "gep_ref");
    var("rhs_region", s.rhs_region());
    var("rhs", s.rhs());
    expr("offset", s.offset());
    var("lhs_region", s.lhs_region());
    var("lhs", s.lhs());
    write_debug_info(s.get_debug_info());
    m_w.end_object();
  }

  virtual void visit(const assume_ref_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "assume_ref");
    ref_cst("cond", s.constraint());
    m_w.end_object();
  }

  virtual void visit(const assert_ref_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "assert_ref");
    ref_cst("cond", s.constraint());
    write_debug_info(s.get_debug_info());
    m_w.end_object();
  }

  virtual void visit(const select_ref_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "select_ref");
    var("cond", s.cond());
    // The region of an operand is absent when that operand is the null
    // constant, which has no region to speak of.
    opt_var("left_region", s.left_rgn());
    var_or_cst("left", s.left_ref());
    opt_var("right_region", s.right_rgn());
    var_or_cst("right", s.right_ref());
    var("lhs_region", s.lhs_rgn());
    var("lhs", s.lhs_ref());
    write_debug_info(s.get_debug_info());
    m_w.end_object();
  }

  virtual void visit(const int_to_ref_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "int_to_ref");
    var("rhs", s.int_var());
    var("lhs_region", s.region());
    var("lhs_ref", s.ref_var());
    write_debug_info(s.get_debug_info());
    m_w.end_object();
  }

  virtual void visit(const ref_to_int_t &s) override {
    m_w.begin_object();
    m_w.kv_string("stmt", "ref_to_int");
    var("rhs_region", s.region());
    var("rhs_ref", s.ref_var());
    var("lhs", s.int_var());
    write_debug_info(s.get_debug_info());
    m_w.end_object();
  }
};

} // end namespace json_impl

namespace json_impl {
/** Default for the `extra` hook below: add nothing to a block. */
struct no_extra_block_fields {
  template <typename Label>
  void operator()(const Label &, crab::json::writer &) const {}
};
} // end namespace json_impl

/**
 * Serialize a single basic block.
 *
 * `extra` is called after "stmts" and before "succs", and may add further keys
 * to the block object -- the invariant holding at the block, say. It must leave
 * the writer inside that same object.
 */
template <typename CFG,
          typename ExtraFields = json_impl::no_extra_block_fields>
void basic_block_to_json(const CFG &cfg,
                         const typename CFG::basic_block_label_t &label,
                         crab::json::writer &w,
                         ExtraFields extra = ExtraFields()) {
  using basic_block_t = typename CFG::basic_block_t;

  w.begin_object();
  w.kv_string("label", basic_block_traits<basic_block_t>::to_string(label));

  w.key("stmts");
  w.begin_array();
  const basic_block_t &bb = cfg.get_node(label);
  for (auto const&s : bb) {
    json_impl::statement_to_json_visitor<CFG> vis(w);
    s.accept(&vis);
  }
  w.end_array();

  extra(label, w);

  w.key("succs");
  w.begin_array(true /*compact*/);
  for (auto const &succ : cfg.next_nodes(label)) {
    w.value_string(basic_block_traits<basic_block_t>::to_string(succ));
  }
  w.end_array();

  w.end_object();
}

/**
 * Serialize a CFG as JSON onto `w`.
 *
 * Emits one object; the caller decides what document it lives in. `extra` is
 * forwarded to every block: see basic_block_to_json.
 */
template <typename CFG,
          typename ExtraFields = json_impl::no_extra_block_fields>
void cfg_to_json(const CFG &cfg, crab::json::writer &w,
                 ExtraFields extra = ExtraFields()) {
  using basic_block_t = typename CFG::basic_block_t;

  w.begin_object();

  if (cfg.has_func_decl()) {
    auto const &decl = cfg.get_func_decl();
    w.kv_string("name", decl.get_func_name());
    w.key("declaration");
    w.begin_object();
    w.key("inputs");
    w.begin_array();
    for (auto const &v : decl.get_inputs()) {
      crab::json::write(w, v);
    }
    w.end_array();
    w.key("outputs");
    w.begin_array();
    for (auto const &v : decl.get_outputs()) {
      crab::json::write(w, v);
    }
    w.end_array();
    w.end_object();
  } else {
    w.kv_null("name");
    w.kv_null("declaration");
  }

  w.kv_string("entry",
              basic_block_traits<basic_block_t>::to_string(cfg.entry()));
  if (cfg.has_exit()) {
    w.kv_string("exit",
                basic_block_traits<basic_block_t>::to_string(cfg.exit()));
  } else {
    w.kv_null("exit");
  }

  // Sorted by label: see the note on stability in the file comment.
  std::vector<std::pair<std::string, typename CFG::basic_block_label_t>> labels;
  for (auto it = cfg.label_begin(), et = cfg.label_end(); it != et; ++it) {
    labels.emplace_back(basic_block_traits<basic_block_t>::to_string(*it), *it);
  }
  std::sort(labels.begin(), labels.end(),
            [](const std::pair<std::string,
                               typename CFG::basic_block_label_t> &x,
               const std::pair<std::string,
                               typename CFG::basic_block_label_t> &y) {
              return x.first < y.first;
            });

  w.key("blocks");
  w.begin_array();
  for (auto const &l : labels) {
    basic_block_to_json(cfg, l.second, w, extra);
  }
  w.end_array();

  w.end_object();
}

/** Convenience: serialize a CFG as a standalone JSON document. */
template <typename CFG> void cfg_to_json(const CFG &cfg, crab::crab_os &os) {
  crab::json::writer w(os);
  cfg_to_json(cfg, w);
  w.finish();
}

} // end namespace cfg
} // end namespace crab
