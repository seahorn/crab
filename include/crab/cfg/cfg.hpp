#pragma once

/*
 * Build a CFG to interface with the abstract domains and fixpoint
 * iterators.
 *
 * === Syntax ===
 *
 * The CFG is intended to be expressive enough for representing the
 * semantics of a variety of programming languages and with different
 * levels of abstractions.  All the CFG statements are strongly
 * typed. However, only variables need to be typed. The types of
 * constants can be inferred from the context since they always appear
 * together with at least one variable. Types form a **flat** lattice
 * consisting of:
 *
 * - booleans,
 * - integers,
 * - reals,
 * - references,
 * - array of booleans/integers/reals,
 * - regions of booleans/integers/reals
 * - regions of references.
 * - regions of arrays of booleans/integers/reals
 *
 * Crab CFG supports the modeling of:
 *
 *   - arithmetic operations over integers or reals,
 *   - boolean operations,
 *   - C-like pointers (through regions and references),
 *   - uni-dimensional arrays of booleans, integers or reals
 *     (useful for C-like arrays and heap abstractions).
 *   - and functions declarations and callsites
 *
 * Since the CFG is strongly typed, it provides several methods to
 * convert betwen types:
 * 
 *   - between integers and booleans
 *   - between integers and references
 *   - between integers with different bitwidths
 * 
 * === Semantics ===
 *
 * The semantics of boolean and numerical is similar to a language
 * like C. For arrays, the main feature is that arrays are guaranteed
 * to be disjoint from each other (i.e., a write in array A cannot
 * affect the contents of another array B). The semantics of arrays
 * reads and writes are similar to SMT arrays.
 *
 * The semantics for regions and references is more novel and it is
 * described in the wiki (https://github.com/seahorn/crab/wiki).
 *
 * === Function calls ===
 *
 * In Crab, a function is represented by a CFG with a special function
 * declaration statement that describes the formal parameters INPUTS
 * and return parameters OUTPUTS.  Very importantly, Crab will report
 * an error if INPUTS and OUTPUTS are not disjoint. This assumption
 * about disjointness makes easier the inter-procedural analyses.
 *
 * Crab allows to have multiple CFGs in memory and it combines them in
 * a call graph where each node is a Crab CFG and there is an edge
 * from n1 to n2 if n1 has a callsite that calls the function at n2.
 *
 * If you want to use Crab's interprocedural analysis is important to
 * know that the analysis will compute summaries that only express
 * relationships over INPUT and OUTPUT variables. This affects in how
 * you should encode your program into the Crab CFG language.
 *
 * For instance, if you have a function `foo(p)` where `p` is a
 * reference. If foo updates the content of p then you should write the
 * CFG for `foo` as follows:
 *
 *       function_decl(foo, p_in, p_out);
 *       gep_ref(p_out, p_in);
 *       ... // rest of the CFG reads or writes p_out and keeps intact p_in
 *
 * The same things should be done with arrays:
 *
 *       function_decl(foo, a_in, a_out);
 *       array_assign(a_out, a_in);
 *       ... // rest of the CFG reads or writes a_out and keeps intact a_in
 *
 * === Important note ===
 *
 * Objects of the class cfg are not copyable. Instead, we provide a
 * class cfg_ref that wraps cfg references into copyable and
 * assignable objects.
 *
 * === CFG Limitations ===
 *
 * - The CFG language does not allow to express floating point
 *   operations.
 * - Integer and reals cannot be mixed.
 */

#include <crab/cfg/basic_block_traits.hpp>
#include <crab/cfg/cfg_operators.hpp>
#include <crab/domains/discrete_domains.hpp>
#include <crab/domains/interval.hpp>
#include <crab/numbers/bignums.hpp>
#include <crab/types/linear_constraints.hpp>
#include <crab/types/reference_constraints.hpp>
#include <crab/types/variable.hpp>
#include <crab/types/tag.hpp>

#include <boost/iterator/indirect_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <unordered_map>
#include <unordered_set>

namespace crab {
namespace cfg {

enum stmt_code {
  UNDEF = 0,
  // numerical
  BIN_OP = 20,
  ASSIGN = 21,
  ASSUME = 22,
  UNREACH = 23,
  SELECT = 24,
  ASSERT = 25,
  // arrays
  ARR_INIT = 30,
  ARR_STORE = 31,
  ARR_LOAD = 32,
  ARR_ASSIGN = 33,
  // regions and references
  REF_MAKE = 40,
  REF_LOAD = 41,
  REF_STORE = 42,
  REF_GEP = 43,
  REF_ASSUME = 44,
  REF_ASSERT = 45,
  REF_SELECT = 46,
  REGION_INIT = 47,
  REF_TO_INT = 48,
  INT_TO_REF = 49,
  REGION_COPY = 50,
  REGION_CAST = 51,  
  REF_REMOVE = 52,
  // functions calls
  CALLSITE = 60,
  CRAB_INTRINSIC = 61,
  // integers/arrays/pointers/boolean
  HAVOC = 70,
  // boolean
  BOOL_BIN_OP = 80,
  BOOL_ASSIGN_CST = 81,
  BOOL_ASSIGN_VAR = 82,
  BOOL_ASSUME = 83,
  BOOL_SELECT = 84,
  BOOL_ASSERT = 85,
  // casts
  INT_CAST = 90
};

template <typename Number, typename VariableName>
class live {
public:
  using variable_t = variable<Number, VariableName>;

private:
  using live_set_t = std::vector<variable_t>;

public:
  using const_use_iterator = typename live_set_t::const_iterator;
  using const_def_iterator = typename live_set_t::const_iterator;

private:
  live_set_t m_uses;
  live_set_t m_defs;

  void add(live_set_t &s, variable_t v) {
    auto it = find(s.begin(), s.end(), v);
    if (it == s.end())
      s.push_back(v);
  }

public:
  live() {}

  void add_use(const variable_t v) { add(m_uses, v); }
  void add_def(const variable_t v) { add(m_defs, v); }

  const_use_iterator uses_begin() const { return m_uses.begin(); }
  const_use_iterator uses_end() const { return m_uses.end(); }
  const_use_iterator defs_begin() const { return m_defs.begin(); }
  const_use_iterator defs_end() const { return m_defs.end(); }

  size_t num_uses() const { return m_uses.size(); }
  size_t num_defs() const { return m_defs.size(); }

  friend crab_os &operator<<(crab_os &o,
                             const live<Number, VariableName> &live) {
    o << "Use={";
    for (auto const &v :
         boost::make_iterator_range(live.uses_begin(), live.uses_end()))
      o << v << ",";
    o << "} Def={";
    for (auto const &v :
         boost::make_iterator_range(live.defs_begin(), live.defs_end()))
      o << v << ",";
    o << "}";
    return o;
  }
};

class debug_info {
private:
  std::string m_file;
  int m_line;
  int m_col;
  int64_t m_id; 
public:
  debug_info()
    : m_file(""), m_line(-1), m_col(-1), m_id(-1) {}
  debug_info(std::string file, unsigned line, unsigned col, int64_t num_id)
    : m_file(file), m_line(line), m_col(col), m_id(num_id) {}
  debug_info(int64_t num_id)
    : m_file(""), m_line(-1), m_col(-1), m_id(num_id) {}
  
  bool operator<(const debug_info &other) const {
    if (m_id < other.m_id) {
      return true;
    } else if (m_id > other.m_id) {
      return false;
    }
    if (m_file < other.m_file) {
      return true;
    } else if (m_file > other.m_file) {
      return false;
    }
    if (m_line < other.m_line) {
      return true;
    } else if (m_line > other.m_line) {
      return false;
    }
    return (m_col < other.m_col);
  }

  bool operator==(const debug_info &other) const {
    return (m_id == other.m_id &&
	    m_file == other.m_file &&
	    m_line == other.m_line &&
            m_col == other.m_col);
  }

  bool has_debug() const {
    return ((m_id >= 0) ||
	    ((m_file != "") && (m_line >= 0) && (m_col >= 0)));
  }

  std::string get_file() const { return m_file; }
  int get_line() const { return m_line; }
  int get_column() const { return m_col; }
  int64_t get_id() const { return m_id; }
  
  void write(crab_os &o) const {
    o << "loc("
      << "file=";
    if (m_file == "") {
      o << "none";
    } else {
      o << m_file;
    }
    o << " "
      << "line=" << m_line << " "
      << "col=" << m_col << ") ";
    if (m_id >= 0) {
      o << "id=" << m_id;
    }
  }
};

inline crab_os &operator<<(crab_os &o, const debug_info &l) {
  l.write(o);
  return o;
}
  
template <class BasicBlockLabel, class VariableName, class Number>
class basic_block;

template <class BasicBlockLabel, class Number, class VariableName>
struct statement_visitor;

template <class BasicBlockLabel, class Number, class VariableName>
class statement {
public:
  using live_t = live<Number, VariableName>;
  using basic_block_t = basic_block<BasicBlockLabel, VariableName, Number>;

protected:
  live_t m_live;
  stmt_code m_stmt_code;
  basic_block_t *m_parent;
  debug_info m_dbg_info;

  statement(stmt_code code, basic_block_t *parent,
            debug_info dbg_info = debug_info())
      : m_stmt_code(code), m_parent(parent), m_dbg_info(dbg_info) {}

public:
  virtual ~statement() = default;

  bool is_bin_op() const { return (m_stmt_code == BIN_OP); }
  bool is_assign() const { return (m_stmt_code == ASSIGN); }
  bool is_assume() const { return (m_stmt_code == ASSUME); }
  bool is_select() const { return (m_stmt_code == SELECT); }
  bool is_assert() const { return (m_stmt_code == ASSERT); }
  bool is_int_cast() const { return (m_stmt_code == INT_CAST); }
  bool is_havoc() const { return m_stmt_code == HAVOC; }
  bool is_unreachable() const { return m_stmt_code == UNREACH; }
  bool is_arr_init() const { return (m_stmt_code == ARR_INIT); }
  bool is_arr_read() const { return (m_stmt_code == ARR_LOAD); }
  bool is_arr_write() const { return (m_stmt_code == ARR_STORE); }
  bool is_arr_assign() const { return (m_stmt_code == ARR_ASSIGN); }
  bool is_ref_make() const { return m_stmt_code == REF_MAKE; }
  bool is_ref_remove() const { return m_stmt_code == REF_REMOVE;}
  bool is_ref_load() const { return m_stmt_code == REF_LOAD; }
  bool is_ref_store() const { return m_stmt_code == REF_STORE; }
  bool is_ref_gep() const { return m_stmt_code == REF_GEP; }
  bool is_ref_assume() const { return (m_stmt_code == REF_ASSUME); }
  bool is_ref_assert() const { return (m_stmt_code == REF_ASSERT); }
  bool is_ref_select() const { return (m_stmt_code == REF_SELECT); }
  bool is_ref_to_int() const { return (m_stmt_code == REF_TO_INT); }
  bool is_int_to_ref() const { return (m_stmt_code == INT_TO_REF); }
  bool is_region_init() const { return (m_stmt_code == REGION_INIT); }
  bool is_region_copy() const { return (m_stmt_code == REGION_COPY); }
  bool is_region_cast() const { return (m_stmt_code == REGION_CAST); }  
  bool is_bool_bin_op() const { return (m_stmt_code == BOOL_BIN_OP); }
  bool is_bool_assign_cst() const { return (m_stmt_code == BOOL_ASSIGN_CST); }
  bool is_bool_assign_var() const { return (m_stmt_code == BOOL_ASSIGN_VAR); }
  bool is_bool_assume() const { return (m_stmt_code == BOOL_ASSUME); }
  bool is_bool_assert() const { return (m_stmt_code == BOOL_ASSERT); }
  bool is_bool_select() const { return (m_stmt_code == BOOL_SELECT); }
  bool is_callsite() const { return (m_stmt_code == CALLSITE); }
  bool is_intrinsic() const { return (m_stmt_code == CRAB_INTRINSIC); }

  const live_t &get_live() const { return m_live; }

  const debug_info &get_debug_info() const { return m_dbg_info; }

  basic_block_t *get_parent() { return m_parent; }

  const basic_block_t *get_parent() const { return m_parent; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *) = 0;

  virtual void write(crab_os &o) const = 0;

  virtual statement<BasicBlockLabel, Number, VariableName> *
  clone(basic_block_t *parent) const = 0;

  // for gdb
  void dump() const { write(crab::errs()); }

  friend crab_os &
  operator<<(crab_os &o,
             const statement<BasicBlockLabel, Number, VariableName> &s) {
    s.write(o);
    return o;
  }
};

/*
 *  Numerical statements
 */

template <class BasicBlockLabel, class Number, class VariableName>
class binary_op : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = binary_op<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using linear_expression_t = ikos::linear_expression<Number, VariableName>;

  binary_op(variable_t lhs, binary_operation_t op, linear_expression_t op1,
            linear_expression_t op2, basic_block_t *parent,
            debug_info dbg_info = debug_info())
      : statement_t(BIN_OP, parent, dbg_info), m_lhs(lhs), m_op(op), m_op1(op1),
        m_op2(op2) {
    this->m_live.add_def(m_lhs);
    for (auto const &v : m_op1.variables()) {
      this->m_live.add_use(v);
    }
    for (auto const &v : m_op2.variables()) {
      this->m_live.add_use(v);
    }
  }

  const variable_t &lhs() const { return m_lhs; }

  binary_operation_t op() const { return m_op; }

  const linear_expression_t &left() const { return m_op1; }

  const linear_expression_t &right() const { return m_op2; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_op, m_op1, m_op2, parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << m_lhs << " = " << m_op1 << m_op << m_op2;
  }

private:
  variable_t m_lhs;
  binary_operation_t m_op;
  linear_expression_t m_op1;
  linear_expression_t m_op2;
};

template <class BasicBlockLabel, class Number, class VariableName>
class assignment : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = assignment<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using linear_expression_t = ikos::linear_expression<Number, VariableName>;

  assignment(variable_t lhs, linear_expression_t rhs, basic_block_t *parent)
      : statement_t(ASSIGN, parent), m_lhs(lhs), m_rhs(rhs) {
    this->m_live.add_def(m_lhs);
    for (auto const &v : m_rhs.variables())
      this->m_live.add_use(v);
  }

  const variable_t &lhs() const { return m_lhs; }

  const linear_expression_t &rhs() const { return m_rhs; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_rhs, parent);
  }

  virtual void write(crab_os &o) const override {
    o << m_lhs << " = " << m_rhs; // << " " << this->m_live;
  }

private:
  variable_t m_lhs;
  linear_expression_t m_rhs;
};

template <class BasicBlockLabel, class Number, class VariableName>
class assume_stmt : public statement<BasicBlockLabel, Number, VariableName> {

  using this_type = assume_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using linear_constraint_t = ikos::linear_constraint<Number, VariableName>;

  assume_stmt(linear_constraint_t cst, basic_block_t *parent)
      : statement_t(ASSUME, parent), m_cst(cst) {
    for (auto const &v : cst.variables())
      this->m_live.add_use(v);
  }

  const linear_constraint_t &constraint() const { return m_cst; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_cst, parent);
  }

  virtual void write(crab_os &o) const override {
    o << "assume(" << m_cst << ")"; //  << " " << this->m_live;
  }

private:
  linear_constraint_t m_cst;
};

template <class BasicBlockLabel, class Number, class VariableName>
class unreachable_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = unreachable_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;

  unreachable_stmt(basic_block_t *parent) : statement_t(UNREACH, parent) {}

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(parent);
  }

  virtual void write(crab_os &o) const override  {
    o << "unreachable";
  }
};

template <class BasicBlockLabel, class Number, class VariableName>
class havoc_stmt : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = havoc_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;

  havoc_stmt(variable_t lhs, std::string comment, basic_block_t *parent)
      : statement_t(HAVOC, parent), m_lhs(lhs), m_comment(comment) {
    this->m_live.add_def(m_lhs);
  }

  const variable_t &get_variable() const { return m_lhs; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_comment, parent);
  }

  virtual void write(crab_os &o) const override {
    o << "havoc(" << m_lhs << ")";
    if (!m_comment.empty()) {
      o << " /* " << m_comment << "*/";
    }
  }

private:
  variable_t m_lhs;
  std::string m_comment;
};

// select x, c, e1, e2:
//    if c > 0 then x=e1 else x=e2
//
// Note that a select instruction is not strictly needed and can be
// simulated by splitting blocks. However, frontends like LLVM can
// generate many select instructions so we prefer to support
// natively to avoid a blow up in the size of the CFG.
template <class BasicBlockLabel, class Number, class VariableName>
class select_stmt : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = select_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using linear_expression_t = ikos::linear_expression<Number, VariableName>;
  using linear_constraint_t = ikos::linear_constraint<Number, VariableName>;

  select_stmt(variable_t lhs, linear_constraint_t cond, linear_expression_t e1,
              linear_expression_t e2, basic_block_t *parent)
      : statement_t(SELECT, parent), m_lhs(lhs), m_cond(cond), m_e1(e1),
        m_e2(e2) {
    this->m_live.add_def(m_lhs);
    for (auto const &v : m_cond.variables())
      this->m_live.add_use(v);
    for (auto const &v : m_e1.variables())
      this->m_live.add_use(v);
    for (auto const &v : m_e2.variables())
      this->m_live.add_use(v);
  }

  const variable_t &lhs() const { return m_lhs; }

  const linear_constraint_t &cond() const { return m_cond; }

  const linear_expression_t &left() const { return m_e1; }

  const linear_expression_t &right() const { return m_e2; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_cond, m_e1, m_e2, parent);
  }

  virtual void write(crab_os &o) const override {
    o << m_lhs << " = "
      << "ite(" << m_cond << "," << m_e1 << "," << m_e2 << ")";
  }

private:
  variable_t m_lhs;
  linear_constraint_t m_cond;
  linear_expression_t m_e1;
  linear_expression_t m_e2;
};

template <class BasicBlockLabel, class Number, class VariableName>
class assert_stmt : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = assert_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using linear_constraint_t = ikos::linear_constraint<Number, VariableName>;

  assert_stmt(linear_constraint_t cst, basic_block_t *parent,
              debug_info dbg_info = debug_info())
      : statement_t(ASSERT, parent, dbg_info), m_cst(cst) {
    for (auto const &v : cst.variables())
      this->m_live.add_use(v);
  }

  const linear_constraint_t &constraint() const { return m_cst; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_cst, parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << "assert(" << m_cst << ")";
    if (this->m_dbg_info.has_debug()) {
      o << "   /* " << this->m_dbg_info <<  " */";
    }
  }

private:
  linear_constraint_t m_cst;
};

template <class BasicBlockLabel, class Number, class VariableName>
class int_cast_stmt : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = int_cast_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using bitwidth_t = typename variable_t::bitwidth_t;

  int_cast_stmt(cast_operation_t op, variable_t src, variable_t dst,
                basic_block_t *parent, debug_info dbg_info = debug_info())
      : statement_t(INT_CAST, parent, dbg_info), m_op(op), m_src(src),
        m_dst(dst) {
    this->m_live.add_use(m_src);
    this->m_live.add_def(m_dst);
  }

  cast_operation_t op() const { return m_op; }
  const variable_t &src() const { return m_src; }
  bitwidth_t src_width() const { return get_bitwidth(m_src); }
  const variable_t &dst() const { return m_dst; }
  bitwidth_t dst_width() const { return get_bitwidth(m_dst); }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_op, m_src, m_dst, parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    // bitwidths are casted to int, otherwise operator<< may try
    // to print them as characters if bitwidth_t = uint8_t
    o << m_op << " " << m_src << ":" << (int)src_width() << " to " << m_dst
      << ":" << (int)dst_width();
  }

private:
  static bitwidth_t get_bitwidth(const variable_t &v) {
    auto ty = v.get_type();
    assert(ty.is_integer() || ty.is_bool());
    return (ty.is_integer() ? ty.get_integer_bitwidth() : 1);
  }

  cast_operation_t m_op;
  variable_t m_src;
  variable_t m_dst;
};

/*
 *  Array statements
 */

// XXX: the array statements array_init_stmt, array_store_stmt,
// and array_load_stmt take as one of their template parameters
// Number which is the type for the array indexes. Although array
// indexes should be always integers we keep it as a template
// parameter in case an analysis over a type different from
// integers (e.g., reals) is done. Note that we don't allow mixing
// non-integers and integers so we cannot have analysis where all
// variables are non-integers except array indexes.
//
// Each of these statements requires an element size, that is, the
// number of bytes that are being accessed. If the front-end is
// LLVM, then the element size is always known at compilation
// time. However, with other front-ends (e.g., BPF programs) the
// element size is stored in a variable so that's why the type of
// the element size is not just a constant integer but it can also
// be a variable.

//! Initialize all array elements to some variable or number.
//  The semantics is similar to constant arrays in SMT.
template <class BasicBlockLabel, class Number, class VariableName>
class array_init_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = array_init_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using linear_expression_t = ikos::linear_expression<Number, VariableName>;
  using variable_t = variable<Number, VariableName>;
  using type_t = typename variable_t::type_t;

  array_init_stmt(variable_t arr, linear_expression_t elem_size,
                  linear_expression_t lb, linear_expression_t ub,
                  linear_expression_t val, basic_block_t *parent)
      : statement_t(ARR_INIT, parent), m_arr(arr), m_elem_size(elem_size),
        m_lb(lb), m_ub(ub), m_val(val) {

    this->m_live.add_def(m_arr);
    for (auto const &v : m_elem_size.variables()) {
      this->m_live.add_use(v);
    }
    for (auto const &v : m_lb.variables()) {
      this->m_live.add_use(v);
    }
    for (auto const &v : m_ub.variables()) {
      this->m_live.add_use(v);
    }
    for (auto const &v : m_val.variables()) {
      this->m_live.add_use(v);
    }
  }

  const variable_t &array() const { return m_arr; }

  type_t array_type() const { return m_arr.get_type(); }

  const linear_expression_t &elem_size() const { return m_elem_size; }

  const linear_expression_t &lb_index() const { return m_lb; }

  const linear_expression_t &ub_index() const { return m_ub; }

  const linear_expression_t &val() const { return m_val; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_arr, m_elem_size, m_lb, m_ub, m_val, parent);
  }

  void write(crab_os &o) const override {
    o << m_arr << "[" << m_lb << "..." << m_ub << "] := " << m_val;
  }

private:
  // forall i \in [lb,ub] % elem_size :: arr[i] := val and
  // forall j < lb or j > ub :: arr[j] is undefined.
  variable_t m_arr;
  linear_expression_t m_elem_size; //! size in bytes
  linear_expression_t m_lb;
  linear_expression_t m_ub;
  linear_expression_t m_val;
};

template <class BasicBlockLabel, class Number, class VariableName>
class array_store_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = array_store_stmt<BasicBlockLabel, Number, VariableName>;

public:
  // forall i \in [lb,ub] % elem_size :: arr[i] := val

  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using linear_expression_t = ikos::linear_expression<Number, VariableName>;
  using variable_t = variable<Number, VariableName>;
  using type_t = typename variable_t::type_t;

  // Constructor for in-place array stores
  array_store_stmt(variable_t arr, linear_expression_t elem_size,
                   linear_expression_t lb, linear_expression_t ub,
                   linear_expression_t value, bool is_strong_update,
                   basic_block_t *parent)
      : statement_t(ARR_STORE, parent), m_arr(arr), m_elem_size(elem_size),
        m_lb(lb), m_ub(ub), m_value(value),
        m_is_strong_update(is_strong_update) {

    this->m_live.add_def(m_arr);
    this->m_live.add_use(m_arr);
    for (auto const &v : m_elem_size.variables()) {
      this->m_live.add_use(v);
    }
    for (auto const &v : m_lb.variables()) {
      this->m_live.add_use(v);
    }
    for (auto const &v : m_ub.variables()) {
      this->m_live.add_use(v);
    }
    for (auto const &v : m_value.variables()) {
      this->m_live.add_use(v);
    }
  }

  // Return the array variable.
  const variable_t &array() const { return m_arr; }

  const linear_expression_t &lb_index() const { return m_lb; }

  const linear_expression_t &ub_index() const { return m_ub; }

  const linear_expression_t &value() const { return m_value; }

  type_t array_type() const { return m_arr.get_type(); }

  const linear_expression_t &elem_size() const { return m_elem_size; }

  bool is_strong_update() const { return m_is_strong_update; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_arr, m_elem_size, m_lb, m_ub, m_value,
                         m_is_strong_update, parent);
  }

  virtual void write(crab_os &o) const override {
    o << "array_store(";
    if (m_lb.equal(m_ub)) {
      o << m_arr << "," << m_lb << "," << m_value << ",sz=" << elem_size();
    } else {
      o << m_arr << "," << m_lb << ".." << m_ub << "," << m_value
        << ",sz=" << elem_size();
    }
    o << ")";
  }

private:
  // -- The array variable.
  variable_t m_arr;
  // the size of the array element
  linear_expression_t m_elem_size; //! size in bytes
  // smallest array index
  linear_expression_t m_lb;
  // largest array index
  linear_expression_t m_ub;
  // value stored in the array
  linear_expression_t m_value;
  // whether the store is a strong update. This might help the
  // abstract domain. If unknown set to false.  Only makes sense
  // if m_lb is equal to m_ub.
  bool m_is_strong_update;
};

template <class BasicBlockLabel, class Number, class VariableName>
class array_load_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = array_load_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using linear_expression_t = ikos::linear_expression<Number, VariableName>;
  using variable_t = variable<Number, VariableName>;
  using type_t = typename variable_t::type_t;

  array_load_stmt(variable_t lhs, variable_t arr, linear_expression_t elem_size,
                  linear_expression_t index, basic_block_t *parent)
      : statement_t(ARR_LOAD, parent), m_lhs(lhs), m_array(arr),
        m_elem_size(elem_size), m_index(index) {

    this->m_live.add_def(lhs);
    this->m_live.add_use(m_array);
    for (auto const &v : m_elem_size.variables()) {
      this->m_live.add_use(v);
    }
    for (auto const &v : m_index.variables()) {
      this->m_live.add_use(v);
    }
  }

  const variable_t &lhs() const { return m_lhs; }

  const variable_t &array() const { return m_array; }

  type_t array_type() const { return m_array.get_type(); }

  const linear_expression_t &index() const { return m_index; }

  const linear_expression_t &elem_size() const { return m_elem_size; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_array, m_elem_size, m_index, parent);
  }

  virtual void write(crab_os &o) const override {
    o << m_lhs << " = "
      << "array_load(" << m_array << "," << m_index << ",sz=" << elem_size()
      << ")";
  }

private:
  variable_t m_lhs;
  variable_t m_array;
  linear_expression_t m_elem_size; //! size in bytes
  linear_expression_t m_index;
};

template <class BasicBlockLabel, class Number, class VariableName>
class array_assign_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  //! a = b
  using this_type = array_assign_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using type_t = typename variable_t::type_t;

  array_assign_stmt(variable_t lhs, variable_t rhs, basic_block_t *parent)
      : statement_t(ARR_ASSIGN, parent), m_lhs(lhs), m_rhs(rhs) {
    this->m_live.add_def(lhs);
    this->m_live.add_use(rhs);
  }

  const variable_t &lhs() const { return m_lhs; }

  const variable_t &rhs() const { return m_rhs; }

  type_t array_type() const { return m_lhs.get_type(); }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_rhs, parent);
  }

  virtual void write(crab_os &o) const override {
    o << "array_assign(" << m_lhs << ", " << m_rhs << ")";
  }

private:
  variable_t m_lhs;
  variable_t m_rhs;
};

/*
 * References and Regions
 */
  
template <class BasicBlockLabel, class Number, class VariableName>
class region_init_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  // region_init(region);
  using this_type = region_init_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;

  region_init_stmt(variable_t region, basic_block_t *parent,
                   debug_info dbg_info = debug_info())
      : statement_t(REGION_INIT, parent, dbg_info), m_region(region) {
    this->m_live.add_def(m_region);
  }

  const variable_t &region() const { return m_region; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_region, parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << "region_init(" << m_region << ":" << m_region.get_type() << ")";
  }

private:
  variable_t m_region;
};

template <class BasicBlockLabel, class Number, class VariableName>
class region_copy_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  //! a = b
  using this_type = region_copy_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using type_t = typename variable_t::type_t;

  region_copy_stmt(variable_t lhs_region, variable_t rhs_region,
                   basic_block_t *parent)
      : statement_t(REGION_COPY, parent), m_lhs_region(lhs_region),
        m_rhs_region(rhs_region) {
    this->m_live.add_def(m_lhs_region);
    this->m_live.add_use(m_rhs_region);
  }

  const variable_t &lhs_region() const { return m_lhs_region; }

  const variable_t &rhs_region() const { return m_rhs_region; }

  type_t region_type() const { return m_lhs_region.get_type(); }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs_region, m_rhs_region, parent);
  }

  virtual void write(crab_os &o) const override {
    o << "region_copy(" << m_lhs_region << ":" << m_lhs_region.get_type()
      << ", " << m_rhs_region << ":" << m_rhs_region.get_type() << ")";
  }

private:
  variable_t m_lhs_region;
  variable_t m_rhs_region;
};


template <class BasicBlockLabel, class Number, class VariableName>
class region_cast_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = region_cast_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using type_t = typename variable_t::type_t;

  region_cast_stmt(variable_t src_region, variable_t dst_region,
                   basic_block_t *parent)
      : statement_t(REGION_CAST, parent),
	m_src_region(src_region),
        m_dst_region(dst_region) {
    this->m_live.add_use(m_src_region);    
    this->m_live.add_def(m_dst_region);
  }

  const variable_t &src() const { return m_src_region; }

  const variable_t &dst() const { return m_dst_region; }

  type_t src_type() const { return m_src_region.get_type(); }
  type_t dst_type() const { return m_dst_region.get_type(); }  

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_src_region, m_dst_region, parent);
  }

  virtual void write(crab_os &o) const override {
    o << "region_cast "
      << src() << ":" << src_type() << " to "
      << dst() << ":" << dst_type();
  }

private:
  variable_t m_src_region;
  variable_t m_dst_region;
};
  
template <class BasicBlockLabel, class Number, class VariableName>
class make_ref_stmt : public statement<BasicBlockLabel, Number, VariableName> {
  // lhs := make_ref(region, as)
  using this_type = make_ref_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using variable_or_constant_t = variable_or_constant<Number, VariableName>;
  
  make_ref_stmt(variable_t lhs, variable_t region,
		variable_or_constant_t size,
		crab::tag as, basic_block_t *parent,
                debug_info dbg_info = debug_info())
      : statement_t(REF_MAKE, parent, dbg_info),
	m_lhs(lhs), m_region(region), m_size(size), m_alloc_site(as) {
    this->m_live.add_def(m_lhs);
    this->m_live.add_use(m_region);
  }

  const variable_t &lhs() const { return m_lhs; }

  const variable_t &region() const { return m_region; }

  const crab::tag &alloc_site() const { return m_alloc_site;}

  const variable_or_constant_t &size() const { return m_size;}
  
  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_region, m_size, m_alloc_site, parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << m_lhs << " := "
      << "make_ref(" << m_region << ":" << m_region.get_type() << ","
      << m_size << ","
      << m_alloc_site << ")";
  }

private:
  variable_t m_lhs;
  variable_t m_region;
  variable_or_constant_t m_size;
  crab::tag m_alloc_site;
};


template <class BasicBlockLabel, class Number, class VariableName>
class remove_ref_stmt : public statement<BasicBlockLabel, Number, VariableName> {
  // remove_ref(region, ref)
  using this_type = remove_ref_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;

  remove_ref_stmt(variable_t region, variable_t ref,
		  basic_block_t *parent,
		  debug_info dbg_info = debug_info())
    : statement_t(REF_REMOVE, parent, dbg_info),
      m_region(region), m_ref(ref) {
    this->m_live.add_use(m_region);    
    this->m_live.add_use(m_ref);
  }

  const variable_t &region() const { return m_region; }

  const variable_t &ref() const { return m_ref; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_region, m_ref, parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << "remove_ref(" << m_region << ":" << m_region.get_type() << ","
      << m_ref << ":" << m_ref.get_type() << ")";
  }

private:
  variable_t m_region;  
  variable_t m_ref;
};
  
template <class BasicBlockLabel, class Number, class VariableName>
class load_from_ref_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  // lhs := load_from_ref(ref, region);
  using this_type = load_from_ref_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;

  load_from_ref_stmt(variable_t lhs, variable_t ref, variable_t region,
                     basic_block_t *parent, debug_info dbg_info = debug_info())
      : statement_t(REF_LOAD, parent, dbg_info), m_lhs(lhs), m_ref(ref),
        m_region(region) {
    this->m_live.add_def(m_lhs);
    this->m_live.add_use(m_region);
    this->m_live.add_use(m_ref);
  }

  const variable_t &lhs() const { return m_lhs; }

  const variable_t &ref() const { return m_ref; }

  const variable_t &region() const { return m_region; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_ref, m_region, parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << m_lhs << ":" << m_lhs.get_type() << " := "
      << "load_from_ref(" << m_region << ":" << m_region.get_type() << ","
      << m_ref << ":" << m_ref.get_type() << ")";
  }

private:
  variable_t m_lhs;
  variable_t m_ref;
  variable_t m_region;
};

template <class BasicBlockLabel, class Number, class VariableName>
class store_to_ref_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  // store_to_ref(ref, region, val);
  using this_type = store_to_ref_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using variable_or_constant_t = variable_or_constant<Number, VariableName>;

  store_to_ref_stmt(variable_t ref, variable_t region,
                    variable_or_constant_t val, basic_block_t *parent,
                    debug_info dbg_info = debug_info())
      : statement_t(REF_STORE, parent, dbg_info), m_ref(ref), m_region(region),
        m_val(val) {
    this->m_live.add_def(m_region);    
    this->m_live.add_use(m_region);
    this->m_live.add_use(m_ref);
    if (m_val.is_variable()) {
      this->m_live.add_use(m_val.get_variable());
    }
  }

  const variable_t &ref() const { return m_ref; }

  const variable_or_constant_t &val() const { return m_val; }

  const variable_t &region() const { return m_region; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_ref, m_region, m_val, parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << "store_to_ref(" << m_region << ":" << m_region.get_type() << ","
      << m_ref << ":" << m_ref.get_type() << "," << m_val << ":"
      << m_val.get_type() << ")";
  }

private:
  variable_t m_ref;
  variable_t m_region;
  variable_or_constant_t m_val;
};

template <class BasicBlockLabel, class Number, class VariableName>
class gep_ref_stmt : public statement<BasicBlockLabel, Number, VariableName> {
  // (lhs, lhs_region) := gep_ref(rhs, rhs_region, offset)
  using this_type = gep_ref_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using linear_expression_t = ikos::linear_expression<Number, VariableName>;
  gep_ref_stmt(variable_t lhs, variable_t lhs_region, variable_t rhs,
               variable_t rhs_region, linear_expression_t offset,
               basic_block_t *parent, debug_info dbg_info = debug_info())
      : statement_t(REF_GEP, parent, dbg_info), m_lhs(lhs),
        m_lhs_region(lhs_region), m_rhs(rhs), m_rhs_region(rhs_region),
        m_offset(offset) {
    // JN: not sure whether lhs_region should be a definition or an
    // use.
    // 
    // Although the abstract state of lhs_region might change, for
    // dead code/variable elimination purposes, this statement can be
    // removed if lhs is not used.
    this->m_live.add_def(m_lhs);    
    //this->m_live.add_def(m_lhs_region);
    this->m_live.add_use(m_lhs_region);
    this->m_live.add_use(m_rhs);
    this->m_live.add_use(m_rhs_region);
    for (auto const &v : m_offset.variables()) {
      this->m_live.add_use(v);
    }
  }

  const variable_t &lhs() const { return m_lhs; }

  const variable_t &rhs() const { return m_rhs; }

  const linear_expression_t &offset() const { return m_offset; }

  const variable_t &lhs_region() const { return m_lhs_region; }

  const variable_t &rhs_region() const { return m_rhs_region; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_lhs_region, m_rhs, m_rhs_region, m_offset,
                         parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << "(" << m_lhs_region << ":" << m_lhs_region.get_type() << "," << m_lhs
      << ":" << m_lhs.get_type() << ") := ";
    if (m_offset.equal(Number(0))) {
      o << "gep_ref(" << m_rhs_region << ":" << m_rhs_region.get_type() << ","
        << m_rhs << ":" << m_rhs.get_type() << ")";
    } else {
      o << "gep_ref(" << m_rhs_region << ":" << m_rhs_region.get_type() << ","
        << m_rhs << ":" << m_rhs.get_type() << " + " << m_offset << ")";
    }
    return;
  }

private:
  variable_t m_lhs;
  variable_t m_lhs_region;
  variable_t m_rhs;
  variable_t m_rhs_region;
  linear_expression_t m_offset;
};

template <class BasicBlockLabel, class Number, class VariableName>
class assume_ref_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = assume_ref_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using reference_constraint_t = reference_constraint<Number, VariableName>;

  assume_ref_stmt(reference_constraint_t cst, basic_block_t *parent)
      : statement_t(REF_ASSUME, parent), m_cst(cst) {
    if (!m_cst.is_tautology() && !m_cst.is_contradiction()) {
      if (m_cst.is_unary()) {
        this->m_live.add_use(m_cst.lhs());
      } else {
        this->m_live.add_use(m_cst.lhs());
        this->m_live.add_use(m_cst.rhs());
      }
    }
  }

  const reference_constraint_t &constraint() const { return m_cst; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_cst, parent);
  }

  virtual void write(crab_os &o) const override {
    o << "assume(" << m_cst << ")";
  }

private:
  reference_constraint_t m_cst;
};

template <class BasicBlockLabel, class Number, class VariableName>
class assert_ref_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = assert_ref_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using reference_constraint_t = reference_constraint<Number, VariableName>;

  assert_ref_stmt(reference_constraint_t cst, basic_block_t *parent,
                  debug_info dbg_info = debug_info())
      : statement_t(REF_ASSERT, parent, dbg_info), m_cst(cst) {
    if (!m_cst.is_tautology() && !m_cst.is_contradiction()) {
      if (m_cst.is_unary()) {
        this->m_live.add_use(m_cst.lhs());
      } else {
        this->m_live.add_use(m_cst.lhs());
        this->m_live.add_use(m_cst.rhs());
      }
    }
  }

  const reference_constraint_t &constraint() const { return m_cst; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_cst, parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << "assert(" << m_cst << ")";
    if (this->m_dbg_info.has_debug()) {
      o << "   /* " << this->m_dbg_info <<  " */";      
    }
  }

private:
  reference_constraint_t m_cst;
};

// ref_select lhs, b, op1, op2:
//    if b then lhs=op1 else lhs=op2
template <class BasicBlockLabel, class Number, class VariableName>
class ref_select_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = ref_select_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using variable_or_constant_t = variable_or_constant<Number, VariableName>;

  /* opX_rgn can be none if opX_ref is a null constant */
  ref_select_stmt(variable_t lhs_ref, variable_t lhs_rgn, variable_t cond,
                  variable_or_constant_t op1_ref,
                  boost::optional<variable_t> op1_rgn,
                  variable_or_constant_t op2_ref,
                  boost::optional<variable_t> op2_rgn, basic_block_t *parent)
      : statement_t(REF_SELECT, parent), m_lhs_ref(lhs_ref), m_lhs_rgn(lhs_rgn),
        m_cond(cond), m_op1_ref(op1_ref), m_op1_rgn(op1_rgn),
        m_op2_ref(op2_ref), m_op2_rgn(op2_rgn) {

    this->m_live.add_def(m_lhs_ref);
    this->m_live.add_def(m_lhs_rgn);
    this->m_live.add_use(m_cond);
    if (m_op1_ref.is_variable())
      this->m_live.add_use(m_op1_ref.get_variable());
    if (m_op2_ref.is_variable())
      this->m_live.add_use(m_op2_ref.get_variable());
    if (m_op1_rgn)
      this->m_live.add_use(*m_op1_rgn);
    if (m_op2_rgn)
      this->m_live.add_use(*m_op2_rgn);
  }

  const variable_t &lhs_ref() const { return m_lhs_ref; }

  const variable_t &lhs_rgn() const { return m_lhs_rgn; }

  const variable_t &cond() const { return m_cond; }

  const variable_or_constant_t &left_ref() const { return m_op1_ref; }

  // if left_ref() is a null constant then left_rgn() returns none
  boost::optional<variable_t> left_rgn() const { return m_op1_rgn; }

  const variable_or_constant_t &right_ref() const { return m_op2_ref; }

  // if right_ref() is a null constant then right_rgn() returns none
  boost::optional<variable_t> right_rgn() const { return m_op2_rgn; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs_ref, m_lhs_rgn, m_cond, m_op1_ref, m_op1_rgn,
                         m_op2_ref, m_op2_rgn, parent);
  }

  virtual void write(crab_os &o) const override {
    o << "(" << m_lhs_ref << "," << m_lhs_rgn << ")"
      << " = "
      << "ite(" << m_cond << ",";
    if (m_op1_rgn) {
      o << "(" << m_op1_ref << "," << *m_op1_rgn << ")";
    } else {
      o << m_op1_ref;
    }
    o << ",";
    if (m_op2_rgn) {
      o << "(" << m_op2_ref << "," << *m_op2_rgn << ")";
    } else {
      o << m_op2_ref;
    }
    o << ")";
  }

private:
  variable_t m_lhs_ref; // pre: reference type
  variable_t m_lhs_rgn;
  variable_t m_cond;                // pre: boolean type
  variable_or_constant_t m_op1_ref; // pre: reference type
  boost::optional<variable_t> m_op1_rgn;
  variable_or_constant_t m_op2_ref; // pre: reference type
  boost::optional<variable_t> m_op2_rgn;
};

template <class BasicBlockLabel, class Number, class VariableName>
class ref_to_int_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  // inv_var := ref_to_int(region, ref_var)
  using this_type = ref_to_int_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;

  ref_to_int_stmt(variable_t region, variable_t ref_var, variable_t int_var,
                  basic_block_t *parent, debug_info dbg_info = debug_info())
      : statement_t(REF_TO_INT, parent, dbg_info), m_region(region),
        m_ref_var(ref_var), m_int_var(int_var) {
    this->m_live.add_use(m_ref_var);
    this->m_live.add_def(m_int_var);
  }

  const variable_t &region() const { return m_region; }

  const variable_t &ref_var() const { return m_ref_var; }

  const variable_t &int_var() const { return m_int_var; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_region, m_ref_var, m_int_var, parent,
                         this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << m_int_var << ":" << m_int_var.get_type() << " := "
      << "ref_to_int(" << m_region << ":" << m_region.get_type() << ","
      << m_ref_var << ":" << m_ref_var.get_type() << ")";
  }

private:
  variable_t m_region;
  variable_t m_ref_var;
  variable_t m_int_var;
};

template <class BasicBlockLabel, class Number, class VariableName>
class int_to_ref_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  // ref_var := int_to_ref(int_var, region)
  using this_type = int_to_ref_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;

  int_to_ref_stmt(variable_t int_var, variable_t region, variable_t ref_var,
                  basic_block_t *parent, debug_info dbg_info = debug_info())
      : statement_t(INT_TO_REF, parent, dbg_info), m_int_var(int_var),
        m_region(region), m_ref_var(ref_var) {
    this->m_live.add_use(m_int_var);
    this->m_live.add_def(m_ref_var);
  }

  const variable_t &int_var() const { return m_int_var; }

  const variable_t &region() const { return m_region; }

  const variable_t &ref_var() const { return m_ref_var; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_int_var, m_region, m_ref_var, parent,
                         this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << m_ref_var << ":" << m_ref_var.get_type() << " := "
      << "int_to_ref(" << m_region << ":" << m_region.get_type() << ","
      << m_int_var << ":" << m_int_var.get_type() << ")";
  }

private:
  variable_t m_int_var;
  variable_t m_region;
  variable_t m_ref_var;
};

/*
  Function calls
*/

template <class BasicBlockLabel, class Number, class VariableName>
class callsite_stmt : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = callsite_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using type_t = typename variable_t::type_t;

  callsite_stmt(std::string func_name, const std::vector<variable_t> &args,
                basic_block_t *parent)
      : statement_t(CALLSITE, parent), m_func_name(func_name) {

    std::copy(args.begin(), args.end(), std::back_inserter(m_args));
    for (auto arg : m_args) {
      this->m_live.add_use(arg);
    }
  }

  callsite_stmt(std::string func_name, const std::vector<variable_t> &lhs,
                const std::vector<variable_t> &args, basic_block_t *parent)
      : statement_t(CALLSITE, parent), m_func_name(func_name) {

    std::copy(args.begin(), args.end(), std::back_inserter(m_args));
    for (auto arg : m_args) {
      this->m_live.add_use(arg);
    }

    std::copy(lhs.begin(), lhs.end(), std::back_inserter(m_lhs));
    for (auto arg : m_lhs) {
      this->m_live.add_def(arg);
    }
  }

  const std::string &get_func_name() const { return m_func_name; }
  
  const std::vector<variable_t> &get_lhs() const { return m_lhs; }

  unsigned get_num_lhs() const { return m_lhs.size(); }

  const variable_t &get_lhs_name(unsigned idx) const {
    if (idx >= m_lhs.size()) {
      CRAB_ERROR("Out-of-bound access to call site parameter");
    }
    return m_lhs[idx];
  }

  type_t get_lhs_type(unsigned idx) const {
    if (idx >= m_lhs.size()) {
      CRAB_ERROR("Out-of-bound access to call site parameter");
    }
    return m_lhs[idx].get_type();
  }

  const std::vector<variable_t> &get_args() const { return m_args; }

  unsigned get_num_args() const { return m_args.size(); }

  const variable_t &get_arg_name(unsigned idx) const {
    if (idx >= m_args.size()) {
      CRAB_ERROR("Out-of-bound access to call site parameter");
    }
    return m_args[idx];
  }

  type_t get_arg_type(unsigned idx) const {
    if (idx >= m_args.size()) {
      CRAB_ERROR("Out-of-bound access to call site parameter");
    }
    return m_args[idx].get_type();
  }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_func_name, m_lhs, m_args, parent);
  }

  virtual void write(crab_os &o) const override {
    if (m_lhs.empty()) {
      // do nothing
    } else if (m_lhs.size() == 1) {
      o << (*m_lhs.begin()) << ":" << (*m_lhs.begin()).get_type() << " =";
    } else {
      o << "(";
      for (auto It = m_lhs.begin(), Et = m_lhs.end(); It != Et;) {
        o << (*It) << ":" << (*It).get_type();
        ++It;
        if (It != Et)
          o << ",";
      }
      o << ")=";
    }
    o << " call " << m_func_name << "(";
    for (auto It = m_args.begin(), Et = m_args.end(); It != Et;) {
      o << *It << ":" << (*It).get_type();
      ++It;
      if (It != Et)
        o << ",";
    }
    o << ")";
  }

private:
  std::string m_func_name;
  std::vector<variable_t> m_lhs;
  std::vector<variable_t> m_args;
};


/* An intrinsic function is an "internal" function with semantics
   defined by the Crab domains */
template <class BasicBlockLabel, class Number, class VariableName>
class intrinsic_stmt : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = intrinsic_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using variable_or_constant_t = variable_or_constant<Number, VariableName>;
  using type_t = typename variable_t::type_t;

  intrinsic_stmt(std::string intrinsic_name,
                 const std::vector<variable_or_constant_t> &args,
		 basic_block_t *parent,
		 debug_info dbg_info = debug_info())
    : statement_t(CRAB_INTRINSIC, parent, dbg_info), m_intrinsic_name(intrinsic_name) {

    std::copy(args.begin(), args.end(), std::back_inserter(m_args));
    for (auto arg : m_args) {
      if (arg.is_variable()) {
	this->m_live.add_use(arg.get_variable());
      }
    }
  }

  intrinsic_stmt(std::string intrinsic_name,
		 const std::vector<variable_t> &lhs,
                 const std::vector<variable_or_constant_t> &args,
		 basic_block_t *parent,
		 debug_info dbg_info = debug_info())
    : statement_t(CRAB_INTRINSIC, parent, dbg_info), m_intrinsic_name(intrinsic_name) {

    std::copy(args.begin(), args.end(), std::back_inserter(m_args));
    for (auto arg : m_args) {
      if (arg.is_variable()) {
	this->m_live.add_use(arg.get_variable());
      }
    }

    std::copy(lhs.begin(), lhs.end(), std::back_inserter(m_lhs));
    for (auto arg : m_lhs) {
      this->m_live.add_def(arg);
    }
  }
  
  const std::vector<variable_t> &get_lhs() const { return m_lhs; }

  const std::string &get_intrinsic_name() const { return m_intrinsic_name; }

  const std::vector<variable_or_constant_t> &get_args() const { return m_args; }

  unsigned get_num_args() const { return m_args.size(); }

  const variable_t &get_arg_name(unsigned idx) const {
    if (idx >= m_args.size())
      CRAB_ERROR("Out-of-bound access to intrinsic parameter");

    return m_args[idx];
  }

  type_t get_arg_type(unsigned idx) const {
    if (idx >= m_args.size())
      CRAB_ERROR("Out-of-bound access to intrinsic parameter");

    return m_args[idx].get_type();
  }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_intrinsic_name, m_lhs, m_args, parent);
  }

  virtual void write(crab_os &o) const override {
    if (m_lhs.empty()) {
      // do nothing
    } else if (m_lhs.size() == 1) {
      o << (*m_lhs.begin()) << " = ";
    } else {
      o << "(";
      for (auto It = m_lhs.begin(), Et = m_lhs.end(); It != Et;) {
        o << (*It);
        ++It;
        if (It != Et)
          o << ",";
      }
      o << ")= ";
    }
    o << "crab_intrinsic(" << m_intrinsic_name;
    if (!m_args.empty()) {
      o << ",";
    }
    for (auto It = m_args.begin(), Et = m_args.end(); It != Et;) {
      o << *It << ":" << (*It).get_type();
      ++It;
      if (It != Et)
        o << ",";
    }
    o << ")";
    if (this->m_dbg_info.has_debug()) {
      o << "   /* " << this->m_dbg_info <<  " */";
    }
  }

private:
  std::string m_intrinsic_name;
  std::vector<variable_t> m_lhs;
  std::vector<variable_or_constant_t> m_args;
};

/*
   Boolean statements
*/

template <class BasicBlockLabel, class Number, class VariableName>
class bool_assign_cst
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = bool_assign_cst<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;
  using linear_constraint_t = ikos::linear_constraint<Number, VariableName>;
  using reference_constraint_t = reference_constraint<Number, VariableName>;
  using linear_or_reference_constraint_t =
      boost::variant<linear_constraint_t, reference_constraint_t>;

  bool_assign_cst(variable_t lhs, linear_constraint_t rhs,
                  basic_block_t *parent)
      : statement_t(BOOL_ASSIGN_CST, parent), m_is_rhs_linear_constraint(true),
        m_lhs(lhs), m_rhs(rhs) {
    this->m_live.add_def(m_lhs);
    for (auto const &v : rhs_as_linear_constraint().variables())
      this->m_live.add_use(v);
  }

  bool_assign_cst(variable_t lhs, reference_constraint_t rhs,
                  basic_block_t *parent)
      : statement_t(BOOL_ASSIGN_CST, parent), m_is_rhs_linear_constraint(false),
        m_lhs(lhs), m_rhs(rhs) {
    this->m_live.add_def(m_lhs);
    for (auto const &v : rhs_as_reference_constraint().variables())
      this->m_live.add_use(v);
  }

  const variable_t &lhs() const { return m_lhs; }

  const linear_constraint_t &rhs_as_linear_constraint() const {
    if (m_is_rhs_linear_constraint) {
      return boost::get<linear_constraint_t>(m_rhs);
    } else {
      CRAB_ERROR("rhs of bool_assign_cst is not a linear constraint");
    }
  }

  const reference_constraint_t &rhs_as_reference_constraint() const {
    if (!m_is_rhs_linear_constraint) {
      return boost::get<reference_constraint_t>(m_rhs);
    } else {
      CRAB_ERROR("rhs of bool_assign_cst is not a reference constraint");
    }
  }

  bool is_rhs_linear_constraint() const { return m_is_rhs_linear_constraint; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    if (is_rhs_linear_constraint()) {
      return new this_type(m_lhs, rhs_as_linear_constraint(), parent);
    } else {
      return new this_type(m_lhs, rhs_as_reference_constraint(), parent);
    }
  }

  virtual void write(crab_os &o) const override {
    o << m_lhs << " = ";
    if (is_rhs_linear_constraint()) {
      auto const &rhs = rhs_as_linear_constraint();
      if (rhs.is_tautology()) {
        o << "true";
      } else if (rhs.is_contradiction()) {
        o << "false";
      } else {
        o << "(" << rhs << ")";
      }
    } else {
      auto const &rhs = rhs_as_reference_constraint();
      if (rhs.is_tautology()) {
        o << "true";
      } else if (rhs.is_contradiction()) {
        o << "false";
      } else {
        o << "(" << rhs << ")";
      }
    }
  }

private:
  bool m_is_rhs_linear_constraint;
  variable_t m_lhs; // pre: boolean type
  linear_or_reference_constraint_t m_rhs;
};

template <class BasicBlockLabel, class Number, class VariableName>
class bool_assign_var
    : public statement<BasicBlockLabel, Number, VariableName> {
  // Note that this can be simulated with bool_binary_op (e.g.,
  // b1 := b2 ----> b1 := b2 or false). However, we create a
  // special statement to assign a variable to another because it
  // is a very common operation.

  using this_type = bool_assign_var<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;

  bool_assign_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                  basic_block_t *parent)
      : statement_t(BOOL_ASSIGN_VAR, parent), m_lhs(lhs), m_rhs(rhs),
        m_is_rhs_negated(is_not_rhs) {
    this->m_live.add_def(m_lhs);
    this->m_live.add_use(m_rhs);
  }

  const variable_t &lhs() const { return m_lhs; }

  const variable_t &rhs() const { return m_rhs; }

  bool is_rhs_negated() const { return m_is_rhs_negated; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_rhs, m_is_rhs_negated, parent);
  }

  virtual void write(crab_os &o) const override {
    o << m_lhs << " = ";
    if (is_rhs_negated()) {
      o << "not(" << m_rhs << ")";
    } else {
      o << m_rhs;
    }
  }

private:
  variable_t m_lhs; // pre: boolean type
  variable_t m_rhs; // pre: boolean type

  // if m_is_rhs_negated is true then lhs := not(rhs).
  bool m_is_rhs_negated;
};

template <class BasicBlockLabel, class Number, class VariableName>
class bool_binary_op : public statement<BasicBlockLabel, Number, VariableName> {
  // b1:= b2 and b3
  // b1:= b2 or b3
  // b1:= b2 xor b3
  using this_type = bool_binary_op<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;

  bool_binary_op(variable_t lhs, bool_binary_operation_t op, variable_t op1,
                 variable_t op2, basic_block_t *parent,
                 debug_info dbg_info = debug_info())
      : statement_t(BOOL_BIN_OP, parent, dbg_info), m_lhs(lhs), m_op(op),
        m_op1(op1), m_op2(op2) {
    this->m_live.add_def(m_lhs);
    this->m_live.add_use(m_op1);
    this->m_live.add_use(m_op2);
  }

  const variable_t &lhs() const { return m_lhs; }

  bool_binary_operation_t op() const { return m_op; }

  const variable_t &left() const { return m_op1; }

  const variable_t &right() const { return m_op2; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_op, m_op1, m_op2, parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << m_lhs << " = " << m_op1 << m_op << m_op2;
  }

private:
  variable_t m_lhs; // pre: boolean type
  bool_binary_operation_t m_op;
  variable_t m_op1; // pre: boolean type
  variable_t m_op2; // pre: boolean type
};

template <class BasicBlockLabel, class Number, class VariableName>
class bool_assume_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = bool_assume_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;

  bool_assume_stmt(variable_t v, bool is_negated, basic_block_t *parent)
      : statement_t(BOOL_ASSUME, parent), m_var(v), m_is_negated(is_negated) {
    this->m_live.add_use(v);
  }

  const variable_t &cond() const { return m_var; }

  bool is_negated() const { return m_is_negated; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_var, m_is_negated, parent);
  }

  virtual void write(crab_os &o) const override {
    if (m_is_negated) {
      o << "assume(not(" << m_var << "))";
    } else {
      o << "assume(" << m_var << ")";
    }
  }

private:
  variable_t m_var; // pre: boolean type
  bool m_is_negated;
};

// select b1, b2, b3, b4:
//    if b2 then b1=b3 else b1=b4
template <class BasicBlockLabel, class Number, class VariableName>
class bool_select_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = bool_select_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;

  bool_select_stmt(variable_t lhs, variable_t cond, variable_t b1,
                   variable_t b2, basic_block_t *parent)
      : statement_t(BOOL_SELECT, parent), m_lhs(lhs), m_cond(cond), m_b1(b1),
        m_b2(b2) {
    this->m_live.add_def(m_lhs);
    this->m_live.add_use(m_cond);
    this->m_live.add_use(m_b1);
    this->m_live.add_use(m_b2);
  }

  const variable_t &lhs() const { return m_lhs; }

  const variable_t &cond() const { return m_cond; }

  const variable_t &left() const { return m_b1; }

  const variable_t &right() const { return m_b2; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_lhs, m_cond, m_b1, m_b2, parent);
  }

  virtual void write(crab_os &o) const override {
    o << m_lhs << " = "
      << "ite(" << m_cond << "," << m_b1 << "," << m_b2 << ")";
  }

private:
  variable_t m_lhs;  // pre: boolean type
  variable_t m_cond; // pre: boolean type
  variable_t m_b1;   // pre: boolean type
  variable_t m_b2;   // pre: boolean type
};

template <class BasicBlockLabel, class Number, class VariableName>
class bool_assert_stmt
    : public statement<BasicBlockLabel, Number, VariableName> {
  using this_type = bool_assert_stmt<BasicBlockLabel, Number, VariableName>;

public:
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = typename statement_t::basic_block_t;
  using variable_t = variable<Number, VariableName>;

  bool_assert_stmt(variable_t v, basic_block_t *parent,
                   debug_info dbg_info = debug_info())
      : statement_t(BOOL_ASSERT, parent, dbg_info), m_var(v) {
    this->m_live.add_use(v);
  }

  const variable_t &cond() const { return m_var; }

  virtual void
  accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) override {
    v->visit(*this);
  }

  virtual statement_t *clone(basic_block_t *parent) const override {
    return new this_type(m_var, parent, this->m_dbg_info);
  }

  virtual void write(crab_os &o) const override {
    o << "assert(" << m_var << ")";
    if (this->m_dbg_info.has_debug()) {
      o << "   /* " << this->m_dbg_info <<  " */";
    }
  }

private:
  variable_t m_var; // pre: boolean type
};

template <class BasicBlockLabel, class VariableName, class Number>
class cfg;

template <class BasicBlockLabel, class VariableName, class Number>
class basic_block {
  friend class cfg<BasicBlockLabel, VariableName, Number>;

public:
  using number_t = Number;
  using varname_t = VariableName;
  using basic_block_label_t = BasicBlockLabel;

  // helper types to build statements
  using variable_t = variable<Number, VariableName>;
  using variable_or_constant_t = variable_or_constant<Number, VariableName>;
  using lin_exp_t = ikos::linear_expression<Number, VariableName>;
  using lin_cst_t = ikos::linear_constraint<Number, VariableName>;
  using ref_cst_t = reference_constraint<Number, VariableName>;
  using statement_t = statement<BasicBlockLabel, Number, VariableName>;
  using basic_block_t = basic_block<BasicBlockLabel, VariableName, Number>;
  using interval_t = ikos::interval<Number>;

private:
  using bb_id_set_t = std::vector<BasicBlockLabel>;
  using stmt_list_t = std::vector<statement_t *>;

public:
  // -- iterators

  using succ_iterator = typename bb_id_set_t::iterator;
  using const_succ_iterator = typename bb_id_set_t::const_iterator;
  using pred_iterator = succ_iterator;
  using const_pred_iterator = const_succ_iterator;
  using iterator = boost::indirect_iterator<typename stmt_list_t::iterator>;
  using const_iterator =
      boost::indirect_iterator<typename stmt_list_t::const_iterator>;
  using reverse_iterator =
      boost::indirect_iterator<typename stmt_list_t::reverse_iterator>;
  using const_reverse_iterator =
      boost::indirect_iterator<typename stmt_list_t::const_reverse_iterator>;
  using live_domain_t = ikos::discrete_domain<variable_t>;

  // -- statements

  using havoc_t = havoc_stmt<BasicBlockLabel, Number, VariableName>;
  using unreach_t = unreachable_stmt<BasicBlockLabel, Number, VariableName>;
  // Numerical
  using bin_op_t = binary_op<BasicBlockLabel, Number, VariableName>;
  using assign_t = assignment<BasicBlockLabel, Number, VariableName>;
  using assume_t = assume_stmt<BasicBlockLabel, Number, VariableName>;
  using select_t = select_stmt<BasicBlockLabel, Number, VariableName>;
  using assert_t = assert_stmt<BasicBlockLabel, Number, VariableName>;
  using int_cast_t = int_cast_stmt<BasicBlockLabel, Number, VariableName>;
  // Functions
  using callsite_t = callsite_stmt<BasicBlockLabel, Number, VariableName>;
  using intrinsic_t = intrinsic_stmt<BasicBlockLabel, Number, VariableName>;
  // Arrays
  using arr_init_t = array_init_stmt<BasicBlockLabel, Number, VariableName>;
  using arr_store_t = array_store_stmt<BasicBlockLabel, Number, VariableName>;
  using arr_load_t = array_load_stmt<BasicBlockLabel, Number, VariableName>;
  using arr_assign_t = array_assign_stmt<BasicBlockLabel, Number, VariableName>;
  // Regions and references
  using region_init_t = region_init_stmt<BasicBlockLabel, Number, VariableName>;
  using region_copy_t = region_copy_stmt<BasicBlockLabel, Number, VariableName>;
  using region_cast_t = region_cast_stmt<BasicBlockLabel, Number, VariableName>;  
  using make_ref_t = make_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using remove_ref_t = remove_ref_stmt<BasicBlockLabel, Number, VariableName>;  
  using load_from_ref_t =
      load_from_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using store_to_ref_t =
      store_to_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using gep_ref_t = gep_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using assume_ref_t = assume_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using assert_ref_t = assert_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using ref_select_t = ref_select_stmt<BasicBlockLabel, Number, VariableName>;
  using int_to_ref_t = int_to_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using ref_to_int_t = ref_to_int_stmt<BasicBlockLabel, Number, VariableName>;
  // Boolean
  using bool_bin_op_t = bool_binary_op<BasicBlockLabel, Number, VariableName>;
  using bool_assign_cst_t =
      bool_assign_cst<BasicBlockLabel, Number, VariableName>;
  using bool_assign_var_t =
      bool_assign_var<BasicBlockLabel, Number, VariableName>;
  using bool_assume_t = bool_assume_stmt<BasicBlockLabel, Number, VariableName>;
  using bool_select_t = bool_select_stmt<BasicBlockLabel, Number, VariableName>;
  using bool_assert_t = bool_assert_stmt<BasicBlockLabel, Number, VariableName>;

private:
  BasicBlockLabel m_bb_id;
  stmt_list_t m_stmts;
  bb_id_set_t m_prev, m_next;
  // Ideally it should be size_t to indicate any position within the
  // block. For now, we only allow to insert either at front or at
  // the back (default). Note that if insertions at the front are
  // very common we should replace stmt_list_t from a vector to a
  // deque.
  bool m_insert_point_at_front;
  // set of used/def variables
  live_domain_t m_live;

  void insert_adjacent(bb_id_set_t &c, BasicBlockLabel e) {
    if (std::find(c.begin(), c.end(), e) == c.end()) {
      c.push_back(e);
    }
  }

  void remove_adjacent(bb_id_set_t &c, BasicBlockLabel e) {
    if (std::find(c.begin(), c.end(), e) != c.end()) {
      c.erase(std::remove(c.begin(), c.end(), e), c.end());
    }
  }

  basic_block(BasicBlockLabel bb_id)
      : m_bb_id(bb_id), m_insert_point_at_front(false),
        m_live(live_domain_t::bottom()) {}

  static basic_block_t *create(BasicBlockLabel bb_id) {
    return new basic_block_t(bb_id);
  }

  void update_uses_and_defs(const statement_t &s) {
    auto const &ls = s.get_live();
    for (auto &v : boost::make_iterator_range(ls.uses_begin(), ls.uses_end())) {
      m_live += v;
    }
    for (auto &v : boost::make_iterator_range(ls.defs_begin(), ls.defs_end())) {
      m_live += v;
    }
  }

  const statement_t *insert(statement_t *stmt) {
    if (m_insert_point_at_front) {
      m_stmts.insert(m_stmts.begin(), stmt);
      m_insert_point_at_front = false;
    } else {
      m_stmts.push_back(stmt);
    }
    update_uses_and_defs(*stmt);
    return stmt;
  }

public:
  basic_block(const basic_block_t &o) = delete;

  basic_block_t &operator=(const basic_block_t &o) = delete;

  // The basic block owns the statements
  ~basic_block() {
    for (unsigned i = 0, e = m_stmts.size(); i < e; ++i) {
      delete m_stmts[i];
    }
  }

  // it will be set to false after the first insertion
  void set_insert_point_front() { m_insert_point_at_front = true; }

  basic_block_t *clone() const {
    // The basic block labels (i.e., identifiers) are not cloned.

    basic_block_t *b = new basic_block_t(label());
    for (auto &s : boost::make_iterator_range(begin(), end())) {
      b->m_stmts.push_back(s.clone(b));
    }

    for (auto id : boost::make_iterator_range(prev_blocks())) {
      b->m_prev.push_back(id);
    }
    for (auto id : boost::make_iterator_range(next_blocks())) {
      b->m_next.push_back(id);
    }

    b->m_live = m_live;
    return b;
  }

  const BasicBlockLabel &label() const { return m_bb_id; }

  std::string name() const {
    return basic_block_traits<basic_block_t>::to_string(m_bb_id);
  }

  iterator begin() { return boost::make_indirect_iterator(m_stmts.begin()); }
  iterator end() { return boost::make_indirect_iterator(m_stmts.end()); }
  const_iterator begin() const {
    return boost::make_indirect_iterator(m_stmts.begin());
  }
  const_iterator end() const {
    return boost::make_indirect_iterator(m_stmts.end());
  }

  reverse_iterator rbegin() {
    return boost::make_indirect_iterator(m_stmts.rbegin());
  }
  reverse_iterator rend() {
    return boost::make_indirect_iterator(m_stmts.rend());
  }
  const_reverse_iterator rbegin() const {
    return boost::make_indirect_iterator(m_stmts.rbegin());
  }
  const_reverse_iterator rend() const {
    return boost::make_indirect_iterator(m_stmts.rend());
  }

  // Return a mutable reference to allow modifying the CFG
  statement_t& operator[](unsigned int i) {
    if (i >= size()) {
      CRAB_ERROR("Out-of-bound access in basic_block::operator[]");
    }
    return *m_stmts[i];
  }
  
  size_t size() const { return m_stmts.size(); }

  live_domain_t &live() { return m_live; }

  const live_domain_t &live() const { return m_live; }

  // Collect the set of uses and definitions of the basic block
  void update_uses_and_defs() {
    for (const statement_t *s : m_stmts) {
      update_uses_and_defs(*s);
    }
  }

  void accept(statement_visitor<BasicBlockLabel, Number, VariableName> *v) {
    v->visit(*this);
  }

  std::pair<succ_iterator, succ_iterator> next_blocks() {
    return std::make_pair(m_next.begin(), m_next.end());
  }

  std::pair<pred_iterator, pred_iterator> prev_blocks() {
    return std::make_pair(m_prev.begin(), m_prev.end());
  }

  std::pair<const_succ_iterator, const_succ_iterator> next_blocks() const {
    return std::make_pair(m_next.begin(), m_next.end());
  }

  std::pair<const_pred_iterator, const_pred_iterator> prev_blocks() const {
    return std::make_pair(m_prev.begin(), m_prev.end());
  }

  size_t in_degree() const {
    return m_prev.size();
  }

  size_t out_degree() const {
    return m_next.size();
  }

  // Add a cfg edge from *this to b  
  void operator>>(basic_block_t &b) {
    add_succ(b);
  }

  // Add a cfg edge from *this to b
  void add_succ(basic_block_t &b) {
    insert_adjacent(m_next, b.m_bb_id);
    insert_adjacent(b.m_prev, m_bb_id);
  }

  // Remove a cfg edge from *this to b
  void operator-=(basic_block_t &b) {
    remove_adjacent(m_next, b.m_bb_id);
    remove_adjacent(b.m_prev, m_bb_id);
  }

  // insert all statements of other at the front
  void copy_front(const basic_block_t &other) {
    std::vector<statement_t *> cloned_stmts;
    cloned_stmts.reserve(other.size());
    std::transform(other.m_stmts.begin(), other.m_stmts.end(),
                   std::back_inserter(cloned_stmts),
                   [this](const statement_t *s) { return s->clone(this); });

    m_stmts.insert(m_stmts.begin(), cloned_stmts.begin(), cloned_stmts.end());
    m_live = m_live | other.m_live;
  }

  // insert all statements of other at the back
  void copy_back(const basic_block_t &other) {
    std::vector<statement_t *> cloned_stmts;
    cloned_stmts.reserve(other.size());
    std::transform(other.m_stmts.begin(), other.m_stmts.end(),
                   std::back_inserter(cloned_stmts),
                   [this](const statement_t *s) { return s->clone(this); });

    m_stmts.insert(m_stmts.end(), cloned_stmts.begin(), cloned_stmts.end());
    m_live = m_live | other.m_live;
  }

  // insert all statements of other at the back
  void move_back(basic_block_t &other) {
    m_stmts.reserve(m_stmts.size() + other.m_stmts.size());
    std::move(other.m_stmts.begin(), other.m_stmts.end(),
              std::back_inserter(m_stmts));
  }

  // Remove s (and free) from this
  void remove(const statement_t *s, bool must_update_uses_and_defs = true) {
    // remove statement using the remove-erase idiom
    m_stmts.erase(
        std::remove_if(m_stmts.begin(), m_stmts.end(),
                       [s](const statement_t *o) { return (o == s); }),
        m_stmts.end());

    if (must_update_uses_and_defs) {
      update_uses_and_defs();
    }
    delete s;
  }

  // Pre: old is a statement at position i in this
  //      new is not part of any basic block.
  // Post: old is deleted (and freed) and new is at position i in
  //       the basic block.
  void replace(statement_t *old_s, statement_t *new_s) {
    std::replace(m_stmts.begin(), m_stmts.end(), old_s, new_s);
    delete old_s;
  }

  void write(crab_os &o) const {
    o << basic_block_traits<basic_block_t>::to_string(m_bb_id) << ":\n";
    for (auto const &s : *this) {
      o << "  " << s << ";\n";
    }
    std::pair<const_succ_iterator, const_succ_iterator> p = next_blocks();
    const_succ_iterator it = p.first;
    const_succ_iterator et = p.second;
    if (it != et) {
      o << "  "
        << "goto ";
      for (; it != et;) {
        o << basic_block_traits<basic_block_t>::to_string(*it);
        ++it;
        if (it == et) {
          o << ";";
        } else {
          o << ",";
        }
      }
    }
    o << "\n";
  }

  // for gdb
  void dump() const { write(crab::errs()); }

  /**********************************************************************/
  /*                  Begin BUILD CFG STATEMENTS                        */
  /**********************************************************************/
  const statement_t *add(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_ADD, op1, op2, this));
  }

  const statement_t *add(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_ADD, op1, op2, this));
  }

  const statement_t *sub(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_SUB, op1, op2, this));
  }

  const statement_t *sub(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_SUB, op1, op2, this));
  }

  const statement_t *mul(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_MUL, op1, op2, this));
  }

  const statement_t *mul(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_MUL, op1, op2, this));
  }

  // signed division
  const statement_t *div(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_SDIV, op1, op2, this));
  }

  const statement_t *div(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_SDIV, op1, op2, this));
  }

  // unsigned division
  const statement_t *udiv(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_UDIV, op1, op2, this));
  }

  const statement_t *udiv(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_UDIV, op1, op2, this));
  }

  // signed rem
  const statement_t *rem(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_SREM, op1, op2, this));
  }

  const statement_t *rem(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_SREM, op1, op2, this));
  }

  // unsigned rem
  const statement_t *urem(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_UREM, op1, op2, this));
  }

  const statement_t *urem(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_UREM, op1, op2, this));
  }

  const statement_t *bitwise_and(variable_t lhs, variable_t op1,
                                 variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_AND, op1, op2, this));
  }

  const statement_t *bitwise_and(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_AND, op1, op2, this));
  }

  const statement_t *bitwise_or(variable_t lhs, variable_t op1,
                                variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_OR, op1, op2, this));
  }

  const statement_t *bitwise_or(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_OR, op1, op2, this));
  }

  const statement_t *bitwise_xor(variable_t lhs, variable_t op1,
                                 variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_XOR, op1, op2, this));
  }

  const statement_t *bitwise_xor(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_XOR, op1, op2, this));
  }

  const statement_t *shl(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_SHL, op1, op2, this));
  }

  const statement_t *shl(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_SHL, op1, op2, this));
  }

  const statement_t *lshr(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_LSHR, op1, op2, this));
  }

  const statement_t *lshr(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_LSHR, op1, op2, this));
  }

  const statement_t *ashr(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bin_op_t(lhs, BINOP_ASHR, op1, op2, this));
  }

  const statement_t *ashr(variable_t lhs, variable_t op1, number_t op2) {
    return insert(new bin_op_t(lhs, BINOP_ASHR, op1, op2, this));
  }

  // numerical assignment
  const statement_t *assign(variable_t lhs, lin_exp_t rhs) {
    return insert(new assign_t(lhs, rhs, this));
  }

  const statement_t *assume(lin_cst_t cst) {
    return insert(new assume_t(cst, this));
  }

  const statement_t *havoc(variable_t lhs, std::string comment = "") {
    return insert(new havoc_t(lhs, comment, this));
  }

  const statement_t *unreachable() { return insert(new unreach_t(this)); }

  const statement_t *select(variable_t lhs, variable_t v, lin_exp_t e1,
                            lin_exp_t e2) {
    lin_cst_t cond = (v >= number_t(1));
    return insert(new select_t(lhs, cond, e1, e2, this));
  }

  const statement_t *select(variable_t lhs, lin_cst_t cond, lin_exp_t e1,
                            lin_exp_t e2) {
    return insert(new select_t(lhs, cond, e1, e2, this));
  }

  const statement_t *assertion(lin_cst_t cst, debug_info di = debug_info()) {
    return insert(new assert_t(cst, this, di));
  }

  const statement_t *truncate(variable_t src, variable_t dst) {
    return insert(new int_cast_t(CAST_TRUNC, src, dst, this));
  }

  const statement_t *sext(variable_t src, variable_t dst) {
    return insert(new int_cast_t(CAST_SEXT, src, dst, this));
  }

  const statement_t *zext(variable_t src, variable_t dst) {
    return insert(new int_cast_t(CAST_ZEXT, src, dst, this));
  }

  const statement_t *callsite(std::string func,
                              const std::vector<variable_t> &lhs,
                              const std::vector<variable_t> &args) {
    return insert(new callsite_t(func, lhs, args, this));
  }

  const statement_t *intrinsic(std::string name,
                               const std::vector<variable_t> &lhs,
                               const std::vector<variable_or_constant_t> &args,
			       debug_info di = debug_info()) {
    return insert(new intrinsic_t(name, lhs, args, this, di));
  }

  const statement_t *array_init(variable_t a, lin_exp_t lb_idx,
                                lin_exp_t ub_idx, lin_exp_t v,
                                lin_exp_t elem_size) {
    return insert(new arr_init_t(a, elem_size, lb_idx, ub_idx, v, this));
  }

  const statement_t *array_store(variable_t arr, lin_exp_t idx, lin_exp_t v,
                                 lin_exp_t elem_size) {
    return insert(new arr_store_t(arr, elem_size, idx, idx, v, false, this));
  }

  const statement_t *array_store(variable_t arr, lin_exp_t idx, lin_exp_t v,
                                 lin_exp_t elem_size, bool is_strong_update) {
    return insert(
        new arr_store_t(arr, elem_size, idx, idx, v, is_strong_update, this));
  }

  const statement_t *array_store_range(variable_t arr, lin_exp_t lb_idx,
                                       lin_exp_t ub_idx, lin_exp_t v,
                                       lin_exp_t elem_size) {
    return insert(
        new arr_store_t(arr, elem_size, lb_idx, ub_idx, v, false, this));
  }

  const statement_t *array_load(variable_t lhs, variable_t arr, lin_exp_t idx,
                                lin_exp_t elem_size) {
    return insert(new arr_load_t(lhs, arr, elem_size, idx, this));
  }

  const statement_t *array_assign(variable_t lhs, variable_t rhs) {
    return insert(new arr_assign_t(lhs, rhs, this));
  }

  const statement_t *region_init(variable_t region) {
    return insert(new region_init_t(region, this));
  }

  const statement_t *region_copy(variable_t lhs, variable_t rhs) {
    return insert(new region_copy_t(lhs, rhs, this));
  }

  const statement_t *region_cast(variable_t src, variable_t dst) {
    return insert(new region_cast_t(src, dst, this));
  }
  
  const statement_t *make_ref(variable_t lhs_ref, variable_t region,
			      variable_or_constant_t size, crab::tag as) {
    return insert(new make_ref_t(lhs_ref, region, size, as, this));
  }

  const statement_t *remove_ref(variable_t region, variable_t ref) {
    return insert(new remove_ref_t(region, ref, this));
  }
  
  const statement_t *load_from_ref(variable_t lhs, variable_t ref,
                                   variable_t region) {
    return insert(new load_from_ref_t(lhs, ref, region, this));
  }

  const statement_t *store_to_ref(variable_t ref, variable_t region,
                                  variable_or_constant_t val) {
    return insert(new store_to_ref_t(ref, region, val, this));
  }

  const statement_t *gep_ref(variable_t lhs_ref, variable_t lhs_region,
                             variable_t rhs_ref, variable_t rhs_region,
                             lin_exp_t offset = number_t(0)) {
    return insert(
        new gep_ref_t(lhs_ref, lhs_region, rhs_ref, rhs_region, offset, this));
  }

  const statement_t *assume_ref(ref_cst_t cst) {
    return insert(new assume_ref_t(cst, this));
  }

  const statement_t *assert_ref(ref_cst_t cst, debug_info di = debug_info()) {
    return insert(new assert_ref_t(cst, this, di));
  }

  // (lhs_ref, lhs_rgn) := select_ref(cond, (op1_ref, op1_rgn), (op2_ref,
  // op2_rgn))
  const statement_t *select_ref(variable_t lhs_ref, variable_t lhs_rgn,
                                variable_t cond, variable_t op1_ref,
                                variable_t op1_rgn, variable_t op2_ref,
                                variable_t op2_rgn) {
    return insert(new ref_select_t(lhs_ref, lhs_rgn, cond, op1_ref, op1_rgn,
                                   op2_ref, op2_rgn, this));
  }

  // (lhs_ref, lhs_rgn) := select_ref(cond, null, (op_ref, op_rgn))
  const statement_t *select_ref_null_true_value(variable_t lhs_ref,
                                                variable_t lhs_rgn,
                                                variable_t cond,
                                                variable_t op_ref,
                                                variable_t op_rgn) {
    return insert(new ref_select_t(
        lhs_ref, lhs_rgn, cond, variable_or_constant_t::make_reference_null(),
        boost::none, op_ref, op_rgn, this));
  }

  // (lhs_ref, lhs_rgn) := select_ref(cond, (op_ref, op_rgn), null)
  const statement_t *select_ref_null_false_value(variable_t lhs_ref,
                                                 variable_t lhs_rgn,
                                                 variable_t cond,
                                                 variable_t op_ref,
                                                 variable_t op_rgn) {
    return insert(new ref_select_t(
        lhs_ref, lhs_rgn, cond, op_ref, op_rgn,
        variable_or_constant_t::make_reference_null(), boost::none, this));
  }

  const statement_t *int_to_ref(variable_t int_var, variable_t region,
                                variable_t ref_var) {
    return insert(new int_to_ref_t(int_var, region, ref_var, this));
  }

  const statement_t *ref_to_int(variable_t region, variable_t ref_var,
                                variable_t int_var) {
    return insert(new ref_to_int_t(region, ref_var, int_var, this));
  }

  const statement_t *bool_assign(variable_t lhs, lin_cst_t rhs) {
    return insert(new bool_assign_cst_t(lhs, rhs, this));
  }

  const statement_t *bool_assign(variable_t lhs, ref_cst_t rhs) {
    return insert(new bool_assign_cst_t(lhs, rhs, this));
  }

  const statement_t *bool_assign(variable_t lhs, variable_t rhs,
                                 bool is_not_rhs = false) {
    return insert(new bool_assign_var_t(lhs, rhs, is_not_rhs, this));
  }

  const statement_t *bool_not_assign(variable_t lhs, variable_t rhs) {
    return insert(new bool_assign_var_t(lhs, rhs, true, this));
  }
  
  const statement_t *bool_assume(variable_t c) {
    return insert(new bool_assume_t(c, false, this));
  }

  const statement_t *bool_not_assume(variable_t c) {
    return insert(new bool_assume_t(c, true, this));
  }

  const statement_t *bool_assert(variable_t c, debug_info di = debug_info()) {
    return insert(new bool_assert_t(c, this, di));
  }

  const statement_t *bool_select(variable_t lhs, variable_t cond, variable_t b1,
                                 variable_t b2) {
    return insert(new bool_select_t(lhs, cond, b1, b2, this));
  }

  const statement_t *bool_and(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bool_bin_op_t(lhs, BINOP_BAND, op1, op2, this));
  }

  const statement_t *bool_or(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bool_bin_op_t(lhs, BINOP_BOR, op1, op2, this));
  }

  const statement_t *bool_xor(variable_t lhs, variable_t op1, variable_t op2) {
    return insert(new bool_bin_op_t(lhs, BINOP_BXOR, op1, op2, this));
  }
  /**********************************************************************/
  /*                  End BUILD CFG STATEMENTS                          */
  /**********************************************************************/

  friend crab_os &operator<<(crab_os &o, const basic_block_t &b) {
    b.write(o);
    return o;
  }
};

// Viewing a BasicBlock with all statements reversed. Useful for
// backward analysis.
template <class BasicBlock>
class basic_block_rev {
public:
  using number_t = typename BasicBlock::number_t;
  using varname_t = typename BasicBlock::varname_t;
  using variable_t = typename BasicBlock::variable_t;
  using basic_block_label_t = typename BasicBlock::basic_block_label_t;

  using basic_block_rev_t = basic_block_rev<BasicBlock>;

  using succ_iterator = typename BasicBlock::succ_iterator;
  using const_succ_iterator = typename BasicBlock::const_succ_iterator;
  using pred_iterator = succ_iterator;
  using const_pred_iterator = const_succ_iterator;

  using iterator = typename BasicBlock::reverse_iterator;
  using const_iterator = typename BasicBlock::const_reverse_iterator;
  using live_domain_t = ikos::discrete_domain<variable_t>;

private:
  BasicBlock &_bb;

public:
  basic_block_rev(BasicBlock &bb) : _bb(bb) {}

  const basic_block_label_t &label() const { return _bb.label(); }

  std::string name() const { return _bb.name(); }

  iterator begin() { return _bb.rbegin(); }

  iterator end() { return _bb.rend(); }

  const_iterator begin() const { return _bb.rbegin(); }

  const_iterator end() const { return _bb.rend(); }

  std::size_t size() const { return _bb.size(); }

  void accept(statement_visitor<basic_block_label_t, number_t, varname_t> *v) {
    v->visit(*this);
  }

  live_domain_t &live() { return _bb.live(); }

  const live_domain_t &live() const { return _bb.live(); }

  std::pair<succ_iterator, succ_iterator> next_blocks() {
    return _bb.prev_blocks();
  }

  std::pair<pred_iterator, pred_iterator> prev_blocks() {
    return _bb.next_blocks();
  }

  std::pair<const_succ_iterator, const_succ_iterator> next_blocks() const {
    return _bb.prev_blocks();
  }

  std::pair<const_pred_iterator, const_pred_iterator> prev_blocks() const {
    return _bb.next_blocks();
  }

  void write(crab_os &o) const {
    o << name() << ":\n";
    for (auto const &s : *this) {
      o << "  " << s << ";\n";
    }
    o << "--> [";
    for (auto const &n : boost::make_iterator_range(next_blocks())) {
      o << n << ";";
    }
    o << "]\n";
  }

  // for gdb
  void dump() const { write(crab::errs()); }

  friend crab_os &operator<<(crab_os &o, const basic_block_rev_t &b) {
    b.write(o);
    return o;
  }
};

/** Visitor class for statements **/
template <class BasicBlockLabel, class Number, class VariableName>
struct statement_visitor {
  using bin_op_t = binary_op<BasicBlockLabel, Number, VariableName>;
  using assign_t = assignment<BasicBlockLabel, Number, VariableName>;
  using assume_t = assume_stmt<BasicBlockLabel, Number, VariableName>;
  using select_t = select_stmt<BasicBlockLabel, Number, VariableName>;
  using assert_t = assert_stmt<BasicBlockLabel, Number, VariableName>;
  using int_cast_t = int_cast_stmt<BasicBlockLabel, Number, VariableName>;
  using havoc_t = havoc_stmt<BasicBlockLabel, Number, VariableName>;
  using unreach_t = unreachable_stmt<BasicBlockLabel, Number, VariableName>;
  using callsite_t = callsite_stmt<BasicBlockLabel, Number, VariableName>;
  using intrinsic_t = intrinsic_stmt<BasicBlockLabel, Number, VariableName>;
  using arr_init_t = array_init_stmt<BasicBlockLabel, Number, VariableName>;
  using arr_store_t = array_store_stmt<BasicBlockLabel, Number, VariableName>;
  using arr_load_t = array_load_stmt<BasicBlockLabel, Number, VariableName>;
  using arr_assign_t = array_assign_stmt<BasicBlockLabel, Number, VariableName>;
  using make_ref_t = make_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using remove_ref_t = remove_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using region_init_t = region_init_stmt<BasicBlockLabel, Number, VariableName>;
  using region_copy_t = region_copy_stmt<BasicBlockLabel, Number, VariableName>;
  using region_cast_t = region_cast_stmt<BasicBlockLabel, Number, VariableName>;  
  using load_from_ref_t =
      load_from_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using store_to_ref_t =
      store_to_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using gep_ref_t = gep_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using assume_ref_t = assume_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using assert_ref_t = assert_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using select_ref_t = ref_select_stmt<BasicBlockLabel, Number, VariableName>;
  using int_to_ref_t = int_to_ref_stmt<BasicBlockLabel, Number, VariableName>;
  using ref_to_int_t = ref_to_int_stmt<BasicBlockLabel, Number, VariableName>;
  using bool_bin_op_t = bool_binary_op<BasicBlockLabel, Number, VariableName>;
  using bool_assign_cst_t =
      bool_assign_cst<BasicBlockLabel, Number, VariableName>;
  using bool_assign_var_t =
      bool_assign_var<BasicBlockLabel, Number, VariableName>;
  using bool_assume_t = bool_assume_stmt<BasicBlockLabel, Number, VariableName>;
  using bool_select_t = bool_select_stmt<BasicBlockLabel, Number, VariableName>;
  using bool_assert_t = bool_assert_stmt<BasicBlockLabel, Number, VariableName>;

  virtual void visit(bin_op_t &){};
  virtual void visit(assign_t &){};
  virtual void visit(assume_t &){};
  virtual void visit(select_t &){};
  virtual void visit(assert_t &){};
  virtual void visit(int_cast_t &){};
  virtual void visit(unreach_t &){};
  virtual void visit(havoc_t &){};
  virtual void visit(callsite_t &){};
  virtual void visit(intrinsic_t &){};
  virtual void visit(arr_init_t &){};
  virtual void visit(arr_store_t &){};
  virtual void visit(arr_load_t &){};
  virtual void visit(arr_assign_t &){};
  virtual void visit(region_init_t &) {}
  virtual void visit(region_copy_t &) {}
  virtual void visit(region_cast_t &) {}  
  virtual void visit(make_ref_t &) {}
  virtual void visit(remove_ref_t &) {}  
  virtual void visit(load_from_ref_t &) {}
  virtual void visit(store_to_ref_t &) {}
  virtual void visit(gep_ref_t &) {}
  virtual void visit(assume_ref_t &){};
  virtual void visit(assert_ref_t &){};
  virtual void visit(select_ref_t &){};
  virtual void visit(int_to_ref_t &){};
  virtual void visit(ref_to_int_t &){};
  virtual void visit(bool_bin_op_t &){};
  virtual void visit(bool_assign_cst_t &){};
  virtual void visit(bool_assign_var_t &){};
  virtual void visit(bool_assume_t &){};
  virtual void visit(bool_select_t &){};
  virtual void visit(bool_assert_t &){};

  void visit(basic_block<BasicBlockLabel, VariableName, Number> &b) {
    for (auto &s : b) {
      s.accept(this);
    }
  }

  template <typename BasicBlock> void visit(basic_block_rev<BasicBlock> &b) {
    for (auto &s : b) {
      s.accept(this);
    }
  }

  virtual ~statement_visitor() {}
};

template <class Number, class VariableName>
class function_decl {
public:
  using variable_t = variable<Number, VariableName>;
  using type_t = typename variable_t::type_t;

private:
  std::string m_func_name;
  std::vector<variable_t> m_inputs;
  std::vector<variable_t> m_outputs;

  using param_iterator = typename std::vector<variable_t>::iterator;
  using const_param_iterator = typename std::vector<variable_t>::const_iterator;
  using this_type = function_decl<Number, VariableName>;

public:
  function_decl() : m_func_name("") {}

  function_decl(std::string func_name, std::vector<variable_t> inputs,
                std::vector<variable_t> outputs)
      : m_func_name(func_name), m_inputs(inputs), m_outputs(outputs) {

    // CFG restriction: inputs and outputs must be disjoint,
    // otherwise we cannot produce meaningful input-output
    // relations.
    std::set<variable_t> s;
    for (auto &tv : m_inputs) {
      s.insert(tv);
    }
    for (auto &tv : m_outputs) {
      s.insert(tv);
    }
    if (s.size() != (m_inputs.size() + m_outputs.size())) {
      CRAB_ERROR("interprocedural analysis requires that for each function ",
                 "its set of inputs and outputs must be disjoint.");
    }
  }

  function_decl(const this_type &o)
      : m_func_name(o.m_func_name), m_inputs(o.m_inputs),
        m_outputs(o.m_outputs) {}

  function_decl(const this_type &&o)
      : m_func_name(std::move(o.m_func_name)), m_inputs(std::move(o.m_inputs)),
        m_outputs(std::move(o.m_outputs)) {}

  this_type &operator=(const this_type &o) {
    if (this != &o) {
      m_func_name = o.m_func_name;
      m_inputs = o.m_inputs;
      m_outputs = o.m_outputs;
    }
    return *this;
  }

  this_type &operator=(const this_type &&o) {
    m_func_name = std::move(o.m_func_name);
    m_inputs = std::move(o.m_inputs);
    m_outputs = std::move(o.m_outputs);
    return *this;
  }
  
  bool operator==(const this_type &o) const {
    if (m_func_name != o.m_func_name) {
      return false;
    }

    unsigned ninputs = get_num_inputs();
    unsigned noutputs = get_num_outputs();

    if (ninputs != o.get_num_inputs()) {
      return false;
    }

    if (noutputs != o.get_num_outputs()) {
      return false;
    }

    for (unsigned i = 0, e = ninputs; i < e; ++i) {
      if (get_input_type(i) != o.get_input_type(i)) {
        return false;
      }
    }

    for (unsigned i = 0, e = noutputs; i < e; ++i) {
      if (get_output_type(i) != o.get_output_type(i)) {
        return false;
      }
    }

    return true;
  }

  const std::string &get_func_name() const { return m_func_name; }

  const std::vector<variable_t> &get_inputs() const { return m_inputs; }

  const std::vector<variable_t> &get_outputs() const { return m_outputs; }

  unsigned get_num_inputs() const { return m_inputs.size(); }

  unsigned get_num_outputs() const { return m_outputs.size(); }

  const variable_t &get_input_name(unsigned idx) const {
    if (idx >= m_inputs.size())
      CRAB_ERROR("Out-of-bound access to function input parameter");
    return m_inputs[idx];
  }

  type_t get_input_type(unsigned idx) const {
    if (idx >= m_inputs.size())
      CRAB_ERROR("Out-of-bound access to function output parameter");
    return m_inputs[idx].get_type();
  }

  const variable_t &get_output_name(unsigned idx) const {
    if (idx >= m_outputs.size())
      CRAB_ERROR("Out-of-bound access to function input parameter");
    return m_outputs[idx];
  }

  type_t get_output_type(unsigned idx) const {
    if (idx >= m_outputs.size())
      CRAB_ERROR("Out-of-bound access to function output parameter");
    return m_outputs[idx].get_type();
  }

  void write(crab_os &o) const {

    if (m_outputs.empty()) {
      o << "void";
    } else if (m_outputs.size() == 1) {
      auto out = *(m_outputs.begin());
      o << out << ":" << out.get_type();
    } else {
      o << "(";
      for (auto It = m_outputs.begin(), Et = m_outputs.end(); It != Et;) {
        auto out = *It;
        o << out << ":" << out.get_type();
        ++It;
        if (It != Et)
          o << ",";
      }
      o << ")";
    }

    o << " declare " << m_func_name << "(";
    for (const_param_iterator It = m_inputs.begin(), Et = m_inputs.end();
         It != Et;) {
      o << (*It) << ":" << (*It).get_type();
      ++It;
      if (It != Et)
        o << ",";
    }
    o << ")";
  }

  friend crab_os &operator<<(crab_os &o,
                             const function_decl<Number, VariableName> &decl) {
    decl.write(o);
    return o;
  }
};

// forward declarations
template <class Any> class cfg_rev;
template <class Any> class cfg_ref;

template <class BasicBlockLabel, class VariableName, class Number>
class cfg {
public:
  using number_t = Number;
  using basic_block_label_t = BasicBlockLabel;
  using node_t = basic_block_label_t; // for Bgl graphs
  using varname_t = VariableName;
  using variable_t = variable<number_t, varname_t>;
  using fdecl_t = function_decl<number_t, varname_t>;
  using basic_block_t = basic_block<BasicBlockLabel, VariableName, number_t>;
  using statement_t = statement<BasicBlockLabel, number_t, VariableName>;

  using succ_iterator = typename basic_block_t::succ_iterator;
  using pred_iterator = typename basic_block_t::pred_iterator;
  using const_succ_iterator = typename basic_block_t::const_succ_iterator;
  using const_pred_iterator = typename basic_block_t::const_pred_iterator;

  using succ_range = boost::iterator_range<succ_iterator>;
  using pred_range = boost::iterator_range<pred_iterator>;
  using const_succ_range = boost::iterator_range<const_succ_iterator>;
  using const_pred_range = boost::iterator_range<const_pred_iterator>;

private:
  using cfg_t = cfg<BasicBlockLabel, VariableName, Number>;
  using basic_block_map_t =
      std::unordered_map<BasicBlockLabel, basic_block_t *>;
  using binding_t = typename basic_block_map_t::value_type;
  using live_domain_t = typename basic_block_t::live_domain_t;

  struct get_ref : public std::unary_function<binding_t, basic_block_t> {
    get_ref() {}
    basic_block_t &operator()(const binding_t &p) const { return *(p.second); }
  };

  struct get_label : public std::unary_function<binding_t, BasicBlockLabel> {
    get_label() {}
    BasicBlockLabel operator()(const binding_t &p) const {
      return p.second->label();
    }
  };

public:
  using iterator =
      boost::transform_iterator<get_ref, typename basic_block_map_t::iterator>;
  using const_iterator =
      boost::transform_iterator<get_ref,
                                typename basic_block_map_t::const_iterator>;
  using label_iterator =
      boost::transform_iterator<get_label,
                                typename basic_block_map_t::iterator>;
  using const_label_iterator =
      boost::transform_iterator<get_label,
                                typename basic_block_map_t::const_iterator>;

  using var_iterator = typename std::vector<varname_t>::iterator;
  using const_var_iterator = typename std::vector<varname_t>::const_iterator;

private:
  BasicBlockLabel m_entry;
  boost::optional<BasicBlockLabel> m_exit;
  basic_block_map_t m_blocks;
  fdecl_t m_func_decl;

  using visited_t = std::unordered_set<BasicBlockLabel>;
  template <typename T>
  void dfs_rec(BasicBlockLabel curId, visited_t &visited, T f) const {
    if (!visited.insert(curId).second)
      return;

    const basic_block_t &cur = get_node(curId);
    f(cur);
    for (auto const &n : boost::make_iterator_range(cur.next_blocks())) {
      dfs_rec(n, visited, f);
    }
  }

  template <typename T> void dfs(T f) const {
    visited_t visited;
    dfs_rec(m_entry, visited, f);
  }

  struct print_block {
    crab_os &m_o;
    print_block(crab_os &o) : m_o(o) {}
    void operator()(const basic_block_t &B) { B.write(m_o); }
  };

public:
  // --- needed by crab::cg::call_graph<CFG>::cg_node
  cfg() {}

  cfg(const cfg_t &o) = delete;

  cfg_t &operator=(const cfg_t &o) = delete;

  cfg(BasicBlockLabel entry) : m_entry(entry), m_exit(boost::none) {
    m_blocks.insert(binding_t(m_entry, basic_block_t::create(m_entry)));
  }

  cfg(BasicBlockLabel entry, BasicBlockLabel exit)
      : m_entry(entry), m_exit(exit) {
    m_blocks.insert(binding_t(m_entry, basic_block_t::create(m_entry)));
  }

  cfg(BasicBlockLabel entry, BasicBlockLabel exit, fdecl_t func_decl)
      : m_entry(entry), m_exit(exit), m_func_decl(func_decl) {
    m_blocks.insert(binding_t(m_entry, basic_block_t::create(m_entry)));
  }

  // The cfg owns the basic blocks
  ~cfg() {
    for (auto &kv : m_blocks) {
      delete kv.second;
    }
  }

  cfg_t *clone() const {
    cfg_t *copy_cfg = new cfg_t();
    copy_cfg->m_entry = m_entry;
    copy_cfg->m_exit = m_exit;
    copy_cfg->m_func_decl = m_func_decl;

    for (auto const &bb : boost::make_iterator_range(begin(), end())) {
      basic_block_t *copy_bb = bb.clone();
      copy_cfg->m_blocks.insert(binding_t(copy_bb->label(), copy_bb));
    }
    return copy_cfg;
  }

  bool has_func_decl() const {
    return (!(m_func_decl.get_func_name() == "" &&
              m_func_decl.get_num_inputs() == 0 &&
              m_func_decl.get_num_outputs() == 0));
  }

  const fdecl_t &get_func_decl() const { return m_func_decl; }

  bool has_exit() const { return (bool)m_exit; }

  const BasicBlockLabel &exit() const {
    if (has_exit())
      return *m_exit;
    CRAB_ERROR("cfg does not have an exit block");
  }

  //! set method to mark the exit block after the cfg has been
  //! created.
  void set_exit(BasicBlockLabel exit) { m_exit = exit; }

  //! set method to add the function declaration after the cfg has
  //! been created.
  void set_func_decl(fdecl_t decl) { m_func_decl = decl; }

  // --- Begin ikos fixpoint API

  const BasicBlockLabel &entry() const { return m_entry; }

  const_succ_range next_nodes(BasicBlockLabel bb_id) const {
    const basic_block_t &b = get_node(bb_id);
    return boost::make_iterator_range(b.next_blocks());
  }

  const_pred_range prev_nodes(BasicBlockLabel bb_id) const {
    const basic_block_t &b = get_node(bb_id);
    return boost::make_iterator_range(b.prev_blocks());
  }

  succ_range next_nodes(BasicBlockLabel bb_id) {
    basic_block_t &b = get_node(bb_id);
    return boost::make_iterator_range(b.next_blocks());
  }

  pred_range prev_nodes(BasicBlockLabel bb_id) {
    basic_block_t &b = get_node(bb_id);
    return boost::make_iterator_range(b.prev_blocks());
  }

  basic_block_t &get_node(BasicBlockLabel bb_id) {
    auto it = m_blocks.find(bb_id);
    if (it == m_blocks.end()) {
      CRAB_ERROR("Basic block ", bb_id, " not found in the CFG: ", __LINE__);
    }

    return *(it->second);
  }

  const basic_block_t &get_node(BasicBlockLabel bb_id) const {
    auto it = m_blocks.find(bb_id);
    if (it == m_blocks.end()) {
      CRAB_ERROR("Basic block ", bb_id, " not found in the CFG: ", __LINE__);
    }

    return *(it->second);
  }

  // --- End ikos fixpoint API
  
  basic_block_t &insert(BasicBlockLabel bb_id) {
    auto it = m_blocks.find(bb_id);
    if (it != m_blocks.end())
      return *(it->second);

    basic_block_t *block = basic_block_t::create(bb_id);
    m_blocks.insert(binding_t(bb_id, block));
    return *block;
  }

  void remove(BasicBlockLabel bb_id) {
    if (bb_id == m_entry) {
      CRAB_ERROR("Cannot remove entry block");
    }

    if (m_exit && *m_exit == bb_id) {
      CRAB_ERROR("Cannot remove exit block");
    }

    std::vector<std::pair<basic_block_t *, basic_block_t *>> dead_edges;
    basic_block_t *bb = &(get_node(bb_id));

    for (auto id : boost::make_iterator_range(bb->prev_blocks())) {
      if (bb_id != id) {
        basic_block_t &p = get_node(id);
        dead_edges.push_back({&p, bb});
      }
    }

    for (auto id : boost::make_iterator_range(bb->next_blocks())) {
      if (bb_id != id) {
        basic_block_t &s = get_node(id);
        dead_edges.push_back({bb, &s});
      }
    }

    for (auto p : dead_edges) {
      (*p.first) -= (*p.second);
    }

    m_blocks.erase(bb_id);
    delete bb;
  }

  // Return a fresh basic block that is not part of the current CFG.
  std::unique_ptr<basic_block_t> create_unlinked_block(BasicBlockLabel bb) {
    return std::unique_ptr<basic_block_t>(new basic_block_t(bb));
  }
  
  // Return all variables (either used or defined) in the cfg.
  //
  // This operation is linear on the size of the cfg to still keep
  // a valid set in case a block is removed.
  std::vector<varname_t> get_vars() const {
    live_domain_t ls = live_domain_t::bottom();
    for (auto const &b : boost::make_iterator_range(begin(), end()))
      ls = ls | b.live();
    // std::vector<varname_t> vars(ls.size());
    // vars.insert(vars.end(), ls.begin(), ls.end());
    std::vector<varname_t> vars;
    for (auto v : ls)
      vars.push_back(v);
    return vars;
  }

  //! return a begin iterator of BasicBlock's
  iterator begin() {
    return boost::make_transform_iterator(m_blocks.begin(), get_ref());
  }

  //! return an end iterator of BasicBlock's
  iterator end() {
    return boost::make_transform_iterator(m_blocks.end(), get_ref());
  }

  const_iterator begin() const {
    return boost::make_transform_iterator(m_blocks.begin(), get_ref());
  }

  const_iterator end() const {
    return boost::make_transform_iterator(m_blocks.end(), get_ref());
  }

  //! return a begin iterator of BasicBlockLabel's
  label_iterator label_begin() {
    return boost::make_transform_iterator(m_blocks.begin(), get_label());
  }

  //! return an end iterator of BasicBlockLabel's
  label_iterator label_end() {
    return boost::make_transform_iterator(m_blocks.end(), get_label());
  }

  const_label_iterator label_begin() const {
    return boost::make_transform_iterator(m_blocks.begin(), get_label());
  }

  const_label_iterator label_end() const {
    return boost::make_transform_iterator(m_blocks.end(), get_label());
  }

  size_t size() const { return m_blocks.size(); }

  void write(crab_os &o) const {
    if (has_func_decl()) {
      o << m_func_decl << "\n";
    }

    print_block f(o);
    dfs(f);
  }

  // for gdb
  void dump() const {
    crab::errs() << "Number of basic blocks=" << size() << "\n";
    for (auto &bb : boost::make_iterator_range(begin(), end())) {
      bb.dump();
    }
  }

  friend crab_os &operator<<(crab_os &o, const cfg_t &cfg) {
    cfg.write(o);
    return o;
  }

  void simplify() {
    merge_blocks();
    remove_unreachable_blocks();
    remove_useless_blocks();
    // after removing useless blocks there can be opportunities to
    // merge more blocks.
    merge_blocks();
    // merge_blocks();
  }

private:
  ////
  // Trivial cfg simplifications
  // TODO: move to transform directory
  ////

  // Helpers

  basic_block_t &get_child(BasicBlockLabel b) {
    // assert(has_one_child(b));
    auto rng = next_nodes(b);
    return get_node(*(rng.begin()));
  }

  basic_block_t &get_parent(BasicBlockLabel b) {
    // assert(has_one_parent(b));
    auto rng = prev_nodes(b);
    return get_node(*(rng.begin()));
  }

#if 1
  bool has_one_child(BasicBlockLabel b) const {
    auto rng = next_nodes(b);
    return (std::distance(rng.begin(), rng.end()) == 1);
  }

  bool has_one_parent(BasicBlockLabel b) const {
    auto rng = prev_nodes(b);
    return (std::distance(rng.begin(), rng.end()) == 1);
  }

  void merge_blocks_rec(BasicBlockLabel curId, visited_t &visited) {
    if (!visited.insert(curId).second)
      return;

    basic_block_t &cur = get_node(curId);

    if (has_one_child(curId) && has_one_parent(curId)) {
      basic_block_t &parent = get_parent(curId);
      basic_block_t &child = get_child(curId);

      // Merge with its parent if it's its only child.
      if (has_one_child(parent.label())) {
        // fold cur into parent
        parent.copy_back(cur);
        visited.erase(curId);
	if (has_exit() && exit() == curId) {
	  set_exit(parent.label());
	}
        remove(curId);
        parent >> child;
        merge_blocks_rec(child.label(), visited);
        return;
      }
    }

    for (auto n : boost::make_iterator_range(cur.next_blocks())) {
      merge_blocks_rec(n, visited);
    }
  }

  // Merges a basic block into its predecessor if there is only one
  // and the predecessor only has one successor.
  void merge_blocks() {
    visited_t visited;
    merge_blocks_rec(entry(), visited);
  }
#else
  // Non-recursive version thanks to Prevail
  void merge_blocks() {
    std::set<basic_block_label_t> worklist(this->label_begin(),
                                           this->label_end());
    while (!worklist.empty()) {
      auto label = *worklist.begin();
      worklist.erase(label);

      basic_block_t &bb = get_node(label);

      // skip bb for now but it will be folded into its parent when
      // the parent is processed.
      if (bb.in_degree() == 1 && get_parent(label).out_degree() == 1) {
        continue;
      }
      while (bb.out_degree() == 1) {
        basic_block_t &next_bb = get_child(label);

        // skip self-loops or if next_bb has in-degree > 1
        if (&next_bb == &bb || next_bb.in_degree() != 1) {
          break;
        }

        worklist.erase(next_bb.label());

        if (next_bb.label() == m_exit) {
          m_exit = label;
        }

        // fold next_bb into bb
        bb.copy_back(next_bb);
        bb -= next_bb;
        auto children = next_bb.m_next;
        for (auto const &next_next_label : children) {
          basic_block_t &next_next_bb = get_node(next_next_label);
          bb >> next_next_bb;
        }
	if (has_exit() && exit() == next_bb.label()) {
	  set_exit(bb.label());
	}
	remove(next_bb.label());
      }
    }
  }
#endif

  // mark reachable blocks from entry
  template <class AnyCfg>
  void mark_alive_blocks(AnyCfg &cfg, visited_t &alive_set) {
    std::vector<typename AnyCfg::basic_block_label_t> worklist;
    worklist.push_back(cfg.entry());

    while (!worklist.empty()) {
      auto bb_label = worklist.back();
      worklist.pop_back();
      alive_set.insert(bb_label);

      for (auto child : cfg.next_nodes(bb_label)) {
        if (!(alive_set.count(child) > 0)) {
          worklist.push_back(child);
        }
      }
    }
  }

  // remove unreachable blocks
  void remove_unreachable_blocks() {
    visited_t alive, dead;
    mark_alive_blocks(*this, alive);
    
    for (auto const &bb : *this) {
      if (!(alive.count(bb.label()) > 0)) {
        dead.insert(bb.label());
      }
    }

    for (auto bb_id : dead) {
      if (has_exit() && exit() == bb_id) {
	// If the exit block is unreachable then we still keep it because any
	// well-formed CFG should have an exit.
	continue;
      }
      remove(bb_id);
    }
  }

  // remove blocks that cannot reach the exit block
  void remove_useless_blocks() {
    if (!has_exit())
      return;

    cfg_rev<cfg_ref<cfg_t>> rev_cfg(*this);
    visited_t useful;
    std::vector<basic_block_label_t> useless;
    mark_alive_blocks(rev_cfg, useful);
    

    for (auto const &bb : *this) {
      if (!(useful.count(bb.label()) > 0)) {
        useless.push_back(bb.label());
      }
    }

    for (auto bb_id : useless) {
      if (bb_id == entry()) {
	// if the entry block cannot reach the exit then we still keep it 
	// because any well-formed CFG must have an entry block.
	continue;
      }
      remove(bb_id);
    }
  }
};

// A lightweight object that wraps a pointer a CFG into a copyable,
// assignable object.
template <class CFG>
class cfg_ref {
public:
  // CFG's typedefs
  using basic_block_label_t = typename CFG::basic_block_label_t;
  using node_t = typename CFG::node_t;
  using varname_t = typename CFG::varname_t;
  using number_t = typename CFG::number_t;
  using variable_t = typename CFG::variable_t;
  using fdecl_t = typename CFG::fdecl_t;
  using basic_block_t = typename CFG::basic_block_t;
  using statement_t = typename CFG::statement_t;

  using succ_iterator = typename CFG::succ_iterator;
  using pred_iterator = typename CFG::pred_iterator;
  using const_succ_iterator = typename CFG::const_succ_iterator;
  using const_pred_iterator = typename CFG::const_pred_iterator;
  using succ_range = typename CFG::succ_range;
  using pred_range = typename CFG::pred_range;
  using const_succ_range = typename CFG::const_succ_range;
  using const_pred_range = typename CFG::const_pred_range;
  using iterator = typename CFG::iterator;
  using const_iterator = typename CFG::const_iterator;
  using label_iterator = typename CFG::label_iterator;
  using const_label_iterator = typename CFG::const_label_iterator;
  using var_iterator = typename CFG::var_iterator;
  using const_var_iterator = typename CFG::const_var_iterator;

private:
  CFG* m_ref;

public:
  // default constructor needed by crab::cg::CallGraph<CFG>::CgNode
  cfg_ref(): m_ref(nullptr) {}
  
  cfg_ref(CFG &cfg) : m_ref(&cfg) {}
  
  cfg_ref(const cfg_ref<CFG> &o) = default;
  cfg_ref(cfg_ref<CFG> &&o) = default;
  cfg_ref<CFG> &operator=(const cfg_ref<CFG> &o) = default;
  cfg_ref<CFG> &operator=(cfg_ref<CFG> &&o) = default;
  
  const CFG &get() const {
    assert(m_ref);
    return *m_ref;
  }

  CFG &get() {
    assert(m_ref);
    return *m_ref;
  }

  const basic_block_label_t &entry() const {
    return get().entry();
  }

  const_succ_range next_nodes(basic_block_label_t bb) const {
    return get().next_nodes(bb);
  }

  const_pred_range prev_nodes(basic_block_label_t bb) const {
    return get().prev_nodes(bb);
  }

  succ_range next_nodes(basic_block_label_t bb) {
    return get().next_nodes(bb);
  }

  pred_range prev_nodes(basic_block_label_t bb) {
    return get().prev_nodes(bb);
  }

  basic_block_t &get_node(basic_block_label_t bb) {
    return get().get_node(bb);
  }

  const basic_block_t &get_node(basic_block_label_t bb) const {
    return get().get_node(bb);
  }

  size_t size() const {
    return get().size();
  }

  iterator begin() {
    return get().begin();
  }

  iterator end() {
    return get().end();
  }

  const_iterator begin() const {
    return get().begin();
  }

  const_iterator end() const {
    return get().end();
  }

  label_iterator label_begin() {
    return get().label_begin();
  }

  label_iterator label_end() {
    return get().label_end();
  }

  const_label_iterator label_begin() const {
    return get().label_begin();
  }

  const_label_iterator label_end() const {
    return get().label_end();
  }

  bool has_func_decl() const {
    return get().has_func_decl();
  }

  const fdecl_t &get_func_decl() const {
    return get().get_func_decl();
  }

  bool has_exit() const {
    return get().has_exit();
  }

  const basic_block_label_t &exit() const {
    return get().exit();
  }

  friend crab_os &operator<<(crab_os &o, const cfg_ref<CFG> &cfg) {
    o << cfg.get();
    return o;
  }

  // for gdb
  void dump() const {
    get().dump();
  }

  void simplify() {
    get().simplify();
  }
};

// Viewing a CFG with all edges and block statements
// reversed. Useful for backward analysis.
template <class CFGRef> // CFGRef must be copyable!
class cfg_rev {
public:
  using basic_block_label_t = typename CFGRef::basic_block_label_t;
  using basic_block_t = basic_block_rev<typename CFGRef::basic_block_t>;
  using node_t = basic_block_label_t; // for Bgl graphs
  using varname_t = typename CFGRef::varname_t;
  using number_t = typename CFGRef::number_t;
  using variable_t = typename CFGRef::variable_t;
  using fdecl_t = typename CFGRef::fdecl_t;
  using statement_t = typename CFGRef::statement_t;

  using pred_range = typename CFGRef::succ_range;
  using succ_range = typename CFGRef::pred_range;
  using const_pred_range = typename CFGRef::const_succ_range;
  using const_succ_range = typename CFGRef::const_pred_range;

  // For BGL
  using succ_iterator = typename basic_block_t::succ_iterator;
  using pred_iterator = typename basic_block_t::pred_iterator;
  using const_succ_iterator = typename basic_block_t::const_succ_iterator;
  using const_pred_iterator = typename basic_block_t::const_pred_iterator;

  using cfg_rev_t = cfg_rev<CFGRef>;

private:
  struct getRev : public std::unary_function<typename CFGRef::basic_block_t,
                                             basic_block_t> {
    const std::unordered_map<basic_block_label_t, basic_block_t> &_rev_bbs;

    getRev(
        const std::unordered_map<basic_block_label_t, basic_block_t> &rev_bbs)
        : _rev_bbs(rev_bbs) {}

    const basic_block_t &operator()(typename CFGRef::basic_block_t &bb) const {
      auto it = _rev_bbs.find(bb.label());
      if (it != _rev_bbs.end())
        return it->second;
      CRAB_ERROR("Basic block ", bb.label(),
                 " not found in the CFG: ", __LINE__);
    }
  };

  using visited_t = std::unordered_set<basic_block_label_t>;

  template <typename T>
  void dfs_rec(basic_block_label_t curId, visited_t &visited, T f) const {
    if (!visited.insert(curId).second)
      return;
    const basic_block_t &cur = get_node(curId);
    f(cur);
    for (auto const &n : next_nodes(curId)) {
      dfs_rec(n, visited, f);
    }
  }

  template <typename T> void dfs(T f) const {
    visited_t visited;
    dfs_rec(entry(), visited, f);
  }

  struct print_block {
    crab_os &m_o;
    print_block(crab_os &o) : m_o(o) {}
    void operator()(const basic_block_t &B) { B.write(m_o); }
  };

public:
  using iterator = boost::transform_iterator<getRev, typename CFGRef::iterator>;
  using const_iterator =
      boost::transform_iterator<getRev, typename CFGRef::const_iterator>;
  using label_iterator = typename CFGRef::label_iterator;
  using const_label_iterator = typename CFGRef::const_label_iterator;
  using var_iterator = typename CFGRef::var_iterator;
  using const_var_iterator = typename CFGRef::const_var_iterator;

private:
  CFGRef _cfg;
  std::unordered_map<basic_block_label_t, basic_block_t> _rev_bbs;

public:
  // --- hook needed by crab::cg::CallGraph<CFGRef>::CgNode
  cfg_rev() {}

  cfg_rev(CFGRef cfg) : _cfg(cfg) {
    // Create basic_block_rev from BasicBlock objects
    // Note that basic_block_rev is also a view of BasicBlock so it
    // doesn't modify BasicBlock objects.
    for (auto &bb : cfg) {
      basic_block_t rev(bb);
      _rev_bbs.insert(std::make_pair(bb.label(), rev));
    }
  }

  cfg_rev(const cfg_rev_t &o) : _cfg(o._cfg), _rev_bbs(o._rev_bbs) {}

  cfg_rev(cfg_rev_t &&o)
      : _cfg(std::move(o._cfg)), _rev_bbs(std::move(o._rev_bbs)) {}

  cfg_rev_t &operator=(const cfg_rev_t &o) {
    if (this != &o) {
      _cfg = o._cfg;
      _rev_bbs = o._rev_bbs;
    }
    return *this;
  }

  cfg_rev_t &operator=(cfg_rev_t &&o) {
    _cfg = std::move(o._cfg);
    _rev_bbs = std::move(o._rev_bbs);
    return *this;
  }

  const basic_block_label_t &entry() const {
    if (!_cfg.has_exit())
      CRAB_ERROR("Entry not found!");
    return _cfg.exit();
  }

  const_succ_range next_nodes(basic_block_label_t bb) const {
    return _cfg.prev_nodes(bb);
  }

  const_pred_range prev_nodes(basic_block_label_t bb) const {
    return _cfg.next_nodes(bb);
  }

  succ_range next_nodes(basic_block_label_t bb) { return _cfg.prev_nodes(bb); }

  pred_range prev_nodes(basic_block_label_t bb) { return _cfg.next_nodes(bb); }

  basic_block_t &get_node(basic_block_label_t bb_id) {
    auto it = _rev_bbs.find(bb_id);
    if (it == _rev_bbs.end())
      CRAB_ERROR("Basic block ", bb_id, " not found in the CFG: ", __LINE__);
    return it->second;
  }

  const basic_block_t &get_node(basic_block_label_t bb_id) const {
    auto it = _rev_bbs.find(bb_id);
    if (it == _rev_bbs.end())
      CRAB_ERROR("Basic block ", bb_id, " not found in the CFG: ", __LINE__);
    return it->second;
  }

  iterator begin() {
    return boost::make_transform_iterator(_cfg.begin(), getRev(_rev_bbs));
  }

  iterator end() {
    return boost::make_transform_iterator(_cfg.end(), getRev(_rev_bbs));
  }

  const_iterator begin() const {
    return boost::make_transform_iterator(_cfg.begin(), getRev(_rev_bbs));
  }

  const_iterator end() const {
    return boost::make_transform_iterator(_cfg.end(), getRev(_rev_bbs));
  }

  label_iterator label_begin() { return _cfg.label_begin(); }

  label_iterator label_end() { return _cfg.label_end(); }

  const_label_iterator label_begin() const { return _cfg.label_begin(); }

  const_label_iterator label_end() const { return _cfg.label_end(); }

  bool has_func_decl() const { return _cfg.has_func_decl(); }

  const fdecl_t &get_func_decl() const { return _cfg.get_func_decl(); }

  bool has_exit() const { return true; }

  const basic_block_label_t &exit() const { return _cfg.entry(); }

  void write(crab_os &o) const {
    if (has_func_decl()) {
      o << get_func_decl() << "\n";
    }
    print_block f(o);
    dfs(f);
  }

  friend crab_os &operator<<(crab_os &o, const cfg_rev<CFGRef> &cfg) {
    cfg.write(o);
    return o;
  }

  void simplify() {}
};

// Wrapper to allow compare callsites and function declarations. This
// class can be avoided, if a callsite keeps a pointer to the callee's
// CFG function and function declarations keep a pointer to its
// CFG. That would require to build CFGs in two phases.  For now, we
// have this wrapper that it's less efficient but hopefully it's not a
// big bottleneck.
template<typename CFG>  
class callsite_or_fdecl {
public:
  using callsite_t = typename CFG::basic_block_t::callsite_t;
  using fdecl_t = typename CFG::fdecl_t;
private:
  const callsite_t *m_cs;
  const fdecl_t *m_fdecl;

  static bool match_callee_and_function(const callsite_t &cs, const fdecl_t &fdecl) {
    if (cs.get_func_name() != fdecl.get_func_name() ||
	cs.get_num_args() != fdecl.get_num_inputs() ||
	cs.get_num_lhs() != fdecl.get_num_outputs()) {
      return false;
    }
    for (unsigned i=0, args=fdecl.get_num_inputs(); i<args;++i) {
      if (cs.get_arg_type(i) != fdecl.get_input_type(i)) {
	return false;
      }
    }
    for (unsigned i=0, args=fdecl.get_num_outputs(); i<args;++i) {
      if (cs.get_lhs_type(i) != fdecl.get_output_type(i)) {
	return false;
      }
    }
    return true;
  }

  static bool match_callees(const callsite_t &cs1, const callsite_t &cs2) {
    if (cs1.get_func_name() != cs2.get_func_name() ||
	cs1.get_num_args() != cs2.get_num_args() ||
	cs1.get_num_lhs() != cs2.get_num_lhs()) {
      return false;
    }
    for (unsigned i=0, args=cs1.get_num_args(); i<args;++i) {
      if (cs1.get_arg_type(i) != cs2.get_arg_type(i)) {
	return false;
      }
    }
    for (unsigned i=0, args=cs1.get_num_lhs(); i<args;++i) {
      if (cs1.get_lhs_type(i) != cs2.get_lhs_type(i)) {
	return false;
      }
    }
    return true;
  }

  static bool match_functions(const fdecl_t &fd1, const fdecl_t &fd2) {
    if (fd1.get_func_name() != fd2.get_func_name() ||
	fd1.get_num_inputs() != fd2.get_num_inputs() ||
	fd1.get_num_outputs() != fd2.get_num_outputs()) {
      return false;
    }
    for (unsigned i=0, args=fd1.get_num_inputs(); i<args;++i) {
      if (fd1.get_input_type(i) != fd2.get_input_type(i)) {
	return false;
      }
    }
    for (unsigned i=0, args=fd1.get_num_outputs(); i<args;++i) {
      if (fd1.get_output_type(i) != fd2.get_output_type(i)) {
	return false;
      }
    }
    return true;
  }
  
  static size_t combine(size_t seed, size_t hash_val) {
    // Similar to boost::hash_combine
    seed ^= hash_val + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }

  static size_t compute_hash(const callsite_t &cs) {
    size_t res = boost::hash_value(cs.get_func_name());
    for (unsigned i=0, args=cs.get_num_args(); i < args; ++i) {
      combine(res, cs.get_arg_type(i).hash());
    }
    for (unsigned i=0, args=cs.get_num_lhs(); i < args; ++i) {
      combine(res, cs.get_lhs_type(i).hash());
    }
    return res;
  }

  static size_t compute_hash(const fdecl_t &d) {
    size_t res = boost::hash_value(d.get_func_name());
    for (unsigned i=0, args=d.get_num_inputs(); i < args; ++i) {
      combine(res, d.get_input_type(i).hash());
    }
    for (unsigned i=0, args=d.get_num_outputs(); i < args; ++i) {
      combine(res, d.get_output_type(i).hash());
    }
    return res;
  }
  
  public:
  callsite_or_fdecl(const callsite_t *cs): m_cs(cs), m_fdecl(nullptr) {}
  callsite_or_fdecl(const fdecl_t *fdecl): m_cs(nullptr), m_fdecl(fdecl) {}
  callsite_or_fdecl(const callsite_or_fdecl &o) = default;
  callsite_or_fdecl(callsite_or_fdecl &&o) = default;
  callsite_or_fdecl& operator=(const callsite_or_fdecl &o) = default;
  callsite_or_fdecl& operator=(callsite_or_fdecl &&o) = default;

  // non-standard equality: two callsites or function declarations are
  // equal if they refer to the same function. For callsites, we mean
  // the called function.
  bool operator==(const callsite_or_fdecl &o) const {
    if (m_cs && o.m_cs) {
      // return true iff the two callsites call the same function
      return match_callees(*m_cs, *o.m_cs);
    } else if (m_fdecl && o.m_fdecl) {
      // return true iff the two function declarations refer to the
      // same function.
      return match_functions(*m_fdecl, *o.m_fdecl);
    } else if (m_cs && o.m_fdecl) {
      // return true iff the called function is the same than the one
      // referred by the function declaration
      return match_callee_and_function(*m_cs, *(o.m_fdecl));
    } else {
      // as above
      assert(m_fdecl && o.m_cs);
      return match_callee_and_function(*(o.m_cs), *m_fdecl);
    } 
  }

  size_t hash() const {
    return (m_cs ? compute_hash(*m_cs): compute_hash(*m_fdecl));
  }

  void write(crab_os &o) const {
    if (m_cs) {
      m_cs->write(o);
    } else {
      m_fdecl->write(o);
    }
  }
};

namespace callsite_or_fdecl_details {
template<class CFG>
struct hasher {
  size_t operator()(const callsite_or_fdecl<CFG> &o) const {
    return o.hash();
  }
};
template<class CFG>  
struct compare {
  bool operator()(const callsite_or_fdecl<CFG> &o1,
		  const callsite_or_fdecl<CFG> &o2) const {
    return o1 == o2;
  }
};
} // end namespace callsite_or_fdecl_details

// Specialized unordered map when the key if callsite_or_fdecl.
template<class CFG, class Value>  
using callsite_or_fdecl_map =
  std::unordered_map<callsite_or_fdecl<CFG>, Value,
		     callsite_or_fdecl_details::hasher<CFG>,
		     callsite_or_fdecl_details::compare<CFG>>;
    
} // end namespace cfg
} // end namespace crab

namespace crab {
namespace cfg {
template <class B, class V, class N>
bool operator==(cfg<B, V, N> const &a, cfg<B, V, N> const &b) {
  if (!a.has_func_decl() || !b.has_func_decl()) {
    CRAB_ERROR("cannot call operator== of a cfg because function declaration "
               "is missing");
  }
  return (a.get_func_decl() == b.get_func_decl());
}

template <class CFG>
bool operator==(cfg_ref<CFG> const &a, cfg_ref<CFG> const &b) {
  if (!a.has_func_decl() || !b.has_func_decl()) {
    CRAB_ERROR("cannot call operator== of a cfg because function declaration "
               "is missing");
  }
  return (a.get_func_decl() == b.get_func_decl());
}

template <class CFGRef>
bool operator==(cfg_rev<CFGRef> const &a, cfg_rev<CFGRef> const &b) {
  if (!a.has_func_decl() || !b.has_func_decl()) {
    CRAB_ERROR("cannot call operator== of a cfg because function declaration "
               "is missing");
  }
  return (a.get_func_decl() == b.get_func_decl());
}
} // namespace cfg

template <typename BasicBlock>
class basic_block_traits<cfg::basic_block_rev<BasicBlock>> {
public:
  using basic_block_label_t = typename BasicBlock::basic_block_label_t;
  static std::string to_string(const basic_block_label_t &bb) {
    return basic_block_traits<BasicBlock>::to_string(bb);
  }
};

} // namespace crab

namespace std {
/**  specialization of std::hash for cfg variants **/
template <class B, class V, class N> struct hash<crab::cfg::cfg<B, V, N>> {
  using cfg_t = crab::cfg::cfg<B, V, N>;
  size_t operator()(const cfg_t &_cfg) const {
    if (!_cfg.has_func_decl()) {
      CRAB_ERROR("cannot hash a cfg because function declaration is missing");
    }
    crab::cfg::callsite_or_fdecl<cfg_t> sig(_cfg.get_func_decl());
    return sig.hash();
  }
};

template <class CFG> struct hash<crab::cfg::cfg_ref<CFG>> {
  using cfg_ref_t = crab::cfg::cfg_ref<CFG>;
  size_t operator()(const cfg_ref_t &_cfg) const {
    if (!_cfg.has_func_decl()) {
      CRAB_ERROR("cannot hash a cfg because function declaration is missing");
    }
    crab::cfg::callsite_or_fdecl<cfg_ref_t> sig(&(_cfg.get_func_decl()));
    return sig.hash();
  }
};

template <class CFGRef> struct hash<crab::cfg::cfg_rev<CFGRef>> {
  using cfg_rev_t = crab::cfg::cfg_rev<CFGRef>;
  size_t operator()(const cfg_rev_t &_cfg) const {
    if (!_cfg.has_func_decl()) {
      CRAB_ERROR("cannot hash a cfg because function declaration is missing");
    }
    crab::cfg::callsite_or_fdecl<cfg_rev_t> sig(&(_cfg.get_func_decl()));
    return sig.hash();
  }
};
} // end namespace std
