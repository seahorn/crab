#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/cfg/cfg_to_json.hpp>
#include <crab/support/json.hpp>
#include <crab/types/linear_constraints_to_json.hpp>
#include <crab/types/reference_constraints_to_json.hpp>
#include <crab/types/variable_to_json.hpp>

#include <boost/core/lightweight_test.hpp>

#include <algorithm>
#include <string>
#include <vector>

/**
 * Tests for the JSON exporters: support/json.hpp, types/variable_to_json.hpp,
 * types/linear_constraints_to_json.hpp and cfg/cfg_to_json.hpp.
 *
 * Uses Boost.Core lightweight test: checks are silent on success and report
 * failures to stderr, so this test contributes no stdout and needs no entry in
 * expected_results.out.
 *
 * Expectations are written as single-line JSON and compared against the
 * exporter's pretty-printed output run through `compact()`, which drops the
 * whitespace *between* tokens but keeps whitespace inside strings. That keeps
 * the expectations readable without pinning the indentation.
 */

using namespace crab::cfg_impl;

using z_lin_cst_sys_t =
    ikos::linear_constraint_system<ikos::z_number, varname_t>;
using z_dis_lin_cst_sys_t =
    ikos::disjunctive_linear_constraint_system<ikos::z_number, varname_t>;

namespace {

// Drop whitespace that sits between JSON tokens, keeping whatever is inside
// string literals (including escapes, so that a \" does not end the string).
std::string compact(const std::string &s) {
  std::string res;
  bool in_string = false;
  for (size_t i = 0; i < s.size(); ++i) {
    char c = s[i];
    if (in_string) {
      res += c;
      if (c == '\\' && i + 1 < s.size()) {
        res += s[++i];
      } else if (c == '"') {
        in_string = false;
      }
      continue;
    }
    if (c == '"') {
      in_string = true;
      res += c;
    } else if (c != ' ' && c != '\n' && c != '\t' && c != '\r') {
      res += c;
    }
  }
  return res;
}

// Brackets are balanced and every string literal is closed. Catches a writer
// that forgets an end_object/end_array or emits an unescaped quote.
bool well_formed(const std::string &s) {
  std::vector<char> stack;
  bool in_string = false;
  for (size_t i = 0; i < s.size(); ++i) {
    char c = s[i];
    if (in_string) {
      if (c == '\\') {
        ++i;
      } else if (c == '"') {
        in_string = false;
      }
      continue;
    }
    switch (c) {
    case '"':
      in_string = true;
      break;
    case '{':
    case '[':
      stack.push_back(c);
      break;
    case '}':
      if (stack.empty() || stack.back() != '{')
        return false;
      stack.pop_back();
      break;
    case ']':
      if (stack.empty() || stack.back() != '[')
        return false;
      stack.pop_back();
      break;
    default:
      break;
    }
  }
  return stack.empty() && !in_string;
}

template <typename T> std::string to_json(const T &x) {
  crab::crab_string_os os;
  crab::json::writer w(os);
  crab::json::write(w, x);
  w.finish();
  return compact(os.str());
}

std::string cfg_to_json_string(const z_cfg_t &cfg) {
  crab::crab_string_os os;
  crab::cfg::cfg_to_json(cfg, os);
  return os.str();
}

bool contains(const std::string &haystack, const std::string &needle) {
  return haystack.find(needle) != std::string::npos;
}

} // namespace

// ---------------------------------------------------------------------------
// support/json.hpp
// ---------------------------------------------------------------------------

void test_escaping() {
  BOOST_TEST_EQ(crab::json::escape("plain"), std::string("plain"));
  BOOST_TEST_EQ(crab::json::escape("a\"b"), std::string("a\\\"b"));
  BOOST_TEST_EQ(crab::json::escape("a\\b"), std::string("a\\\\b"));
  BOOST_TEST_EQ(crab::json::escape("a\nb"), std::string("a\\nb"));
  BOOST_TEST_EQ(crab::json::escape("a\tb"), std::string("a\\tb"));
  BOOST_TEST_EQ(crab::json::escape("a\rb"), std::string("a\\rb"));
  BOOST_TEST_EQ(crab::json::escape("a\bb"), std::string("a\\bb"));
  BOOST_TEST_EQ(crab::json::escape("a\fb"), std::string("a\\fb"));
  // Other control characters use the \u00XX form.
  BOOST_TEST_EQ(crab::json::escape(std::string("a\x01"
                                               "b")),
                std::string("a\\u0001b"));
  BOOST_TEST_EQ(crab::json::escape(std::string("\x1f")), std::string("\\u001f"));
}

void test_writer_layout() {
  // A pretty-printed object with a nested compact array.
  crab::crab_string_os os;
  crab::json::writer w(os);
  w.begin_object();
  w.kv_string("name", "foo");
  w.key("succs");
  w.begin_array(true /*compact*/);
  w.value_string("a");
  w.value_string("b");
  w.end_array();
  w.end_object();
  w.finish();
  BOOST_TEST_EQ(os.str(), std::string("{\n  \"name\": \"foo\",\n  \"succs\": "
                                      "[\"a\", \"b\"]\n}\n"));
}

void test_writer_empty_containers() {
  // An empty array/object keeps its brackets on one line: the closing bracket
  // is only put on a fresh line when the frame had items.
  crab::crab_string_os os;
  crab::json::writer w(os);
  w.begin_object();
  w.key("stmts");
  w.begin_array();
  w.end_array();
  w.key("decl");
  w.begin_object();
  w.end_object();
  w.end_object();
  w.finish();
  BOOST_TEST_EQ(os.str(),
                std::string("{\n  \"stmts\": [],\n  \"decl\": {}\n}\n"));
}

void test_writer_scalars() {
  crab::crab_string_os os;
  crab::json::writer w(os);
  w.begin_object(true /*compact*/);
  w.kv_bool("b", true);
  w.kv_unsigned("u", 32);
  w.kv_int("i", -7);
  w.kv_null("n");
  // Numbers are strings so that arbitrary precision survives.
  w.kv_number("big", "170141183460469231731687303715884105727");
  w.end_object();
  w.finish();
  BOOST_TEST_EQ(os.str(),
                std::string("{\"b\": true, \"u\": 32, \"i\": -7, \"n\": null, "
                            "\"big\": "
                            "\"170141183460469231731687303715884105727\"}\n"));
}

// ---------------------------------------------------------------------------
// types/variable_to_json.hpp
// ---------------------------------------------------------------------------

void test_variable_json() {
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var b(vfac["b"], crab::BOOL_TYPE);
  z_var a(vfac["a"], crab::ARR_INT_TYPE);

  BOOST_TEST_EQ(to_json(x),
                std::string("{\"name\":\"x\",\"type\":{\"kind\":\"int\","
                            "\"bitwidth\":32}}"));
  BOOST_TEST_EQ(to_json(b),
                std::string("{\"name\":\"b\",\"type\":{\"kind\":\"bool\"}}"));
  BOOST_TEST_EQ(
      to_json(a),
      std::string("{\"name\":\"a\",\"type\":{\"kind\":\"int_array\"}}"));

  // The bitwidth is part of the type, not of the name.
  z_var x8(vfac["x8"], crab::INT_TYPE, 8);
  BOOST_TEST(contains(to_json(x8), "\"bitwidth\":8"));
}

// ---------------------------------------------------------------------------
// types/linear_constraints_to_json.hpp
// ---------------------------------------------------------------------------

void test_linear_expression_json() {
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  // 2x - y + 3
  z_lin_exp_t e = 2 * x - y + 3;
  BOOST_TEST_EQ(to_json(e),
                std::string("{\"type\":{\"kind\":\"int\",\"bitwidth\":32},"
                            "\"terms\":[[\"2\",\"x\"],[\"-1\",\"y\"]],"
                            "\"const\":\"3\"}"));

  // A constant expression has no variables, hence no type.
  z_lin_exp_t c(ikos::z_number(42));
  BOOST_TEST_EQ(to_json(c),
                std::string("{\"type\":null,\"terms\":[],\"const\":\"42\"}"));
}

void test_linear_constraint_json() {
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  // Crab stores `expr <op> 0`; the exporter moves the constant to the rhs.
  BOOST_TEST_EQ(to_json(z_lin_cst_t(x <= 9)),
                std::string("{\"op\":\"<=\",\"type\":{\"kind\":\"int\","
                            "\"bitwidth\":32},\"terms\":[[\"1\",\"x\"]],"
                            "\"const\":\"9\"}"));
  BOOST_TEST(contains(to_json(z_lin_cst_t(x - y == 0)), "\"op\":\"=\""));
  BOOST_TEST(contains(to_json(z_lin_cst_t(x != 3)), "\"op\":\"!=\""));
  BOOST_TEST(contains(to_json(z_lin_cst_t(x < 5)), "\"op\":\"<\""));

  // 2x - y <= 3
  BOOST_TEST_EQ(to_json(z_lin_cst_t(2 * x - y <= 3)),
                std::string("{\"op\":\"<=\",\"type\":{\"kind\":\"int\","
                            "\"bitwidth\":32},\"terms\":[[\"2\",\"x\"],"
                            "[\"-1\",\"y\"]],\"const\":\"3\"}"));

  // Tautologies and contradictions are nullary: no terms, no const, no type.
  BOOST_TEST_EQ(to_json(z_lin_cst_t::get_true()),
                std::string("{\"op\":\"true\"}"));
  BOOST_TEST_EQ(to_json(z_lin_cst_t::get_false()),
                std::string("{\"op\":\"false\"}"));

  // Arbitrary-precision coefficients survive as decimal strings.
  ikos::z_number big("170141183460469231731687303715884105727");
  BOOST_TEST(contains(to_json(z_lin_cst_t(x <= big)),
                      "\"const\":\"170141183460469231731687303715884105727\""));
}

void test_constraint_system_json() {
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);

  z_lin_cst_sys_t csts;
  csts += (x >= 0);
  csts += (x - y <= 10);
  std::string out = to_json(csts);
  BOOST_TEST(well_formed(out));
  // A conjunction is a bare array of constraints.
  BOOST_TEST(out.front() == '[');
  BOOST_TEST_EQ(std::count(out.begin(), out.end(), '{'), std::ptrdiff_t(4));
}

void test_disjunctive_constraint_system_json() {
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);

  // Bottom and top are tagged, so that they cannot be confused with each other
  // or with an empty conjunction.
  z_dis_lin_cst_sys_t bot(true /*is_false*/);
  BOOST_TEST_EQ(to_json(bot), std::string("{\"kind\":\"false\"}"));

  z_dis_lin_cst_sys_t top(false /*is_false*/);
  BOOST_TEST_EQ(to_json(top), std::string("{\"kind\":\"true\"}"));

  z_lin_cst_sys_t d1;
  d1 += (x <= 0);
  z_lin_cst_sys_t d2;
  d2 += (x >= 10);
  z_dis_lin_cst_sys_t disj(false);
  disj += d1;
  disj += d2;
  std::string out = to_json(disj);
  BOOST_TEST(well_formed(out));
  BOOST_TEST(contains(out, "\"kind\":\"disj\""));
  BOOST_TEST(contains(out, "\"disjuncts\":[["));
  // Two disjuncts, one constraint each.
  BOOST_TEST(contains(out, "]],\"const\":\"0\"}],[{"));
}

// ---------------------------------------------------------------------------
// cfg/cfg_to_json.hpp
// ---------------------------------------------------------------------------

void test_cfg_json_header() {
  variable_factory_t vfac;
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var r(vfac["r"], crab::INT_TYPE, 32);

  crab::cfg::function_decl<ikos::z_number, varname_t> decl("foo", {i}, {r});
  z_cfg_t cfg("entry", "ret", decl);
  BB((&cfg), entry);
  BB((&cfg), ret);
  entry >> ret;
  entry.assign(r, i);

  std::string out = compact(cfg_to_json_string(cfg));
  BOOST_TEST(well_formed(out));
  BOOST_TEST(contains(out, "\"name\":\"foo\""));
  BOOST_TEST(contains(out, "\"entry\":\"entry\""));
  BOOST_TEST(contains(out, "\"exit\":\"ret\""));
  BOOST_TEST(contains(out, "\"declaration\":{\"inputs\":[{\"name\":\"i\""));
  BOOST_TEST(contains(out, "\"outputs\":[{\"name\":\"r\""));
}

void test_cfg_json_without_declaration_or_exit() {
  // Single-argument constructor => no exit block, and no function declaration.
  z_cfg_t cfg("entry");
  BB((&cfg), entry);
  BB((&cfg), sink);
  entry >> sink;

  std::string out = compact(cfg_to_json_string(cfg));
  BOOST_TEST(well_formed(out));
  BOOST_TEST(contains(out, "\"name\":null"));
  BOOST_TEST(contains(out, "\"declaration\":null"));
  BOOST_TEST(contains(out, "\"exit\":null"));
}

void test_cfg_json_block_order_is_sorted() {
  // Blocks are stored in an unordered_map, so the exporter sorts them by label
  // to keep the output stable. Insert them out of order on purpose.
  z_cfg_t cfg("m_entry");
  BB((&cfg), m_entry);
  BB((&cfg), zz);
  BB((&cfg), aa);
  BB((&cfg), kk);
  m_entry >> zz;
  m_entry >> aa;
  m_entry >> kk;

  std::string out = cfg_to_json_string(cfg);
  size_t aa_pos = out.find("\"label\": \"aa\"");
  size_t kk_pos = out.find("\"label\": \"kk\"");
  size_t me_pos = out.find("\"label\": \"m_entry\"");
  size_t zz_pos = out.find("\"label\": \"zz\"");
  BOOST_TEST(aa_pos != std::string::npos);
  BOOST_TEST(kk_pos != std::string::npos);
  BOOST_TEST(me_pos != std::string::npos);
  BOOST_TEST(zz_pos != std::string::npos);
  BOOST_TEST(aa_pos < kk_pos);
  BOOST_TEST(kk_pos < me_pos);
  BOOST_TEST(me_pos < zz_pos);

  // Successors keep the order in which the edges were added.
  BOOST_TEST(contains(out, "\"succs\": [\"zz\", \"aa\", \"kk\"]"));

  // Same CFG, same bytes: the output must not depend on iteration order.
  BOOST_TEST_EQ(out, cfg_to_json_string(cfg));
}

// Every statement kind the exporter supports, in one CFG.
void test_cfg_json_statements() {
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 8);
  z_var b1(vfac["b1"], crab::BOOL_TYPE);
  z_var b2(vfac["b2"], crab::BOOL_TYPE);
  z_var a(vfac["a"], crab::ARR_INT_TYPE);
  z_var a2(vfac["a2"], crab::ARR_INT_TYPE);

  z_cfg_t cfg("entry", "ret");
  BB((&cfg), entry);
  BB((&cfg), ret);
  entry >> ret;

  // integers
  entry.assign(x, 0);
  entry.add(x, x, 1);
  entry.udiv(y, x, 2);
  entry.ashr(y, y, 1);
  entry.assume(x <= 100);
  entry.select(y, x <= 0, 1, 2);
  entry.havoc(y);
  entry.sext(w, x);
  entry.assertion(x != 7, crab::cfg::debug_info("t.c", 10, 4, 99));
  // booleans
  entry.bool_assign(b1, x >= 3);
  entry.bool_assign(b2, b1, true /*is_negated*/);
  entry.bool_and(b1, b1, b2);
  entry.bool_assume(b1);
  entry.bool_select(b2, b1, b1, b2);
  entry.bool_assert(b2);
  // arrays
  entry.array_init(a, 0, 10, 0, 4);
  entry.array_store(a, x, 3, 4, true /*is_strong_update*/);
  entry.array_load(y, a, x, 4);
  entry.array_assign(a2, a);
  // calls
  entry.callsite("callee", {y}, {x});
  ret.unreachable();

  std::string pretty = cfg_to_json_string(cfg);
  BOOST_TEST(well_formed(pretty));
  std::string out = compact(pretty);

  BOOST_TEST(contains(out, "\"stmt\":\"bool_assign_cst\",\"cst_kind\":"
                           "\"linear\","));

  const char *kinds[] = {
      "assign",     "binop",       "assume",          "select",
      "havoc",      "cast",        "assert",          "bool_assign_cst",
      "bool_binop", "bool_assume", "bool_select",     "bool_assert",
      "array_init", "array_store", "array_load",      "array_assign",
      "callsite",   "unreachable", "bool_assign_var"};
  for (auto k : kinds) {
    BOOST_TEST(contains(out, std::string("\"stmt\":\"") + k + "\""));
  }

  // Operator names are spelled out, not printed as symbols.
  BOOST_TEST(contains(out, "\"op\":\"add\""));
  BOOST_TEST(contains(out, "\"op\":\"udiv\""));
  BOOST_TEST(contains(out, "\"op\":\"ashr\""));
  BOOST_TEST(contains(out, "\"op\":\"and\""));
  BOOST_TEST(contains(out, "\"op\":\"sext\""));
  BOOST_TEST(!contains(out, "\"op\":\">>_r\""));

  // A cast records both widths, each next to the operand it belongs to.
  BOOST_TEST(contains(out, "\"rhs_width\":8"));
  BOOST_TEST(contains(out, "\"lhs_width\":32"));

  // Debug info is attached when present and null otherwise.
  BOOST_TEST(contains(
      out, "\"loc\":{\"file\":\"t.c\",\"line\":10,\"col\":4,\"id\":99}"));
  BOOST_TEST(contains(out, "\"loc\":null"));

  // Flags that change the meaning of a statement are exported.
  BOOST_TEST(contains(out, "\"negated\":true"));

  // array_store does not export strong_update, but it does export elem_size,
  // which comes last in both array_init and array_store.
  BOOST_TEST(!contains(out, "\"strong_update\""));
  BOOST_TEST(contains(out, "\"val\":{\"type\":null,\"terms\":[],\"const\":"
                           "\"0\"},\"elem_size\":"));
  BOOST_TEST(contains(out, "\"value\":{\"type\":null,\"terms\":[],\"const\":"
                           "\"3\"},\"elem_size\":"));
  // In array_load it comes right after lhs instead.
  BOOST_TEST(contains(out, "\"stmt\":\"array_load\",\"array\":{\"name\":"
                           "\"a\","));
  BOOST_TEST(contains(out, "\"elem_size\":{\"type\":null,\"terms\":[],"
                           "\"const\":\"4\"},\"lhs\":{\"name\":\"y\","));

  // A callsite records callee, lhs and args.
  BOOST_TEST(contains(out, "\"callee\":\"callee\""));
}

// Region and reference statements, added one kind at a time. Every kind the
// exporter does not handle yet still raises CRAB_ERROR, so this test grows as
// support does.
void test_cfg_json_region_statements() {
  variable_factory_t vfac;
  z_var mem(vfac["region_0"], crab::REG_INT_TYPE, 32);
  z_var mem2(vfac["region_1"], crab::REG_INT_TYPE, 32);
  z_var unk(vfac["region_2"], crab::REG_UNKNOWN_TYPE);

  z_cfg_t cfg("entry");
  BB((&cfg), entry);
  entry.region_init(mem);
  entry.region_init(mem2);
  entry.region_copy(mem2, mem);
  entry.region_cast(unk, mem);

  std::string out = compact(cfg_to_json_string(cfg));
  BOOST_TEST(well_formed(out));
  // The written region is the lhs, like every other statement.
  BOOST_TEST(contains(out, "\"stmt\":\"region_init\",\"lhs\":{\"name\":"
                           "\"region_0\",\"type\":{\"kind\":\"int_region\","
                           "\"bitwidth\":32}}"));
  BOOST_TEST(contains(out, "\"stmt\":\"region_copy\",\"rhs\":{\"name\":"
                           "\"region_0\",\"type\":{\"kind\":\"int_region\","
                           "\"bitwidth\":32}},\"lhs\":{\"name\":\"region_1\""));
  // region_cast repeats the operand types: a cast is precisely the statement
  // where the two differ.
  BOOST_TEST(contains(out, "\"stmt\":\"region_cast\",\"rhs\":{\"name\":"
                           "\"region_2\",\"type\":{\"kind\":\"unknown_region\"}"
                           "},\"rhs_type\":{\"kind\":\"unknown_region\"},"
                           "\"lhs\":{\"name\":\"region_0\",\"type\":{\"kind\":"
                           "\"int_region\",\"bitwidth\":32}},\"lhs_type\":"
                           "{\"kind\":\"int_region\",\"bitwidth\":32}"));
}

// Reference statements. The convention across all of them is that source
// operands are printed before destination operands.
void test_cfg_json_reference_statements() {
  variable_factory_t vfac;
  crab::tag_manager as_man;
  z_var mem(vfac["region_0"], crab::REG_INT_TYPE, 32);
  z_var mem2(vfac["region_1"], crab::REG_INT_TYPE, 32);
  z_var p(vfac["p"], crab::REF_TYPE);
  z_var q(vfac["q"], crab::REF_TYPE);
  z_var r(vfac["r"], crab::REF_TYPE);
  z_var v(vfac["v"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  z_var c(vfac["c"], crab::BOOL_TYPE);
  z_var_or_cst_t size4(ikos::z_number(4), crab::variable_type(crab::INT_TYPE, 32));

  z_cfg_t cfg("entry");
  BB((&cfg), entry);
  entry.region_init(mem);
  entry.region_init(mem2);
  entry.make_ref(p, mem, size4, as_man.mk_tag());
  entry.store_to_ref(p, mem, z_var_or_cst_t(v));
  entry.load_from_ref(v, p, mem);
  entry.gep_ref(q, mem2, p, mem, 8);
  entry.select_ref(r, mem2, c, p, mem, q, mem2);
  entry.assume_ref(z_ref_cst_t::mk_not_null(p));
  entry.assert_ref(z_ref_cst_t::mk_eq(p, q, ikos::z_number(4)));
  entry.int_to_ref(n, mem, r);
  entry.ref_to_int(mem, r, n);
  entry.remove_ref(mem, p);
  entry.bool_assign(c, z_ref_cst_t::mk_not_null(q));

  std::string out = compact(cfg_to_json_string(cfg));
  BOOST_TEST(well_formed(out));

  // make_ref: size and allocation site are inputs; the new reference and the
  // region it lives in are the outputs, so both take the lhs prefix.
  BOOST_TEST(contains(out, "\"stmt\":\"make_ref\",\"size\":{\"kind\":"
                           "\"const\",\"value\":\"4\",\"type\":{\"kind\":"
                           "\"int\",\"bitwidth\":32}},\"alloc_site\":"));
  BOOST_TEST(contains(out, "\"lhs_region\":{\"name\":\"region_0\","));
  // remove_ref and load/store put the region before the reference.
  BOOST_TEST(contains(out, "\"stmt\":\"remove_ref\",\"region\":{\"name\":"
                           "\"region_0\",\"type\":{\"kind\":\"int_region\","
                           "\"bitwidth\":32}},\"ref\":{\"name\":\"p\","
                           "\"type\":{\"kind\":\"ref\"}}"));
  BOOST_TEST(contains(out, "\"stmt\":\"store_to_ref\",\"region\":{\"name\":"
                           "\"region_0\","));
  BOOST_TEST(contains(out, "\"value\":{\"kind\":\"var\",\"name\":\"v\","));
  BOOST_TEST(contains(out, "\"stmt\":\"load_from_ref\",\"region\":{\"name\":"
                           "\"region_0\","));
  // gep_ref: rhs side first, then the offset, then the lhs side.
  BOOST_TEST(contains(out, "\"stmt\":\"gep_ref\",\"rhs_region\":{\"name\":"
                           "\"region_0\","));
  BOOST_TEST(contains(out, "\"rhs\":{\"name\":\"p\",\"type\":{\"kind\":\"ref\"}"
                           "},\"offset\":"));
  BOOST_TEST(contains(out, "\"lhs_region\":{\"name\":\"region_1\","));
  // select_ref: cond and both branches first, destination last.
  BOOST_TEST(contains(out, "\"stmt\":\"select_ref\",\"cond\":{\"name\":\"c\","));
  BOOST_TEST(contains(out, "\"left_region\":{\"name\":\"region_0\","));
  BOOST_TEST(contains(out, "\"right_region\":{\"name\":\"region_1\","));
  // int_to_ref / ref_to_int.
  // The two conversions are mirrors: the reference operand carries the _ref
  // suffix, on whichever side the statement writes it.
  BOOST_TEST(contains(out, "\"stmt\":\"int_to_ref\",\"rhs\":{\"name\":"
                           "\"n\",\"type\":{\"kind\":\"int\",\"bitwidth\":"
                           "32}},\"lhs_region\":{\"name\":\"region_0\","));
  BOOST_TEST(contains(out, "\"lhs_ref\":{\"name\":\"r\",\"type\":{\"kind\":"
                           "\"ref\"}}"));
  BOOST_TEST(contains(out, "\"stmt\":\"ref_to_int\",\"rhs_region\":{\"name\":"
                           "\"region_0\","));
  BOOST_TEST(contains(out, "\"rhs_ref\":{\"name\":\"r\",\"type\":{\"kind\":"
                           "\"ref\"}},\"lhs\":{\"name\":\"n\","));
  // bool_assign_cst tags which flavour of constraint its rhs is.
  BOOST_TEST(contains(out, "\"stmt\":\"bool_assign_cst\",\"cst_kind\":"
                           "\"reference\",\"rhs\":{\"op\":\"!=\","
                           "\"lhs\":{\"name\":\"q\","));
}

// Reference constraints, as carried by assume_ref and assert_ref.
void test_reference_constraint_json() {
  variable_factory_t vfac;
  z_var p(vfac["p"], crab::REF_TYPE);
  z_var q(vfac["q"], crab::REF_TYPE);

  // A unary constraint compares against null, which is why there is no "rhs":
  // its absence is what distinguishes `p != null` from `p != q`.
  BOOST_TEST_EQ(to_json(z_ref_cst_t::mk_not_null(p)),
                std::string("{\"op\":\"!=\",\"lhs\":{\"name\":\"p\",\"type\":"
                            "{\"kind\":\"ref\"}}}"));
  BOOST_TEST_EQ(to_json(z_ref_cst_t::mk_null(p)),
                std::string("{\"op\":\"==\",\"lhs\":{\"name\":\"p\",\"type\":"
                            "{\"kind\":\"ref\"}}}"));
  // A binary constraint carries the other reference and the offset.
  BOOST_TEST_EQ(to_json(z_ref_cst_t::mk_eq(p, q, ikos::z_number(4))),
                std::string("{\"op\":\"==\",\"lhs\":{\"name\":\"p\",\"type\":"
                            "{\"kind\":\"ref\"}},\"rhs\":{\"name\":\"q\","
                            "\"type\":{\"kind\":\"ref\"}},\"offset\":\"4\"}"));
  BOOST_TEST(contains(to_json(z_ref_cst_t::mk_lt_null(p)), "\"op\":\"<\""));
  BOOST_TEST(contains(to_json(z_ref_cst_t::mk_ge_null(p)), "\"op\":\">=\""));
  // Tautologies and contradictions are nullary, like their linear counterpart.
  BOOST_TEST_EQ(to_json(z_ref_cst_t::mk_true()),
                std::string("{\"op\":\"true\"}"));
  BOOST_TEST_EQ(to_json(z_ref_cst_t::mk_false()),
                std::string("{\"op\":\"false\"}"));
}

// The value partitioning intrinsics are directives with no concrete semantics
// and are the only statements the exporter is allowed to drop.
void test_cfg_json_skips_value_partition_intrinsics() {
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);

  z_cfg_t cfg("entry");
  BB((&cfg), entry);
  entry.intrinsic("value_partition_start", {}, {z_var_or_cst_t(x)});
  entry.havoc(x);
  entry.intrinsic("value_partition_end", {}, {z_var_or_cst_t(x)});

  std::string out = compact(cfg_to_json_string(cfg));
  BOOST_TEST(well_formed(out));
  BOOST_TEST(!contains(out, "value_partition"));
  // Exactly one statement survives, and the array is still well formed: the
  // skipped statements must not leave a dangling separator behind.
  BOOST_TEST(contains(out, "\"stmts\":[{\"stmt\":\"havoc\""));
  BOOST_TEST_EQ(std::count(out.begin(), out.end(), '['),
                std::count(out.begin(), out.end(), ']'));
}

// A block with no statements and no successors still emits both keys.
void test_cfg_json_empty_block() {
  z_cfg_t cfg("entry");
  cfg.insert("entry");
  std::string out = compact(cfg_to_json_string(cfg));
  BOOST_TEST(well_formed(out));
  BOOST_TEST(contains(out, "\"stmts\":[]"));
  BOOST_TEST(contains(out, "\"succs\":[]"));
}

// Labels are JSON strings, so a label with characters that must be escaped
// has to come back out escaped.
void test_cfg_json_escapes_labels() {
  z_cfg_t cfg("entry");
  BB((&cfg), entry);
  auto &odd = cfg.insert("a\"b\\c");
  entry >> odd;

  std::string out = cfg_to_json_string(cfg);
  BOOST_TEST(well_formed(out));
  BOOST_TEST(contains(out, "\"label\": \"a\\\"b\\\\c\""));
  BOOST_TEST(contains(out, "\"succs\": [\"a\\\"b\\\\c\"]"));
}

// The exporter must see the CFG as the analyses do, i.e. after simplify().
void test_cfg_json_after_simplify() {
  variable_factory_t vfac;
  z_var x(vfac["x"], crab::INT_TYPE, 32);

  z_cfg_t cfg("entry", "ret");
  BB((&cfg), entry);
  BB((&cfg), middle);
  BB((&cfg), ret);
  entry >> middle;
  middle >> ret;
  entry.assign(x, 0);
  middle.add(x, x, 1);

  size_t before = cfg.size();
  cfg.simplify();
  BOOST_TEST(cfg.size() < before);

  std::string out = compact(cfg_to_json_string(cfg));
  BOOST_TEST(well_formed(out));
  // "middle" was merged away, so it must not appear as a block.
  BOOST_TEST(!contains(out, "\"label\":\"middle\""));
  // Its statement survives in the merged block.
  BOOST_TEST(contains(out, "\"stmt\":\"binop\""));
}

int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool /*stats_enabled*/) -> int {
    test_escaping();
    test_writer_layout();
    test_writer_empty_containers();
    test_writer_scalars();

    test_variable_json();

    test_linear_expression_json();
    test_linear_constraint_json();
    test_constraint_system_json();
    test_disjunctive_constraint_system_json();

    test_cfg_json_header();
    test_cfg_json_without_declaration_or_exit();
    test_cfg_json_block_order_is_sorted();
    test_cfg_json_statements();
    test_cfg_json_region_statements();
    test_cfg_json_reference_statements();
    test_reference_constraint_json();
    test_cfg_json_skips_value_partition_intrinsics();
    test_cfg_json_empty_block();
    test_cfg_json_escapes_labels();
    test_cfg_json_after_simplify();

    return boost::report_errors();
  });
}
