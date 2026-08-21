#pragma once

/**
 * A minimal streaming JSON writer.
 *
 * Crab has no JSON dependency and this is deliberately not a JSON *parser* nor
 * a DOM: it is just enough to emit well-formed, pretty-printed JSON from the
 * various `*_to_json` writers (see cfg_to_json.hpp,
 * linear_constraints_to_json.hpp).
 *
 * Usage:
 *
 *   crab::json::writer w(crab::outs());
 *   w.begin_object();
 *   w.key("name"); w.value_string("foo");
 *   w.key("blocks"); w.begin_array();
 *      w.value_string("entry");
 *   w.end_array();
 *   w.end_object();
 *
 * Arrays and objects may be opened in *compact* mode, in which case no newlines
 * or indentation are emitted inside them. That keeps small tuples such as
 * `["1", "x"]` on a single line.
 *
 * Numbers: JSON numbers are IEEE doubles for most consumers, which cannot
 * represent Crab's arbitrary-precision numbers. Every writer in Crab therefore
 * emits numeric values as *strings* (see value_number). Consumers must parse
 * them with a bignum type.
 */

#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <string>
#include <vector>

namespace crab {
namespace json {

/** Escape a string as a JSON string body (without the surrounding quotes). */
inline std::string escape(const std::string &s) {
  std::string res;
  res.reserve(s.size());
  for (char c : s) {
    switch (c) {
    case '"':
      res += "\\\"";
      break;
    case '\\':
      res += "\\\\";
      break;
    case '\b':
      res += "\\b";
      break;
    case '\f':
      res += "\\f";
      break;
    case '\n':
      res += "\\n";
      break;
    case '\r':
      res += "\\r";
      break;
    case '\t':
      res += "\\t";
      break;
    default:
      if (static_cast<unsigned char>(c) < 0x20) {
        static const char *digits = "0123456789abcdef";
        res += "\\u00";
        res += digits[(static_cast<unsigned char>(c) >> 4) & 0xf];
        res += digits[static_cast<unsigned char>(c) & 0xf];
      } else {
        res += c;
      }
    }
  }
  return res;
}

class writer {
  struct frame {
    bool is_object;
    bool has_items;
    bool compact;
  };

  crab::crab_os &m_os;
  unsigned m_indent_width;
  std::vector<frame> m_stack;
  // True right after key(): the next value is emitted in place rather than on
  // a fresh line.
  bool m_pending_key;
  // Number of enclosing compact frames; while > 0 no newlines are emitted.
  unsigned m_compact_depth;

  bool compact() const { return m_compact_depth > 0; }

  void indent(size_t depth) {
    for (size_t i = 0; i < depth * m_indent_width; ++i) {
      m_os << " ";
    }
  }

  /** Emit the separator and indentation preceding the next item. */
  void next_item() {
    if (m_pending_key) {
      m_pending_key = false;
      return;
    }
    if (m_stack.empty()) {
      return;
    }
    if (m_stack.back().has_items) {
      m_os << ",";
      if (compact()) {
        m_os << " ";
      }
    }
    if (!compact()) {
      m_os << "\n";
      indent(m_stack.size());
    }
    m_stack.back().has_items = true;
  }

  void open(bool is_object, bool is_compact) {
    next_item();
    m_os << (is_object ? "{" : "[");
    m_stack.push_back(frame{is_object, false, is_compact});
    if (is_compact) {
      m_compact_depth++;
    }
  }

  void close(bool is_object) {
    if (m_stack.empty()) {
      CRAB_ERROR("json::writer: unbalanced close");
    }
    frame f = m_stack.back();
    if (f.is_object != is_object) {
      CRAB_ERROR("json::writer: mismatched close");
    }
    // Whether *this* frame was laid out compactly, which is what decides
    // where its closing bracket goes. Must be read before popping the frame's
    // own contribution to m_compact_depth.
    const bool was_compact = compact();
    m_stack.pop_back();
    if (f.compact) {
      m_compact_depth--;
    }
    if (f.has_items && !was_compact) {
      m_os << "\n";
      indent(m_stack.size());
    }
    m_os << (is_object ? "}" : "]");
  }

public:
  explicit writer(crab::crab_os &os, unsigned indent_width = 2)
      : m_os(os), m_indent_width(indent_width), m_pending_key(false),
        m_compact_depth(0) {}

  writer(const writer &) = delete;
  writer &operator=(const writer &) = delete;

  void begin_object(bool is_compact = false) { open(true, is_compact); }
  void end_object() { close(true); }
  void begin_array(bool is_compact = false) { open(false, is_compact); }
  void end_array() { close(false); }

  void key(const std::string &k) {
    next_item();
    m_os << "\"" << escape(k) << "\": ";
    m_pending_key = true;
  }

  void value_string(const std::string &s) {
    next_item();
    m_os << "\"" << escape(s) << "\"";
  }

  /**
   * Numeric values are emitted as JSON *strings* so that arbitrary-precision
   * numbers survive consumers whose JSON numbers are doubles.
   */
  void value_number(const std::string &digits) { value_string(digits); }

  void value_bool(bool b) {
    next_item();
    m_os << (b ? "true" : "false");
  }

  void value_unsigned(unsigned long long n) {
    next_item();
    m_os << n;
  }

  void value_int(long long n) {
    next_item();
    m_os << n;
  }

  void value_null() {
    next_item();
    m_os << "null";
  }

  /** Emit an already-formatted JSON fragment verbatim. */
  void value_raw(const std::string &json) {
    next_item();
    m_os << json;
  }

  // -- convenience -----------------------------------------------------------

  void kv_string(const std::string &k, const std::string &v) {
    key(k);
    value_string(v);
  }
  void kv_number(const std::string &k, const std::string &digits) {
    key(k);
    value_number(digits);
  }
  void kv_bool(const std::string &k, bool v) {
    key(k);
    value_bool(v);
  }
  void kv_unsigned(const std::string &k, unsigned long long v) {
    key(k);
    value_unsigned(v);
  }
  void kv_int(const std::string &k, long long v) {
    key(k);
    value_int(v);
  }
  void kv_null(const std::string &k) {
    key(k);
    value_null();
  }

  /** Finish the document with a trailing newline. */
  void finish() {
    if (!m_stack.empty()) {
      CRAB_ERROR("json::writer: unbalanced document");
    }
    m_os << "\n";
  }
};

/** Stringify anything with a crab_os operator<<. */
template <typename T> inline std::string to_string(const T &x) {
  crab::crab_string_os os;
  os << x;
  return os.str();
}

} // end namespace json
} // end namespace crab
