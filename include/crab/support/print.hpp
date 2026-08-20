#pragma once

/*
 * A small facility for printing containers (and optional-like values) to a
 * crab_os, so that write()/operator<< implementations do not need to
 * hand-roll separator loops.
 *
 * Everything lives in namespace crab::print. Code inside namespace crab
 * writes print::seq(v); downstream code writes crab::print::seq(v), or
 * abbreviates with `namespace cp = crab::print;`.
 *
 * Quick reference (o is a crab_os):
 *
 *   o << print::seq(v);                      // [a, b, c]
 *   o << print::seq(v, print::fmt_set());    // {a, b, c}
 *   o << print::seq(v).sep("; ").bare();     // a; b; c
 *   o << print::kv(m);                       // {k1 -> v1; k2 -> v2}
 *   o << print::opt(x);                      // value of *x, or "none"
 *   o << print::opt(x, "<absent>");          // custom empty marker
 *
 *   print::print_range(o, v, print::fmt_seq_tight());   // [a,b,c]
 *   print::print_range_with(o, v, fn);       // fn(crab_os&, const elem&)
 *   print::to_string(x);                     // print anything to a std::string
 *
 *   print::separator sep(o, ", ");           // for loops that filter elements
 *   for (auto &e : v) {                      // or span several sibling loops
 *     if (skip(e)) continue;
 *     sep.next() << e;
 *   }
 *
 * Before/after -- the typical hand-rolled separator loop
 *
 *   o << "{";
 *   for (auto it = m.begin(); it != m.end();) {
 *     o << it->first << " -> " << it->second;
 *     if (++it != m.end()) { o << "; "; }
 *   }
 *   o << "}";
 *
 * becomes the one-liner
 *
 *   o << print::kv(m);
 *
 * Elements are printed by trying, in order: operator<< (so anything already
 * streamable prints exactly as before), std::pair (first -> second),
 * container (recurse), optional-like (dereference or print the empty
 * marker), then write(crab_os&). See print_element below for the rationale
 * of that order.
 *
 * Design notes:
 *
 * - There is deliberately NO generic operator<<(crab_os&, const T&).
 *   crab_os's operator<< overloads are non-template members, and bool, float,
 *   short, enums and T* reach them only by promotion/conversion; a generic
 *   const T& template would be an identity match and silently steal them.
 *   The seq/kv/opt wrappers sidestep this: their operator<< is a plain
 *   non-template hidden friend that can only ever match the wrapper.
 *
 * - Raw pointers are NEVER dereferenced automatically. A raw T* keeps
 *   binding to crab_os::operator<<(const void*) and prints as an address,
 *   as it always has. If a caller wants the pointee they must say so
 *   explicitly with print::opt(p).
 *
 * - This header is C++14: no if constexpr, no std::void_t, no fold
 *   expressions, hence the make_void / priority_tag machinery below.
 */

#include <crab/support/os.hpp>

#include <string>
#include <type_traits>
#include <utility>

namespace crab {
namespace print {

/*
 * Describes ONE level of container printing: brackets, separator, the arrow
 * between the members of a pair, what to print for an empty container, the
 * marker for an empty optional-like element, and whether the separator also
 * trails the last element (the "{a;b;c;}" debug style).
 *
 * Nested elements do not inherit the enclosing format: a nested container is
 * printed with its own type default (fmt_map() if its elements are pairs,
 * fmt_seq() otherwise). So map<K, vector<V>> prints {k -> [a, b], ...}.
 *
 * The setters return a modified copy, so formats chain:
 *   format().brackets("{", "}").sep("; ")
 */
class format {
public:
  const char *m_open = "[";
  const char *m_close = "]";
  const char *m_sep = ", ";
  const char *m_arrow = " -> ";
  // printed instead of m_open/m_close when the range is empty (nullptr means
  // print m_open followed by m_close)
  const char *m_empty = nullptr;
  // marker for an empty optional-like element. "none" matches the existing
  // convention (e.g. cfg.hpp); it is not "_|_", which in crab means bottom,
  // a lattice value -- an absent value is not a bottom value.
  const char *m_none = "none";
  bool m_trailing = false;

  format sep(const char *s) const {
    format f(*this);
    f.m_sep = s;
    return f;
  }
  format brackets(const char *open, const char *close) const {
    format f(*this);
    f.m_open = open;
    f.m_close = close;
    return f;
  }
  format bare() const { return brackets("", ""); }
  format arrow(const char *a) const {
    format f(*this);
    f.m_arrow = a;
    return f;
  }
  format empty(const char *e) const {
    format f(*this);
    f.m_empty = e;
    return f;
  }
  format none(const char *n) const {
    format f(*this);
    f.m_none = n;
    return f;
  }
  format trailing(bool b = true) const {
    format f(*this);
    f.m_trailing = b;
    return f;
  }
};

/*
 * Presets, named after the conventions they reproduce so that migrating an
 * existing printer is a lookup rather than a judgement call.
 */
// [a, b, c]
inline format fmt_seq() { return format(); }
// [a,b,c]
inline format fmt_seq_tight() { return format().sep(","); }
// {a, b, c}
inline format fmt_set() { return format().brackets("{", "}"); }
// {a,b,c}
inline format fmt_set_tight() { return fmt_set().sep(","); }
// {k1 -> v1; k2 -> v2} -- the patricia-tree map printers
inline format fmt_map() { return format().brackets("{", "}").sep("; "); }
// same bytes as fmt_map(); named separately so patricia call sites document
// which convention they came from. NOTE: patricia-tree iterators yield a
// binding type that is NOT std::pair (it only has first/second members), so
// the automatic pair rank does not fire for them -- print a patricia map
// with print_range_with and an explicit "k -> v" callback, passing this
// format (see region_domain.hpp for an example).
inline format fmt_patricia() { return fmt_map(); }
// [|v1, v2|] -- the graph printers
inline format fmt_graph() { return format().brackets("[|", "|]"); }
// {a;b;c;} -- the CRAB_LOG debug loops
inline format fmt_debug() {
  return format().brackets("{", "}").sep(";").trailing();
}

/*
 * A stateful joiner for the loops that cannot become a single range call:
 * loops that `continue`-filter elements, or two sibling loops sharing one
 * "first" flag (e.g. the split_dbm printer). The first call to next() prints
 * nothing; every later call prints the separator. Either way it returns the
 * stream, so the natural use is:
 *
 *   separator sep(o, ", ");
 *   for (...) { if (filtered) continue; sep.next() << elem; }
 */
class separator {
  crab_os &m_o;
  const char *m_sep;
  bool m_first;

public:
  separator(crab_os &o, const char *sep = ", ")
      : m_o(o), m_sep(sep), m_first(true) {}
  crab_os &next() {
    if (!m_first) {
      m_o << m_sep;
    }
    m_first = false;
    return m_o;
  }
  // true iff next() has not fired yet -- for "print a footer only if we
  // printed anything" call sites
  bool first() const { return m_first; }
  void reset() { m_first = true; }
};

template <typename Range>
void print_range(crab_os &o, const Range &r, const format &fmt);

namespace detail {

// C++14 has no std::void_t. Hand-rolled through a struct rather than a
// direct alias: with a direct alias some compilers do not guarantee that
// unused alias arguments trigger SFINAE (CWG 1558).
template <typename...> struct make_void { using type = void; };
template <typename... Ts> using void_t = typename make_void<Ts...>::type;

// Rank ordering for the dispatch ladder below: an argument of type
// priority_tag<N> prefers an overload taking priority_tag<N> over one taking
// priority_tag<N-1> (derived-to-base conversion), so ranks are tried
// highest-first and SFINAE'd-out ranks fall through.
template <unsigned N> struct priority_tag : priority_tag<N - 1> {};
template <> struct priority_tag<0> {};

template <typename T> struct always_false : std::false_type {};

template <typename T> struct is_pair : std::false_type {};
template <typename A, typename B>
struct is_pair<std::pair<A, B>> : std::true_type {};

// Ranges are detected through MEMBER begin()/end() only, never std::begin:
// with std::begin every string literal (const char[N]) would count as a
// range. std::string is excluded explicitly, and the reason is nesting, not
// ambiguity: it has member begin()/end(), so without the exclusion
// vector<string> would print [[h,i],[t,h,e,r,e]] instead of [hi, there].
template <typename T, typename = void>
struct has_member_begin_end : std::false_type {};
template <typename T>
struct has_member_begin_end<T,
                            void_t<decltype(std::declval<const T &>().begin()),
                                   decltype(std::declval<const T &>().end())>>
    : std::true_type {};

template <typename T>
struct is_range
    : std::integral_constant<bool, has_member_begin_end<T>::value &&
                                       !std::is_same<T, std::string>::value> {};

template <typename T, typename = void>
struct has_os_insert : std::false_type {};
template <typename T>
struct has_os_insert<
    T, void_t<decltype(std::declval<crab_os &>() << std::declval<const T &>())>>
    : std::true_type {};

// Probes T&, not const T&: a few types declare write(crab_os&) non-const
// (e.g. fwd_analyzer), and a const-only probe would miss them. A const
// write() is still callable through a T&, so this probe accepts both.
template <typename T, typename = void> struct has_write : std::false_type {};
template <typename T>
struct has_write<
    T, void_t<decltype(std::declval<T &>().write(std::declval<crab_os &>()))>>
    : std::true_type {};

template <typename T, typename = void> struct has_deref : std::false_type {};
template <typename T>
struct has_deref<T, void_t<decltype(*std::declval<const T &>())>>
    : std::true_type {};

template <typename T, typename = void>
struct is_bool_testable : std::false_type {};
template <typename T>
struct is_bool_testable<
    T, void_t<decltype(static_cast<bool>(std::declval<const T &>()))>>
    : std::true_type {};

// Structural detection of boost::optional, std::shared_ptr, std::unique_ptr
// (and std::optional under C++17) without including any of their headers:
// dereferenceable and contextually convertible to bool. Ranges keep range
// semantics, and raw pointers are excluded on purpose -- see the "never
// auto-dereference a raw pointer" note at the top of this file.
template <typename T>
struct is_optional_like
    : std::integral_constant<
          bool, has_deref<T>::value && is_bool_testable<T>::value &&
                    !is_range<T>::value && !std::is_pointer<T>::value> {};

// For a nested range printed as an element: its own type default.
template <typename Range> struct range_element {
  using type = typename std::decay<
      decltype(*std::declval<const Range &>().begin())>::type;
};

template <typename Range> inline format default_format() {
  return is_pair<typename range_element<Range>::type>::value ? fmt_map()
                                                             : fmt_seq();
}

/*
 * The element dispatch ladder. print_element(o, x, fmt) prints one element
 * by trying, highest rank first:
 *
 *   operator<< > pair > range > optional-like > write() > static_assert
 *
 * operator<< is the TOP rank on purpose: an element that can already be
 * streamed keeps streaming exactly as it always has, so routing it through
 * the facility can never change its printed bytes. This matters because
 * several crab types have BOTH iterators and an operator<< (e.g.
 * linear_constraint_t iterates its terms); ranking ranges higher would
 * silently print those as bracketed element lists.
 *
 * operator<< is also preferred over write(): the two agree for crab's own
 * types, but templates like separate_domains are instantiated by downstream
 * clients with their types, where they can differ, and a type's own
 * operator<< is the author's stated intent.
 *
 * Caveat of the top rank: has_os_insert matches through implicit
 * conversions, so an optional-like type with an IMPLICIT operator bool
 * would stream as 0/1 instead of dereferencing. boost::optional,
 * shared_ptr and unique_ptr all declare theirs explicit, so they take the
 * optional-like rank as intended.
 *
 * The fmt argument supplies this level's pair arrow and optional-empty
 * marker only; a nested range starts over with its own type default.
 */
template <typename T>
void print_element(crab_os &o, const T &x, const format &fmt);

template <typename T,
          typename std::enable_if<has_os_insert<T>::value, int>::type = 0>
void print_element_impl(crab_os &o, const T &x, const format &,
                        priority_tag<4>) {
  o << x;
}

template <typename A, typename B>
void print_element_impl(crab_os &o, const std::pair<A, B> &p, const format &fmt,
                        priority_tag<3>) {
  print_element(o, p.first, fmt);
  o << fmt.m_arrow;
  print_element(o, p.second, fmt);
}

template <typename T,
          typename std::enable_if<is_range<T>::value, int>::type = 0>
void print_element_impl(crab_os &o, const T &r, const format &,
                        priority_tag<2>) {
  print_range(o, r, default_format<T>());
}

template <typename T,
          typename std::enable_if<is_optional_like<T>::value, int>::type = 0>
void print_element_impl(crab_os &o, const T &x, const format &fmt,
                        priority_tag<1>) {
  if (x) {
    print_element(o, *x, fmt);
  } else {
    o << fmt.m_none;
  }
}

template <typename T,
          typename std::enable_if<has_write<T>::value, int>::type = 0>
void print_element_impl(crab_os &o, const T &x, const format &,
                        priority_tag<0>) {
  // const_cast because a few write(crab_os&) methods are missing a const
  // qualifier (see has_write above). This rank REQUIRES that write() does
  // not actually mutate: calling a genuinely mutating write() on an object
  // that was defined const would be undefined behavior. A write() that
  // mutates its object is broken for printing anyway; fix its signature
  // rather than relying on this rank.
  const_cast<T &>(x).write(o);
}

template <typename T> struct is_printable_element {
  static constexpr bool value = is_pair<T>::value || is_range<T>::value ||
                                is_optional_like<T>::value ||
                                has_os_insert<T>::value || has_write<T>::value;
};

template <typename T>
void print_element(crab_os &o, const T &x, const format &fmt) {
  static_assert(is_printable_element<T>::value,
                "crab::print: this element type is not printable. It needs "
                "an operator<<(crab_os&, const T&) or a write(crab_os&) "
                "method; alternatively use print_range_with with an explicit "
                "callback.");
  print_element_impl(o, x, fmt, priority_tag<4>{});
}

} // namespace detail

template <typename Range>
void print_range(crab_os &o, const Range &r, const format &fmt) {
  auto it = r.begin();
  auto end = r.end();
  if (it == end && fmt.m_empty) {
    o << fmt.m_empty;
    return;
  }
  o << fmt.m_open;
  for (bool first = true; it != end; ++it, first = false) {
    // in the trailing style the separator is emitted after every element
    // instead of between elements
    if (!first && !fmt.m_trailing) {
      o << fmt.m_sep;
    }
    detail::print_element(o, *it, fmt);
    // the "{a;b;c;}" style: the separator also trails the last element
    if (fmt.m_trailing) {
      o << fmt.m_sep;
    }
  }
  o << fmt.m_close;
}

template <typename Range> void print_range(crab_os &o, const Range &r) {
  print_range(o, r, fmt_seq());
}

template <typename It>
void print_range(crab_os &o, It first, It last, const format &fmt) {
  if (first == last && fmt.m_empty) {
    o << fmt.m_empty;
    return;
  }
  o << fmt.m_open;
  for (It it = first; it != last; ++it) {
    if (it != first && !fmt.m_trailing) {
      o << fmt.m_sep;
    }
    detail::print_element(o, *it, fmt);
    if (fmt.m_trailing) {
      o << fmt.m_sep;
    }
  }
  o << fmt.m_close;
}

template <typename It> void print_range(crab_os &o, It first, It last) {
  print_range(o, first, last, fmt_seq());
}

/*
 * Like print_range, but each element is printed by fn(o, elem) instead of
 * the dispatch ladder. fn takes the stream (matching the write(crab_os&)
 * convention), so the same callback works against a crab_string_os.
 *
 * Use this where the element printing must be pinned down explicitly, e.g.
 * converting a printer whose Key/Value are client-supplied template
 * parameters: an explicit write-calling callback is byte-exact by
 * construction, independent of what operator<<s the client type has.
 */
template <typename Range, typename Fn>
void print_range_with(crab_os &o, const Range &r, Fn fn, const format &fmt) {
  auto it = r.begin();
  auto end = r.end();
  if (it == end && fmt.m_empty) {
    o << fmt.m_empty;
    return;
  }
  o << fmt.m_open;
  for (bool first = true; it != end; ++it, first = false) {
    // in the trailing style the separator is emitted after every element
    // instead of between elements
    if (!first && !fmt.m_trailing) {
      o << fmt.m_sep;
    }
    fn(o, *it);
    if (fmt.m_trailing) {
      o << fmt.m_sep;
    }
  }
  o << fmt.m_close;
}

template <typename Range, typename Fn>
void print_range_with(crab_os &o, const Range &r, Fn fn) {
  print_range_with(o, r, fn, fmt_seq());
}

template <typename It, typename Fn>
void print_range_with(crab_os &o, It first, It last, Fn fn, const format &fmt) {
  if (first == last && fmt.m_empty) {
    o << fmt.m_empty;
    return;
  }
  o << fmt.m_open;
  for (It it = first; it != last; ++it) {
    if (it != first && !fmt.m_trailing) {
      o << fmt.m_sep;
    }
    fn(o, *it);
    if (fmt.m_trailing) {
      o << fmt.m_sep;
    }
  }
  o << fmt.m_close;
}

template <typename It, typename Fn>
void print_range_with(crab_os &o, It first, It last, Fn fn) {
  print_range_with(o, first, last, fn, fmt_seq());
}

// Print anything the ladder accepts (including the seq/kv/opt wrappers)
// to a std::string, folding the usual crab_string_os dance.
template <typename T> std::string to_string(const T &x) {
  crab_string_os os;
  detail::print_element(os, x, format());
  return os.str();
}

/*
 * The ergonomic wrappers: o << print::seq(v), o << print::kv(m),
 * o << print::opt(x).
 *
 * seq_ref/kv_ref/opt_ref versus the print_range/print_range_with functions
 * above: one implementation, two entry styles. seq_ref::write() simply
 * delegates to print_range(), so both print the same bytes. Use the wrapper
 * form in operator<< chains and as CRAB_ERROR/CRAB_WARN arguments -- a
 * wrapper is a streamable value. Use the function forms inside larger
 * write() bodies: they print immediately, and print_range_with takes an
 * explicit per-element callback for when element printing must be pinned
 * down (see its comment above).
 *
 * Each holds a POINTER to the wrapped object, not a reference: CRAB_ERROR /
 * CRAB_WARN funnel their arguments through ___print___ by value, and a
 * reference member would make those copies dangle-prone; a pointer member
 * copies safely. The wrappers are still view types -- they must not outlive
 * the wrapped object, so use them directly in a print expression rather
 * than storing them.
 *
 * operator<< is a hidden friend: a non-template function that only ADL can
 * find and that can only ever match the wrapper, so it cannot interfere
 * with crab_os's own non-template operator<< overloads.
 */

// A range printed as a sequence: o << print::seq(v) -> [a, b, c]
template <typename Range> class seq_ref {
  const Range *m_r;
  format m_fmt;

public:
  seq_ref(const Range &r, const format &fmt) : m_r(&r), m_fmt(fmt) {}
  seq_ref sep(const char *s) const { return seq_ref(*m_r, m_fmt.sep(s)); }
  seq_ref brackets(const char *open, const char *close) const {
    return seq_ref(*m_r, m_fmt.brackets(open, close));
  }
  seq_ref bare() const { return seq_ref(*m_r, m_fmt.bare()); }
  seq_ref empty(const char *e) const { return seq_ref(*m_r, m_fmt.empty(e)); }
  seq_ref trailing(bool b = true) const {
    return seq_ref(*m_r, m_fmt.trailing(b));
  }
  void write(crab_os &o) const { print_range(o, *m_r, m_fmt); }
  friend crab_os &operator<<(crab_os &o, const seq_ref &s) {
    s.write(o);
    return o;
  }
};

template <typename Range> seq_ref<Range> seq(const Range &r) {
  return seq_ref<Range>(r, fmt_seq());
}
template <typename Range>
seq_ref<Range> seq(const Range &r, const format &fmt) {
  return seq_ref<Range>(r, fmt);
}

// A range of pairs printed as a map: o << print::kv(m) -> {k1 -> v1; k2 -> v2}
template <typename Map> class kv_ref {
  const Map *m_m;
  format m_fmt;

public:
  kv_ref(const Map &m, const format &fmt) : m_m(&m), m_fmt(fmt) {}
  kv_ref sep(const char *s) const { return kv_ref(*m_m, m_fmt.sep(s)); }
  kv_ref brackets(const char *open, const char *close) const {
    return kv_ref(*m_m, m_fmt.brackets(open, close));
  }
  kv_ref bare() const { return kv_ref(*m_m, m_fmt.bare()); }
  kv_ref arrow(const char *a) const { return kv_ref(*m_m, m_fmt.arrow(a)); }
  kv_ref empty(const char *e) const { return kv_ref(*m_m, m_fmt.empty(e)); }
  void write(crab_os &o) const { print_range(o, *m_m, m_fmt); }
  friend crab_os &operator<<(crab_os &o, const kv_ref &m) {
    m.write(o);
    return o;
  }
};

template <typename Map> kv_ref<Map> kv(const Map &m) {
  return kv_ref<Map>(m, fmt_map());
}
template <typename Map> kv_ref<Map> kv(const Map &m, const format &fmt) {
  return kv_ref<Map>(m, fmt);
}

// An optional-like (or raw pointer, made explicit here) printed as its
// value, or a marker when empty: o << print::opt(x) -> *x or "none"
template <typename T> class opt_ref {
  const T *m_x;
  const char *m_none;

public:
  opt_ref(const T &x, const char *none) : m_x(&x), m_none(none) {}
  void write(crab_os &o) const {
    static_assert(detail::has_deref<T>::value &&
                      detail::is_bool_testable<T>::value,
                  "crab::print::opt requires an optional-like type "
                  "(dereferenceable and contextually convertible to bool) or "
                  "a raw pointer");
    if (static_cast<bool>(*m_x)) {
      detail::print_element(o, **m_x, format().none(m_none));
    } else {
      o << m_none;
    }
  }
  friend crab_os &operator<<(crab_os &o, const opt_ref &x) {
    x.write(o);
    return o;
  }
};

template <typename T> opt_ref<T> opt(const T &x, const char *none = "none") {
  return opt_ref<T>(x, none);
}

} // namespace print
} // namespace crab
