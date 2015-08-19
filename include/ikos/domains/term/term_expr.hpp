#ifndef __IKOS_TERM_EXPR_H__
#define __IKOS_TERM_EXPR_H__
#include <map>
#include <memory>
#include <boost/optional.hpp>

using namespace std;
using namespace boost;

// ERGH, make it properly flexible later.
namespace ikos {
  namespace term {
    typedef int term_id;
    typedef int var_id;

    enum term_kind { TERM_CONST, TERM_VAR, TERM_APP };

    template<class Num, class Ftor>
    class term : public writeable {
    public:
      typedef term<Num, Ftor> term_t;
      typedef std::shared_ptr<term_t> term_ptr;

      virtual term_kind kind(void) = 0;

      bool operator<(term_t& other);

      virtual ostream& write(ostream& o) = 0;

      int depth;
    };

    template<class Num, class Ftor>
    class term_ref {
    public:
      typedef term<Num, Ftor> term_t;
      typedef std::shared_ptr<term_t> term_ptr;

      term_ref(term_ptr _p)
        : p(_p)
      { }

      term_ref(term_t* _p)
        : p(_p)
      { }

      bool operator<(const term_ref& other) const
      {
        return *p < *(other.p);
      }

      term_ptr p;
    };

    // Dispatch for different kinds of terms.
    template<class Num, class Ftor>
    class const_term : public term<Num, Ftor> {
    public:
      const_term(Num _val)
        : val(_val)
      { }
      term_kind kind(void) { return TERM_CONST; }

      ostream& write(ostream& o)
      {
        return (o << "c(" << val << ")");
      }
      Num val;
    };

    template<class Num, class Ftor>
    class var_term : public term<Num, Ftor> {
    public:
      var_term(var_id _var)
        : var(_var)
      { }
      term_kind kind(void) { return TERM_VAR; }
      ostream& write(ostream& o)
      {
        return (o << "v(" << var << ")");
      }
      var_id var;
    };
    
    template<class Num, class Ftor>
    class ftor_term : public term<Num, Ftor> {
    public:
      ftor_term(Ftor _f, std::vector<term_id>& _args)
        : ftor(_f), args(_args)
      { }
      term_kind kind(void) { return TERM_APP; }
      ostream& write(ostream& o)
      {
        o << "<op>(";
        bool first = true;
        for(term_id c : args)
        {
          if(first)
            first = false;
          else
            o << ", ";
          o << c;
        }
        return (o << ")");
      }

      Ftor ftor;
      std::vector<term_id> args;
    };

    // These functions are super-unsafe.
    template<class Num, class Ftor>
    var_id term_var(term<Num, Ftor>* term)
    {
      return static_cast< var_term<Num, Ftor>* >(term)->var;
    }
    template<class Num, class Ftor>
    Num term_const(term<Num, Ftor>* term)
    {
      return static_cast< const_term<Num, Ftor>* >(term)->val;
    }

    template<class Num, class Ftor>
    Ftor term_ftor(term<Num, Ftor>* term)
    {
      return static_cast< ftor_term<Num, Ftor>* >(term)->ftor;
    }

    template<class Num, class Ftor>
    std::vector<term_id>& term_args(term<Num, Ftor>* term)
    {
      return static_cast< ftor_term<Num, Ftor>* >(term)->args;
    }

    // Term ordering
    template<class Num, class Ftor> 
    bool term<Num, Ftor>::operator<(term_t& other)
    {
      if(kind() != other.kind())
        return kind() < other.kind();

      switch(kind())
      {
        case TERM_CONST:
          {
            return term_const(this) < term_const(&other);
          }
          break;
        case TERM_VAR:
          {
            return term_var(this) < term_var(&other);
          }
        case TERM_APP:
          {
            if(term_ftor(this) != term_ftor(&other))
              return term_ftor(this) < term_ftor(&other);

            std::vector<int>& xargs(term_args(this));
            std::vector<int>& yargs(term_args(&other));

            if(xargs.size() != yargs.size())
              return xargs.size() < yargs.size();

            for(unsigned int ii = 0; ii< xargs.size(); ii++)
            {
              if(xargs[ii] != yargs[ii])
                return xargs[ii] < yargs[ii];
            }
            return false;
          }
        default:
          assert(0 && "Unsupported term kind.");
          return false;
      }
    }

    template<class Num, class Ftor>
    class term_table : public writeable {
    public:
      typedef term_table<Num, Ftor> term_table_t;
      typedef term<Num, Ftor> term_t;
      typedef std::shared_ptr<term_t> term_ptr;
      typedef term_ref<Num, Ftor> term_ref_t;
      typedef term_id term_id_t;

      typedef const_term<Num, Ftor> const_term_t;
      typedef var_term<Num, Ftor> var_term_t;
      typedef ftor_term<Num, Ftor> ftor_term_t;

      // For establishing a mapping between tables
      typedef std::map<term_id, term_id> term_map_t;
      typedef std::map<std::pair<term_id, term_id>, term_id> gener_map_t;
               
      term_table(void)
        : free_var(0)
      { }

      term_table(const term_table_t& o)
        : free_var(o.free_var), _map(o._map),
          terms(o.terms), _parents(o._parents), // term_refs(o.term_refs),
          _depth(o._depth),
          free_terms(o.free_terms)
      { }

      term_table_t& operator=(const term_table_t& o)
      {
        free_var = o.free_var;
        _map = o._map;
        terms = o.terms;
//        term_refs = o.term_refs;
        _parents = o._parents;
        _depth = o._depth;
        free_terms = o.free_terms;
        return *this;
      }

      term_id make_const(const Num& n) {
        term_ref_t ref(new const_term_t(n));
        return add_term(ref);
      };
      optional<term_id> find_const(Num& n) {
        term_ref_t ref(new const_term_t(n));
        return find_term(ref);
      }

      term_id make_var(var_id& n) {
        free_var = std::max(free_var, n+1);
        term_ref_t ref(new var_term_t(n));
        return add_term(ref);
      };

      term_id fresh_var(void) {
        var_id v = free_var++;
        return make_var(v); 
      }

      template <typename T, typename ... Types>
      void collect_args(std::vector<T>& vec, T& x, Types ... rest)
      {
        vec.push_back(x); collect_args(vec, rest...);
      }

      template <typename T>
      void collect_args(std::vector<T>& vec) { }

      term_id apply_ftor(const Ftor& f, std::vector<term_id>& ids)
      {
        term_ref_t ref(new ftor_term_t(f, ids));
        return add_term(ref);
      }

      template <typename ... Types>
      term_id apply_ftor(const Ftor& f, Types ... args)
      {
        std::vector<term_id> ids;
        collect_args(ids, args...);
        return apply_ftor(f, ids);
      }

      optional<term_id> find_ftor(Ftor& f, std::vector<term_id>& ids)
      {
        term_ref_t ref(new ftor_term_t(f, ids));
        return find_term(ref); 
      }

      template<typename ... Types>
      optional<term_id> find_ftor(Ftor& f, Types ... args)
      {
        std::vector<term_id> ids;
        collect_args(ids, args...);
        return find_ftor(f, ids);
      }


      term_t* get_term_ptr(term_id t)
      {
        return terms[t].p.get();
      }

      void add_ref(term_id t)
      {
        // term_refs[t]++;
      }

      void deref(term_id t, std::vector<term_id>& forgotten)
      {
//        fprintf(stdout, "WARNING: term_table::deref not properly implemented.");
        return;

        /*
        assert(term_refs[t]);
        term_refs[t]--;
        if(!term_refs[t])
        {
          forgotten.push_back(t);
          free_terms.push_back(t);
          term_ref_t ref(terms[t]);
          _map.erase(ref);
          if(ref.p.get()->kind() == TERM_APP)
          {
            for(term_id c : term_args(ref.p.get()))
              deref(c, forgotten);
          }
        }
        */
      }

      // Check if a tx is a generalization of ty, given an existing context.
      bool map_leq(term_table_t& y, term_id tx, term_id ty, term_map_t& map)
      {
        auto it = map.find(ty);
        if(it != map.end())
          return (*it).second == tx;

        term_ref_t ry(y.terms[ty]);
        term_ref_t rx(terms[tx]);
        switch(ry.p.get()->kind())
        {
          case TERM_CONST:
            {
              if(rx.p.get()->kind() != TERM_CONST || term_const(rx.p.get()) != term_const(ry.p.get()))
                return false;
            }
            break;
          case TERM_APP:
            {
              if(rx.p.get()->kind() != TERM_APP || term_ftor(rx.p.get()) != term_ftor(ry.p.get()))
                return false;

              std::vector<int>& xargs(term_args(rx.p.get()));
              std::vector<int>& yargs(term_args(ry.p.get()));

              if(xargs.size() != yargs.size())
                return false;

              for(unsigned int ii = 0; ii< xargs.size(); ii++)
              {
                if(!map_leq(y, xargs[ii], yargs[ii], map))
                  return false;
              }
            }
            break;
         case TERM_VAR:
            break;
        }
        map[ty] = tx;
        return true;
      }

      term_id_t generalize(term_table_t& y, term_id tx, term_id ty, term_table_t& out, gener_map_t& g_map)
      {
        auto txy = std::make_pair(tx, ty); 
        auto it = g_map.find(txy);
        if(it != g_map.end())
        {
          return (*it).second;
        } else {
          // Haven't found this pair yet.
          // Generalize the arguments.
          auto px(terms[tx].p.get());
          auto py(y.terms[ty].p.get());

          term_kind kx = px->kind();
          term_kind ky = py->kind();

          // If either is a variable, we just create a fresh variable
          term_id_t ret;
          do
          {
            if(kx == TERM_VAR || ky == TERM_VAR || kx != ky)
            {
              ret = out.fresh_var(); break;
            }

            if(kx == TERM_CONST)
            {
              if(term_const(px) != term_const(py))
              {
                ret = out.fresh_var();
              } else {
                ret = out.make_const(term_const(px));
              }
            } else {
              assert(kx == TERM_APP);
              if(term_ftor(px) != term_ftor(py))
              {
                ret = out.fresh_var(); break;
              }

              std::vector<term_id>& xargs(term_args(px));
              std::vector<term_id>& yargs(term_args(py));
              if(xargs.size() != yargs.size())
              {
                ret = out.fresh_var(); break;
              }

              std::vector<term_id> xyargs;
              for(unsigned int ii = 0; ii < xargs.size(); ii++)
                xyargs.push_back(generalize(y, xargs[ii], yargs[ii], out, g_map));
              ret = out.apply_ftor(term_ftor(px), xyargs);
            }
          } while(0);
          
          g_map[txy] = ret;
          return ret;
        }
      }
      
      std::vector<term_id_t>& parents(term_id_t id)
      {
        return _parents[id];
      }

      ostream& write(ostream& o) { 
        bool first = true;
        for(int ti = 0; ti < terms.size(); ti++)
        {
          if(first)
            first = false;
          else
            o << ", ";
          term_t* p(terms[ti].p.get());
          o << ti << " -> " << *p;
        }
        return o;
      }

      int size() { return terms.size(); }

      int depth(term_id t) { return _depth[t]; }

    protected:
      optional<term_id> find_term(term_ref_t ref)
      {
        auto it = _map.find(ref);
        if(it != _map.end())
          return optional<term_id>(it->second);
        else
          return optional<term_id>();
      }

      // When we know that a germ doesn't already exist.
      term_id fresh_term(term_ref_t ref)
      {
        if(free_terms.size() > 0)
        {
          term_id t = free_terms.back();
          free_terms.pop_back();
          terms[t] = ref;
          _parents[t].clear();
          _depth[t] = 0;
//          term_refs[t] = 0;
          return t;
        } else {
//          term_id t = term_refs.size();
//          term_refs.push_back(0);
          term_id t = terms.size();
          terms.push_back(ref);
          _parents.push_back(std::vector<term_id>());
          _depth.push_back(0);
          return t;
        }
      }

      term_id add_term(term_ref_t ref)
      {
        auto it = _map.find(ref);
        if(it != _map.end())
        {
          return (*it).second;
        } else {
          term_id id = fresh_term(ref);
          _map[ref] = id;
          if(ref.p.get()->kind() == TERM_APP)
          {
            unsigned int c_depth = 0;
            for(term_id c : term_args(ref.p.get()))
            {
              // add_ref(c);
              _parents[c].push_back(id);
              c_depth = std::max(c_depth, _depth[c]);
            }
            _depth[id] = 1+c_depth;
          }
          assert(_map.size() == id+1);
          return id;
        }
      }

      int free_var;
      std::map<term_ref_t, term_id> _map;
      std::vector<term_ref_t> terms;
      std::vector< std::vector<term_id_t> > _parents;
      std::vector<unsigned int> _depth;
      std::vector<term_id> free_terms;
    };
  }
}
#endif
