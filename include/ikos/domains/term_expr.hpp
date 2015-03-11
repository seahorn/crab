#ifndef __IKOS_TERM_EXPR_H__
#define __IKOS_TERM_EXPR_H__
#include <map>
#include <memory>
#include <boost/optional.hpp>

// ERGH, make it properly flexible later.
namespace ikos {
  namespace term {
    typedef int term_id;
    typedef int var_id;

    enum term_kind { TERM_CONST, TERM_VAR, TERM_APP };

    template<class Num, class Ftor>
    class term {
    public:
      typedef term<Num, Ftor> term_t;
      typedef std::shared_ptr<term_t> term_ptr;

      virtual term_kind kind(void) = 0;
      bool operator<(term_t& other);
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

      /*
      term_t& operator*(void)
      {
        return *p; 
      }
      */

      /*
      term_t* operator->(void) {
        return (term_t*) p;
      }
      operator term_t*(void) { return (term_t*) p; }
      */

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
      Num val;
    };

    template<class Num, class Ftor>
    class var_term : public term<Num, Ftor> {
    public:
      var_term(var_id _var)
        : var(_var)
      { }
      term_kind kind(void) { return TERM_VAR; }
      var_id var;
    };
    
    template<class Num, class Ftor>
    class ftor_term : public term<Num, Ftor> {
    public:
      ftor_term(Ftor _f, std::vector<term_id>& _args)
        : ftor(_f), args(_args)
      { }
      term_kind kind(void) { return TERM_APP; }
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
    class term_table {
    public:
      typedef term<Num, Ftor> term_t;
      typedef std::shared_ptr<term_t> term_ptr;
      typedef term_ref<Num, Ftor> term_ref_t;
      typedef term_id term_id_t;

      typedef const_term<Num, Ftor> const_term_t;
      typedef var_term<Num, Ftor> var_term_t;
      typedef ftor_term<Num, Ftor> ftor_term_t;

      term_id term_const(Num& n) {
        term_ref_t ref(new const_term_t(n));
        return add_term(ref);
      };

      term_id term_var(var_id& n) {
        free_var = std::max(free_var, n+1);
        term_ref_t ref(new var_term_t(n));
        return add_term(ref);
      };

      term_id fresh_var(void) {
        var_id v = free_var++;
        return term_var(v); 
      }
      
      template <typename T>
      void collect_args(std::vector<T>& vec) { }

      template <typename T, typename ... Types>
      void collect_args(std::vector<T>& vec, T& x, Types ... rest)
      {
        vec.push_back(x); collect_args(vec, rest...);
      }

      term_id apply_ftor(Ftor& f, std::vector<term_id>& ids)
      {
        term_ref_t ref(new ftor_term_t(f, ids));
        return add_term(ref);
      }

      template <typename ... Types>
      term_id apply_ftor(Ftor& f, Types ... args)
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

      void add_ref(term_id t)
      {
        term_refs[t]++;
      }

      void deref(term_id t, std::vector<term_id>& forgotten)
      {
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
      }
    protected:
      optional<term_id> find_term(term_ref_t ref)
      {
        auto it = _map.find(ref);
        if(it != _map.end())
          return optional<term_id>(*it);
        else
          return optional<term_id>();
      }

      term_id fresh_term(void)
      {
        if(free_terms.size() > 0)
        {
          term_id t = free_terms.back();
          free_terms.pop_back();
          return t;
        } else {
          term_id t = term_refs.size();
          term_refs.push_back(0);
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
          term_id id = fresh_term();
          _map[ref] = id;
          if(ref.p.get()->kind() == TERM_APP)
          {
            for(term_id c : term_args(ref.p.get()))
              add_ref(c);
          }
          assert(_map.size() == id+1);
          return id;
        }
      }

      int free_var;
      std::map<term_ref_t, term_id> _map;
      std::vector<term_ref_t> terms;
      std::vector<unsigned int> term_refs;
      std::vector<term_id> free_terms;
    };
  }
}
#endif
