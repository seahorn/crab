#pragma once

namespace crab {
// This is pretty similar to std::reference_wrapper but it allows
// default constructor.
template<class T> class reference_wrapper {
public:
  using type = reference_wrapper<T>;
  
  reference_wrapper(): m_p(nullptr) {}
  reference_wrapper(T &p): m_p(&p) {}
  
  reference_wrapper(const type &other) = default;    
  type &operator=(const type &other) = default;
  reference_wrapper(type &&other) = default;    
  type &operator=(type &&other) = default;
  
  T& get() const {
    assert(m_p);
    return *m_p;
  }
private:
  T* m_p;
};
} // end namespace
