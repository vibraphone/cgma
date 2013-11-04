#ifndef MANAGED_PTR_VECTOR
#define MANAGED_PTR_VECTOR

#include <vector>
#include <memory>
#include <cstddef>

//! A vector of pointers to objects that are automatically deleted when the container itself is deleted.
/*! Loosely
modeled after the boost ptr containers.
*/
template <typename X>
class ManagedPtrVector
{
public:
  typedef X value_type;
  typedef X& reference_type;
  typedef size_t size_type;
  typedef std::vector<X*> container_type;
  
  class iterator
  {
  public:
    typedef iterator this_type;
    typedef X& reference;
    typedef size_t size_type;
    typedef size_t difference_type;
    
    iterator()
      {}
    iterator(const this_type& rhs)
        :mIter(rhs.mIter)
      {}
    iterator(const typename ManagedPtrVector<X>::container_type::iterator& rhs)
        : mIter(rhs)
      {}
    ~iterator()
      {}
    
    X* operator->()
      { return *mIter; }
    reference operator*() const
      { return **mIter; }
    
    bool operator==(const this_type& rhs) const
      { return this->mIter == rhs.mIter; }
    bool operator!=(const this_type& rhs) const
      { return this->mIter != rhs.mIter; }
    
    this_type& operator++()
      {
        ++mIter;
        return *this;
      }
    
    this_type operator++(int)
      {
        this_type rv = *this;
        ++mIter;
        return rv;
      }
    
    this_type& operator--()
      {
        --mIter;
        return *this;
      }
    
    this_type operator--(int)
      {
        this_type rv = *this;
        --mIter;
        return rv;
      }

    this_type operator+(difference_type n)
    {
      this_type rv = *this;
      rv += n;
      return rv;
    }
    
    this_type& operator+=(difference_type n)
      {
        mIter += n;
        return *this;
      }
    
    this_type& operator-=(difference_type n)
      {
        mIter -= n;
        return *this;
      }
    
    reference operator[](difference_type i) const
      {
        return reference(*(mIter+i));
      }
    
  private:
    typename ManagedPtrVector<X>::container_type::iterator mIter;
  };
  
  
  ManagedPtrVector()
  {}
  ~ManagedPtrVector()
  { clear(); }
  
  // Add an object to the container.
  // The container becomes "owned" by the container.
  void push_back(X* obj)
  { mContainer.push_back(obj); }

  iterator begin()
  { return iterator(mContainer.begin()); }
  iterator end()
  { return iterator(mContainer.end()); }

  //! Refer by index
  reference_type operator[](size_type i)
  { return *(mContainer[i]); }

  //! Remove an item from the list without deleting it.
  /*! The returned auto_ptr now owns the pointer.*/
  std::auto_ptr<X> release(iterator to_release)
  {
    // save a raw pointer.
    X* rv = to_release.operator->();
    // find the iterator in the container.  We don't have access
    // to the internal iterator, so we have to loop through to find it.
    // This could probably be optimized.
    typename container_type::iterator i = mContainer.begin();
    for (;
      i != mContainer.end() && iterator(i) != to_release;
      ++i)
    {}
    // we either found the iterator, or the end.  erase it.
    mContainer.erase(i);
    return std::auto_ptr<X>(rv);
  }

  //! Replace one object with another.
  /*! The returned auto_ptr now owns the pointer removed from
      the container, and the container owns the object pointed to
      by \a new_val.
  */
  std::auto_ptr<X> replace(iterator to_replace, std::auto_ptr<X> new_val)
  {
    // save a raw pointer.
    X* rv = to_replace.operator->();
    // find the iterator in the container.  We don't have access
    // to the internal iterator, so we have to loop through to find it.
    // This could probably be optimized.
    typename container_type::iterator i = mContainer.begin();
    for (;
      i != mContainer.end();
      ++i)
    {
      if (iterator(i) == to_replace)
      {
        *i = new_val;
      }
    }
    return std::auto_ptr<X>(NULL);
  }

  // Delete all objects contained in this container.
  void clear()
  {
    // delete all the pointers
    for (typename container_type::iterator i = mContainer.begin();
      i != mContainer.end();
      i++)
    {
      delete *i;
    }
    mContainer.clear();
  }

  size_type size() const
  { return mContainer.size(); }

private:
  // Make both of these illegal
  ManagedPtrVector(const ManagedPtrVector<X>&);
  ManagedPtrVector& operator=(const ManagedPtrVector<X>&);

  container_type mContainer;
};



#endif
