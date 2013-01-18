
#ifndef CUBIT_HASHER
#define CUBIT_HASHER

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>

// Hasher class optimized for:
//   simple integer keys
//   fairly contiguous keys (not memory efficient for non-contiguous keys)
// provides constant time lookups, insertions, and removals
//   with emphasis on fast lookups
// T is sub container type, B is number of hashing bits
template <class Y, unsigned int B>
class Hasher
{
public:

  struct IterTag : public std::forward_iterator_tag {};
  template <class T> struct BaseIter
  {
    typedef IterTag iterator_category;
    typedef int difference_type;
    typedef T value_type;
    typedef T* pointer;
    typedef T& reference;
  };
  
  class iterator;
  typedef iterator const_iterator;

  typedef unsigned int Key;
  typedef typename Y::value_type value_type;

  // number of bits we represent
  enum { NumberOfBits = Y::NumberOfBits + B };

  // number of bits our sub containers manage
  enum { SubSize = Y::NumberOfBits };

  // number of sub holders we have
  enum { MySize = 1 << B };

  // helper mask
  enum { MyMask = ( 1 << SubSize ) - 1 };

  value_type Null;

  // constructor
  Hasher() : Null(0)
  {
    for(Key i=0; i<MySize; i++)
      mChildren[i] = 0;
    mCount=0;
  }
  
  // destructor
  ~Hasher()
  {
    for(Key i=0; i<MySize; i++)
      if(mChildren[i])
        delete mChildren[i];
  }


  // add a value for a key
  bool add(Key idx, value_type ptr)
  {
    assert(ptr != 0);
    Key myidx = my_index(idx);
    Y* y = mChildren[myidx];
    if(!y)
    {
      y = mChildren[myidx] = new Y();
    }
    bool inserted = y->add(sub_index(idx), ptr);
    if(inserted)
    {
      this->mCount++;
    }
    return inserted;
  }

  // remove a key
  bool remove(Key idx)
  {
    Key myidx = my_index(idx);
    Y* y = mChildren[myidx];
    bool removed = y ? y->remove(sub_index(idx)) : false;
    if(removed)
    {
      this->mCount--;
      if(y->count() == 0)
      {
        delete mChildren[myidx];
        mChildren[myidx] = 0;
      }
    }
    return removed;
  }

  // lookup a value
  const value_type& operator[](Key idx) const
  {
    Key myidx = my_index(idx);
    if(mChildren[myidx])
      return (*mChildren[myidx])[sub_index(idx)];
    return this->Null;
  }
  
  value_type& operator[](Key idx)
  {
    Key myidx = my_index(idx);
    if(mChildren[myidx])
      return (*mChildren[myidx])[sub_index(idx)];
    return this->Null;
  }
  
  Key count() const
  {
    return mCount;
  }

  iterator begin()
  {
    return iterator(this->mChildren, true);
  }
  
  iterator end()
  {
    return iterator(this->mChildren, false);
  }

  class iterator : public BaseIter< typename Y::value_type >
  {
      friend class Hasher;
      typedef typename Y::value_type& reference;
      enum { MySize = 1 << B };
    public:

      iterator() : mChildren(NULL), index(MySize) {}
      iterator(const iterator& copy)
      {
        this->operator=(copy);
      }

      // prefix incrementer
      iterator& operator++()
      {
        if(index != MySize)
        {
          ++iter;
          typename Y::iterator e = mChildren[index]->end();
          for(; this->iter != e && !*iter; ++iter) {}
          if(iter == e)
          {
            for(++index; index < MySize && !mChildren[index]; ++index) {}
            if(index != MySize)
            {
              e = mChildren[index]->end();
              for(iter = mChildren[index]->begin(); this->iter != e && !*iter; ++iter) {}
            }
          }
        }
        return *this;
      }

      // postfix incrementer
      iterator operator++(int)
      {
        iterator tmp(*this);
        this->operator++();
        return tmp;
      }
      
      bool operator==(const iterator& other) const
      {
        if( MySize == this->index && MySize == other.index)
          return true;
        return this->index == other.index && this->iter == other.iter;
      }

      bool operator!=(const iterator& other) const
      {
        return !operator==(other);
      }

      reference operator*()
      {
        return *iter;
      }

      iterator& operator=(const iterator& other)
      {
        this->mChildren = other.mChildren;
        this->index = other.index;
        if(other.mChildren)
          this->iter = other.iter;
        return *this;
      }
      
    protected:
      iterator(Y** c, bool b)
      {
        mChildren = c;
        index = MySize;
        if(b)
        {
          Y** s = c;
          for(; s<c+MySize && !*s; s++){}
          if(s<c+MySize && *s)
          {
            index = s - c;
            iter = (*s)->begin();
            for(; 0 == *iter; ++iter) {}
          }
        }
      }

      Y** mChildren;
      int index;
      typename Y::iterator iter;
  };


  // given key, get index for which subcontainer
  Key my_index(Key idx) const
  {
    return idx >> SubSize;
  }

  // given key, get index for subcontainer
  Key sub_index(Key idx) const
  {
    return idx & MyMask;
  }

private:

  Key mCount;
  Y* mChildren[MySize];
};


// vector class used as the real storage for the hasher
template <class Y, unsigned int S>
class HasherVector : public std::vector<Y>
{
  public:
    typedef unsigned int Key;
    typedef typename std::vector<Y>::value_type value_type;
    enum { NumberOfBits = S };

    HasherVector() 
        : std::vector<Y>(1<<S, (Y)0), mCount(0) 
      {
      }

    bool add(Key idx, value_type& ptr)
    {
      Y& p = (*this)[idx];
      if(!p)
      {
        p = ptr;
        ++this->mCount;
        return true;
      }
      return false;
    }

    bool remove(Key idx)
    {
      Y& p = (*this)[idx];
      if(p)
      {
        p = 0;
        --this->mCount;
        return true;
      }
      return false;
    }

    Key count() const
    {
      return mCount;
    }

  private:
    Key mCount;
};

#if 0
int main(int, char**)
{

  Hasher <10, HasherVector<int> > myHash;

  for(int i = 1; i < 5000000; i++)
  {
    myHash.add(i,i);
  }

  //printf("finished...\n");
}


  typedef Hasher<5, Hasher<5, HasherVector<int*> > > MyHash;
  MyHash mytree(15);

  //mytree.add(1,(int*)(0x01));
  
  //std::vector<int*> vec;


  //HasherVector<int*> myhasher(5, (int*)0x055);
  //HasherVector<int*>::iterator jter = myhasher.begin();
  //for(; jter != myhasher.end(); ++jter)
  //{
  //  int* val = *jter;
  //}

  for(int i=0; i<10000; i++)
  {
    //unsigned int key = rand();
     unsigned int key = i;
    int* value = (int*)rand();
    bool added = mytree.add(key, value);
    if(!added)
      printf("wasn't added\n");
    for(int j=0; j<1000; j++)
      mytree[key];
    if(value != mytree[key])
      printf("error %i != %i\n", (int)value, (int)mytree[key]);
    //mytree.remove(key);
  }
  
  MyHash::iterator iter = mytree.begin();
  for(; iter != mytree.end(); ++iter)
  {
    int* val = *iter;
  }
}
#endif

#endif // #ifndef CUBIT_HASHER

