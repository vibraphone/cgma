//-------------------------------------------------------------------------
// Purpose       : Simple array template
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------

#ifndef VG_ARRAY_HPP
#define VG_ARRAY_HPP

#include <assert.h>

template <class T> class VGArray
{
  public:
  
    VGArray( int initial_size = 0 );
    
    ~VGArray();
    
    inline int size() const;
    // get current array size
    
    void size( int new_size );
    // increase or decrease size of array
    
    void size_end( int new_size );
    // increase or decrease size, populating the 
    // new array begining with the last element at the
    // previous size in the last slot of the new size,
    // and filling the array in decreasing index order.
    
    inline void push( const T& entry );
    // increase array size by 1 and add entry at last index
    
    inline const T& pop();
    // decrease array size by 1 and return removed entry
    
    inline T& operator[]( int index );
    
    inline const T& operator[]( int index ) const;
    
    void remove( int index );
    // remove entry at index and shift all higher indices down
    
    void insert( const T& entry, int index );
    // insert entry at specified index, shifting all higher indices up
    
    int find( const T& entry, int search_from = 0 ) const;
    // find the index of an occurance of the passed entry
    // in the array.  search starts from search_from and
    // ends at the last index.  -1 is returned if no match
    // is found.
    
    void reverse();
    // Reverse order of list.
    
  private:
  
    const VGArray& operator=( const VGArray<T>& );
    VGArray( const VGArray<T>& );
      //- do not allow assignment
  
    static int storage( int size );
      //- calculate storage to allocate for specified size
      //- (smallest power-2 greater than size)
  
    T* data_;
    int storage_;
    int size_;
};

#include "VGArray.cpp"

template <class T> inline
int VGArray<T>::size() const
  { return size_; }

template <class T> inline
T& VGArray<T>::operator[]( int index )
{
  assert( index >= 0 && index < size_ );
  return data_[index];
}

template <class T> inline
const T& VGArray<T>::operator[]( int index ) const
{
  assert( index >= 0 && index < size_ );
  return data_[index];
}
  
template <class T> inline
void VGArray<T>::push( const T& entry )
{
  int s = size_;
  size( s + 1 );
  data_[s] = entry;
}

#endif
