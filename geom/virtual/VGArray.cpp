
#include <algorithm>
#include "VGArray.hpp"
#include "CubitDefines.h"

template <class T>
VGArray<T>::VGArray( int init_size )
{
  if( init_size <= 0 )
  {
    storage_ = 0;
    size_ = 0;
    data_ = 0;
  }
  else
  {
    size_ = init_size;
    storage_ = storage( init_size );
    data_ = new T[storage_];
  }
}

template <class T>
VGArray<T>::~VGArray( )
{
  delete [] data_;
  storage_ = size_ = 0;
  data_ = 0;
}

template <class T>
void VGArray<T>::size( int new_size )
{
  if( new_size <= 0 ) 
  {
    new_size = 0;
  }
  else if( new_size > storage_ )
  {
    int new_stor = storage(new_size);
    assert( new_stor > storage_ );
    T* new_array = new T[new_stor];
    if( data_ )
    {
      memcpy( new_array, data_, size_ * sizeof(T) );
      delete [] data_;
    }
    data_ = new_array;
    storage_ = new_stor;
  }
  size_ = new_size;
}

template <class T>
void VGArray<T>::size_end( int new_size )
{
  if( new_size <= 0 ) 
  {
    new_size = 0;
  }
  else if( new_size > storage_ )
  {
    int new_stor = storage(new_size);
    assert( new_stor > storage_ );
    T* new_array = new T[new_stor];
    if( data_ )
    {
      memcpy( new_array + new_size - size_, data_, size_ * sizeof(T) );
      delete [] data_;
    }
    data_ = new_array;
    storage_ = new_stor;
  }
  else if (new_size > size_)
  {
    T* read = data_ + size_;
    T* write = data_ + new_size;
    while ( read > data_ )
      *(--write) = *(--read);
  }
  else if (new_size < size_)
  {
    T* read = data_ + size_ - new_size;
    T* write = data_;
    T* end = data_ + size_;
    while ( read < end )
      *(write--) = *(read--);
  }
  size_ = new_size;
}

template <class T>
int VGArray<T>::storage( int size )
{
  assert( size <= 0x40000000 ); //no power of 2 greater than size (overflow)

  int result;
  for( result = 1; result < size; result = result << 1 );
  assert(result>=size);
  return result;
}

template <class T>
void VGArray<T>::remove( int index )
{
  assert( index >= 0 && index < size_ );
  size_--;
  T* ptr = data_ + index;
  T* end = data_ + size_;
  while( ptr < end )
  {
    T* next = ptr + 1;
    *ptr = *next;
    ptr = next;
  }
}

template <class T>
void VGArray<T>::insert( const T& entry, int index )
{
  assert( index >= 0 && index <= size_ );

  size( size_ + 1 );
  T* ptr = data_ + size_ - 1;
  T* end = data_ + index;

  while( ptr > end )
  {
    T* next = ptr - 1;
    *ptr = *next;
    ptr = next;
  }

  data_[index] = entry;
}

template <class T>
int VGArray<T>::find( const T& entry, int search_from ) const
{
  assert( search_from >= 0 );
  
  if( data_ )
  {
    T* ptr = data_ + search_from;
    T* end = data_ + size_;
    for(; ptr < end; ptr++ )
    {
      if( *ptr == entry )
        return ptr - data_;
    }
  }
  return -1;
}

template <class T>
void VGArray<T>::reverse()
{
  T* start = data_;
  T* end   = data_ + size_ - 1;
  while (start < end)
  {
    std::swap(*start,*end);
    start++;
    end--;
  }
}

