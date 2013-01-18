//- Class: DLIList
//-
//- Assumptions: All data are stored contiguously in the array with 
//-              empty slots at the end of the array.
//-
//- Owner: Darryl Melander

#ifndef DLILIST_HPP
#define DLILIST_HPP

#include "CubitDefines.h"
#include "CubitMessage.hpp"
#include <cstring>
#include <cassert>
#include <vector>
#include <set>
#include <stdexcept>

#define DLI_COUNT_INCREMENT 8
#define DLI_COUNT_FACTOR 1.5
template<class X> class DLIListIterator;

//! A list class, similar to a std::vector<>.
/*! DLIList is implemented as an array that is grown by a
specified amount when the list becomes full.
Most insertions and deletions at any point other than the end
of the list are handled inefficiently 
since all data from that point to the end is bubbled
down/up to fill/create the void. Operators {+} and {+=} 
are provided as an efficient means of assigning one list
to another or appending one list onto another.  The list has
a current position, which is a zero-based index into the 
list.  Many of the member functions operate on the current 
item (the item at the current position) or move the current
position.
*/
template<class X> class DLIList
{
public:
  friend class DLIListIterator<X>;

  typedef typename std::vector<X>::reference reference;
  typedef typename std::vector<X>::const_reference const_reference;
  typedef typename std::vector<X>::pointer pointer;
  typedef typename std::vector<X>::const_pointer const_pointer;
  
  //! Constructor: Create a list, allocating enough storage for \a size elements.
  /*! Although enough space is allocated for \a size elements, the list starts
  out empty (the size() function will return 0).  If the requested memory size
  is 0, then storage is not allocated until the first element is added to the list. 
  Additional storage space will be allocated if the list grows beyond its
  original size.
  \param size The amount of storage to pre-allocate for this list.
  */
  explicit DLIList (int size = 0);
  
  //! Copy constructor.
  /*! \param from The list to be copied. */
  DLIList(const DLIList<X>& from);
  
  //! Copy constructor for std::vector
  /*! \param from The list to be copied. */
  explicit DLIList(const std::vector<X>& from);

  //! Destructor: Free all resources used by this list.
  /*! The list and its storage space are freed.  Note that if this is a list
  of pointers, this destructor will \a NOT delete the objects whose pointers
  are stored in the list.
  */
  ~DLIList();
  
  //! Move the pointer to the next element in the list.
  /*! If pointer reaches the end of the list, wrap around to the beginning */
  void step();
  
  //! Move the pointer \a n positions forward in the list.
  /*! The pointer will wrap around to the beginning of the list if the end 
  is reached.  If \a n is less than zero, move backward.
  \param n The number of positions to move forward or backward.
  */
  void step(int n);
  
  //! Move the pointer to the previous element in the list. 
  /*! If it reaches the beginning of the list, wrap around to the end. */
  void back();
  
  //! Move the pointer \a n positions backward in the list.
  /*! The pointer will wrap around to the end of the list if the beginning
  is reached. If \a n is less than zero, move forward.
  */
  void back(int n);
  
  //! Set the pointer to the beginning of the list.
  void reset();
  //! Set the pointer to the end of the list.
  void last();
  
  //! Delete all elements in the list.
  /*! This function does not release memory already allocated for the list. This
  call is more efficient than creating a new list repeatedly within a loop.
  */
  void clean_out();
  
  //! Reduces the size of the list by \a k.
  /*! No list storage is freed.  Items in the array are not changed.
  If the current position is beyond the new end of the list, then the 
  current position is moved to the new end of the list.
  */
  void shrink(int k);
  
  //! Returns CUBIT_TRUE if the current position is at the beginning of the list.
  CubitBoolean is_at_beginning() const;
  
  //! Returns CUBIT_TRUE if the current position is at the end of the list.
  CubitBoolean is_at_end() const;

  //! Reverse the items in the list.
  void reverse();
  
  //! Create a copy of a list.
  /*! This is the most efficient way to do this. 
  \return A reference to this list.*/
  DLIList<X>& operator=(const DLIList<X>& from);

  //! Create a DLIList from a std::vector.
  /*! This is the most efficient way to do this. 
  \return A reference to this list.*/
  DLIList<X>& operator=(const std::vector<X>& from);

  //! Create a copy of the list an iterator was obtained from.
  /*! If \a from_iterator is not associated with a list,
  then this list becomes empty.
  \param from_iterator An iterator to another list.
  \return A reference to this list.
  */
  DLIList<X>& operator=(const DLIListIterator<X>& from_iterator);
  
  //! Append one list to another list.
  /*!  This is the most efficient way to append two lists together. 
  \param from The list to be appended to this list.
  \return A reference to this list.
  */
  DLIList<X>& operator+=(const DLIList<X>& from);

  //! Append a std::vector to a DLIList.
  /*!  This is the most efficient way to append two lists together. 
  \param from The list to be appended to this list.
  \return A reference to this list.
  */
  DLIList<X>& operator+=(const std::vector<X>& from);

  //! Subtract one list from another list.
  /*! Any element in \from which is also found in this list
  will be removed from this list.  Element ordered is retained.
  The \a from list is not modified.
  \param from The list to be subtracted from this list.
  \return A reference to this list.*/
  DLIList<X>& operator-=(const DLIList<X>& from);

  //! Compare two lists for equality.
  /*! Two lists are considered equal if they have the same contents.
  The elements do not need to be in the same order.  If values are
  repeated, they do have to appear an equal number of times in each list.
  \param from The list to compare with this list.
  \return True if the lists are equal, false otherwise.
  */
  int operator==(const DLIList<X>& from);

  //! Compare two lists for inequality.
  /*! Two lists are considered equal if they have the same contents.
  The elements do not need to be in the same order.  If values are
  repeated, they do have to appear an equal number of times in each list.
  \param from The list to compare with this list.
  \return False if the lists are equal, true otherwise.
  */
  int operator!=(const DLIList<X>& from);

  //! Gets a reference to the element with the given index.
  /*! The index is zero-based.
  \param index The index to the desired element.
  \return A reference to the indicated element.
  */
  const_reference operator[](int index) const;
  reference operator[](int index);

  //! Gets a reference to the last element in the list.
  reference last_item( void );
  const_reference last_item( void ) const;

  //! Merges the contents of another list into this list.
  /*! Each element in \a merge_list is added to this list, but only if
  it does not already appear in this list.  The result is a list where
  no value appears twice.
  
  This function runs faster if you know that no value is repeated in the
  \a merge_list.  This is indicated with the \a merge_list_unique parameter.
  If \a merge_list_unique is
  \a true, then elements of \merge_list will be checked against the original
  contents of this list, but not against the other elements of \a merge_list.
  If \a merge_list_unique is \a false, then each element will also be checked
  against the other elements of \merge_list.

  \param merge_list The list whose elements will be incorporated into this list.
  \param merge_list_unique A flag indicating whether to skip a check for 
  uniqueness between elements of \a merge_list.
  */
  void merge_unique(const DLIList<X>& merge_list, 
                    bool merge_list_unique = false);

  //! Merges the contents of a list of a different type into this list.
  /*! This function is like merge_unique(), except that the type of object
  stored by \a merge_list is not the same as this list's type.  The type of 
  object stored in the other list must be able to be static_cast<> to this 
  list's type.  
  \param merge_list The list whose elements will be incorporated into this list.
  \param merge_list_unique A flag indicating whether to skip a check for 
  uniqueness between elements of \a merge_list.
  \sa merge_unique()
  */
  template<typename Y> inline void casting_merge_unique(const DLIList<Y>& merge_list,
                                                        bool merge_list_unique = false)
    {
        // Save the current index of the merge_list
      int old_size = size();   
      int i, j, check_index;
      
      X new_item;

      // The resulting list will be at least as large as the larger of the two lists.
      // Reserve space so we don't have to reallocate so often.  Note that if
      // this list is already bigger than merge_list, the reserve won't
      // make the list shorter.
      reserve(merge_list.size());

      for ( i = 0; i < merge_list.size(); i++)
      {
          // Get the item from the merge_list and insert it into "this"
          // list if it doesn't already exist there.
        new_item = static_cast<X>(merge_list[i]);
        check_index = merge_list_unique ? old_size : size();

        // Append the new item and then remove it if necessary.
        append(new_item);
        for ( j = 0; j < check_index; j++ )
        {
          if ( listArray[j] == new_item )
          {
            listArray.resize(listArray.size()-1);
            break;
          }
        }
      }
    }

  //! Remove all elements that are not also in \a merge_list, preserving order.
  /*! This version is O(n^2).  If the order of elements in the list does 
  not matter, then consider using intersect_undordered(), which is O(n + sort_time).
  \param merge_list The list containing values to keep.
  \sa intersect_unordered(const DLIList<X>& merge_list)
  */
  void intersect(const DLIList<X>& merge_list);
  //! Remove all elements that are not also in \a merge_list, not preserving order.
  /*! This version is O(n + sort_time).  If the order of elements in the list 
  is significant, then use intersect(), which is O(n^2) but preserves list order.
  \param merge_list The list containing values to keep.
  \sa intersect(const DLIList<X>& merge_list)
  */
  void intersect_unordered(const DLIList<X>& merge_list);

  //! Appends the new item to the list, but only if it isn't already in the list.
  /*! In either case, the current position is not changed.
  \return CUBIT_TRUE if the item was added, otherwise CUBIT_FALSE.
  */
  CubitBoolean append_unique(const_reference new_item);
  
  //! Ensure each element of the list only appears once, not preserving order.
  /*! The list is first sorted for speed, so the order of elements may change.
  The current position is set to 0.
  \sa uniquify_ordered()
  */
  void uniquify_unordered();
    
  //! Ensure each element of the list only appears once, preserving order.
  /*! The order of elements is preserved, but at a speed cost compared to 
  uniquify_unordered().  The current position is set to 0.
  \sa uniquify_unordered()
  */
  void uniquify_ordered();
  
  //! Removes the current item from the list.
  /*! Remaining items with an index higher than the current item
  are moved up in the list.  The current position is set to the next item
  in the list (i.e., the index does not change) unless the removed item was
  the last item in the list, in which case the current position is set to
  the beginning of the list.
  \return The removed item.
  \sa remove(X val)
  \sa remove_all_with_value(X val)
  \sa omit(X val)
  */
  X remove ();
  
  //! Removes the next item with value \a val from the list.
  /*! Only one element is removed from the list, even if the specified
  value occurs multiple times in the list.  The list is searched from 
  the current position to the end of the list, then from the beginning of 
  the list to the current position.  The current position continues to 
  point at the same element, unless it is the element which was removed,
  in which case the behavior is identical to remove().
  \param val The value of the item to remove.
  \return Whether the item was found and removed.
  \sa remove()
  \sa remove_all_with_value(X val)
  \sa omit(X val)
  */
  bool remove (const_reference val);

  //! Remove all instances of a given value from the list.
  /*! This function can be used to remove multiple items efficiently.
  First, change the value of any elements you want to remove to \a val.
  Next, call this function with that same value.  This is more efficient
  than removing each element one at a time.  After this function call,
  the current position is set to the start of the list.
  \sa remove()
  \sa remove(X val)
  \sa omit(X val)
  */
  void remove_all_with_value(const_reference val);

  //! Removes all instances of a value from the list.
  /*! The current position of the list is not changed, unless the 
  current item is removed from the list.  If the current item is
  removed, then the current position is moved back one item, looping
  to the back of the list if necessary.
  \param val The value to remove from the list.
  \return CUBIT_TRUE if at least one instance the item was removed, CUBIT_FALSE if not.
  \sa remove()
  \sa remove(X val)
  \sa remove_all_with_value(X val)
  */
  CubitBoolean omit (const_reference val);
  
  //! Returns the value of the current item.
  /*! \return The value of the current item. If the list is empty, X(0) is returned.
  */
  const_reference get () const;
  
  //! Returns the value at next position in the list.
  /*! The current position is not changed.  If the current position is at the
  end of the list, the function wraps to the front of the list and returns the
  value of the first item.
  \return The value of the next item. If the list is empty, X(0) is returned.
  */
  X next () const;
  
  //! Returns the value at \a n positions forward in list.
  /*! The current position is not changed.  If \a n positions beyond the 
  current position would move past the end of the list, the function wraps
  around to the start of the list.
  \param n The number of positions beyond the current position to return.
  \return The value of the item \a n positions beyond the current position. 
  If the list is empty, X(0) is returned.
  */
  X next (int n) const;
  
  //! Returns the value at the previous position in list.
  /*! The current position is not changed.  If the current position is at the
  beginning of the list, the function wraps to the end of the list and returns the
  value of the last item.
  \return The value of the previous item. If the list is empty, X(0) is returned.
  */
  X prev () const;
  
  //! Returns the value at \a n positions back in list.
  /*! The current position is not changed.  If \a n positions before the 
  current position would move before the beginning of the list, the function wraps
  around to the end of the list.
  \param n The number of positions behind the current position to return.
  \return The value of the item \a n positions before the current position. 
  If the list is empty, X(0) is returned.
  */
  X prev (int n) const;
  
  //! Returns the current value, then advances the current position by one.
  /*! If the current position is at the end of the list, the function will
  wrap to the beginning of the list.
  \return The value at the current position, before stepping.
  */
  const_reference get_and_step ();
  
  //! Returns the current value, then moves the current position back one.
  /*! If the current position is at the beginning of the list, the function will
  wrap to the end of the list.
  \return The value at the current position, before stepping back.
  If the list is empty, X(0) is returned.
  */
  const_reference get_and_back ();
  
  //! Advances the current position by one, then returns the current value.
  /*! If the current position is at the end of the list, the function will
  wrap to the beginning of the list.
  \return The value at the current position, after stepping.
  */
  const_reference step_and_get ();
  
  //! Moves to the next instance of this value, wrapping if necessary.
  /*! \return Returns \a true if the item was found in the list, otherwise returns \a false.
  */
  CubitBoolean move_to(const_reference item);

  //! Return \a true if the item is in the list, \a false otherwise.
  /*! \return Returns \a true if the item was found in the list, otherwise returns \a false.
  */
  CubitBoolean is_in_list(const_reference item) const;

  //! Returns and removes last value in list.
  /*! \return Returns the last value in the list.  If the list is empty, returns X(0).
  */
  X pop();

  //! Puts a value on the end of the list.
  /*! The current position is then set to the end of the list.
  \param val The value to place at the end of the list.
  \return The value placed on the end of the list.
  \sa push_back(X val)
  */
  X push(X val);

  //! Insert an item into the list, after current item.
  /*! The newly inserted item becomes the current item.  This function is
  inefficient due to bubbling.
  \param val The item to insert.
  */
  void insert (X val);

  //! Add a value at start of the list. 
  /*! The current position is set to the beginning of the list.
  Inefficient due to bubbling.
  \param val The value to place at the start of the list.
  */
  void insert_first(X val);

  //! Place an item at the end of the list.
  /*! This is the most efficient way to insert an item to the list.
  The current position is unchanged.
  \param val The value to place at the end of the list.  
  */
  void append(const_reference new_item);

  //! Remove the current value and put last value in the list in its place.
  /*!  If list order isn't important, this is much more efficient than remove().
  \return The value removed from the list, or X(0) if the list was empty.
  */
  X extract();
  
  //! Change the value of the current item.
  /*! Because this function does not actually remove or insert an element,
  it is quite efficient.  If the list is empty, the new value will not be
  inserted into the list.
  \return The former value of the current item, or X(0) if the list was empty.
  */
  X change_to(X val);
  
  //! Orders the list elements from lowest to highest.
  /*! The sort order is determined by operator>= for the
  stored element type.  
  
  Use reverse after this call if you want to go the
  other direction.
  */
  void sort();

  //! A function which determines the relative order of objects.
  /*!
  The SortFunction should return a negative number if \a a should
  be before \a b, a positive number if \a b comes before \a a, and
  zero if they are equal, or relative order doesn't matter.
  \param a The first object in the comparison
  \param b The second object in the comparison
  \sa sort(SortFunction f)
  */
  typedef int (*SortFunction)(X& a, X& b);

  //! Orders the list elements from lowest to highest, as defined by a function.
  /*! The sort order is determined by the passed in SortFunction.
  
  Use reverse after this call if you want to go the
  other direction.
  \param f The function which determines the sort order.
  \sa SortFunction
  */
  void sort(SortFunction f);

  //! Allocate enough space for at least \a min_size elements.
  /*! If there is already enough space allocated, the function does nothing; this
  function will never reduce the amount of memory allocated to the list.
  \param min_size The minimum number of elements to be prepared to store.
  */
  void reserve(int min_size);
  
  //! Returns the number of items in the list
  /*! \return The number of items current in the list. */
  int size() const
    { return listArray.size(); }

  //! Returns current index of the current position
  /*! \return the current position */
  int get_index()
     { return index; }
  
  //! Return the index of an item in list, or -1 if the item is not in the list.
  /*! The location of the first instance of \a val is returned as a zero-based
  index from the start of the list.  The list is searched starting from the 
  current position, wrapping to the front of the list if necessary.
  \param val The value to search for.
  \return The index of the first instance of \a val, or -1 if the value is not found.
  */
  int where_is_item(const_reference val) const;

  //! Returns the number of bytes allocated for this list's storage space.
  int memory_use(CubitBoolean verbose_boolean);
   
  //! Copy this list's contents into an array.
  /*! It is assumed that \a other_array is big enough to hold all of this
  list's elements.  No check is made to verify this.
  \param other_array The array into which this list's contents will be copied.
  */
  void copy_to(X *other_array);

  //! Copy items from an array into this list.
  /*! Any prior contents of this list are removed.  The list allocates additional storage if necessary to hold all of
  \a other_array's contents.
  \param other_array The array from which items will be copied.
  \param other_size The number of items to be copied from \a other_array.
  */
  void copy_from(X *other_array, const int other_size);

  //! Moves to the nearest instance of \a val, searching both forward and backward.
  /*! 
  \return True if an item with the specified value was found, false otherwise.
  */
  CubitBoolean move_to_nearby(const X val);

  //! Returns the distance to the nearest element with a given value.
  /*! The returned value is always positive, regardless of whether the
  closest element is before or after the current position.
  \param val The value to search for.
  \return The distance to the closest element with value \a val.
  */
  int distance_to_nearby(const X val);

  //! Set the current position to point at one of the two items, where the previous item is the other item.
  /*! This function looks for a place in the list where these two items are adjacent 
  to each other, in either order.  The current position is then set to the later
  of the two items.  The function acts in a wrapped manner, such that the last
  item in the list is considered adjacent to the first item in the list.
  \param val1 One of the values to look for.
  \param val2 The other value to look for.
  \return \a true if the two items are found adjacent to each other, \a false otherwise.
  */
  CubitBoolean move_between(X val1, X val2);

  //! Return a std::vector with the same contents as this list.
  /*! This function creates a std::vector from the DLIList and returns it.
  */
  const std::vector<X>& as_vector();

private:
  void lengthen_list(int by_how_much = DLI_COUNT_INCREMENT,
                     double by_what_factor = DLI_COUNT_FACTOR );
    //- Makes the array bigger. Multiply current size 
    //- "by_what_factor" and add 'by_how_much'.
  
  int  index;      // Index of current item in {listArray}
  std::vector<X> listArray;
};

// *** inline function definitions *******************************************
template <class X> inline 
DLIList<X>& DLIList<X>::operator=(const DLIListIterator<X>& from_iter)
{ 
  if ( from_iter.dlIList ) 
    *this = *(from_iter.dlIList); 
  else
    clean_out();
  return *this; 
}

template <class X> inline
typename DLIList<X>::const_reference DLIList<X>::operator[](int index) const
{
  if(index < 0 || index >= (int)listArray.size())
    throw std::out_of_range("Index out of Bounds\n");
  return listArray[index];
}

template <class X> inline
typename DLIList<X>::reference DLIList<X>::operator[](int index)
{
  if(index < 0 || index >= (int)listArray.size())
    throw std::out_of_range("Index out of Bounds\n");
  return listArray[index];
}

template <class X> inline
typename DLIList<X>::reference DLIList<X>::last_item(void)
{
    assert( listArray.size() > 0 );
    return listArray[listArray.size()-1];
}

template <class X> inline
typename DLIList<X>::const_reference DLIList<X>::last_item(void) const
{
    assert( listArray.size() > 0 );
    return listArray[listArray.size()-1];
}

template <class X> inline void DLIList<X>::step()
{
  if (!listArray.empty())
    index++; 
  if (index >= (int)listArray.size())
    index = 0;
}

template <class X> inline void DLIList<X>::step(int n)
{
  const int itemCount = listArray.size();
  if (itemCount)
  {
    index += n;
    if (index < 0)
      index = itemCount - (-index)%itemCount;
    if (index >= itemCount)
      index %= itemCount;
  }
}

// Chop off the last (k) items in the list.
template <class X> inline void DLIList<X>::shrink(int k)
{
  const size_t itemCount = listArray.size();
  if (k > 0)
  {
    if (itemCount > k)
      listArray.resize(itemCount - k);
    else
      listArray.clear();
  }
  if (index >= itemCount)
  {
    if (itemCount)
      index = itemCount-1;
    else
      index = 0;
  }
}

template <class X> inline void DLIList<X>::back()
{
  const size_t itemCount = listArray.size();
  if (itemCount)
  {
    index--;
    if (index < 0)
      index = itemCount - 1;
  }
}

template <class X> inline void DLIList<X>::back(int n)
{
  step(-n);
}

template <class X> inline void DLIList<X>::reset()
{
  index = 0;
}

template <class X> inline void DLIList<X>::last()
{
  const size_t itemCount = listArray.size();
  if (itemCount)
    index = itemCount - 1;
}

template <class X> inline CubitBoolean DLIList<X>::is_at_beginning() const
{
  if (!index)
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}

template <class X> inline CubitBoolean DLIList<X>::is_at_end() const
{
  const size_t itemCount = listArray.size();
  if (itemCount == 0 || index == (int)itemCount-1)
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}

template <class X> inline void DLIList<X>::clean_out()
{
  listArray.clear();
  index = 0;
}

template <class X> inline CubitBoolean DLIList<X>::is_in_list(const_reference item) const
{
  return where_is_item(item) >= 0 ? CUBIT_TRUE : CUBIT_FALSE;
}

//- Add item to end of list, do not change current index value.
//- This function used to reduce time required to do an insert 
//- followed by a move_to back to where we were.
template <class X> inline void DLIList<X>::append(const_reference new_item)
{
    // see if the list must be lengthened
  if (listArray.capacity() == listArray.size())
    lengthen_list();

  listArray.push_back(new_item);
}

template <class X> inline typename DLIList<X>::const_reference DLIList<X>::get() const
{
  if ( listArray.empty() )
  {
    throw std::out_of_range("Attempted get of empty DLIList\n");
    /*PRINT_WARNING("Attempted get of empty DLIList\n");*/
  }

  return listArray[index];
}

template <class X> inline typename DLIList<X>::const_reference DLIList<X>::get_and_step()
{
  if ( listArray.empty() )
  {
    throw std::out_of_range("Attempted get_and_step from empty DLIList\n");
    /*PRINT_WARNING("Attempted get_and_step from empty DLIList\n");*/
  }
  const_reference temp = listArray[index++];
  if (index == (int)listArray.size())
    index=0;
  return temp;
}

template <class X> inline typename DLIList<X>::const_reference DLIList<X>::get_and_back ()
{
   if ( listArray.empty() )
   {
     throw std::out_of_range("Attempted get_and_back from empty DLIList\n");
      /*PRINT_WARNING("Attempted get_and_back from empty DLIList\n");*/
   }
   const_reference temp = listArray[index--];
   if (index < 0)
      index=listArray.size()-1;
   return temp;
}

template <class X> inline X DLIList<X>::next() const
{
  if (listArray.empty())
  {
    throw std::out_of_range("Attempted next of empty DLIList\n");
    /*PRINT_WARNING("Attempted next of empty DLIList\n");*/
  }
  else if (index == listArray.size()-1)
    return listArray[0];
  else
    return listArray[index+1];
}

template <class X> inline X DLIList<X>::next(int n) const
{
  if ( listArray.empty() )
  {
    throw std::out_of_range("Attempted next of empty DLIList\n");
    /*PRINT_WARNING("Attempted next of empty DLIList\n");*/
  }
  else
  {
      // return the proper index
      // beware of negative n leading to negative new_index
    int new_index = index+n;
    while (new_index < 0)
      new_index += listArray.size();
    if (new_index >= (int)listArray.size())
      new_index %= listArray.size();
    
    return listArray[new_index];
  }
}

template <class X> inline X DLIList<X>::prev() const
{
  return this->next(-1);
}

template <class X> inline X DLIList<X>::prev(int n) const
{
  return this->next(-n);
}

//- put new item in list at index 0 and make it current
template <class X> inline void DLIList<X>::insert_first(X new_item)
{
  // set index to -1 ... index will be set to 0 in insert
  index = -1;
  insert(new_item);
}

template <class X> inline CubitBoolean DLIList<X>::move_to (const_reference item)
{
  int item_index = where_is_item(item);
  CubitBoolean rv = CUBIT_FALSE;
  if (item_index >= 0)
  {
    index = item_index;
    rv = CUBIT_TRUE;
  }
  return rv;
}

template <class X> inline X DLIList<X>::change_to (X new_value)
{
  if ( listArray.empty() )
  {
    throw std::out_of_range("DLIList: Attempted update of empty list\n");
    //PRINT_WARNING("DLIList: Attempted update of empty list\n");
  }

  X temp = listArray[index];
  listArray[index] = new_value;
  return temp;
}

// Removes every instance of 'val' in list.
// Set to beginning of list after this call.
template <class X> inline void DLIList<X>::remove_all_with_value(const_reference val)
{
  int j = 0;
  int i = 0;
  const size_t itemCount = listArray.size();

  for ( ; i < (int)itemCount; i++)
    if (listArray[i] != val && j++ != i)
      listArray[j-1] = listArray[i];

  listArray.resize(j);
  index = 0;
}

// Searches for item and removes the next instance from the list.
// If the item was found and removed it is returned, else 0 is returned.
template <class X> inline bool DLIList<X>::remove (const_reference item)
{
  if (move_to(item))
  {
    remove();
    return true;
  }
  return false;
}

template <class X> inline int DLIList<X>::where_is_item (const_reference item) const
{
  if (listArray.empty())
    return -1;

  const int itemCount = listArray.size();
  
    // loop through list searching for item ...
    // if found return index
  
    // Search from current index to end of array
  int i;
  for (i=index; i < itemCount; i++)
    if (listArray[i] == item)
      return i;
  
    // Now search from beginning of array to index...
  for (i = 0; i < index && i < itemCount; i++)
    if (listArray[i] == item)
      return i;
  
    // item is not in array, return -1
  return -1;
}
//- Append this new_item to the list if it doesn't already exist in it.
//- Make sure the index is unchanged.
template <class X> inline CubitBoolean DLIList<X>::append_unique(const_reference new_item)
{
    // Append new_item, if it isn't already there.
  if( where_is_item(new_item) < 0 ) 
  {
    append (new_item);
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

template <class X> inline X DLIList<X>::push(X value)
{
  //- swapped last and append; current position is set new end
   append(value);
   last();
   return value;
}

template <class X> inline X DLIList<X>::pop()
{
  last();
  return remove();
}

//- Constructor: Create a list with initial size 0. The list
//- will be grown by {increment} each time it is filled. Memory for the
//- list is not allocated until the first element is inserted using
//- {insertLink}.
template <class X> inline DLIList<X>::DLIList (int size)
{
   index      = 0;
   listArray.reserve(size);
}

//- Copy Constructor
template <class X> inline DLIList<X>::DLIList(const DLIList<X>& from)
{
   if (&from != this)
   {
     index = from.index;
     listArray = from.listArray;
   }
}

//- Copy constructor for std::vector
template <class X> inline DLIList<X>::DLIList(const std::vector<X>& from)
{
    // Setup the variables
    index = 0;
    listArray = from;
}

// Destructor
template <class X> inline DLIList<X>::~DLIList()
{
}

template <class X> inline void DLIList<X>::reserve(int min_size)
{
  listArray.reserve(min_size);
}

template <class X> inline void DLIList<X>::lengthen_list(int by_how_much,
                                                  double by_what_factor)
{
    // Make a new array
  int new_size = (int) ((double)listArray.capacity() * by_what_factor) + by_how_much;
  reserve(new_size);
}

//- put new item in list after current item and make it current
template <class X> inline void DLIList<X>::insert(X new_item)
{
     // see if the list must be lengthened
   if ( listArray.size() == listArray.capacity())
   {
      lengthen_list();
   }

     // set new index
   if ( !listArray.empty() )
   {
      index++;
   }
   else
   {
      index = 0;
   }

   listArray.insert(listArray.begin()+index, new_item);
}

//- merge the input list, merge_list, with the current list, making
//- sure not to add duplicate items into the current list
template <class X> inline
void DLIList<X>::merge_unique ( const DLIList<X>& merge_list,
                                bool  merge_list_unique )
{
  this->casting_merge_unique(merge_list, merge_list_unique);
}

template <class X> inline void DLIList<X>::intersect_unordered(
  const DLIList<X>& merge_list )
{
  if ( &merge_list == this )
     return;

  DLIList <X> intersect_list(merge_list);   // copy input list so can sort
  sort();
  intersect_list.sort();

  typename std::vector<X>::iterator iter1 = listArray.begin();                      // iterator for this array
  typename std::vector<X>::iterator end1 = listArray.end();               // end of this array
  typename std::vector<X>::iterator iter2 = intersect_list.listArray.begin();       // iterstor for other array
  typename std::vector<X>::iterator end2 = intersect_list.listArray.end();// end of other array
  typename std::vector<X>::iterator insert;                     // location of last insert
  bool first = true;

  for ( ; iter1 < end1; ++iter1 )
  {
    while (iter2 < end2 && *iter2 < *iter1)
      ++iter2;

    if (iter2 == end2)
      break;

    // items are the same and ...
    if (*iter2 == *iter1)
    {
      // is the first item or ...
      if(first)
      {
        first = false;
        insert = listArray.begin();
        *insert = *iter1;
      }
      // is not the same as the previous item
      else if(*iter1 != *insert)
      {
        *++insert = *iter1;
      }
    }
  }

  if(first)
  {
    // no intersections
    listArray.clear();
  }
  else
  {
    listArray.resize(insert - listArray.begin() + 1);
  }
  reset();
}

template <class X> inline void DLIList<X>::intersect ( const DLIList<X>& merge_list )
{
  if ( &merge_list == this )
     return;

  const int itemCount = listArray.size();
  std::vector<X> tmp;

  for ( int i=0; i<itemCount; i++ )
  {
    if (merge_list.is_in_list(listArray[i]))
    {
      tmp.push_back(listArray[i]);
    }
  }

  this->listArray.swap(tmp);
  index = 0;
}

//template <class X> inline void DLIList<X>::intersect ( void* merge_list )
//{
//  intersect( *(DLIList<X>*)merge_list );
//}


//- remove the item at the current location and return a pointer to it.
//- The next node becomes the current node
//- Returns 0 if there are no items in list
template <class X> inline X DLIList<X>::remove ()
{
   if ( listArray.empty() )
   {
     throw std::out_of_range("Attempted link removal from empty DLIList\n");
      /*PRINT_WARNING("Attempted link removal from empty DLIList\n");*/
   }

     // save the current value
   X temp = listArray[index];

   listArray.erase(listArray.begin()+index);

   if ( index == (int)listArray.size() )
   {
      index = 0;
   }

   return temp;
}


//- remove the item at the current location and return a pointer to it.
//- used for efficiency in cases where preservation of list order is not
//- important.  moves last list item (itemCount - 1) to current index and
//- decrements itemCount.  eliminates the need to perform the list bubble
//- down (i.e. cut_link) but sacrifices list order in the process.  this
//- function should not be used when up-stream order from the removed node is
//- important.  when processing a list using this function the user should
//- reset the list to the head (index = 0) before beginning to ensure all
//- list nodes are processed properly.
//- Returns 0 if there are no items in list
template <class X> inline X DLIList<X>::extract ()
{
   if ( listArray.empty() )
   {
     throw std::out_of_range("Attempted link removal from empty DLIList\n");
      /*PRINT_WARNING("Attempted link removal from empty DLIList\n");*/
   }

     // save the current value
   X temp = listArray[index];

     // assign last node to the current index
   listArray[index] = listArray[listArray.size() - 1];

     // decrement
   listArray.resize(listArray.size()-1);
   if ( index == listArray.size() && index )
        // The choices here are index at beginning or end.
        // End seems to work better when 'extract' is being used
        // with calls to 'move_to_item'.
      index--;

   return temp;
}

//+//Added so list removals don't disturb current position. PRK 05-23-94
//+//Corrected for omitting the last item in the list. PRK 09-16-94
//+//Corrected for omitting before the current position. PRK 10-07-94
//- Finds instance of item by matching value and delets first instance
//- of it from the list. The current position of the list is not changed.
//- Returns CUBIT_TRUE if the item was found and deleted, CUBIT_FALSE if not.
template <class X> inline CubitBoolean DLIList<X>::omit(const_reference old_val)
{
   int scan_index;
   int squeeze_index = 0;
   CubitBoolean found = CUBIT_FALSE;
   for(scan_index = 0; scan_index < listArray.size(); ++scan_index)
   {
      if(listArray[scan_index] == old_val)
      {
         found = CUBIT_TRUE;
         if(index == scan_index)
            index = squeeze_index - 1;
      }
      else
      {
         if(scan_index != squeeze_index)
         {
            listArray[squeeze_index] = listArray[scan_index];
            if(index == scan_index) index = squeeze_index;
         }
         ++squeeze_index;
      }
   }

   if(found)
   {
      listArray.resize(squeeze_index);
//+//   If the first item was deleted and index pointed to it, make an
//+//   adjustment here. If itemCount is zero, don't assign -1 again.
      if(index < 0)
         index = listArray.empty() ? 0 : listArray.size()-1;
   }

   return found;
}

template <class X> inline typename DLIList<X>::const_reference DLIList<X>::step_and_get ()
{
   if ( listArray.empty() )
   {
     throw std::out_of_range("Attempted step_and_get from empty DLIList\n");
      /*PRINT_WARNING("Attempted step_and_get from empty DLIList\n");*/
   }

   if (++index == (int) listArray.size())
      index=0;
   return listArray[index];
}

template <class X> inline DLIList<X>& DLIList<X>::
                       operator=(const DLIList<X>& from)
{
   if (this != &from)
   {
      index = from.index;
      listArray = from.listArray;
   }
   return *this;
}

template <class X> inline DLIList<X>& DLIList<X>::
                       operator=(const std::vector<X>& from)
{
  index = 0;
  listArray = from;
  return *this;
}

template <class X> inline DLIList<X>& DLIList<X>::
                       operator+=(const DLIList<X>& from)
{
   listArray.insert(listArray.end(), from.listArray.begin(), from.listArray.end());
   return *this;
}

template <class X> inline DLIList<X>& DLIList<X>::
                       operator+=(const std::vector<X>& from)
{
  listArray.insert(listArray.end(), from.listArray.begin(), from.listArray.end());
  return *this;
}

template <class X> inline DLIList<X>& DLIList<X>::
                       operator-=(const DLIList<X>& from)
{
    // step through items in from list, removing them from this list.
   for (int i = from.listArray.size(); i--; )
   {
     // quit early if this list is empty.
     if (listArray.empty())
       break;
     remove_all_with_value(from.listArray[i]);
   }
   return *this;
}

template <class X> inline int DLIList<X>::operator==(const DLIList<X> &from)
{
  if(listArray.size() != from.listArray.size())
     return CUBIT_FALSE;
  DLIList<X> temp_list = from;
  for( int i = 0; i < listArray.size(); i++)
     if(temp_list.move_to(listArray[i]))
        temp_list.remove();
  return temp_list.listArray.size() == 0;
}
template <class X> inline int DLIList<X>::operator!=(const DLIList<X> &from)
{
  return !( *this == from );
}

// Sorts list from low to high value, according to the > operator
// of the list type.
// List is sorted using a standard Heap Sort algorithm ("Algorithms in C++",
// Robert Sedgewick).
template <class X> inline void DLIList<X>::sort(typename DLIList<X>::SortFunction f)
{
    // Only sort it if there is more than one
    // item in the list
  const size_t itemCount = listArray.size();
  if (itemCount > 1)
  {
    int mid = (itemCount >> 1) + 1;
    int ir = itemCount;
    X temp_element;

      // You loop until ir is 1 (ir = iterations remaining)
    while(CUBIT_TRUE)
    {
      if (mid > 1)
      {
        mid--;
        temp_element = listArray[mid - 1];
      }
      else
      {
        ir--;
        temp_element = listArray[ir];
        listArray[ir] = listArray[0];
        if (ir == 1)
        {
          listArray[0] = temp_element;
          return;
        }
      }

      int i = mid;
      int j = mid + mid;

      while (j <= ir)
      {
        if (j < ir)
        {
          if (f(listArray[j], listArray[j - 1]) >= 0)
            j++;
        }
        if (f(listArray[j - 1], temp_element) >= 0)
        {
          listArray[i - 1] = listArray[j - 1];
          i = j;
          j += j;
        }
        else
        {
          j = ir + 1;
        }
      }
      listArray[i - 1] = temp_element;
    }
  }
}

template <class X> inline void DLIList<X>::sort()
{
  const size_t itemCount = listArray.size();
  // Only sort it if there is more than one
  // item in the list
  if (itemCount > 1)
  {
    int mid = (itemCount >> 1) + 1;
    int ir = itemCount;
    X temp_element;

      // You loop until ir is 1 (ir = iterations remaining)
    while(CUBIT_TRUE)
    {
      if (mid > 1)
      {
        mid--;
        temp_element = listArray[mid - 1];
      }
      else
      {
        ir--;
        temp_element = listArray[ir];
        listArray[ir] = listArray[0];
        if (ir == 1)
        {
          listArray[0] = temp_element;
          return;
        }
      }

      int i = mid;
      int j = mid + mid;

      while (j <= ir)
      {
        if (j < ir)
        {
          if (listArray[j] >= listArray[j - 1])
            j++;
        }
        if (listArray[j - 1] >= temp_element)
        {
          listArray[i - 1] = listArray[j - 1];
          i = j;
          j += j;
        }
        else
        {
          j = ir + 1;
        }
      }
      listArray[i - 1] = temp_element;
    }
  }

}

template <class X> inline void DLIList<X>::reverse()
{
  int front = 0;
  int tail  = listArray.size()-1;
  X temp;

  while (front < tail)
  {
     temp             = listArray[front];
     listArray[front] = listArray[tail];
     listArray[tail]  = temp;
     tail--;
     front++;
  }
}

template <class X> inline int DLIList<X>::memory_use(CubitBoolean verbose_boolean)
{
   // report amount of memory allocated

   int size = listArray.capacity() * sizeof(X);

   if (verbose_boolean)
   {
      PRINT_INFO("      DLIList: %d bytes\n",size);
   }

   return size;
}

template <class X> inline void DLIList<X>::copy_to(X *other_array)
    //- copy this list's listArray into other_array
{
  if(other_array == 0)
    throw std::invalid_argument("Array is NULL");

  if (listArray.size())
  {
    const size_t itemCount = listArray.size();
    for(size_t i=0; i<itemCount; i++)
    {
      other_array[i] = listArray[i];
    }
  }
}

template <class X> inline void DLIList<X>::copy_from(X *other_array, const int other_size)
  //- copy other_array into listArray
{
  if(other_array == 0)
    throw std::invalid_argument("Array is NULL");

  listArray.clear();
  listArray.insert(listArray.end(), other_array, other_array+other_size);
  index = 0;
}

template <class X> inline CubitBoolean DLIList<X>::move_to_nearby(const X item)
{
//  if (nullItem && (X)nullItem == item)
//     return CUBIT_FALSE;
//  else
  if (listArray.empty())
     return CUBIT_FALSE;
  else

  {
    const size_t itemCount = listArray.size();
    typename std::vector<X>::iterator ptr_up = listArray.begin() + index;
    if (*ptr_up == item)
       return CUBIT_TRUE;

    int i_up, i_down;
    i_up = i_down = index;
    typename std::vector<X>::iterator ptr_down = ptr_up;
    while (1)
    {
        // check forward in the list increment
      if ( ++i_up < itemCount )
         ptr_up++;
      else
      {
        i_up = 0;
        ptr_up = listArray.begin();
      }
        // check
      if ( *ptr_up == item )
      {
        index = i_up;
        return CUBIT_TRUE;
      }
      if ( i_up == i_down )
      {
        return CUBIT_FALSE;
      }

        // check backward in the list
        //  increment
      if ( --i_down >= 0 )
         ptr_down--;
      else
      {
        i_down = itemCount-1;
        ptr_down = listArray.begin()+ i_down;
      }
        // check
      if ( *ptr_down == item )
      {
        index = i_down;
        return CUBIT_TRUE;
      }
      if ( i_up == i_down )
      {
        return CUBIT_FALSE;
      }

    }
  }
}

template <class X> inline int DLIList<X>::distance_to_nearby(const X body)
{
//  if (nullItem && (X)nullItem == body)
//     return CUBIT_FALSE;
//  else
  {
    int old_index = index;
    move_to_nearby( body );
    int distance = abs(index - old_index);
    if (distance > listArray.size() / 2)
       distance = listArray.size() - distance;
    index= old_index;
    return distance;
  }
}

template <class X> inline CubitBoolean DLIList<X>::move_between(X item1, X item2)
{
  {
    if ( listArray.empty() )
       return CUBIT_FALSE;

    const size_t itemCount = listArray.size();

    for ( int i = 0; i < itemCount; i++ )
    {
      if ( listArray[i] == item1 )
      {
          //  dance around looking for item2
        if ( i+1 < itemCount ) {
          if ( listArray[i+1] == item2 ) {
            index = i+1;
            return CUBIT_TRUE;
          }
        }
        if ( i > 0 ) {
          if ( listArray[i-1] == item2 ) {
            index = i;
            return CUBIT_TRUE;
          }
        }
        if ( ( i+1 == itemCount && listArray[0] == item2 )
             || ( i == 0 && listArray[itemCount-1] == item2 ) )
        {
          index = 0;
          return CUBIT_TRUE;
        }
      }
    }
    return CUBIT_FALSE;
  }
}

// Removes every instance of 'val' in list.
// Set to beginning of list after this call.
template <class X> inline void DLIList<X>::uniquify_unordered()
{
  const size_t itemCount = listArray.size();
  if ( itemCount < 2 )
    return;

  sort();

  int j = 1;
  int i = 0;

  while (j < itemCount)
  {
    if (listArray[i] != listArray[j])
      listArray[++i] = listArray[j++];
    else
      j++;
  }

  listArray.resize(i + 1);
  index = 0;
}

template <class X> inline void DLIList<X>::uniquify_ordered()
{
  std::set<X> encountered;
  typename std::vector<X>::iterator read, write, end = listArray.end();

    // To avoid copying the array onto itself in the case
    // where the list contains no duplicates, loop twice.

    // Find first duplicate entry.
  for ( read = listArray.begin(); read != end; ++read )
    if ( ! encountered.insert(*read).second )
      break;

    // Now compact array, removing duplicates.  If there
    // are no duplicates, this loop will not run (read == end).
  for ( write = read; read != end; ++read )
    if ( encountered.insert(*read).second )
      *(write++) = *read;

  listArray.resize(write - listArray.begin());
  index = 0;
}

template <class X> inline const std::vector<X>& DLIList<X>::as_vector()
{
  return listArray;
}

template <class X> class DLIListIterator
{
  friend class DLIList<X>;
public:
    // constructor
  DLIListIterator() :
      mIndex(0),
      dlIList(NULL)
    {}
  
    // destructor
  virtual ~DLIListIterator()
    {}
  
  void watch( const DLIList <X> *dl_i_list )
    {
      dlIList = dl_i_list;
      mIndex = dlIList->index;
    }
  
  void reset()
    { mIndex = 0; }
  
  void step( int n = 1 ) 
    {
      if (!size())
        mIndex = 0;
      else
      {
        mIndex += n;
        while (mIndex < 0)
          mIndex += dlIList->size();
        mIndex %= dlIList->size();
      }
    }
  
  int size() const
    { return dlIList ? dlIList->size() : 0; }
  
  X get() const
    { return next(0); }
  X next( int n = 1 ) const;
  X get_and_step( int n = 1 ) 
    {
      X temp = get();
      step( n );
      return temp;
    }
  
    //- Moves to the next instance of this value, wrapping if necessary
  CubitBoolean move_to(X item)
    { 
      mIndex = where_is_item( item );
      return mIndex >= 0 ? CUBIT_TRUE : CUBIT_FALSE; 
    }
  
    //- Return True if item is in the list, false otherwise.
  CubitBoolean is_in_list(X item) const
    { return where_is_item(item) >= 0 ? CUBIT_TRUE : CUBIT_FALSE; }
  
    // Returns CUBIT_TRUE if we're pointing at the last item in the list.
  CubitBoolean is_at_end() const
    {
      if (size() == 0 || mIndex == size() - 1)
        return CUBIT_TRUE;
      return CUBIT_FALSE;
    }
  
    //- returns CUBIT_TRUE if we're pointing at the firt item in the list.
  CubitBoolean is_at_beginning() const
    {
      if (!mIndex || mIndex == size())
        return CUBIT_TRUE;
      return CUBIT_FALSE;
    }
  
protected:
    // utility functions
  int where_is_item( X item ) const
    { return dlIList ? dlIList->where_is_item( item ) : -1; }
  
    // data
  int mIndex;
  const DLIList<X> *dlIList;
};

template <class X> 
inline X DLIListIterator<X>::next( int n ) const
{
  int size_loc = size();
  if ( size_loc == 0 )
    return NULL;
  int index_loc = mIndex + n;
  while ( index_loc < 0 )
    index_loc += size_loc;
  while ( index_loc >= size_loc )
    index_loc -= size_loc;
  
    // if dlIList == NULL, then size == 0 and returned already 
  return dlIList->listArray[index_loc];
}

#endif //- DLILIST_HPP

