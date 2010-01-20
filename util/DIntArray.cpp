//- Class: DIntArray
//- Description:
//- Owner: Kirk Walton
//- Checked by:

#include "DIntArray.hpp"
#include <cassert>

void DIntArray::append(int item) 
{
  intArray.push_back(item);
}

int DIntArray::size() 
{
  return (intArray.size());
}

int DIntArray::operator[](int index)
{
  return (intArray[index]);
}
