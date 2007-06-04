//-------------------------------------------------------------------------
// Filename      : CompositeAttrib.cpp
//
// Purpose       : Container for attribute data placed on composite geometry.
//
// Special Notes : This object is intended for internal use by CompositeGeom
//                 exclusively.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/01/03
//-------------------------------------------------------------------------

#include "CompositeAttrib.hpp"
#include "CubitString.hpp"
#include "CubitSimpleAttrib.hpp"

//-------------------------------------------------------------------------
// Purpose       : Construct composite attribute object from CSA
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/30/03
//-------------------------------------------------------------------------
CompositeAttrib::CompositeAttrib( CubitSimpleAttrib* attrib, CompositeAttrib* n ) 
  : int_array(0), real_array(0), string_array(0), next(n)
{
  int_count = attrib->int_data_list()->size();
  if (int_count) 
  {
    int_array = new int[int_count];
    attrib->int_data_list()->reset();
    for (int i = 0; i < int_count; i++)
      int_array[i] = *attrib->int_data_list()->get_and_step();
  }

  real_count = attrib->double_data_list()->size();
  if (real_count) 
  {
    real_array = new double[real_count];
    attrib->double_data_list()->reset();
    for (int i = 0; i < real_count; i++)
      real_array[i] = *attrib->double_data_list()->get_and_step();
  }

  string_count = attrib->string_data_list()->size();
  if (string_count) 
  {
    string_array = new CubitString[string_count];
    attrib->string_data_list()->reset();
    for (int i = 0; i < string_count; i++)
      string_array[i] = *attrib->string_data_list()->get_and_step();
  }
}  
  
//-------------------------------------------------------------------------
// Purpose       : Copy constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/30/03
//-------------------------------------------------------------------------
CompositeAttrib::CompositeAttrib( const CompositeAttrib& copy) 
  : int_count(copy.int_count), int_array(0), 
    real_count(copy.real_count), real_array(0),
    string_count(copy.string_count), string_array(0), 
    next(0)
{
  if (int_count)
  {
    int_array = new int[int_count];
    memcpy(int_array, copy.int_array, int_count * sizeof(int));
  }
  if (real_count)
  {
    real_array = new double[real_count];
    memcpy(real_array, copy.real_array, real_count * sizeof(double));
  }
  if (string_count)
  {
    string_array = new CubitString[string_count];
    CubitString *read = copy.string_array;
    CubitString *write = string_array;
    CubitString* end = read + string_count;
    while (read < end)
      *(write++) = *(read++);
  }
}
    
//-------------------------------------------------------------------------
// Purpose       : Populate attribute for saving
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/30/03
//-------------------------------------------------------------------------
void CompositeAttrib::append_to_csa( CubitSimpleAttrib* attrib ) const
{
  append_to_lists(*attrib->string_data_list(),
                  *attrib->int_data_list(),
                  *attrib->double_data_list());
}

//-------------------------------------------------------------------------
// Purpose       : Common functionality for append_to_csa() and csa()
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/30/03
//-------------------------------------------------------------------------
void CompositeAttrib::append_to_lists( DLIList<CubitString*>& string_list,
                                       DLIList<int*>& int_list,
                                       DLIList<double*>& real_list ) const
{
  int* iitor = int_array;
  int *const iend = int_array + int_count;
  while( iitor < iend )
    int_list.append(iitor++);
  
  double* ritor = real_array;
  double *const rend = real_array + real_count;
  while( ritor < rend )
    real_list.append(ritor++);
  
  CubitString* sitor = string_array;
  CubitString *const send = string_array + string_count;
  while( sitor < send )
    string_list.append(sitor++);
}
  
//-------------------------------------------------------------------------
// Purpose       : Copy data into CSA for return to higher level
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/30/03
//-------------------------------------------------------------------------
CubitSimpleAttrib* CompositeAttrib::csa() const
{
  DLIList<CubitString*> string_list(string_count);
  DLIList<int*> int_list(int_count);
  DLIList<double*> real_list(real_count);
  append_to_lists( string_list, int_list, real_list );
  CubitSimpleAttrib* result = new CubitSimpleAttrib;
  result->initialize_from_lists_of_ptrs( &string_list,
                                real_count ? &real_list : 0, 
                                int_count ? &int_list : 0 );
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : CompositeAttrib destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/30/03
//-------------------------------------------------------------------------
CompositeAttrib::~CompositeAttrib()
{
  delete [] int_array;
  delete [] real_array;
  delete [] string_array;
}

//-------------------------------------------------------------------------
// Purpose       : Compare to a csa
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/30/03
//-------------------------------------------------------------------------
bool CompositeAttrib::equals( CubitSimpleAttrib* attrib ) const
{
  if (attrib->string_data_list()->size() != string_count ||
      attrib->int_data_list()->size() != int_count ||
      attrib->double_data_list()->size() != real_count)
    return false;
  
  int i;
  CubitString* s_itor = string_array;
  attrib->string_data_list()->reset();
  for (i = 0; i < attrib->string_data_list()->size(); i++)
    if (*(attrib->string_data_list()->next(i)) != *(s_itor++))
      return false;
  
  int *i_itor = int_array;
  attrib->int_data_list()->reset();
  for (i = 0; i < attrib->int_data_list()->size(); i++)
    if (*(attrib->int_data_list()->next(i)) != *(i_itor++))
      return false;
  
  double *r_itor = real_array;
  attrib->double_data_list()->reset();
  for (i = 0; i < attrib->double_data_list()->size(); i++)
    if (*(attrib->double_data_list()->next(i)) != *(r_itor++))
      return false;
  
  return true;
}
