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
CompositeAttrib::CompositeAttrib( const CubitSimpleAttrib& attrib, CompositeAttrib* n )
  : int_array(0), real_array(0), string_array(0), next(n)
{
  int_array = attrib.int_data_list();
  real_array = attrib.double_data_list();
  string_array = attrib.string_data_list();
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
  : next(0)
{
  int_array = copy.int_array;
  real_array = copy.real_array;
  string_array = copy.string_array;
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
void CompositeAttrib::append_to_csa( CubitSimpleAttrib& attrib ) const
{
  append_to_lists(attrib.string_data_list(),
                  attrib.int_data_list(),
                  attrib.double_data_list());
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
void CompositeAttrib::append_to_lists( std::vector<CubitString>& string_list,
                                       std::vector<int>& int_list,
                                       std::vector<double>& real_list ) const
{
  int_list.insert(int_list.end(), int_array.begin(), int_array.end());
  real_list.insert(real_list.end(), real_array.begin(), real_array.end());
  string_list.insert(string_list.end(), string_array.begin(), string_array.end());
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
CubitSimpleAttrib CompositeAttrib::csa() const
{
  std::vector<CubitString> string_list;
  std::vector<int> int_list;
  std::vector<double> real_list;
  append_to_lists( string_list, int_list, real_list );
  return CubitSimpleAttrib(&string_list, &real_list, &int_list);
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
bool CompositeAttrib::equals( const CubitSimpleAttrib& attrib ) const
{
  if (attrib.string_data_list().size() != string_array.size() ||
      attrib.int_data_list().size() != int_array.size() ||
      attrib.double_data_list().size() != real_array.size())
    return false;
  
  return std::equal(string_array.begin(), string_array.end(), attrib.string_data_list().begin()) &&
      std::equal(int_array.begin(), int_array.end(), attrib.int_data_list().begin()) &&
      std::equal(real_array.begin(), real_array.end(), attrib.double_data_list().begin());
}
