//-------------------------------------------------------------------------
// Filename      : CompositeAttrib.hpp
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

#ifndef COMPOSITE_ATTRIB_HPP
#define COMPOSITE_ATTRIB_HPP

class CubitSimpleAttrib;
class CubitString;
#include <vector>

class CompositeAttrib
{
  private:
  
    std::vector<int> int_array;
    std::vector<double> real_array;
    std::vector<CubitString> string_array;
  
    void append_to_lists( std::vector<CubitString>& strings,
                          std::vector<int>& ints,
                          std::vector<double>& reals ) const;
  
  public:
  
    const CubitString& name() const { return string_array[0]; }
    
    CompositeAttrib* next;
  
    CompositeAttrib( const CubitSimpleAttrib& attrib, CompositeAttrib* next_ptr );
    CompositeAttrib( const CompositeAttrib& copy );
    ~CompositeAttrib();
    
    void append_to_csa( CubitSimpleAttrib& attrib ) const;
    
    bool equals( const CubitSimpleAttrib& attrib ) const;
    
    CubitSimpleAttrib csa( ) const;
    
};
    

#endif
