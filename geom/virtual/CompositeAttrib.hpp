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
template <class X> class DLIList; 

class CompositeAttrib
{
  private:
  
    int int_count;
    int* int_array;
    
    int real_count;
    double* real_array;
    
    int string_count;
    CubitString* string_array;
  
    void append_to_lists( DLIList<CubitString*>& strings,
                          DLIList<int*>& ints,
                          DLIList<double*>& reals ) const;
  
  public:
  
    const CubitString& name() const { return *string_array; }
    
    CompositeAttrib* next;
  
    CompositeAttrib( CubitSimpleAttrib* attrib, CompositeAttrib* next_ptr );
    CompositeAttrib( const CompositeAttrib& copy );
    ~CompositeAttrib();
    
    void append_to_csa( CubitSimpleAttrib* attrib ) const;
    
#ifdef BOYD15
    void merge_strings( const CompositeAttrib& dead );
#endif
    
    bool equals( CubitSimpleAttrib* attrib ) const;
    
    CubitSimpleAttrib* csa( ) const;
    
};
    

#endif
