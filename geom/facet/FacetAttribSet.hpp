//-------------------------------------------------------------------------
// Filename      : FacetAttribSet.hpp
//
// Purpose       : Common attrib functions for MBG
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/01/03
//-------------------------------------------------------------------------

#ifndef FACET_BRIDGE_HPP
#define FACET_BRIDGE_HPP

#include <DLIList.hpp>

class CubitSimpleAttrib;
class FacetAttrib;
class CubitString;

class FacetAttribSet 
{

  public:
  
    FacetAttribSet() : listHead(0) {}
    
    ~FacetAttribSet() { remove_all_attributes(); }
    
    void append_attribute( const CubitSimpleAttrib& );
    
    void remove_attribute( const CubitSimpleAttrib& );
    
    void remove_all_attributes();
    
    CubitStatus get_attributes( DLIList<CubitSimpleAttrib>& ) const;
    
    CubitStatus get_attributes( const CubitString& name,
                                DLIList<CubitSimpleAttrib>& ) const;
    
    CubitStatus save_attributes( FILE* file ) const;
    
    CubitStatus restore_attributes( FILE* file, unsigned int endian );
    
    int attribute_count() const;
    
  private:
  
    FacetAttrib* listHead;
};

#endif
