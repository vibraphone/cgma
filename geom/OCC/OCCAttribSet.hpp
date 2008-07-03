//-------------------------------------------------------------------------
// Filename      : OCCAttribSet.hpp
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
#include "TDF_Label.hxx"
#include "TDF_TagSource.hxx"
#include "OCCQueryEngine.hpp"

class CubitSimpleAttrib;
class OCCAttrib;
class CubitString;
class TopoDS_Shape;

class OCCAttribSet 
{

  public:
  
    OCCAttribSet() : listHead(0) 
    {myLabel = TDF_TagSource::NewChild(OCCQueryEngine::mainLabel);}
    
    ~OCCAttribSet() { remove_all_attributes(); }
    
    void append_attribute( CubitSimpleAttrib*, TopoDS_Shape& shape );
    
    void remove_attribute( CubitSimpleAttrib* );
    
    void remove_all_attributes();
    
    CubitStatus get_attributes( DLIList<CubitSimpleAttrib*>& ) const;
    
    CubitStatus get_attributes( const CubitString& name,
                                DLIList<CubitSimpleAttrib*>& ) const;
    
    CubitStatus save_attributes( FILE* file ) const;
    
    CubitStatus restore_attributes( FILE* file, unsigned int endian );
    
    int attribute_count() const;
    
  private:
  
    OCCAttrib* listHead;
    TDF_Label   myLabel;
};

#endif
