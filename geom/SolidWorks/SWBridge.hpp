//-------------------------------------------------------------------------
// Filename      : SWBridge.hpp
//
// Purpose       : Many functions are identical for each SW-specific
//                 TopologyBridge.  This class implements those functions.
//
// Creator       : Joel Kopp, Darryl Melander
//
// Creation Date : 9/25/00
//
// Owner         : Joel Kopp, Darryl Melander
//-------------------------------------------------------------------------

#ifndef SW_BRIDGE_HPP
#define SW_BRIDGE_HPP

#include "CubitDefines.h"
#include "CubitEntity.hpp"

//struct IEntity;
class CubitSimpleAttrib;

#include "DLIList.hpp"


class SWQueryEngine;
class SWPart;

class SWBridge
{
public:
  SWBridge();
  virtual ~SWBridge();
    // Constructor/Destructor.
    // SWBridge is responsible for connecting/unhooking
    // SolidWorks and Cubit
  
  //void unhook_from_SW();
    // Removes the links between Cubit and SW.
    // ENTITY_ptr() is NULL after this function.
  
  void append_simple_attribute_virt(CubitSimpleAttrib* attrib_ptr);
    //- The purpose of this function is to append an 
    //- attribute to the OSME. The  is attached to each of the 
    //- underlying solid model entities this one points to.
  
  void remove_simple_attribute_virt(CubitSimpleAttrib* attrib_ptr);
    //R void
    //I CubitSimpleAttrib*
    //I- A reference to a CubitSimpleAttrib object which is the object
    //I- that is to be removed to this OSME object.
    //- The purpose of this function is to remove a simple
    //- attribute from the OSME. The attribute is attached to each of the
    //- underlying solid model entities this one points to.
  
  void remove_all_simple_attribute_virt();
    //R void
    //I-
    //- The purpose of this function is to remove all simple
    //- attributes from the OSME. 
  
  CubitStatus get_simple_attribute(const CubitString& name,
                                   DLIList<CubitSimpleAttrib*>& );
  CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>& cubit_simple_attrib_list);
  //static CubitStatus get_simple_attribute(LPENTITY entity,
  //                                        DLCubitSimpleAttribList& attrib_list);
    //R CubitSimpleAttrib*
    //R- the returned cubit simple attribute.
    //- The purpose of this function is to get the attributes
    //- of the geometry entity. The name is attached to the underlying solid
    //- model entity(ies) this one points to.
    //- MJP Note:
    //- This is the code that implements the requirement that names
    //- of VGI Entities propagate across solid model boolean
    //- operations.  The success of this relies, of course, on the underlying
    //- solid modeler being able to propagate attributes across
    //- such operations on its entities. If it cannot, then "names"
    //- of VGI entities will not propagate.
  
  
public:
  
  //static void print_attribs(LPDISPATCH);
    //- print all the attributes info for this entity

private:

  DLIList<CubitSimpleAttrib *> m_simpleAttribList;

protected:
};

#endif
