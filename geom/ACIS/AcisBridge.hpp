//-------------------------------------------------------------------------
// Filename      : AcisBridge.hpp
//
// Purpose       : Many functions are identical for each Acis-specific
//                 TopologyBridge.  This class implements those functions.
//
// Creator       : Darryl Melander
//
// Creation Date : 03/16/99
//
// Owner         : Darryl Melander
//-------------------------------------------------------------------------

#ifndef ACIS_BRIDGE_HPP
#define ACIS_BRIDGE_HPP

#include "CubitDefines.h"

class ENTITY;
class CubitSimpleAttrib;
class CubitSimpleAttrib;
template <class X> class DLIList;
class Lump;
class ShellSM;
class Surface;
class LoopSM;
class Curve;
class CoEdgeSM;
class Point;
class AcisQueryEngine;
class BodySM;
class CubitString;

class AcisBridge
{
public:
  static unsigned bridgeCount; // Track number of bridges to catch mem leaks.

  AcisBridge(ENTITY* entity);
  virtual ~AcisBridge();
    // Constructor/Destructor.
    // AcisBridge is responsible for connecting/unhooking
    // Acis and Cubit
  
  void unhook_from_ACIS();
    // Removes the links between Cubit and ACIS.
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
  
  CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>& attrib_list);
  static CubitStatus get_simple_attribute(ENTITY *entity,
                                          DLIList<CubitSimpleAttrib*>& attrib_list);
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
    
  CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib*>& result_list );
  static 
  CubitStatus get_simple_attribute( ENTITY* entity, const CubitString& name,
                                    DLIList<CubitSimpleAttrib*>& result_list );
  
  ENTITY* ENTITY_ptr()
    { return myENTITY; }
  const ENTITY* ENTITY_ptr() const
    { return myENTITY; }

protected:
#ifdef BOYD14
  void get_bodies(DLIList<BodySM*> &bodies);
  void get_lumps(DLIList<Lump*> &lumps);
  void get_shells(DLIList<ShellSM*> &shellsms);
  void get_surfaces(DLIList<Surface*> &surfaces);
  void get_loops(DLIList<LoopSM*> &loopsms);
  void get_curves(DLIList<Curve*> &curves);
  void get_coedges(DLIList<CoEdgeSM*> &coedgesms);
  void get_points(DLIList<Point*> &points);
    //- Topology traversal of TBs.  Common implementation for
    //- each type of AcisBridge.
#endif
public:
  
  void ENTITY_ptr(ENTITY* entity);
  
#ifdef BOYD14
  static void print_attribs(ENTITY *);
    //- print all the attributes info for this entity
#endif

  AcisQueryEngine *get_acis_query_engine() const;
    //- return the query engine for this bridge type

private:
  ENTITY* myENTITY;
};

inline CubitStatus AcisBridge::get_simple_attribute(
  DLIList<CubitSimpleAttrib*>& attrib_list)
{
  return get_simple_attribute(myENTITY, attrib_list);
}

inline CubitStatus AcisBridge::get_simple_attribute(
  const CubitString& name, DLIList<CubitSimpleAttrib*>& attrib_list)
{
  return get_simple_attribute(myENTITY, name, attrib_list);
}

#endif
