//-------------------------------------------------------------------------
// Filename      : OCCLoop.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Steven J. Owen
//
// Creation Date : 12/06/00
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

#ifndef LOOP_Facet_HPP
#define LOOP_Facet_HPP

#include "CubitDefines.h"
#include "CubitEntity.hpp"
#include "LoopSM.hpp"
#include "DLIList.hpp"

#include <TopoDS_Wire.hxx>

class GeometryEntity;

class OCCBody;
class OCCLump;
class OCCShell;
class OCCSurface;
class OCCCoEdge;
class OCCCurve;
class OCCPoint;

class OCCLoop : public LoopSM
{
public :
 
  OCCLoop( TopoDS_Wire *theWire );	
  OCCLoop( Surface *surf_ptr,
             DLIList<CoEdgeSM*> &coedge_list );
    //I- surf_ptr
    //I- A pointer to the set of CoEdges that bound this loop
  
  OCCLoop( DLIList<CoEdgeSM*> &coedge_list );
    //I- A pointer to the set of CoEdges that bound this loop

  virtual ~OCCLoop() ;
    //- The destructor
  
  virtual void append_simple_attribute_virt(CubitSimpleAttrib*);
    //R void
    //I 
    //I- 
    //I- that is to be appended to this OSME object.
    //- The purpose of this function is to append a 
    //- attribute to the OSME. The  is attached to each of the 
    //- underlying solid model entities this one points to.
  
  virtual void remove_simple_attribute_virt(CubitSimpleAttrib*);
    //R void
    //I CubitSimpleAttrib*
    //I- A reference to a CubitSimpleAttrib object which is the object
    //I- that is to be removed to this OSME object.
    //- The purpose of this function is to remove a simple
    //- attribute from the OSME. The attribute is attached to each of the
    //- underlying solid model entities this one points to.
  
  virtual void remove_all_simple_attribute_virt();
    //R void
    //I-
    //- The purpose of this function is to remove all simple
    //- attributes from the OSME. 
  
  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>&);
  virtual CubitStatus get_simple_attribute(const CubitString& name,
                                           DLIList<CubitSimpleAttrib*>&);
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
  
  virtual CubitBox bounding_box() const ;
    // see comments in GeometryEntity.hpp
  
  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
#ifdef BOYD14
  void get_bodies  ( DLIList<OCCBody   *>& bodies   );
#endif
  void get_lumps   ( DLIList<OCCLump   *>& lumps    );
  void get_shells  ( DLIList<OCCShell  *>& shells   );
#ifdef BOYD14
  void get_surfaces( DLIList<OCCSurface*>& surfaces );
  void get_loops   ( DLIList<OCCLoop   *>& loops    );
#endif
  void get_coedges ( DLIList<OCCCoEdge *>& coedges  );
  void get_curves  ( DLIList<OCCCurve  *>& curves   );
#ifdef BOYD14
  void get_points  ( DLIList<OCCPoint  *>& points   );
#endif

  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );

  void add_surface( Surface *new_surface_ptr )
    { mySurface = new_surface_ptr; }

  void reverse()
    { myCoEdges.reverse(); }
    
  inline Surface* get_surface() const { return mySurface; }
  
  inline void remove_surface() { mySurface = 0; }
  
  void disconnect_all_coedges();

protected: 

private:
  Surface *mySurface;
  DLIList<CoEdgeSM*> myCoEdges;
  TopoDS_Wire *myTopoDSWire;

};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

