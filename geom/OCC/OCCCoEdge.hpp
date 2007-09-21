//-------------------------------------------------------------------------
// Filename      : OCCCoEdge.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Steven J. Owen
//
// Creation Date : 07/23/00
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

#ifndef FACETCOEDGE_HPP
#define FACETCOEDGE_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitEntity.hpp"
#include "CoEdgeSM.hpp"
#include <TopoDS_Edge.hxx>
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class TopologyEntity;
class LoopSM;

class OCCBody;
class OCCLump;
class OCCShell;
class OCCSurface;
class OCCLoop;
class OCCCurve;
class OCCPoint;

// ********** END FORWARD DECLARATIONS     **********

class OCCCoEdge : public CoEdgeSM
{
public:
  
  OCCCoEdge(TopoDS_Edge *theEdge, Curve *curv_ptr);
  OCCCoEdge(Curve *curv_ptr, LoopSM *loop_ptr, CubitSense sense);
    //- A constructor
  //
  OCCCoEdge(Curve *curv_ptr, CubitSense sense);
    //- A constructor (for save/restore)
  
  virtual ~OCCCoEdge() ;
    //- The destructor
    
  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
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
  
  CubitSense sense();
    //- returns the sense of the underlying coedge wrt the underlying edge

#ifdef BOYD14
  void get_bodies  ( DLIList<OCCBody   *>& bodies   );
#endif
  void get_lumps   ( DLIList<OCCLump   *>& lumps    );
  void get_shells  ( DLIList<OCCShell  *>& shells   );
#ifdef BOYD14
  void get_surfaces( DLIList<OCCSurface*>& surfaces );
  void get_loops   ( DLIList<OCCLoop   *>& loops    );
  void get_coedges ( DLIList<OCCCoEdge *>& coedges  );
  void get_points  ( DLIList<OCCPoint  *>& points   );
#endif
  void get_curves  ( DLIList<OCCCurve  *>& curves   );

  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );

  void add_loop( LoopSM *loop_ptr )
    { myLoop = loop_ptr; }
    //- set the loop pointer that this coedge is asociated

#ifdef BOYD14
  void set_sense( CubitSense sense );
    //- set the sense of the coedge wrt the surface loop
#endif
  CubitSense get_sense()
    {return edgeSense;}
    //- get the sense of the coedge wrt the surface loop
  void reverse_sense();

  Curve *curve()
    {return myCurve;}
    //- get the curve associated with this coedge
    
  inline LoopSM* get_loop() const { return myLoop; }
  
  inline void remove_loop() { myLoop = 0; }
  inline void remove_curve() { myCurve = 0; }

protected: 
  
private:
  LoopSM *myLoop;
  Curve *myCurve;
  CubitSense edgeSense;
  TopoDS_Edge *myTopoDSEdge;
};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

