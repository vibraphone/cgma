//-------------------------------------------------------------------------
// Filename      : FacetLoop.hpp
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

class GeometryEntity;

class FacetBody;
class FacetLump;
class FacetShell;
class FacetSurface;
class FacetCoEdge;
class FacetCurve;
class FacetPoint;

class FacetLoop : public LoopSM
{
public :
  
  FacetLoop( Surface *surf_ptr,
             DLIList<CoEdgeSM*> &coedge_list );
    //I- surf_ptr
    //I- A pointer to the set of CoEdges that bound this loop
  
  FacetLoop( DLIList<CoEdgeSM*> &coedge_list );
    //I- A pointer to the set of CoEdges that bound this loop

  virtual ~FacetLoop() ;
    //- The destructor
  
  virtual CubitBoolean is_external() ;
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- Returns CUBIT_TRUE if the Loop is an external Loop and CUBIT_FALSE
    //- otherwise.
  
  virtual LoopType loop_type() ;

  virtual void append_simple_attribute_virt(const CubitSimpleAttrib&);
    //R void
    //I 
    //I- 
    //I- that is to be appended to this OSME object.
    //- The purpose of this function is to append a 
    //- attribute to the OSME. The  is attached to each of the 
    //- underlying solid model entities this one points to.
  
  virtual void remove_simple_attribute_virt(const CubitSimpleAttrib&);
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
  
  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib>&);
  virtual CubitStatus get_simple_attribute(const CubitString& name,
                                           DLIList<CubitSimpleAttrib>&);
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
  
  void get_lumps   ( DLIList<FacetLump   *>& lumps    );
  void get_shells  ( DLIList<FacetShell  *>& shells   );
  void get_coedges ( DLIList<FacetCoEdge *>& coedges  );
  void get_curves  ( DLIList<FacetCurve  *>& curves   );

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

};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

