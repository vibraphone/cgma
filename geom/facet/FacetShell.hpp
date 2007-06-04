//-------------------------------------------------------------------------
// Filename      : FacetShell.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef FACET_SHELL_HPP
#define FACET_SHELL_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "ShellSM.hpp"

class FacetBody;
class FacetLump;
class FacetSurface;
class FacetLoop;
class FacetCoEdge;
class FacetCurve;
class FacetPoint;

// ********** END CUBIT INCLUDES           **********

class FacetShell : public ShellSM
{
public:
  
  FacetShell(Lump* my_lump,
             DLIList<Surface*> &my_surfs );
    //- Constructor with lists of attached lumps and surfaces.
  
  FacetShell( DLIList<Surface*> &my_surfs );
    //- Constructor with lists of attached surfaces.
  
  virtual ~FacetShell() ;
    //- Destructor.
  void add_lump(Lump* lump_ptr);
      
  
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

#ifdef BOYD14
  void get_bodies  ( DLIList<FacetBody   *>& bodies   );
#endif
  void get_lumps   ( DLIList<FacetLump   *>& lumps    );
#ifdef BOYD14
  void get_shells  ( DLIList<FacetShell  *>& shells   );
#endif
  void get_surfaces( DLIList<FacetSurface*>& surfaces );
#ifdef BOYD14
  void get_loops   ( DLIList<FacetLoop   *>& loops    );
#endif
  void get_coedges ( DLIList<FacetCoEdge *>& coedges  );
  void get_curves  ( DLIList<FacetCurve  *>& curves   );
#ifdef BOYD14
  void get_points  ( DLIList<FacetPoint  *>& points   );
#endif

  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );
  
  inline Lump* get_lump() const { return myLump; }
  
  inline void remove_lump() { myLump = 0; }
 
  void disconnect_surfaces( DLIList<FacetSurface*> &surfs_to_disconnect );
  void disconnect_all_surfaces();
  
  void reverse(); // invert sense of each surface as used in this shell.
  void reverse_surfaces(); //Actually flip the surface... do not change sense.

  CubitPointContainment point_containment( const CubitVector &point );
  
    //determine whether this is a sheet shell or not.
    // This function may have problems with certain non-manifold geometries
    // It is looking for facets that aren't attached to another facet for
    // one or more of its edges.
    // NOTE (mbrewer): this can probably be improved by going to the
    // curves and checking for the number of co-edges on each curve.
  CubitBoolean is_sheet();

protected: 
  
private:
  Lump* myLump;
  DLIList<Surface*> mySurfs;
};



#endif

