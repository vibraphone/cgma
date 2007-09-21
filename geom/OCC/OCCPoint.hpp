//-------------------------------------------------------------------------
// Filename      : OCCPoint.hpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Steven J. Owen
//
// Creation Date : 08/02/96
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

#ifndef POINT_FACET_HPP
#define POINT_FACET_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "Point.hpp"
#include "OCCAttribSet.hpp"
#include <TopoDS_Vertex.hxx>
#include <gp_Pnt.hxx>

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class TopologyEntity;
class CubitSimpleAttrib;
class RefVertex;
class RefVolume;
class RefVolume;
class CubitPoint;
class OCCAttrib;

class OCCBody;
class OCCLump;
class OCCShell;
class OCCSurface;
class OCCLoop;
class OCCCoEdge;
class OCCCurve;

// ********** END FORWARD DECLARATIONS     **********

class OCCPoint : public Point
{
private:

  DLIList<Curve*> myCurves;
  CubitPoint *myPoint;
  CubitBoolean iCreated;
  TopoDS_Vertex *myTopoDSVertex;

  OCCAttribSet attribSet;
    //List of OCCAttrib*'s instead of CubitSimpleAttribs 

public :
  
  OCCPoint(CubitVector &location, DLIList<Curve*> &curves );
    //I- CubitVector &location
    //I- location of point (creates a CubiPoint).
    //I- DLIList<Curve*> curves
    //I- curves attaced to point

  OCCPoint(CubitPoint *thePoint, DLIList<Curve*> &curves );
    //I- CubitPoint *thePoint
    //I- pointer to the CubitPoint associated with OCCPoint
    //I- DLIList<Curve*> curves
    //I- curves attaced to point
 
  OCCPoint(CubitPoint *thePoint ); 
    //I- CubitPoint *thePoint
    //I- pointer to the CubitPoint associated with OCCPoint

  OCCPoint(TopoDS_Vertex *thePoint, DLIList<Curve*> &curves );
    //I- TopoDS_Vertex *thePoint
    //I- pointer to the TopoDS_Vertex associated with OCCPoint
    //I- DLIList<Curve*> curves
    //I- curves attaced to point
 
  OCCPoint(TopoDS_Vertex *thePoint ); 
    //I- TopoDS_Vertex *thePoint
    //I- pointer to the TopoDS_Vertex associated with OCCPoint

  virtual ~OCCPoint();
    //- The destructor

#ifdef BOYD14
  OCCPoint *copy();
    // make a new copy of this point and return it
#endif
      
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
  
  virtual CubitVector coordinates() const;
    //R CubitVector
    //R- Contains the coordinate values {x y z} of this Point
    //- Returns the spatial coordinates of this Point.
  
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
  CubitStatus move( CubitVector &delta );
    //- Move the first vertex's point by the delta.
#endif

  void add_curve( Curve* curv_ptr )
    { myCurves.append_unique( curv_ptr ); }
    // associate this point with a curve

  CubitPoint *get_cubit_point()
    { return myPoint; }
    // return the CubitPoint associated with this facet point

  CubitStatus save_attribs( FILE* file_ptr );
    // Write FactAttribs out to file

  CubitStatus restore_attribs( FILE* file_ptr, unsigned int endian );
    // Read FactAttribs from file

  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );

#ifdef BOYD14
  void get_bodies  ( DLIList<OCCBody   *>& bodies   );
#endif
  void get_lumps   ( DLIList<OCCLump   *>& lumps    );
  void get_shells  ( DLIList<OCCShell  *>& shells   );
  void get_surfaces( DLIList<OCCSurface*>& surfaces );
  void get_loops   ( DLIList<OCCLoop   *>& loops    );
  void get_coedges ( DLIList<OCCCoEdge *>& coedges  );
  void get_curves  ( DLIList<OCCCurve  *>& curves   );
#ifdef BOYD14
  void get_points  ( DLIList<OCCPoint  *>& points   );
#endif
  
  CubitStatus disconnect_curve( OCCCurve* curve );
  
  inline bool has_parent_curve() const { return myCurves.size() > 0; }
};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif
