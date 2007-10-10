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

#ifndef POINT_OCCPOINT_HPP
#define POINT_OCCPOINT_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "Point.hpp"
#include <stdio.h>
#include "gp_Pnt.hxx"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class CubitSimpleAttrib;
class OCCAttrib;
// ********** END FORWARD DECLARATIONS     **********

class OCCPoint : public Point
{
private:

  gp_Pnt myPoint;

public :
  
  OCCPoint(CubitVector &location );
    //I- CubitVector &location
    //I- location of point (creates a CubiPoint).
    //I- DLIList<Curve*> curves
    //I- curves attaced to point

  OCCPoint(gp_Pnt& thePoint ):myPoint(thePoint){};
    //I- gp_Pnt *thePoint
    //I- pointer to the TopoDS_Vertex associated with OCCPoint
    //I- DLIList<Curve*> curves
    //I- curves attaced to point
 
  gp_Pnt get_gp()const;
  void set_gp( const 
gp_Pnt gp_ptr);
  // - get/set the gp_Pnt associated with this object
  virtual ~OCCPoint();
    //- The destructor

#ifdef BOYD14
  OCCPoint *copy();
    // make a new copy of this point and return it
#endif
      
  virtual CubitVector coordinates() const;
    //R CubitVector
    //R- Contains the coordinate values {x y z} of this Point
    //- Returns the spatial coordinates of this Point.
  
  virtual void set_coord(CubitVector &location );
  virtual void set_coord(double x, double y, double z);

  CubitBoolean is_equal(const OCCPoint & other, double Tol);

  double distance(const OCCPoint & other);
  double SquareDistance (const OCCPoint & other);

  virtual CubitBox bounding_box() const ;
    // see comments in GeometryEntity.hpp
  
  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );

  virtual void append_simple_attribute_virt(CubitSimpleAttrib*);
  virtual void remove_simple_attribute_virt(CubitSimpleAttrib*); 
  virtual void remove_all_simple_attribute_virt();
  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>&);
  virtual CubitStatus get_simple_attribute(const CubitString&, 
                                           DLIList<CubitSimpleAttrib*>&);
  
  CubitStatus save_attribs( FILE *file_ptr );

  CubitStatus restore_attribs( FILE *file_ptr, unsigned int endian );

};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif
