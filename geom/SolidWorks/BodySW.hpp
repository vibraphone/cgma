//-------------------------------------------------------------------------
// Filename      : BodySW.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Joel Kopp, Xuechen Liu
//
// Creation Date : 8/9/00
//
// Owner         : Joel Kopp, Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef BODY_SW_HPP
#define BODY_SW_HPP


#include "CubitDefines.h"
#include "CubitEntity.hpp"

#include "BodySM.hpp"
#include "SWPart.hpp"

class CubitString;
template <class X> class DLIList;


class BodySW : public BodySM
{
public:
  
  BodySW(SWPart *pPart);
    //- Constructor with a pointer to a SW BODY.
  
  virtual ~BodySW() ;
    //- The destructor.
  
  IBody2 *get_BODY_ptr() const;
    //R LPBODY
    //R- Pointer to a BODY
    //- Returns a pointer to the SW BODY associated with this BodySW
    //- object.
  
  virtual GeometryQueryEngine* get_geometry_query_engine() const;
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
  
  virtual CubitStatus get_simple_attribute(const CubitString& name,
                                           DLIList<CubitSimpleAttrib*>&);
  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>&);
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
  
  
  virtual CubitStatus move(double x_offset, double y_offset, double z_offset);
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I dx, dy, dz
    //I- Offset values in each of the 3 Cartesian coordinate directions
    //- Move the SW BODY by dx, dy and dz
  
  virtual CubitStatus rotate( double theta_x, double theta_y, double theta_z,
							  double theta);
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I theta_x, theta_y, theta_z
    //I- Angle of rotation in radians
    //- Rotate the SW BODY theta_x around the x-axis, theta_y around the
    //- y -axis, and theta_z around the z-axis.  Corresponds to a 3-1-3 
    //- rotation matrix.
  
  virtual CubitStatus scale(double scaling_factor);

  virtual CubitStatus scale(double scaling_factor_x,
                            double scaling_factor_y,
                            double scaling_factor_z);
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I scaling_factor
    //I- Scaling factor
    //- Scale the SW BODY by the factor, scaling_factor
  
  virtual CubitStatus restore();
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I- SW Component pointer
    //- Restore the SW BODY by replacing the rotation matrix with its
    //- inverse, the translation with their negative values, and the 
    //- scaling factor with its reciprocal
  
  CubitStatus reverse() ;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //- Reverse the curve orientations on this body
  
  static CubitStatus reverse(IBody2 *body);
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I BODYPtr
    //- Reverse the curve orientations on this body
  
  virtual CubitStatus get_transforms( CubitTransformMatrix &tfm );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I BODYPtr
    //- return the transformation matrix
  
  void set_BODY_ptr(IBody2 *bodyDisp);
    //I LPBODY
    //I- A BODY pointer to be associated with this object.
    // This will replace the BODY previously associated with this object.

  //int validate(const CubitString &user_name);
    //- does an api_entity_check for the body.
  
  void bodysms(DLIList<BodySM*> &bodies);
  void lumps(DLIList<Lump*> &lumps);
  void shellsms(DLIList<ShellSM*> &shellsms);
  void surfaces(DLIList<Surface*> &surfaces);
  void loopsms(DLIList<LoopSM*> &loopsms);
  void curves(DLIList<Curve*> &curves);
  void coedgesms(DLIList<CoEdgeSM*> &coedgesms);
  void points(DLIList<Point*> &points);
    //- topology traversal of TB's; need to implement at this level 'cuz
    //- don't know how many SW entities per TB entity

//    void get_mass_props(CubitVector &cofg);
  virtual CubitStatus mass_properties( CubitVector& centroid, double& volume );
  virtual CubitPointContainment point_containment( const CubitVector& pos );


  virtual void get_parents_virt(DLIList<TopologyBridge*> &parents );
  virtual void get_children_virt(DLIList<TopologyBridge*> &children );

private:
	IBody2 *m_pSWBody;
    SWPart *m_pSWPart;
};


#endif

