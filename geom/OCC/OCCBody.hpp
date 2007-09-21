//-------------------------------------------------------------------------
// Filename      : OCCBody.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 7/18/00
//
//-------------------------------------------------------------------------

#ifndef FACET_BODY_HPP
#define FACET_BODY_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "BodySM.hpp"
#include "CubitTransformMatrix.hpp"
#include "OCCAttribSet.hpp"

#include <TopoDS_CompSolid.hxx>
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class Body;
class TopologyEntity;
class CubitString;
class OCCAttrib;

class OCCLump;
class OCCShell;
class OCCSurface;
class OCCLoop;
class OCCCoEdge;
class OCCCurve;
class OCCPoint;

// ********** END FORWARD DECLARATIONS     **********

class OCCBody : public BodySM
{
public:
  
  OCCBody(TopoDS_Shape *theShape);
  OCCBody(DLIList<Lump*> &myLumps);
    //- Constructor with a pointer to a ACIS BODY.
  virtual ~OCCBody() ;
    //- The destructor.

  CubitBoolean can_be_deleted( DLIList <Body*> &body_list );
  
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
  
  virtual BodySM* copy();
    //R OCCBody*
    //R- Pointer to a OCCBody object
    //- Copies this OCCBody object (including the ACIS BODY that it
    //- contains) and returns a pointer to a new OCCBody object.
  
  virtual CubitStatus move(double , double , double );
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I dx, dy, dz
    //I- Offset values in each of the 3 Cartesian coordinate directions
    //- Move the ACIS BODY by dx, dy and dz
  
  virtual CubitStatus rotate( double , double , double , 
                              double );
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I x, y, z
    //I- Axis of rotation
    //I angle_in_degrees
    //I- Angle of rotation in degrees
    //- Rotate the ACIS BODY angle degrees about a vector defined by 
    //- x, y and z
  
  virtual CubitStatus scale(double, double, double);
  
  virtual CubitStatus scale(double);
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I scaling_factor
    //I- Scaling factor
    //- Scale the ACIS BODY by the factor, scaling_factor

  CubitStatus reflect(double,double,double);
    //- reflect about an axis

  virtual CubitStatus restore();
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //- Restore the ACIS BODY by replacing the transformation matrix 
    //- associated with it with a unit matrix
  
#ifdef BOYD14
  CubitStatus reverse() ;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I BODYPtr
    //- Reverse the face orientations on this body
#endif

  CubitStatus get_transforms( CubitTransformMatrix &tfm );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I BODYPtr
    //- return the transformation matrix for this body
  
  CubitStatus set_transforms( CubitTransformMatrix tfm );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE that myTransforms was 
    // set correctly
    //I BODYPtr

  int validate(const CubitString &, DLIList <TopologyEntity*>&);
    //- does an api_entity_check for the body.
  
  CubitStatus save_attribs( FILE* file_ptr );
    // Write FactAttribs out to file

  CubitStatus restore_attribs( FILE* file_ptr, unsigned int endian );
    // Read FactAttribs from file

#ifdef BOYD14
  void get_bodies  ( DLIList<OCCBody   *>& bodies   );
#endif
  void get_lumps   ( DLIList<OCCLump   *>& lumps    );
  void get_shells  ( DLIList<OCCShell  *>& shells   );
  void get_surfaces( DLIList<OCCSurface*>& surfaces );
  void get_loops   ( DLIList<OCCLoop   *>& loops    );
  void get_coedges ( DLIList<OCCCoEdge *>& coedges  );
  void get_curves  ( DLIList<OCCCurve  *>& curves   );
  void get_points  ( DLIList<OCCPoint  *>& points   );

  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );
  
  void disconnect_all_lumps();
  void add_lump( OCCLump *lump_to_add );
  void remove_lump( OCCLump *lump_to_remove );

  virtual CubitStatus mass_properties( CubitVector& result, double& volume );
  
  virtual CubitPointContainment point_containment( const CubitVector& pos );

protected: 
  
private:
  CubitStatus transform( CubitTransformMatrix &tfmat, CubitBoolean is_rotation );
    // main function for applying transforms to facet-based bodies

  void init_edge_flags( DLIList<Surface *>&surf_list, int flag );
    // set the flags on the facet edges

  CubitTransformMatrix myTransforms;
  DLIList<Lump*> myLumps;
    //List of the attached lumps for the traversal functions.
  OCCAttribSet attribSet;
    //List of OCCAttrib*'s instead of CubitSimpleAttribs 
  TopoDS_Shape *myTopoDSShape;
};



#endif

