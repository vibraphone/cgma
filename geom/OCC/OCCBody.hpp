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

#ifndef OCC_BODY_HPP
#define OCC_BODY_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "BodySM.hpp"
#include "CubitTransformMatrix.hpp"
#include "CubitSimpleAttrib.hpp"
#include "CubitBox.hpp"
#include <TopAbs_ShapeEnum.hxx>
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class Body;
class TopologyEntity;
class CubitString;
class OCCLump;
class OCCShell;
class OCCSurface;
class OCCLoop;
class OCCCurve;
class OCCPoint;
class BRepBuilderAPI_Transform;
class BRepAlgoAPI_BooleanOperation;
class BRepBuilderAPI_MakeShape;
class LocOpe_SplitShape;
class TopoDS_Shape;
class TopoDS_Compound;
// ********** END FORWARD DECLARATIONS     **********

class OCCBody : public BodySM
{
public:
  
  // Currently, the compound should only consists solids.
  OCCBody(TopoDS_Compound *theShape, 
          OCCSurface* surface = NULL, OCCShell* shell = NULL, Lump* lump = NULL);

  OCCBody(DLIList<Lump*>& my_lumps, DLIList<OCCShell*>& shells, DLIList<OCCSurface*>& surfaces);
  void lumps(DLIList<Lump*>& my_lumps); //add lump list to myLumps
  DLIList<Lump*> lumps(){return myLumps;} 
  
  void shells(DLIList<OCCShell*> shells) {myShells = shells;}
  void shells(OCCShell* shell) {myShells.append(shell);}
  DLIList<OCCShell*> shells() {return myShells;}

  virtual ~OCCBody() ;
    //- The destructor.

  virtual GeometryQueryEngine* get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  TopoDS_Compound *get_TopoDS_Shape(); 
  void set_TopoDS_Shape( TopoDS_Compound& theshape);

  virtual CubitStatus get_transforms( CubitTransformMatrix &tfm );
  //R CubitStatus
  //R- CUBIT_SUCCESS/CUBIT_FAILURE
  //I BODYPtr
  //- return the transformation matrix for this body

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
  
  void update_bounding_box();
    // calculate bounding box.

  CubitBox get_bounding_box();
    // return bounding box.

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

  virtual void get_parents_virt( DLIList<TopologyBridge*>& parents );
  virtual void get_children_virt( DLIList<TopologyBridge*>& children );
  
  void get_all_surfaces(DLIList<OCCSurface*> &surfaces);
  void get_all_curves(DLIList<OCCCurve*> &curves);
  void get_all_points(DLIList<OCCPoint*> &points);

  virtual CubitStatus mass_properties( CubitVector& result, double& volume );
  
  virtual CubitPointContainment point_containment( const CubitVector& pos );

  //update the underlining CompSolid
  CubitStatus update_OCC_entity( BRepBuilderAPI_Transform *aBRepTrsf,
                                 BRepAlgoAPI_BooleanOperation *op = NULL);

  static CubitStatus update_OCC_entity(TopoDS_Shape& old_shape,
                                       TopoDS_Shape& new_shape,
                                       BRepBuilderAPI_MakeShape *op,
                                       LocOpe_SplitShape* sp = NULL );

  DLIList<OCCSurface*> my_sheet_surfaces() {return mySheetSurfaces;}

  void set_sheet_surfaces(DLIList<OCCSurface*> surfaces); 

  void set_sheet_surfaces(OCCSurface* surf) {mySheetSurfaces.append(surf);}
  virtual CubitBoolean is_sheet_body();

  static TopoDS_Compound* make_Compound(DLIList<Lump*>& my_lumps,
                                        DLIList<OCCShell*>& shells,
                                        DLIList<OCCSurface*>& surfaces);

  CubitStatus transform(BRepBuilderAPI_Transform& aBRepTrsf); 
protected: 
private:

  DLIList<Lump*> myLumps;
    //List of the attached lumps for the traversal functions.
  DLIList<CubitSimpleAttrib*> csa_list;
    //List of CubitSimpleAttribs 
  TopoDS_Compound *myTopoDSShape;

  CubitBox boundingbox;

  DLIList<OCCSurface*> mySheetSurfaces;

  DLIList<OCCShell*>  myShells; 
};


#endif

