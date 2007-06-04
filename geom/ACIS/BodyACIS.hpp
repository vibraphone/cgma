//-------------------------------------------------------------------------
// Filename      : BodyACIS.hpp
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

#ifndef BODY_ACIS_HPP
#define BODY_ACIS_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
#if CUBIT_ACIS_VERSION < 1100
#include "baseutil/vector/transf.hxx"
#include "kernel/kerndata/top/body.hxx"
#else
#include "transf.hxx"
#include "body.hxx"
#endif
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"

#include "BodySM.hpp"
#include "AcisBridge.hpp"
#include "AcisTypes.h"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class Body;
class BODY;
class TopologyEntity;
class CubitString;
class insanity_list;
// ********** END FORWARD DECLARATIONS     **********

class BodyACIS : public BodySM, public AcisBridge
{
public:
  
  BodyACIS(BODY* BODYPtr = NULL);
    //- Constructor with a pointer to a ACIS BODY.
  
  virtual ~BodyACIS() ;
    //- The destructor.
  
  BODY* get_BODY_ptr() const;
    //R BODY*
    //R- Pointer to a BODY
    //- Returns a pointer to the ACIS BODY associated with this BodyACIS
    //- object.
  
  static BODY* get_BODY_ptr(Body *body_ptr, CubitBoolean verbose=CUBIT_TRUE);
  static BODY* get_BODY_ptr(const Body *body_ptr );
    //R BODY*
    //R- Pointer to a BODY
    //- Returns a pointer to the ACIS BODY associated with this BodyACIS
    //- object.
  
  virtual GeometryQueryEngine* get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  bool is_sheet_body() const;
  static bool is_sheet_body( BODY *BODY_ptr );
  
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
  

  virtual CubitStatus get_transforms( CubitTransformMatrix &tfm );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I BODYPtr
    //- return the transformation matrix
  
  void set_BODY_ptr(BODY* BODY_ptr);
    //I BODYPtr
    //I- A BODY pointer to be associated with this object.
    // This will replace the BODY previously associated with this object.

  int validate(const CubitString &user_name,
               DLIList <TopologyEntity*> &bad_entities);
    //- does an api_entity_check for the body.
  
  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );

  virtual CubitStatus mass_properties( CubitVector& centroid, double& volumem );
  
  virtual CubitPointContainment point_containment( const CubitVector& pos );

protected: 
  
private:
  
#ifdef BOYD14
  CubitStatus transform_BODY ( SPAtransf transformation_matrix );
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I transformation_matrix
    //I- ACIS transf object representing a transformation matrix
    //- This function applies the input transformation to the transformation
    //- that already exists on the BODY_.  It then "converts"
    //- the transformation applied to the BODY to an identity matrix --
    //- ACIS takes care of modifying the actual local coordinates of
    //- the BODY to "absorb" the transformation.
#endif


  CubitStatus show_bad_geom( DLIList <TopologyEntity*> &ent_list, 
                             const CubitString &user_name);
  // This was put in only to display bad geometry from a validate.  Needs to 
  // be encapsulated in another class...

  void convert_entity_list( insanity_list *ent_list,
                            DLIList <TopologyEntity*> &topology_entities);
    // Given the acis entities, gets the corrisponding topology entities in
    // cubit by which they are represented.
  

};


// ********** BEGIN INLINE FUNCTIONS       **********

inline BODY* BodyACIS::get_BODY_ptr() const
{
  return (BODY*)ENTITY_ptr();
}

// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

