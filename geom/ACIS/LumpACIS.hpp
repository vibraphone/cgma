//-------------------------------------------------------------------------
// Filename      : LumpACIS.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef LUMP_ACIS_HPP
#define LUMP_ACIS_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "Lump.hpp"
#include "AcisBridge.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class LUMP ;
class ENTITY_LIST;
class TopologyEntity;
// ********** END FORWARD DECLARATIONS     **********

class LumpACIS : public Lump, public AcisBridge
{
public:
  
  LumpACIS(LUMP* LUMPPtr = NULL);
    //I- LUMP*
    //I- A pointer to the LUMP which the object will be associated with.
    //- This constructor takes a pointer to a LUMP to which it will be
    //- connected.
  
  virtual ~LumpACIS();
    //- The destructor
  
  LUMP* get_LUMP_ptr() const;
  void set_LUMP_ptr(LUMP* LUMP_ptr);
    // Get/set the LUMP associated with this object.
  
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
  
  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I GEPtr
    //I- A pointer to a GeometryEntity object which *should* be of type
    //I- LumpACIS.
    //- This function merges "this" LumpACIS object with the input GEPtr
    //- object.  Checks the type of the input object to make sure it is
    //- another LumpACIS object.
  
  virtual double measure();
    //R double
    //R- The numeric value of the measure (its units depend on the dimension
    //R- of the RefEntity being "measured")
    //- A generic geometric extent function.
    //- Returns volume for Lump, area for Surface, length for Curve and 
    //- 1.0 for Point

  virtual CubitStatus mass_properties( CubitVector &centroid, double &volume );

  virtual int validate(const CubitString &entity_name,
                       DLIList <TopologyEntity*> &bad_entities);
    //R int
    //R- 1 if there are bad entities, 0 if the entity is good.
    //R- Populates bad_entities list with bad entities.
    //I const CubitString
    //I- name of entity (lump) in cubit.
  
  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );

protected: 
  
private:
  void convert_entity_list(ENTITY_LIST &ent_list,
                           DLIList <TopologyEntity*> &topo_entity_list);
    //- From the ent_list, it populates the topo_entity_list with the
    //- related cubit entities from their corrisponding acis entity.
} ;


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

