//-------------------------------------------------------------------------
// Filename      : attrib_cubit_owner.hpp
//
// Purpose       : This attribute represents a pointer from an
//                 ACIS ENTITY to a Cubit AcisBridge
//
// Creator       : Greg Neilson
//
// Owner         : Tim Tautges
//-------------------------------------------------------------------------

#ifndef ATTRIB_CUBIT_OWNER_HPP
#define ATTRIB_CUBIT_OWNER_HPP

// ********** BEGIN ACIS INCLUDES             **********
// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********
#include "CubitDefines.h"
#include "attrib_snl.hpp"
#include "AcisTypes.h"
// ********** END CUBIT INCLUDES              **********

// ********** BEGIN MACRO DEFINITIONS         **********
extern int ATTRIB_CUBIT_OWNER_TYPE;
#define ATTRIB_CUBIT_OWNER_LEVEL (ATTRIB_SNL_LEVEL + 1)
// ********** END MACRO DEFINITIONS           **********

// ********** BEGIN FORWARD DECLARATIONS      **********
class AcisBridge;
class TopologyEntity;
template <class X> class DLIList;
class TopologyBridge;
class ENTITY_LIST;

// ********** END FORWARD DECLARATIONS        **********

class ATTRIB_CUBIT_OWNER: public ATTRIB_SNL
{
public:
  
  ATTRIB_CUBIT_OWNER( ENTITY* = NULL, AcisBridge* = NULL);
  
  void set_cubit_owner( AcisBridge* );
  
  AcisBridge* cubit_owner() { return cubitOwnerData; }
  
  AcisBridge* cubit_owner() const { return cubitOwnerData; }
  
  static void cubit_owner(ENTITY_LIST &entity_list,
                          DLIList<AcisBridge*> &tb_list);
  
  static void cubit_owner(ENTITY_LIST &entity_list,
                          DLIList<TopologyBridge*> &tb_list);
  
  static AcisBridge* cubit_owner(ENTITY *ENTITY_ptr);
  
  static TopologyEntity* get_topology_entity(ENTITY* ENTITY_ptr);

  static void set_cubit_owner(ENTITY* entity, AcisBridge *cubit_entity);
  
  static void remove_cubit_owner(ENTITY* ENTITY_ptr,
                                 CubitBoolean children_too = CUBIT_FALSE);

  virtual void split_owner( ENTITY *entity);
  
  virtual void merge_owner( ENTITY *entity, logical delete_this);
  
  virtual void trans_owner( SPAtransf const& );
  
  virtual logical savable() const;
  virtual logical copyable() const;
  
  virtual void to_tolerant_owner( ENTITY *tol_ent );

  static void set_copyable( bool copying_attribs );
  
  ATTRIB_FUNCTIONS(ATTRIB_CUBIT_OWNER, NONE)
    
    private:
  AcisBridge *cubitOwnerData;
  static bool copyingAttribs;
};

#endif


