//-------------------------------------------------------------------------
// Filename      : attrib_cubit_owner.hpp
//
// Purpose       : COM implementation. This attribute represents a pointer 
//                 from an SW ENTITY to a Cubit SWBridge
//
// Creator       : Joel Kopp, Greg Neilson
//
// Owner         : Tim Tautges
//-------------------------------------------------------------------------

#ifndef SW_ATTRIB_CUBIT_OWNER_HPP
#define SW_ATTRIB_CUBIT_OWNER_HPP

//#include "amapp.h"
#import "C:\Program Files\SolidWorks\sldworks.tlb" raw_interfaces_only, raw_native_types, no_namespace, named_guids     //Import the SolidWorks type library

#include "CubitDefines.h"
//#include "attrib_snl.hpp"

//extern int ATTRIB_CUBIT_OWNER_TYPE;
//#define ATTRIB_CUBIT_OWNER_LEVEL (ATTRIB_SNL_LEVEL + 1)

class SWBridge;
class TopologyEntity;
class TopologyBridge;
#include "DLIList.hpp"


class SW_ATTRIB_CUBIT_OWNER //: public ATTRIB_SNL
{
public:
  static HRESULT initialize();

  static void cubit_owner( DLIList<IUnknown *> &entity_list,
                           DLIList<SWBridge*> &tb_list );

  static void cubit_owner(DLIList<IUnknown *> &entity_list,
                          DLIList<TopologyBridge*> &tb_list);

  static SWBridge* cubit_owner(IUnknown *entity);

  static void cubit_owner( enum EntityType ownerType,
                           DLIList<TopologyBridge*> &tb_list );

  static TopologyEntity* get_topology_entity(IUnknown *entity);

  static void set_cubit_owner(IUnknown *entity, SWBridge *cubit_entity);

  static void remove_cubit_owner(IUnknown *entity,
                                 CubitBoolean children_too = CUBIT_FALSE);

private:
  SW_ATTRIB_CUBIT_OWNER( IUnknown* = NULL, SWBridge* = NULL);
  virtual ~SW_ATTRIB_CUBIT_OWNER();

  void set_cubit_owner( SWBridge* );

  SWBridge* cubit_owner() { return cubitOwnerData; }
  SWBridge* cubit_owner() const { return cubitOwnerData; }

  //virtual void split_owner(LPDISPATCH entityDisp);

  virtual void merge_owner(IUnknown *entity, VARIANT_BOOL delete_this);

  virtual void trans_owner(double transform[16]);

  //ATTRIB_FUNCTIONS(ATTRIB_CUBIT_OWNER, NONE)


private:
  
  static IAttributeDef *attDef;
  static CubitBoolean initialized; // = CUBIT_FALSE;
  static int ownerAttCount; // = 0;

  IAttribute *attInstance;

  SWBridge *cubitOwnerData;

  BSTR bsAttribName;//  CString attrib_name;
  //TCHAR bufferAtt[1000];
};

#endif


