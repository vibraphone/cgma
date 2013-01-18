//- Class:          CAEntityName
//- Owner:          Greg Nielson
//- Description:    Cubit attribute for entity names.
//- Checked by:
//- Version:

#ifndef CA_ENTITY_NAME_HPP
#define CA_ENTITY_NAME_HPP

#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CubitDefines.h"
#include "CADefines.hpp"

class RefEntity;

class CUBIT_GEOM_EXPORT CAEntityName: public CubitAttrib
{
private:
 
  DLIList<CubitString> entityNames;

public:

  virtual ~CAEntityName();

  CAEntityName(RefEntity*, const CubitSimpleAttrib&);
    //- create a CAEN from a simple attribute

  CubitStatus actuate();

  CubitStatus update();

  CubitStatus reset();
    //- reset function, cleans out name lists

  CubitSimpleAttrib split_owner();

  void merge_owner(CubitAttrib *deletable_attrib);

  CubitSimpleAttrib cubit_simple_attrib();

  int int_attrib_type() {return CA_ENTITY_NAME;}

  void print();
  
};

CubitAttrib* CAEntityName_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa);

#endif

