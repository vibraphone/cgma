//- Class:          CAEntityId
//- Owner:          Dong Zhu
//- Description:    Cubit attribute for entity Ids.
//- Checked by:
//- Version:

#ifndef CA_ENTITY_ID_HPP
#define CA_ENTITY_ID_HPP

#include "CubitAttrib.hpp"
#include "CubitDefines.h"
#include "CADefines.hpp"

class RefEntity;

class CUBIT_GEOM_EXPORT CAEntityId : public CubitAttrib
{
private:
 
  int entityId;
    //- entity id

  int boundingUid;
    //- uid of bounding entity; for vertices, this entity is the
    //- start vertex

  CubitSense boundingSense;
    //- for surfaces, the sense of the surface wrt the bounding entity

  CubitVector *boundingXYZ;
    //- for single vertex curves, the xyz value of s=1/3; needed to
    //- determine sense for these curves

public:

  virtual ~CAEntityId();

  CAEntityId (RefEntity*, const CubitSimpleAttrib &);
    //- create a CAEID from a simple attribute

  CubitStatus actuate();

  CubitStatus update();

  CubitStatus reset() {return CUBIT_SUCCESS;};
    //- don't need an active reset function since there aren't any
    //- lists on this CA
  
  void merge_owner(CubitAttrib *deletable_attrib);

  CubitSimpleAttrib cubit_simple_attrib();

  int int_attrib_type() {return CA_ENTITY_ID;}

  int id() {return entityId;}

  void print();
  
};

CubitAttrib* CAEntityId_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa);

#endif

