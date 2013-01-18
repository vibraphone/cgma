//- Class:          CABodies
//- Description:    Cubit attribute for bodies entity is part of.
//- Author: Hong-Jun Kim
//- Checked by: Tim Tautges
//- Version:

#ifndef CA_BODIES_HPP
#define CA_BODIES_HPP

#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CubitDefines.h"
#include "CADefines.hpp"

class RefEntity;

class CUBIT_GEOM_EXPORT CABodies: public CubitAttrib
{
private:
 
  DLIList<int> bodyUniqueId;
  //- body Unique ids containing attribOwnerEntity

public:

  virtual ~CABodies();

  CABodies(RefEntity*);

  CABodies(RefEntity*, CubitSimpleAttrib *);
    //- create a CAB from a simple attribute

  virtual const type_info& entity_type_info() const
     { return typeid(CABodies);}
  //R- The geometric modeler type
  //- This function returns the type of the geometric modeler.

  CubitStatus actuate();

  CubitStatus update();

  CubitSimpleAttrib* cubit_simple_attrib();

  CubitStatus reset();
  //- reset this attribute

  int int_attrib_type() {return CA_BODIES;}
  //- returns the enumerated attribute type

  
};

CubitAttrib* CABodies_creator(RefEntity* entity, CubitSimpleAttrib *p_csa);

#endif



