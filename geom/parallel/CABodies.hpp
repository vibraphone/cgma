//- Class:          CABodies
//- Description:    Cubit attribute for bodies entity is part of.
//- Author: Hong-Jun Kim
//- Version:

#ifndef CA_BODIES_HPP
#define CA_BODIES_HPP

#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CubitDefines.h"
#include "CADefines.hpp"

#include <typeinfo>
#if !defined(_MSC_VER)
using std::type_info;
#endif

class RefEntity;

class CUBIT_GEOM_EXPORT CABodies: public CubitAttrib
{
private:
 
  int m_interface, m_uniqueID;

  DLIList<int> m_sharedBodies; // shared bodies

  DLIList<int> m_sharedProcs; // shared processors

  DLIList<int> m_ghostProcs; // ghost processors

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

  CubitSimpleAttrib cubit_simple_attrib();

  CubitStatus reset();
  //- reset this attribute

  int int_attrib_type() {return CA_BODIES;}
  //- returns the enumerated attribute type
};

CubitAttrib* CABodies_creator(RefEntity* entity, const CubitSimpleAttrib& p_csa);

#endif



