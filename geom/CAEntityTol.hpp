//- Class:          CAEntityTol 
//- Owner:          W. Roshan Quadros
//- Description:    Cubit attribute for tolerance at entity 
//- Checked by:
//- Version:

#ifndef CA_ENTITY_TOL_HPP
#define CA_ENTITY_TOL_HPP

#include "CubitAttrib.hpp"
#include "CubitDefines.h"
#include "CADefines.hpp"

class RefEntity;

class CUBIT_GEOM_EXPORT CAEntityTol : public CubitAttrib
{
private:
 
  double entityTol;
  //- entity tolerance

 public:

  virtual ~CAEntityTol();

  CAEntityTol (RefEntity*);

  CAEntityTol (RefEntity*, CubitSimpleAttrib *);
    //- create a CAEID from a simple attribute

  //HEADER- RTTI and safe casting functions.
  virtual const type_info& entity_type_info() const
     { return typeid(CAEntityTol);}
  //R- The geometric modeler type
  //- This function returns the type of the geometric modeler.

  CubitStatus actuate();

  CubitStatus update();

  CubitStatus reset(); 
    //- don't need an active reset function since there aren't any
    //- lists on this CA
  
  void merge_owner(CubitAttrib *deletable_attrib);
  CubitSimpleAttrib *split_owner();

  CubitSimpleAttrib* cubit_simple_attrib();

  int int_attrib_type() {return CA_ENTITY_TOL;}

  double tolerance() {return entityTol;}

  void print();
  
};

CubitAttrib* CAEntityTol_creator(RefEntity* entity, CubitSimpleAttrib *p_csa);

#endif

