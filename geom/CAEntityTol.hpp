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

  CAEntityTol (RefEntity*, const CubitSimpleAttrib&);
    //- create a CAEID from a simple attribute

  CubitStatus actuate();

  CubitStatus update();

  CubitStatus reset(); 
    //- don't need an active reset function since there aren't any
    //- lists on this CA
  
  void merge_owner(CubitAttrib *deletable_attrib);
  CubitSimpleAttrib split_owner();

  CubitSimpleAttrib cubit_simple_attrib();

  int int_attrib_type() {return CA_ENTITY_TOL;}

  double tolerance() {return entityTol;}

  void print();
  
};

CubitAttrib* CAEntityTol_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa);

#endif

