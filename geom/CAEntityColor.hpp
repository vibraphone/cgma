//- Class:          CAEntityColor
//- Owner:          Joel Kopp
//- Description:    Cubit attribute for entity colors.
//- Checked by:
//- Version:

#ifndef CA_ENTITY_COLOR_HPP
#define CA_ENTITY_COLOR_HPP

#include "CubitAttrib.hpp"
#include "CubitDefines.h"
#include "CADefines.hpp"

class RefEntity;

class CUBIT_GEOM_EXPORT CAEntityColor : public CubitAttrib
{
private:
 
  int entityColor;
    //- entity color

public:

  virtual ~CAEntityColor();

  CAEntityColor (RefEntity*, const CubitSimpleAttrib& );
    //- create a CAEID from a simple attribute

  CubitStatus actuate();

  CubitStatus update();

  CubitStatus reset() {return CUBIT_SUCCESS;};
    //- don't need an active reset function since there aren't any
    //- lists on this CA
  
  void merge_owner(CubitAttrib *deletable_attrib);

  CubitSimpleAttrib cubit_simple_attrib();

  int int_attrib_type() {return CA_ENTITY_COLOR;}

  int color() {return entityColor;}

  void print();
  
};

CubitAttrib* CAEntityColor_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa);

#endif

