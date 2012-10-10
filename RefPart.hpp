//- Class:       RefPart
//- Description: The RefPart class contains implementation of the RefPart Class
//-
//- Owner:       Madhan Narayanan

#ifndef REFPART_HPP
#define REFPART_HPP

#include "RefGroup.hpp"
#include "CubitDefines.h"

template <class X> class DLIList;
class CollectionEntity;
class GeometryQueryEngine;
class PartSM;

class RefPart : public RefGroup
{ 	
public:

  RefPart(CollectionEntity*);
  //- Class contructors

  virtual ~RefPart();
  //- Class destructor

  virtual const type_info& entity_type_info() const
  { return typeid(RefPart); }
  //- return the type for this geometryEntity

  static const char* get_class_name()
     {return "Part";}
  virtual const char* class_name() const
     {return get_class_name();}
  
  CubitStatus remove();
  
  PartSM* get_part_sm_ptr() const;
  
  virtual GeometryQueryEngine* get_geometry_query_engine() const;
  //R GeometryQueryEngine*
  //R- A pointer to the geometric modeling engine associated with
  //R- the object.
  //- This function returns a pointer to the geometric modeling engine
  //- associated with the object.
    
};

#endif // REFPART_HPP

