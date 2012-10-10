//- Class:       RefAssembly
//- Description: The RefAssembly class contains implementation of the RefAssembly Class
//-
//- Owner:       Madhan Narayanan

#ifndef REFASSEMBLY_HPP
#define REFASSEMBLY_HPP

#include "RefGroup.hpp"

template <class X> class DLIList;
class CollectionEntity;
class GeometryQueryEngine;
class AssemblySM;

class RefAssembly : public RefGroup
{ 	
public:

  RefAssembly(CollectionEntity*);
  //- Class contructors

  virtual ~RefAssembly();
  //- Class destructor

  virtual const type_info& entity_type_info() const
  { return typeid(RefAssembly); }
  //- return the type for this geometryEntity

  static const char* get_class_name()
     {return "Assembly";}
  virtual const char* class_name() const
     {return get_class_name();}
  
  CubitStatus remove();
  
  AssemblySM* get_assembly_sm_ptr() const;
  
  virtual GeometryQueryEngine* get_geometry_query_engine() const;
  //R GeometryQueryEngine*
  //R- A pointer to the geometric modeling engine associated with
  //R- the object.
  //- This function returns a pointer to the geometric modeling engine
  //- associated with the object.

};

#endif // REFASSEMBLY_HPP

