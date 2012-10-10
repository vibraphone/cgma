//-------------------------------------------------------------------------
// Filename      : RefCoordSys.hpp
//
// Purpose       : Implementation of the RefCsys class.
//                 
//		   
//
// Special Notes :
//
// Creator       : Eric Nielsen
//
// Creation Date : 03/12/98
//
// Owner         : Eric Nielsen
//-------------------------------------------------------------------------


#ifndef REFCOORDSYS_HPP
#define REFCOORDSYS_HPP

#include "DLIList.hpp"
#include "OtherRefEntity.hpp"
#include "CubitBox.hpp"

class Csys;
class GeometryQueryEngine;

class RefCoordSys : public OtherRefEntity
{ 	
public:

  RefCoordSys(Csys*);
  //- Class contructors

  virtual ~RefCoordSys();
  //- Class destructor

  virtual const type_info& entity_type_info() const
  { return typeid(RefCoordSys); }
  //- return the type for this geometryEntity

  CubitStatus remove();
  
  virtual GeometryQueryEngine* get_geometry_query_engine() const;
  //R GeometryQueryEngine*
  //R- A pointer to the geometric modeling engine associated with
  //R- the object.
  //- This function returns a pointer to the geometric modeling engine
  //- associated with the object.

  CubitVector origin () const;
  //R CubitVector
  //R- Contains the origin values {x y z} of this RefCoordSys
  //- Returns the spatial coordinates of this RefCoordSys.
  

  virtual CubitBox bounding_box();
  //- Returns the bounding box of this entity

  virtual const char* class_name() const;
  //- return class name string.

  virtual DagType dag_type() const;

  
  private :
  
        void initialize();
      //- Initializes all member data


  
};

#endif // REFCOORDSYS_HPP

