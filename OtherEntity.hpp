//-------------------------------------------------------------------------
// Filename      : OtherEntity.hpp
//
// Purpose       : This file contains the declarations of the base class 
//                 OtherEntity.  This is the root of the class hierarchy
//                 the encapsulates all entities which do NOT fit into 
//                 the GeometryEntity Class such as coordinate systems.
//
// Special Notes : 
//
// Creator       : Eric Nielsen 
//
// Creation Date : 03/012/98
//
// Owner         : Eric Nielsen
//-------------------------------------------------------------------------

#ifndef OTHER_ENTITY_HPP
#define OTHER_ENTITY_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN MOTIF INCLUDES         **********
// ********** END MOTIF INCLUDES           **********

// ********** BEGIN OPEN INVENTOR INCLUDES **********
// ********** END OPEN INVENTOR INCLUDES   **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitMessage.hpp"
#include "CubitSimpleAttrib.hpp"
#include "DLIList.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN MACROS DEFINITIONS     **********
// ********** END MACROS DEFINITIONS       **********

// ********** BEGIN FORWARD DECLARATIONS   **********

class OtherRefEntity;
class GeometryQueryEngine ;

// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN ENUM DEFINITIONS       **********



// ********** END ENUM DEFINITIONS         **********

class OtherEntity 
{
   public :

      OtherEntity() ;
      //- The default constructor

      virtual ~OtherEntity();
      //- The destructor
      //- The destructor is made pure virtual to prevent instantiation
      //- of this class.

      virtual void append_name_attribute(CubitString const& name) = 0;
      //R void
      //I name
      //I- A reference to a constant CubitString object which is the name
      //I- that is to be appended to the OtherEntity object.
      //- The purpose of this function is to append a string
      //- attribute to the entity. 

      virtual void append_simple_attribute(CubitSimpleAttrib*) = 0;
      //R void
      //I name
      //I- A reference to a constant CubitString object which is the name
      //I- that is to be appended to the OtherEntity object.
      //- The purpose of this function is to append a string
      //- attribute to the OE.
      
      virtual void remove_name_attribute( CubitString const& name ) = 0;
      //R void
      //I name
      //I- A reference to a constant CubitString object which is the name
      //I- that is to be removed to this OSME object.
      //- The purpose of this function is to remove a string (name)
      //- attribute from the OE. 
      
      virtual void remove_simple_attribute(CubitSimpleAttrib*) = 0;
      //R void
      //I CubitSimpleAttrib*
      //I- A reference to a CubitSimpleAttrib object which is the object
      //I- that is to be removed to this OE object.
      //- The purpose of this function is to remove a simple
      //- attribute from the OE. 
      
      virtual CubitString get_name_attribute()
      { return CubitString(); }
      //R CubitString
      //R- The returned name
      //- The purpose of this function is to get the "name" attribute
      //- of the other entity. 
      
      //- The default implementation returns an empty string.
      //- MJP Note:
      //- This is the code that implements the requirement that names
      //- of VGI Entities propagate across solid model boolean
      //- operations.  The success of this relies, of course, on the underlying
      //- solid modeler being able to propagate attributes across
      //- such operations on its entities. If it cannot, then "names"
      //- of VGI entities will not propagate.

      virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>&) = 0;
      //R CubitSimpleAttrib*
      //R- the returned cubit simple attribute.
      //- The purpose of this function is to get the attributes
      //- of the other entity. 
      //- The default implementation returns a NULL.
      //- MJP Note:
      //- This is the code that implements the requirement that names
      //- of VGI Entities propagate across solid model boolean
      //- operations.  The success of this relies, of course, on the underlying
      //- solid modeler being able to propagate attributes across
      //- such operations on its entities. If it cannot, then "names"
      //- of VGI entities will not propagate.
      
      virtual CubitStatus set_owner_attribute(
                           OtherRefEntity* otherRefEntityPtr) = 0;
      //R CubitStatus
      //R- CUBIT_FAILURE/SUCCESS
      //I OtherRefEntityPtr
      //I- Pointer to a OtherRefEntity
      //- This (pure virtual) function attaches the input OtherRefEntity 
      //- pointer, as an OWNER attribute, to the underlying solid model 
      //- entity(s) of this OtherEntity object.
      
      virtual OtherRefEntity* get_owner_attribute() = 0;
      //R OtherRefEntity*
      //R- Pointer to the OtherRefEntity using this OtherEntity
      //- This (pure virtual) function returns a pointer to the
      //- OtherRefEntity that is using this object, i.e. owner of this
      //- OtherEntity. It is implied that a OtherEntity is not
      //- used by more than one OtherRefEntity.


      
  
      //*******MJP NOTE: NOT YET SURE IF THIS IS REQUIRED**********
      //virtual void remove_cubit_owner_attributes() {};
      //- Removes all the "Cubit Owner" attributes attached to each of the  
      //- associated solid model entities of this GeometryEntity.  A Cubit 
      //- Owner attribute is a pointer to the Cubit TopologyEntity associated 
      //- with the solid model entities.
      //- Note: A do-nothing implementation is inlined here for those
      //-       types of GeometryEntities that do not have underlying
      //-       solid model representations.
      //*******MJP NOTE: NOT YET SURE IF THIS IS REQUIRED**********

      virtual GeometryQueryEngine* 
                 get_geometry_query_engine() const = 0;
      //R GeometryQueryEngine*
      //R- A pointer to the geometric modeling engine associated with
      //R- the object.
      //- This function returns a pointer to the geometric modeling engine
      //- associated with the object.

      virtual int validate(const CubitString &user_name);
      //- Check that entity is valid. Returns number of problems detected.
      
      
      

      
      
      

   protected:

   
       CubitStatus merge_owner_attribute(OtherEntity* deadOE) ;
      //R CubitStatus
      //R- CUBIT_SUCCESS/CUBIT_FAILURE
      //I deadGE
      //I- The OtherEntity whose owner attribute is about to be merged 
      //I- with this.
      //- This function merges the owner attributes of two OtherEntities.
      //- At the end of this function, the owner of the deadOE will be
      //- the same as owner of this OE. The owner attribute of deadOE is
      //- changed to the owner attribute of this OE.

   
   private:

};

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********
 
// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********
 
// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********
 
#endif

