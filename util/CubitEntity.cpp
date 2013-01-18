//- Class:       CubitEntity
//- Description: Implementation code for the CubitEntity base class.
//-              Once code re-design is complete, this should be made
//-              an abstract class.
//- Owner:       Bill Bohnhoff
//- Checked By:
//- Version: $Id: 


#include <cassert>
#include <cstring>

#include "CubitDefines.h"
#include "CubitVector.hpp"
#include "CubitBox.hpp"
#include "CubitString.hpp"
#include "CubitUtil.hpp"
#include "CubitEntity.hpp"
#include "CubitMessage.hpp"
#include "CastTo.hpp"

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes : 
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/20/96
//-------------------------------------------------------------------------

CubitEntity::~CubitEntity() 
{ 
}

//- Need the following definitions for mesh entities; they should NEVER
//- be called unless there is a logic error

void
CubitEntity::is_visible(int )
{ assert (0);}

int 
CubitEntity::is_visible() const
{ assert (0); return 0; }

void
CubitEntity::is_transparent( int )
{ assert (0); }

int
CubitEntity::is_transparent() const
{ assert (0); return 0; }

int  CubitEntity::id() const { return entityId; }

void CubitEntity::set_id(int i) { entityId = i; }
    //- set the id of this entity to i

int CubitEntity::validate()
{
  // No way to validate a generic CubitEntity.
  // Override in subclasses.
  return 0;
}

void
CubitEntity::color(int)
{ 
    // assert (0);
}

int
CubitEntity::color() const
{ 
    // assert (0); 
  return -1;
}

CubitVector CubitEntity::center_point()
{
  return bounding_box().center();
}
