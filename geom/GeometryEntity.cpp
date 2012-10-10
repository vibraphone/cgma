//-------------------------------------------------------------------------
// Filename      : GeometryEntity.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/31/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitMessage.hpp"
#include "DLIList.hpp"
#include "GeometryEntity.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : Destructor. Does not do anything right now.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/05/96
//-------------------------------------------------------------------------

GeometryEntity::~GeometryEntity()
{}

void GeometryEntity::get_saved_names( DLIList<CubitString*> &names )
{
  std::vector<CubitString>::iterator iter = myNames.begin();
  for(; iter != myNames.end(); iter++ )
  {
    CubitString *tmp_name = &(*iter);    
    names.append( tmp_name );
  }
}

void GeometryEntity::set_saved_names( DLIList<CubitString*> names )
{
  myNames.clear();
  int i;
  for( i=names.size(); i--; )   
    myNames.push_back( *(names.get_and_step()) );  
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
