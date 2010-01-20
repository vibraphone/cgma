//- Class: ToolData
//- Owner: Scott Mitchell
//- Description: This is tool specific data that should chain off of CubitNode.
//-              This is the generic form, all instances should be
//-              more specific.
//- Checked By: 
//- Version: $Id: 

#ifndef TOOL_DATA_HPP
#define TOOL_DATA_HPP

#include "CubitDefines.h"
#include "CubitUtilConfigure.h"
#include <cassert>

class ToolDataUser;

class CUBIT_UTIL_EXPORT ToolData
{
  private:

    ToolData *nextToolData;
    //- next tool_data in the chain. Last one has this member NULL.

  protected:

    ToolData() {nextToolData = NULL;}
    //- protected constructor to eliminate creation by other than
    //- derived classes

  public:

    virtual ~ToolData();
    //- destructor
 
    ToolData* next_tool_data()               { return nextToolData; }
    void      next_tool_data(ToolData *data) { assert( data != this );
    nextToolData = data; }
    //- access to next tooldata in chain
  
    // handy for ToolDataUser::delete_TD, see e.g. DoubletPillower.cc

        
    virtual ToolData* propogate(ToolDataUser* new_td_user);
    //- propogate() receives the ToolData User that has been copied or split off from the 
    //- ToolDataUser this TD is on and returns the ToolData that should be put on the new 
    //- ToolDataUser.  If no new TD should be created, it returns NULL.

    virtual ToolData* merge(ToolDataUser* other_td_user);
    //- merge() receives a ToolDataUser that is about to be merged with the ToolDataUser that this
    //- TD is on.  It should process what should happen to this and any similar tooldata on the
    //- other ToolDataUser and return the TD that should be assigned to the merged entity.
    //- Note: ToolDataUser deletes any TD that is on it when itself is deleted.  The calling function
    //- should probably remove any TD returned by this function from the entities before deleting them.
  
};
// ********** BEGIN INLINE FUNCTIONS       **********

// ********** END INLINE FUNCTIONS       **********


#endif // TOOL_DATA_HPP








