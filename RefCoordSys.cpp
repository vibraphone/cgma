//-------------------------------------------------------------------------
// Filename      : RefCoordSys.cc
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

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "RefEntityName.hpp"
#include "DLIList.hpp"
#include "CastTo.hpp"
#include "RefCoordSys.hpp"
#include "Csys.hpp"
#include "RefEntityFactory.hpp"

#include "ProeDefines.h"

RefCoordSys::RefCoordSys(Csys* csys_ptr) 
{
   // Set the OtherEntity pointer
   if (csys_ptr != NULL)
   {
      set_other_entity_ptr(csys_ptr) ;
   }
   else
   {
      PRINT_ERROR("ERROR: In the RefCoordSys(Csys*) constructor\n");
      PRINT_ERROR("       Input Csys pointer is NULL\n");
      assert(CUBIT_FALSE);
   }
   
   // Initialize the member data
   initialize();

   // Notify Model about the creation of this object
   CubitObserver::notify_static_observers(this, MODEL_ENTITY_CONSTRUCTED) ;

}

RefCoordSys::~RefCoordSys ()
{
     // remove this entity from the global list
   RefEntityFactory::instance()->remove(this);

   // Notify Model about the destruction of this object
   CubitObserver::notify_static_observers(this, MODEL_ENTITY_DESTRUCTED);
}


CubitStatus RefCoordSys::remove()
{
 
   delete otherEntityPtr;
   
   return CUBIT_SUCCESS;
}

GeometryQueryEngine* RefCoordSys::get_geometry_query_engine() const
{
   return otherEntityPtr->get_geometry_query_engine();   
}

//-------------------------------------------------------------------------
// Purpose       : Returns the global origin of this RefCoordSys.
//
// Special Notes :
//
// Creator       : Eric Nielsen
//
// Creation Date : 03/16/98
//-------------------------------------------------------------------------
CubitVector RefCoordSys::origin() const
{
   Csys* csys_ptr = CAST_TO(otherEntityPtr, Csys);
   
   if (csys_ptr != NULL)
   {
      return csys_ptr->origin();
   }
   
   else
   {
      PRINT_ERROR("ERROR: In RefCoordSys::origin\n");
      PRINT_ERROR("       RefCoordSys %d has no OtherEntity attached to it.\n", id());
      PRINT_ERROR("       THIS IS A BUG - PLEASE REPORT IT.\n");
      assert( csys_ptr != NULL );
      return CubitVector();
   }
}


// ********** BEGIN PRIVATE FUNCTIONS      **********

//-------------------------------------------------------------------------
// Purpose       : Initializes member data
//
// Special Notes :
//
// Creator       : Eric Nielsen
//
// Creation Date : 03/16/98
//-------------------------------------------------------------------------
void RefCoordSys::initialize()
{
   // Set the Entity ID for this new RefCoordSys
   entityId = RefEntityFactory::instance()->next_ref_vertex_id();
}



CubitBox RefCoordSys::bounding_box()
{

   CubitBox dummy;
   PRINT_ERROR("Error:  RefCoordSys::bounding_box()\n");
   PRINT_ERROR("        function not implemented yet\n");
   
   return dummy;      
      
   
}

const char* RefCoordSys::class_name() const
{
	PRINT_ERROR(STUB_ERROR_MSG);
    return CUBIT_FALSE;
}

DagType RefCoordSys::dag_type() const
{
	PRINT_ERROR(STUB_ERROR_MSG);
	return DagType::invalid_type();
}
