//-------------------------------------------------------------------------
// Filename      : GeometryModifyTool.cpp
//
// Purpose       : The tool used to create and modify geometry
//
// Special Notes :
//
// Creator       : Tim Tautges
//
// Creation Date : 2/01
//
//-------------------------------------------------------------------------

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
//#include <iostream>

#include "GeometryModifyTool.hpp"
#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "GeometryEntity.hpp"
#include "GeometryQueryTool.hpp"
#include "MergeTool.hpp"
#include "SurfaceOverlapTool.hpp"
#include "GeometryQueryEngine.hpp"
#include "GeometryModifyEngine.hpp"
#include "AnalyticGeometryTool.hpp"
#include "DAG.hpp"
#include "TopologyBridge.hpp"
#include "ModelQueryEngine.hpp"
#include "CADefines.hpp"
#include "GeomMeasureTool.hpp"

#include "RefEntity.hpp"
#include "RefEntityFactory.hpp"
#include "RefEntityName.hpp"
#include "BasicTopologyEntity.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"

#include "CoVertex.hpp"
#include "CoEdge.hpp"
#include "CoFace.hpp"
#include "CoVolume.hpp"

#include "Chain.hpp"
#include "Loop.hpp"
#include "Shell.hpp"
#include "Body.hpp"

#include "Lump.hpp"
#include "Surface.hpp"
#include "Curve.hpp"
#include "Point.hpp"
#include "BodySM.hpp"
#include "ShellSM.hpp"
#include "LoopSM.hpp"
#include "CoEdgeSM.hpp"

#include "CubitAttrib.hpp"
#include "CubitVector.hpp"
#include "CubitPlane.hpp"
#include "CubitTransformMatrix.hpp"
#include "GfxPreview.hpp"

#include "DLIList.hpp"

#include "DRefFaceArray.hpp"
#include "DRefEdgeArray.hpp"
#include "DRefVertexArray.hpp"


#include "CubitMessage.hpp"
#include "SettingHandler.hpp"

#include "CastTo.hpp"
#include "CpuTimer.hpp"

#include "BridgeManager.hpp"
#include "SplitSurfaceTool.hpp"
#include "OffsetSplitTool.hpp"
#include "AutoMidsurfaceTool.hpp"
#include "TDSurfaceOverlap.hpp"
#include "GfxDebug.hpp"

#include "CubitUndo.hpp"

/* Work around stupid #define hz equal to HZ on IBM */
#ifdef hz
#  undef hz
#endif
#ifdef PROE
#include "CompositeTool.hpp"
#endif

GeometryModifyTool* GeometryModifyTool::instance_ = 0;
CubitBoolean GeometryModifyTool::allEdgesImprint = CUBIT_TRUE;
CubitBoolean GeometryModifyTool::groupImprint = CUBIT_TRUE;
CubitBoolean GeometryModifyTool::newIds = CUBIT_FALSE;
CubitBoolean GeometryModifyTool::sepAfterWebcut = CUBIT_TRUE;
CubitBoolean GeometryModifyTool::booleanRegularize = CUBIT_TRUE;
CubitBoolean GeometryModifyTool::uniteMixedModels = CUBIT_TRUE;
CubitBoolean GeometryModifyTool::oldNames = CUBIT_FALSE;
CubitBoolean GeometryModifyTool::meshAutodelete = CUBIT_TRUE;
CubitBoolean GeometryModifyTool::meshAutodeleteRemesh = CUBIT_FALSE;
RefEntity* GeometryModifyTool::copyingEntity = NULL;

//-------------------------------------------------------------------------
// Purpose       : Controls access and creation of the sole instance of this
//                 class.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------

GeometryModifyTool* GeometryModifyTool::instance(GeometryModifyEngine *GMEPtr)
{
     // Check to see if we have created an instance of the class
     // If not, proceed to create one.

   if (instance_ == 0)
   {
        // When creating the instance, we should always have a valid
        // SMEPtr. If not, complain.

     instance_ = new GeometryModifyTool (GMEPtr) ;

        // check to make sure there's a ref entity factory extant
     //RefEntityFactory *factory =
     RefEntityFactory::instance();
   }

     // If there is an existing instance of the class, check if there
     // was a request to set default solid modeling engine. If so, be nice
     // to the calling routine :) :) and kindly set the default solid
     // modeling engine.

   else if ( GMEPtr != NULL && instance_->gmeList.move_to(GMEPtr)) {
     delete instance_->gmeList.remove();
     instance_->gmeList.insert(GMEPtr);
   }

     // Return the a pointer to the instance of the class.

   return instance_ ;
}

//-------------------------------------------------------------------------
// Purpose       : Destructor.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------
GeometryModifyTool::~GeometryModifyTool ()
{
   instance_ = NULL;
     //Kill the geometry modify engine(s)
   int i;
   for (i = gmeList.size(); i > 0; i--)
     delete gmeList.get_and_step();

   gmeList.clean_out();
}

//-------------------------------------------------------------------------
// Purpose       : Deletes instance variable
//
// Special Notes :
//
// Creator       : Corey Ernst
//
// Creation Date : 12/31/07
//-------------------------------------------------------------------------
void GeometryModifyTool::delete_instance()
{
  if( NULL != instance_ )
  {
    delete instance_;
    instance_ = NULL;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a sphere with the
//                 given radius
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::sphere(double radius,
                                 int x_shift,
                                 int y_shift,
                                 int z_shift,
                                 double inner_radius,
                                 bool delete_side )
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // First make sure the radius is not zero or less
   if ( radius <= 0.0 )
   {
      PRINT_ERROR("In GeometryModifyTool::sphere\n"
                  "       Cannot make a sphere of radius %f\n",
                  radius);
      return NULL ;
   }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state();

  GeometryModifyEngine *gme = gmeList.get();

  BodySM* body_sm = gme->sphere(radius);

  if (x_shift != 0 || y_shift != 0 || z_shift != 0)
  {
    BodySM *brick = gme->brick( 2.0 * radius,
                                2.0 * radius,
                                2.0 * radius);

    CubitVector delta( x_shift * radius,
                       y_shift * radius,
                       z_shift * radius );

    gme->get_gqe()->translate( brick, delta );

    DLIList<BodySM*> new_sms;
    DLIList<BodySM*> from_bodies(1);
    from_bodies.append( body_sm );
    bool keep_old = true;
    CubitStatus bool_status;
    if( delete_side == false )
       bool_status = gme->intersect( brick, from_bodies, new_sms, keep_old );
    else
    {
      DLIList<BodySM*> tool_bodies(1);
      tool_bodies.append( brick );
      bool imprint = false;
      bool_status = gme->subtract( tool_bodies, from_bodies, new_sms, imprint, keep_old );
    }

    //delete the old bodies
    gme->get_gqe()->delete_solid_model_entities( body_sm );
    gme->get_gqe()->delete_solid_model_entities( brick );

    if( bool_status == CUBIT_FAILURE || new_sms.size() == 0 )
    {
      PRINT_ERROR("Problems creating sphere.\n" );
      return NULL ;
    }

    body_sm = new_sms.get();
  }

  Body *new_body = NULL;
  BodySM* inner_body_sm = NULL;
  if( inner_radius )
  {
    inner_body_sm = gme->sphere(inner_radius);
    DLIList<BodySM*> new_sms;
    DLIList<BodySM*> from_bodies(1);
    DLIList<BodySM*> tool_bodies(1);
    from_bodies.append( body_sm );
    tool_bodies.append( inner_body_sm );
    bool imprint = false;
    bool keep_old = true;
    CubitStatus subtract_status =
       gme->subtract( tool_bodies, from_bodies, new_sms, imprint, keep_old );

     //delete the old bodies
    gme->get_gqe()->delete_solid_model_entities( body_sm );
    gme->get_gqe()->delete_solid_model_entities( inner_body_sm );

    if( subtract_status == CUBIT_FAILURE || new_sms.size() == 0 )
    {
      PRINT_ERROR("Problems creating sphere with inner radius.\n" );
      return NULL ;
    }

    body_sm = new_sms.get();
  }

  if (!body_sm)
  {
     PRINT_ERROR("In GeometryModifyTool::sphere\n"
                 "       Problems building a volume from the sphere.\n");
  }
  else
    new_body = GeometryQueryTool::instance()->make_Body(body_sm);

  if( CubitUndo::get_undo_enabled() )
  {
    if( new_body )
      CubitUndo::note_result_body( new_body );
    else
      CubitUndo::remove_last_undo();
  }

  return new_body;
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a brick with given width,
//                 depth, and height.
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::brick(double width, double depth, double height)
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // First make sure that none of the input values are zero or less
   if ( width <= GEOMETRY_RESABS || depth <= GEOMETRY_RESABS || height <= GEOMETRY_RESABS )
   {
      PRINT_ERROR("In GeometryModifyTool::brick\n"
                  "    Cannot make a cuboid of size %f x %f x %f\n"
                  "       All dimensions must be > 0.0\n",
                  width, depth, height);
      return NULL;
   }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state();

   Body *new_body = NULL;
   BodySM* bodyPtr = gmeList.get()->brick(width, depth, height);

   if (bodyPtr == NULL)
   {
      PRINT_ERROR("In GeometryModifyTool::brick\n"
                  "       Problem creating a brick.\n") ;
   }
   else
     new_body = GeometryQueryTool::instance()->make_Body(bodyPtr);

  if( CubitUndo::get_undo_enabled() )
  {
    if( new_body )
      CubitUndo::note_result_body( new_body );
    else
      CubitUndo::remove_last_undo();
  }

   return new_body;
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a brick at given location,
//                 orientation, and size.
//
// Special Notes : If one of the dimensions is zero, a planar sheet is created.
//
// Creator       : Steve Storm
//
// Creation Date : 10/09/2000
//-------------------------------------------------------------------------
Body* GeometryModifyTool::brick(const CubitVector &center,
                                const CubitVector axes[3],
                                const CubitVector &extension)
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

   double width = 2.0*extension.x();
   double height = 2.0*extension.y();
   double depth = 2.0*extension.z();
     // First make sure that entity creation is possible.  Allow at most
     // only one of the dimensions to be zero - in which case a planar
     // sheet is created.
   if ( width < 0.0 || height < 0.0 || depth < 0.0 )
   {
      PRINT_ERROR( "Cannot make a brick of size %f x %f x %f\n"
                   "     Negative dimensions are not allowed.\n",
                   width, height, depth );
      return NULL ;
   }

   int wz = width < GEOMETRY_RESABS;
   int hz = height < GEOMETRY_RESABS;
   int dz = depth < GEOMETRY_RESABS;
   if( wz || hz || dz )
   {
      if( (wz + hz + dz) > 1 )
      {
         PRINT_ERROR( "Cannot make a sheet of size %f x %f x %f\n"
            "     At least two dimensions must be nonzero.\n",
            width, height, depth );
         return NULL ;
      }

      // Make a sheet body instead of a cuboid
      CubitVector sheet_axes[2];
      if( wz )
      {
         sheet_axes[0] = axes[1];
         sheet_axes[1] = axes[2];
         return planar_sheet( center, sheet_axes, height, depth );
      }
      else if( hz )
      {
         sheet_axes[0] = axes[2];
         sheet_axes[1] = axes[0];
         return planar_sheet( center, sheet_axes, depth, width );
      }
      else
      {
         sheet_axes[0] = axes[0];
         sheet_axes[1] = axes[1];
         return planar_sheet( center, sheet_axes, width, height );
      }
   }
   else
   {
      if( CubitUndo::get_undo_enabled() )
        CubitUndo::save_state();

      BodySM* bodyPtr = gmeList.get()->brick(center, axes, extension) ;
      Body *new_body = NULL;
      if (bodyPtr == NULL)
      {
         PRINT_ERROR("In GeometryTool::brick\n"
            "       Problem creating a brick.\n") ;
      }
      else
        new_body = GeometryQueryTool::instance()->make_Body(bodyPtr);

      if( CubitUndo::get_undo_enabled() )
      {
        if( new_body )
          CubitUndo::note_result_body( new_body );
        else
          CubitUndo::remove_last_undo();
      }

      return new_body;
   }
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a prism with height,
//                 sides, major, minor radii
//
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::prism( double height, int sides, double major, double minor)
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // First make sure that the input values make sense
   if ( major < minor ||
        major <= GEOMETRY_RESABS || minor <= GEOMETRY_RESABS ||
        height <= GEOMETRY_RESABS || sides < 3 )
   {
      PRINT_ERROR("In GeometryModifyTool::prism\n"
                  "   Cannot make a prism of major-radius = %f,"
                  " minor-radius = %f, height = %f, and"
                  " number of sides = %d\n"
                  "   All dimensions must be > 0.0, the major-radius\n"
                  " should be greater than or equal to the minor-"
                  "radius,\n and the number of sides should be greater"
                  " than 2.\n", major, minor, height, sides);

      return NULL;
   }

   if( CubitUndo::get_undo_enabled() )
     CubitUndo::save_state();

     // Create a Body that represents the prism
   BodySM* bodyPtr = gmeList.get()->prism(height,  sides, major, minor) ;
   Body *new_body = NULL;
   if (bodyPtr == NULL)
   {
      PRINT_ERROR("In GeometryModifyTool::prism\n"
                  "       Problems building a volume from the prism.\n") ;
   }
   else
     new_body = GeometryQueryTool::instance()->make_Body(bodyPtr);

  if( CubitUndo::get_undo_enabled() )
  {
    if( new_body )
      CubitUndo::note_result_body( new_body );
    else
      CubitUndo::remove_last_undo();
  }

   return new_body;
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a pyramid with given
//                 height, major, minor,top, sides
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::pyramid  ( double height, int sides, double major,
                             double minor,  double top)
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // First make sure that the input values make sense
   if ( major < minor ||
        major <= GEOMETRY_RESABS || minor <= GEOMETRY_RESABS ||
        height <= GEOMETRY_RESABS || top < 0.0 || sides < 3  )
   {
      PRINT_ERROR("In GeometryModifyTool::pyramid\n"
                  "    Cannot make a pyramid of major-radius = %f,"
                  " minor-radius = %f, height = %f,"
                  " top %f, and number of sides = %d\n"
                  "     All dimensions must be > 0.0, the major-radius"
                  " should be greater than the minor-radius, the "
                  "number of sides should be greater than 2.\n",
                  major, minor, height, top, sides);
      return NULL;
   }

   if( CubitUndo::get_undo_enabled() )
     CubitUndo::save_state();

     // Create a Body that represents the prism
   BodySM* bodyPtr = gmeList.get()->pyramid ( height, sides, major, minor, top);
   Body *new_body = NULL;
   if (bodyPtr == NULL)
   {
      PRINT_ERROR("In GeometryModifyTool::pyramid\n"
                  "      Problems building a volume from the pyramid.\n");
   }
   else
     new_body = GeometryQueryTool::instance()->make_Body(bodyPtr);

  if( CubitUndo::get_undo_enabled() )
  {
    if( new_body )
      CubitUndo::note_result_body( new_body );
    else
      CubitUndo::remove_last_undo();
  }

   return new_body;
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a frustum of a cone with
//                 height,
//                 input radius in x-direction at base
//                 input radius in y-direction at base
//                 input radius in x-direction at top
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::cylinder ( double hi, double r1, double r2, double r3 )
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // First make sure that the input values make sense
   if ( hi <= GEOMETRY_RESABS || r1 <= GEOMETRY_RESABS ||
        r2 <= GEOMETRY_RESABS || r3 < 0.0 )
   {
      PRINT_ERROR("In GeometryModifyTool::cylinder\n"
                  " Cannot make a frustum of a cone with height = %f,"
                  " lower x-direction radius = %f,"
                  " lower y-direction radius = %f, and "
                  "top radius = %f\n"
                  "       All dimensions must be > 0.0\n",
                  hi, r1, r2, r3);
      return NULL;
   }

   if( CubitUndo::get_undo_enabled() )
     CubitUndo::save_state();

     // Create a Body that represents the prism
   BodySM* bodyPtr = gmeList.get()->cylinder( hi, r1, r2, r3);
   Body *new_body = NULL;

   if (bodyPtr == NULL)
   {
      PRINT_ERROR("In GeometryModifyTool::cylinder\n"
                  "       Problems building a volume from the conical frustum.\n");
   }
   else
     new_body = GeometryQueryTool::instance()->make_Body(bodyPtr);

  if( CubitUndo::get_undo_enabled() )
  {
    if( new_body )
      CubitUndo::note_result_body( new_body );
    else
      CubitUndo::remove_last_undo();
  }

   return new_body;
}

//-------------------------------------------------------------------------
// Purpose       : Creates a Body corresponding to a torus with major_radius
//                 and minor_radius
//
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 09/18/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::torus( double r1, double r2 )
{   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }


     // First make sure that the input values make sense
   if ( r1 <= r2 || r1 <= GEOMETRY_RESABS || r2 <= 0.0)
   {
      PRINT_ERROR("In GeometryModifyTool::torus\n"
                  "  Cannot make a torus of major-radius = %f and"
                  " minor-radius = %f\n",
                  r1, r2);
      return NULL;
   }

   if( CubitUndo::get_undo_enabled() )
     CubitUndo::save_state();

     // Create a Body that represents the torus
   BodySM* bodyPtr = gmeList.get()->torus(r1, r2) ;
   Body *new_body = NULL;

   if (bodyPtr == NULL)
   {
      PRINT_ERROR("In GeometryModifyTool::torus\n"
                  "       Problems building a volume from the torus.\n") ;
   }
   else
     new_body = GeometryQueryTool::instance()->make_Body(bodyPtr);

  if( CubitUndo::get_undo_enabled() )
  {
    if( new_body )
      CubitUndo::note_result_body( new_body );
    else
      CubitUndo::remove_last_undo();
  }

   return new_body;
}

//-------------------------------------------------------------------------
// Purpose       : Creates a planar sheet with the given corner locations.
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 05/10/05
//-------------------------------------------------------------------------
Body* GeometryModifyTool::planar_sheet( const CubitVector& p1,
                                        const CubitVector& p2,
                                        const CubitVector& p3,
                                        const CubitVector& p4 )
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

   if( CubitUndo::get_undo_enabled() )
     CubitUndo::save_state();

  // Create a Body that represents the sheet
  BodySM* body_ptr = gmeList.get()->planar_sheet(p1, p2, p3, p4) ;
  Body *new_body = NULL;
  if( body_ptr == NULL )
  {
    PRINT_ERROR("In GeometryTool::planar_sheet\n"
      "       Problems building a volume from the sheet.\n") ;
  }
  else
    new_body = GeometryQueryTool::instance()->make_Body(body_ptr);

  if( CubitUndo::get_undo_enabled() )
  {
    if( new_body )
      CubitUndo::note_result_body( new_body );
    else
      CubitUndo::remove_last_undo();
  }

  return new_body;
}

//-------------------------------------------------------------------------
// Purpose       : Creates a planar sheet with the given input location,
//                 orientation and size.
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 10/09/00
//-------------------------------------------------------------------------
Body* GeometryModifyTool::planar_sheet ( const CubitVector &center,
                                         const CubitVector axes[2],
                                         double width, double height )
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

   CubitVector p1, p2, p3, p4;

   // Get the corners of the sheet
   center.next_point( axes[0], width/2.0, p1 );
   p1.next_point( axes[1], -height/2.0, p1 );
   p1.next_point( axes[1], height, p2 );
   p2.next_point( axes[0], -width, p3 );
   p3.next_point( axes[1], -height, p4 );

   if( CubitUndo::get_undo_enabled() )
     CubitUndo::save_state();

   // Create a Body that represents the sheet
   BodySM* body_ptr = gmeList.get()->planar_sheet(p1, p2, p3, p4) ;
   Body *new_body = NULL;

   if( body_ptr == NULL )
   {
      PRINT_ERROR("In GeometryTool::planar_sheet\n"
         "       Problems building a volume from the sheet.\n") ;
   }
   else
     new_body = GeometryQueryTool::instance()->make_Body(body_ptr);

  if( CubitUndo::get_undo_enabled() )
  {
    if( new_body )
      CubitUndo::note_result_body( new_body );
    else
      CubitUndo::remove_last_undo();
  }

  return new_body;
}

//-------------------------------------------------------------------------
// Purpose       : Creates a planar sheet with the minimal amount of area
//                 to cut through the given 3D bounding box.  The sheet
//                 can be extended optionally to make it larger than this
//                 minimal size.
//
// Special Notes : extension_types: 0-none, 1-percentage, 2-absolute
//
// Creator       : Steve Storm
//
// Creation Date : 10/19/98
//-------------------------------------------------------------------------
Body* GeometryModifyTool::planar_sheet ( const CubitPlane& plane,
                                         const CubitBox& bounding_box,
                                         int extension_type,
                                         double extension )
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

   CubitVector p1, p2, p3, p4;

   // Get the corners of the sheet
   if( AnalyticGeometryTool::instance()->
      min_pln_box_int_corners( plane, bounding_box, extension_type,
                               extension, p1, p2, p3, p4 ) == CUBIT_FAILURE )
      return NULL;

   if( CubitUndo::get_undo_enabled() )
     CubitUndo::save_state();

   // Create a Body that represents the sheet
   BodySM* body_ptr = gmeList.get()->planar_sheet(p1, p2, p3, p4) ;
   Body *new_body = NULL;

   if( body_ptr == NULL )
   {
      PRINT_ERROR("In GeometryModifyTool::planar_sheet\n"
         "       Problems building a volume from the sheet.\n") ;
   }
   else
     new_body = GeometryQueryTool::instance()->make_Body(body_ptr);

   if( CubitUndo::get_undo_enabled() )
   {
     if( new_body )
       CubitUndo::note_result_body( new_body );
     else
       CubitUndo::remove_last_undo();
   }

   return new_body;
}

RefVertex* GeometryModifyTool::make_RefVertex( RefVertex *vertex ) const
{
  if ( vertex == NULL )
  {
    PRINT_ERROR("Vertex is NULL\n");
    return NULL;
  }

  TopologyBridge* bridge = 0;
  GeometryModifyEngine* engine = get_engine(vertex, &bridge);
  if (engine == NULL)
  {
     PRINT_ERROR( "%s (vertex %d) does not have a modify engine.\n",
                  vertex->entity_name().c_str(),
                  vertex->id() );
     return 0;
  }
  
  CubitVector point = vertex->coordinates();

  // Call the default GeometryModifyEngine to create a new Point
  Point* point_ptr = engine->make_Point(point);
  
  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state();

  // Use the Point to create a RefVertex
  RefVertex* ref_vertex_ptr = RefEntityFactory::instance()->construct_RefVertex(point_ptr) ;

  if( CubitUndo::get_undo_enabled() )
  {
    if( ref_vertex_ptr )
      CubitUndo::note_result_entity( ref_vertex_ptr );
    else
      CubitUndo::remove_last_undo();
  }

  //transfer the names
  DLIList<CubitString*> names;
  vertex->entity_names( names );

  int i;
  for( i=names.size(); i--; )
  {
    CubitString *tmp_name = names.get_and_step();
    ref_vertex_ptr->entity_name( *tmp_name ); 
  }

  // Send a message to the model indicating the vertex was created
  CubitObserver::notify_static_observers(ref_vertex_ptr, FREE_REF_ENTITY_GENERATED);

  // Return the newly created RefVertex
  return ref_vertex_ptr ;
}


//-------------------------------------------------------------------------
// Purpose       : This function takes a type information and a location
//                 to create a RefVertex. The underlying representation
//                 of the RefVertex is determined by the default
//                 GeometryModifyEngine.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/27/97
//-------------------------------------------------------------------------

RefVertex* GeometryModifyTool::make_RefVertex(
                                        CubitVector const& point, int color) const
{
   if (0 == gmeList.size())
   {
      PRINT_WARNING("No active geometry engine.\n");
      return NULL;
   }

     // Call the default GeometryModifyEngine to create a new Point
   Point* point_ptr = gmeList.get()->make_Point(point);

     // If we get a NULL pointer, give a warning message and return
     // a NULL pointer.
   if ( point_ptr == NULL )
   {
      PRINT_WARNING("In GeometryModifyTool::make_RefVertex\n"
                    "         Got a NULL pointer to a Point.\n"
                    "         Cannot make a RefVertex.\n") ;
      return (RefVertex *)NULL ;
   }

   if( CubitUndo::get_undo_enabled() )
     CubitUndo::save_state();

     // Use the Point to create a RefVertex
   RefVertex* ref_vertex_ptr = RefEntityFactory::instance()->construct_RefVertex(point_ptr) ;

  if( CubitUndo::get_undo_enabled() )
  {
    if( ref_vertex_ptr )
      CubitUndo::note_result_entity( ref_vertex_ptr );
    else
      CubitUndo::remove_last_undo();
  }

   ref_vertex_ptr->color(color);

     // Send a message to the model indicating the vertex was created
   CubitObserver::notify_static_observers(ref_vertex_ptr, FREE_REF_ENTITY_GENERATED);

     // Return the newly created RefVertex
   return ref_vertex_ptr ;
}

Body *GeometryModifyTool::make_Body(Surface *surface) const
{
  // make a Body from a Surface
  BodySM *body_sm = surface->bodysm();

  GeometryModifyEngine* gme = get_engine(surface);

  if (gme && !body_sm)
    body_sm = gme->make_BodySM(surface);

  assert(body_sm != 0);

  return GeometryQueryTool::instance()->make_Body(body_sm);
}

//-------------------------------------------------------------------------
// Purpose       : This function takes RefEdge type information, too
//
// Creator       : David White
//
// Creation Date : 10/9/97
//-------------------------------------------------------------------------
RefEdge* GeometryModifyTool::make_RefEdge(RefEdge *ref_edge_ptr,
                                          bool copy_attribs ) const
{
  if ( ref_edge_ptr == NULL )
  {
    PRINT_ERROR("curve is NULL\n");
    return (RefEdge*)NULL;
  }

  TopologyBridge* bridge = 0;
  GeometryModifyEngine* engine = get_engine(ref_edge_ptr, &bridge);
  Curve *old_curve = dynamic_cast<Curve*>(bridge);
  if (engine == NULL)
  {
     PRINT_ERROR( "%s (curve %d) does not have a modify engine.\n",
                  ref_edge_ptr->entity_name().c_str(),
                  ref_edge_ptr->id() );
     return 0;
  }

  Curve *tmp_curve;
  if( copy_attribs )
  {
    DLIList<RefEntity*> tmp_list;
    tmp_list.append( ref_edge_ptr );
    TopologyBridge *curve_bridge;
    prepare_for_copy( ref_edge_ptr, curve_bridge );
    tmp_curve = CAST_TO( curve_bridge, Curve);
  }
  else
    tmp_curve = old_curve;

  Curve *new_curve = engine->make_Curve( tmp_curve );
  if (!new_curve)
  {
    if( copy_attribs )
      clean_up_from_copy_failure( old_curve );
    return (RefEdge *)NULL;
  }

  TopologyBridge *new_curve_bridge;
  if( copy_attribs )
  {
    new_curve_bridge = new_curve;
    finish_copy( new_curve_bridge, old_curve );
    new_curve = CAST_TO( new_curve_bridge, Curve);
  }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state();

  // Complete the task of linking this new Curve into the rest of the
  // geometry datastructures and return the new RefEdge.
  RefEdge *new_ref_edge = GeometryQueryTool::instance()->make_free_RefEdge(new_curve);

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::note_result_entity( new_ref_edge );

  return new_ref_edge;
}

CubitStatus GeometryModifyTool::prepare_for_copy( RefEntity *ref_ent,
                                                  TopologyBridge *&top_bridge )
{
  //save attribute settings
  CGMApp::instance()->save_current_attribute_states();

  CGMApp::instance()->attrib_manager()->auto_flag(1);
  // Groups are saved directly in the Cubit file, not with attributes
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_GROUP, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_GROUP, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_GROUP, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_GROUP, CUBIT_FALSE);

  // The mesh container is only used for the Exodus save/resore method
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_MESH_CONTAINER, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_MESH_CONTAINER, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_MESH_CONTAINER, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_MESH_CONTAINER, CUBIT_FALSE);

  // Genesis Entities are saved directly in the Cubit file, not with attributes
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_GENESIS_ENTITY, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_GENESIS_ENTITY, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_GENESIS_ENTITY, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_GENESIS_ENTITY, CUBIT_FALSE);

  //Don't save out entity ids
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_ENTITY_ID, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_ENTITY_ID, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_ENTITY_ID, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_ENTITY_ID, CUBIT_FALSE);

  //Don't save out graphics
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_GRAPHICS_OPTS, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_GRAPHICS_OPTS, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_GRAPHICS_OPTS, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_GRAPHICS_OPTS, CUBIT_FALSE);

  //Flag to say "we're copying" attributes
  copyingEntity = ref_ent;

  //Get all children
  DLIList<RefEntity*> child_list;
  ref_ent->get_all_child_ref_entities( child_list );

  //Puts attributes on all the entities
  child_list.append( ref_ent );
  CubitAttribUser::auto_update_cubit_attrib( child_list );

  //Prepare virtual
  TopologyEntity *topo_ptr = dynamic_cast<TopologyEntity*>( ref_ent );
  top_bridge = topo_ptr->bridge_manager()->topology_bridge();
  DLIList<TopologyBridge*> bridge_list;
  bridge_list.append( top_bridge );
  GeometryQueryTool::instance()->ige_export_geom( bridge_list );

  //Should only have 1 bridge in the list
  assert( bridge_list.size() == 1 );
  TopologyBridge *tmp_bridge_after = bridge_list.get();
  if( top_bridge != tmp_bridge_after)
    top_bridge = tmp_bridge_after;

  //Done copying attributes
  copyingEntity = NULL;

  return CUBIT_SUCCESS;

}

CubitStatus GeometryModifyTool::clean_up_from_copy_failure( TopologyBridge *old_bridge )
{
  //Remove attributes on underlying entities of virtual geometry
  //entities which do not have corresponding RefEntities.
  DLIList<TopologyBridge*> bridge_list;
  bridge_list.append( old_bridge );
  GeometryQueryTool::instance()->ige_remove_attributes( bridge_list );

  //Remove attributes on original entity and children
  DLIList<RefEntity*> child_list;
  RefEntity *ref_ent = CAST_TO( old_bridge->topology_entity(), RefEntity );
  ref_ent->get_all_child_ref_entities( child_list );
  //child_list.append( ref_ent );
  CubitAttribUser::clear_all_simple_attrib( child_list );
  child_list.clean_out();
  child_list.append( ref_ent );
  CubitAttribUser::clear_all_simple_attrib( child_list );

  CGMApp::instance()->restore_previous_attribute_states();

return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::finish_copy( TopologyBridge *&new_bridge,
                                             TopologyBridge *old_bridge )
{
  //Remove attributes on underlying entities of virtual geometry
  //entities which do not have corresponding RefEntities.
  DLIList<TopologyBridge*> bridge_list;
  bridge_list.append( old_bridge );
  GeometryQueryTool::instance()->ige_remove_attributes( bridge_list );

  //Remove attributes on original entity and children
  DLIList<RefEntity*> child_list;
  RefEntity *ref_ent = CAST_TO( old_bridge->topology_entity(), RefEntity );
  ref_ent->get_all_child_ref_entities( child_list );
  //child_list.append( ref_ent );
  CubitAttribUser::clear_all_simple_attrib( child_list );
  child_list.clean_out();
  child_list.append( ref_ent );
  CubitAttribUser::clear_all_simple_attrib( child_list );


  //Restore virtual
  //Could create virtual geometry here so the Topology Bridge can change
  bridge_list.clean_out();
  bridge_list.append( new_bridge );
  TopologyBridge *tmp_bridge_before = new_bridge;
  GeometryQueryTool::instance()->ige_import_geom( bridge_list );
  assert( bridge_list.size() == 1 );
  if( tmp_bridge_before != bridge_list.get() )
    new_bridge = bridge_list.get();

  //make the RefEntities
  Curve *curve = NULL;
  Surface *surface = NULL;
  Lump *lump = NULL;
  BodySM *body = NULL;

  RefEntity *new_ref_ent = NULL;
  if( (curve = CAST_TO( new_bridge, Curve ) ) != NULL )
    new_ref_ent = GeometryQueryTool::instance()->make_RefEdge( curve );
  else if( (surface = CAST_TO( new_bridge, Surface) ) != NULL )
    new_ref_ent = GeometryQueryTool::instance()->make_RefFace( surface );
  else if( (lump = CAST_TO( new_bridge, Lump) ) != NULL )
    new_ref_ent = GeometryQueryTool::instance()->make_Body( lump->bodysm() );
  else if( (body = CAST_TO( new_bridge, BodySM) ) != NULL )
    new_ref_ent = GeometryQueryTool::instance()->make_Body( body );

  //actuate the attributes on everything
  child_list.clean_out();
  new_ref_ent->get_all_child_ref_entities( child_list );
  child_list.append( new_ref_ent );
  GeometryQueryTool::import_actuate( child_list );

  //Remove attributes on new entity and children
  child_list.clean_out();
  new_ref_ent->get_all_child_ref_entities( child_list );
  CubitAttribUser::clear_all_simple_attrib( child_list );
  child_list.clean_out();
  child_list.append( new_ref_ent );
  CubitAttribUser::clear_all_simple_attrib( child_list );

  CGMApp::instance()->restore_previous_attribute_states();
  return CUBIT_SUCCESS;
}
#ifdef PROE
CubitStatus GeometryModifyTool::prepare_for_topology_update( BodySM* old_bodysm )
{
  DLIList<BodySM*> old_bodysms;
  old_bodysms.append(old_bodysm);

  do_attribute_setup();
  push_vg_attributes_before_modify( old_bodysms );

  return CUBIT_SUCCESS;
}
CubitStatus GeometryModifyTool::finish_topology_update( BodySM* new_bodysm,
                                                        Body* old_body )
{
  DLIList<Body*> input_bodies,result_bodies;
  DLIList<BodySM*> new_bodysms, old_bodysms;
  input_bodies.append(old_body);
  new_bodysms.append(new_bodysm);
  old_bodysms.clean_out();

  DLIList<int> merged_surface_ids;
  DLIList<int> merged_curve_ids;


  get_merged_curve_and_surface_ids( input_bodies, merged_surface_ids, merged_curve_ids );
  ///*
  //Store information about what surfaces were originally merged
  //
  //This fixes a problem caused by the lack of support
  //for automatically remerging surfaces that contain any
  //virtual geometry.
  //This should be removed if this issue is ever fixed - AHH
  std::map< int , DLIList<Surface*> > merged_map;
  for(int i=0; i< merged_surface_ids.size(); i++ )
    {
    int refface_id = merged_surface_ids.get_and_step();
    RefFace* old_merged_refface = RefEntityFactory::instance()->get_ref_face( refface_id );
    if( old_merged_refface && old_merged_refface->bridge_manager()->number_of_bridges() > 1 )
      {
      DLIList<Surface*> merged_surfsms;
      DLIList<TopologyBridge*> bridge_list;
      old_merged_refface->bridge_manager()->get_bridge_list( bridge_list );
      for(int j=0; j< bridge_list.size(); j++ )
        {
        Surface* merging_surf = CAST_TO(bridge_list.get_and_step(),Surface);
        if( merging_surf )
          merged_surfsms.append_unique( merging_surf );
        }
      merged_map[ refface_id ] = merged_surfsms;
      }
    }
  //*/
  /*
  //check to see if any curves have been orphaned inside virtual geometry
  DLIList<Point*> point_list;
  new_bodysm->points( point_list );
  CubitBoolean loop = CUBIT_TRUE;
  while( loop )
    {
    for(i=0; i< point_list.size(); i++)
      {
      //loop all the points and see if any will be on only one live curve
      Point* curr_point = point_list.get_and_step();
      DLIList<Curve*> curves_on_point;
      DLIList<Curve*> curves_to_keep;
      curr_point->curves( curves_on_point );
      for(int j=0; j< curves_on_point.size(); j++)
        {
        Curve* curr_curve = curves_on_point.get_and_step();
        DLIList<CubitSimpleAttrib*> attrib_list;
        CubitSimpleAttrib* attrib = CompositeEngine::find_attribute_by_name( curr_curve, "COMPOSITE_GEOM" );
        if( attrib )
          {
          attrib_list.append(attrib);
          }
        else
          curves_to_keep.append_unique( curr_curve );
        //curr_curve->get_simple_attribute( "COMPOSITE_GEOM", attrib_list );
        //if ( !attrib_list.size() )
          //curves_to_keep.append_unique( curr_curve );
        }
      //if all but one curve on this point will be removed we need to remove the other one too
      if( curves_to_keep.size() == 1 )
        {
        CubitString name("COMPOSITE_GEOM");
        DLIList<CubitString*> string_list;
        string_list.append( &name );
        CubitSimpleAttrib geom_attrib( &string_list, 0, 0 );
        curves_to_keep.get()->append_simple_attribute_virt( &geom_attrib );
        curr_point->append_simple_attribute_virt( &geom_attrib );
        loop = CUBIT_FALSE;
        }
      else if( curves_to_keep.size() == 0 )
        {
        DLIList<CubitSimpleAttrib*> attrib_list;
        curr_point->get_simple_attribute( "COMPOSIT_GEOM", attrib_list );
        if( !attrib_list.size() )
          {
          CubitString name("COMPOSITE_GEOM");
          DLIList<CubitString*> string_list;
          string_list.append( &name );
          CubitSimpleAttrib geom_attrib( &string_list, 0, 0 );
          curr_point->append_simple_attribute_virt( &geom_attrib );
          loop = CUBIT_TRUE;
          }
        }
      }
    if( loop )
      loop = CUBIT_FALSE;
    else
      loop = CUBIT_TRUE;
    }
  */
  GeometryModifyEngine* gme = GeometryModifyTool::instance()->get_engine(old_body);

  restore_vg_after_modify(new_bodysms, input_bodies, gme);


  DLIList<RefVolume*> volume_list;

  int i;
  for(i=new_bodysms.size(); i--;)
  {
    BodySM *bsm = new_bodysms.get_and_step();
    Body *body = dynamic_cast<Body*>(bsm->topology_entity());
    if(body)
    {
      // Append to the total list of volumes.
      body->ref_volumes(volume_list);
    }
  }
  // get all child entities (only get entities below volumes)
  DLIList<RefEntity*> child_list, ref_ent_list;
  CAST_LIST_TO_PARENT(volume_list, ref_ent_list);
  RefEntity::get_all_child_ref_entities( ref_ent_list, child_list );

  // Only push the id attributes if we are doing persistent ids.
  for(i=child_list.size(); i--;)
  {
  child_list.get_and_step()->auto_actuate_cubit_attrib(CUBIT_FALSE,CUBIT_TRUE);
  }

  remove_pushed_attributes(new_bodysms, input_bodies);

  finish_sm_op( input_bodies, new_bodysms , result_bodies);

  fixup_merged_entities( merged_surface_ids, merged_curve_ids);

  //Look for RefEdges that are orphaned inside of virtual surfaces
  CubitBoolean loop = CUBIT_TRUE;
  DLIList<RefVertex*> vertex_list;
  while(loop)
    {
    loop = CUBIT_FALSE;
    vertex_list.clean_out();
    for(i=0; i< volume_list.size(); i++)
      {
      RefVolume* curr_volume = volume_list.get_and_step();
      curr_volume->ref_vertices( vertex_list );
      }
    for(i=0; i< vertex_list.size(); i++)
      {
      RefVertex* curr_vertex = vertex_list.get_and_step();
      DLIList<RefEdge*> edges_on_vertex;
      curr_vertex->ref_edges( edges_on_vertex );
      edges_on_vertex.uniquify_unordered();
      if( edges_on_vertex.size() == 1 )
        {
        RefEdge* edge_to_remove = edges_on_vertex.get();
        CompositeTool::instance()->remove_edge(edge_to_remove);
        loop = CUBIT_TRUE;
        break;
        }
      }
    }
  ///*
  //Force merge faces that were missed by finish_sm_op
  //
  //This fixes a problem caused by the lack of support
  //for automatically remerging surfaces that contain any
  //virtual geometry.
  //This should be removed if this issue is ever fixed - AHH
  for(i=0; i< merged_surface_ids.size(); i++)
    {
    int refface_id = merged_surface_ids.get_and_step();
    DLIList<RefFace*> ref_faces_to_remerge;
    if(merged_map.find( refface_id ) != merged_map.end())
      {
      DLIList<Surface*> merged_surfsms = merged_map[ refface_id ];
      for(int k=0; k< merged_surfsms.size(); k++)
        {
        Surface* merging_surf = merged_surfsms.get_and_step();
        RefFace* merging_refface = CAST_TO(merging_surf->topology_entity(),RefFace);
        if(merging_refface)
          ref_faces_to_remerge.append_unique(merging_refface);
        }
      }
    if( ref_faces_to_remerge.size() == 2 )
      {
      RefFace* first_remerge_face = ref_faces_to_remerge.get_and_step();
      RefFace* second_remerge_face = ref_faces_to_remerge.get_and_step();
      int first_id = first_remerge_face->id();
      int second_id = second_remerge_face->id();
      MergeTool::instance()->force_merge( first_remerge_face , second_remerge_face );
      PRINT_WARNING("Remerging Surfaces %i and %i\n",first_id,second_id);
      }
    }
  //*/
  //do_attribute_cleanup();

  return CUBIT_SUCCESS;
}

#endif
GeometryModifyEngine*
GeometryModifyTool::make_RefEdge_common( RefVertex* start_vertex,
                                         RefVertex* end_vertex,
                                         Point*& start_point,
                                         Point*& end_point,
                                         RefFace* ref_face,
                                         Surface** surface ) const
{
  DLIList<TopologyEntity*> entity_list(3);
  DLIList<TopologyBridge*> bridge_list(3);
  TopologyBridge* bridge = 0;
  GeometryModifyEngine* gme = 0;
  start_point = end_point = 0;

  bool new_start_point = start_vertex->get_parents() > 0;
  bool new_end_point = end_vertex->get_parents() > 0;

  if (ref_face)
    entity_list.append( ref_face );
  if (!new_start_point)
    entity_list.append( start_vertex );
  if (!new_end_point)
    entity_list.append( end_vertex );

  if (entity_list.size())
    gme = common_modify_engine( entity_list, bridge_list );

  if (gme)
  {
    bridge_list.reset();
    if (ref_face)
      *surface = dynamic_cast<Surface*>(bridge_list.get_and_step());
    if (!new_start_point)
      start_point = dynamic_cast<Point*>(bridge_list.get_and_step());
    if (!new_end_point)
      end_point = dynamic_cast<Point*>(bridge_list.get_and_step());
  }
  else if (ref_face)
  {
    gme = get_engine( ref_face, &bridge );
    if (!gme)
    {
      PRINT_ERROR("No modify engine for surface %d\n", ref_face->id()) ;
      return 0;
    }
    *surface = dynamic_cast<Surface*>(bridge);
    GeometryQueryEngine* gqe = bridge->get_geometry_query_engine();
    if (!new_start_point)
      start_point = dynamic_cast<Point*>( start_vertex->
                                  bridge_manager()->topology_bridge( gqe ) );
    if (!new_end_point)
      end_point = dynamic_cast<Point*>( end_vertex->
                                  bridge_manager()->topology_bridge( gqe ) );
  }
  else if (!new_start_point && (gme = get_engine( start_vertex, &bridge )))
  {
    start_point = dynamic_cast<Point*>(bridge);
    if (!new_end_point)
      end_point = dynamic_cast<Point*>( end_vertex->bridge_manager()->
                    topology_bridge( bridge->get_geometry_query_engine() ) );
  }
  else if (!new_end_point && (gme = get_engine( end_vertex, &bridge )))
  {
    end_point = dynamic_cast<Point*>(bridge);
  }
  else
  {
    gme = get_gme();
  }


  if (!start_point)
    start_point = gme->make_Point( start_vertex->coordinates() );
  if (!end_point)
    end_point = gme->make_Point( end_vertex->coordinates());

  return gme;
}

//-------------------------------------------------------------------------
// Purpose       : This function uses the vertices and surface to create
//                 a curve along this surface.  The solid modeler will actually
//                 do most of the work.
// Note          : The optional third vertex will not be part of the curve,
//                 it is used for interpolation purposes.
//
// Creator       : David White
//
// Creation Date : 10/9/97
//
// Rewriting for new GeometryQueryTool interface - J.Kraftcheck 9/03
//-------------------------------------------------------------------------
RefEdge* GeometryModifyTool::make_RefEdge(RefVertex *ref_vertex_1,
                                    RefVertex *ref_vertex_2,
                                    RefFace* ref_face_ptr,
                                    RefVertex const* ref_vertex_3 ) const
{
	  //make sure that the vertices are on the ref-face.
   CubitVector vert_1 = ref_vertex_1->coordinates();
   CubitVector vert_2 = ref_vertex_2->coordinates();

   ref_face_ptr->move_to_surface(vert_1);
   ref_face_ptr->move_to_surface(vert_2);

   GeometryQueryEngine *gqe = ref_face_ptr->get_geometry_query_engine();

   if (!ref_vertex_1->coordinates().within_tolerance(vert_1, gqe->get_sme_resabs_tolerance()))
   {
      PRINT_ERROR("vertices must lie within tolerance to the given"
                  " surface.\n"
                  "%s (Vertex %d) does not lie on %s (Surface %d).\n",
                  ref_vertex_1->entity_name().c_str(),
                  ref_vertex_1->id(),
                  ref_face_ptr->entity_name().c_str(),
                  ref_face_ptr->id() );
      return (RefEdge*)NULL;
   }
   else if (!ref_vertex_2->coordinates().within_tolerance(vert_2, gqe->get_sme_resabs_tolerance() ))
   {
      PRINT_ERROR("vertices must lie within tolerance to the given"
                  " surface.\n"
                  "%s (Vertex %d) does not lie on %s (Surface %d).\n",
                  ref_vertex_2->entity_name().c_str(),
                  ref_vertex_2->id(),
                  ref_face_ptr->entity_name().c_str(),
                  ref_face_ptr->id() );
      return (RefEdge*)NULL;
   }
     //Now let us find the points on the surface.  We want to
     //create them as we go for accuracy.

	  // Get the GME of the first RefVertex.
   GeometryModifyEngine* GMEPtr = 0;

     // Extract the end Points to be used to make the RefEdge
   Point* point_ptr1 = NULL;
   Point* point_ptr2 = NULL;
   Surface* surface_ptr = NULL;

   // Look for a common GeometryModifyEngine
  GMEPtr = make_RefEdge_common( ref_vertex_1, ref_vertex_2,
                                point_ptr1, point_ptr2,
                                ref_face_ptr, &surface_ptr );

  //If we did not find a common GeometryModifyEngine, fail
  if( ! GMEPtr )
  {
    PRINT_ERROR("Surface %d, vertex %d and vertex %d do not "
      "belong to the same geometric modeling engine.\n", ref_face_ptr->id(),
      ref_vertex_1->id(), ref_vertex_2->id() );

    return 0;
  }

   CubitVector *third_vector_ptr = NULL;
   CubitVector third_vector;
   if ( ref_vertex_3 != NULL )
   {
      third_vector = ref_vertex_3->coordinates();
      third_vector_ptr = &third_vector;
      ref_face_ptr->move_to_surface( third_vector );
   }

     // Make sure that we get back valid Points
   assert ( point_ptr1 != NULL && point_ptr2 != NULL );




     // Request the GME to create a Curve using the Points
   Curve *curve_ptr;
   curve_ptr = GMEPtr->make_Curve( point_ptr1,
                          point_ptr2,
                          surface_ptr,
                          third_vector_ptr);


     // If we get a NULL pointer, give a warning message and return
     // a NULL pointer.
   if ( curve_ptr == NULL )
   {
      PRINT_WARNING("In GeometryModifyTool::make_RefEdge\n"
                    "\tProblems making a spline curve from Vertex %d and %d\n"
                    "\tand the input list of positions.\n",
                    ref_vertex_1->id(), ref_vertex_2->id());
      return (RefEdge *)NULL ;
   }

   if( CubitUndo::get_undo_enabled() )
   {
     //if endpoints are free vertices, need to save them out
     DLIList<RefVertex*> verts_to_save;
     verts_to_save.append( ref_vertex_1 );
     verts_to_save.append( ref_vertex_2 );
     bool save_only_if_free = true;
     CubitUndo::save_state_with_cubit_file( verts_to_save, save_only_if_free );
   }

     // Complete the task of linking this new Curve into the rest of the
     // geometry datastructures and return the new RefEdge.
   RefEdge *new_edge = GeometryQueryTool::instance()->make_free_RefEdge(curve_ptr);
   if( CubitUndo::get_undo_enabled() )
   {
     if( new_edge )
       CubitUndo::note_result_entity( new_edge );
   }
   return new_edge;
}

//-------------------------------------------------------------------------
// Purpose       : This function takes RefEdge type information, two
//                 RefVertices, and a list of positions in space (represented
//                 by CubitVectors) to create a RefEdge. The RefVertices must
//                 be associated with the same GeometryModifyEngine.
//                 The return value can be a NULL pointer, if the RefEdge
//                 cannot be succesfully made for some reason.
//
// Special Notes : The interpolated curve is a spline curve.
//                 If the input refface_ptr is not NULL, the points are first
//                 moved to the surface before interpolation.
//
// Creator       : Malcolm Panthaki
//
// Creation Date : 05/15/97
//
// Rewriting for new GeometryQueryTool interface - J.Kraftcheck 9/03
//-------------------------------------------------------------------------
RefEdge* GeometryModifyTool::make_RefEdge(GeometryType ref_edge_type,
                                    RefVertex *ref_vertex_1,
                                    RefVertex *ref_vertex_2,
                                    DLIList<CubitVector*>& vector_list,
                                    RefFace* refface_ptr) const
{
   GeometryModifyEngine* GMEPtr = 0;

     // Extract the end Points to be used to make the RefEdge
   Point* point_ptr1 = NULL;
   Point* point_ptr2 = NULL;
   Surface* surface_ptr = NULL;

   // Look for a common GeometryModifyEngine
  GMEPtr = make_RefEdge_common( ref_vertex_1, ref_vertex_2,
                                point_ptr1, point_ptr2,
                                refface_ptr, &surface_ptr );
  if (!GMEPtr)
    return 0;

     // Make sure that we get back valid Points
   assert ( point_ptr1 != NULL && point_ptr2 != NULL );

   if( CubitUndo::get_undo_enabled() )
   {
     //if endpoints are free vertices, need to save them out
     DLIList<RefVertex*> verts_to_save;
     verts_to_save.append( ref_vertex_1 );
     verts_to_save.append( ref_vertex_2 );
     bool save_only_if_free = true;
     CubitUndo::save_state_with_cubit_file( verts_to_save, save_only_if_free );
   }

     // Request the GME to create a Curve using the Points
   Curve* curve_ptr = GMEPtr->make_Curve(ref_edge_type,
                                         point_ptr1,
                                         point_ptr2,
                                         vector_list,
                                         refface_ptr ? refface_ptr->get_surface_ptr() : 0);

     // If we get a NULL pointer, give a warning message and return
     // a NULL pointer.
   if ( curve_ptr == NULL )
   {
      PRINT_WARNING("In GeometryModifyTool::make_RefEdge\n"
                    "         Got a NULL pointer to a Curve.\n"
                    "         Problems making a spline RefEdge from RefVertex %d and %d\n"
                    "         and the input list of positions.\n",
                    ref_vertex_1->id(), ref_vertex_2->id());
      return (RefEdge *)NULL ;
   }

     // Complete the task of linking this new Curve into the rest of the
     // geometry datastructures and return the new RefEdge.
   RefEdge *new_edge = GeometryQueryTool::instance()->make_free_RefEdge(curve_ptr);
   if( CubitUndo::get_undo_enabled() )
   {
     if( new_edge )
       CubitUndo::note_result_entity( new_edge );
   }
   return new_edge;
}

//-------------------------------------------------------------------------
// Purpose       : This function takes a type information, two
//                 RefVertices, an intermediate position, and a sense
//                 information to create a RefEdge. The RefVertices must
//                 be associated with the same GeometryModifyEngine.
//                 The return value can be a NULL pointer, if the RefEdge
//                 cannot be succesfully made for some reason.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/27/97
//-------------------------------------------------------------------------
RefEdge* GeometryModifyTool::make_RefEdge(GeometryType ref_edge_type,
                                    RefVertex *ref_vertex_1,
                                    RefVertex *ref_vertex_2,
                                    CubitVector const* intermediate_point,
                                    CubitSense sense) const
{
   GeometryModifyEngine* GMEPtr = 0;

     // Extract the Points to be used to make the RefEdge
   Point* point_ptr1 = NULL;
   Point* point_ptr2 = NULL;

    // Look for a common geometric modeling engine to use
  GMEPtr = make_RefEdge_common( ref_vertex_1, ref_vertex_2,
                                point_ptr1, point_ptr2 );
  if (!GMEPtr)
    return 0;

     // Make sure that we get back valid Points
   assert ( point_ptr1 != NULL && point_ptr2 != NULL ) ;

   if( CubitUndo::get_undo_enabled() )
   {
     //if endpoints are free vertices, need to save them out
     DLIList<RefVertex*> verts_to_save;
     verts_to_save.append( ref_vertex_1 );
     verts_to_save.append( ref_vertex_2 );
     bool save_only_if_free = true;
     CubitUndo::save_state_with_cubit_file( verts_to_save, save_only_if_free );
   }

     // Request the GME to create a Curve using the Points
   Curve* curve_ptr = GMEPtr->make_Curve(ref_edge_type,
                                         point_ptr1,
                                         point_ptr2,
                                         intermediate_point,
                                         sense) ;

     // If we get a NULL pointer, give a warning message and return
     // a NULL pointer.
   if ( curve_ptr == NULL )
   {
      PRINT_WARNING("In GeometryModifyTool::make_RefEdge\n"
                    "         Got a NULL pointer to a Curve.\n"
                    "         Problems making RefEdge from RefVertex %d and %d\n",
                    ref_vertex_1->id(), ref_vertex_2->id());
      return (RefEdge *)NULL ;
   }

     // Complete the task of linking this new Curve into the rest of the
     // geometry datastructures and return the new RefEdge.
   RefEdge *new_edge = GeometryQueryTool::instance()->make_free_RefEdge( curve_ptr );
   if( CubitUndo::get_undo_enabled() )
   {
     if( new_edge )
       CubitUndo::note_result_entity( new_edge );
   }
   return new_edge;
}

//-------------------------------------------------------------------------
// Purpose       : This function creates a RefFace from an existing RefFace.
//                 The new face is a sheet body, double sided face with no
//                 volume to the body.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/27/97
//-------------------------------------------------------------------------

RefFace* GeometryModifyTool::make_RefFace(RefFace *from_ref_face,
                                    CubitBoolean extended_from) const
{
   GeometryModifyEngine* GME_ptr = get_engine(from_ref_face);

   if (!GME_ptr)
   {
      PRINT_ERROR("Cannot create surface from another virtual surface.\n");
      return NULL;
   }
     //From the surface create a new surface.
   TopologyBridge* bridge = 0;
   GeometryModifyEngine* engine = get_engine(from_ref_face, &bridge);
   Surface* old_surface_ptr = dynamic_cast<Surface*>(bridge);
   if ( engine == NULL )
   {
      PRINT_ERROR("%s (surface %d) does not have a modify engine.\n",
                  from_ref_face->entity_name().c_str(),
                  from_ref_face->id() );
      return (RefFace*)NULL;
   }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state();

  //this list will get all the TB's what we'll be copying
  TopologyBridge *top_bridge = old_surface_ptr;
  if( !extended_from )
  {
    prepare_for_copy( from_ref_face, top_bridge );
  }

  Surface *tmp_surface = CAST_TO( top_bridge, Surface );
  Surface* new_surface_ptr = engine->make_Surface( tmp_surface,
                                                    extended_from );

  if (!new_surface_ptr)
  {
    PRINT_ERROR("Surface copy failed.\n");
    if( !extended_from )
      clean_up_from_copy_failure( old_surface_ptr );

   if( CubitUndo::get_undo_enabled() )
        CubitUndo::remove_last_undo();

    return 0;
  }

  if( !extended_from  )
  {
    TopologyBridge *top_bridge_new = new_surface_ptr;
    finish_copy( top_bridge_new, old_surface_ptr );
    new_surface_ptr = CAST_TO( top_bridge_new, Surface );
  }

  Body *new_Body = make_Body(new_surface_ptr);
  DLIList<RefFace*> ref_faces;
  new_Body->ref_faces(ref_faces);
  assert(ref_faces.size() > 0);

  if( CubitUndo::get_undo_enabled() )
  {
    if( new_Body )
      CubitUndo::note_result_body( new_Body );
    else
      CubitUndo::remove_last_undo();
  }

  return ref_faces.get();
}

//-------------------------------------------------------------------------
// Purpose       : This function creates a sheet body by extending out a
//                 set of surfaces.  The sheet body is a double sided face
//                 with no volume.
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 02/28/08
//-------------------------------------------------------------------------
Body*
GeometryModifyTool::make_extended_sheet( DLIList<RefFace*> &ref_face_list,
                                         CubitBox *clip_box_ptr,
                                         bool preview ) const
{
  if( !ref_face_list.size() )
    return 0;

  // Check for virtual geometry
  DLIList<RefEntity*> ref_ent_list;
  CAST_LIST_TO_PARENT(ref_face_list, ref_ent_list);
  if ( GeometryQueryTool::instance()->contains_intermediate_geometry(ref_ent_list) )
  {
    PRINT_ERROR("EXTENDING surfaces containing virtual geometry is not\n"
      "       allowed. Delete virtual geometry on these surfaces\n"
      "       before operation.\n" );
    return 0;
  }

  // Look for a common GeometryModifyEngine for all of the RefFaces
  int count = ref_face_list.size();
  DLIList<TopologyBridge*> bridge_list(count);
  DLIList<TopologyEntity*> entity_list(count);
  CAST_LIST_TO_PARENT( ref_face_list, entity_list );

  GeometryModifyEngine* GME_ptr =
    common_modify_engine( entity_list, bridge_list );
  if(! GME_ptr )
  {
     PRINT_ERROR("Cannot construct an extended sheet using surfaces that\n"
                 "       do not share a common geometry engine.\n");
     return 0;
  }

  DLIList<Surface*> surface_list(count);
  CAST_LIST( bridge_list, surface_list, Surface );

  BodySM *bodySM_ptr = GME_ptr->make_extended_sheet( surface_list, 
    clip_box_ptr, preview );

  if( bodySM_ptr )
    return GeometryQueryTool::instance()->make_Body(bodySM_ptr);
  else
    return 0;
}

//-------------------------------------------------------------------------
// Purpose       : This function takes a type information and a list of
//                 RefEdges to create a Body with just one RefFace. The
//                 underlying representation of the RefFace is determined
//                 by the GeometryModifyEngine of the RefEdges. All the
//                 RefEdges in the list must be associated with the same
//                 GeometryModifyEngine. The return value can be a
//                 NULL pointer, if the RefFace cannot be succesfully
//                 made for some reason.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 10/23/97
//-------------------------------------------------------------------------

RefFace* GeometryModifyTool::make_RefFace(GeometryType ref_face_type,
                                    DLIList<RefEdge*>& ref_edge_list,
                                    bool is_free_face,
                                    RefFace *ref_face_ptr,
                                    bool check_edges ) const
{
    //Look for a common GeometryModifyEngine for all of
    //the RefEdges.

  const int count = ref_edge_list.size() + (ref_face_ptr ? 1 : 0);
  DLIList<TopologyBridge*> bridge_list(count);
  DLIList<TopologyEntity*> entity_list(count);
  CAST_LIST_TO_PARENT( ref_edge_list, entity_list );
  if( ref_face_ptr ) entity_list.append( ref_face_ptr );

  GeometryModifyEngine* GME_ptr =
    common_modify_engine( entity_list, bridge_list );
  if(! GME_ptr )
  {
     PRINT_ERROR("Cannot construct a Surface using entities that do "
                 "not share a common geometry engine.\n");
     return 0;
  }

   Surface* old_surface_ptr = 0;
   if (ref_face_ptr)
    old_surface_ptr = dynamic_cast<Surface*>(bridge_list.pop());

   //Collect all the names on vertices to propagate after you create
   //the surface
   DLIList<CubitVector> vertex_coordinates;
   DLIList<CubitString*> vertex_names;
   DLIList<RefEdge*> free_ref_edges;
   int kk;
   for( kk=ref_edge_list.size(); kk--; )
   {
     DLIList<CubitString*> tmp_names;
     RefEdge *tmp_edge = ref_edge_list.get_and_step();

     if( tmp_edge->num_parent_ref_entities() == 0 )
       free_ref_edges.append( tmp_edge );

     RefVertex *s_vertex = tmp_edge->start_vertex();
     RefVertex *e_vertex = tmp_edge->end_vertex();
     int jj;

     s_vertex->entity_names( tmp_names );
     for( jj=tmp_names.size(); jj--; )
     {
       CubitVector tmp_vec = tmp_edge->start_vertex()->coordinates();
       CubitString *name = tmp_names.get_and_step();
       vertex_coordinates.append( tmp_vec );
       vertex_names.append( new CubitString(*name) );
     }

     tmp_names.clean_out();
     e_vertex->entity_names( tmp_names );
     for( jj=tmp_names.size(); jj--; )
     {
       CubitVector tmp_vec = tmp_edge->end_vertex()->coordinates();
       CubitString *name = tmp_names.get_and_step();
       vertex_coordinates.append( tmp_vec );
       vertex_names.append( new CubitString(*name) );
     }
   }

   DLIList<Curve*> curve_list(ref_edge_list.size());
   CAST_LIST( bridge_list, curve_list, Curve );

     // Use the Curves to create a Surface
   Surface* surface_ptr = GME_ptr->make_Surface(ref_face_type, curve_list,
                                                old_surface_ptr, check_edges) ;

   if (surface_ptr == NULL) {
     PRINT_ERROR("Couldn't make new RefFace.\n");
     return NULL;
   }

   GeometryQueryTool* const gqt = GeometryQueryTool::instance();

   RefFace* result_face = gqt->make_free_RefFace(surface_ptr, is_free_face);
   gqt->cleanout_deactivated_geometry();

   //send out events for free curves saying that their 'free' status has
   //be changed
  
  for( kk=0; kk<free_ref_edges.size(); kk++ )
  {
    RefEdge *free_edge = free_ref_edges.get_and_step();
    CubitObserver::notify_static_observers( free_edge, TOP_LEVEL_ENTITY_DESTRUCTED );
    CGMHistory::Event evt(CGMHistory::TOP_LEVEL_ENTITY_DELETED, free_edge );
    GeometryQueryTool::instance()->history().add_event(evt);
  }


   //look for a vertex at the same location of the original
   //vertex(s).  Add the name to this new vertex.
   DLIList<RefVertex*> tmp_verts;
   result_face->ref_vertices( tmp_verts);
   for( kk=vertex_coordinates.size(); kk--; )
   {
     CubitVector tmp_coord = vertex_coordinates.get_and_step();
     CubitString *tmp_name = vertex_names.get_and_step();

     int jj;
     for( jj=tmp_verts.size(); jj--; )
     {
       RefVertex *tmp_vert = tmp_verts.get_and_step();
       if( tmp_coord.distance_between( tmp_vert->coordinates() ) < GEOMETRY_RESABS )
       {
         //add the name if it doesn't already exist
         RefEntityName::instance()->add_refentity_name( tmp_vert, *tmp_name );
       }
     }
   }

   return result_face;
}

//-------------------------------------------------------------------------
// Purpose       : This function takes a list of RefVolumes to create a
//                 Body. The underlying representation of the Body is
//                 determined by the GeometryModifyEngine of the
//                 RefVolumes. All the RefVolumes in the list must be
//                 associated with the same GeometryModifyEngine. The
//                 return value can be a NULL pointer, if the RefFace
//                 cannot be succesfully made for some reason.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------

Body* GeometryModifyTool::make_Body(DLIList<RefVolume*>& ref_volume_list) const
{
  DLIList<TopologyEntity*> topo_list;
  DLIList<Lump*> lump_list;
  DLIList<TopologyBridge*> bridge_list;
  CAST_LIST_TO_PARENT( ref_volume_list, topo_list );
  GeometryModifyEngine* GME_ptr =
    common_modify_engine( topo_list, bridge_list, CUBIT_FALSE );
  if( ! GME_ptr )
  {
    PRINT_ERROR("The specified volumes do not belong to the same "
      "geometry engine, and therefore cannot be used to form a "
      "single body.");
    return 0;
  }
  CAST_LIST( bridge_list, lump_list, Lump );
  assert( ref_volume_list.size() == lump_list.size() );


     // Use the Lumps to create a BodySM
   BodySM* bodySM_ptr = GME_ptr->make_BodySM(lump_list) ;

     // If we get a NULL pointer, give a warning message and return
     // a NULL pointer.
   if ( bodySM_ptr == NULL )
   {
      PRINT_WARNING("In GeometryModifyTool::make_Body\n"
                    "         Got a NULL pointer to a BodySM.\n"
                    "         Cannot make a Body.\n") ;
      return (Body *)NULL;
   }

   return GeometryQueryTool::instance()->make_Body(bodySM_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : This function creates a sheet body from an existing
//                 RefFace.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------

Body* GeometryModifyTool::make_Body(RefFace *from_ref_face,
                              CubitBoolean extended_from) const
{
     // Given the arguments, make a RefFace.
   RefFace* new_ref_face = this->make_RefFace(from_ref_face,
                                              extended_from );
   if (!new_ref_face)
   {
      return (Body *)NULL;
   }

   DLIList<Body*> bodies;
   new_ref_face->bodies(bodies);

   if (!bodies.size()) {
     GeometryEntity *ge_ptr = new_ref_face->get_geometry_entity_ptr();
     Surface *surf_ptr = CAST_TO(ge_ptr, Surface);
     assert(surf_ptr != 0);

     BodySM *body_sm = gmeList.get()->make_BodySM(surf_ptr);
     return GeometryQueryTool::instance()->make_Body(body_sm);
   }

   else return bodies.get();
}
//-------------------------------------------------------------------------
// Purpose       : This function takes a type information and a list of
//                 RefEdges to create a Body that has just one RefFace.
//                 The underlying representation of the Body is
//                 determined by the GeometryModifyEngine of the
//                 RefEdges. All the RefEdges in the list must be
//                 associated with the same GeometryModifyEngine. The
//                 return value can be a NULL pointer, if the Body cannot
//                 be succesfully made for some reason.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------

Body* GeometryModifyTool::make_Body(GeometryType ref_face_type,
                              DLIList<RefEdge*>& ref_edge_list,
                              RefFace *ref_face_ptr) const
{
   if( CubitUndo::get_undo_enabled() )
     CubitUndo::save_state_with_cubit_file( ref_edge_list );
  
   bool is_free_face = false;

     // Given the arguments, make a RefFace.
   RefFace* new_ref_face = this->make_RefFace(ref_face_type,
                                              ref_edge_list,
                                              is_free_face,
                                              ref_face_ptr);

   if( new_ref_face == NULL )
      return NULL;

      // If new_ref_face doesn't have a body, create one.
   DLIList<Body*> bodies;
   new_ref_face->bodies(bodies);

   if (!bodies.size()) {
     Surface *surf_ptr = new_ref_face->get_surface_ptr();
     assert(surf_ptr != 0);
     GeometryModifyEngine* engine = get_engine(surf_ptr);
     assert(engine != 0);

     BodySM *body_sm = engine->make_BodySM(surf_ptr);
     Body* body = GeometryQueryTool::instance()->make_Body(body_sm);
     bodies.append( body );
   }

   if( CubitUndo::get_undo_enabled() )
   {
     if( bodies.size() )
       CubitUndo::note_result_entity( bodies.get() ); 
     else
      CubitUndo::remove_last_undo();
   }

     // Return the body containing the new RefFace.
   return bodies.get();
}

//-------------------------------------------------------------------------
// Purpose       : The following functions copy the input Body to create
//                 a new one.  Some transformations may be applied to the
//                 newly copied Body before returning it:
//                    copy
//                    copy_and_move
//                    copy_and_rotate
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/1/96
//-------------------------------------------------------------------------
Body* GeometryModifyTool::copy_body ( Body* bodyPtr )
{
   BodySM* body_sm = bodyPtr->get_body_sm_ptr();
   if (!body_sm)
   {
     PRINT_ERROR("Body %d is invalid -- no attached BodySM.\n", bodyPtr->id());
     return 0;
   }

   //this list will get all the TB's what we'll be copying
#ifdef BOYD17
   DLIList<TopologyBridge*> original_bridges;
#endif
   DLIList<RefEntity*> tmp_list;
   tmp_list.append( bodyPtr );
   TopologyBridge *top_bridge;
   prepare_for_copy( bodyPtr, top_bridge );

   BodySM *tmp_body_sm = CAST_TO( top_bridge, BodySM );

   GeometryModifyEngine* GMEPtr = get_engine(tmp_body_sm);
   if ( GMEPtr == NULL )
   {
     clean_up_from_copy_failure( top_bridge );
     PRINT_ERROR("Cannot copy volume %d\n"
                 "     Copy is supported for bodies based on solid "
                 "models only\n", bodyPtr->ref_volume()->id() );
     return NULL ;
   }

   BodySM* new_body = GMEPtr->copy_body(tmp_body_sm);

   if (!new_body)
   {
     clean_up_from_copy_failure( top_bridge );
     PRINT_ERROR("Failed to copy volume %d\n", bodyPtr->ref_volume()->id());
     return 0;
   }

   TopologyBridge *top_bridge_new = new_body;
   finish_copy( top_bridge_new, body_sm );

   new_body = CAST_TO( top_bridge_new, BodySM );
   return GeometryQueryTool::instance()->make_Body(new_body);
}

//-------------------------------------------------------------------------
// Purpose       : Check if bodies can be webcut
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/10/03
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::okay_to_modify(
                                 DLIList<Body*>& webcut_body_list,
                                 const char* op ) const
{
  CubitStatus ret = CUBIT_SUCCESS;

  Body *b = webcut_body_list.get();
  if(b)
  {
    GeometryModifyEngine *gme = get_engine(b);
    if(gme)
    {
      if(!gme->supports_interoperability() &&
          contains_intermediate_geom(webcut_body_list))
      {
        PRINT_ERROR("Intermixing real and virtual geometry operations using the current solid modeling kernel is not allowed.\n");
        ret = CUBIT_FAILURE;
      }
    }
  }

  if(ret == CUBIT_SUCCESS)
  {
    // If the operation is not one of the ones below...
    if(strcmp(op, "WEBCUT") &&
      strcmp(op, "CHOP") &&
      strcmp(op, "UNITE") &&
      strcmp(op, "TWEAK") &&
      strcmp(op, "IMPRINT") &&
      strcmp(op, "REGULARIZE") &&
      strcmp(op, "SPLIT_SURFACE") &&
      strcmp(op, "REMOVE_TOPOLOGY") &&
      strcmp(op, "SPLIT"))
    {
      if (contains_intermediate_geom(webcut_body_list))
      {
        PRINT_ERROR("Performing %s on volumes containing virtual geometry is not allowed.\n", op);
        ret = CUBIT_FAILURE;
      }
    }
    else
    {
      if(contains_partitions(webcut_body_list))
      {
        PRINT_ERROR("Performing %s on volumes containing virtual partitions is not allowed.\n", op);
        ret = CUBIT_FAILURE;
      }
    }
  }

  return ret;
}


//-------------------------------------------------------------------------
// Purpose       : Common code for finishing webcut operations.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/25/03
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::finish_webcut( DLIList<Body*>& webcut_body_list,
                                               DLIList<BodySM*>& result_body_sms,
                                               CubitBoolean merge,
                                               CubitStatus status,
                                               DLIList<Body*>& result_list,
                                               DLIList<int> *merged_surface_ids,
                                               DLIList<int> *merged_curve_ids,
                                               CubitBoolean print_info) const
{
   if (!finish_sm_op(webcut_body_list, result_body_sms, result_list,(bool)print_info))
     status = CUBIT_FAILURE;

   DLIList<Body*> temp_result_list;
   if (status)
   {
     status = separate_body_after_webcut( result_list, temp_result_list);
     result_list = temp_result_list;
   }

   if( merged_surface_ids && merged_curve_ids )
     fixup_merged_entities( *merged_surface_ids, *merged_curve_ids );

   if (merge && status)
   {
     DLIList<Body*> temp_results(result_list);
     status = MergeTool::instance()->merge_bodies( temp_results );
   }

   return status;
}


//-------------------------------------------------------------------------
// Purpose       : Complete solid modeling operation on bodies
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/25/03
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::finish_sm_op( DLIList<Body*>& input_bodies,
                                              DLIList<BodySM*>& new_bodies,
                                              DLIList<Body*>& result_bodies,
                                              bool print_info ) const
{
  int i;
  GeometryQueryTool* gqt = GeometryQueryTool::instance();

  DLIList<int> updated_ids, created_ids, destroyed_ids;
  DLIList<int> updated_vol_ids, created_vol_ids, destroyed_vol_ids;

  // traverse the body object and remove any meshes from modified objects
  int b;
  for (b = 0; b < input_bodies.size(); b++) {
	  Body* body = input_bodies.get_and_step();
	  body_premodify(body);
  }

    // Remove dead bodies
  input_bodies.reset();
  for (i = input_bodies.size(); i--;)
  {
    Body* body = input_bodies.step_and_get();

    BodySM* bodysm = body->get_body_sm_ptr();
    if (!bodysm)
    {
        // If the body was destroyed, update the list
      destroyed_ids.append( body->id() );

      DLIList<RefVolume*> temp_vols;
      body->ref_volumes( temp_vols );
      for( int nv = temp_vols.size(); nv > 0; nv-- )
      {
        if( !temp_vols.get()->get_lump_ptr() )
          destroyed_vol_ids.append( temp_vols.get_and_step()->id() );
      }

      input_bodies.change_to(0);
      gqt->destroy_dead_entity(body);
    }
    else
    {
      remove_dead_entity_names( body );
    }
  }
  gqt->cleanout_deactivated_geometry();

    // Create new bodies
  new_bodies.last();
  for (i = new_bodies.size(); i--; )
  {
    BodySM* bodysm = new_bodies.step_and_get();
    bool newbody = bodysm->owner() == 0;
    Body* body = gqt->make_Body(bodysm);
    result_bodies.append(body);

    if (newbody)
    {
      created_ids.append(body->id());

      DLIList<RefVolume*> temp_vols;
      body->ref_volumes( temp_vols );
      for( int nv = temp_vols.size(); nv > 0; nv-- )
      {
        created_vol_ids.append( temp_vols.get_and_step()->id() );
      }
    }
    else
    {
      updated_ids.append(body->id());

      DLIList<RefVolume*> temp_vols;
      body->ref_volumes( temp_vols );
      for( int nv = temp_vols.size(); nv > 0; nv-- )
      {
        updated_vol_ids.append( temp_vols.get_and_step()->id() );
      }
    }
  }
  gqt->cleanout_deactivated_geometry();

  if (print_info)
  {
    if( DEBUG_FLAG( 153 ) )
    {
      if (created_ids.size())
         CubitUtil::list_entity_ids( "Created body(s): ", created_ids );
      if (updated_ids.size())
         CubitUtil::list_entity_ids( "Updated body(s): ", updated_ids );
      if (destroyed_ids.size())
         CubitUtil::list_entity_ids( "Destroyed body(s): ", destroyed_ids );
    }

    if (created_vol_ids.size())
       CubitUtil::list_entity_ids( "Created volume(s): ", created_vol_ids );
    if (updated_vol_ids.size())
       CubitUtil::list_entity_ids( "Updated volume(s): ", updated_vol_ids );
    if (destroyed_vol_ids.size())
       CubitUtil::list_entity_ids( "Destroyed volume(s): ", destroyed_vol_ids );
  }

  return CUBIT_SUCCESS;
}



//-------------------------------------------------------------------------
// Purpose       : Webcut a body with a cylinder given the input parameters.
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 10/28/97
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::webcut_with_cylinder(
                                       DLIList<Body*>& webcut_body_list,
                                       double radius,
                                       const CubitVector &axis,
                                       const CubitVector &center,
                                       DLIList<Body*>& results_list,
                                       DLIList<Body*> &neighboring_bodies,
                                       ImprintType imprint_type,
                                       CubitBoolean merge,
                                       CubitBoolean preview)
{
  if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
    return CUBIT_FAILURE;

  CubitStatus rval = CUBIT_SUCCESS;

  if (preview)
  {
    // find the bounding box for the cylinder
    CubitBox bounding_box;
    Body* body_ptr = webcut_body_list.get_and_step();
    bounding_box = body_ptr->bounding_box();

    int i;
    for( i=1; i<webcut_body_list.size(); i++ )
    {
       body_ptr = webcut_body_list.get_and_step();
       bounding_box |= body_ptr->bounding_box();
    }

    int color = CUBIT_BLUE;
    GfxPreview::clear();
    GfxPreview::draw_cylinder(axis, center, bounding_box, (float) radius, color);
    GfxPreview::flush();
    return rval;
  }

  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<Body*> bodies_to_save;
    bodies_to_save += webcut_body_list;
    bodies_to_save += neighboring_bodies;
    CubitUndo::save_state_with_cubit_file( bodies_to_save );
  }

  const int count = webcut_body_list.size();
  DLIList<BodySM*> result_sm_list;
  DLIList<Body*> body_list(webcut_body_list);
  DLIList<BodySM*> engine_body_sms(count);
  DLIList<Body*> engine_bodies(count);
  GeometryModifyEngine* gme = 0;

  if(!preview)
    do_attribute_setup();

  while ( (gme = group_bodies_by_engine(body_list, engine_bodies, engine_body_sms)) )
  {
    //get all bodysms that we might modify in this operation
    DLIList<BodySM*> neighbor_imprint_list;
    DLIList<BodySM*> bodies_sm_to_modify;
    DLIList<Body*> bodies_to_modify;
    bodies_to_modify += engine_bodies;
    int i;
    for( i=neighboring_bodies.size(); i--; )
    {
      Body *tmp_body = neighboring_bodies.get_and_step();
      BodySM *tmp_body_sm = tmp_body->get_body_sm_ptr();
      if( gme == get_engine(tmp_body_sm ) )
      {
        neighbor_imprint_list.append( tmp_body_sm ); 
        bodies_to_modify.append( tmp_body );
      }
    }

    DLIList<int> merged_surface_ids;
    DLIList<int> merged_curve_ids;

    if(!preview)
    {
      bodies_sm_to_modify += neighbor_imprint_list;

      push_vg_attributes_before_modify( bodies_sm_to_modify );
      //get all the child entities that have been merged
      get_merged_curve_and_surface_ids( bodies_to_modify, merged_surface_ids, merged_curve_ids );
    }

    // note that preview actually gets handled before this point
    CubitStatus status = webcut_w_cylinder( engine_body_sms, radius, axis,
              center, neighbor_imprint_list, result_sm_list, imprint_type );

    if ( status != CUBIT_FAILURE )
    {
      if(!preview)
      {
        restore_vg_after_modify(result_sm_list, bodies_to_modify, gme);
        remove_pushed_attributes(result_sm_list, bodies_to_modify);
      }
      status = finish_webcut( engine_bodies, result_sm_list, merge, status,
                              results_list, &merged_surface_ids, &merged_curve_ids );
    }
    else
    {
      if(!preview)
        remove_pushed_attributes(result_sm_list, engine_bodies);
    }


    engine_bodies.clean_out();
    engine_body_sms.clean_out();
    result_sm_list.clean_out();

    if ( status == CUBIT_FAILURE )
    {
      rval = CUBIT_FAILURE;
      break;
    }
  }

  if(!preview)
    do_attribute_cleanup();

  if( CubitUndo::get_undo_enabled() )
  {
    if( rval == CUBIT_SUCCESS )
      CubitUndo::note_result_bodies( results_list );
    else
      CubitUndo::remove_last_undo();
  }

  return rval;
}

CubitStatus GeometryModifyTool::webcut_w_cylinder(
                                        DLIList<BodySM*> &webcut_body_list,
                                        double radius,
                                        const CubitVector &axis,
                                        const CubitVector &center,
                                        DLIList<BodySM*>& neighbor_imprint_list,
                                        DLIList<BodySM*>& results_list,
                                        ImprintType imprint_type )
{
  GeometryModifyEngine* gme = 0;
  gme = get_engine(webcut_body_list.get()); 
  assert(gme);

  double max_size =  0.;
  //lets find the distance to the center for each body and take
  //the max.
  double curr;
  CubitVector cent_bod;
  CubitBox bounding_box;
  BodySM *body_ptr;
  bounding_box = webcut_body_list[0]->bounding_box();
  cent_bod =  bounding_box.center();
  cent_bod = cent_bod - center;
  curr = cent_bod.length();
  if ( curr > max_size )
     max_size = curr;


  for ( int ii = webcut_body_list.size()-1; ii > 0; ii-- )
    {
      body_ptr = webcut_body_list[ii];
      bounding_box |= body_ptr->bounding_box();
      cent_bod = body_ptr->bounding_box().center();
      cent_bod = cent_bod - center;
      curr = cent_bod.length();
      if ( curr > max_size )
        max_size = curr;
    }

  curr = bounding_box.diagonal().length();

  if ( curr > max_size )
     max_size = curr;

  double height = 0.0;
  if ( center.x() > max_size )
    {
      height = 500.0 * center.x();
    }
  else if ( center.y() > max_size )
    {
      height = 500.0 * center.y();
    }
  else if ( center.z() > max_size )
    {
      height = 500.0 * center.z();
    }
  else
    {
      height = 500.0 * max_size;
    }

  //lets make certain we have a valid height..
  if ( height < GEOMETRY_RESABS )
    {
      height = 500.0;
    }

  BodySM *cutting_tool_ptr = gme->cylinder( height, radius, radius, radius );

  if( cutting_tool_ptr == NULL )
    return CUBIT_FAILURE;

  //transform the cyclinder to cernter and axis
  // The current frustum is centered on the z axis.
  CubitVector axis2(0., 0., 1.0 );
  //now find the normal to the current axis and axis we want to be
  //at. This normal is where we will rotate about.
  CubitVector normal_axis = axis2 * axis;
  if ( normal_axis.length() > CUBIT_RESABS )
    {
       //angle in degrees.
       double angle = normal_axis.vector_angle( axis2, axis );
       gme->get_gqe()->rotate(cutting_tool_ptr, normal_axis, angle);
    }
  gme->get_gqe()->translate(cutting_tool_ptr, center);

  CubitStatus stat = gme->webcut( webcut_body_list, cutting_tool_ptr, 
              neighbor_imprint_list, results_list, imprint_type) ;

  // Delete the BodySM that was created to be used as a tool
  gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

  return stat;
}

void GeometryModifyTool::do_attribute_setup(void)
{
  //save attribute settings
  CGMApp::instance()->save_current_attribute_states();
  
  //Turn off all attributes
  CubitAttribManager *attrib_manager = CGMApp::instance()->attrib_manager();
  attrib_manager->set_all_auto_update_flags( CUBIT_FALSE ); 
  attrib_manager->set_all_auto_actuate_flags( CUBIT_FALSE ); 
  attrib_manager->set_all_auto_write_flags( CUBIT_FALSE ); 
  attrib_manager->set_all_auto_read_flags( CUBIT_FALSE );  

  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_ENTITY_NAME, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_ENTITY_NAME, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_ENTITY_NAME, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_ENTITY_NAME, CUBIT_TRUE);

  // enable update, actuate, write, read for composite attributes
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_COMPOSITE_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_COMPOSITE_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_COMPOSITE_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_COMPOSITE_VG, CUBIT_TRUE);
  // enable update, actuate, write, read for partition attributes
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_PARTITION_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_PARTITION_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_PARTITION_VG, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_PARTITION_VG, CUBIT_TRUE);
  // enable update, actuate, write, read for entity ids
  // We will use these ID atts to make sure the ref entities associated with unmodified
  // virtual bridges will maintain the same ids after real operations.
  CGMApp::instance()->attrib_manager()->set_auto_update_flag(CA_ENTITY_ID, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(CA_ENTITY_ID, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_write_flag(CA_ENTITY_ID, CUBIT_TRUE);
  CGMApp::instance()->attrib_manager()->set_auto_read_flag(CA_ENTITY_ID, CUBIT_TRUE);
}

void GeometryModifyTool::do_attribute_cleanup(void)
{
  CGMApp::instance()->restore_previous_attribute_states();
}

// Push the virtual geometry attributes down to the underlying solid model
// topology.
void GeometryModifyTool::push_vg_attributes_before_modify(DLIList<BodySM*> &old_sms)
{
  // Get all of the ref entities involved and push CA_ENTITY_ID attributes
  // on to them.  This will help us maintain ids on virtual that doesn't
  // get modified by the real operation.
  int i;
  DLIList<RefVolume*> volume_list; 

  for(i=old_sms.size(); i--;)
  {
    BodySM *bsm = old_sms.get_and_step();
    Body *body = dynamic_cast<Body*>(bsm->topology_entity());
    if(body)
    {
      // Append to the total list of volumes.
      body->ref_volumes(volume_list);
    }
  }
  // get all child entities (only get entities below volumes)
  DLIList<RefEntity*> child_list, ref_ent_list;
  CAST_LIST_TO_PARENT(volume_list, ref_ent_list);
  RefEntity::get_all_child_ref_entities( ref_ent_list, child_list );  

  // Only push the id attributes if we are doing persistent ids.
  if(!get_new_ids())
    CubitAttribUser::auto_update_cubit_attrib(child_list);

  DLIList<TopologyBridge*> top_bridges;
  CAST_LIST_TO_PARENT(old_sms, top_bridges);
  GeometryQueryTool::instance()->ige_export_geom(top_bridges);
}

// Push imprint specific attributes down to the underlying solid model
// topology prior to imprinting.
void GeometryModifyTool::push_imprint_attributes_before_modify(DLIList<BodySM*> &body_sms)
{
  GeometryQueryTool::instance()->ige_push_imprint_attributes_before_modify(body_sms);
}

void GeometryModifyTool::push_named_attributes_to_curves_and_points(DLIList<TopologyBridge*> &tb_list,
                                                             const char *name_in)
{
  GeometryQueryTool::instance()->ige_push_named_attributes_to_curves_and_points(tb_list, name_in);
}

// Cleanup imprint attributes that were pushed onto the underlying solid
// model topology.
void GeometryModifyTool::remove_imprint_attributes_after_modify(DLIList<BodySM*> &old_sms,
                                                                DLIList<BodySM*> &new_sms)
{
  GeometryQueryTool::instance()->ige_remove_imprint_attributes_after_modify(old_sms, new_sms);
}

#ifdef CAT
CubitStatus GeometryModifyTool::webcut_across_translate(
                                          DLIList<Body*> &webcut_body_list,
                                          RefFace *ref_face_top,
                                          RefFace *ref_face_bottom,
                                          DLIList<Body*> &results_list,
                                          ImprintType imprint_type,
                                          CubitBoolean merge,
                                          CubitBoolean preview)
{
  if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
    return CUBIT_FAILURE;

  DLIList<Body*> original_body_list = webcut_body_list;
  const int count = webcut_body_list.size() + 2;
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);
  CAST_LIST_TO_PARENT(webcut_body_list, entity_list);
  entity_list.append( ref_face_top );
  entity_list.append( ref_face_bottom );
  GeometryModifyEngine* gme = common_modify_engine( entity_list, bridge_list );

  if ( !gme )
  {
     PRINT_ERROR("Performing WEBCUTS on volumes containing geometry from\n"
                 "different modeling engines is not allowed.\n"
                 "Delete uncommon geometry on these volumes before operation.\n\n");
     return CUBIT_FAILURE;
  }

  Surface* surface_bottom = dynamic_cast<Surface*>(bridge_list.pop());
  Surface* surface_top    = dynamic_cast<Surface*>(bridge_list.pop());
  DLIList<BodySM*> result_sm_list, webcut_sm_list(bridge_list.size());
  CAST_LIST(bridge_list, webcut_sm_list, BodySM);

  if (!preview)
  {
    do_attribute_setup();
    push_vg_attributes_before_modify(webcut_sm_list);
  }

  DLIList<BodySM*> neighbor_list;
  CubitStatus result_val = gme->webcut_across_translate (
    webcut_sm_list, surface_top, surface_bottom, result_sm_list, neighbor_list, imprint_type, preview );

  // if we're doing the real thing (not preview) finish the process
  if (!preview)
  {
    restore_vg_after_modify(result_sm_list, original_body_list, gme);
    remove_pushed_attributes(result_sm_list, original_body_list);
    result_val = finish_webcut( webcut_body_list, result_sm_list, merge,
                                result_val, results_list);
    do_attribute_cleanup();
  }

  return result_val;
}
#endif

//-------------------------------------------------------------------------
// Purpose       : This functions webcuts a list of bodies using a loop
//                 of curves
//
// Special Notes :
//
// Creator       : Eric W. Nielsen
//
// Creation Date : 12/20/99
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::webcut_with_curve_loop(
                                         DLIList<Body*>& webcut_body_list,
                                         DLIList<RefEdge*>& refedge_list,
                                         DLIList<Body*>& results_list,
                                         DLIList<Body*> &neighboring_bodies,
                                         ImprintType imprint_type,
                                         CubitBoolean merge,
                                         CubitBoolean preview)

{
  if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
    return CUBIT_FAILURE;

   int i;
   const int count = webcut_body_list.size() + refedge_list.size();
   DLIList<TopologyEntity*> entity_list(count);
   DLIList<TopologyBridge*> bridge_list(count);
   CAST_LIST_TO_PARENT(webcut_body_list, entity_list);
   refedge_list.reset();
   for (i = refedge_list.size(); i--; )
     entity_list.append(refedge_list.get_and_step());
   GeometryModifyEngine* gme = common_modify_engine( entity_list, bridge_list );

   if ( !gme )
   {
      PRINT_ERROR("Performing WEBCUTS on volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
   {
     DLIList<Body*> bodies_to_save;
     bodies_to_save += webcut_body_list;
     bodies_to_save += neighboring_bodies;
     CubitUndo::save_state_with_cubit_file( bodies_to_save );
   }

   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), result_sm_list;
   DLIList<Curve*> curve_list(refedge_list.size());
   CAST_LIST(bridge_list, body_sm_list, BodySM);
   CAST_LIST(bridge_list, curve_list, Curve);
   assert(body_sm_list.size() == webcut_body_list.size());
   assert(curve_list.size() == refedge_list.size());

   DLIList<int> merged_surface_ids;
   DLIList<int> merged_curve_ids;
   DLIList<BodySM*> neighbor_imprint_list;
   DLIList<Body*> bodies_to_modify;

   //make copies of the curves.
   DLIList<Curve*> copied_curves;
   Curve* temp_curve = NULL;
   for (int i = curve_list.size(); i--;)
     {
       temp_curve = gme->make_Curve(curve_list.get_and_step());
       if(temp_curve != NULL)
         copied_curves.append(temp_curve);
     }

   //make a face out of the curves
   Surface * surf = gme->make_Surface(PLANE_SURFACE_TYPE, copied_curves, NULL, false );
   if (surf == NULL)
     {
       PRINT_ERROR("webcut tool surface is not created from acis.\n");
       return CUBIT_FAILURE;
     }


   //get cutting tool BodySM.
   BodySM* cutting_tool_ptr = gme->make_BodySM(surf);
   assert(cutting_tool_ptr );

   if (!preview)
   {
     int i;
     for( i=neighboring_bodies.size(); i--; )
     {
       Body *neighbor_body = neighboring_bodies.get_and_step();
       BodySM *tmp_body = neighbor_body->get_body_sm_ptr(); 
       GeometryModifyEngine *neighbor_gme = get_engine( tmp_body );
       
       if( gme == neighbor_gme )
         neighbor_imprint_list.append( tmp_body ); 
     }

     do_attribute_setup();
     DLIList<BodySM*> bodies_sm_to_modify;
     bodies_sm_to_modify += body_sm_list;
     bodies_sm_to_modify += neighbor_imprint_list;
     push_vg_attributes_before_modify( bodies_sm_to_modify );
     bodies_to_modify += webcut_body_list;
     bodies_to_modify += neighboring_bodies;
     get_merged_curve_and_surface_ids( bodies_to_modify, merged_surface_ids, merged_curve_ids );
   }
  

   CubitStatus result_val = gme->webcut(
                     body_sm_list, cutting_tool_ptr,
                     neighbor_imprint_list,
                     result_sm_list,
                     imprint_type, preview) ;

   // Delete the BodySM that was created to be used as a tool
   gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

   // finish up if we're doing the real thing
   if (!preview)
   {
     restore_vg_after_modify(result_sm_list, bodies_to_modify, gme);
     remove_pushed_attributes(result_sm_list, bodies_to_modify);

     result_val = finish_webcut( webcut_body_list, result_sm_list, merge,
                                 result_val, results_list,
                                 &merged_surface_ids, &merged_curve_ids );
     do_attribute_cleanup();
   }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( result_val == CUBIT_SUCCESS )
      CubitUndo::note_result_bodies( results_list );
    else
      CubitUndo::remove_last_undo();
  }

   return result_val;
}

//-------------------------------------------------------------------------
// Purpose       : This functions webcuts a list of bodies through a plane
//                 defined by a refFace. The newly created bodies are
//                 merged and imprinted depeding on the respective flags.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/16/96
//-------------------------------------------------------------------------

CubitStatus GeometryModifyTool::webcut_with_surface(
                                      DLIList<Body*>& webcut_body_list,
                                      RefFace* refFace,
                                      DLIList<Body*>& results_list,
                                      DLIList<Body*> &neighboring_bodies,
                                      ImprintType imprint_type,
                                      CubitBoolean merge,
                                      CubitBoolean preview)
{
  //make sure that this refface is planar, or you'll get unexpected results
  if( refFace->geometry_type() != PLANE_SURFACE_TYPE )
  {
    PRINT_ERROR("Surface %d is not planar.\n", refFace->id() );
    return CUBIT_FAILURE;
  }

     // Using the face, get three positions (CubitVectors :) :) )
     // The positions can be used by a GeometryModifyEngine to generate its
     // own webcutting tool and perform the webcut.
   CubitVector vector1 = refFace->position_from_u_v(0, 1) ;
   CubitVector vector3 = refFace->position_from_u_v(1, 0) ;
   CubitVector vector2 = refFace->position_from_u_v(0, 0) ;

   CubitStatus result_val = this->webcut_with_plane(webcut_body_list, vector1,
                                  vector2, vector3, results_list, neighboring_bodies,
                                  imprint_type, merge, preview) ;

   return result_val;
}

//-------------------------------------------------------------------------
// Purpose       : This functions webcuts a list of bodies using a plane
//                 defined by a refFace. The newly created bodies are
//                 merged and imprinted depending on the respective flags.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/16/96
//-------------------------------------------------------------------------

CubitStatus GeometryModifyTool::webcut_with_body(
                                   DLIList<Body*>& webcut_body_list,
                                   Body* tool_body,
                                   DLIList<Body*>& results_list,
                                   DLIList<Body*> &neighboring_bodies,
                                   ImprintType imprint_type,
                                   CubitBoolean merge,
                                   CubitBoolean preview)
{

  CubitStatus ret;

  if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
    return CUBIT_FAILURE;

  DLIList<BodySM*> body_sm_list(webcut_body_list.size()), result_sm_list;
  BodySM* tool_sm = tool_body->get_body_sm_ptr();
  GeometryModifyEngine* gme = 0;
  if (!tool_sm || !(gme = get_engine(tool_sm)))
  {
    PRINT_DEBUG_153("No modify engine available for body %d.\n", tool_body->id());
    PRINT_ERROR("No modify engine available for volume %d.\n", tool_body->ref_volume()->id());
    return CUBIT_FAILURE;
  }

  webcut_body_list.reset();
  int i;
  for (i = webcut_body_list.size(); i--; )
  {
    Body* body_ptr = webcut_body_list.get_and_step();
    BodySM* sm_ptr = body_ptr->get_body_sm_ptr();
    if (!sm_ptr || get_engine(sm_ptr) != gme)
    {
      PRINT_DEBUG_153("ERROR: Body %d does not belong to the same geometric engine "
                  "as body %d\n", body_ptr->id(), tool_body->id());
      PRINT_ERROR("Volume %d does not belong to the same geometric engine "
                  "as volume %d\n", body_ptr->ref_volume()->id(), tool_body->ref_volume()->id());
      return CUBIT_FAILURE;
    }

    body_sm_list.append(sm_ptr);
  }

  DLIList<int> merged_surface_ids;
  DLIList<int> merged_curve_ids;
  DLIList<BodySM*> neighbor_imprint_list;
  DLIList<Body*> bodies_to_modify;

  if (!preview)
  {
    if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
    {
      DLIList<Body*> bodies_to_save;
      bodies_to_save += webcut_body_list;
      bodies_to_save += neighboring_bodies;
      CubitUndo::save_state_with_cubit_file( bodies_to_save );
    }

    int i;
    for( i=neighboring_bodies.size(); i--; )
    {
      Body *neighbor_body = neighboring_bodies.get_and_step();
      BodySM *tmp_body = neighbor_body->get_body_sm_ptr(); 
      GeometryModifyEngine *neighbor_gme = get_engine( tmp_body );
      
      if( gme == neighbor_gme )
        neighbor_imprint_list.append( tmp_body ); 
    }

    do_attribute_setup();
    DLIList<BodySM*> bodies_sm_to_modify;
    bodies_sm_to_modify += body_sm_list;
    bodies_sm_to_modify += neighbor_imprint_list;
    push_vg_attributes_before_modify( bodies_sm_to_modify );
    bodies_to_modify += webcut_body_list;
    bodies_to_modify += neighboring_bodies;
    get_merged_curve_and_surface_ids( bodies_to_modify, merged_surface_ids, merged_curve_ids );
  }


  int count = gme->webcut(body_sm_list, tool_sm,
                          neighbor_imprint_list,
                          result_sm_list, imprint_type, preview);

  // finish up if we're doing the real thing
  if (!preview)
  {
    restore_vg_after_modify(result_sm_list, bodies_to_modify, gme);
    remove_pushed_attributes(result_sm_list, bodies_to_modify );

    ret = finish_webcut(webcut_body_list, result_sm_list, merge,
                          count > 0 ? CUBIT_SUCCESS : CUBIT_FAILURE,
                          results_list, &merged_surface_ids, &merged_curve_ids );

    do_attribute_cleanup();
  }
  else
    ret = count > 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;

  if( preview == CUBIT_FALSE && CubitUndo::get_undo_enabled() )
  {
    if( ret == CUBIT_SUCCESS )
      CubitUndo::note_result_bodies( results_list );
    else
      CubitUndo::remove_last_undo();
  }

  return ret;
}

CubitStatus GeometryModifyTool::section(
                                   DLIList<Body*> &section_body_list,
                                   const CubitVector &point_1,
                                   const CubitVector &point_2,
                                   const CubitVector &point_3,
                                   DLIList<Body*> &new_body_list,
                                   CubitBoolean keep_normal_side,
                                   CubitBoolean keep_old )
{
  if (!okay_to_modify( section_body_list, "SECTION" ))
    return CUBIT_FAILURE;
   CubitStatus rval = CUBIT_SUCCESS;

   const int count = section_body_list.size();
   DLIList<BodySM*> result_sm_list;
   DLIList<Body*> body_list(section_body_list);
   DLIList<BodySM*> engine_body_sms(count);
   DLIList<Body*> engine_bodies(count);
   GeometryModifyEngine* gme = 0;

   if( CubitUndo::get_undo_enabled() )
   {
     if( keep_old )
       CubitUndo::save_state();
     else
       CubitUndo::save_state_with_cubit_file( section_body_list );
   }

   while ( (gme = group_bodies_by_engine(body_list, engine_bodies, engine_body_sms)) )
   {
      //get all the child entities that have been merged
      DLIList<int> merged_surface_ids;
      DLIList<int> merged_curve_ids;
      get_merged_curve_and_surface_ids( engine_bodies, merged_surface_ids, merged_curve_ids );

      CubitStatus result = gme->section( engine_body_sms,
                                         point_1, point_2, point_3,
                                         result_sm_list,
                                         keep_normal_side, keep_old );

      if (!finish_sm_op( engine_bodies, result_sm_list, new_body_list ))
        result = CUBIT_FAILURE;
      if (!result)
        rval = CUBIT_FAILURE;

      if( merged_surface_ids.size() || merged_curve_ids.size() )
        fixup_merged_entities( merged_surface_ids, merged_curve_ids);

      engine_body_sms.clean_out();
      engine_bodies.clean_out();
      result_sm_list.clean_out();
   }

  if( CubitUndo::get_undo_enabled() )
  {
    if( new_body_list.size() ) //if there are new bodies...something succeeded
      CubitUndo::note_result_bodies( new_body_list );
    else
      CubitUndo::remove_last_undo();
  }

   return rval;
}

CubitStatus
GeometryModifyTool::offset_curves( DLIList<RefEdge*>& ref_edge_list,
                                   double offset_distance,
                                   const CubitVector& offset_direction,
                                   int gap_type )
{
  // Make sure all curves are from same geometry engine
  DLIList<TopologyEntity*> entity_list(ref_edge_list.size());
  DLIList<TopologyBridge*> bridge_list(ref_edge_list.size());
  CAST_LIST_TO_PARENT(ref_edge_list, entity_list);
  GeometryModifyEngine* gme_ptr = common_modify_engine(entity_list, bridge_list);
  if ( !gme_ptr )
  {
    PRINT_ERROR("Can't create offset curves that don't all originate from same\n"
      "       modeling engine.\n");
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state();

  DLIList<Curve*> curve_list(bridge_list.size()), result_list;
  CAST_LIST(bridge_list, curve_list, Curve);
  assert(curve_list.size() == ref_edge_list.size());

  CubitStatus rval = gme_ptr->offset_curves( curve_list, result_list,
                        offset_distance, offset_direction, gap_type );
  assert( rval || !result_list.size() );
  result_list.reset();

  DLIList<RefEntity*> created_edges;
  for (int i = result_list.size(); i--; )
  {
    RefEdge *new_edge = GeometryQueryTool::instance()->make_free_RefEdge(result_list.get_and_step());
    created_edges.append( new_edge );
  }

  if( CubitUndo::get_undo_enabled() )
  {
    if( created_edges.size() )
      CubitUndo::note_result_entities( created_edges );
    else
      CubitUndo::remove_last_undo();
  }

  return rval;
}

CubitStatus GeometryModifyTool::align_body( Body *body_ptr,
                                            RefFace* my_face,
                                            RefFace *target_face,
                                            CubitVector &my_center,
                                            CubitVector &axis,
                                            CubitVector &target_center,
                                            double &angle)
{
     // Initialize angle to zero, so we are passing back a valid
     // value if we never need to rotate the body
   angle = 0.0;

     //Align the body, which contains my_face, with the target_face.

     //Essentially we want to align the body at my_face with the target_face.
   my_center = my_face->center_point();
   target_center = target_face->center_point();

     //Get the correct volume from my_volume.
   DLIList<RefVolume*> ref_vols;
   body_ptr->ref_volumes(ref_vols);
   DLIList<RefFace*> ref_face_list;
   RefVolume *ref_vol = NULL;

   CubitBoolean found_vol = CUBIT_FALSE;

   for ( int ii = ref_vols.size(); ii > 0; ii-- )
   {
      ref_vol = ref_vols.get_and_step();
      ref_face_list.clean_out();
      ref_vol->ref_faces( ref_face_list );
      if ( ref_face_list.move_to(my_face) )
      {
         found_vol = CUBIT_TRUE;
         break;
      }
   }
   if ( !found_vol )
   {
      PRINT_ERROR( "Cannot find appropriate volume that %s (surface%d) is on.\n",
                   my_face->entity_name().c_str(),
                   my_face->id() );
      PRINT_DEBUG_153( "%s (surface%d) is not on %s (body %d)\n",
                       my_face->entity_name().c_str(),
                       my_face->id(),
                       body_ptr->entity_name().c_str(),
                       body_ptr->id() );
      return CUBIT_FAILURE;
   }
   CubitVector my_vector = my_face->normal_at( my_center, ref_vol );
   CubitVector target_vector = target_face->normal_at( target_center );

   if( CubitUndo::get_undo_enabled() )
   {
     DLIList<Body*> align_list(1);
     align_list.append( body_ptr );
     CubitUndo::save_state_with_cubit_file( align_list );
   }

   axis = my_vector * target_vector;
     //Move the body so that the center of my_face
     //is at the origin.
   GeometryQueryTool::instance()->translate( body_ptr, -my_center );

     //If the axis length is zero we dont want to do a rotation,
     //they are already parallel..
   if (axis.length() > CUBIT_RESABS )
   {
      angle = axis.vector_angle( my_vector, target_vector );

        //Now we have the angle and the axis to rotate about...
      angle = 180.0*angle/CUBIT_PI;
      
      //if the angle is greater 180, we want the normals pointing
      //into one another, not in the same direction
      if( angle > 90 )
        angle = angle - 180;

        //Now rotate the body about the axis.
      GeometryQueryTool::instance()->rotate( body_ptr, axis, angle );
   }
     //Now move the body to the location of the target_center.
   GeometryQueryTool::instance()->translate( body_ptr, target_center );

     //That should do it.
   return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::align_body(Body *body_ptr,
                                           DLIList<RefEntity*>& ref_ent_list,
                                           CubitVector first_vector,
                                           CubitVector second_vector,
                                           CubitVector &my_center_1,
                                           CubitVector &my_center_2,
                                           CubitVector &axis_of_rot,
                                           double &angle,
                                           double &angle_2)
{
   DLIList<RefFace*> ref_face_list;
   DLIList<RefVertex*> ref_vertex_list;

   CAST_LIST(ref_ent_list, ref_face_list, RefFace);
   CAST_LIST(ref_ent_list, ref_vertex_list, RefVertex);

   angle = 0.0;
   angle_2 = 0.0;

   if( ref_face_list.size() == 1 &&
      ref_vertex_list.size() == 0 )
   {

      RefFace *my_face = ref_face_list.get_and_step();
      my_center_1 = my_face->center_point();
      my_center_2 = my_center_1;

      //Get the correct volume from my_volume.
      DLIList<RefVolume*> ref_vols;
      body_ptr->ref_volumes(ref_vols);
      DLIList<RefFace*> face_list;
      RefVolume *ref_vol = NULL;

      CubitBoolean found_vol = CUBIT_FALSE;

      for ( int ii = ref_vols.size(); ii > 0; ii-- )
      {
         ref_vol = ref_vols.get_and_step();
         face_list.clean_out();
         ref_vol->ref_faces( face_list );

         if ( face_list.move_to(my_face) )
         {
            found_vol = CUBIT_TRUE;
            break;
         }
      }

      if ( !found_vol )
      {
         PRINT_ERROR( "Can't find appropriate volume that surface %d is on.\n",
                      my_face->id() );

         PRINT_DEBUG_153( "surface%d is not on body %d\n",
            my_face->id(),
            body_ptr->id() );
         return CUBIT_FAILURE;
      }

      if( CubitUndo::get_undo_enabled() )
      {
        DLIList<Body*> align_list(1);
        align_list.append( body_ptr );
        CubitUndo::save_state_with_cubit_file( align_list );
      }

      //translate body so surface's centroid is at 0,0,0
      GeometryQueryTool::instance()->translate( body_ptr, -my_center_1 );

      //get axis of rotation
      CubitVector my_vector = my_face->normal_at( my_center_1, ref_vol );
      my_vector.normalize();
      axis_of_rot = first_vector * my_vector;

      //get angle of rotation and rotate
      if (axis_of_rot.length() > CUBIT_RESABS )
      {
         angle = axis_of_rot.vector_angle( my_vector, first_vector );
         //Now we have the angle and the axis to rotate about...
         angle = 180.0*angle/CUBIT_PI;
         //Now rotate the body about the axis.
         GeometryQueryTool::instance()->rotate( body_ptr, axis_of_rot, angle );
      }
      else  //is aligned already----flip 180 degrees
      {
        axis_of_rot = second_vector;
        angle = 180.0;
        GeometryQueryTool::instance()->rotate( body_ptr, second_vector, angle );
      }

      if( CubitUndo::get_undo_enabled() )
      {
        DLIList<Body*> align_list(1);
        align_list.append( body_ptr );
        CubitUndo::save_state_with_cubit_file( align_list );
      }

      //translate body back
      GeometryQueryTool::instance()->translate( body_ptr, my_center_1 );

      //That should do it.
      return CUBIT_SUCCESS;
  }
  else if( ref_face_list.size() == 0 &&
     ref_vertex_list.size() == 2 )
  {

     RefVertex *my_vertex_1 = ref_vertex_list.get_and_step();
     RefVertex *my_vertex_2 = ref_vertex_list.get_and_step();
     my_center_1 = my_vertex_1->center_point();
     CubitVector my_center_ver_2 = my_vertex_2->center_point();

     //Get the correct volume from my_volume.
     DLIList<RefVolume*> ref_vols;
     body_ptr->ref_volumes(ref_vols);
     DLIList<RefVertex*> vertex_list;
     RefVolume *ref_vol;

     CubitBoolean found_vol = CUBIT_FALSE;

     for ( int ii = ref_vols.size(); ii > 0; ii-- )
     {
        ref_vol = ref_vols.get_and_step();
        vertex_list.clean_out();
        ref_vol->ref_vertices( vertex_list );

        if ( vertex_list.move_to(my_vertex_1)&&
           vertex_list.move_to(my_vertex_2))
        {
           found_vol = CUBIT_TRUE;
           break;
        }
     }

     if ( !found_vol )
     {
        PRINT_ERROR( "Can't find appropriate volume that vertex %d or %d is on.\n",
                     my_vertex_1->id(),
                     my_vertex_2->id() );

        PRINT_DEBUG_153( "vertex %d or %d is not on body %d\n",
           my_vertex_1->id(),
           my_vertex_2->id(),
           body_ptr->id() );
        return CUBIT_FAILURE;
     }

     CubitVector my_vector = my_center_ver_2 - my_center_1;
     CubitVector target_vector = first_vector;
     axis_of_rot = my_vector * target_vector;

     if( CubitUndo::get_undo_enabled() )
     {
       DLIList<Body*> align_list(1);
       align_list.append( body_ptr );
       CubitUndo::save_state_with_cubit_file( align_list );
     }

     //Move the body so that the center of the vertex
     //is at the origin.
     GeometryQueryTool::instance()->translate( body_ptr, -my_center_1 );

     my_center_2 = my_vertex_1->center_point();
     CubitVector temp_center_ver_2 = my_vertex_2->center_point();
     my_vector = temp_center_ver_2 - my_center_2;

     my_vector.normalize();
     target_vector.normalize();
     CubitVector test_vec = my_vector - target_vector;

     if (axis_of_rot.length() > CUBIT_RESABS )
     {
        angle = axis_of_rot.vector_angle( my_vector, target_vector );
        //Now we have the angle and the axis to rotate about...
        angle = 180.0*angle/CUBIT_PI;
        //Now rotate the body about the axis.
        GeometryQueryTool::instance()->rotate( body_ptr, axis_of_rot, angle );
     }
     else if ( test_vec.normalize() > CUBIT_RESABS)
     {
        axis_of_rot = my_center_1 * target_vector;
        angle += 180.0;
        //Now rotate the body about the axis.
        GeometryQueryTool::instance()->rotate( body_ptr, axis_of_rot, 180.0 );
     }

     CubitVector body_center = body_ptr->center_point();
     body_center.normalize();
     CubitVector temp_vec = first_vector * body_center;

     // room for improvment
     temp_vec.set( temp_vec.x()*temp_vec.x(),
        temp_vec.y()*temp_vec.y(),
        temp_vec.z()*temp_vec.z() );

     body_center += temp_vec;

     if( body_center.x() < -CUBIT_RESABS ||
        body_center.y() < -CUBIT_RESABS ||
        body_center.z() < -CUBIT_RESABS )
     {
        GeometryModifyEngine *gePtr1 = get_engine(body_ptr);
        DLIList<Surface*> bad_face_list;
        body_ptr->get_body_sm_ptr()->surfaces(bad_face_list);
        angle_2= 180.0;
        //Now rotate the body about the axis.
        GeometryQueryTool::instance()->rotate( body_ptr, first_vector, angle_2 );
        if (bad_face_list.size() && body_ptr->is_sheet_body())
        {
           CubitStatus flip_result = gePtr1->flip_normals(bad_face_list);
           PRINT_WARNING("Reversing all surface normals within sheet body %d \n"
                         "to remain in the positive plane\n",body_ptr->id());
           if ( flip_result == CUBIT_FAILURE )
           {
             if( CubitUndo::get_undo_enabled() )
               CubitUndo::remove_last_undo();

             return CUBIT_FAILURE;
           }
        }
     }
     return CUBIT_SUCCESS;
  }
  else
  {
     PRINT_ERROR("Incorrect syntax type 'help align'. Nothing aligned\n");
     return CUBIT_FAILURE;
  }
}


CubitStatus GeometryModifyTool::sweep_setup(
                         const char* name,
                         DLIList<RefEntity*>& input_entity_list,
                         DLIList<Body*>& output_body_list,
                         GeometryModifyEngine*& output_engine,
                         CubitBoolean& changed_new_ids,
                         DLIList<GeometryEntity*>& output_geom_list,
                         DLIList<RefEdge*>* input_edge_list,
                         DLIList<Curve*>* output_curve_list )
{
  int i;

  if (input_entity_list.size() == 0)
    return CUBIT_FAILURE;

  DLIList<RefFace*> ref_face_list;
  DLIList<RefEdge*> ref_edge_list;
  CAST_LIST( input_entity_list, ref_face_list, RefFace );
  CAST_LIST( input_entity_list, ref_edge_list, RefEdge );
  if( ref_face_list.size() && ref_edge_list.size() )
  {
    PRINT_ERROR( "Currently cannot sweep curves and surfaces at the same time\n" );
    return CUBIT_FAILURE;
  }

  if( (ref_face_list.size()+ref_edge_list.size()) != input_entity_list.size() )
  {
    PRINT_ERROR( "Can only sweep surfaces or curves\n" );
    return CUBIT_FAILURE;
  }

  // Currently faces to be swept have to belong to different bodies.  It
  // doesn't work well to sweep multiple faces from the same body (also
  // current implementation in AGE crashes).
  if( ref_face_list.size() )
  {
    int free_surf_cnt = 0; // Number of free surfaces being swept
    DLIList<Body*> body_list;
    RefFace *ref_face_ptr;
    for( i=ref_face_list.size(); i--; )
    {
      ref_face_ptr = ref_face_list.get_and_step();
      DLIList<Body*> tmp_body_list;
      ref_face_ptr->bodies( tmp_body_list );
      if( tmp_body_list.size() > 1 )
      {
        PRINT_ERROR( "Surface %d belongs to more than one volume; cannot sweep.\n",
          ref_face_ptr->id() );
        return CUBIT_FAILURE;
      }

      if( tmp_body_list.size() > 0 )
        body_list.append_unique( tmp_body_list.get() );
      else
        free_surf_cnt++;
    }

    if( body_list.size()+free_surf_cnt != ref_face_list.size() )
    {
      PRINT_ERROR( "The surfaces to be swept cannot belong to the same volume\n" );
      return CUBIT_FAILURE;
    }

    output_body_list = body_list;
  }

  // Make sure all entities are from same modify engine
  DLIList<TopologyEntity*> te_list;
  CAST_LIST(input_entity_list, te_list, TopologyEntity);

  // This is for sweeping along a curve list
  if (input_edge_list)
  {
    input_edge_list->reset();
    for( i=input_edge_list->size(); i--; )
      te_list.append( input_edge_list->get_and_step() );
  }

  DLIList<TopologyBridge*> bridge_list( te_list.size() );
  output_engine = common_modify_engine( te_list, bridge_list );
  if( output_engine == NULL )
  {
     PRINT_ERROR("Can't sweep with entities from different modeling engines.\n");
     return CUBIT_FAILURE;
  }

  // Check for virtual anywhere in participating bodies.  This will catch cases where
  // the result from sweeping would interact with some virtual somewhere. bwc 12/22/05.
  for(i=output_body_list.size(); i--;)
  {
     if(GeometryQueryTool::instance()->contains_intermediate_geometry(
        output_body_list.get_and_step()))
     {
        PRINT_ERROR("Can't sweep faces of bodies containing virtual geometry.\n");
        return CUBIT_FAILURE;
     }
  }

  changed_new_ids = CUBIT_FALSE;

  TopologyBridge* bridge_ptr;
  bridge_list.reset();
  for( i=input_entity_list.size(); i--; )
  {
    bridge_ptr = bridge_list.get_and_step();
    GeometryEntity* geom_ptr = dynamic_cast<GeometryEntity*>(bridge_ptr);
    output_geom_list.append(geom_ptr);
  }

  if( input_edge_list )
  {
    for( i=input_edge_list->size(); i--; )
    {
      bridge_ptr = bridge_list.get_and_step();
      Curve* curve_ptr = dynamic_cast<Curve*>(bridge_ptr);
      output_curve_list->append(curve_ptr);
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::sweep_finish(
                                     const char* const,
                                     DLIList<Body*>& input_body_list,
                                     DLIList<BodySM*>& new_body_list,
                                     DLIList<Body*>& output_body_list,
                                     CubitBoolean changed_new_ids )
{
  int i;
  DLIList<BodySM*> regen_list( new_body_list );
  GeometryQueryTool *gqt = GeometryQueryTool::instance();

  input_body_list.reset();
  for (i = input_body_list.size(); i--; )
  {
    Body* body = input_body_list.get_and_step();
    BodySM* bodysm = body->get_body_sm_ptr();
    if (bodysm)
    {
      regen_list.append_unique(bodysm);
      remove_dead_entity_names(body);
      output_body_list.append( body );
    }
    else
    {
      gqt->destroy_dead_entity(body);
    }
  }
  gqt->cleanout_deactivated_geometry();

  regen_list.reset();
  for (i = regen_list.size(); i--; )
  {
    BodySM* bodysm = regen_list.get_and_step();
    Body* body = gqt->make_Body(bodysm);
    output_body_list.append( body );
    PRINT_INFO("%s volume %d\n",
      new_body_list.is_in_list(bodysm) ? "Created swept" : "Updated",
      body->ref_volume()->id());
  }

  if (changed_new_ids) set_new_ids(CUBIT_FALSE);
  gqt->cleanout_deactivated_geometry();

  return CUBIT_SUCCESS;
}





CubitStatus GeometryModifyTool::sweep_rotational(
                                        DLIList<RefEntity*>& ref_ent_list,
                                        const CubitVector& point,
                                        const CubitVector& sweep_axis,
                                        double angle,
                                        int steps,
                                        double draft_angle,
                                        int draft_type,
                                        CubitBoolean switchside,
                                        CubitBoolean make_solid,
                                        CubitBoolean rigid)
{
  DLIList<Body*> body_list;
  DLIList<GeometryEntity*> geom_list;
  GeometryModifyEngine* gePtr1 = 0;
  CubitBoolean change_newids;
  if (!sweep_setup("rotational", ref_ent_list, body_list, gePtr1,
                   change_newids, geom_list))
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled())
  {
    DLIList<RefEdge*> edge_list;
    DLIList<RefFace*> face_list;

    CAST_LIST( ref_ent_list, edge_list, RefEdge );
    CAST_LIST( ref_ent_list, face_list, RefFace );

    //Edges aren't consumed, so there's nothing to save out
    if( edge_list.size() )
      CubitUndo::save_state();
    else
     //Faces will get consumed so you have to save out original entities
      CubitUndo::save_state_with_cubit_file( face_list );
  }

  DLIList<BodySM*> result_list;
  CubitStatus status = gePtr1->sweep_rotational( geom_list,
                                                 result_list,
                                                 point,
                                                 sweep_axis,
                                                 angle,
                                                 steps,
                                                 draft_angle,
                                                 draft_type,
                                                 switchside,
                                                 make_solid,
                                                 rigid);
  DLIList<Body*> output_body_list;
  if (!sweep_finish("rotational", body_list, result_list, output_body_list, change_newids))
    status = CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled())
  {
    if( status == CUBIT_FAILURE )
      CubitUndo::remove_last_undo();
    else
      CubitUndo::note_result_bodies( output_body_list );
  }

  return status;
}

CubitStatus GeometryModifyTool::sweep_translational(
                                         DLIList<RefEntity*>& ref_ent_list,
                                         const CubitVector& sweep_vector,
                                         double draft_angle,
                                         int draft_type,
                                         CubitBoolean switchside,
                                         CubitBoolean rigid)
{
  DLIList<Body*> body_list;
  DLIList<GeometryEntity*> geom_list;
  GeometryModifyEngine* gePtr1 = 0;
  CubitBoolean change_newids;
  if (!sweep_setup("translational", ref_ent_list, body_list, gePtr1,
                   change_newids, geom_list))
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled())
  {
    DLIList<RefEdge*> edge_list;
    DLIList<RefFace*> face_list;

    CAST_LIST( ref_ent_list, edge_list, RefEdge );
    CAST_LIST( ref_ent_list, face_list, RefFace );

    //Edges aren't consumed, so there's nothing to save out
    if( edge_list.size() )
      CubitUndo::save_state();
    else
     //Faces will get consumed so you have to save out original entities
      CubitUndo::save_state_with_cubit_file( face_list );
  }

  DLIList<BodySM*> result_list;
  CubitStatus status = gePtr1->
    sweep_translational( geom_list,
                         result_list,
                         sweep_vector,
                         draft_angle,
                         draft_type,
                         switchside,
                         rigid);

  DLIList<Body*> output_body_list;
  if (!sweep_finish("translational", body_list, result_list, output_body_list, change_newids))
    status = CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled())
  {
    if( status == CUBIT_FAILURE )
      CubitUndo::remove_last_undo();
    else
      CubitUndo::note_result_bodies( output_body_list );
  }

  return status;
}

//Author:: Jonathan Bugman
//Sept 10, 2006
CubitStatus GeometryModifyTool::sweep_curve_target(CubitPlane ref_plane,
                                             DLIList<RefEntity*>& ref_ent_list)
{
	double distance1;
	double distance2;
	double distance3;
	double temp_counter=0;
	CubitVector begin_point;
	CubitVector target_begin_point;
	CubitVector get_mid_point;
	CubitVector target_mid_point;
	CubitVector end_point;
	CubitVector target_end_point;
	DLIList<RefEdge*> edge_list;
	CAST_LIST(ref_ent_list, edge_list, RefEdge);
	CubitVector result;
	double max_distance_for_edge=0;
	CubitVector max_result;
	//make sure that there are actually edges in list
	if(edge_list.size() > 0)
	{
		//this part of the code steps through all edges that could be
		//selected by user and finds the largest distance of all edges from three
		//points: midpoint, endpoint, and beginpoint
		for (int i=0;i<edge_list.size();i++)
		{
			//find the midpoint vector of the edge
			edge_list[i]->mid_point(get_mid_point);

			//Project the midpoint onto the specified plane
			target_mid_point = ref_plane.project(get_mid_point);

			//Calculate the distance between the mid_point, and target_point
			distance1 = target_mid_point.distance_between(get_mid_point);

			//the next two blocks are just copies of the above three steps
			double fraction_along_curve = 1;
			edge_list[i]->position_from_fraction(fraction_along_curve, end_point);
			target_end_point = ref_plane.project(end_point);
			distance2 = target_end_point.distance_between(end_point);

			fraction_along_curve = 0;
			edge_list[i]->position_from_fraction(fraction_along_curve, begin_point);
			target_begin_point = ref_plane.project(begin_point);
			distance3 = target_begin_point.distance_between(begin_point);

			//see which of the three distances is greater for that edge
			//and compare its distance to the other edges that have already been calculated
			//saving the vector of the largest distance
			if (distance1>distance2 && distance1>distance3)
			{
				result = 2*(target_mid_point - get_mid_point);
				max_distance_for_edge=distance1;
				if (max_distance_for_edge>temp_counter)
				{
					temp_counter=max_distance_for_edge;
					max_result=result;
				}
			}
			else if (distance2>distance3)
			{
				result = 2*(target_end_point - end_point);
				max_distance_for_edge=distance2;
				if (max_distance_for_edge>temp_counter)
				{
					temp_counter=max_distance_for_edge;
					max_result=result;
				}
			}
			else
			{
				result = 2*(target_begin_point - begin_point);
				max_distance_for_edge=distance3;
				if (max_distance_for_edge>temp_counter)
				{
					temp_counter=max_distance_for_edge;
					max_result=result;
				}
			}
		}

		//using the 'facet edges' defined by the drawing geometry take the leading
		//and trailing edge projection vectors of this small edge.  With that,
		//calc the dot product and if negative, the vectors are opposite of each
		//other meaning the curve travels through the plane and the sweep function will fail
		for (int ii=0;ii<edge_list.size();ii++)
		{
			Curve *facet_curve;
			facet_curve=edge_list[ii]->get_curve_ptr();
			CubitStatus response;
			GMem g_mem;

			//get number of points and their locations
			//on the curve as defined by the drawing geometry algorithm
			response = facet_curve->get_geometry_query_engine()->
				get_graphics( facet_curve, &g_mem );
			int num_points = g_mem.pointListCount;

			if (response==CUBIT_FAILURE || num_points == 0)
			{
				PRINT_WARNING("Facet Curve calling function failed\n" );
			}

			GPoint *point_data = g_mem.point_list();

			for (int jj=0; jj<= (num_points-1); jj++)
			{
				//get first points vectors
				CubitVector point_1;
				point_1.x(point_data[jj].x);
				point_1.y(point_data[jj].y);
				point_1.z(point_data[jj].z);
				//get second points vectors
				CubitVector point_2;
				point_2.x(point_data[jj+1].x);
				point_2.y(point_data[jj+1].y);
				point_2.z(point_data[jj+1].z);
				//project the two points onto target plane
				CubitVector target_point_1;
				target_point_1 = ref_plane.project(point_1);
				CubitVector target_point_2;
				target_point_2 = ref_plane.project(point_2);
				//calc vector from point on curve to point on surface
				CubitVector vec_1 = point_1 - target_point_1;
				CubitVector vec_2 = point_2 - target_point_2;
				//make them unit vectors
				vec_1.normalize();
				vec_2.normalize();
				//calculate dot product
				double dot = vec_1.x()*vec_2.x()+vec_1.y()*vec_2.y()+vec_1.z()*vec_2.z();
				//check to see if dot product sign is zero
				//the and statement checks for if the first or last vertex of the edge
				//which may be sitting on the surface, this is alright
				//and should still be able to sweep
				if (dot<0 && (jj+1!=num_points || jj==0))
				{
					PRINT_ERROR( "The edge is traveling through the plane\n" );
					PRINT_ERROR( "and thus forces the sweep direction to switch\n" );
					PRINT_ERROR( "direction. Try splitting the edge.\n" );
					return CUBIT_FAILURE;
				}
			}
		}
		DLIList<BodySM*> webcut_results_list;
		DLIList<Body*> body_list;
		DLIList<BodySM*> sweep_result_list;
		DLIList<GeometryEntity*> geom_list;
		GeometryModifyEngine* gePtr1 = 0;
		CubitBoolean change_newids;

		if (!sweep_setup("translational", ref_ent_list , body_list, gePtr1,
			change_newids, geom_list))
			return CUBIT_FAILURE;

		//below block is default settings to be fed into sweep_translational
		CubitBoolean rigid = CUBIT_FALSE;
		CubitBoolean switchside = CUBIT_FALSE;
		int draft_type = 0;
		double draft_angle=0.0;

		//sweep the curve through and hopefully past the target surface
		CubitStatus status = gePtr1->sweep_translational( geom_list,sweep_result_list,
			max_result,draft_angle, draft_type,switchside,rigid);

		if (status == 0)
		{
			PRINT_ERROR( "Sweep operation failed!\n" );
			//delete sweep_result_list memory
			gePtr1->get_gqe()->delete_solid_model_entities(sweep_result_list);
			return CUBIT_FAILURE;
		}

		//below lines are used to define the ref_plane with three points for webcut input
		//find normal vector to user specified plane
		CubitVector vec1 = ref_plane.normal();
		CubitVector vec2;
		CubitVector vec3;
		//get orthogonal vectors of the normal vector
		vec1.orthogonal_vectors(vec2, vec3);
		//place the orthogonal vectors at a point on the plane as opposed to the origin
		vec2=vec2+target_mid_point;
		vec3=vec3+target_mid_point;

                DLIList<BodySM*> neighbor_imprint_list;
		//do a webcut with a plane created from the three projected points above
		CubitStatus status2 = gePtr1->webcut(sweep_result_list,target_mid_point,
			vec2,vec3, neighbor_imprint_list, webcut_results_list);

		if (status2 == 0)
		{
			PRINT_ERROR( "Sweep operation worked; however, webcut operation failed.\n" );
			//delete memory since it failed
			gePtr1->get_gqe()->delete_solid_model_entities(sweep_result_list);
			return CUBIT_FAILURE;
		}

		if (webcut_results_list.size()==0)
		{
			PRINT_ERROR( "Number of bodies from webcut is zero, unable to perform rest of sweep operation\n" );
			//delete memory since it failed
			gePtr1->get_gqe()->delete_solid_model_entities(sweep_result_list);
			return CUBIT_FAILURE;
		}

		CubitVector created_curve_mid;
		CubitVector user_mid;
		DLIList<BodySM*> keep_bodies_list;
		DLIList<Curve*> curve_list;

		//generate a list of curves from the webcut_results_list
		for (int counter = 0;counter < webcut_results_list.size(); counter++)
		{
			BodySM* OO = webcut_results_list[counter];
			OO->curves(curve_list);
		}

		//step through each of the user specified edges
		for (int edge_counter = 0; edge_counter < edge_list.size(); edge_counter++)
		{
			//find the midpoint of a user specified edge
			edge_list[edge_counter]->position_from_fraction(.5,user_mid);
			double min_dist=DBL_MAX;
			double distance;
			Curve* closest_curve=0;

			//step through all of the curves
			for (int counter2 = 0; counter2 < curve_list.size(); counter2++)
			{
				//find the midpoint of the created curve
				curve_list[counter2]->position_from_fraction(.5,created_curve_mid);
				//calculate the distance between the two midpoints
				distance=created_curve_mid.distance_between(user_mid);

				//find the minimum distance between all the curves
				if (distance<min_dist)
				{
					//reset the min_distance
					min_dist=distance;
					//keep track of which curve is associated with the shortest distance
					closest_curve=curve_list[counter2];
				}
			}
			//find the parent body of the minimum distance curve
			BodySM* body_sm = closest_curve->bodysm();
			//append that body to a keep list
			keep_bodies_list.append(body_sm);
		}
		//delete bodies which repeat (useful only when in granite engine)
		keep_bodies_list.uniquify_unordered();
		//subtract the keep_list from the webcut_results_list
		webcut_results_list -=keep_bodies_list;
		//delete webcut_results_list memory
		gePtr1->get_gqe()->delete_solid_model_entities(webcut_results_list);

                if( CubitUndo::get_undo_enabled() )
                  CubitUndo::save_state();

		//builds ref bodies
    DLIList<Body*> output_body_list;
		if (!sweep_finish("translational", body_list, keep_bodies_list, output_body_list, change_newids))
		{
			gePtr1->get_gqe()->delete_solid_model_entities(keep_bodies_list);
			status = CUBIT_FAILURE;
		}

    if( CubitUndo::get_undo_enabled())
    {
      if( status == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else
        CubitUndo::note_result_bodies( output_body_list );
    }
	}
	else
	{
		PRINT_ERROR( "No edge(s) found - sweep creation failed.\n" );
		return CUBIT_FAILURE;
	}

	return CUBIT_SUCCESS;
}
// Author: Derek Quam
// Dec. 1, 2008
CubitStatus GeometryModifyTool::sweep_surface_target(RefFace *face,
                                                     Body *target_body,
                                                     CubitVector distance,
                                                     CubitPlane stop_plane,
                                                     double magnitude)
{
  //  Declare local variables
  Surface *source_surf;
  bool set_dist = true, set_plane = true, set_mag = true;
  CubitVector direction;
  CubitStatus status;
  DLIList<Body*> body_list;
  DLIList<GeometryEntity*> geom_list;
  GeometryModifyEngine* gePtr1 = 0;
  CubitBoolean change_newids;
  DLIList<RefEntity*> ref_ent_list;
  ref_ent_list.append(face);

  // Check to make sure the target body is a solid
  if (target_body->is_sheet_body())
  {
    PRINT_ERROR("Target body must be a solid body.\n");
    return CUBIT_FAILURE;
  }

  // Set up the sweep
  if (!sweep_setup("target", ref_ent_list, body_list, gePtr1, change_newids, geom_list))
    return CUBIT_FAILURE;

  // Get the Surface * from the RefFace *
  source_surf = face->get_surface_ptr();

  // Check if the direction vector and stop plane were specified
  if (distance.length() < 0.0001)
    set_dist = false;
  if (stop_plane.normal().length() < 0.0001)
    set_plane = false;
  if (magnitude < 0.0001)
    set_mag = false;

  // Calculate the direction of the sweep
  if (!set_dist)
    direction = face->normal_at(face->center_point());
  else
    direction = distance;
  direction.normalize();

  double length = 0.0;
  if (!set_mag && !set_plane)
  {
    CubitVector center_body = target_body->center_point();
    CubitVector center_face = face->center_point();
    length = center_face.distance_between(center_body);
  }
  else if (set_plane)
  {
    length = stop_plane.intersect(face->center_point(), direction).distance_between(face->center_point());
  }
  else
  {
    length = magnitude;
  }
  if (set_mag && length > magnitude)
    length = magnitude;
  direction *= length;

  DLIList<BodySM*> new_bodies;
  DLIList<GeometryEntity*> src_ents;
  src_ents.append( source_surf );
  status = gePtr1->sweep_translational( src_ents, new_bodies, direction, 0, 0, false, false, 0, target_body->get_body_sm_ptr());

  if (status != CUBIT_SUCCESS)
    return status;

  // Make all the new bodies
  DLIList<Body*> output_body_list;
  if (!sweep_finish("target", body_list, new_bodies, output_body_list, change_newids))
    status = CUBIT_FAILURE;
    
  
  /*
  for (int i = 0; i < new_bodies.size(); i++)
    GeometryQueryTool::instance()->make_Body(new_bodies.get_and_step());
    */
    

  return CUBIT_SUCCESS;
}

// Author: Andrew Rout and Derek Quam
// Nov. 14, 2008
CubitStatus GeometryModifyTool::sweep_curve_target(DLIList<RefEdge*>& edge_list,
                                                   Body *target_body,
                                                   DLIList<Body*> &out_bodies,
                                                   CubitVector distance,
                                                   CubitPlane stop_plane,
                                                   bool unite)
{
  DLIList<BodySM*> new_bodies;
  DLIList<Curve*> curve_list;
  bool set_dist = true, set_plane = true;
  double larDist = 1.0;
  CubitVector dir = distance;
  CubitStatus status;

  // Check to make sure the target body is a sheetbody
  if (!target_body->is_sheet_body())
  {
    PRINT_ERROR("Target body must be a sheet body.\n");
    return CUBIT_FAILURE;
  }

  // Get the Curve *'s from the RefEdge *'s
  for (int i = 0; i < edge_list.size(); i++)
  {
    edge_list[i]->get_curve_ptr()->set_saved_id(edge_list[i]->id());
    curve_list.append(edge_list[i]->get_curve_ptr());
  }
  // Check if the direction vector and stop plane were specified
  if (distance.length() < 0.0001)
    set_dist = false;
  if (stop_plane.normal().length() < 0.0001)
    set_plane = false;

  // Check inputs
  if (!set_plane && !set_dist)
  {
    PRINT_ERROR("User must specify a stop plane, a direction, or both.\n");
    return CUBIT_FAILURE;
  }

  // Calculate the direction vector
  double begDist, midDist, endDist;
  for (int i = 0; i < edge_list.size(); i++)
  {
    CubitVector beg, mid, end, begP, midP, endP;
    RefEdge *temp = edge_list.get_and_step();

    // Retrieve the beginning, middle, and end coordinates of the edge
    beg = temp->start_coordinates();
    temp->mid_point(mid);
    end = temp->end_coordinates();

    if (set_plane)
    {
      // Project the start, mid, and end point onto the stop plane
      begP = stop_plane.project(beg);
      midP = stop_plane.project(mid);
      endP = stop_plane.project(end);

      // Calculate the distance between the points
      begDist = beg.distance_between(begP);
      midDist = mid.distance_between(midP);
      endDist = end.distance_between(endP);
    }
    else  // No stop plane specified
    {
      begDist = beg.distance_between(target_body->center_point());
      midDist = mid.distance_between(target_body->center_point());
      endDist = end.distance_between(target_body->center_point());
    }

    // Find the largest distance
    if (begDist > larDist)
      larDist = begDist;
    if (midDist > larDist)
      larDist = midDist;
    if (endDist > larDist)
      larDist = endDist;

    // Make sure the plane normal is pointing the right way
    if (set_plane)
    {
      double planeNorm = stop_plane.normal().interior_angle(beg-begP);
      if (planeNorm <= 90 || planeNorm >= 270)
        stop_plane.reverse();
    }

    if (!set_dist)
      dir += (midP-mid);

  } // End for loop
  if (!set_dist)
    dir /= edge_list.size();

  // Unitize the direction vector
  dir /= dir.length();

  // Add the magnitude to the direction vector and check for intersection with the stop plane
  dir *= larDist;
  if (set_plane)
  {
    double angle = dir.interior_angle(stop_plane.normal());
    if (angle >= 90 && angle <= 270)
      PRINT_WARNING("Direction vector does not intersect stop plane!\n");
  }

  // Call the geometry function
  DLIList<GeometryEntity*> source_ents;
  CAST_LIST_TO_PARENT(curve_list, source_ents);
  status = get_gme()->sweep_translational(source_ents, new_bodies, dir, 0, 0, false, false, 0, target_body->get_body_sm_ptr());
  if (status != CUBIT_SUCCESS)
    return CUBIT_FAILURE;
    
  if (unite) {
    DLIList<BodySM*> tmp_bodies(new_bodies);
    new_bodies.clean_out();
    status = get_gme()->unite( tmp_bodies, new_bodies );
    if (status != CUBIT_SUCCESS)
      return CUBIT_FAILURE;
  }

  // Make all the new bodies
  for (int i = 0; i < new_bodies.size(); i++)
  {
    out_bodies.append(GeometryQueryTool::instance()->make_Body(new_bodies.get_and_step()));
    PRINT_INFO("Created volume in body %d\n", out_bodies[i]->id());
  }

  return CUBIT_SUCCESS;
}

//Author: Jonathan Bugman
//Sept 10, 2006
CubitStatus GeometryModifyTool::sweep_surface_target( CubitPlane ref_plane,
                                    DLIList<RefEntity*>& ref_ent_list)
{
	DLIList<RefFace*> surface_list;
	CAST_LIST(ref_ent_list, surface_list, RefFace);
	double distance1;
	double distance2;
	double distance3;
	double temp_counter=0;
	CubitVector begin_point;
	CubitVector target_begin_point;
	CubitVector get_mid_point;
	CubitVector target_mid_point;
	CubitVector end_point;
	CubitVector target_end_point;
	CubitVector result;
	double max_distance_for_edge=0;
	CubitVector max_result;
	DLIList<RefEdge*> surf_edge_list;

	//make sure that only one surface has been selected
	if(surface_list.size() == 0)
	{
		PRINT_ERROR( "No surface found - sweep surface to target failed.\n" );
		return CUBIT_FAILURE;
	}

	//make sure that there are actually surfaces in list
	else if(surface_list.size() > 1)
	{
		PRINT_ERROR( "You can only use this operation \n with one surface selected at a time\n" );
		return CUBIT_FAILURE;
	}
	else
	{
		//this part of the code steps through all edges that could be
		//selected by user and finds the largest distance of all edges from three
		//points: midpoint, endpoint, and beginpoint
		for (int i=0;i<surface_list.size();i++)
		{
			surface_list[i]->ref_edges(surf_edge_list);
			for (int j=0;j<surf_edge_list.size();j++)
			{
				//get midpoint of edge on surface
				surf_edge_list[j]->mid_point(get_mid_point);
				//Project the midpoint of each surface edge onto the specified plane
				target_mid_point = ref_plane.project(get_mid_point);

				//Calculate the distance between the mid_point, and target_point
				distance1 = target_mid_point.distance_between(get_mid_point);

				//the next two blocks are just copies of the above three steps

				int fraction_along_curve = 1;
				surf_edge_list[j]->position_from_fraction(fraction_along_curve, end_point);
				target_end_point = ref_plane.project(end_point);
				distance2 = target_end_point.distance_between(end_point);

				fraction_along_curve = 0;
				surf_edge_list[j]->position_from_fraction(fraction_along_curve, begin_point);
				target_begin_point = ref_plane.project(begin_point);
				distance3 = target_begin_point.distance_between(begin_point);

				//see which of the three distances is greater of the edge
				//and compare its distance to the other edges
				if (distance1>distance2 && distance1>distance3)
				{
					result = 2*(target_mid_point - get_mid_point);
					max_distance_for_edge=distance1;
					if (max_distance_for_edge>temp_counter)
					{
						temp_counter = max_distance_for_edge;
						max_result=result;

					}
				}
				else if (distance2>distance3)
				{
					result = 2*(target_end_point - end_point);
					max_distance_for_edge=distance2;
					if (max_distance_for_edge>temp_counter)
					{
						temp_counter = max_distance_for_edge;
						max_result=result;

					}
				}
				else
				{
					result = 2*(target_begin_point - begin_point);
					max_distance_for_edge=distance3;
					if (max_distance_for_edge>temp_counter)
					{
						temp_counter = max_distance_for_edge;
						max_result=result;

					}
				}
				//just checking to make sure the user didn't specify a surface as both
				//a target surface and sweeping surface which would result in a CUBIT crash
				if (i==surface_list.size()-1 && j==surf_edge_list.size()-1 &&
					max_distance_for_edge < 0.000000000000001)
				{
					PRINT_ERROR( "The sweep distance is less than the geometry tolerance\n" );
					PRINT_ERROR( "This may be caused by selecting the same\n");
					PRINT_ERROR( "sweep surface and target surface\n" );
					return CUBIT_FAILURE;
				}
			}
		}
		//using the facet edges defined by the drawing geometry code; take the leading
		//and trailing edge vectors of each facet edge;  With that, calc the dot product
		//and if negative (meaning the vector directions switched) the curve travels
		//through the plane and the sweep function will fail
		surface_list.reset();
		surf_edge_list.clean_out();
		for (int kk=0;kk<surface_list.size();kk++)
		{
			surface_list[kk]->ref_edges(surf_edge_list);
			for (int ii=0;ii<surf_edge_list.size();ii++)
			{
				Curve *facet_curve;
				facet_curve=surf_edge_list[ii]->get_curve_ptr();
				CubitStatus response;
				GMem g_mem;

				//get number of points and their locations
				//on the curve as defined by the drawing geometry algorithm
				response = facet_curve->get_geometry_query_engine()->
					get_graphics( facet_curve, &g_mem );

        int num_points = g_mem.pointListCount;

				if (response==CUBIT_FAILURE || num_points == 0)
				{
					PRINT_WARNING("Unable to preview a curve\n" );
				}

				GPoint *point_data = g_mem.point_list();

				for (int jj=0; jj<= (num_points-1); jj++)
				{
					//get first points vectors
					CubitVector point_1;
					point_1.x(point_data[jj].x);
					point_1.y(point_data[jj].y);
					point_1.z(point_data[jj].z);
					//get second points vectors
					CubitVector point_2;
					point_2.x(point_data[jj+1].x);
					point_2.y(point_data[jj+1].y);
					point_2.z(point_data[jj+1].z);
					//project the two points onto target plane
					CubitVector target_point_1;
					target_point_1 = ref_plane.project(point_1);
					CubitVector target_point_2;
					target_point_2 = ref_plane.project(point_2);
					//calc vector from point on curve to point on surface
					CubitVector vec_1 = point_1 - target_point_1;
					CubitVector vec_2 = point_2 - target_point_2;
					//make them unit vectors
					vec_1.normalize();
					vec_2.normalize();
					//double dot = DotProduct(vec_1,vec_2);
					//calculate dot product
					double dot = vec_1.x()*vec_2.x()+vec_1.y()*vec_2.y()+vec_1.z()*vec_2.z();

					//check to see if dot product sign is zero
					//the and statement checks for if the first or last vertex of the edge
					//which may be sitting on the surface, this is alright
					//and should still be able to sweep
					if (dot<0 && (jj+1!=num_points || jj==0))
					{
						PRINT_ERROR( "The surface is traveling through the plane\n" );
						PRINT_ERROR( "and thus forces the sweep direction to switch\n" );
						PRINT_ERROR( "direction. Try splitting the surface.\n" );
						return CUBIT_FAILURE;
					}
				}
			}
		}

		DLIList<Body*> body_list;
		DLIList<GeometryEntity*> geom_list;
		GeometryModifyEngine* gePtr1 = 0;
		CubitBoolean change_newids;

		if (!sweep_setup("translational", ref_ent_list, body_list, gePtr1,
			change_newids, geom_list))
			return CUBIT_FAILURE;

    if( CubitUndo::get_undo_enabled())
      //Faces will get consumed so you have to save out original entities
      CubitUndo::save_state_with_cubit_file( surface_list );

		DLIList<BodySM*> sweep_result_list;

		//below block is default settings to be fed into sweep_translational
		CubitBoolean rigid = CUBIT_FALSE;
		CubitBoolean switchside = CUBIT_FALSE;
		int draft_type = 0;
		double draft_angle=0.0;

		//find normal vector to user specified plane
		CubitVector vec1 = ref_plane.normal();
		CubitVector vec2;
		CubitVector vec3;
		//get orthogonal vectors of the normal vector
		vec1.orthogonal_vectors(vec2, vec3);
		//place the orthogonal vectors at a point on the plane as opposed to the origin
		vec2=vec2+target_mid_point;
		vec3=vec3+target_mid_point;
		DLIList<BodySM*> webcut_results_list;
		//get a point on the first user specified surface before it gets consumed in the sweep_translational function
		CubitVector mid_point_surface;
		mid_point_surface = surface_list[0]->center_point();

		//sweep the curve down through and hopefully past the target surface
		CubitStatus status = gePtr1->sweep_translational( geom_list,sweep_result_list,
			max_result,draft_angle, draft_type,switchside,rigid);

		if (status == 0)
		{
			//If in here, sweep_translational failed so delete result_list and
			//print an error to the screen for the user
			gePtr1->get_gqe()->delete_solid_model_entities(sweep_result_list);
			PRINT_ERROR( "Error occured in the sweep operation.\n" );

                        if( CubitUndo::get_undo_enabled())
                          CubitUndo::remove_last_undo();
			return CUBIT_FAILURE;
		}

                DLIList<BodySM*> neighbor_imprint_list;
		//do a webcut with a plane created from the three projected points
		CubitStatus status2 = gePtr1->webcut(sweep_result_list,target_mid_point,
			vec2,vec3,neighbor_imprint_list, webcut_results_list);

		if (status2 == 0)
		{
			//If in here, webcut operation failed so delete result_list and
			//print an error to the screen for the user
			gePtr1->get_gqe()->delete_solid_model_entities(sweep_result_list);
			gePtr1->get_gqe()->delete_solid_model_entities(webcut_results_list);
			PRINT_ERROR( "Error occured in the webcut operation.\n" );

                        if( CubitUndo::get_undo_enabled())
                          CubitUndo::remove_last_undo();
			return CUBIT_FAILURE;
		}

		DLIList<BodySM*> keep_bodies_list = webcut_results_list;
		int bodies_deleted_counter=0;
		CubitVector target_mid_point_volume;
		CubitVector mid_point_volume;
		CubitVector target_mid_point_surface;
		double volume_size;

		//project the mid_point of user specified surface onto target plane
		target_mid_point_surface = ref_plane.project(mid_point_surface);
		//calculate vector between the two points and then normalize
		CubitVector vec_1 = mid_point_surface - target_mid_point_surface;
		vec_1.normalize();

		//step through the
		for (int counter=0;counter<webcut_results_list.size();counter++)
		{
			//find the geometric midpoint of the body and project that point on the target plane
			webcut_results_list[counter]->mass_properties(mid_point_volume,volume_size);
			target_mid_point_volume = ref_plane.project(mid_point_volume);
			//generate a vector between the two points and then normalize
			CubitVector vec_2 = mid_point_volume - target_mid_point_volume;
			vec_2.normalize();

			//calculate the dot product
			double dot = vec_1.x()*vec_2.x()+vec_1.y()*vec_2.y()+vec_1.z()*vec_2.z();

			//if a negative dot product delete it because it is on the opposite side of the target plane
			if (dot < 0)
			{
				keep_bodies_list.remove(webcut_results_list[counter]);
				bodies_deleted_counter = bodies_deleted_counter + 1;
			}

		}

		//test to see if all bodies have been deleted and if so let the user know
		if (bodies_deleted_counter == webcut_results_list.size())
		{
			PRINT_ERROR( "All sweeped surfaces deleted - sweep_target failed.\n" );
			PRINT_ERROR( "This may be due to granite engine limitations and/or\n" );
			PRINT_ERROR( "angle between curve and target surface\n" );

                        if( CubitUndo::get_undo_enabled())
                          CubitUndo::remove_last_undo();
			return CUBIT_FAILURE;
		}

		//delete webcut_results_list since it is no longer of use
		webcut_results_list -= keep_bodies_list;
		gePtr1->get_gqe()->delete_solid_model_entities(webcut_results_list);

		//builds ref bodies
    DLIList<Body*> output_body_list;
		if (!sweep_finish("translational", body_list, keep_bodies_list, output_body_list, change_newids))
			status = CUBIT_FAILURE;

    if( CubitUndo::get_undo_enabled())
    {
      if( status == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else
        CubitUndo::note_result_bodies( output_body_list );
    }

	}
	return CUBIT_SUCCESS;
}


CubitStatus GeometryModifyTool::sweep_perpendicular( DLIList<RefEntity*>& ref_ent_list,
                                                     double distance,
                                                     double draft_angle,
                                                     int draft_type,
                                                     CubitBoolean switchside,
                                                     CubitBoolean rigid)
{
  DLIList<Body*> body_list;
  DLIList<GeometryEntity*> geom_list;
  GeometryModifyEngine* gePtr1 = 0;
  CubitBoolean change_newids;
  if (!sweep_setup("translational", ref_ent_list, body_list, gePtr1,
                   change_newids, geom_list))
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled())
  {
    DLIList<RefEdge*> edge_list;
    DLIList<RefFace*> face_list;

    CAST_LIST( ref_ent_list, edge_list, RefEdge );
    CAST_LIST( ref_ent_list, face_list, RefFace );

    //Edges aren't consumed, so there's nothing to save out
    if( edge_list.size() )
      CubitUndo::save_state();
    else
     //Faces will get consumed so you have to save out original entities
      CubitUndo::save_state_with_cubit_file( face_list );
  }

  DLIList<BodySM*> result_list;
  CubitStatus status = gePtr1->
    sweep_perpendicular( geom_list,
                         result_list,
                         distance,
                         draft_angle,
                         draft_type,
                         switchside,
                         rigid);

  DLIList<Body*> output_body_list;
  if (!sweep_finish("perpendicular", body_list, result_list, output_body_list, change_newids))
    status = CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled())
  {
    if( status == CUBIT_FAILURE )
      CubitUndo::remove_last_undo();
    else
      CubitUndo::note_result_bodies( output_body_list );
  }

  return status;
}
CubitStatus GeometryModifyTool::sweep_along_curve(DLIList<RefEntity*>& ref_ent_list,
                                                  DLIList<RefEdge*>& ref_edge_list,
                                                  double draft_angle,
                                                  int draft_type,
                                                  CubitBoolean rigid)

{
  DLIList<GeometryEntity*> geom_list(ref_ent_list.size());
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> body_list(ref_ent_list.size());
  GeometryModifyEngine* engine_ptr = 0;
  CubitBoolean changed_new_ids = CUBIT_FALSE;
  CubitStatus status = sweep_setup( "along_curve",
                                    ref_ent_list,
                                    body_list,
                                    engine_ptr,
                                    changed_new_ids,
                                    geom_list,
                                    &ref_edge_list,
                                    &curve_list );
  if (status != CUBIT_SUCCESS)
  return status;

  if( CubitUndo::get_undo_enabled())
  {
    DLIList<RefEdge*> edge_list;
    DLIList<RefFace*> face_list;

    CAST_LIST( ref_ent_list, edge_list, RefEdge );
    CAST_LIST( ref_ent_list, face_list, RefFace );

    //Edges aren't consumed, so there's nothing to save out
    if( edge_list.size() )
      CubitUndo::save_state();
    else
      //Faces will get consumed so you have to save out original entities
      CubitUndo::save_state_with_cubit_file( face_list );
  }

  DLIList<BodySM*> result_list;
  status = engine_ptr->sweep_along_curve( geom_list,
                                          result_list,
                                          curve_list,
                                          draft_angle,
                                          draft_type,
                                          rigid);

  DLIList<Body*> output_body_list;
  if (!sweep_finish("along_curve", body_list, result_list, output_body_list, changed_new_ids))
  status = CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled())
  {
    if( status == CUBIT_FAILURE )
      CubitUndo::remove_last_undo();
    else
      CubitUndo::note_result_bodies( output_body_list );
  }

  return status;
}

void GeometryModifyTool::initialize_settings() {


  SettingHandler::instance()->add_setting("Group Imprint",
					  GeometryModifyTool::set_group_imprint,
					  GeometryModifyTool::get_group_imprint);

  SettingHandler::instance()->add_setting("NonRegImprint",
					  GeometryModifyTool::set_all_edges_imprint,
					  GeometryModifyTool::get_all_edges_imprint);

  SettingHandler::instance()->add_setting("New Ids",
					  GeometryModifyTool::set_new_ids,
					  GeometryModifyTool::get_new_ids);

  SettingHandler::instance()->add_setting("Separate After Webcut",
					  GeometryModifyTool::set_sep_after_webcut_setting,
					  GeometryModifyTool::get_sep_after_webcut_setting);

  SettingHandler::instance()->add_setting("old names",
            GeometryModifyTool::set_old_names,
            GeometryModifyTool::get_old_names);

  SettingHandler::instance()->add_setting("Mesh Auto Delete",
					  GeometryModifyTool::set_mesh_autodelete,
					  GeometryModifyTool::get_mesh_autodelete);

  SettingHandler::instance()->add_setting("Mesh Auto Delete Cache",
					  GeometryModifyTool::set_mesh_autodelete_remesh,
					  GeometryModifyTool::is_mesh_autodelete_remesh);
}



// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********

//-------------------------------------------------------------------------
// Purpose       : Constructor of the GeometryModifyTool class.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------
GeometryModifyTool::GeometryModifyTool(GeometryModifyEngine*gme_ptr)
{
  if (gme_ptr) gmeList.append(gme_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : This functions webcuts a list of bodies through a plane.
//                 The newly created bodies are and merged depeding on the
//                 respective flags.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------

CubitStatus GeometryModifyTool::webcut_with_brick(
                                           DLIList<Body*>& webcut_body_list,
                                           const CubitVector &center,
                                           const CubitVector axes[3],
                                           const CubitVector &extension,
                                           DLIList<Body*> &results_list,
                                           DLIList<Body*> &neighboring_bodies,
                                           ImprintType imprint_type,
                                           CubitBoolean merge,
                                           CubitBoolean preview)
{
   // Make sure that entity creation is possible.  Allow at most
   // only one of the dimensions to be zero - in which case a planar
   // sheet is used to do the cutting.
   double width = 2.0*extension.x();
   double height = 2.0*extension.y();
   double depth = 2.0*extension.z();
   if ( width < 0.0 || height < 0.0 || depth < 0.0 )
   {
      PRINT_ERROR( "Cannot make a brick of size %f x %f x %f\n"
                   "     Negative dimensions are not allowed.\n",
                   width, height, depth );
      return CUBIT_FAILURE ;
   }

   BodySM *cutting_tool_ptr = NULL;
   CubitBoolean is_sheet_body = CUBIT_FALSE;
   CubitVector p1, p2, p3, p4;
   int wz = width < GEOMETRY_RESABS;
   int hz = height < GEOMETRY_RESABS;
   int dz = depth < GEOMETRY_RESABS;
   int num_zero_dim = wz + hz + dz;
   if( num_zero_dim > 0 )
   {
      if( num_zero_dim > 1 )
      {
         PRINT_ERROR( "Cannot make a sheet of size %f x %f x %f\n"
            "     At least two dimensions must be nonzero.\n",
            width, height, depth );
         return CUBIT_FAILURE ;
      }

      // Make a sheet body instead of a cuboid
      CubitVector sheet_axes[2];
      if( wz )
      {
         sheet_axes[0] = axes[1];
         sheet_axes[1] = axes[2];
         width = depth;
      }
      else if( hz )
      {
         sheet_axes[0] = axes[2];
         sheet_axes[1] = axes[0];
         height = depth;
      }
      else
      {
         sheet_axes[0] = axes[0];
         sheet_axes[1] = axes[1];
      }

      // Create the planar sheet to cut with
      // Get the corners of the sheet
      center.next_point( axes[0], width/2.0, p1 );
      p1.next_point( axes[1], -height/2.0, p1 );
      p1.next_point( axes[1], height, p2 );
      p2.next_point( axes[0], -width, p3 );
      p3.next_point( axes[1], -height, p4 );
   }

   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;  

   if( CubitUndo::get_undo_enabled() )
   {
     DLIList<Body*> bodies_to_save;
     bodies_to_save += webcut_body_list;
     bodies_to_save += neighboring_bodies;
     CubitUndo::save_state_with_cubit_file( bodies_to_save );
   }

   CubitStatus rval = CUBIT_SUCCESS;

   const int count = webcut_body_list.size();
   DLIList<BodySM*> result_sm_list;
   DLIList<Body*> body_list(webcut_body_list);
   DLIList<BodySM*> engine_body_sms(count);
   DLIList<Body*> engine_bodies(count);
   GeometryModifyEngine* gme = 0;

   while ( (gme = group_bodies_by_engine(body_list, engine_bodies, engine_body_sms)) )
   {
     DLIList<int> merged_surface_ids;
     DLIList<int> merged_curve_ids;
     DLIList<BodySM*> neighbor_imprint_list;
     if (!preview)
     {
       int i;
       for( i=neighboring_bodies.size(); i--; )
       {
         Body *neighbor_body = neighboring_bodies.get_and_step();
         BodySM *tmp_body = neighbor_body->get_body_sm_ptr(); 
         GeometryModifyEngine *neighbor_gme = get_engine( tmp_body );
         
         if( gme == neighbor_gme )
         {
           neighbor_imprint_list.append( tmp_body ); 
           engine_bodies.append( neighbor_body );
         }
       }

       do_attribute_setup();
       DLIList<BodySM*> bodies_to_modify;
       bodies_to_modify += engine_body_sms;
       bodies_to_modify += neighbor_imprint_list;
       push_vg_attributes_before_modify( bodies_to_modify );
       get_merged_curve_and_surface_ids( engine_bodies, merged_surface_ids, merged_curve_ids );
     }

      // Create the brick to cut with
      if (is_sheet_body)
	cutting_tool_ptr = gme->planar_sheet(p1,p2,p3,p4);
      else
        cutting_tool_ptr = gme->brick( center, axes, extension );
      if( cutting_tool_ptr == NULL )
         return CUBIT_FAILURE;

      CubitStatus status = gme->webcut (
        engine_body_sms, cutting_tool_ptr,
        neighbor_imprint_list,
        result_sm_list, imprint_type, preview );

      // Delete the BodySM that was created to be used as a tool
      gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

      // just continue the loop if previewing
      if (preview)
      {
        rval = status;
        continue;
      }

      restore_vg_after_modify(result_sm_list, engine_bodies, gme);
      remove_pushed_attributes(result_sm_list, engine_bodies);

      status = finish_webcut( engine_bodies, result_sm_list, merge, status, results_list,
                              &merged_surface_ids, &merged_curve_ids );

      if (!status)
        rval = CUBIT_FAILURE;

      engine_bodies.clean_out();
      engine_body_sms.clean_out();
      result_sm_list.clean_out();
   }

   if (!preview)
     do_attribute_cleanup();

   if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
   {
     if( results_list.size() ) //if there are new bodies...something succeeded
       CubitUndo::note_result_bodies( results_list );
     else
       CubitUndo::remove_last_undo();
   }

   return rval;
}

CubitStatus GeometryModifyTool::webcut_with_planar_sheet(
                                           DLIList<Body*>& webcut_body_list,
                                           const CubitVector &center,
                                           const CubitVector axes[2],
                                           double width, double height,
                                           DLIList<Body*> &results_list,
                                           DLIList<Body*> &neighboring_bodies,
                                           ImprintType imprint_type,
                                           CubitBoolean merge,
                                           CubitBoolean preview)
{
   if ( width <= GEOMETRY_RESABS || height <= GEOMETRY_RESABS )
   {
      PRINT_ERROR( "Cannot webcut with a sheet of size %f x %f\n"
                   "     Negative or zero dimensions are not allowed.\n",
                   width, height );
      return CUBIT_FAILURE ;
   }

   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   CubitStatus rval = CUBIT_SUCCESS;

  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<Body*> bodies_to_save;
    bodies_to_save += webcut_body_list;
    bodies_to_save += neighboring_bodies;
    CubitUndo::save_state_with_cubit_file( bodies_to_save );
  }

   const int count = webcut_body_list.size();
   DLIList<BodySM*> result_sm_list;
   DLIList<Body*> body_list(webcut_body_list);
   DLIList<BodySM*> engine_body_sms(count);
   DLIList<Body*> engine_bodies(count);
   GeometryModifyEngine* gme = 0;

   while ( (gme = group_bodies_by_engine(body_list, engine_bodies, engine_body_sms)) )
   {
     DLIList<int> merged_surface_ids;
     DLIList<int> merged_curve_ids;
     DLIList<BodySM*> neighbor_imprint_list;
     if (!preview)
     {
       int i;
       for( i=neighboring_bodies.size(); i--; )
       {
         Body *neighbor_body = neighboring_bodies.get_and_step();
         BodySM *tmp_body = neighbor_body->get_body_sm_ptr(); 
         GeometryModifyEngine *neighbor_gme = get_engine( tmp_body );
         
         if( gme == neighbor_gme )
         {
           neighbor_imprint_list.append( tmp_body ); 
           engine_bodies.append( neighbor_body );
         }
       }

       do_attribute_setup();
       DLIList<BodySM*> bodies_sms_to_modify;
       bodies_sms_to_modify += engine_body_sms;
       bodies_sms_to_modify += neighbor_imprint_list;
       push_vg_attributes_before_modify( bodies_sms_to_modify );
       get_merged_curve_and_surface_ids( engine_bodies, merged_surface_ids, merged_curve_ids );
     }

      // Create the planar sheet to cut with
      CubitVector p1, p2, p3, p4;

      // Get the corners of the sheet
      center.next_point( axes[0], width/2.0, p1 );
      p1.next_point( axes[1], -height/2.0, p1 );
      p1.next_point( axes[1], height, p2 );
      p2.next_point( axes[0], -width, p3 );
      p3.next_point( axes[1], -height, p4 );

      BodySM *cutting_tool_ptr = gme->planar_sheet(p1,p2,p3,p4);
      if( cutting_tool_ptr == NULL )
         return CUBIT_FAILURE;

      CubitStatus status = gme->webcut (
        engine_body_sms, cutting_tool_ptr,
        neighbor_imprint_list,
        result_sm_list, imprint_type, preview );
 
      // Delete the BodySM that was created to be used as a tool
      gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;

      // just continue the loop if previewing
      if (preview)
      {
        rval = status;
        continue;
      }

      restore_vg_after_modify(result_sm_list, engine_bodies, gme);
      remove_pushed_attributes(result_sm_list, engine_bodies);

      status = finish_webcut( engine_bodies, result_sm_list, merge, status, results_list,
                              &merged_surface_ids, &merged_curve_ids );
      if (!status)
        rval = CUBIT_FAILURE;

      engine_bodies.clean_out();
      engine_body_sms.clean_out();
      result_sm_list.clean_out();
   }

   if (!preview)
     do_attribute_cleanup();

  if( preview == CUBIT_FALSE && CubitUndo::get_undo_enabled() )
  {
    if( rval == CUBIT_SUCCESS )
      CubitUndo::note_result_bodies( results_list );
    else
      CubitUndo::remove_last_undo();
  }

   return rval;
}

CubitStatus GeometryModifyTool::webcut_with_plane(
                                    DLIList<Body*>& webcut_body_list,
                                    const CubitVector &vector1,
                                    const CubitVector &vector2,
                                    const CubitVector &vector3,
                                    DLIList<Body*>& results_list,
                                    DLIList<Body*> &neighboring_bodies,
                                    ImprintType imprint_type,
                                    CubitBoolean merge,
                                    CubitBoolean preview)
{
  if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
    return CUBIT_FAILURE;

  CubitStatus rval = CUBIT_SUCCESS;
  if (preview)
  {
    GeometryModifyTool::plane_preview(webcut_body_list, vector1, vector2, vector3);
    return rval;
  }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    DLIList<Body*> bodies_to_save;
    bodies_to_save += webcut_body_list;
    bodies_to_save += neighboring_bodies;
    CubitUndo::save_state_with_cubit_file( bodies_to_save );
  }

  const int count = webcut_body_list.size();
  DLIList<BodySM*> temp_sm_list(webcut_body_list.size());
  DLIList<BodySM*> result_sm_list;
  DLIList<Body*> body_list(webcut_body_list);
  DLIList<BodySM*> engine_body_sms(count);
  DLIList<Body*> engine_bodies(count);
  GeometryModifyEngine* gme = 0;

  // all preview stuff handled before this point
  if(!preview)
    do_attribute_setup();

  while ( (gme = group_bodies_by_engine(body_list, engine_bodies, engine_body_sms)) )
  {

    //get all the child entities that have been merged
    DLIList<int> merged_surface_ids;
    DLIList<int> merged_curve_ids;
    DLIList<BodySM*> neighbor_imprint_list;

    if (!preview)
    {
      int i;
      for( i=neighboring_bodies.size(); i--; )
      {
        Body *neighbor_body = neighboring_bodies.get_and_step();
        BodySM *tmp_body = neighbor_body->get_body_sm_ptr(); 
        GeometryModifyEngine *neighbor_gme = get_engine( tmp_body );
        
        if( gme == neighbor_gme )
        {
          neighbor_imprint_list.append( tmp_body ); 
          engine_bodies.append( neighbor_body );
        }
      }

      DLIList<BodySM*> bodies_sms_to_modify;
      bodies_sms_to_modify += engine_body_sms;
      bodies_sms_to_modify += neighbor_imprint_list;
      push_vg_attributes_before_modify( bodies_sms_to_modify );
      get_merged_curve_and_surface_ids( engine_bodies, merged_surface_ids, merged_curve_ids );
    }

    CubitStatus status = gme->webcut(engine_body_sms, vector1, vector2,
              vector3, neighbor_imprint_list, result_sm_list, imprint_type, preview );

    if ( status != CUBIT_FAILURE )
    {
      if(!preview)
      {
        restore_vg_after_modify(result_sm_list, engine_bodies, gme);
        remove_pushed_attributes(result_sm_list, engine_bodies);
      }
      status = finish_webcut( engine_bodies, result_sm_list, merge, status, results_list,
                              &merged_surface_ids, &merged_curve_ids );
    }
    else
    {
      if(!preview)
        remove_pushed_attributes(result_sm_list, engine_bodies);
    }

    engine_bodies.clean_out();
    engine_body_sms.clean_out();
    result_sm_list.clean_out();

    if ( status == CUBIT_FAILURE )
    {
      rval = CUBIT_FAILURE;
      break;
    }
  }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( rval == CUBIT_SUCCESS )
      CubitUndo::note_result_bodies( results_list );
    else
      CubitUndo::remove_last_undo();
  }

  if(!preview)
    do_attribute_cleanup();

  return rval;
}

void GeometryModifyTool::remove_pushed_attributes(DLIList<BodySM*> &new_sms,
                                                        DLIList<Body*> &old_bodies)
{
  DLIList<TopologyBridge*> old_bridges(old_bodies.size());
  DLIList<TopologyBridge*> new_bridges(new_sms.size());
  CAST_LIST(new_sms, new_bridges, TopologyBridge);

  // Get bridges for all of the old Bodies.
  int k;
  for(k = old_bodies.size(); k>0; k--)
  {
    Body *body = old_bodies.get_and_step();
    TopologyBridge *tb = body->bridge_manager()->topology_bridge();
    if(tb)
    {
      old_bridges.append(tb);
      DLIList<TopologyBridge*> bridge_list;
      bridge_list.append(tb);
      // Add any bodies with composites to the new_sms list so that
      // make_Body gets called on them.  This will make sure that the
      // virtual gets ref entities properly built.
      if(this->contains_composites(bridge_list))
      {
        BodySM *bsm = dynamic_cast<BodySM*>(tb);
        if(bsm)
          new_sms.append_unique(bsm);
      }
    }
  }

  // Make a list including all of the bridges passed in.
  DLIList<TopologyBridge*> all_bridges;
  all_bridges = new_bridges;
  for(k=old_bridges.size(); k--;)
    all_bridges.append_unique(old_bridges.get_and_step());

  // At this point we don't need any more attributes on the underlying
  // entities so make sure they are cleaned up.
  GeometryQueryTool::instance()->ige_remove_attributes( all_bridges );
}

CubitStatus GeometryModifyTool::restore_vg_after_modify(DLIList<BodySM*> &new_sms,
                                                        DLIList<Body*> &old_bodies,
                                                        GeometryModifyEngine *gme)
{
  DLIList<TopologyBridge*> old_bridges(old_bodies.size());
  DLIList<TopologyBridge*> new_bridges(new_sms.size());
  CAST_LIST(new_sms, new_bridges, TopologyBridge);

  // Get bridges for all of the old Bodies.
  int k;
  for(k = old_bodies.size(); k>0; k--)
  {
    Body *body = old_bodies.get_and_step();
    TopologyBridge *tb = body->bridge_manager()->topology_bridge();
    if(tb)
    {
      old_bridges.append(tb);
      DLIList<TopologyBridge*> bridge_list;
      bridge_list.append(tb);
      // Add any bodies with composites to the new_sms list so that
      // make_Body gets called on them.  This will make sure that the
      // virtual gets ref entities properly built.
      if(this->contains_composites(bridge_list))
      {
        BodySM *bsm = dynamic_cast<BodySM*>(tb);
        if(bsm)
          new_sms.append_unique(bsm);
      }
    }
  }

  // Make a list including all of the bridges passed in.
  DLIList<TopologyBridge*> all_bridges;
  all_bridges = new_bridges;
  for(k=old_bridges.size(); k--;)
    all_bridges.append_unique(old_bridges.get_and_step());

  DLIList<TopologyBridge*> tbs_to_check;
  if(gme)
    gme->get_possible_invalid_tbs(all_bridges, tbs_to_check);

  DLIList<Surface*> all_surfs;
  DLIList<Curve*> all_curves;
  DLIList<Point*> all_points;
  if(tbs_to_check.size() > 0)
  {
    for(k=tbs_to_check.size(); k--;)
    {
      TopologyBridge *tb = tbs_to_check.get_and_step();
      Surface *surf = dynamic_cast<Surface*>(tb);
      if(surf)
        all_surfs.append(surf);
      else
      {
        Curve *cur = dynamic_cast<Curve*>(tb);
        if(cur)
          all_curves.append(cur);
        else
        {
          Point *pt = dynamic_cast<Point*>(tb);
          if(pt)
            all_points.append(pt);
        }
      }
    }
  }

  // This function has been changed to blown away any virtual (really only doing
  // composites right now).  The virtual will rebuilt from the attributes stored
  // on the ACIS entities.
  GeometryQueryTool::instance()->ige_remove_modified(all_surfs, all_curves, all_points);

  //Restore virtual
  GeometryQueryTool::instance()->ige_import_geom( all_bridges );

  return CUBIT_SUCCESS;
}

GeometryModifyEngine* GeometryModifyTool::group_bodies_by_engine(
                                  DLIList<Body*>& remaining_bodies,
                                  DLIList<Body*>& engine_bodies,
                                  DLIList<BodySM*>& engine_body_sms ) const
{
  int i = remaining_bodies.size();
  remaining_bodies.reset();
  GeometryModifyEngine* engine = 0;

  if (i == 0)
    return 0;

    // Skip over any bodies that don't have a modify engine.
  while (i--)
  {
    Body* body = remaining_bodies.get();
    TopologyBridge* bridge = 0;
    engine = get_engine(body, &bridge);
    if (engine)
    {
      remaining_bodies.change_to(0);
      engine_bodies.append(body);
      engine_body_sms.append(dynamic_cast<BodySM*>(bridge));
      remaining_bodies.step();
      break;
    }
    remaining_bodies.step();
  }

  // catch case where no engine was found
  if (0 == engine)
  {
    PRINT_WARNING("No geometry modify engine found for this operation.");
    return engine;
  }


    // Get remaining bodies with same modify engine.
  while (i--)
  {
    Body* body = remaining_bodies.get();
    TopologyBridge* bridge = 0;
    if (get_engine(body, &bridge) == engine)
    {
      remaining_bodies.change_to(0);
      engine_bodies.append(body);
      engine_body_sms.append(dynamic_cast<BodySM*>(bridge));
    }
    remaining_bodies.step();
  }

  remaining_bodies.remove_all_with_value(0);
  return engine;
}




// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

CubitStatus GeometryModifyTool::unite( DLIList<BodySM*> &body_sm_list,
                                       DLIList<BodySM*> &new_body_sm_list,
                                       bool keep_old)
{
  //this assumes that all bodies have the same modify engine
  GeometryModifyEngine *gme = get_engine( body_sm_list.get() );
  CubitStatus result = gme->unite(body_sm_list, new_body_sm_list, keep_old);
  return result;
}

CubitStatus
GeometryModifyTool::unite( DLIList<Body*> &bodies,
                           DLIList<Body*> &new_body_list,
                           bool keep_old )
{
  if( bodies.size() <= 1 )
  {
    PRINT_WARNING("There is only one volume in the list. Nothing modified.\n");
    return CUBIT_FAILURE;
  }
  if (!okay_to_modify( bodies, "UNITE" ))
    return CUBIT_FAILURE;

  int i;
  const int count = bodies.size();
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);
  bodies.reset();
  for (i = bodies.size(); i--; )
    entity_list.append_unique(bodies.get_and_step());
  GeometryModifyEngine* gme = common_modify_engine( entity_list, bridge_list );

  if( !gme )
  {
    PRINT_ERROR("Performing UNITE with volumes containing geometry from\n"
      "different modeling engines is not allowed.\n"
      "Delete uncommon geometry on these volumes before operation.\n\n");
    return CUBIT_FAILURE;
  }

  // Cubit can't mesh mixed sheet/solid bodies that are united together. If
  // required, separate the unite between these types.
  CubitStatus result;
  if( GeometryModifyTool::instance()->unite_mixed_models() )
    result = unite_all( gme, bodies, new_body_list, keep_old );
  else
    result = unite_separately( gme, bodies, new_body_list, keep_old );

  if( result == CUBIT_FAILURE )
    PRINT_ERROR("UNITE failed\n");

  return result;
}

CubitStatus
GeometryModifyTool::unite_separately( GeometryModifyEngine *gme_ptr,
                                      DLIList<Body*> &bodies,
                                      DLIList<Body*> &new_body_list,
                                      bool keep_old )
{
  // Cubit can't mesh mixed sheet/solid bodies that are united together. Sort
  // based on these types.
  int i;
  Body *body_ptr;
  DLIList<Body*> solid_body_list;
  DLIList<Body*> sheet_body_list;
  bodies.reset();
  for( i=bodies.size(); i--; )
  {
    body_ptr = bodies.get_and_step();

    if( body_ptr->is_sheet_body() )
      sheet_body_list.append( body_ptr );
    else
      solid_body_list.append( body_ptr );
  }

  if( sheet_body_list.size() == 1 && solid_body_list.size() == 1 )
  {
    PRINT_ERROR( "Cannot unite solid and sheet bodies together\n" );
    return CUBIT_FAILURE;
  }

  // Setup undo
  if( CubitUndo::get_undo_enabled() )
  {
    if( keep_old )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( bodies );
  }

  // Unite solids with each other, and sheets with each other separately
  CubitStatus result1 = CUBIT_SUCCESS;
  CubitStatus result2 = CUBIT_SUCCESS;
  if( solid_body_list.size() > 1 )
    result1 = unite_private( gme_ptr, solid_body_list, new_body_list, keep_old );
  if( sheet_body_list.size() > 1 )
    result2 = unite_private( gme_ptr, sheet_body_list, new_body_list, keep_old );

  // Finish undo
  if( CubitUndo::get_undo_enabled() )
  {
    if( new_body_list.size() )
      CubitUndo::note_result_bodies( new_body_list );
    else
      CubitUndo::remove_last_undo();
  }

  // Return success if both unites successful
  if( result1 == CUBIT_SUCCESS && result2 == CUBIT_SUCCESS )
    return CUBIT_SUCCESS;

  // Return success if either unite was successful
  if( (solid_body_list.size() > 1 && result1 == CUBIT_SUCCESS) ||
      (sheet_body_list.size() > 1 && result2 == CUBIT_SUCCESS) )
  {
    // Give warning if one or the other failed
    if( result1 == CUBIT_FAILURE )
      PRINT_WARNING( "Unite of solid volumes failed\n" );
    if( result2 == CUBIT_FAILURE )
      PRINT_WARNING( "Unite of sheet bodies failed\n" );

    return CUBIT_SUCCESS;
  }

  return CUBIT_FAILURE;
}

CubitStatus
GeometryModifyTool::unite_all( GeometryModifyEngine *gme_ptr,
                               DLIList<Body*> &bodies,
                               DLIList<Body*> &new_body_list,
                               bool keep_old )
{
  // Setup undo
  if( CubitUndo::get_undo_enabled() )
  {
    if( keep_old )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( bodies );
  }

  // Unite solids with each other, and sheets with each other separately
  CubitStatus result = unite_private( gme_ptr, bodies, new_body_list, keep_old );

  // Finish undo
  if( CubitUndo::get_undo_enabled() )
  {
    if( new_body_list.size() )
      CubitUndo::note_result_bodies( new_body_list );
    else
      CubitUndo::remove_last_undo();
  }
  
  return result;
}

// Private workhorse function for unite
CubitStatus
GeometryModifyTool::unite_private( GeometryModifyEngine *gme_ptr,
                                   DLIList<Body*> &body_list,
                                   DLIList<Body*> &new_body_list,
                                   bool keep_old )
{
  if( !body_list.size() )
    return CUBIT_SUCCESS;

  int i, j;
  Body *body_ptr;
  CubitStatus result;

  // Give 1st body all the names of all bodies being united
  std::list<CubitString> names_list;
  DLIList<CubitString*> entity_names;

  body_list.reset();
  for( i=body_list.size(); i--; )
  {
    body_ptr = body_list.get_and_step();

    // See if body has names
    if( body_ptr->num_names() )
    {
      // Put the names in a list
      body_ptr->entity_names( entity_names );
      entity_names.reset();

      // Loop through names
      for( j=entity_names.size(); j--; )
        names_list.push_back( *entity_names.get_and_step() );

      entity_names.clean_out();
      body_ptr->remove_entity_names();
    }
  }

  do_attribute_setup();

  DLIList<TopologyEntity*> entity_list(body_list.size());
  DLIList<TopologyBridge*> bridge_list(body_list.size());
  body_list.reset();
  for( i=body_list.size(); i--; )
    entity_list.append_unique(body_list.get_and_step());
  common_modify_engine( entity_list, bridge_list );

  DLIList<BodySM*> body_sm_list(body_list.size());
  CAST_LIST(bridge_list, body_sm_list, BodySM);

  push_vg_attributes_before_modify(body_sm_list);

  DLIList<int> merged_surface_ids;
  DLIList<int> merged_curve_ids;

  get_merged_curve_and_surface_ids( body_list, merged_surface_ids, merged_curve_ids );

  DLIList<BodySM*> new_body_sm_list;
  result = unite( body_sm_list, new_body_sm_list, keep_old );

  restore_vg_after_modify( new_body_sm_list, body_list, gme_ptr );
  remove_pushed_attributes( new_body_sm_list, body_list );

  DLIList<Body*> result_list;
  if( !finish_sm_op(body_list, new_body_sm_list, result_list) )
    result = CUBIT_FAILURE;

  if( keep_old == CUBIT_FALSE )
    fixup_merged_entities( merged_surface_ids, merged_curve_ids );

  do_attribute_cleanup();

  if( result )
  {
    new_body_list += result_list;

    for( j=result_list.size(); j--; )
    {
      //Add names to 1st body
      std::list<CubitString>::iterator iter, end = names_list.end();
      for (iter = names_list.begin(); iter != end; ++iter)
        result_list.get_and_step()->entity_name( *iter );
    }
  }

  return result;
}

CubitStatus GeometryModifyTool::chop( DLIList<Body*> &bodies,
                                      DLIList<Body*> &intersectBodies,
                                      DLIList<Body*> &outsideBodies,
                                      Body*& leftoversBody,
                                      bool keep_old,
                                      bool nonreg )
{
   leftoversBody = 0;
   if (bodies.size() <= 1)
   {
      PRINT_WARNING("There is only one body in the list. Nothing modified\n");
      return CUBIT_FAILURE;
   }
   if (!okay_to_modify( bodies, "CHOP" ))
     return CUBIT_FAILURE;

  int i;
  const int count = bodies.size();
  DLIList<Body*> original_body_list = bodies;
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);
  bodies.reset();
  for (i = bodies.size(); i--; )
    entity_list.append_unique(bodies.get_and_step());
  GeometryModifyEngine* gme = common_modify_engine( entity_list, bridge_list );

  if ( !gme )
  {
      PRINT_ERROR("Performing CHOP with volumes containing geometry\n"
                  " from different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
  }

   DLIList<BodySM*> body_sm_list;
   CAST_LIST( bridge_list, body_sm_list, BodySM );
   DLIList<BodySM*> intersect_bodies, outside_bodies;
   BodySM *leftovers_body = 0;

   if( CubitUndo::get_undo_enabled() )
   {
     if( keep_old )
       CubitUndo::save_state();
     else
     {
       //Get all the bodies associated with the vertex
       CubitUndo::save_state_with_cubit_file( bodies );
     }
   }

   DLIList<int> merged_surface_ids;
   DLIList<int> merged_curve_ids;
   DLIList<Body*> tmp_bodies(1);
   tmp_bodies.append( bodies.get() );
   if( keep_old == CUBIT_FALSE )
     get_merged_curve_and_surface_ids( tmp_bodies, merged_surface_ids, merged_curve_ids );

   do_attribute_setup();

   // Push attributes down onto the blank body (first one in list).
   DLIList<BodySM*> tmp_body_sm_list;
   body_sm_list.reset();
   tmp_body_sm_list.append(body_sm_list.get());
   push_vg_attributes_before_modify(tmp_body_sm_list);

   CubitStatus result = gme->chop( body_sm_list, intersect_bodies,
                          outside_bodies, leftovers_body, keep_old, nonreg );

   if( result == CUBIT_FAILURE )
   {
     if( CubitUndo::get_undo_enabled() )
       CubitUndo::remove_last_undo();

     PRINT_ERROR("CHOP failed\n");
     remove_pushed_attributes(tmp_body_sm_list, tmp_bodies);
     do_attribute_cleanup();
     return CUBIT_FAILURE;
   }

   DLIList<BodySM*> all_sms = intersect_bodies;
   all_sms += outside_bodies;

   restore_vg_after_modify(all_sms, tmp_bodies, gme);
   remove_pushed_attributes(all_sms, tmp_bodies);

   DLIList<Body*> result_bodies;

   body_sm_list.clean_out();
   body_sm_list += intersect_bodies;

   CubitStatus stat = finish_sm_op(bodies, body_sm_list, result_bodies);

   if( keep_old == CUBIT_FALSE )
     fixup_merged_entities( merged_surface_ids, merged_curve_ids);

   if( CubitUndo::get_undo_enabled() )
   {
     if( stat == CUBIT_SUCCESS )
       CubitUndo::note_result_bodies( result_bodies );
     else
       CubitUndo::remove_last_undo();
   }

   if( stat == CUBIT_FAILURE )
   {
     if( CubitUndo::get_undo_enabled() )
       CubitUndo::remove_last_undo();

     PRINT_ERROR("CHOP failed\n");
     do_attribute_cleanup();
     return CUBIT_FAILURE;
   }
   intersectBodies += result_bodies;
   bodies.clean_out();

   body_sm_list.clean_out();
   body_sm_list += outside_bodies;
   result_bodies.clean_out();
   if (!finish_sm_op(bodies, body_sm_list, result_bodies))
   {
     if( CubitUndo::get_undo_enabled() )
       CubitUndo::remove_last_undo();

     PRINT_ERROR("CHOP failed\n");
     do_attribute_cleanup();
     return CUBIT_FAILURE;
   }
   outsideBodies += result_bodies;

   if( leftovers_body )
   {
     body_sm_list.clean_out();
     body_sm_list.append( leftovers_body );
     result_bodies.clean_out();
     if (!finish_sm_op(bodies, body_sm_list, result_bodies))
     {
       if( CubitUndo::get_undo_enabled() )
         CubitUndo::remove_last_undo();

       PRINT_ERROR("CHOP failed\n");
       do_attribute_cleanup();
       return CUBIT_FAILURE;
     }
     leftoversBody = result_bodies.get();

   }

  if( CubitUndo::get_undo_enabled() )
  {
    if( leftoversBody )
    {
      DLIList<Body*> tmp_list(1);
      tmp_list.append( leftoversBody );
      CubitUndo::note_result_bodies( tmp_list );
    }
    CubitUndo::note_result_bodies( intersectBodies );
    CubitUndo::note_result_bodies( outsideBodies );
  }

  do_attribute_cleanup();
  return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::thicken( DLIList<Body*>& bodies,
                                         DLIList<Body*>& new_bodies,
                                         double depth,
                                         bool both )
{
  if (bodies.size() <= 0)
  {
     PRINT_WARNING("There is only one body in the list. Nothing modified\n");
     return CUBIT_FAILURE;
  }

  if (!okay_to_modify( bodies, "THICKEN" ))
    return CUBIT_FAILURE;

  // Get the GeometryEngine for each Body of the list to check
  // if they are the same and if they are GeometryModifyEngine

  const int count = bodies.size();
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);
  CAST_LIST_TO_PARENT(bodies, entity_list);
  GeometryModifyEngine* gme = common_modify_engine( entity_list, bridge_list );

  if (!gme)
  {
     PRINT_ERROR("Performing THICKEN with volumes containing geometry\n"
                 " from different modeling engines is not allowed.\n"
                 "Delete uncommon geometry on these volumes before operation.\n\n");
     return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state_with_cubit_file( bodies );

  DLIList<BodySM*> new_sms(count);
  DLIList<BodySM*> body_sms(count);
  CAST_LIST(bridge_list, body_sms, BodySM);
  CubitStatus result = gme->thicken( body_sms, new_sms, depth, both);

  // check for resued entities, they have been moved and we need to notify observers
  DLIList<RefEntity*> entities_to_update;
  int i;
  for(i=0; i<new_sms.size(); i++)
  {
    BodySM* bodysm = new_sms.get_and_step();
    DLIList<TopologyBridge*> to_check;
    DLIList<TopologyBridge*> tmp;
    DLIList<Surface*> surfs;
    bodysm->surfaces(surfs);
    DLIList<Curve*> curves;
    bodysm->curves(curves);
    DLIList<Point*> points;
    bodysm->points(points);
    to_check.append(bodysm);
    to_check.append(bodysm->lump());
    CAST_LIST_TO_PARENT(surfs, tmp);
    to_check += tmp;
    CAST_LIST_TO_PARENT(curves, tmp);
    to_check += tmp;
    CAST_LIST_TO_PARENT(points, tmp);
    to_check += tmp;

    int k;
    for(k=0; k<to_check.size(); k++)
      if(BridgeManager* m = to_check.get_and_step()->bridge_manager())
        if(TopologyEntity* t = m->topology_entity())
          entities_to_update.append(CAST_TO(t, RefEntity));

  }

  if (!finish_sm_op(bodies, new_sms, new_bodies))
    result = CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled() )
  {
    if( new_bodies.size() )
      CubitUndo::note_result_bodies( new_bodies );
    else
      CubitUndo::remove_last_undo();
  }

  // Update graphics
  while (entities_to_update.size())
    entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );

  return result;
}


CubitStatus GeometryModifyTool::validate_normals(DLIList<Body*>& bodies,
                                                 RefFace *surf_ref,
                                                 bool reverse)
{
   if (bodies.size() <= 0)
   {
      PRINT_WARNING("There are no entities in the list. Nothing modified\n");
      return CUBIT_FAILURE;
   }

   DLIList<RefEntity*> temp;
   CAST_LIST_TO_PARENT(bodies, temp);
   if ( !same_modify_engine(temp, CUBIT_TRUE))
   {
      PRINT_ERROR("Performing VALIDATE NORMALS with volumes containing geometry\n"
            "      from different modeling engines is not allowed. Delete uncommon\n"
            "      geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   // Get the GeometryEngine for each Body of the list to check
   // if they are the same and if they are GeometryModifyEngine
   GeometryModifyEngine* gePtr1 = get_engine(bodies.get());
   GeometryModifyEngine* gePtr2;
   Body* Body_ptr = NULL;
   bodies.reset();
   for( int i = bodies.size(); i > 0; i--)
   {
      Body_ptr = bodies.get_and_step();
      gePtr2 = get_engine(Body_ptr);
      if (gePtr1 != gePtr2)
      {
         PRINT_ERROR("In GeometryModifyTool::validate_normals\n"
            "  Volumes are associated with different GMEs. \n");
         return CUBIT_FAILURE;
      }
      if ( gePtr2 == NULL)
      {
         PRINT_ERROR("In GeometryModifyTool::validate_normals\n"
            "       Body %d is not associated with a SME.\n",
            Body_ptr->id());
         return CUBIT_FAILURE;
      }
   }


   DLIList<RefFace*> ref_face_list;
   DLIList<RefFace*> free_face_list;
   DLIList<RefFace*> bad_face_list;

   //gather all the faces from all the bodies
   bodies.reset();
   for( int i=bodies.size(); i--; )
   {
     Body* BodyPtr = bodies.get_and_step();
     BodyPtr->ref_faces(free_face_list);
   }

   RefFace *ref_face_ptr;
   RefFace *inter_face_ptr;

   if(surf_ref != NULL)   // getting the starting surface
   {
     ref_face_ptr = surf_ref;
   }
   else
   {
     ref_face_ptr = free_face_list.get_and_step();
   }

   ref_face_list.append(ref_face_ptr);
   free_face_list.remove(ref_face_ptr);
   free_face_list.reset();

   while(free_face_list.size())
   {
      DLIList<RefEdge*> curve_list;
      ref_face_ptr->ref_edges(curve_list);
      free_face_list.reset();
      for(int jj=free_face_list.size(); jj > 0; jj--)  // getting a new searching surface
      {
         inter_face_ptr = free_face_list.get_and_step();
         DLIList<RefEdge*> inter_curve_list;
         inter_face_ptr->ref_edges(inter_curve_list);
         curve_list.reset();
        // PRINT_INFO("base face %d working on face %d\n", ref_face_ptr->id(),inter_face_ptr->id());

         for (int k= curve_list.size(); k > 0; k--)  // looping through all of the surface curves
         {
            RefEdge *ref_check_curve = curve_list.step_and_get();
            inter_curve_list.reset();

            for (int kk = inter_curve_list.size(); kk > 0; kk--)
            {
               RefEdge *check_curve = inter_curve_list.step_and_get();

               if(ref_check_curve == check_curve)  // finding if a surface is connected
               {
                  DLIList<CoEdge*> coedge_list;
                  check_curve->co_edges(coedge_list);

                  CoEdge* first_coedge = coedge_list.get_and_step();
                  CoEdge* second_coedge = coedge_list.get_and_step();

                  if((first_coedge->get_sense() == second_coedge->get_sense() &&
                     !bad_face_list.is_in_list(ref_face_ptr) ) ||
                     (first_coedge->get_sense() != second_coedge->get_sense() &&
                     bad_face_list.is_in_list(ref_face_ptr) ))    // finding if a surface has a fliped normal
                  {
                     bad_face_list.append(inter_face_ptr);
                     PRINT_INFO("Surface %d is not consistent\n", inter_face_ptr->id());
                  }

                  ref_face_list.append(inter_face_ptr);  // adding to the searched list
                  free_face_list.remove(inter_face_ptr); // removing from the unsearched list
               }

            }

         }

      }

      ref_face_list.remove(ref_face_ptr);     // removeing from the searched list
      ref_face_list.last();
      if(ref_face_list.size() <= 0)
      {
         PRINT_ERROR("In GeometryModifyTool::validate_normals\n"
            "       all surfaces must be connected\n");
         return CUBIT_FAILURE;
      }
      ref_face_ptr = ref_face_list.get();
  }

  if (reverse && bad_face_list.size())
  {
     //CubitStatus result = gePtr1->flip_normals(bad_face_list);
     CubitStatus result = GeometryModifyTool::instance()->reverse(bad_face_list);
     if ( result == CUBIT_FAILURE )
     {
       return CUBIT_FAILURE;
    }
  }
  else if(!bad_face_list.size())
  {
    PRINT_INFO("All surfaces are consistent\n");
  }
  return CUBIT_SUCCESS;
}



CubitStatus GeometryModifyTool::subtract( Body* tool_body, DLIList<Body*> &from_bodies,
                                          DLIList<Body*> &new_bodies,
                                          bool imprint,
                                          bool keep_old )
{
  DLIList<Body*> temp_body_list;
  temp_body_list.append(tool_body);
  return subtract(temp_body_list, from_bodies, new_bodies, imprint, keep_old);
}


CubitStatus GeometryModifyTool::subtract( DLIList<Body*>  &tool_body_list,
                                          DLIList<Body*> &from_bodies,
                                          DLIList<Body*> &new_bodies,
                                          bool imprint,
                                          bool keep_old )
{
   if(tool_body_list.size() == 0 )
       return CUBIT_FAILURE;
   if(from_bodies.size() == 0 )
       return CUBIT_FAILURE;
   DLIList<Body*> tem_bodies(tool_body_list);

   // cannot subtract from self
   int old_size = from_bodies.size();
   from_bodies -= tool_body_list;
   if (from_bodies.size() != old_size)
   {
     PRINT_WARNING("Cannot subtract body from itself.  Ignoring \"from\" body.\n");
     if (!from_bodies.size())
       return CUBIT_FAILURE;
   }

   tem_bodies += from_bodies;
   if (!okay_to_modify( tem_bodies, "SUBTRACT" ))
     return CUBIT_FAILURE;

   DLIList<BodySM*> tool_sms(tool_body_list.size());
   DLIList<BodySM*> from_sms(from_bodies.size());
   GeometryModifyEngine* gme = common_modify_engine( tool_body_list, tool_sms );
   GeometryModifyEngine* gme2 = common_modify_engine( from_bodies, from_sms );
   if (!gme || gme != gme2)
   {
      PRINT_ERROR("Performing SUBTRACTION with volumes containing geometry\n"
                  "from different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   if( CubitUndo::get_undo_enabled() )
   {
     if( keep_old )
       CubitUndo::save_state();
     else
     {
       //Get all the bodies associated with the vertex
       DLIList<Body*> bodies;
       bodies += tool_body_list;
       bodies += from_bodies;
       CubitUndo::save_state_with_cubit_file( bodies );
     }
   }

   DLIList<int> merged_surface_ids;
   DLIList<int> merged_curve_ids;
   if( keep_old == CUBIT_FALSE )
     get_merged_curve_and_surface_ids( from_bodies, merged_surface_ids, merged_curve_ids );

     // Do the subtract operation
   DLIList<BodySM*> new_sms;
   CubitStatus result = gme->subtract(tool_sms, from_sms, new_sms, imprint, keep_old );

   if( CubitUndo::get_undo_enabled() && result == CUBIT_FAILURE )
     CubitUndo::remove_last_undo();

   result = finish_sm_op(tem_bodies, new_sms, new_bodies);

   if( CubitUndo::get_undo_enabled() )
   {
     if( result == CUBIT_SUCCESS )
       CubitUndo::note_result_bodies( new_bodies );
     else
       CubitUndo::remove_last_undo();
   }

   if ( result == CUBIT_FAILURE )
   {
     PRINT_ERROR("Subtract FAILED\n" );
     return CUBIT_FAILURE;
   }

   if( keep_old == CUBIT_FALSE )
     fixup_merged_entities( merged_surface_ids, merged_curve_ids);

   return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::intersect( DLIList<Body*> &from_bodies,
                                           DLIList<Body*> &new_bodies,
                                           bool keep_old )
{
  DLIList<Body*> tem_bodies = from_bodies;
  if (!okay_to_modify( tem_bodies, "INTERSECT" ))
    return CUBIT_FAILURE;

  DLIList<BodySM*> from_sm_list(tem_bodies.size());
  GeometryModifyEngine* engine = common_modify_engine(tem_bodies, from_sm_list);
  if ( NULL == engine )
  {
    PRINT_ERROR("Performing INTERSECTION with volumes containing geometry\n"
                "from different modeling engines is not allowed.\n"
                "Delete uncommon geometry on these volumes before operation.\n\n");
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
  {
    if( keep_old )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( from_bodies );
  }

  DLIList<int> merged_surface_ids;
  DLIList<int> merged_curve_ids;
  if( keep_old == CUBIT_FALSE )
    get_merged_curve_and_surface_ids( from_bodies, merged_surface_ids, merged_curve_ids );

  GeometryModifyEngine* gme_ptr = get_engine(from_sm_list.get());
  GeometryQueryEngine* gqe_ptr = gme_ptr->get_gqe();

  DLIList<BodySM*> all_new_bodysms;
  int i,j;
  for( i=0; i<from_sm_list.size(); i++ )
  {
    from_sm_list.reset();
    from_sm_list.step(i);
    BodySM *body1 = from_sm_list.get_and_step();

    for(j=i+1; j<from_sm_list.size(); j++ )
    {
      BodySM *body2 = from_sm_list.get_and_step();

      if( body1 == body2 )
        continue;

      //copy the bodies
      BodySM *body1_copy = gme_ptr->copy_body( body1 );
      BodySM *body2_copy = gme_ptr->copy_body( body2 );

      DLIList<BodySM*> tmp_sm_list(1);
      tmp_sm_list.append( body2_copy );
      DLIList<BodySM*> new_sms;

      CubitStatus result =
          engine->intersect(body1_copy, tmp_sm_list, new_sms, true );

      //delete the copies
      gqe_ptr->delete_solid_model_entities( body1_copy );
      gqe_ptr->delete_solid_model_entities( body2_copy );

      if ( result == CUBIT_FAILURE || new_sms.size() == 0 )
      {
        RefEntity* ref_ent1 = dynamic_cast<RefEntity*>(body1->topology_entity());
        RefEntity* ref_ent2 = dynamic_cast<RefEntity*>(body2->topology_entity());

        PRINT_WARNING("INTERSECTION of %s with %s failed\n",
          ref_ent1->entity_name().c_str(),
          ref_ent2->entity_name().c_str() );
        continue;

      }

      all_new_bodysms += new_sms;
    }
  }

  //now make all the RefEntities
  all_new_bodysms.reset();
  for( i=all_new_bodysms.size(); i--; )
  {
    Body *new_body = GeometryQueryTool::instance()->make_Body(all_new_bodysms.get_and_step());
    if( new_body )
      new_bodies.append( new_body );
  }

  if( CubitUndo::get_undo_enabled() )
  {
    if( all_new_bodysms.size() )
      CubitUndo::note_result_bodies( new_bodies );
    else
      CubitUndo::remove_last_undo();
  }

  if( keep_old == CUBIT_FALSE )
    fixup_merged_entities( merged_surface_ids, merged_curve_ids);

  return CUBIT_SUCCESS;
}



CubitStatus GeometryModifyTool::intersect( Body *tool_body_ptr,
                                           DLIList<Body*> &from_bodies,
                                           DLIList<Body*> &new_bodies,
                                           bool keep_old )
{
  if(tool_body_ptr == NULL )
    return CUBIT_FAILURE;
  if(from_bodies.size() == 0 || from_bodies.get() == NULL )
    return CUBIT_FAILURE;

  DLIList<Body*> tem_bodies = from_bodies;
  tem_bodies.append( tool_body_ptr );
  if (!okay_to_modify( tem_bodies, "INTERSECT" ))
    return CUBIT_FAILURE;

  DLIList<BodySM*> from_sm_list(tem_bodies.size());
  GeometryModifyEngine* engine = common_modify_engine(tem_bodies, from_sm_list);
  if ( NULL == engine )
  {
    PRINT_ERROR("Performing INTERSECTION with volumes containing geometry\n"
                "from different modeling engines is not allowed.\n"
                "Delete uncommon geometry on these volumes before operation.\n\n");
    return CUBIT_FAILURE;
  }

  BodySM* tool_sm = from_sm_list.pop();


  //cannot intersect tool with itself
  from_sm_list.remove_all_with_value( tool_sm );
  if( from_sm_list.size() == 0 )
  {
    PRINT_ERROR("Cannot intersect volume %d from itself\n",
                 tool_body_ptr->ref_volume()->id() );
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
  {
    if( keep_old )
      CubitUndo::save_state();
    else
    {
      DLIList<Body*> bodies;
      bodies.append( tool_body_ptr );
      bodies += from_bodies;
      CubitUndo::save_state_with_cubit_file( bodies );
    }
  }

  DLIList<int> merged_surface_ids;
  DLIList<int> merged_curve_ids;
  if( keep_old == CUBIT_FALSE )
    get_merged_curve_and_surface_ids( from_bodies, merged_surface_ids, merged_curve_ids );

     // Do the intersect operation
  DLIList<BodySM*> new_sms;
  CubitStatus result =
      engine->intersect(tool_sm, from_sm_list, new_sms, keep_old );
  result = finish_sm_op(tem_bodies, new_sms, new_bodies);

  if( CubitUndo::get_undo_enabled() )
  {
    if( result == CUBIT_SUCCESS )
      CubitUndo::note_result_bodies( new_bodies );
    else
      CubitUndo::remove_last_undo();
  }

  if ( result == CUBIT_FAILURE )
  {
    PRINT_ERROR("Intersect FAILED\n" );
    return CUBIT_FAILURE;
  }

  if( keep_old == CUBIT_FALSE )
    fixup_merged_entities( merged_surface_ids, merged_curve_ids);

  return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::imprint( DLIList<Body*> &from_body_list,
                                         DLIList<Body*> &new_body_list,
                                         CubitBoolean keep_old )
{
  if( from_body_list.size() == 1 )
  {
    PRINT_WARNING("Need more than 1 body or volume to imprint.\n");
    return CUBIT_FAILURE;
  }

   //if (get_group_imprint() == CUBIT_FALSE)
   //{
   //  CubitStatus result = imprint_singly( from_body_list, new_body_list, keep_old );
   //  return result;
   //}

     // Check the GeometryEngine for each of the Body's; check to
     // make sure they're all the same
   from_body_list.reset();
   if (!okay_to_modify( from_body_list, "IMPRINT" ))
     return CUBIT_FAILURE;

     //Check for repeats in each individual list and for overlap
     //between the two lists.
   from_body_list.uniquify_ordered();

   DLIList<BodySM*> from_sms(from_body_list.size()), new_sms;
   GeometryModifyEngine* gePtr1 = common_modify_engine(from_body_list, from_sms);
   if ( !gePtr1 )
   {
      PRINT_ERROR("Performing IMPRINT with volumes containing geometry\n"
                  "from different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   if( CubitUndo::get_undo_enabled() )
   {
     if( keep_old )
       CubitUndo::save_state();
     else
       CubitUndo::save_state_with_cubit_file( from_body_list );
   }

   int process_composites = 0;
   if(contains_composites(from_body_list))
     process_composites = 1;

   if(process_composites)
   {
      // Push virtual attributes down to solid model topology before
      // doing the imprint.
      do_attribute_setup();
      push_vg_attributes_before_modify(from_sms);
      // This must be done after pushing the vg atts because it uses them.
      DLIList<TopologyBridge*> tb_list;
      CAST_LIST(from_sms, tb_list, TopologyBridge);
      push_named_attributes_to_curves_and_points(tb_list, "IMPRINTER");
      push_named_attributes_to_curves_and_points(tb_list, "ORIGINAL");
   }

   DLIList<TopologyBridge*> new_tbs, att_tbs;
   CubitStatus result =  gePtr1->imprint(from_sms, new_sms, keep_old, &new_tbs,
     &att_tbs);

   int i, j;
   if(process_composites)
   {
     if(result == CUBIT_SUCCESS)
     {
        // Analyze the results and adjust virtual attributes as necessary.
       DLIList<TopologyBridge*> tb_list;
       CAST_LIST(new_sms, tb_list, TopologyBridge);
        GeometryQueryTool::instance()->ige_attribute_after_imprinting(new_tbs, att_tbs,
          tb_list, from_body_list);

        // Clean up attributes.
        remove_imprint_attributes_after_modify(from_sms, new_sms);

        // Restore the virtual geometry.
        restore_vg_after_modify(new_sms, from_body_list, gePtr1);
     }
     remove_pushed_attributes(new_sms, from_body_list);
   }

   if (get_old_names() == CUBIT_FALSE)
   {
     if (!finish_sm_op(from_body_list, new_sms, new_body_list))
       result = CUBIT_FAILURE;

     if( CubitUndo::get_undo_enabled() )
     {
       if( result == CUBIT_SUCCESS )
         CubitUndo::note_result_bodies( new_body_list );
       else
         CubitUndo::remove_last_undo();
     }

     if(process_composites)
       do_attribute_cleanup();

     return result;
   }

   if(process_composites)
     do_attribute_cleanup();

   // If old_names is true, need to make sure things are deleted in
   // the correct order so that entities get the same @A type extension
   // on their names.

   // Update existing bodies.
   from_body_list.reset();
   for (i = from_body_list.size(); i--; )
   {
     Body* body = from_body_list.get();
     BodySM* body_sm = body->get_body_sm_ptr();
     if (!body_sm)
     {
       GeometryQueryTool::instance()->destroy_dead_entity(body);
       from_body_list.change_to(0);
     }
     else
     {
       remove_dead_entity_names(body);
       GeometryQueryTool::instance()->make_Body(body_sm);
     }
     from_body_list.step();
   }

      // Construct new bodies
   new_sms.reset();
   for (j = new_sms.size(); j--; )
   {
     BodySM* body_sm = new_sms.get_and_step();
     Body* body = GeometryQueryTool::instance()->make_Body(body_sm);
     new_body_list.append(body);
   }
   GeometryQueryTool::instance()->cleanout_deactivated_geometry();

   if( CubitUndo::get_undo_enabled() )
   {
     if( result == CUBIT_SUCCESS )
       CubitUndo::note_result_bodies( new_body_list );
     else
       CubitUndo::remove_last_undo();
   }

   return result;
}

CubitStatus GeometryModifyTool::scale( Body *&body,
                                       const CubitVector& factors,
                                       bool check_to_transform,
                                       bool preview /*=false*/)
{
  if( check_to_transform )
    if (!GeometryQueryTool::instance()->okay_to_transform( body ))
      return CUBIT_FAILURE;

  if (preview)
  {
    GfxPreview::clear();
    DLIList<RefEdge*> edges;
    body->ref_edges(edges);
    for (int i = 0; i < edges.size(); i++)
    {
      GMem poly, prev;
      edges[i]->get_graphics(poly);
      GPoint *prev_points = NULL;
      prev_points = new GPoint[poly.point_list_size()];
      for (int j = 0; j < poly.point_list_size(); j++)
      {
        CubitVector tempV(poly.point_list()[j].x, poly.point_list()[j].y, poly.point_list()[j].z);
        tempV.x(tempV.x()*factors.x());
        tempV.y(tempV.y()*factors.y());
        tempV.z(tempV.z()*factors.z());
        prev_points[j].x = tempV.x();
        prev_points[j].y = tempV.y();
        prev_points[j].z = tempV.z();
      }
      GfxPreview::draw_polyline(prev_points, poly.point_list_size(), CUBIT_BLUE);
      delete [] prev_points;
      prev_points = NULL;
    }
    GfxPreview::flush();
    return CUBIT_SUCCESS;
  }

  BodySM* bodysm = body->get_body_sm_ptr();
  GeometryModifyEngine* engine = get_engine( bodysm );
  CubitStatus result;
  if( !engine )
  {
    GeometryQueryEngine* tmp_engine = bodysm->get_geometry_query_engine();
    result = tmp_engine->scale( bodysm, factors );
  }
  else
    result = engine->scale( bodysm, factors );

  //for non-uniform scaling, topology can change...need to update stuff
  if( factors.x() != factors.y() ||
      factors.y() != factors.z() ||
      factors.z() != factors.z() )
    body = GeometryQueryTool::instance()->make_Body(bodysm);

  if (result)
  {
    CubitTransformMatrix xform;
    xform.scale_about_origin( factors );
    GeometryQueryTool::instance()->notify_intermediate_of_transform( body, xform );
    GeometryQueryTool::instance()->notify_observers_of_transform( body );
  }
  else
    PRINT_ERROR("Scale of %s (%s %d) failed.\n",
      body->entity_name().c_str(), body->class_name(), body->id() );

  return result;
}


/*
CubitStatus GeometryModifyTool::imprint_singly( DLIList<Body*> &from_body_list,
                                                DLIList<Body*> &new_body_list,
                                                CubitBoolean keep_old )
{

     // Check the GeometryEngine for each of the Body's; check to
     // make sure they're all the same
   from_body_list.reset();
   if (!okay_to_modify( from_body_list, "IMPRINT" ))
     return CUBIT_FAILURE;

     //Check for repeats in each individual list and for overlap
     //between the two lists.
   from_body_list.uniquify_ordered();

   DLIList<BodySM*> from_sms(from_body_list.size());
   GeometryModifyEngine* gePtr1 = common_modify_engine(from_body_list, from_sms);
   if ( !gePtr1 )
   {
      PRINT_ERROR("Performing IMPRINT with volumes containing geometry\n"
                  "from different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   CubitBoolean reset_new_ids = CUBIT_FALSE;
   if (get_new_ids() == CUBIT_FALSE) {
     PRINT_WARNING("New ids must be TRUE when group imprint FALSE; "
                   "setting new ids TRUE for this operation only.\n");
     set_new_ids(CUBIT_TRUE);
     reset_new_ids = CUBIT_TRUE;
   }

     // do the imprinting; bounding box should be checked in
     // SMEEngine function
   new_body_list = from_body_list;
   GeometryQueryTool* gqt = GeometryQueryTool::instance();

   if( CubitUndo::get_undo_enabled() )
   {
     if( keep_old )
       CubitUndo::save_state();
     else
       CubitUndo::save_state_with_cubit_file( from_body_list );
   }

   int i;
   for (i = 0; i < new_body_list.size(); i++) {
     for (int j = 1; j < new_body_list.size()-i; j++) {
       new_body_list.reset();
       new_body_list.step(i);
       Body *body_1 = new_body_list.get();
       Body *body_2 = new_body_list.next(j);
       Body *newBodyPtr1 = NULL;
       Body *newBodyPtr2 = NULL;

       if (body_1 == body_2) {
         PRINT_WARNING("Can't imprint a volume with itself.\n");

         if (reset_new_ids == CUBIT_TRUE) set_new_ids(CUBIT_FALSE);

         return CUBIT_FAILURE;
       }

       BodySM *new_sm_1, *new_sm_2;
       CubitStatus status = gePtr1->imprint(body_1->get_body_sm_ptr(),
                                            body_2->get_body_sm_ptr(),
                                            new_sm_1, new_sm_2,
                                            keep_old);

       if ( status != CUBIT_FAILURE &&
            (new_sm_1 != NULL ||
             new_sm_2 != NULL) )
       {
         from_body_list.reset();
         from_body_list.step(i);

           //put the new ones in the new list and
           //remove the olds ones.
         if ( new_sm_1 != NULL )
         {
           if (!body_1->get_body_sm_ptr())
             gqt->destroy_dead_entity(body_1);
           newBodyPtr1 = gqt->make_Body(new_sm_1);
           new_body_list.change_to(newBodyPtr1);
           from_body_list.change_to(0);
         }
         else
           gqt->make_Body(body_1->get_body_sm_ptr());

         new_body_list.step(j);
         from_body_list.step(j);
         if ( new_sm_2 != NULL )
         {
           if (!body_2->get_body_sm_ptr())
             gqt->destroy_dead_entity(body_2);
           newBodyPtr2 = gqt->make_Body(new_sm_2);
           new_body_list.change_to(newBodyPtr2);
           from_body_list.change_to(NULL);
         }
         else
           gqt->make_Body(body_2->get_body_sm_ptr());
       }
       gqt->cleanout_deactivated_geometry();
     }
   }

   from_body_list.remove_all_with_value(NULL);
   Body *temp_body;
   for (i = from_body_list.size(); i > 0; i--)
   {
     temp_body = from_body_list.get_and_step();
     while (new_body_list.move_to(temp_body)) new_body_list.remove();
   }

   if( CubitUndo::get_undo_enabled() )
   {
     if( new_body_list.size() )
       CubitUndo::note_result_bodies( new_body_list );
     else
       CubitUndo::remove_last_undo();
   }

   PRINT_INFO("\n");
   if( DEBUG_FLAG( 153 ) )
   {
     PRINT_INFO( "  New Body(ies) created:");
     new_body_list.reset();
     for (i = 0; i < new_body_list.size(); i++)
     {
       if (i != 0) PRINT_INFO( ",");
       PRINT_INFO( " %d", new_body_list.get_and_step()->id());
     }
     PRINT_INFO("\n");
     PRINT_INFO( "  Original Body(ies) retained:");
     from_body_list.reset();
     for (i = 0; i < from_body_list.size(); i++)
     {
       if (i != 0) PRINT_INFO( ",");
       PRINT_INFO( " %d", from_body_list.get_and_step()->id());
     }
     PRINT_INFO("\n");
   }

   PRINT_INFO( "  New Volume(s) created:");
   new_body_list.reset();
   DLIList<RefVolume*> new_vol_list;
   for (i = 0; i < new_body_list.size(); i++)
   {
     DLIList<RefVolume*> t2;
     new_body_list.get_and_step()->ref_volumes( t2 );
     new_vol_list += t2;
   }

   for( i = 0; i < new_vol_list.size(); i++ )
   {
     if (i != 0) PRINT_INFO( ",");
     PRINT_INFO( " %d", new_vol_list.get_and_step()->id());
   }
   PRINT_INFO("\n");

   PRINT_INFO( "  Original Volume(s) retained:");
   from_body_list.reset();
   DLIList<RefVolume*> from_vol_list;
   for (i = 0; i < from_body_list.size(); i++)
   {
     DLIList<RefVolume*> t2;
     from_body_list.get_and_step()->ref_volumes( t2 );
     from_vol_list += t2;
   }

   for( i = 0; i < from_vol_list.size(); i++ )
   {
     if (i != 0) PRINT_INFO( ",");
     PRINT_INFO( " %d", from_vol_list.get_and_step()->id());
   }
   PRINT_INFO("\n");

   if (reset_new_ids) set_new_ids(CUBIT_FALSE);

   return CUBIT_SUCCESS;
}
*/

CubitStatus GeometryModifyTool::imprint( DLIList<Body*> &body_list,
                                         DLIList<RefEdge*> &ref_edge_list,
                                         DLIList<Body*>& new_body_list,
                                         CubitBoolean keep_old_body,
                                         CubitBoolean show_messages)
{
  // Check the GeometryEngine for each of the bodies; check to
  // make sure they're all the same
  body_list.reset();
  int i;

  if (!okay_to_modify( body_list, "IMPRINT" ))
    return CUBIT_FAILURE;

  const int count = body_list.size() + ref_edge_list.size();
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);
  CAST_LIST_TO_PARENT(body_list, entity_list);
  ref_edge_list.reset();
  for (i = ref_edge_list.size(); i--;)
    entity_list.append(ref_edge_list.get_and_step());

  GeometryModifyEngine* gePtr1 = common_modify_engine(entity_list, bridge_list);

  DLIList<BodySM*> body_sm_list(body_list.size());
  DLIList<Curve*> curve_list(ref_edge_list.size());
  CAST_LIST(bridge_list, body_sm_list, BodySM);
  CAST_LIST(bridge_list, curve_list, Curve);

  if ( !gePtr1 ||
      body_sm_list.size() != body_list.size() ||
      curve_list.size() != ref_edge_list.size() )
  {
    PRINT_ERROR("Performing IMPRINT with volumes containing geometry from\n"
                "different modeling engines is not allowed.\n"
                "Delete uncommon geometry on these volumes before operation.\n\n");
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( body_list );
  }

  int process_composites = 0;
  if(contains_composites(body_list))
    process_composites = 1;

  if(process_composites)
  {
    // Turn certain attributes on.
    do_attribute_setup();
    // Push virtual attributes down to solid model topology before
    // doing the imprint.
    push_vg_attributes_before_modify(body_sm_list);
    // Put "ORIGINAL" attributes on the bodies being imprinted and
    // the curves as these originally existed.
    DLIList<TopologyBridge*> tb_list;
    CAST_LIST(body_sm_list, tb_list, TopologyBridge);
    push_named_attributes_to_curves_and_points(tb_list, "ORIGINAL");
    tb_list.clean_out();
    CAST_LIST(curve_list, tb_list, TopologyBridge);
    push_named_attributes_to_curves_and_points(tb_list, "ORIGINAL");
  }

  DLIList<BodySM*> new_sm_list;
  // The bridges doing the imprinting often get split during the process but
  // because of the way we are making copies, the IMPRINTER attribute doesn't
  // get propagated to them.  temporary_bridges will be filled in with any
  // additional IMPRINTER bridges we need to consider below when deciding whether to
  // keep composite attributes.
  DLIList<TopologyBridge*> temporary_bridges;
  CubitStatus status = gePtr1->imprint( body_sm_list, curve_list,
            new_sm_list, temporary_bridges, keep_old_body, show_messages);

  temporary_bridges.uniquify_ordered();

  if(status == CUBIT_FAILURE)
  {
    if(process_composites)
    {
      remove_pushed_attributes(new_sm_list, body_list);
      do_attribute_cleanup();
    }

    while(temporary_bridges.size())
      delete temporary_bridges.pop();

    return status;
  }
  else
  {
    if(process_composites)
    {
      DLIList<TopologyBridge*> tb_list, new_tbs, att_tbs;
      // Analyze the results and adjust virtual attributes as necessary.
      CAST_LIST(new_sm_list, tb_list, TopologyBridge);
      // The bridges coming back in temporary_bridges may not have IMPRINTER
      // attributes on them becuase of the way they were generated below.  Make
      // sure they get IMPRINTER attributes.
      push_named_attributes_to_curves_and_points(temporary_bridges, "IMPRINTER");
      tb_list += temporary_bridges;
      GeometryQueryTool::instance()->ige_attribute_after_imprinting(new_tbs, att_tbs,
                                                                    tb_list, body_list);

      // Clean up attributes.
      remove_imprint_attributes_after_modify(body_sm_list, new_sm_list);

      // Restore the virtual geometry.
      restore_vg_after_modify(new_sm_list, body_list, gePtr1);
      remove_pushed_attributes(new_sm_list, body_list);
    }
  }

  while(temporary_bridges.size())
    delete temporary_bridges.pop();

   status = finish_sm_op(body_list, new_sm_list, new_body_list);

  if(process_composites)
    do_attribute_cleanup();

   if( CubitUndo::get_undo_enabled() )
   {
     if( status == CUBIT_SUCCESS )
       CubitUndo::note_result_bodies( new_body_list );
     else
       CubitUndo::remove_last_undo();
   }

   return status;
}

CubitStatus GeometryModifyTool::imprint( DLIList<RefFace*> &ref_face_list,
                                         DLIList<RefEdge*> &ref_edge_list,
                                         DLIList<Body*>& new_body_list,
                                         CubitBoolean keep_old_body )
{
  //get the owning bodies of the faces and edges
  DLIList<Body*> body_list;
  int j;
  for(j=ref_face_list.size(); j--;)
    ref_face_list.get_and_step()->bodies( body_list );
  for(j=ref_edge_list.size(); j--;)
    ref_edge_list.get_and_step()->bodies( body_list );
  body_list.uniquify_unordered();
  if (!okay_to_modify( body_list, "IMPRINT" ))
    return CUBIT_FAILURE;

  DLIList<ModelEntity*> temp_list, temp_list_2, body_me_list;
  CAST_LIST_TO_PARENT(ref_face_list, temp_list);
  CAST_LIST_TO_PARENT(ref_edge_list, temp_list_2);
  temp_list += temp_list_2;
  ModelQueryEngine::instance()->query_model(temp_list, DagType::body_type(), body_me_list );

  DLIList<Surface*> surf_list(ref_face_list.size());
  DLIList<Curve*> curve_list(ref_edge_list.size());
  GeometryModifyEngine* gePtr1 = common_modify_engine( ref_face_list,ref_edge_list,surf_list,
                             curve_list,true);
  if ( !gePtr1 )
  {
    PRINT_ERROR("Performing IMPRINT with volumes containing geometry from\n"
            "different modeling engines is not allowed.\n"
            "Delete uncommon geometry on these volumes before operation.\n\n");
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_face_list );
  }

  int process_composites = 0;
  if(contains_composites(body_list))
    process_composites = 1;

  int i;
  DLIList<BodySM*> body_sm_list;
  if(process_composites)
  {
    // Turn certain attributes on.
    do_attribute_setup();
    for(i=body_list.size(); i--;)
      body_sm_list.append_unique(body_list.get_and_step()->get_body_sm_ptr());
    // Push virtual attributes down to solid model topology before
    // doing the imprint.
    push_vg_attributes_before_modify(body_sm_list);
    // Put "ORIGINAL" attributes on the bodies being imprinted and
    // the curves as these originally existed.
    DLIList<TopologyBridge*> tmp_tb_list;
    CAST_LIST(surf_list, tmp_tb_list, TopologyBridge);
    push_named_attributes_to_curves_and_points(tmp_tb_list, "ORIGINAL");
  }

  DLIList<BodySM*> new_sm_list;
  // The bridges doing the imprinting often get split during the process but
  // because of the way we are making copies, the IMPRINTER attribute doesn't
  // get propagated to them.  temporary_bridges will be filled in with any
  // additional IMPRINTER bridges we need to consider below when deciding whether to
  // keep composite attributes.
  DLIList<TopologyBridge*> temporary_bridges;
  CubitStatus status = gePtr1->imprint( surf_list, curve_list, temporary_bridges,
                                         new_sm_list, keep_old_body );
  temporary_bridges.uniquify_ordered();

  if(process_composites)
  {
    // Analyze the results and adjust virtual attributes as necessary.
    DLIList<TopologyBridge*> tb_list;
    CAST_LIST(new_sm_list, tb_list, TopologyBridge);
    // The bridges coming back in temporary_bridges may not have IMPRINTER
    // attributes on them becuase of the way they were generated below.  Make
    // sure they get IMPRINTER attributes.
    push_named_attributes_to_curves_and_points(temporary_bridges, "IMPRINTER");
    tb_list += temporary_bridges;
    DLIList<TopologyBridge*> new_tbs, att_tbs;
    GeometryQueryTool::instance()->ige_attribute_after_imprinting(new_tbs, att_tbs,
                    tb_list, body_list);

    // Clean up attributes.
    remove_imprint_attributes_after_modify(body_sm_list, new_sm_list);

    restore_vg_after_modify(body_sm_list, body_list, gePtr1);
    remove_pushed_attributes(body_sm_list, body_list);
  }

  while(temporary_bridges.size())
    delete temporary_bridges.pop();

   status = finish_sm_op(body_list, new_sm_list, new_body_list);

  if(process_composites)
    do_attribute_cleanup();

  if( CubitUndo::get_undo_enabled() )
   {
     if( status == CUBIT_SUCCESS )
       CubitUndo::note_result_bodies( new_body_list );
     else
       CubitUndo::remove_last_undo();
   }

   return status;
}

CubitStatus GeometryModifyTool::imprint( DLIList<Surface*> &surface_list,
                                         DLIList<DLIList<Curve*>*> &curve_lists_list,
                                         Body*& /*new_body*/,
                                         CubitBoolean keep_old_body,
                                         CubitBoolean expand)
{
  int i;
  DLIList<Curve*> *curve_list_ptr;

  // Check to see if any curves exist - if none, just exit
  int have_curves = 0;
  for( i=curve_lists_list.size(); i--; )
  {
    curve_list_ptr = curve_lists_list.get_and_step();
    if( curve_list_ptr->size() )
    {
      have_curves++;
      break;
    }
  }
  if( !have_curves )
    return CUBIT_SUCCESS;

  // Get parent bodies
  DLIList<Body*> old_body_list;
  surface_list.reset();
  for( i=surface_list.size(); i--; )
  {
    Surface *surf_ptr = surface_list.get_and_step();
    RefEntity* ref_ent = dynamic_cast<RefEntity*>(surf_ptr->topology_entity());
    RefFace *ref_face_ptr = CAST_TO( ref_ent, RefFace );

    if( ref_face_ptr )
    {
      DLIList<Body*> body_list;
      ref_face_ptr->bodies( body_list );
      old_body_list.merge_unique( body_list );
    }
  }

  if( old_body_list.size() > 1 )
  {
    PRINT_ERROR( "This operation requires all surfaces to be from the same volume\n" );
    // Note: this restriction could be pretty easily lifted by sorting the
    //       input lists and calling the GeometryModifyEngine for each body
    //       separately, or having engines handle this.
    return CUBIT_FAILURE;
  }


  // In order to support imprinting on composite surfaces we will
  // get any surfaces underlying the surfaces passed in.  For now
  // we will only do this if a single surface is coming in but
  // it could be extended for multiple surfaces as well.
  DLIList<Surface*> new_surface_list;
  if(surface_list.size() == 1)
  {
    GeometryQueryEngine *gqe = surface_list.get()->get_geometry_query_engine();
    DLIList<TopologyBridge*> tbs;
    gqe->get_underlying_surfaces(surface_list.get(), tbs);
    if(tbs.size() > 0)
    {
      for(int k=tbs.size(); k--;)
        new_surface_list.append(dynamic_cast<Surface*>(tbs.get_and_step()));
    }
    else
      new_surface_list.append(surface_list.get());
  }
  else
    new_surface_list = surface_list;

  // Check engines - must all be the same
  GeometryModifyEngine* gme;
  new_surface_list.reset();
  gme = get_engine( new_surface_list.get() );
  for( i=new_surface_list.size(); i--; )
  {
    Surface *surf_ptr = new_surface_list.get_and_step();
    GeometryModifyEngine* gme2 = get_engine( surf_ptr );
    if( gme != gme2 )
    {
      PRINT_ERROR( "All surfaces being imprinted must be from the same geometry engine\n" );
      return CUBIT_FAILURE;
    }
  }

  int j;
  for( i=curve_lists_list.size(); i--; )
  {
    curve_list_ptr = curve_lists_list.get_and_step();
    for( j=curve_list_ptr->size(); j--; )
    {
      Curve *curve_ptr = curve_list_ptr->get_and_step();
      GeometryModifyEngine* gme2 = get_engine( curve_ptr );
      if( gme != gme2 )
      {
        PRINT_ERROR( "Curves used to imprint must be from same geometry engine as Surface\n" );
        return CUBIT_FAILURE;
      }
    }
  }

   if( CubitUndo::get_undo_enabled() )
   {
     if( keep_old_body )
       CubitUndo::save_state();
     else
       CubitUndo::save_state_with_cubit_file( old_body_list );
   }

  int process_composites = 0;
  if(contains_composites(old_body_list))
    process_composites = 1;

  BodySM* new_body_sm = 0;
  DLIList<BodySM*> body_sm_list;

  DLIList<TopologyBridge*> tb_list;
  if(process_composites)
  {
    do_attribute_setup();
    for(i=old_body_list.size(); i--;)
      body_sm_list.append_unique(old_body_list.get_and_step()->get_body_sm_ptr());
    push_vg_attributes_before_modify(body_sm_list);
 //   push_imprint_attributes_before_modify(body_sm_list);
    DLIList<TopologyBridge*> tmp_tb_list;
    CAST_LIST(new_surface_list, tmp_tb_list, TopologyBridge);
    push_named_attributes_to_curves_and_points(tmp_tb_list, "ORIGINAL");

    for(i=curve_lists_list.size(); i>0; i--)
    {
      DLIList<Curve*> *cur_list = curve_lists_list.get_and_step();
      for(j=cur_list->size(); j>0; j--)
      {
        Curve *cur_curve = cur_list->get_and_step();
        tb_list.append(cur_curve);
      }
    }
    push_named_attributes_to_curves_and_points(tb_list, "IMPRINTER");
  }

  DLIList<TopologyBridge*> new_tbs, att_tbs;
  CubitStatus status = gme->imprint( new_surface_list, curve_lists_list,
    new_body_sm, keep_old_body, expand, &new_tbs, &att_tbs );

  DLIList<Body*> new_body_list;
  DLIList<BodySM*> new_sm_list;
  new_sm_list.append( new_body_sm );

  if(process_composites)
  {
    // Analyze the results and adjust virtual attributes as necessary.
       DLIList<TopologyBridge*> tmp_tb_list;
       CAST_LIST(new_sm_list, tmp_tb_list, TopologyBridge);
       tb_list += tmp_tb_list;
    GeometryQueryTool::instance()->ige_attribute_after_imprinting(new_tbs, att_tbs,
                    tb_list, old_body_list);

    // Clean up attributes.
    remove_imprint_attributes_after_modify(body_sm_list, new_sm_list);

    restore_vg_after_modify(body_sm_list, old_body_list, gme);
    remove_pushed_attributes(body_sm_list, old_body_list);
  }

  status = finish_sm_op(old_body_list, new_sm_list, new_body_list);

  if(process_composites)
    do_attribute_cleanup();

  if( CubitUndo::get_undo_enabled() )
  {
    if( status == CUBIT_SUCCESS )
      CubitUndo::note_result_bodies( new_body_list );
    else
      CubitUndo::remove_last_undo();
  }

  return status;
}

CubitStatus GeometryModifyTool::imprint( DLIList<Body*> &body_list,
                                         DLIList<CubitVector*> &vector_list,
                                         DLIList<Body*>& new_body_list,
                                         CubitBoolean keep_old_body,
                                         CubitBoolean merge )
{
  // Check the GeometryEngine for each of the RefEdges; check to
  // make sure they're all the same
  body_list.reset();

  if (!okay_to_modify( body_list, "IMPRINT" ))
    return CUBIT_FAILURE;

  DLIList<BodySM*> body_sm_list(body_list.size()), new_sm_list;
  GeometryModifyEngine* gePtr1 = common_modify_engine(body_list, body_sm_list);
  if ( !gePtr1 )
  {
    PRINT_ERROR("Performing IMPRINT with volumes containing geometry from\n"
                "different modeling engines is not allowed.\n"
                "Delete uncommon geometry on these volumes before operation.\n\n");
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( body_list );
  }

  int process_composites = 0;
  if(contains_composites(body_list))
    process_composites = 1;

  int i;
  DLIList<TopologyBridge*> temporary_bridges;
  if(process_composites)
  {
    // Turn certain attributes on.
    do_attribute_setup();
    // Push virtual attributes down to solid model topology before
    // doing the imprint.
    push_vg_attributes_before_modify(body_sm_list);
    // Create temporary bridges for the vector positions.  We do
    // this so that we can put an IMPRINTER attribute on them
    // and use them later for deciding whether to keep composite
    // attributes or not.
    for(i=vector_list.size(); i>0; i--)
    {
      CubitVector *vec = vector_list.get_and_step();
      Point *pt = gePtr1->make_Point(*vec);
      temporary_bridges.append(pt);
    }
    push_named_attributes_to_curves_and_points(temporary_bridges, "IMPRINTER");
    DLIList<TopologyBridge*> tmp_tb_list;
    CAST_LIST(body_sm_list, tmp_tb_list, TopologyBridge);
    // Put "ORIGINAL" attributes on the bridges that originally existed.
    push_named_attributes_to_curves_and_points(tmp_tb_list, "ORIGINAL");
  }

  DLIList<TopologyBridge*> new_tbs, att_tbs;
  CubitStatus status = gePtr1->imprint( body_sm_list, vector_list,new_sm_list,
    keep_old_body, &new_tbs,&att_tbs );

  temporary_bridges.uniquify_ordered();

  if(process_composites)
  {
    // Analyze the results and adjust virtual attributes as necessary.
    DLIList<TopologyBridge*> tb_list;
    CAST_LIST(new_sm_list, tb_list, TopologyBridge);
    tb_list += temporary_bridges;
    GeometryQueryTool::instance()->ige_attribute_after_imprinting(new_tbs, att_tbs,
                      tb_list, body_list);

    while(temporary_bridges.size())
      delete temporary_bridges.pop();

    // Clean up attributes.
    remove_imprint_attributes_after_modify(body_sm_list, new_sm_list);

    // Restore the virtual geometry.
    restore_vg_after_modify(new_sm_list, body_list, gePtr1);
    remove_pushed_attributes(new_sm_list, body_list);
  }

  status = finish_sm_op(body_list, new_sm_list, new_body_list);

  if(process_composites)
    do_attribute_cleanup();

   if( merge )
   {
     DLIList<Body*> bodies_to_merge;
     int i;
     for( i=new_body_list.size(); i--; )
     {
       Body *tmp_body = new_body_list.get_and_step();
       DLIList<RefEdge*> ref_edge_list;
       tmp_body->ref_edges( ref_edge_list );

       int j;
       for( j=ref_edge_list.size(); j--; )
       {
         RefEdge *tmp_edge = ref_edge_list.get_and_step();
         DLIList<Body*> body_list;
         tmp_edge->bodies( body_list );
         bodies_to_merge.merge_unique( body_list );
       }
     }
     MergeTool::instance()->merge_bodies( bodies_to_merge );
   }

   if( CubitUndo::get_undo_enabled() )
   {
     if( status == CUBIT_SUCCESS )
       CubitUndo::note_result_bodies( new_body_list );
     else
       CubitUndo::remove_last_undo();
   }

   return status;
}
CubitStatus GeometryModifyTool::project_edges( DLIList<RefFace*> &ref_face_list,
                                               DLIList<RefEdge*> &ref_edge_list_in,
                                               DLIList<RefEdge*> &ref_edge_list_new,
                                               CubitBoolean trim_projected)
{
  int i, j;
  
   // Check the GeometryEngine for each of the RefEdges; check to
   // make sure they're all the same
  DLIList<Surface*> surface_list(ref_face_list.size());
  DLIList<Curve*> curve_list_in(ref_edge_list_in.size()), curve_list_new;
  GeometryModifyEngine* gme = common_modify_engine( ref_face_list,
                                                    ref_edge_list_in,
                                                    surface_list,
                                                    curve_list_in );

  if ( !gme ) {
    PRINT_ERROR("In GeometryTool::create_blend\n"
                "       Curves have different modify engines.\n");
    return CUBIT_FAILURE;
  }

   CubitStatus status = gme->
     project_edges( surface_list, curve_list_in, curve_list_new);

  if( CubitUndo::get_undo_enabled() && status == CUBIT_SUCCESS )
    CubitUndo::save_state();


   curve_list_new.reset();
   
   if(trim_projected){
     DLIList<Curve*> tmp_curves, all_new_curves;
     Curve* tmp_curve;
     Surface* tmp_surface;
     
     for(i = 0; i< surface_list.size(); i++){
       tmp_curves.clean_out();
       tmp_surface = surface_list.get_and_step();
       for(j=0; j<curve_list_new.size(); j++){
         tmp_curve = curve_list_new.get_and_step();
         status = gme->curve_surface_intersection( tmp_surface, tmp_curve, tmp_curves);
       }
       all_new_curves += tmp_curves;

     }

     if(!all_new_curves.size()){
       if(!curve_list_new.size()){
         PRINT_ERROR("Projection resulted in no curves.\n");
         return CUBIT_FAILURE;
       }
       else{
         PRINT_WARNING("No curve remained after trimming operation.  \n \tCurve projection may lie completely outside of trimmed surface.\n");
       }
     }
     
       //fix this...
       //can we just cleanout this list or do we need to delete the entities in it.
     for( i = 0; i< curve_list_new.size(); i++ )
     {
       Curve *tmp_curve = curve_list_new.get_and_step();
       gme->get_gqe()->delete_solid_model_entities( tmp_curve );
     }
     
     curve_list_new.clean_out();
     curve_list_new = all_new_curves;
     
     if( CubitUndo::get_undo_enabled() && status == CUBIT_SUCCESS )
       CubitUndo::save_state();
   }
   
  

   curve_list_new.reset();
   for (i = curve_list_new.size(); i--; )
   {
     Curve* curve = curve_list_new.get_and_step();
     RefEdge* new_edge = GeometryQueryTool::instance()->make_free_RefEdge(curve);
     PRINT_INFO("Created Curve %d\n", new_edge->id());
     ref_edge_list_new.append(new_edge);
     if( CubitUndo::get_undo_enabled() && new_edge )
       CubitUndo::note_result_entity( new_edge );
   }

   if( CubitUndo::get_undo_enabled() )
   {
     if( ref_edge_list_new.size() == 0 )
       CubitUndo::remove_last_undo();
   }

   return status;
}

CubitStatus
GeometryModifyTool::imprint_projected_edges(DLIList<RefFace*> &ref_face_list,
                                            DLIList<RefEdge*> &ref_edge_list_in,
                                            DLIList<Body*>& new_body_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean keep_free_edges)
{

   DLIList<Surface*> surface_list(ref_face_list.size());
   DLIList<Curve*> curve_list(ref_edge_list_in.size());
   GeometryModifyEngine* gme = common_modify_engine( ref_face_list,
                                                     ref_edge_list_in,
                                                     surface_list,
                                                     curve_list );
   if ( !gme )
   {
      PRINT_ERROR("Performing IMPRINT with volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   DLIList<ModelEntity*> query_output, query_input(ref_face_list.size());
   CAST_LIST_TO_PARENT(ref_face_list, query_input);
   ModelQueryEngine::instance()->query_model(query_input, DagType::body_type(), query_output);
   DLIList<Body*> body_list(query_output.size());
   CAST_LIST(query_output, body_list, Body);
   DLIList<BodySM*> new_sm_list;

   if( CubitUndo::get_undo_enabled() )
   {
     if( keep_old_body )
       CubitUndo::save_state();
     else
       CubitUndo::save_state_with_cubit_file( ref_face_list );
   }

   CubitStatus status = gme->imprint_projected_edges( surface_list, curve_list,
                                                      new_sm_list, keep_old_body,keep_free_edges);

   if (!finish_sm_op(body_list, new_sm_list, new_body_list))
     status = CUBIT_FAILURE;

   if( CubitUndo::get_undo_enabled() )
   {
     if( status == CUBIT_FAILURE )
       CubitUndo::remove_last_undo();
     else
       CubitUndo::note_result_bodies( new_body_list );
   }

   return status;
}

CubitStatus
GeometryModifyTool::imprint_projected_edges(DLIList<RefFace*> &ref_face_list,
                                            DLIList<Body*> &body_list,
                                            DLIList<RefEdge*> &ref_edge_list_in,
                                            DLIList<Body*>& new_body_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean keep_free_edges)
{

   DLIList<Surface*> surface_list(ref_face_list.size());
   DLIList<Curve*>     curve_list(ref_edge_list_in.size());
   DLIList<BodySM*>  body_sm_list(body_list.size()), new_sm_list;
   GeometryModifyEngine* gme = common_modify_engine(ref_face_list,
                                                    ref_edge_list_in,
                                                    surface_list,
                                                    curve_list);
   if (!gme || gme != common_modify_engine(body_list, body_sm_list))
   {
      PRINT_ERROR("Performing IMPRINT with volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

     // Get RefFace bodies
   DLIList<ModelEntity*> query_output, query_input(ref_face_list.size());
   CAST_LIST_TO_PARENT(ref_face_list, query_input);
   ModelQueryEngine::instance()->
     query_model( query_input, DagType::body_type(), query_output );
   DLIList<Body*> face_body_list;
   CAST_LIST(query_output, face_body_list, Body);

   CubitStatus status = gme->imprint_projected_edges( surface_list,
                                                      body_sm_list,
                                                      curve_list,
                                                      new_sm_list,
                                                      keep_old_body,
                                                      keep_free_edges);
   face_body_list += body_list;
   if (!finish_sm_op(face_body_list, new_sm_list, new_body_list))
     status = CUBIT_FAILURE;

   return status;

}


RefEntity* GeometryModifyTool::copy_refentity( RefEntity *old_entity )
{
   RefEntity *new_entity = NULL;

   if( old_entity == NULL )
       return (RefEntity*) NULL;

   if( CAST_TO( old_entity, RefVolume ) )
       return (RefEntity*) NULL;

   RefFace *old_face = CAST_TO( old_entity, RefFace );
   if( old_face )
   {
      RefFace *new_face = instance()->make_RefFace( old_face );
      new_entity = CAST_TO( new_face, RefEntity );
   }

   RefEdge *old_edge = CAST_TO( old_entity, RefEdge );
   if( old_edge )
   {
      RefEdge *new_edge = instance()->make_RefEdge( old_edge );
      new_entity = CAST_TO( new_edge, RefEntity );
   }

   RefVertex *old_vert = CAST_TO( old_entity, RefVertex );
   if( old_vert )
   {
      RefVertex *new_vert = instance()->make_RefVertex( old_vert->coordinates() );
      new_entity = CAST_TO( new_vert, RefEntity );
   }

   return new_entity;
}

//  Calls solid modeling engine to webcut with a sheet body.
CubitStatus GeometryModifyTool::webcut_with_sheet( DLIList<Body*> &webcut_body_list,
                                                   Body *sheet_body,
                                                   DLIList<Body*> &new_bodies,
                                                   DLIList<Body*> &neighboring_bodies,
                                                   ImprintType imprint_type,
                                                   CubitBoolean merge,
                                                   CubitBoolean preview)
{
   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<Body*> original_body_list = webcut_body_list;
   webcut_body_list.append(sheet_body);
   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);
   if ( !gme )
   {
      PRINT_ERROR("Performing WEBCUTS with volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }
   BodySM* tool_sm = body_sm_list.pop();
 
   DLIList<int> merged_surface_ids;
   DLIList<int> merged_curve_ids;
   DLIList<BodySM*> neighbor_imprint_list;
   DLIList<Body*> bodies_to_modify;

   if (!preview)
   {
     if( CubitUndo::get_undo_enabled() )
     {
       DLIList<Body*> bodies_to_save;
       bodies_to_save += webcut_body_list;
       bodies_to_save += neighboring_bodies;
       CubitUndo::save_state_with_cubit_file( bodies_to_save );
     }

     int i;
     for( i=neighboring_bodies.size(); i--; )
     {
       Body *neighbor_body = neighboring_bodies.get_and_step();
       BodySM *tmp_body = neighbor_body->get_body_sm_ptr(); 
       GeometryModifyEngine *neighbor_gme = get_engine( tmp_body );
       
       if( gme == neighbor_gme )
         neighbor_imprint_list.append( tmp_body ); 
     }

     do_attribute_setup();
     DLIList<BodySM*> bodies_sm_to_modify;
     bodies_sm_to_modify += body_sm_list;
     bodies_sm_to_modify += neighbor_imprint_list;
     push_vg_attributes_before_modify( bodies_sm_to_modify );
     bodies_to_modify += webcut_body_list;
     bodies_to_modify += neighboring_bodies;
     get_merged_curve_and_surface_ids( bodies_to_modify, merged_surface_ids, merged_curve_ids );
   }

   CubitStatus stat = gme->webcut (
       body_sm_list, tool_sm, neighbor_imprint_list,
       new_sms, imprint_type, preview  );

   if (!preview)
   {
     restore_vg_after_modify(new_sms, bodies_to_modify, gme);
     remove_pushed_attributes(new_sms, bodies_to_modify );
   }

   if( stat == CUBIT_FAILURE )
   {
     do_attribute_cleanup();

     if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
     {
       if( stat == CUBIT_SUCCESS )
         CubitUndo::note_result_bodies( new_bodies );
       else
         CubitUndo::remove_last_undo();
     }
     return CUBIT_FAILURE;
   }

   // finish up if not a preview
   if (!preview)
   {
     stat = finish_webcut( webcut_body_list, new_sms, merge, stat, new_bodies,
                           &merged_surface_ids, &merged_curve_ids );
      // leave webcut_body_list as we found it -- remove appended tool body
     webcut_body_list.pop();

     do_attribute_cleanup();

     if( CubitUndo::get_undo_enabled() ) 
     {
       if( stat == CUBIT_SUCCESS )
         CubitUndo::note_result_bodies( new_bodies );
       else
         CubitUndo::remove_last_undo();
     }
   }


   return stat;
}
//  Calls solid modeling engine to webcut with a sheet body.
CubitStatus GeometryModifyTool::webcut_with_extended_sheet( DLIList<Body*> &webcut_body_list,
                                                           DLIList<RefFace*> &ref_face_list,
                                                           DLIList<Body*> &new_bodies,
                                                           DLIList<Body*> &neighboring_bodies,
                                                           int &num_cut,
                                                           ImprintType imprint_type,
                                                           CubitBoolean merge,
                                                           CubitBoolean preview)
{
   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<RefEntity*> ref_entity_list;
   CAST_LIST( ref_face_list, ref_entity_list, RefEntity );
   if( !same_modify_engine( ref_entity_list ) )
   {
     PRINT_ERROR( "Extending surfaces from different geometry engines is\n"
       "       not allowed.\n\n" );
     return CUBIT_FAILURE;
   }

   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);

   if( !gme )
   {
     PRINT_ERROR("Performing WEBCUTS on volumes containing geometry from\n"
                  "different modeling engines is not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   Surface* surf_ptr = 0;
   TopologyBridge* bridge = ref_face_list.get()->bridge_manager()->topology_bridge(gme->get_gqe());
   surf_ptr = dynamic_cast<Surface*>(bridge);
   
   if ( !surf_ptr )
   {
      PRINT_ERROR("Performing WEBCUTS on volumes containing geometry from a\n"
                  "different modeling engine than the extended surfaces is\n"
                  "not allowed.\n"
                  "Delete uncommon geometry on these volumes before operation.\n\n");
      return CUBIT_FAILURE;
   }

   GeometryQueryEngine *gqe = gme->get_gqe();

   int i;
   RefFace *ref_face_ptr;
   DLIList<Surface*> surf_list;
   for( i=ref_face_list.size(); i--; )
   {
     ref_face_ptr = ref_face_list.get_and_step();
     bridge = ref_face_ptr->bridge_manager()->topology_bridge(gqe);
     surf_ptr = dynamic_cast<Surface*>(bridge);
     surf_list.append( surf_ptr );
   }   

   DLIList<int> merged_surface_ids;
   DLIList<int> merged_curve_ids;
   DLIList<BodySM*> neighbor_imprint_list;
   DLIList<Body*> bodies_to_modify;

   if (!preview)
   {
     if( CubitUndo::get_undo_enabled() )
     {
       DLIList<Body*> bodies_to_save;
       bodies_to_save += webcut_body_list;
       bodies_to_save += neighboring_bodies;
       CubitUndo::save_state_with_cubit_file( bodies_to_save );
     }

     int i;
     for( i=neighboring_bodies.size(); i--; )
     {
       Body *neighbor_body = neighboring_bodies.get_and_step();
       BodySM *tmp_body = neighbor_body->get_body_sm_ptr(); 
       GeometryModifyEngine *neighbor_gme = get_engine( tmp_body );
       
       if( gme == neighbor_gme )
         neighbor_imprint_list.append( tmp_body ); 
     }

     do_attribute_setup();
     DLIList<BodySM*> bodies_sm_to_modify;
     bodies_sm_to_modify += body_sm_list;
     bodies_sm_to_modify += neighbor_imprint_list;
     push_vg_attributes_before_modify( bodies_sm_to_modify );
     bodies_to_modify += webcut_body_list;
     bodies_to_modify += neighboring_bodies;
     get_merged_curve_and_surface_ids( bodies_to_modify, merged_surface_ids, merged_curve_ids );
   }

   // FIXME: this was not implemented correctly in CGMA prior to merge with CGM 12.0
   CubitStatus stat = CUBIT_FAILURE;
   //CubitStatus stat = gme->webcut_with_extended_sheet( body_sm_list,
   //  surf_list, neighbor_imprint_list, new_sms, num_cut, imprint_type,
   //  preview );

   if (!preview)
   {
     restore_vg_after_modify(new_sms, bodies_to_modify, gme);
     remove_pushed_attributes(new_sms, bodies_to_modify);
     stat = finish_webcut(webcut_body_list, new_sms, merge, stat, new_bodies,
                          &merged_surface_ids, &merged_curve_ids );
     do_attribute_cleanup();

     if( CubitUndo::get_undo_enabled() ) 
     {
       if( stat == CUBIT_SUCCESS )
         CubitUndo::note_result_bodies( new_bodies );
       else
         CubitUndo::remove_last_undo();
     }
   }

   return stat;
}


CubitStatus GeometryModifyTool::webcut_with_sweep_surfaces_rotated(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefFace*> &tool_faces,
                            CubitVector &point,
                            CubitVector &sweep_axis,
                            double angle,
                            RefFace* stop_surf,
                            bool up_to_next,
                            DLIList<Body*> &new_bodies,
                            DLIList<Body*> &neighboring_bodies,
                            ImprintType imprint_type,
                            CubitBoolean merge,
                            CubitBoolean preview)
{

   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);
   DLIList<Surface*> surfaces_to_sweep;
   if (gme)
   {
     GeometryModifyEngine* surf_gme = common_modify_engine( tool_faces, surfaces_to_sweep );
     if( gme != surf_gme )
       PRINT_ERROR("All geometry not from the same modeling engine.\n");
   }

   Surface *stop_surface = NULL;
   if( stop_surf )
   {
     GeometryModifyEngine *tmp_gme = get_engine(stop_surf);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Stop surface %d is not of same modeling engine as other input geometry.\n",
                    stop_surf->id() );
       return CUBIT_FAILURE;
     }

     stop_surface = stop_surf->get_surface_ptr();
   } 

   DLIList<int> merged_surface_ids;
   DLIList<int> merged_curve_ids;
   DLIList<BodySM*> neighbor_imprint_list;
   DLIList<Body*> bodies_to_modify;

   if (!preview)
   {
     if( CubitUndo::get_undo_enabled() )
     {
       DLIList<Body*> bodies_to_save;
       bodies_to_save += webcut_body_list;
       bodies_to_save += neighboring_bodies;
       CubitUndo::save_state_with_cubit_file( bodies_to_save );
     }

     int i;
     for( i=neighboring_bodies.size(); i--; )
     {
       Body *neighbor_body = neighboring_bodies.get_and_step();
       BodySM *tmp_body = neighbor_body->get_body_sm_ptr(); 
       GeometryModifyEngine *neighbor_gme = get_engine( tmp_body );
       
       if( gme == neighbor_gme )
         neighbor_imprint_list.append( tmp_body ); 
     }

     do_attribute_setup();
     DLIList<BodySM*> bodies_sm_to_modify;
     bodies_sm_to_modify += body_sm_list;
     bodies_sm_to_modify += neighbor_imprint_list;
     push_vg_attributes_before_modify( bodies_sm_to_modify );
     bodies_to_modify += webcut_body_list;
     bodies_to_modify += neighboring_bodies;
     get_merged_curve_and_surface_ids( bodies_to_modify, merged_surface_ids, merged_curve_ids );
   }

   BodySM* cutting_tool_ptr = NULL;
   CubitStatus stat = prepare_surface_sweep(body_sm_list,surfaces_to_sweep,
                           sweep_axis,false,false,false,
                           up_to_next,stop_surface, NULL, cutting_tool_ptr, &point, &angle);
   if(stat == CUBIT_SUCCESS )
   {
     stat = gme->webcut( body_sm_list, cutting_tool_ptr, neighbor_imprint_list,
                         new_sms, imprint_type, preview );

     // Delete the BodySM that was created to be used as a tool
     gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;
   }

   // if not previewing do the rest of the creation
   if (!preview)
   {
     restore_vg_after_modify(new_sms, bodies_to_modify, gme);
     remove_pushed_attributes(new_sms, bodies_to_modify);
     stat = finish_webcut(webcut_body_list, new_sms, merge, stat, new_bodies,
                          &merged_surface_ids, &merged_curve_ids );
     do_attribute_cleanup();

     if( CubitUndo::get_undo_enabled() ) 
     {
       if( stat == CUBIT_SUCCESS )
         CubitUndo::note_result_bodies( new_bodies );
       else
         CubitUndo::remove_last_undo();
     }
   }

   return stat;
}

CubitStatus GeometryModifyTool::webcut_with_sweep_curves_rotated(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefEdge*> &tool_curves,
                            CubitVector &point,
                            CubitVector &sweep_axis,
                            double angle,
                            RefFace* stop_surf,
                            DLIList<Body*> &new_bodies,
                            DLIList<Body*> &neighboring_bodies,
                            ImprintType imprint_type,
                            CubitBoolean merge,
                            CubitBoolean preview)
{

   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);
   DLIList<Curve*> curves_to_sweep;
   if (gme)
   {
     DLIList<RefFace*> dummy1;
     DLIList<Surface*> dummy2;
     GeometryModifyEngine* surf_gme = common_modify_engine( dummy1, tool_curves,
                                                            dummy2, curves_to_sweep );
     if( gme != surf_gme )
       PRINT_ERROR("All geometry not from the same modeling engine.\n");
   }

   Surface *stop_surface = NULL;
   if( stop_surf )
   {
     GeometryModifyEngine *tmp_gme = get_engine(stop_surf);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Stop surface %d is not of same modeling engine as other input geometry.\n",
                    stop_surf->id() );
       return CUBIT_FAILURE;
     }

     stop_surface = stop_surf->get_surface_ptr();
   }

   DLIList<int> merged_surface_ids;
   DLIList<int> merged_curve_ids;
   DLIList<BodySM*> neighbor_imprint_list;
   DLIList<Body*> bodies_to_modify;

   if (!preview)
   {
     if( CubitUndo::get_undo_enabled() )
     {
       DLIList<Body*> bodies_to_save;
       bodies_to_save += webcut_body_list;
       bodies_to_save += neighboring_bodies;
       CubitUndo::save_state_with_cubit_file( bodies_to_save );
     }

     int i;
     for( i=neighboring_bodies.size(); i--; )
     {
       Body *neighbor_body = neighboring_bodies.get_and_step();
       BodySM *tmp_body = neighbor_body->get_body_sm_ptr(); 
       GeometryModifyEngine *neighbor_gme = get_engine( tmp_body );
       
       if( gme == neighbor_gme )
         neighbor_imprint_list.append( tmp_body ); 
     }

     do_attribute_setup();
     DLIList<BodySM*> bodies_sm_to_modify;
     bodies_sm_to_modify += body_sm_list;
     bodies_sm_to_modify += neighbor_imprint_list;
     push_vg_attributes_before_modify( bodies_sm_to_modify );
     bodies_to_modify += webcut_body_list;
     bodies_to_modify += neighboring_bodies;
     get_merged_curve_and_surface_ids( bodies_to_modify, merged_surface_ids, merged_curve_ids );
   }


   //sweep the curves.
   DLIList<GeometryEntity*> ref_ent_list;
   for(int i = 0; i < curves_to_sweep.size(); i++)
     ref_ent_list.append((GeometryEntity*)(curves_to_sweep.get_and_step())); 

   DLIList<BodySM*> swept_bodies;
   CubitStatus stat = gme->sweep_rotational(ref_ent_list,swept_bodies,point,
                                      sweep_axis,angle,0, 0.0, 0,false,false,
                                      false,stop_surface);
   if(stat == CUBIT_FAILURE  && swept_bodies.size() == 0)
     return stat;

   //stitch faces together
   BodySM* cutting_tool_ptr = NULL;
   stat = gme->stitch_surfs(swept_bodies, cutting_tool_ptr);

   if(stat == CUBIT_FAILURE)
     {
       //delete all swept faces
       for(int i = 0; i < swept_bodies.size(); i++)
         gme->get_gqe()->delete_solid_model_entities(swept_bodies.get_and_step()) ;
       return stat;
     }
   else
   {
     stat = gme->webcut( body_sm_list, cutting_tool_ptr, neighbor_imprint_list,
                         new_sms, imprint_type, preview);
    
     // Delete the BodySM that was created to be used as a tool
     gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;
   }  

   if (!preview)
   {
     restore_vg_after_modify(new_sms, bodies_to_modify, gme);
     remove_pushed_attributes(new_sms, bodies_to_modify);
     stat = finish_webcut(webcut_body_list, new_sms, merge, stat, new_bodies,
                          &merged_surface_ids, &merged_curve_ids );
     do_attribute_cleanup();

     if( CubitUndo::get_undo_enabled() ) 
     {
       if( stat == CUBIT_SUCCESS )
         CubitUndo::note_result_bodies( new_bodies );
       else
         CubitUndo::remove_last_undo();
     }
   }

   return stat;
}

CubitStatus GeometryModifyTool::webcut_with_sweep_curves(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefEdge*> &tool_curves,
                            const CubitVector sweep_vector,
                            bool through_all,
                            RefFace *stop_surf,
                            RefEdge* edge_to_sweep_along,
                            DLIList<Body*> &new_bodies,
                            DLIList<Body*> &neighboring_bodies,
                            ImprintType imprint_type,
                            CubitBoolean merge,
                            CubitBoolean preview)
{
   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);
   DLIList<Curve*> curves_to_sweep;
   if (gme)
   {
     DLIList<RefFace*> dummy1;
     DLIList<Surface*> dummy2;
     GeometryModifyEngine* surf_gme = common_modify_engine( dummy1, tool_curves,
                                                            dummy2, curves_to_sweep );
     if( gme != surf_gme ){
       PRINT_ERROR("All geometry not from the same modeling engine.\n");
       return CUBIT_FAILURE;
     }
   }
   else{
     PRINT_ERROR("All volumes not from the same modeling engine.\n");
     return CUBIT_FAILURE;
   }

   Surface *stop_surface = NULL;
   if( stop_surf )
   {
     GeometryModifyEngine *tmp_gme = get_engine(stop_surf);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Stop surface %d is not of same modeling engine as other input geometry.\n",
                    stop_surf->id() );
       return CUBIT_FAILURE;
     }

     stop_surface = stop_surf->get_surface_ptr();
   }

   Curve *curve_to_sweep_along = NULL;
   if( edge_to_sweep_along )
   {
     GeometryModifyEngine *tmp_gme = get_engine(edge_to_sweep_along);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Curve %d does is not of same modeling as other input geometry.\n",
                    edge_to_sweep_along->id() );
       return CUBIT_FAILURE;
     }

     curve_to_sweep_along = edge_to_sweep_along->get_curve_ptr();
   }

   DLIList<int> merged_surface_ids;
   DLIList<int> merged_curve_ids;
   DLIList<BodySM*> neighbor_imprint_list;
   DLIList<Body*> bodies_to_modify;

   if (!preview)
   {
     if( CubitUndo::get_undo_enabled() )
     {
       DLIList<Body*> bodies_to_save;
       bodies_to_save += webcut_body_list;
       bodies_to_save += neighboring_bodies;
       CubitUndo::save_state_with_cubit_file( bodies_to_save );
     }

     int i;
     for( i=neighboring_bodies.size(); i--; )
     {
       Body *neighbor_body = neighboring_bodies.get_and_step();
       BodySM *tmp_body = neighbor_body->get_body_sm_ptr(); 
       GeometryModifyEngine *neighbor_gme = get_engine( tmp_body );
       
       if( gme == neighbor_gme )
         neighbor_imprint_list.append( tmp_body ); 
     }

     do_attribute_setup();
     DLIList<BodySM*> bodies_sm_to_modify;
     bodies_sm_to_modify += body_sm_list;
     bodies_sm_to_modify += neighbor_imprint_list;
     push_vg_attributes_before_modify( bodies_sm_to_modify );
     bodies_to_modify += webcut_body_list;
     bodies_to_modify += neighboring_bodies;
     get_merged_curve_and_surface_ids( bodies_to_modify, merged_surface_ids, merged_curve_ids );
   }

   CubitVector tmp_sweep_vector = sweep_vector;
   //get model bbox info...will scale sweep vector by its diagonal
   //so that we go far enough
   if( through_all || stop_surf )
     {
       CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
       tmp_sweep_vector.normalize();
       tmp_sweep_vector*=(2*bounding_box.diagonal().length());
     }

   DLIList<GeometryEntity*> ref_ent_list;
   for(int i = 0; i < curves_to_sweep.size(); i++)
     ref_ent_list.append((GeometryEntity*)(curves_to_sweep.get_and_step()));
   DLIList<BodySM*> swept_bodies;
   CubitStatus stat;

   if( curve_to_sweep_along )
   {
     DLIList<Curve*> curves_to_sweep_along;
     curves_to_sweep_along.append(curve_to_sweep_along);
     stat = gme->sweep_along_curve(ref_ent_list, swept_bodies,
                              curves_to_sweep_along,0.0,0,false,stop_surface);

     if (!stat && swept_bodies.size() == 0)
       return stat;
   }

   else
   {
     stat = gme->sweep_translational(ref_ent_list, swept_bodies,
                            tmp_sweep_vector,0.0,0, false, false, stop_surface);

     if (!stat && swept_bodies.size() == 0)
       return stat;
   }

   //stitch faces together
   BodySM* cutting_tool_ptr = NULL;
   stat = gme->stitch_surfs(swept_bodies, cutting_tool_ptr);

   if(stat == CUBIT_FAILURE)
   {
       //delete all swept faces
       for(int i = 0; i < swept_bodies.size(); i++)
         gme->get_gqe()->delete_solid_model_entities(swept_bodies.get_and_step()) ;
       return stat;
   }
   else 
   {
       stat = gme->webcut( body_sm_list, cutting_tool_ptr, neighbor_imprint_list,
                           new_sms, imprint_type, preview );

       // Delete the BodySM that was created to be used as a tool
       gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;
   }
   
   if (!preview)
   {
     restore_vg_after_modify(new_sms, bodies_to_modify, gme);
     remove_pushed_attributes(new_sms, bodies_to_modify);
     stat = finish_webcut(webcut_body_list, new_sms, merge, stat, new_bodies,
                          &merged_surface_ids, &merged_curve_ids );
     do_attribute_cleanup();

     if( CubitUndo::get_undo_enabled() ) 
     {
       if( stat == CUBIT_SUCCESS )
         CubitUndo::note_result_bodies( new_bodies );
       else
         CubitUndo::remove_last_undo();
     }
   }

   return stat;
}

CubitStatus GeometryModifyTool::webcut_with_sweep_surfaces(
                            DLIList<Body*> &webcut_body_list,
                            DLIList<RefFace*> &tool_faces,
                            const CubitVector sweep_vector,
                            bool sweep_perp,
                            bool through_all,
                            bool outward,
                            bool up_to_next,
                            RefFace* stop_surf,
                            RefEdge* edge_to_sweep_along,
                            DLIList<Body*> &new_bodies,
                            DLIList<Body*> &neighboring_bodies,
                            ImprintType imprint_type,
                            CubitBoolean merge,
                            CubitBoolean preview)
{
   if (!okay_to_modify( webcut_body_list, "WEBCUT" ))
     return CUBIT_FAILURE;

   DLIList<BodySM*> body_sm_list(webcut_body_list.size()), new_sms;
   GeometryModifyEngine* gme = common_modify_engine(webcut_body_list, body_sm_list);
   DLIList<Surface*> surfaces_to_sweep;
   if (gme)
   {
     GeometryModifyEngine* surf_gme = common_modify_engine( tool_faces, surfaces_to_sweep );
     if( gme != surf_gme ){
       PRINT_ERROR("All geometry not from the same modeling engine.\n");
       return CUBIT_FAILURE;
     }
   }
   else{
     PRINT_ERROR("All volumes not from the same modeling engine.\n");
     return CUBIT_FAILURE;
   }

   Surface *stop_surface = NULL;
   if( stop_surf )
   {
     GeometryModifyEngine *tmp_gme = get_engine(stop_surf);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Stop surface %d is not of same modeling engine as other input geometry.\n",
                    stop_surf->id() );
       return CUBIT_FAILURE;
     }

     stop_surface = stop_surf->get_surface_ptr();
   }

   Curve *curve_to_sweep_along = NULL;
   if( edge_to_sweep_along )
   {
     GeometryModifyEngine *tmp_gme = get_engine(edge_to_sweep_along);

     if( gme != tmp_gme )
     {
       PRINT_ERROR("Curve %d does is not of same modeling as other input geometry.\n",
                    edge_to_sweep_along->id() );
       return CUBIT_FAILURE;
     }

     //make sure that the curve is not part of the surface(s) being swept
     DLIList<RefFace*> faces_of_edge;
     edge_to_sweep_along->ref_faces( faces_of_edge );
      
     int kk;
     for( kk=faces_of_edge.size(); kk--; ) 
     {
       if( tool_faces.is_in_list( faces_of_edge.get_and_step() ) )
       {
         faces_of_edge.back();
         PRINT_ERROR("Cannot perform sweep.  Curve %d is in Surface %d\n",
                      edge_to_sweep_along->id(), faces_of_edge.get()->id() );
         return CUBIT_FAILURE;
       }
     }

     curve_to_sweep_along = edge_to_sweep_along->get_curve_ptr();
   }

   DLIList<BodySM*> neighbor_imprint_list;
   DLIList<int> merged_surface_ids;
   DLIList<int> merged_curve_ids;
   DLIList<Body*> bodies_to_modify;

   if (!preview)
   {
     if( CubitUndo::get_undo_enabled() )
     {
       DLIList<Body*> bodies_to_save;
       bodies_to_save += webcut_body_list;
       bodies_to_save += neighboring_bodies;
       CubitUndo::save_state_with_cubit_file( bodies_to_save );
     }

     int i;
     for( i=neighboring_bodies.size(); i--; )
     {
       Body *neighbor_body = neighboring_bodies.get_and_step();
       BodySM *tmp_body = neighbor_body->get_body_sm_ptr(); 
       GeometryModifyEngine *neighbor_gme = get_engine( tmp_body );
       
       if( gme == neighbor_gme )
         neighbor_imprint_list.append( tmp_body ); 
     }

     do_attribute_setup();
     DLIList<BodySM*> bodies_sm_to_modify;
     bodies_sm_to_modify += body_sm_list;
     bodies_sm_to_modify += neighbor_imprint_list;
     push_vg_attributes_before_modify( bodies_sm_to_modify );
     bodies_to_modify += webcut_body_list;
     bodies_to_modify += neighboring_bodies;
     get_merged_curve_and_surface_ids( bodies_to_modify, merged_surface_ids, merged_curve_ids );
   }

   BodySM* cutting_tool_ptr = NULL;
   CubitStatus stat = prepare_surface_sweep(
      body_sm_list, surfaces_to_sweep, sweep_vector, sweep_perp, through_all, outward,
      up_to_next, stop_surface, curve_to_sweep_along, cutting_tool_ptr );
   if (stat == CUBIT_SUCCESS)
   {
      stat = gme->webcut( body_sm_list, cutting_tool_ptr, neighbor_imprint_list,
                          new_sms, imprint_type, preview );

      // Delete the BodySM that was created to be used as a tool
      gme->get_gqe()->delete_solid_model_entities(cutting_tool_ptr) ;
   }
   
   if (!preview)
   {
     restore_vg_after_modify(new_sms, bodies_to_modify, gme);
     remove_pushed_attributes(new_sms, bodies_to_modify);
     stat = finish_webcut(webcut_body_list, new_sms, merge, stat, new_bodies,
                          &merged_surface_ids, &merged_curve_ids );
     do_attribute_cleanup();

     if( CubitUndo::get_undo_enabled() ) 
     {
       if( stat  == CUBIT_SUCCESS )
         CubitUndo::note_result_bodies( new_bodies );
       else
         CubitUndo::remove_last_undo();
     }
   }

   return stat;
}

CubitStatus GeometryModifyTool::split_free_curve( RefEdge *ref_edge,
                                                  CubitVector &split_location )
{
  TopologyBridge* bridge = 0;
  GeometryModifyEngine* gme_ptr = get_engine(ref_edge, &bridge);
  Curve *curve = dynamic_cast<Curve*>(bridge);
  
  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<RefEdge*> tmp_ents(1);
    tmp_ents.append( ref_edge );
    CubitUndo::save_state_with_cubit_file( tmp_ents );
  }

  DLIList<Curve*> new_curves;
  gme_ptr->split_free_curve( curve, split_location, new_curves );
 
  if (!new_curves.size())
  {
    if( CubitUndo::get_undo_enabled() )
      CubitUndo::remove_last_undo();
    return CUBIT_FAILURE;
  }

  GeometryQueryTool::instance()->delete_RefEdge( ref_edge );
  
  int i;
  for( i=0; i<new_curves.size(); i++ )
  {
    Curve *new_curve = new_curves.get_and_step();
    RefEdge* new_edge = GeometryQueryTool::instance()->make_free_RefEdge(new_curve);
    if( CubitUndo::get_undo_enabled() )
      CubitUndo::note_result_entity( new_edge );
  }

  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : split a multiple volume body into several bodies
//                 each containing a single volume.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 09/26/97
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::split_body( Body *body_ptr,
                                            DLIList<Body*> &new_bodies ) const
{
  int i;
  DLIList<RefVolume*> ref_vols;
  body_ptr->ref_volumes(ref_vols);
  if ( ref_vols.size() < 2 )
  {
    //no need to split...
    new_bodies.append(body_ptr);
    return CUBIT_SUCCESS;
  }

  DLIList<Body*> b_list;
  b_list.append(body_ptr);
  if (!okay_to_modify( b_list, "SPLIT" ))
    return CUBIT_FAILURE;

  // Call the owning GeometryModifyEngine to split the body
  // so that there is one volume per body.
  BodySM* bodysm_ptr = body_ptr->get_body_sm_ptr();
  GeometryModifyEngine* engine_ptr = get_engine(bodysm_ptr);
  if (!engine_ptr) {
    PRINT_ERROR("There is no modify engine available for this volume."
      " Volume cannot be split.\n");
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<Body*> bodies(1);
    bodies.append( body_ptr );
    CubitUndo::save_state_with_cubit_file( bodies );
  }

  DLIList<BodySM*> new_sm_list;
  CubitStatus stat = engine_ptr->split_body( bodysm_ptr, new_sm_list );
  if ( new_sm_list.size() == 0 )
  {
    PRINT_ERROR("failed in splitting volumes, orginal was lost.\n");
    GeometryQueryTool::instance()->destroy_dead_entity(body_ptr);
    return CUBIT_FAILURE;
  }

  bodysm_ptr = body_ptr->get_body_sm_ptr();
  if (bodysm_ptr)
  {
    remove_dead_entity_names(body_ptr);
    body_ptr = GeometryQueryTool::instance()->make_Body(bodysm_ptr);
  }
  else
  {
    GeometryQueryTool::instance()->destroy_dead_entity(body_ptr);
  }

  new_sm_list.reset();
  for (i = new_sm_list.size(); i--; )
  {
    bodysm_ptr = new_sm_list.get_and_step();
    Body* body = GeometryQueryTool::instance()->make_Body(bodysm_ptr);
    new_bodies.append(body);
  }

  if ( !GeometryModifyTool::instance()->get_new_ids())
  {
    //Now reuse the body ids.
    DLIList<RefVolume*> new_ref_vols;
    for ( int ii = new_bodies.size(); ii > 0; ii-- )
    {
      Body *temp_body = new_bodies.get_and_step();
      new_ref_vols.clean_out();
      temp_body->ref_volumes(new_ref_vols);
      int vol_id = new_ref_vols.get()->id();
      if ( temp_body->id() != vol_id )
      {
        //Check to see if this id is being used by
        //some other body.
        if ( RefEntityFactory::instance()->get_body(vol_id) == NULL )
        {
          temp_body->set_id(vol_id);
        }
        //if it is in use, then we shouldn't mess around with it...
      }
    }
  }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::note_result_bodies( new_bodies );

  GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  return stat;
}

//-------------------------------------------------------------------------
// Purpose       : Separate surfaces from a sheet body into separate bodies.
//                 Connected surfaces will remain united but be placed in
//                 a new body.
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 02/23/2008
//-------------------------------------------------------------------------
CubitStatus
GeometryModifyTool::separate_surfaces( DLIList<RefFace*> &ref_face_list,
                                       DLIList<Body*> &new_bodies )
{
  int i;
  RefFace *ref_face_ptr;

  for( i=ref_face_list.size(); i--; )
  {
    ref_face_ptr = ref_face_list.get();

    // Check for no body
    DLIList<Body*> body_list;
    ref_face_ptr->bodies( body_list );
    if( body_list.size()==0 )
    {
      PRINT_ERROR( "Surface %d is not contained within a parent body.\n"
        "       It cannot be separated.\n", ref_face_ptr->id() );
      return CUBIT_FAILURE;
    }
  }

  // Check for virtual geometry
  DLIList<RefEntity*> ref_ent_list;
  CAST_LIST_TO_PARENT(ref_face_list, ref_ent_list);
  if ( GeometryQueryTool::instance()->contains_intermediate_geometry(ref_ent_list) )
  {
    PRINT_ERROR("SEPARATING surfaces containing virtual geometry is not\n"
      "       allowed. Delete virtual geometry on these surfaces\n"
      "       before operation.\n" );
    return CUBIT_FAILURE;
  }

  // Prep for undo
  DLIList<Body*> body_list;
  for( i=ref_face_list.size(); i--; )
    ref_face_list.get_and_step()->bodies( body_list );
  body_list.uniquify_unordered();

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state_with_cubit_file( body_list );

  // Handle surfaces from different modify engines.  Copy the input list since
  // we will be pulling RefFaces out of it.
  DLIList<RefFace*> copied_ref_face_list = ref_face_list;

  // Keep track of errors and success
  int error_occurred = 0;
  int successful_case = 0;

  GeometryModifyEngine *gme_ptr;
  while( copied_ref_face_list.size() )
  {
    // Need to send in Surfaces
    DLIList<Surface*> surface_list;
    DLIList<RefFace*> gme_face_list;

    gme_ptr = pull_common_surfs( copied_ref_face_list, gme_face_list,
                                 surface_list );
    if( !gme_ptr )
      return CUBIT_FAILURE;

    // Get the owning bodies of the faces...needed for finish_sm_op
    DLIList<Body*> gme_body_list;
    for( i=gme_face_list.size(); i--; )
      gme_face_list.get_and_step()->bodies( gme_body_list );
    gme_body_list.uniquify_unordered();
    
    DLIList<BodySM*> new_sm_list;
    DLIList<Body*> new_body_list;
    if( gme_ptr->separate_surfaces( surface_list, new_sm_list ) == CUBIT_FAILURE ||
        finish_sm_op(gme_body_list, new_sm_list, new_body_list ) == CUBIT_FAILURE )
    {
      error_occurred++;
      continue;
    }

    new_bodies += new_body_list;
    successful_case++;
  }

  if( error_occurred && !successful_case )
    return CUBIT_FAILURE;

  // Following is copied from split_body - to keep same behavior. When all
  // surfaces are separated from a given body, separate_surfaces just separates
  // the body. Without the code below we get new body ids.
  if ( !GeometryModifyTool::instance()->get_new_ids())
  {
    //Now reuse the body ids.
    DLIList<RefVolume*> new_ref_vols;
    for ( int ii = new_bodies.size(); ii > 0; ii-- )
    {
      Body *temp_body = new_bodies.get_and_step();
      new_ref_vols.clean_out();
      temp_body->ref_volumes(new_ref_vols);
      int vol_id = new_ref_vols.get()->id();
      if ( temp_body->id() != vol_id )
      {
        //Check to see if this id is being used by
        //some other body.
        if ( RefEntityFactory::instance()->get_body(vol_id) == NULL )
        {
          temp_body->set_id(vol_id);
        }
        //if it is in use, then we shouldn't mess around with it...
      }
    }
  }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::note_result_bodies( new_bodies );

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Reverse a body
//
// Special Notes : Moved from Body
//
// Creator       : Jason Kraftcheck
//
// Creation Date :
//-------------------------------------------------------------------------
CubitStatus GeometryModifyTool::reverse( DLIList<Body*> &body_list )
{
   DLIList<Body*> reversed_bodies;
   int i;
   for( i=body_list.size(); i--; )
   {
     Body *body = body_list.get_and_step();
     BodySM* body_sm = body->get_body_sm_ptr();
     GeometryModifyEngine* gme = get_engine( body_sm );
     if (!gme) {
        PRINT_ERROR("Body %d does not have a modify engine.\n", body->id());
        continue;
     }

     CubitStatus stat = gme->reverse_body( body_sm );

     if ( CUBIT_SUCCESS != stat )
     {
       PRINT_ERROR("Reverse failed.\n");
       continue;
     }
     else
     {
       reversed_bodies.append( body );
       GeometryQueryTool::instance()->make_Body( body_sm );
       body->notify_all_observers( GEOMETRY_MODIFIED );
       PRINT_INFO("Reversed body %d.\n", body->id());
       continue;
     }
   }

   if( CubitUndo::get_undo_enabled() )
   {
     CubitString undo_command("reverse body ");
     for( i=reversed_bodies.size(); i--; )
     {
       undo_command += reversed_bodies.get_and_step()->id();
       undo_command += " ";
     }
     CubitUndo::set_undo_by_command( undo_command );
   }

   return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : reverse
// Member Type: PUBLIC
// Description: Reverse given surfaces (flip normals)
// Author     : Steve Storm (CAT)
// Date       : 4/3/2007
//===============================================================================
CubitStatus GeometryModifyTool::reverse( DLIList<RefFace*> &ref_face_list )
{
  // Check for virtual geometry
  DLIList<RefEntity*> ref_ent_list;
  CAST_LIST_TO_PARENT(ref_face_list, ref_ent_list);
  if( GeometryQueryTool::instance()->contains_intermediate_geometry(ref_ent_list) )
  {
    PRINT_ERROR("REVERSING surfaces containing virtual geometry is not\n"
      "       allowed. Delete virtual geometry on surfaces before\n"
      "       operation.\n" );
    return CUBIT_FAILURE;
  }

  // Get the owning bodies of the faces
  int i;
  DLIList<Body*> body_list;
  for( i=ref_face_list.size(); i--; )
    ref_face_list.get_and_step()->bodies( body_list );
  body_list.uniquify_unordered();

  // Handle surfaces from different modify engines.  Copy the input list since
  // we will be pulling RefFaces out of it.
  DLIList<RefFace*> copied_ref_face_list = ref_face_list;

  // Keep track of overall errors and successes
  int error_occurred = 0;
  int successful_case = 0;
  GeometryModifyEngine *gme_ptr;

  while( copied_ref_face_list.size() )
  {
    // Need to send in Surfaces
    DLIList<Surface*> surface_list;
    DLIList<RefFace*> common_face_list;

    gme_ptr = pull_common_surfs( copied_ref_face_list, common_face_list,
                                 surface_list );
    if( !gme_ptr )
      return CUBIT_FAILURE;

    if( gme_ptr->flip_normals( surface_list ) == CUBIT_FAILURE )
      error_occurred = 1;
    else
      successful_case = 1;
  }

  if( error_occurred && !successful_case )
    return CUBIT_FAILURE;

  // Update the impacted bodies
  Body *body_ptr;
  BodySM* body_sm_ptr;
  DLIList<BodySM*> new_sm_list;
  for( i=body_list.size(); i--; )
  {
    body_ptr = body_list.get_and_step();
    body_ptr->notify_all_observers( GEOMETRY_MODIFIED );
    body_sm_ptr = body_ptr->get_body_sm_ptr();
    new_sm_list.append( body_sm_ptr );
  }

  DLIList<Body*> new_body_list;
  return finish_sm_op(body_list, new_sm_list, new_body_list );
}

//-------------------------------------------------------------------------
// Purpose       : split a periodic surface.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 09/26/97
//-------------------------------------------------------------------------

CubitStatus GeometryModifyTool::split_periodic(Body *body_ptr,
                                               Body *&new_body_ptr )
{
   DLIList<Body*> b_list(1);
   b_list.append(body_ptr);
   if (!okay_to_modify( b_list, "SPLIT" ))
     return CUBIT_FAILURE;

   BodySM* body_sm = body_ptr->get_body_sm_ptr();
   GeometryModifyEngine* gme = get_engine(body_sm);
   if (!gme) {
      PRINT_ERROR("Volume %d does not have a modify engine.\n", body_ptr->id());
      return CUBIT_FAILURE;
   }

   BodySM* new_sm = 0;

   if (CubitUndo::get_undo_enabled())
     CubitUndo::save_state_with_cubit_file(b_list);

   DLIList<BodySM*> body_sms;
   body_sms.append(body_sm);

   do_attribute_setup();
   push_vg_attributes_before_modify(body_sms);

   // Call the default GeometryModifyEngine to create a new Point
   CubitStatus stat = gme->split_periodic( body_sm, new_sm );

   DLIList<BodySM*> new_bodysm_list;
   if(new_sm)
     new_bodysm_list.append(new_sm);
   DLIList<Body*> old_body_list;
   old_body_list.append(body_ptr);
   restore_vg_after_modify(new_bodysm_list, old_body_list, gme);
   remove_pushed_attributes(new_bodysm_list, old_body_list);
   do_attribute_cleanup();

   update_body(body_ptr);

   new_body_ptr = 0;
   if (new_sm)
     new_body_ptr = GeometryQueryTool::instance()->make_Body(new_sm);

   if (new_body_ptr)
   {
     if (CubitUndo::get_undo_enabled())
       CubitUndo::note_result_body(new_body_ptr);
   }
   else
   {
     if (CubitUndo::get_undo_enabled())
       CubitUndo::remove_last_undo();
   }

   GeometryQueryTool::instance()->cleanout_deactivated_geometry();

   return stat;
}

//===============================================================================
// Function   : split_surface
// Member Type: PUBLIC
// Description: Split a single surface into one or more pieces
// Author     : Steve Storm (CAT)
// Date       : 07/04
//===============================================================================
CubitStatus
GeometryModifyTool::split_surface( RefFace *ref_face_ptr,
                                   DLIList<CubitVector*> &locations,
                                   DLIList<DLIList<CubitVector*>*> &vec_lists,
                                   CubitBoolean preview_flg,
                                   CubitBoolean create_ref_edges_flg,
                                   CubitBoolean clear_previous_previews )
{
  // Check for virtual geometry
  DLIList<RefFace*> ref_face_list;
  ref_face_list.append( ref_face_ptr );

  DLIList<Body*> b_list;
  Body* body = ref_face_ptr->body();
  if (body)
	  b_list.append(body);
  else
  {
	  PRINT_ERROR( "Surface %d is not contained within a parent body.\n"
		  "      It cannot be split.\n", ref_face_ptr->id() );
	  return CUBIT_FAILURE;
  }

  if (!okay_to_modify( b_list, "SPLIT_SURFACE" ))
    return CUBIT_FAILURE;

  int i;
  DLIList<Body*> old_body_list;
  for( i=ref_face_list.size(); i--; )
  {
    ref_face_ptr = ref_face_list.get_and_step();

    DLIList<Body*> body_list;
    ref_face_ptr->bodies( body_list );
    old_body_list.merge_unique( body_list );
  }

  //bad geom with no body -- dont try to imprint this...
  //quick and dirty fix by (aga@cat|1/7/04)
  if( old_body_list.size() < 1 )
  {
    PRINT_ERROR( "Surface %d is not contained within a parent body\n."
      "      It cannot be split.\n", ref_face_ptr->id() );
    return CUBIT_FAILURE;
  }

  SplitSurfaceTool sst;
  if( preview_flg == CUBIT_TRUE )
    return sst.preview( ref_face_ptr, locations, vec_lists, create_ref_edges_flg,
                        clear_previous_previews );
  else
    return sst.split_surface( ref_face_ptr, vec_lists );
}

//===============================================================================
// Function   : split_surface
// Member Type: PUBLIC
// Description: Split multiple surfaces
// Author     : Greg Schebler (CAT)
// Date       : 06/08
//===============================================================================

CubitStatus
GeometryModifyTool::split_surface( DLIList<RefFace*> &ref_face_list,
                                  DLIList<CubitVector*> &locations,
                                  DLIList<DLIList<DLIList<CubitVector*>*>*> &list_of_vec_lists,
                                  CubitBoolean preview_flg,
                                  CubitBoolean create_ref_edges_flg,
                                  CubitBoolean clear_previous_previews )
{
   // Check for virtual geometry
   DLIList<Body*> b_list;
   int gg;
   for( gg = ref_face_list.size() ; gg > 0 ; gg--)
   {
      RefFace* ref_face_ptr = ref_face_list.get_and_step();
      Body* body = ref_face_ptr->body();
      if (body)
         b_list.append(body);
      else
      {
         PRINT_ERROR( "Surface %d is not contained within a parent body.\n"
            "      It cannot be split.\n", ref_face_ptr->id() );
         return CUBIT_FAILURE;
      }
   }
   if (!okay_to_modify( b_list, "SPLIT_SURFACE" ))
      return CUBIT_FAILURE;

   int ii;

   for( ii = ref_face_list.size(); ii > 0 ; ii--)
   {
      DLIList<Body*> old_body_list;
      RefFace* ref_face_ptr = ref_face_list.get_and_step();

      DLIList<Body*> body_list;
      ref_face_ptr->bodies( body_list );
      old_body_list.merge_unique( body_list );

      if( old_body_list.size() < 1 )
      {
         PRINT_ERROR( "Surface %d is not contained within a parent body\n."
            "      It cannot be split.\n", ref_face_ptr->id() );
         return CUBIT_FAILURE;
      }
   }

   SplitSurfaceTool sst;
   if( preview_flg == CUBIT_TRUE )
      return sst.preview( ref_face_list, locations, list_of_vec_lists, create_ref_edges_flg,
      clear_previous_previews );
   else
      return sst.split_surface( ref_face_list, list_of_vec_lists );
}

//===============================================================================
// Function   : split_surfaces_extend
// Member Type: PUBLIC
// Description: Split surfaces by extending hardline curves on the surface
// Author     : Steve Storm (CAT)
// Date       : 10/07
//===============================================================================
CubitStatus
GeometryModifyTool::split_surfaces_extend( DLIList<RefFace*> &ref_face_list,
                                           DLIList<RefVertex*> &ref_vertex_list,
                                           CubitBoolean preview_flg,
                                           CubitBoolean create_ref_edges_flg )
{
  int i;
  RefFace *ref_face_ptr;

  for( i=ref_face_list.size(); i--; )
  {
    ref_face_ptr = ref_face_list.get();

    // Check for no body
    DLIList<Body*> body_list;
    ref_face_ptr->bodies( body_list );
    if( body_list.size()==0 )
    {
      PRINT_ERROR( "Surface %d is not contained within a parent body.\n"
        "       It cannot be split.\n", ref_face_ptr->id() );
      return CUBIT_FAILURE;
    }
  }

  // Check for virtual geometry
  DLIList<RefEntity*> ref_ent_list;
  CAST_LIST_TO_PARENT(ref_face_list, ref_ent_list);
  if ( GeometryQueryTool::instance()->contains_intermediate_geometry(ref_ent_list) )
  {
    PRINT_ERROR("SPLITTING surfaces containing virtual geometry is not\n"
      "       allowed. Delete virtual geometry on these surfaces\n"
      "       before operation.\n" );
    return CUBIT_FAILURE;
  }

  // Make sure all surfaces are from same geometry engine
  if ( !same_modify_engine(ref_ent_list, CUBIT_TRUE) )
  {
    PRINT_ERROR("Performing SPLIT with surfaces containing geometry from\n"
      "different modeling engines is not allowed.\n"
      "Delete uncommon geometry on these surfaces before operation.\n\n");
    return CUBIT_FAILURE;
  }

  SplitSurfaceTool sst;
  return sst.split_surfaces_extend( ref_face_list, ref_vertex_list,
    preview_flg, create_ref_edges_flg );
}

//===============================================================================
// Function   : split_surfaces
// Member Type: PUBLIC
// Description: Split a chain of surfaces into one or more pieces
// Author     : Steve Storm (CAT)
// Date       : 01/04
//===============================================================================
CubitStatus
GeometryModifyTool::split_surfaces( DLIList<RefFace*> &ref_face_list,
                                    int num_segs,
                                    double fraction,
                                    double distance,
                                    RefEdge *from_curve_ptr,
                                    DLIList<RefVertex*> &corner_vertex_list,
                                    DLIList<RefVertex*> &through_vertex_list,
                                    RefEdge *curve_dir_ptr,
                                    CubitBoolean preview_flg,
                                    CubitBoolean create_ref_edges_flg )
{
  // Get parent bodies - all surfs must be from same body
  int i;
  DLIList<Body*> old_body_list;
  RefFace *ref_face_ptr;
  ref_face_list.reset();
  for( i=ref_face_list.size(); i--; )
  {
    ref_face_ptr = ref_face_list.get_and_step();

    DLIList<Body*> body_list;
    ref_face_ptr->bodies( body_list );
    old_body_list.merge_unique( body_list );
  }

  if( old_body_list.size() > 1 )
  {
    PRINT_ERROR( "This operation requires all surfaces to be from the same volume\n" );
    // Note: this restriction could be pretty easily lifted by sorting the
    //       input lists and calling the SplitSurfaceTool separately for each set of
    //       surfaces on each body.
    return CUBIT_FAILURE;
  }

  //bad geom with no body -- dont try to imprint this...
  //quick and dirty fix by (aga@cat|1/7/04)
  if( old_body_list.size() < 1 )
  {
    PRINT_ERROR( "A surface is not contained within a parent body.\n"
      "       It cannot be split.\n");
    return CUBIT_FAILURE;
  }

  // Check for virtual geometry
  DLIList<RefEntity*> ref_ent_list;
  CAST_LIST_TO_PARENT(ref_face_list, ref_ent_list);
  if ( GeometryQueryTool::instance()->contains_intermediate_geometry(ref_ent_list) )
  {
    PRINT_ERROR("SPLITTING surfaces containing virtual geometry is not\n"
      "       allowed. Delete virtual geometry on these surfaces\n"
      "       before operation.\n" );
    return CUBIT_FAILURE;
  }

  // Make sure all surfaces are from same geometry engine
  if ( !same_modify_engine(ref_ent_list, CUBIT_TRUE) )
  {
    PRINT_ERROR("Performing SPLIT with surfaces containing geometry from\n"
      "different modeling engines is not allowed.\n"
      "Delete uncommon geometry on these surfaces before operation.\n\n");
    return CUBIT_FAILURE;
  }

  SplitSurfaceTool split_tool;
  return split_tool.split_surfaces( ref_face_list, num_segs, fraction,
    distance, from_curve_ptr, corner_vertex_list, through_vertex_list,
    curve_dir_ptr, preview_flg, create_ref_edges_flg );
}

//===============================================================================
// Function   : split_surfaces_offset
// Member Type: PUBLIC
// Description: Split a list of surface offset from a curve
// Author     : Sam Showman (CAT)
// Date       : 05/05
//===============================================================================
CubitStatus
GeometryModifyTool::split_surfaces_offset(DLIList<RefFace*> &ref_face_list,
                                          DLIList<RefEdge*> &edge_list,
                                          int num_segs,
                                          double distance,
                                          CubitBoolean partition_flg,
                                          CubitBoolean blunt_flg,
                                          CubitBoolean preview_flg,
                                          CubitBoolean create_ref_edges_flg)
{
  int i;
  DLIList<Body*> old_body_list;
  RefFace *ref_face_ptr;
  ref_face_list.reset();
  for( i=ref_face_list.size(); i--; )
  {
    ref_face_ptr = ref_face_list.get_and_step();

    DLIList<Body*> body_list;
    ref_face_ptr->bodies( body_list );
    old_body_list.merge_unique( body_list );
  }

  //check for geometry with no body
  if( old_body_list.size() < 1 )
  {
    PRINT_ERROR( "A surface is not contained within a parent body.\n"
      "       It cannot be split.\n");
    return CUBIT_FAILURE;
  }

  // Check for virtual geometry
  DLIList<RefEntity*> ref_ent_list;
  CAST_LIST_TO_PARENT(ref_face_list, ref_ent_list);
  if ( GeometryQueryTool::instance()->contains_intermediate_geometry(ref_ent_list) )
  {
    PRINT_ERROR("SPLITTING surfaces containing virtual geometry is not\n"
      "       allowed. Delete virtual geometry on these surfaces\n"
      "       before operation.\n" );
    return CUBIT_FAILURE;
  }

  // Make sure all surfaces are from same geometry engine
  // Make sure all curves are from same geometry engine
  CAST_LIST_TO_PARENT(edge_list, ref_ent_list);
  if ( !same_modify_engine(ref_ent_list, CUBIT_TRUE) )
  {
    PRINT_ERROR("Performing SPLIT with geometry from\n"
      "different modeling engines is not allowed.\n");
    return CUBIT_FAILURE;
  }

  OffsetSplitTool split_tool;
  return split_tool.split_surfaces_offset(
      ref_face_list,
      edge_list,
      num_segs,
      distance,
      partition_flg,
      blunt_flg,
      preview_flg,
      create_ref_edges_flg);
}

//===============================================================================
// Function   : auto_mid_surface
// Member Type: PUBLIC
// Description: Automatically midsurface a volume
// Author     : Sam Showman (CAT)
// Date       : 12/07
//===============================================================================
CubitStatus
GeometryModifyTool::auto_mid_surface(DLIList<Body*> &body_list_in,
                                     DLIList<Body*> &body_list_out,
                                     DLIList<Body*> &old_bodies_midsurfaced,
                                     DLIList<double> &thickness_list,
                                     double lower_tol,
                                     double upper_tol,
                                     CubitBoolean delete_midsurfaced,
 
                                     CubitBoolean preview)
{
	// Check for virtual geometry
	DLIList<RefEntity*> ref_ent_list;
	CAST_LIST_TO_PARENT(body_list_in, ref_ent_list);
	if ( GeometryQueryTool::instance()->contains_intermediate_geometry(ref_ent_list) )
	{
		PRINT_ERROR("Midsurfacing volumes containing virtual geometry is not\n"
			"       allowed. Delete virtual geometry on these surfaces\n"
			"       before operation.\n" );
		return CUBIT_FAILURE;
	}

	// Make sure all surfaces are from same geometry engine
	if ( !same_modify_engine(ref_ent_list, CUBIT_TRUE) )
	{
		PRINT_ERROR("Performing Midsurface with geometry from\n"
			"different modeling engines is not allowed.\n");
		return CUBIT_FAILURE;
	}

	AutoMidsurfaceTool mid_tool;
	DLIList<BodySM*> bodysm_list_out;

    if (CubitUndo::get_undo_enabled() && delete_midsurfaced && !preview)
        CubitUndo::save_state_with_cubit_file(body_list_in);

    CubitStatus result = mid_tool.midsurface(
        body_list_in,
        bodysm_list_out,
        old_bodies_midsurfaced,
        thickness_list,
        lower_tol,
        upper_tol,
        delete_midsurfaced,
        preview);


	if (result == CUBIT_SUCCESS &&
		bodysm_list_out.size() > 0 &&
		!preview)
	{
		if(CubitUndo::get_undo_enabled())
			CubitUndo::save_state();

		for( int i=0; i<bodysm_list_out.size(); i++ )
		{
			body_list_out.append(GeometryQueryTool::instance()->make_Body(bodysm_list_out[i]));
		}

		if( CubitUndo::get_undo_enabled() )
		{
			if( body_list_out.size() )
				CubitUndo::note_result_bodies( body_list_out );
			else
				CubitUndo::remove_last_undo();
		}
	}

	return result;
}

//-------------------------------------------------------------------------
// Purpose       : clean RefEntity
//
// Special Notes :
//
// Creator       : Lingyun Pan (CAT)
//
// Creation Date : 07/15/01
//-------------------------------------------------------------------------

CubitStatus
GeometryModifyTool::regularize_refentity(RefEntity *old_entity_ptr, Body *&new_body_ptr)
{
  BasicTopologyEntity* bte_ptr = dynamic_cast<BasicTopologyEntity*>(old_entity_ptr);
  if (!bte_ptr)
  {
    PRINT_ERROR("Invalid entity passed to regularize_refentity(..)\n");
    return CUBIT_FAILURE;
  }

  DLIList<Body*> body_list;
  bte_ptr->bodies(body_list);

  //ignore free entities
  if( body_list.size() == 0 )
  {
    PRINT_ERROR("%s %d is a free entity.  Cannot regularize it.\n", old_entity_ptr->class_name(), 
                                                                    old_entity_ptr->id() );
    new_body_ptr = NULL;
    return CUBIT_FAILURE;
  }

  if (!okay_to_modify( body_list, "REGULARIZE" ))
    return CUBIT_FAILURE;

  DLIList<TopologyBridge*> bridge_list;
  bte_ptr->bridge_manager()->get_bridge_list(bridge_list);
  bridge_list.reset();

  DLIList<BodySM*> new_sm_list;
  DLIList<Body*> new_body_list;
  CubitStatus stat = CUBIT_SUCCESS;

  do_attribute_setup();

  GeometryModifyEngine *save_gme = NULL;
  for (int i = bridge_list.size(); i--; )
  {
    TopologyBridge* bridge = bridge_list.get_and_step();
    GeometryEntity* geom_ptr = dynamic_cast<GeometryEntity*>(bridge);
    GeometryModifyEngine* gme = get_engine(geom_ptr);
    if(!save_gme)
      save_gme = gme;
    if (!gme) continue;

    DLIList<BodySM*> body_sm_list;
    geom_ptr->bodysms(body_sm_list);
    push_vg_attributes_before_modify(body_sm_list);

    BodySM *new_body_sm = NULL;
    if (!gme->regularize_entity( geom_ptr, new_body_sm ))
      stat = CUBIT_FAILURE;

    if (new_body_sm)
      new_sm_list.append(new_body_sm);
  }

  // This is bad in that it assumes all of the gmes will be the
  // same but I don't know that we truly support the other case
  // anyway.
  if(new_sm_list.size())
  {
    restore_vg_after_modify(new_sm_list, body_list, save_gme);
    remove_pushed_attributes(new_sm_list, body_list);
  }

  if (!finish_sm_op(body_list, new_sm_list, new_body_list))
    stat = CUBIT_FAILURE;

  do_attribute_cleanup();

  new_body_ptr = new_body_list.size() ? new_body_list.get() : 0;
  return stat;
}

CubitStatus GeometryModifyTool::test_regularize_refentity(RefEntity *old_entity_ptr)
{
   DLIList<RefEntity*> tmp_ref_ent_list(1);
   tmp_ref_ent_list.append( old_entity_ptr );
   if( GeometryQueryTool::instance()->contains_intermediate_geometry(tmp_ref_ent_list) )
   {
     return CUBIT_FAILURE;
   }

   BasicTopologyEntity* bte_ptr = dynamic_cast<BasicTopologyEntity*>(old_entity_ptr);
   if (!bte_ptr)
   {
     return CUBIT_FAILURE;
   }

   DLIList<TopologyBridge*> bridge_list;
   bte_ptr->bridge_manager()->get_bridge_list(bridge_list);
   bridge_list.reset();

   CubitStatus stat = CUBIT_SUCCESS;

   for (int i = bridge_list.size(); i--; )
   {
     TopologyBridge* bridge = bridge_list.get_and_step();
     GeometryEntity* geom_ptr = dynamic_cast<GeometryEntity*>(bridge);
     GeometryModifyEngine* gme = get_engine(geom_ptr);
     if (!gme) continue;

     if (!gme->test_regularize_entity( geom_ptr ))
       stat = CUBIT_FAILURE;
   }

   return stat;
}

//-------------------------------------------------------------------------
// Purpose       : split a multiple volume body into several bodies
//                 each containing a single volume.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 09/26/97
//-------------------------------------------------------------------------

CubitStatus GeometryModifyTool::regularize_body(Body *body_ptr,
                                                Body *&new_body )
{
   DLIList<Body*> b_list;
   b_list.append(body_ptr);
   if (!okay_to_modify( b_list, "REGULARIZE" ))
     return CUBIT_FAILURE;

   BodySM* body_sm = body_ptr->get_body_sm_ptr(), *new_sm = 0;
   GeometryModifyEngine* gme = get_engine(body_sm);
   if (!gme) {
     PRINT_ERROR("Volume does not have a modify engine.\n" );
     return CUBIT_FAILURE;
  }

   do_attribute_setup();
   DLIList<BodySM*> body_sm_list, new_bodysm_list;
   body_sm_list.append(body_sm);
   push_vg_attributes_before_modify(body_sm_list);

   CubitStatus stat = gme->regularize_body( body_sm, new_sm );
   if ( new_sm == NULL )
   {
      PRINT_ERROR("REGULARIZATION failure.\n");
      return CUBIT_FAILURE;
   }

   // remove mesh from modified body
   body_premodify(body_ptr);

   new_bodysm_list.append(new_sm);
   restore_vg_after_modify(new_bodysm_list, b_list, gme);
   remove_pushed_attributes(new_bodysm_list, b_list);

   body_sm = body_ptr->get_body_sm_ptr();
   update_body(body_ptr);

   new_body = GeometryQueryTool::instance()->make_Body(new_sm);
   GeometryQueryTool::instance()->cleanout_deactivated_geometry();

   do_attribute_cleanup();

   return stat;
}

CubitStatus
GeometryModifyTool::create_solid_bodies_from_surfs( DLIList<RefFace*> &ref_face_list,
                                            DLIList<Body*> &new_bodies,
                                            CubitBoolean keep_old,
                                            CubitBoolean heal,
                                            CubitBoolean sheet ) const
{
    //First check to make sure the data is all here.
  for ( int ii = ref_face_list.size(); ii > 0; ii-- )
  {
    RefFace *ref_face = ref_face_list.get_and_step();
    DLIList<Body*> bodies;
    ref_face->bodies(bodies);
    if ( bodies.size() > 1 )
    {
      PRINT_ERROR("Can't create a volume with %s (surface %d) is attached to more\n"
                  "than one volume, or if the attached body is not a sheet body.\n",
                  ref_face->entity_name().c_str(),
                  ref_face->id());
      return CUBIT_FAILURE;
    }
    else if ( bodies.size() == 1 )
    {
      if (!bodies.get()->is_sheet_body())
      {
        PRINT_ERROR("Can't create a volume with %s (surface %d), it is\n"
                    "attached to a body that is not a sheet body.\n",
                    ref_face->entity_name().c_str(),
                    ref_face->id());
        return CUBIT_FAILURE;
      }
    }
  }

  DLIList<TopologyEntity*> entity_list(ref_face_list.size());
  DLIList<TopologyBridge*> bridge_list(ref_face_list.size());
  DLIList<Surface*>       surface_list(ref_face_list.size());
  GeometryModifyEngine* gme;

  CAST_LIST_TO_PARENT(ref_face_list, entity_list);
  gme = common_modify_engine(entity_list, bridge_list);
  CAST_LIST(bridge_list, surface_list, Surface);

  if (!gme) {
    PRINT_ERROR("Surfaces do not share a common modify engine.\n");
    return CUBIT_FAILURE;
  }

  DLIList<ModelEntity*> query_output, query_input(ref_face_list.size());
  DLIList<Body*> body_list;

  CAST_LIST_TO_PARENT(ref_face_list, query_input);
  ModelQueryEngine::instance()->
    query_model( query_input, DagType::body_type(), query_output );
  CAST_LIST( query_output, body_list, Body );
  
  int i;
  DLIList<RefFace*> free_face_list;
  for ( i=ref_face_list.size(); i--; )
  {
    RefFace* ref_face = ref_face_list.get_and_step();
    query_output.clean_out();
    ModelQueryEngine::instance()->
      query_model( *ref_face, DagType::body_type(), query_output );
    if (!query_output.size())
      free_face_list.append(ref_face);
  }

  if (CubitUndo::get_undo_enabled())
    CubitUndo::save_state_with_cubit_file(ref_face_list);
  
  //get all the bodysm's 
  DLIList<BodySM*> old_body_sm_list;
  for( i=body_list.size(); i--; )
  {
    Body *tmp_body = body_list.get_and_step();
    TopologyBridge *tb = tmp_body->bridge_manager()->topology_bridge();
    BodySM *tmp_body_sm = CAST_TO(tb, BodySM);
    if( tmp_body_sm )
      old_body_sm_list.append( tmp_body_sm );    
  }  

  // TODO: do I need a clear and a flush here? --KGM
  GeometryModifyTool::instance()->do_attribute_setup();
  GeometryModifyTool::instance()->push_vg_attributes_before_modify(old_body_sm_list);

  DLIList<BodySM*> new_bodies_sm;
  CubitStatus stat = gme->create_solid_bodies_from_surfs( surface_list, new_bodies_sm, keep_old, heal, sheet );
  DLIList<BodySM*> body_sm_list;
  for ( i=new_bodies_sm.size(); i--; )
    body_sm_list.append( new_bodies_sm.get_and_step() );

  GeometryModifyTool::instance()->restore_vg_after_modify(new_bodies_sm, body_list, gme);
  GeometryModifyTool::instance()->remove_pushed_attributes(new_bodies_sm, body_list);

  if (!finish_sm_op( body_list, body_sm_list, new_bodies))
  {
    stat = CUBIT_FAILURE;
    if (CubitUndo::get_undo_enabled())
      CubitUndo::remove_last_undo();
  }
  else
  {
    if (CubitUndo::get_undo_enabled())
      CubitUndo::note_result_bodies(new_bodies);
  }

  GeometryModifyTool::instance()->do_attribute_cleanup();

  DLIList<int> id_list (free_face_list.size());
  while (free_face_list.size())
  {
    RefFace* ref_face = free_face_list.pop();
    if (!ref_face->get_surface_ptr())
    {
      id_list.append(ref_face->id());
      GeometryQueryTool::instance()->
        destroy_dead_entity( ref_face );
    }
  }

  GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  if (id_list.size())
    CubitUtil::list_entity_ids( "Destroyed surface(s) ", id_list );

  //new_body = new_body_list.size() ? new_body_list.get() : 0;
  return stat;
}

//-------------------------------------------------------------------------
// Purpose       : Check for virtual geometry
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/10/03
//-------------------------------------------------------------------------
bool GeometryModifyTool::contains_intermediate_geom( DLIList<Body*>& list ) const
{
  DLIList<TopologyBridge*> bridges;
  list.reset();
  for (int i = list.size(); i--; )
    bridges.append(list.get_and_step()->get_body_sm_ptr());
  return contains_intermediate_geom(bridges);
}

//-------------------------------------------------------------------------
// Purpose       : Check for virtual geometry
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/10/03
//-------------------------------------------------------------------------
bool GeometryModifyTool::contains_intermediate_geom(
                           DLIList<TopologyBridge*>& bridge_list ) const
{
  DLIList<TopologyBridge*> child_list;
  while (bridge_list.size())
  {
    TopologyBridge* bridge = bridge_list.pop();
    TBOwner* bridge_owner = bridge->owner();
    if (bridge_owner && !dynamic_cast<BridgeManager*>(bridge_owner))
      return true;

    child_list.clean_out();
    bridge->get_children_virt(child_list);
    bridge_list += child_list;
  }
  return false;
}

bool GeometryModifyTool::contains_partitions( DLIList<Body*>& list ) const
{
  DLIList<TopologyBridge*> bridges;
  list.reset();
  for (int i = list.size(); i--; )
    bridges.append(list.get_and_step()->get_body_sm_ptr());
  return contains_partitions(bridges);
}

bool GeometryModifyTool::contains_composites( DLIList<Body*>& list ) const
{
  DLIList<TopologyBridge*> bridges;
  list.reset();
  for (int i = list.size(); i--; )
    bridges.append(list.get_and_step()->get_body_sm_ptr());
  return contains_composites(bridges);
}

bool GeometryModifyTool::contains_composites(
                           DLIList<TopologyBridge*>& bridge_list ) const
{
  DLIList<TopologyBridge*> child_list;
  while (bridge_list.size())
  {
    TopologyBridge* bridge = bridge_list.pop();
    TBOwner* bridge_owner = bridge->owner();
    if (bridge_owner && !dynamic_cast<BridgeManager*>(bridge_owner) &&
      GeometryQueryTool::instance()->ige_is_composite(bridge_owner))
      return true;

    child_list.clean_out();
    bridge->get_children_virt(child_list);
    bridge_list += child_list;
  }
  return false;
}

bool GeometryModifyTool::contains_partitions(
                           DLIList<TopologyBridge*>& bridge_list ) const
{
  DLIList<TopologyBridge*> child_list;
  while (bridge_list.size())
  {
    TopologyBridge* bridge = bridge_list.pop();
    TBOwner* bridge_owner = bridge->owner();
    if (bridge_owner && !dynamic_cast<BridgeManager*>(bridge_owner) &&
      GeometryQueryTool::instance()->ige_is_partition(bridge_owner))
      return true;

    child_list.clean_out();
    bridge->get_children_virt(child_list);
    bridge_list += child_list;
  }
  return false;
}


CubitBoolean GeometryModifyTool::same_modify_engine(DLIList<TopologyEntity*> &topo_list) const
{
  GeometryModifyEngine *gePtr1 = get_engine(topo_list.get_and_step());
  GeometryModifyEngine *gePtr2;
  for( int i = 1; i < topo_list.size(); i++)
  {
    gePtr2 = get_engine(topo_list.get_and_step());
    if(gePtr1 != gePtr2)
    {
      return CUBIT_FALSE;
    }
  }
  return CUBIT_TRUE;
}

CubitBoolean GeometryModifyTool::same_modify_engine(DLIList<RefEntity*> &ref_entity_list,
						    CubitBoolean check_children) const
{

  DLIList<RefEntity*> complete_entity_list;

  //Check the check_children option and check all the children if necessary
  if (check_children)
  {
    //Make a complete list of all the RefEntitys and their children
    DLIList<RefEntity*> temp = ref_entity_list;
    RefEntity* ref_entity_ptr;
    int i;
    for (i = 0; i < ref_entity_list.size(); i++)
    {
      ref_entity_ptr = ref_entity_list.get_and_step();
      complete_entity_list.clean_out();
      ref_entity_ptr->get_all_child_ref_entities(complete_entity_list);
      temp += complete_entity_list;
    }
    complete_entity_list.clean_out();
    complete_entity_list.merge_unique(temp);
  }
  else
    complete_entity_list = ref_entity_list;

  //Now make sure all the RefEntitys are from the same geometry engine
  DLIList<TopologyEntity*> te_list;
  CAST_LIST(complete_entity_list, te_list, TopologyEntity);
  return same_modify_engine(te_list);

}

CubitStatus
GeometryModifyTool::trim_curve( RefEdge* trim_curve,
                                const CubitVector& trim_vector,
                                const CubitVector& keep_vector )
{
  // Use geometry engine to find intersections
  TopologyBridge* bridge = 0;
  GeometryModifyEngine* gme_ptr = get_engine(trim_curve, &bridge);
  Curve *new_curve = 0, *curve = dynamic_cast<Curve*>(bridge);

  if (!gme_ptr || !curve)
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<RefEdge*> tmp_ents(1);
    tmp_ents.append( trim_curve );
    CubitUndo::save_state_with_cubit_file( tmp_ents );
  }

  new_curve = gme_ptr->trim_curve( curve, trim_vector, keep_vector );
  if (!new_curve)
  {
    if( CubitUndo::get_undo_enabled() )
      CubitUndo::remove_last_undo();
    return CUBIT_FAILURE;
  }

  GeometryQueryTool::instance()->destroy_dead_entity( trim_curve );

  RefEdge* new_edge = GeometryQueryTool::instance()->make_free_RefEdge(new_curve);

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::note_result_entity( new_edge );

  PRINT_INFO("Created curve %d\n", new_edge->id());

  return CUBIT_SUCCESS;
}

CubitStatus
GeometryModifyTool::surface_intersection( RefFace *ref_face1,
                                          RefFace *ref_face2,
                                          DLIList<RefEdge*> &ref_edge_list )
{
  DLIList<TopologyEntity*> entity_list(2);
  DLIList<TopologyBridge*> bridge_list(2);
  entity_list.append(ref_face1);
  entity_list.append(ref_face2);
  GeometryModifyEngine* gme = common_modify_engine(entity_list, bridge_list);
  if (!gme)
  {
    PRINT_ERROR("Surfaces %d and %d do not share a common solid modeller.\n",
      ref_face1->id(), ref_face2->id() );
    return CUBIT_FAILURE;
  }

  bridge_list.reset();
  Surface* surf0 = dynamic_cast<Surface*>(bridge_list.get_and_step());
  Surface* surf1 = dynamic_cast<Surface*>(bridge_list.get());

  GeometryQueryEngine* gqe = surf0->get_geometry_query_engine();
  // Note the user can set the following value through
  //  set geometry accuracy <val>
  double resabs = gqe->get_sme_resabs_tolerance();

  DLIList<Curve*> curve_list;
  CubitStatus status = gme->surface_intersection( surf0, surf1, curve_list, resabs );

  if( status == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled() && curve_list.size() )
    CubitUndo::save_state();

  int i;
  Curve *curve_ptr;
  RefEdge *ref_edge_ptr;
  curve_list.reset();
  for( i=curve_list.size(); i--; )
  {
    curve_ptr = curve_list.get_and_step();
    ref_edge_ptr = GeometryQueryTool::instance()->
      make_free_RefEdge( curve_ptr );
    if( ref_edge_ptr )
      ref_edge_list.append( ref_edge_ptr );
    else
      delete curve_ptr;
  }

  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<RefEntity*> tmp_list;
    for( i=ref_edge_list.size(); i--; )
      tmp_list.append( ref_edge_list.get_and_step() );
    CubitUndo::note_result_entities( tmp_list );
  }

  return CUBIT_SUCCESS;
}

RefEdge*
GeometryModifyTool::create_arc_three( RefVertex* ref_vertex1,
                                      RefVertex* ref_vertex2,
                                      RefVertex* ref_vertex3,
                                      CubitBoolean full )
{
  DLIList<TopologyEntity*> entity_list(3);
  DLIList<TopologyBridge*> bridge_list(3);
  entity_list.append(ref_vertex1);
  entity_list.append(ref_vertex2);
  entity_list.append(ref_vertex3);
  GeometryModifyEngine* gme = common_modify_engine(entity_list, bridge_list);
  if (!gme)
  {
    PRINT_ERROR("Vertices do not share a common solid modeller.\n");
    return 0;
  }

  //if we can reuse vertices, we decide here
  bool need_new_start_point = false;
  bool need_new_end_point = false;
  if( full )
  {
    need_new_start_point = ref_vertex1->get_parents() > 0;
    if (need_new_start_point)
    {
      bridge_list.reset();
      Point *start_point = gme->make_Point( ref_vertex1->coordinates() );
      bridge_list.change_to( start_point );
    }
  }
  else
  {
    need_new_start_point = ref_vertex1->get_parents() > 0;
    need_new_end_point = ref_vertex3->get_parents() > 0;

    if (need_new_start_point)
    {
      bridge_list.reset();
      Point *start_point = gme->make_Point( ref_vertex1->coordinates() );
      bridge_list.change_to( start_point );
    }
    if (need_new_end_point)
    {
      bridge_list.reset();
      bridge_list.last();
      Point *end_point = gme->make_Point( ref_vertex3->coordinates());
      bridge_list.change_to( end_point );
    }
  }

  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<RefVertex*> vertices_to_save;
    if( !need_new_start_point )
      vertices_to_save.append( ref_vertex1 );
    if( !need_new_end_point )
      vertices_to_save.append( ref_vertex3 );

    if( vertices_to_save.size() )
      CubitUndo::save_state_with_cubit_file( vertices_to_save, true );
    else
      CubitUndo::save_state();
  }

  bridge_list.reset();
  Point* point0 = dynamic_cast<Point*>(bridge_list.next(0));
  Point* point1 = dynamic_cast<Point*>(bridge_list.next(1));
  Point* point2 = dynamic_cast<Point*>(bridge_list.next(2));
  Curve* curve = gme->create_arc_three( point0, point1, point2, full );
  if (!curve)
  {
    if( CubitUndo::get_undo_enabled() )
      CubitUndo::remove_last_undo();
    return 0;
  }

  RefEdge* result = GeometryQueryTool::instance()->make_free_RefEdge(curve);
  PRINT_INFO("Created curve %d\n", result->id());

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::note_result_entity( result );

  return result;
}

RefEdge*
GeometryModifyTool::create_arc_three( RefEdge* ref_edge1,
                                      RefEdge* ref_edge2,
                                      RefEdge* ref_edge3,
                                      CubitBoolean full )
{
  DLIList<TopologyEntity*> entity_list(3);
  DLIList<TopologyBridge*> bridge_list(3);
  entity_list.append(ref_edge1);
  entity_list.append(ref_edge2);
  entity_list.append(ref_edge3);
  GeometryModifyEngine* gme = common_modify_engine(entity_list, bridge_list);
  if (!gme)
  {
    PRINT_ERROR("The input curves must be associated with a common solid modeller.\n");
    return NULL;
  }

  bridge_list.reset();
  Curve* curve0 = dynamic_cast<Curve*>(bridge_list.next(0));
  Curve* curve1 = dynamic_cast<Curve*>(bridge_list.next(1));
  Curve* curve2 = dynamic_cast<Curve*>(bridge_list.next(2));
  Curve* curve = gme->create_arc_three( curve0, curve1, curve2, full );
  if (!curve)
    return 0;

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state();

  RefEdge* result = GeometryQueryTool::instance()->make_free_RefEdge(curve);

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::note_result_entity( result );

  PRINT_INFO("Created curve %d\n", result->id());
  return result;
}

RefEdge*
GeometryModifyTool::create_arc_center_edge( RefVertex* ref_vertex1,
                                            RefVertex* ref_vertex2,
                                            RefVertex* ref_vertex3,
                                            const CubitVector &normal,
                                            double radius,
                                            CubitBoolean full )
{
  DLIList<TopologyEntity*> entity_list(3);
  DLIList<TopologyBridge*> bridge_list(3);
  entity_list.append(ref_vertex1);
  entity_list.append(ref_vertex2);
  entity_list.append(ref_vertex3);
  GeometryModifyEngine* gme = common_modify_engine(entity_list, bridge_list);
  if (!gme)
  {
    PRINT_ERROR("Vertices do not share a common modify engine.\n");
    return 0;
  }

  //if we can reuse vertices, we decide here
  bool need_new_start_point, need_new_end_point = false;
  if( full )
  {
    need_new_start_point = ref_vertex2->get_parents() > 0;
    if (need_new_start_point)
    {
      bridge_list.reset();
      bridge_list.step(1);
      Point *start_point = gme->make_Point( ref_vertex2->coordinates() );
      bridge_list.change_to( start_point );
    }
  }
  else
  {
    need_new_start_point = ref_vertex2->get_parents() > 0;
    need_new_end_point = ref_vertex3->get_parents() > 0;

    if (need_new_start_point)
    {
      bridge_list.reset();
      bridge_list.step();
      Point *start_point = gme->make_Point( ref_vertex2->coordinates() );
      bridge_list.change_to( start_point );
    }
    if (need_new_end_point)
    {
      bridge_list.last();
      Point *end_point = gme->make_Point( ref_vertex3->coordinates());
      bridge_list.change_to( end_point );
    }
  }

  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<RefVertex*> vertices_to_save;
    if( !need_new_start_point )
      vertices_to_save.append( ref_vertex2 );
    if( !need_new_end_point )
      vertices_to_save.append( ref_vertex3 );

    if( vertices_to_save.size() )
      CubitUndo::save_state_with_cubit_file( vertices_to_save, true );
    else
      CubitUndo::save_state();
  }

  bridge_list.reset();
  Point* point0 = dynamic_cast<Point*>(bridge_list.next(0));
  Point* point1 = dynamic_cast<Point*>(bridge_list.next(1));
  Point* point2 = dynamic_cast<Point*>(bridge_list.next(2));
  Curve* curve = gme->create_arc_center_edge( point0, point1, point2,
                                              normal, radius, full );
  if (!curve)
  {
    if( CubitUndo::get_undo_enabled() )
      CubitUndo::remove_last_undo();
    return 0;
  }

  RefEdge* result = GeometryQueryTool::instance()->make_free_RefEdge(curve);
  PRINT_INFO("Created curve %d\n", result->id());

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::note_result_entity( result );

  return result;
}

CubitStatus
GeometryModifyTool::create_curve_combine( DLIList<RefEdge*>& ref_edge_list,
                                    RefEdge *&new_ref_edge_ptr )
{
  DLIList<TopologyBridge*> bridge_list(ref_edge_list.size());
  CubitStatus result = CUBIT_FAILURE;

  int count = ref_edge_list.size();
  if (count == 0)
  {
    PRINT_ERROR("No edges to combine.\n");
    return result;
  }

  DLIList<TopologyEntity*> entity_list(count);
  CAST_LIST_TO_PARENT( ref_edge_list, entity_list );
  GeometryModifyEngine* gme = common_modify_engine(entity_list, bridge_list);
  if (!gme)
  {
    PRINT_ERROR("Edges do not share a common modify engine.\n");
    return result;
  }

  Curve* new_curve_ptr;
  DLIList<Curve*> curve_list(count);
  CAST_LIST( bridge_list, curve_list, Curve );
  result = gme->create_curve_combine(curve_list, new_curve_ptr);
  if (new_curve_ptr)
      new_ref_edge_ptr = CAST_TO(new_curve_ptr->topology_entity(), RefEdge);
  

  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Find a common geometry engine
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/28/00
//-------------------------------------------------------------------------
GeometryModifyEngine*
GeometryModifyTool::common_modify_engine( DLIList<TopologyEntity*>& topology_list,
                                          DLIList<TopologyBridge*>& engine_bridges,
                                          CubitBoolean allow_composites ) const
{
  GeometryModifyEngine* gme_ptr = 0;
  if(allow_composites)
  {
    int i;
    engine_bridges.clean_out();
    for(i=topology_list.size(); i--;)
    {
      TopologyEntity* topo_ptr = topology_list.get_and_step();
      if (!topo_ptr)
        return (GeometryModifyEngine*)NULL;
      TopologyBridge *first_bridge = topo_ptr->bridge_manager()->topology_bridge();
      GeometryQueryEngine *gqe = first_bridge->get_geometry_query_engine();
      DLIList<TopologyBridge*> underlying_bridge_list;
      gqe->get_underlying_bridges(first_bridge, underlying_bridge_list);
      if(underlying_bridge_list.size() == 0)
        underlying_bridge_list.append(first_bridge);
      int k;
      for(k=underlying_bridge_list.size(); k--;)
      {
        TopologyBridge *bridge_ptr = underlying_bridge_list.get_and_step();
        engine_bridges.append( bridge_ptr );
        GeometryModifyEngine *cur_gme = get_engine(bridge_ptr);
        if(!gme_ptr)
          gme_ptr = cur_gme;
        else
        {
          if(gme_ptr != cur_gme)
          {
            gme_ptr = NULL;
            k=0;
            i=0;
          }
        }
      }
    }
  }
  else
  {
    topology_list.reset();

    TopologyEntity* topo_ptr = topology_list.get_and_step();
    if (!topo_ptr)
      return (GeometryModifyEngine*)NULL;
    DLIList<TopologyBridge*> first_bridge_list;
    topo_ptr->bridge_manager()->get_bridge_list( first_bridge_list );

    first_bridge_list.reset();
    for( int i = first_bridge_list.size(); i > 0; i-- )
    {
      TopologyBridge* bridge_ptr = first_bridge_list.get_and_step();
      engine_bridges.clean_out();
      engine_bridges.append( bridge_ptr );
      gme_ptr = get_engine(bridge_ptr);

      if( !gme_ptr )
      return (GeometryModifyEngine*)NULL;

      topology_list.reset();
      topology_list.step(); //skip first entry
      for( int j = topology_list.size(); j > 1; j-- )
      {
        topo_ptr = topology_list.get_and_step();
        bridge_ptr = topo_ptr->bridge_manager()->topology_bridge(gme_ptr->get_gqe());
        if( bridge_ptr ) engine_bridges.append( bridge_ptr );
        else break;
      }

      if( engine_bridges.size() == topology_list.size() )
        break;

      gme_ptr = 0;
    }
  }
  if( !gme_ptr )
  {
    engine_bridges.clean_out();
    PRINT_ERROR("Entities do not belong to the same geometry engine.\n");
  }
  return gme_ptr;
}

//-------------------------------------------------------------------------
// Purpose       : Find common modify engine for Bodies
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/26/03
//-------------------------------------------------------------------------
GeometryModifyEngine*
GeometryModifyTool::common_modify_engine( DLIList<Body*>& input,
                                          DLIList<BodySM*>& output) const
{
  input.reset();
  Body* body_ptr = input.get();
  BodySM* body_sm_ptr = body_ptr->get_body_sm_ptr();
  GeometryModifyEngine* engine = get_engine(body_sm_ptr);

  for (int i = input.size(); i--; )
  {
    body_ptr = input.get_and_step();
    body_sm_ptr = body_ptr->get_body_sm_ptr();
    output.append(body_sm_ptr);

    if (!body_sm_ptr)
    {
      PRINT_ERROR("Body %d is invalid -- no attached BodySM.\n", body_ptr->id());
      output.clean_out();
      return 0;
    }

    if (get_engine(body_sm_ptr) != engine)
    {
      output.clean_out();
      return 0;
    }
  }

  return engine;
}

//-------------------------------------------------------------------------
// Purpose       : Get modify engine common to input RefFaces and RefEdges.
//
// Special Notes : Wrapper around the TopologyEntity/TopologyBridge form.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/26/03
//-------------------------------------------------------------------------
GeometryModifyEngine*
GeometryModifyTool::common_modify_engine( DLIList<RefFace*>& face_list,
                                          DLIList<RefEdge*>& edge_list,
                                          DLIList<Surface*>& surf_list,
                                          DLIList<Curve*>& curv_list,
                                          CubitBoolean allow_composites) const
{
  int i;
  const int count = face_list.size() + edge_list.size();
  DLIList<TopologyEntity*> entity_list(count);
  DLIList<TopologyBridge*> bridge_list(count);

  face_list.reset();
  edge_list.reset();
  for (i = face_list.size(); i--; )
    entity_list.append(face_list.get_and_step());
  for (i = edge_list.size(); i--; )
    entity_list.append(edge_list.get_and_step());

  GeometryModifyEngine* engine = common_modify_engine(entity_list, bridge_list, allow_composites);
  if (!engine)
    return 0;

  entity_list.reset();
  CAST_LIST(bridge_list, surf_list, Surface);
  CAST_LIST(bridge_list, curv_list, Curve  );
  if(allow_composites)
  {
    if(surf_list.size() >= face_list.size() && curv_list.size() >= edge_list.size())
      return engine;
    else
      return 0;
  }
  else
  {
    if (surf_list.size() != face_list.size() || curv_list.size() != edge_list.size())
      return 0;
    else
      return engine;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Get modify engine common to input RefFaces.
//
// Special Notes : Wrapper around the TopologyEntity/TopologyBridge form.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/05/03
//-------------------------------------------------------------------------
GeometryModifyEngine*
GeometryModifyTool::common_modify_engine( DLIList<RefFace*>& face_list,
                                          DLIList<Surface*>& surf_list,
                                          CubitBoolean allow_composites) const
{
  const int size = face_list.size();
  DLIList<TopologyEntity*> topo_list(size);
  DLIList<TopologyBridge*> geom_list(size);
  GeometryModifyEngine* result;

  CAST_LIST_TO_PARENT( face_list, topo_list );
  result = common_modify_engine( topo_list, geom_list, allow_composites );

  CAST_LIST( geom_list, surf_list, Surface );
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Get modify engine common to input RefEdges.
//
// Special Notes : Wrapper around the TopologyEntity/TopologyBridge form.
//
// Creator       : Steve Storm
//
// Creation Date : 03/09/05
//-------------------------------------------------------------------------
GeometryModifyEngine* GeometryModifyTool::common_modify_engine(
                                      DLIList<RefEdge*>& edge_list,
                                      DLIList<Curve*>& curve_list,
                                      CubitBoolean allow_composites) const
{
  const int size = edge_list.size();
  DLIList<TopologyEntity*> topo_list(size);
  DLIList<TopologyBridge*> geom_list(size);
  GeometryModifyEngine* result;

  CAST_LIST_TO_PARENT( edge_list, topo_list );
  result = common_modify_engine( topo_list, geom_list, allow_composites );

  CAST_LIST( geom_list, curve_list, Curve );
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Get modify engine common to input RefVertices.
//
// Special Notes : Wrapper around the TopologyEntity/TopologyBridge form.
//
// Creator       : Steve Storm
//
// Creation Date : 03/25/05
//-------------------------------------------------------------------------
GeometryModifyEngine*
GeometryModifyTool::common_modify_engine( DLIList<RefVertex*>& vertex_list,
                                          DLIList<Point*>& point_list,
                                          CubitBoolean allow_composites) const
{
  const int size = vertex_list.size();
  DLIList<TopologyEntity*> topo_list(size);
  DLIList<TopologyBridge*> geom_list(size);
  GeometryModifyEngine* result;

  CAST_LIST_TO_PARENT( vertex_list, topo_list );
  result = common_modify_engine( topo_list, geom_list, allow_composites );

  CAST_LIST( geom_list, point_list, Point );
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Pull RefFaces with a common GeometryModifyEngine out of
//                 the input ref_face_list.  Place their surfaces in the
//                 output surf_list, and return the common modify engine.
//
// Special Notes : the function returns a NULL pointer if a RefFace without
//                 a modify engine is found in the input list.
//
// Creator       : Steve Storm
//
// Creation Date : 03/02/08
//-------------------------------------------------------------------------
GeometryModifyEngine*
GeometryModifyTool::pull_common_surfs( DLIList<RefFace*> &ref_face_list,
                                       DLIList<RefFace*> &common_face_list,
                                       DLIList<Surface*> &common_surf_list )
{
  int i;
  RefFace *ref_face_ptr;
  Surface *surf_ptr;

  GeometryModifyEngine *gme1 = 0, *gme2 = 0;

  ref_face_list.reset();
  for( i=0; i<ref_face_list.size(); i++ )
  {
    ref_face_ptr = ref_face_list.get();
    surf_ptr = ref_face_ptr->get_surface_ptr();

    if( i==0 )
    {
      common_face_list.append( ref_face_ptr );
      common_surf_list.append( surf_ptr );
      gme1 = get_engine( surf_ptr );
      if (!gme1)
      {
        PRINT_ERROR("Surface %d does not have a modify engine.\n", ref_face_ptr->id());
        return 0;
      }
      ref_face_list.change_to( NULL );
      ref_face_list.step();
      continue;
    }

    gme2 = get_engine( surf_ptr );
    if (!gme2)
    {
      PRINT_ERROR("Surface %d does not have a modify engine.\n", ref_face_ptr->id());
      return 0;
    }

    if( gme2 == gme1 )
    {
      common_face_list.append( ref_face_ptr );
      common_surf_list.append( surf_ptr );
      ref_face_list.change_to( NULL );
    }

    ref_face_list.step();
  }

  ref_face_list.remove_all_with_value( NULL );

  return gme1;
}

// Separates the list of bodies so that there is one body per volume
// after a webcut.  Checks the sepAfterWebcut flag.
CubitStatus GeometryModifyTool::separate_body_after_webcut( DLIList<Body*> &input_list,
                                              DLIList<Body*> &output_list) const
{
    //First see if we should spearate.
  if ( !sepAfterWebcut )
  {
    output_list = input_list;
    return CUBIT_SUCCESS;
  }
  DLIList<Body*> temp_body_list;
  for ( int ii = input_list.size(); ii > 0; ii-- )
  {
    Body *body_ptr = input_list.get_and_step();
    temp_body_list.clean_out();

    bool undo_setting = CubitUndo::get_undo_enabled();
    if( undo_setting == true )
      CubitUndo::set_undo_enabled( false );

    split_body(body_ptr, temp_body_list);

    if( undo_setting == true )
      CubitUndo::set_undo_enabled( true );

    output_list += temp_body_list;
  }
  return CUBIT_SUCCESS;
}

GeometryModifyEngine *GeometryModifyTool::get_engine(TopologyBridge *tb_ptr) const
{
  int i;
  GeometryModifyEngine *gme;
  for (i = 0; i < gmeList.size(); i++) {
    gme = gmeList.next(i);
    if (gme->is_modify_engine(tb_ptr)) return gme;
  }
  return NULL;
}

GeometryModifyEngine *GeometryModifyTool::get_engine(TopologyEntity *te_ptr,
                                                     TopologyBridge **bridge) const
{
  int i;
  GeometryModifyEngine *gme = 0;
  TopologyBridge* tb_ptr = NULL;
  BridgeManager* bm = te_ptr->bridge_manager();
  DLIList<TopologyBridge*> bridges(bm->number_of_bridges());
  bm->get_bridge_list(bridges);
  bridges.reset();
  for (i = bridges.size(); !gme && i--; )
  {
    tb_ptr = bridges.get_and_step();
    gme = get_engine(tb_ptr);
  }

  if (bridge && gme)
    *bridge = tb_ptr;

  return gme;
}

CubitStatus GeometryModifyTool::get_offset_intersections( RefEdge* ref_edge1, RefEdge* ref_edge2,
                                                          DLIList<CubitVector*>& intersection_list,
                                                          double offset, CubitBoolean ext_first )
{
  DLIList<TopologyEntity*> entity_list(2);
  DLIList<TopologyBridge*> bridge_list(2);
  entity_list.append(ref_edge1);
  entity_list.append(ref_edge2);
  GeometryModifyEngine* gme_ptr1 = common_modify_engine(entity_list,bridge_list);

  if( gme_ptr1 == 0 )
  {
      PRINT_ERROR( "Curves %d and %d do not have the same underlying geometry modeling engine\n"
                   "      For intersection calculations, they must be the same\n",
         ref_edge1->id(), ref_edge2->id() );
      return CUBIT_FAILURE;
   }

   bridge_list.reset();
   Curve* curve0 = dynamic_cast<Curve*>(bridge_list.next(0));
   Curve* curve1 = dynamic_cast<Curve*>(bridge_list.next(1));

   return gme_ptr1->get_offset_intersections( curve0, curve1, intersection_list,
                                              offset, ext_first );
}

CubitStatus GeometryModifyTool::get_offset_intersections( RefEdge* ref_edge_ptr,
                                                          RefFace* ref_face_ptr,
                                                          DLIList<CubitVector*>& intersection_list,
                                                          double offset,
                                                          CubitBoolean ext_surf )
{
  // If curve is straight and surface is planar, compute their intersection;
  // otherwise use the geometry engine to do it.

  Curve* curve_ptr = ref_edge_ptr->get_curve_ptr();
  Surface* surface_ptr = ref_face_ptr->get_surface_ptr();

  if( curve_ptr == NULL )
  {
      PRINT_ERROR("Unable to retrieve underlying geometric entity of Curve %d\n"
      "       This is a bug - please report it\n", ref_edge_ptr->id() );
    return CUBIT_FAILURE;
  }

  if( surface_ptr == NULL )
  {
      PRINT_ERROR("Unable to retrieve underlying geometric entity of Surface %d\n"
      "       This is a bug - please report it\n", ref_face_ptr->id() );
    return CUBIT_FAILURE;
  }

  // If straight line and plane, find location right here analytically.
  CubitVector pln_origin, pln_normal;
  CubitVector crv_origin, crv_direction;
  if( ref_face_ptr->get_point_normal( pln_origin, pln_normal ) &&
     ref_edge_ptr->get_point_direction( crv_origin, crv_direction ) )
  {
     double pln_orig[3], pln_norm[3];
     pln_orig[0]=pln_origin.x(); pln_orig[1]=pln_origin.y(); pln_orig[2]=pln_origin.z();
     pln_norm[0]=pln_normal.x(); pln_norm[1]=pln_normal.y(); pln_norm[2]=pln_normal.z();

     double crv_orig[3], crv_dir[3];
     crv_orig[0]=crv_origin.x(); crv_orig[1]=crv_origin.y(); crv_orig[2]=crv_origin.z();
     crv_dir[0]=crv_direction.x(); crv_dir[1]=crv_direction.y(); crv_dir[2]=crv_direction.z();

     AnalyticGeometryTool *agt = AnalyticGeometryTool::instance();

     if( agt->is_vec_perp( pln_norm, crv_dir ) )
     {
        PRINT_ERROR( "Line is parallel to the plane - intersection not possible\n" );
        return CUBIT_FAILURE;
     }

     double ang = agt->angle_vec_vec( pln_norm, crv_dir );
     if( ang > AGT_PI_DIV_2 )
        ang = AGT_PI - ang;

     double int_pnt[3];
     agt->int_ln_pln( crv_orig, crv_dir, pln_orig, pln_norm, int_pnt );

     // ang2 is angle between line and plane
     double ang2 = AGT_PI_DIV_2 - ang;

     double hypotenuse = offset/sin(ang2);

     double final_pnt[3];

     agt->next_pnt( int_pnt, crv_dir, hypotenuse, final_pnt );

     double end1[3], end2[3];
     CubitVector start, end;
     start = ref_edge_ptr->start_coordinates();
     end = ref_edge_ptr->end_coordinates();
     end1[0]=start.x(); end1[1]=start.y(); end1[2]=start.z();
     end2[0]=end.x(); end2[1]=end.y(); end2[2]=end.z();

     CubitVector curve_position;
     if( agt->is_pnt_on_ln_seg( final_pnt, end1, end2 ) )
     {
        curve_position.x( final_pnt[0] );
        curve_position.y( final_pnt[1] );
        curve_position.z( final_pnt[2] );
     }
     else
     {
        agt->reverse_vec( crv_dir, crv_dir );

        agt->next_pnt( int_pnt, crv_dir, hypotenuse, final_pnt );

        if( agt->is_pnt_on_ln_seg( final_pnt, end1, end2 ) )
        {
           curve_position.x( final_pnt[0] );
           curve_position.y( final_pnt[1] );
           curve_position.z( final_pnt[2] );
        }
        else
        {
           PRINT_ERROR( "Resultant point does not lie on bounded curve\n" );
           return CUBIT_FAILURE;
        }
     }
     CubitVector *vec_ptr = new CubitVector( curve_position );
     intersection_list.append( vec_ptr );
  }
  else
  {
     DLIList<TopologyEntity*> entity_list(2);
     DLIList<TopologyBridge*> bridge_list(2);
     entity_list.append(ref_edge_ptr);
     entity_list.append(ref_face_ptr);
     GeometryModifyEngine* gme_ptr1 = common_modify_engine(entity_list,bridge_list);

     if( gme_ptr1 == 0 )
     {
        PRINT_ERROR( "Curve %d and Surface %d do not have the same underlying geometry modeling engine\n"
           "       For intersection calculations, they must be the same\n",
           ref_edge_ptr->id(), ref_face_ptr->id() );
        return CUBIT_FAILURE;
     }

     bridge_list.reset();
     curve_ptr = dynamic_cast<Curve*>(bridge_list.next(0));
     surface_ptr = dynamic_cast<Surface*>(bridge_list.next(1));

     // Use geometry engine to find intersections
     return gme_ptr1->get_offset_intersections( curve_ptr, surface_ptr,
                                                intersection_list, offset, ext_surf );

  }

  return CUBIT_FAILURE;
}

CubitStatus GeometryModifyTool::get_mid_plane( RefFace *ref_face1,
                                               RefFace *ref_face2,
                                               Body *body_to_trim_to,
                                               DLIList<RefFace*> &mid_plane_surfs ) const
{
  if( ref_face1 == ref_face2 )
  {
    PRINT_ERROR("Cannot create midplane between the same surface.\n"
                "       Surface %d was entered twice\n",  ref_face1->id() );
    return CUBIT_FAILURE;
  }

  BodySM* body_sm_to_trim_to = body_to_trim_to->get_body_sm_ptr();
  GeometryModifyEngine *gme1_ptr = get_engine(body_sm_to_trim_to);
  if ( !gme1_ptr )
  {
    PRINT_ERROR("Geometry can't be modified, no associated modify engine.\n");
    return CUBIT_FAILURE;
  }

  CubitVector normal_1, normal_2, point_1, point_2, point_3;
  CubitPlane plane_1, plane_2;
  CubitVector p_mid, n_mid;

  point_1 = ref_face1->center_point();
  point_2 = ref_face2->center_point();

  normal_1 = ref_face1->normal_at(point_1);
  normal_2 = ref_face2->normal_at(point_2);

  plane_1 = CubitPlane(normal_1,point_1);
  plane_2 = CubitPlane(normal_2,point_2);

  if(point_1 == point_2)
  {
    PRINT_ERROR("Since both surfaces share the same point, the midplane is not well-defined\n");
    return CUBIT_FAILURE;
  }
  else
  {
    CubitVector temp1 = point_2;
    temp1 = plane_1.project(temp1);
    temp1 -= point_2;
    if ( temp1.length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
    {
      PRINT_ERROR("Since both planes are the same, the mid-plane is not well-defined.\n");
      return CUBIT_FAILURE;
    }
  }

  if(plane_1.normal()==plane_2.normal() || plane_1.normal()==-plane_2.normal())
  {
    p_mid = (point_1+point_2)/2;
    n_mid = plane_1.normal();
  }
  else
  {

    CubitVector direction_of_line;
    plane_1.intersect(plane_2,p_mid,direction_of_line);
    direction_of_line.normalize();

    // Find if point_1 and point_2 are on the line of intersection
    // If they are, then the mid-plane is not well-defined
    CubitVector p1 = point_1-p_mid;
    CubitVector p2 = point_2-p_mid;
    p1.normalize();
    p2.normalize();

    if(p1==direction_of_line || p1==-direction_of_line)
    {
      PRINT_ERROR("P1 is on the line of intersection.\n");
      return CUBIT_FAILURE;
    }

    if(p2==direction_of_line || p2==-direction_of_line)
    {
      PRINT_ERROR("P2 is on the line of intersection.\n");
      return CUBIT_FAILURE;
    }

    CubitVector v1 = p1 - (p1%direction_of_line)*direction_of_line;
    v1.normalize();

    CubitVector v2 = p2 - (p2%direction_of_line)*direction_of_line;
    v2.normalize();

    n_mid = v1 - v2;
    n_mid.normalize();
  }

  CubitPlane mid_plane(n_mid, p_mid);
  point_1 = p_mid;
    //find three points that will define the infinite plane from the
    //mid plane.
  CubitVector test1(1,0,0), test1n(-1,0,0),test2(0,1,0);
  CubitVector direction1;
    //through the point in any direction just not along the
    //normal direction.
  if(n_mid != test1 && n_mid != test1n )
    direction1 = test1 + n_mid;
  else
    direction1 = test2 + n_mid;

  point_2 = p_mid + direction1;
  point_2 = mid_plane.project(point_2);

  direction1 = point_2-point_1;
  CubitVector direction2 = direction1*n_mid;
  point_3 = point_1 + direction2;

  BodySM* midplane_body_sm = NULL;
  CubitStatus ret = gme1_ptr->get_mid_plane(point_1, point_2, point_3,
                                            body_sm_to_trim_to, midplane_body_sm );

  if (midplane_body_sm)
  {
#ifdef BOYD17
    DLIList<Body*> bodies;
    DLIList<Surface*> surfs;
#endif

    Body *midplane_body;

    midplane_body = GeometryQueryTool::instance()->make_Body(midplane_body_sm);

    DLIList<RefFace*> ref_faces;
    midplane_body->ref_faces( ref_faces );

    //make each surface of the body into its own body
    int i;
    for( i=0; i<ref_faces.size(); i++ )
    {
      RefEntity *new_entity_ptr;
      new_entity_ptr = GeometryModifyTool::instance()->copy_refentity(ref_faces.get_and_step());
      RefFace *ref_face_ptr = CAST_TO(new_entity_ptr, RefFace);
      mid_plane_surfs.append( ref_face_ptr );
    }
    GeometryQueryTool::instance()->delete_Body( midplane_body );
  }
  else
    return CUBIT_FAILURE;

  return ret;
}

//=============================================================================
// Function   : get_planar_mid_surface
// Member Type: LOCAL
// Description: Calculates a mid-surface between 2 planar surfaces.
// Author     : Philippe Pebay
// Date       : 05/15/06
//=============================================================================
CubitStatus get_planar_mid_surface( RefFace* ref_face1,
				    RefFace* ref_face2,
				    BodySM* body_sm_to_trim_to,
				    BodySM*& midsurface_body_sm,
				    GeometryModifyEngine *gme_ptr )
{
    CubitVector normal_1, normal_2, point_1, point_2, point_3;
    CubitPlane plane_1, plane_2;
    CubitVector p_mid, n_mid;

    point_1 = ref_face1->center_point();
    point_2 = ref_face2->center_point();

    normal_1 = ref_face1->normal_at(point_1);
    normal_2 = ref_face2->normal_at(point_2);

    plane_1 = CubitPlane(normal_1,point_1);
    plane_2 = CubitPlane(normal_2,point_2);

    if(point_1 == point_2)
    {
      PRINT_ERROR( "In GeometryModifyTool:: get_planar_mid_surface\n"
		   "       Since both surfaces share the same point, the midsurface is not well-defined\n");
      return CUBIT_FAILURE;
    }
    else
    {
      CubitVector temp1 = point_2;
      temp1 = plane_1.project(temp1);
      temp1 -= point_2;
      if ( temp1.length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
	{
	PRINT_ERROR("In GeometryModifyTool:: get_planar_mid_surface\n"
		    "       Since both planes are the same, the midsurface is not well-defined.\n");
	return CUBIT_FAILURE;
	}
    }

    if ( ( normal_1.about_equal( normal_2 ) ) || ( (-normal_1).about_equal( normal_2 ) ) )
    {
      p_mid = (point_1+point_2)/2;
      n_mid = plane_1.normal();
    }
    else
    {
      CubitVector direction_of_line;
      plane_1.intersect(plane_2,p_mid,direction_of_line);
      direction_of_line.normalize();

      // Find if point_1 and point_2 are on the line of intersection
      // If they are, then the mid-plane is not well-defined
      CubitVector p1 = point_1-p_mid;
      CubitVector p2 = point_2-p_mid;
      p1.normalize();
      p2.normalize();

      if(p1==direction_of_line || p1==-direction_of_line)
      {
	PRINT_ERROR("In GeometryModifyTool:: get_planar_mid_surface\n"
		    "       P1 is on the line of intersection.\n");
	return CUBIT_FAILURE;
      }

      if(p2==direction_of_line || p2==-direction_of_line)
      {
	PRINT_ERROR("In GeometryModifyTool:: get_planar_mid_surface\n"
		    "       P2 is on the line of intersection.\n");
	return CUBIT_FAILURE;
      }

      CubitVector v1 = p1 - (p1%direction_of_line)*direction_of_line;
      v1.normalize();

      CubitVector v2 = p2 - (p2%direction_of_line)*direction_of_line;
      v2.normalize();

      n_mid = v1 - v2;
      n_mid.normalize();
    }

    CubitPlane mid_plane(n_mid, p_mid);
    point_1 = p_mid;

    //find three points that will define the infinite plane from the
    //mid plane.through the point in any direction just not along the
    //normal direction
    CubitVector Xdir(1,0,0), Ydir(0,1,0);
    CubitVector direction1;

    if ( ( ! n_mid.about_equal( Xdir ) ) && ( ! (-n_mid).about_equal( Xdir ) ) )
      direction1 = Xdir + n_mid;
    else
      direction1 = Ydir + n_mid;

    point_2 = p_mid + direction1;
    point_2 = mid_plane.project(point_2);

    direction1 = point_2-point_1;
    CubitVector direction2 = direction1*n_mid;
    point_3 = point_1 + direction2;

    CubitStatus ret = gme_ptr->get_mid_plane(point_1, point_2, point_3,
					      body_sm_to_trim_to, midsurface_body_sm );
    return ret;
}

//=============================================================================
// Function   : get_mid_surface
// Member Type: PUBLIC
// Description: Calculates a mid-surface between 2 surfaces (calls appropriate
//              specialized methods if needed).
// Author     : Philippe Pebay
// Date       : 03/07/06
//=============================================================================
CubitStatus GeometryModifyTool::get_mid_surface( RefFace *ref_face1,
                                               RefFace *ref_face2,
                                               Body *body_to_trim_to,
                                               DLIList<RefFace*> &mid_surface_surfs ) const
{
  if( ref_face1 == ref_face2 )
  {
    PRINT_ERROR("Cannot create midplane between the same surface.\n"
                "       Surface %d was entered twice\n",  ref_face1->id() );
    return CUBIT_FAILURE;
  }

  BodySM *body_sm_to_trim_to = body_to_trim_to->get_body_sm_ptr();
  GeometryModifyEngine *gme1_ptr = get_engine(body_sm_to_trim_to);
  if ( !gme1_ptr )
  {
    PRINT_ERROR("In GeometryModifyTool::get_mid_surface\n"
		"       Geometry can't be modified, no associated modify engine.\n");
    return CUBIT_FAILURE;
  }

  bool found_case = false;
  CubitStatus ret = CUBIT_FAILURE;
  BodySM* midsurface_body_sm = NULL;

  // Plane to plane case
  if ( ( ref_face1->geometry_type() == PLANE_SURFACE_TYPE ) && ( ref_face2->geometry_type() == PLANE_SURFACE_TYPE ) )
  {
    found_case = true;
    ret = get_planar_mid_surface( ref_face1, ref_face2, body_sm_to_trim_to, midsurface_body_sm, gme1_ptr );
  }

  // Quadric to quadric cases
  if ( ( ( ref_face1->geometry_type() == SPHERE_SURFACE_TYPE ) && ( ref_face2->geometry_type() == SPHERE_SURFACE_TYPE ) )
       ||
       ( ( ref_face1->geometry_type() == CONE_SURFACE_TYPE ) && ( ref_face2->geometry_type() == CONE_SURFACE_TYPE ) )
       ||
       ( ( ref_face1->geometry_type() == TORUS_SURFACE_TYPE ) && ( ref_face2->geometry_type() == TORUS_SURFACE_TYPE ) ) )
  {
    found_case = true;
    DLIList<TopologyEntity*> entity_list(2);
    DLIList<TopologyBridge*> bridge_list(2);

    entity_list.append(ref_face1);
    entity_list.append(ref_face2);
    GeometryModifyEngine* gme2_ptr = common_modify_engine(entity_list,bridge_list);

    if( gme2_ptr == 0 )
    {
      PRINT_ERROR( "In GeometryModifyTool::get_mid_surface\n"
		   "       Surfaces %d and %d do not have the same underlying geometry modeling engine.\n"
		   "       For midsurface calculations, they must be the same\n",
		   ref_face1->id(), ref_face2->id() );
      return CUBIT_FAILURE;
    }

    if( gme1_ptr != gme2_ptr )
    {
      PRINT_ERROR( "In GeometryModifyTool::get_mid_surface\n"
		   "       Body and surfaces do not have the same underlying geometry modeling engine.\n"
		   "       For midsurface calculations, they must be the same\n");
      return CUBIT_FAILURE;
    }

    bridge_list.reset();
    Surface* surface1_ptr = dynamic_cast<Surface*>(bridge_list.next(0));
    Surface* surface2_ptr = dynamic_cast<Surface*>(bridge_list.next(1));

    // Sphere to sphere case
    if ( ( ref_face1->geometry_type() == SPHERE_SURFACE_TYPE ) && ( ref_face2->geometry_type() == SPHERE_SURFACE_TYPE ) )
    {
      ret = gme2_ptr->get_spheric_mid_surface( surface1_ptr, surface2_ptr, body_sm_to_trim_to, midsurface_body_sm );
    }

    // Cone to cone case
    if ( ( ref_face1->geometry_type() == CONE_SURFACE_TYPE ) && ( ref_face2->geometry_type() == CONE_SURFACE_TYPE ) )
    {
      ret = gme2_ptr->get_conic_mid_surface( surface1_ptr, surface2_ptr, body_sm_to_trim_to, midsurface_body_sm );
    }

    // Torus to torus case
    if ( ( ref_face1->geometry_type() == TORUS_SURFACE_TYPE ) && ( ref_face2->geometry_type() == TORUS_SURFACE_TYPE ) )
    {
      ret = gme2_ptr->get_toric_mid_surface( surface1_ptr, surface2_ptr, body_sm_to_trim_to, midsurface_body_sm );
    }
  }

  // Unsupported pair of surfaces
  if ( ! found_case )
  {
    PRINT_ERROR("In GeometryModifyTool::get_mid_surface\n"
		"       Midsurface calculation not yet supported for such a pair of surfaces.\n");
    return CUBIT_FAILURE;
  }

  if ( midsurface_body_sm )
  {
    if(CubitUndo::get_undo_enabled())
      CubitUndo::save_state();

    DLIList<Surface*> mid_surfaces;
    DLIList<Body*> new_bodies;
    midsurface_body_sm->surfaces( mid_surfaces);
    //make each surface of the body into its own body
    int i;
    for( i=0; i<mid_surfaces.size(); i++ )
    {
      Surface *tmp_surface = mid_surfaces.get_and_step();
      bool extended_from = false;
      Surface* new_surface_ptr = gme1_ptr->make_Surface( tmp_surface, extended_from );

      Body *new_Body = make_Body(new_surface_ptr);
      new_bodies.append( new_Body );
      DLIList<RefFace*> ref_faces;
      new_Body->ref_faces(ref_faces);
      RefFace *ref_face_ptr = ref_faces.get();
      mid_surface_surfs.append( ref_face_ptr );
    }
    gme1_ptr->get_gqe()->delete_solid_model_entities( midsurface_body_sm );

    if( CubitUndo::get_undo_enabled() )
    {
      if( new_bodies.size() )
        CubitUndo::note_result_bodies( new_bodies );
      else
        CubitUndo::remove_last_undo();
    }

    return ret;
  }
  else
    return CUBIT_FAILURE;
}

CubitStatus GeometryModifyTool::set_default_gme(GeometryModifyEngine* GMEPtr)
{
    if(!GMEPtr)
        return CUBIT_FAILURE;

    int i;
    for (i = 0; i < gmeList.size(); i++)
    {
        if(GMEPtr == gmeList.get())
            break;

        gmeList.step();
    }

    if(i == gmeList.size())
        return CUBIT_FAILURE;

    GeometryModifyEngine* temp_ptr = gmeList.get();
    gmeList.remove();
    gmeList.insert_first(temp_ptr);
    return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Consolidate common setup code for tweak surface operations
// Author     : Jason Kraftcheck
// Date       : 11/05/03
//=============================================================================
GeometryModifyEngine*
GeometryModifyTool::tweak_setup( DLIList<RefFace*> &input_faces,
                                 const char* name,
                                 DLIList<Body*> &output_bodies,
                                 DLIList<Surface*> &output_surfaces,
                                 CubitBoolean allow_composites)
{
  if( input_faces.size() == 0 )
  {
    PRINT_ERROR( "No surfaces specified\n" );
    return 0;
  }

  //are any of the input faces merged? ...disallow that
  int i;
  for( i=input_faces.size(); i--; )
  {
    RefFace *tmp_face = input_faces.get();
    if( tmp_face->is_merged() )
    {
      input_faces.change_to( NULL );
      PRINT_ERROR("Surface %d is a merged surface.  Cannot perform %s operation on it.\n"
                  "       Unmerge it first\n",  tmp_face->id(), name );
    }
    input_faces.step();
  }
  input_faces.remove_all_with_value( NULL );
  if( input_faces.size() == 0 )
    return 0;

  // Get parent bodies
  DLIList<ModelEntity*> query_input(input_faces.size()), query_output;
  CAST_LIST_TO_PARENT(input_faces, query_input);
  ModelQueryEngine::instance()
    ->query_model( query_input, DagType::body_type(), query_output );
  CAST_LIST(query_output, output_bodies, Body);

  if(allow_composites)
  {
    if (!okay_to_modify( output_bodies, "TWEAK" ))
      return 0;
  }
  else
  {
    if ( contains_intermediate_geom(output_bodies))
    {
      PRINT_ERROR("%s surfaces on volumes containing virtual geometry\n"
        "       is not allowed.\n"
        "       Delete virtual geometry on these volumes before operation.\n",
        name);
      return 0;
    }
  }

  // Get engine and corresponding geom entities
  GeometryModifyEngine* gme_ptr;
  gme_ptr = common_modify_engine( input_faces, output_surfaces, allow_composites );
  if (!gme_ptr)
  {
    PRINT_ERROR("%s surfaces on volumes containing surfaces from different\n"
      "       geometry engines is not allowed.\n", name);
  }

  return gme_ptr;
}

//=============================================================================
// Description: Consolidate common setup code for tweak curve operations
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
GeometryModifyEngine*
GeometryModifyTool::tweak_setup( DLIList<RefEdge*> &input_curves,
                                 const char* name,
                                 DLIList<Body*> &output_bodies,
                                 DLIList<Curve*> &output_curves )
{
  if( input_curves.size() == 0 )
  {
    PRINT_ERROR( "No curves specified\n" );
    return 0;
  }

  //are any of the input curves merged? ...disallow that
  int i;
  for( i=input_curves.size(); i--; )
  {
    RefEdge *tmp_edge = input_curves.get();
    if( tmp_edge->is_merged() )
    {
      input_curves.change_to( NULL );
      PRINT_ERROR("Curve %d is a merged curve.  Cannot perform %s operation on it.\n"
                  "       Unmerge it first\n",  tmp_edge->id(), name );
    }
    input_curves.step();
  }
  input_curves.remove_all_with_value( NULL );
  if( input_curves.size() == 0 )
    return 0;

  // Get parent bodies
  DLIList<ModelEntity*> query_input(input_curves.size()), query_output;
  CAST_LIST_TO_PARENT(input_curves, query_input);
  ModelQueryEngine::instance()
    ->query_model( query_input, DagType::body_type(), query_output );
  CAST_LIST(query_output, output_bodies, Body);

  // Check for virtual geometry
  if ( contains_intermediate_geom(output_bodies))
  {
    PRINT_ERROR("%s curves on surfaces containing virtual geometry\n"
      "       is not allowed.\n"
      "       Delete virtual geometry on these entities before operation.\n",
      name);
    return 0;
  }

  // Get engine and corresponding geom entities
  GeometryModifyEngine* gme_ptr;
  gme_ptr = common_modify_engine( input_curves, output_curves );
  if (!gme_ptr)
  {
    PRINT_ERROR("%s curves on entities containing surfaces from different\n"
      "       geometry engines is not allowed.\n", name);
  }

  return gme_ptr;
}

//=============================================================================
// Description: Consolidate common setup code for tweak vertex operations
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
GeometryModifyEngine*
GeometryModifyTool::tweak_setup( DLIList<RefVertex*> &input_vertices,
                                 const char* name,
                                 DLIList<Body*> &output_bodies,
                                 DLIList<Point*> &output_points )
{
  if( input_vertices.size() == 0 )
  {
    PRINT_ERROR( "No vertices specified\n" );
    return 0;
  }

  //are any of the input vertices merged? ...disallow that
  int i;
  for( i=input_vertices.size(); i--; )
  {
    RefVertex *tmp_vert = input_vertices.get();
    if( tmp_vert->is_merged() )
    {
      input_vertices.change_to( NULL );
      PRINT_ERROR("Vertex %d is a merged vertex.  Cannot perform %s operation on it.\n"
                  "       Unmerge it first\n",  tmp_vert->id(), name );
    }
    input_vertices.step();
  }
  input_vertices.remove_all_with_value( NULL );
  if( input_vertices.size() == 0 )
    return 0;

  // Get parent bodies
  DLIList<ModelEntity*> query_input(input_vertices.size()), query_output;
  CAST_LIST_TO_PARENT(input_vertices, query_input);
  ModelQueryEngine::instance()
    ->query_model( query_input, DagType::body_type(), query_output );
  CAST_LIST(query_output, output_bodies, Body);

  // Check for virtual geometry
  if ( contains_intermediate_geom(output_bodies))
  {
    PRINT_ERROR("%s vertices on surfaces containing virtual geometry\n"
      "       is not allowed.\n"
      "       Delete virtual geometry on these entities before operation.\n",
      name);
    return 0;
  }

  // Get engine and corresponding geom entities
  GeometryModifyEngine* gme_ptr;
  gme_ptr = common_modify_engine( input_vertices, output_points );
  if (!gme_ptr)
  {
    PRINT_ERROR("%s vertices on entities containing surfaces from different\n"
      "       geometry engines is not allowed.\n", name);
  }

  return gme_ptr;
}

//=============================================================================
// Description: The user selects a surface they would like to idealize and also selects a radius
//				size for fillets.  The user also specifies whether to consider internal and/or external fillets.
//              The program will identify fillets which meet the users criteria and tweak remove them automatically.
//              There is also a preview and exclude curve capability.
// Author     : Jonathan Bugman
// Date       : 10/23/2008
//=============================================================================
CubitStatus GeometryModifyTool::idealize_fillet_geometry(DLIList<RefEntity*> idealize_entity,
                                                         DLIList<RefEntity*> exclude_entity,
                                                         double fillet_rad,
                                                         CubitBoolean internal_flg,
                                                         CubitBoolean external_flg,
													     CubitBoolean preview)
{
    //cast the DLIList<RefEntity> to a DLIList<RefFace>
    DLIList<RefFace*> face_to_idealize;
    CAST_LIST(idealize_entity, face_to_idealize, RefFace);

    //grabbing geometry tolerance
    double geo_tol = GeometryQueryTool::instance()->get_sme_resabs_tolerance(),temp_fillet_radius;

    //grabbing all the curve loops ONLY from surfaces which are a sheet body
    int y, i, j, z;
    DLIList<RefFace*> sheet_body_idealize_face;
    
    for(i=0; i<face_to_idealize.size(); i++)
    {
        RefFace* target_face = face_to_idealize.get_and_step();
        DLIList<Shell*> shell_list;
        target_face->shells(shell_list);
        for(j=0; j<shell_list.size(); j++)
        {
            Shell* target_shell = shell_list.get_and_step();
            if(target_face->is_nonmanifold( (GroupingEntity*)target_shell ) )
            {         
                sheet_body_idealize_face.append(face_to_idealize[i]);
            }
        }
    }

    face_to_idealize.clean_out();

    //this section is going to remove all excluded curves loopsm from the 'master' loopsm list
    DLIList<Curve*> exclude_cuves;
    DLIList <Body*> old_body_list;
    if(exclude_entity.size()>0)
    {
        //cast the exclude DLIList<RefEntity> to DLIList<RefEdge>
        DLIList<RefEdge*> exclude_edge;
        CAST_LIST(exclude_entity, exclude_edge, RefEdge);

        //switching the DLIList<RefEdge> to DLIList<Curve>
        GeometryModifyEngine* gme_ptr1;
        gme_ptr1 = tweak_setup(exclude_edge,"idealize",old_body_list,exclude_cuves);
        exclude_edge.clean_out();
    }

    //switching the DLIList<RefFace> to DLIList<Surface>
    GeometryModifyEngine* gme_ptr;
    DLIList<Surface*> sheet_body_idealize_surface;
    gme_ptr = tweak_setup(sheet_body_idealize_face,"idealize",old_body_list,sheet_body_idealize_surface);
    sheet_body_idealize_face.clean_out();

    //grab all the loops from each sheet body surface
    DLIList <LoopSM*> idealize_loopSM_list;
    for(y=0;y<sheet_body_idealize_surface.size();y++)
    {
        sheet_body_idealize_surface[y]->loopsms(idealize_loopSM_list);
    }
    sheet_body_idealize_surface.clean_out();

    //search through possible fillet curves filtering only for curves of type arc
    //if it is an arc, does it have a straight line on both sides of it, if so'
    //check the radius of the arc and if it passes test add it to the list of curves to be tweaked removed
    DLIList <Curve*> master_curve_remove_list,possible_fillet_arcs,potential_fillet,internal_fillet, external_fillet, attached_curves;
    CubitVector fillet_center_point, intersection_pt,arc_mid,test_point;
    DLIList <Point*> arc_vertices;
    for(y=0;y<idealize_loopSM_list.size();y++)
    {
        idealize_loopSM_list[y]->curves(possible_fillet_arcs);
        //doing this as a performance boost, it'll keep the code out of the next couple of for loops for situations
        //where there is no fillet possible, for instance, a hole with 2 curves will never be a fillet
        if(possible_fillet_arcs.size()>3)
        {
            for(i=0;i<possible_fillet_arcs.size();i++)
            {
                if(possible_fillet_arcs[i]->geometry_type() == ARC_CURVE_TYPE &&
                    exclude_cuves.is_in_list(possible_fillet_arcs[i])==CUBIT_FALSE)
                {
                    possible_fillet_arcs[i]->points(arc_vertices);
                    //don't need to check for one point as in a hole because I have a check that there needs to be
                    //at least 3 curves in the loop

                    //this is to check that there is only one curve attached to the arc
                    for(z=0;z<arc_vertices.size();z++)
                    {
                        arc_vertices[z]->curves(attached_curves);
                        if(attached_curves.size()!=2)
                        {
                            //I dont' think this break point is going to kick me far enough out of the for loop
                            break;
                        }
                    }

                    possible_fillet_arcs[i]->mid_point(arc_mid);
                    possible_fillet_arcs[i]->get_center_radius(fillet_center_point,temp_fillet_radius);
                    test_point = arc_mid + geo_tol * (fillet_center_point-arc_mid);
                    DLIList<Surface*> test_surf;
                    idealize_loopSM_list[y]->surfaces(test_surf);

                    //this may be dangerous but I'm assuming that a loop is on only one surface
                    CubitPointContainment cpc = test_surf[0]->point_containment(test_point);
                    if(temp_fillet_radius <= fillet_rad && cpc==CUBIT_PNT_INSIDE)
                    {
                        external_fillet.append(possible_fillet_arcs[i]);
                    }
                    else if(temp_fillet_radius <= fillet_rad && cpc==CUBIT_PNT_OUTSIDE)
                    {
                        internal_fillet.append(possible_fillet_arcs[i]);
                    }
                }
            }
        }
        possible_fillet_arcs.clean_out();
    }

    if(internal_flg==CUBIT_TRUE)
    {
        master_curve_remove_list+=internal_fillet;
    }
    if(external_flg==CUBIT_TRUE)
    {
        master_curve_remove_list+=external_fillet;
    }

    //if no arcs are found to be removed, warn the user.
    if(master_curve_remove_list.size()==0)
    {
        PRINT_INFO( "Failed to find any fillets which met users requirements\n\n" );
        //I'm returning success here even though no curves were found
        return CUBIT_SUCCESS;
    }
    else if(preview == CUBIT_TRUE)
    {
        DLIList<BodySM*> new_bodysm_list;
        bool old_error_flag = GET_ERROR_FLAG();
        SET_ERROR_FLAG(false); // don't throw any tweak_remove errors

        CubitStatus stat = gme_ptr->tweak_remove(master_curve_remove_list, new_bodysm_list,CUBIT_FALSE, CUBIT_TRUE );

        SET_ERROR_FLAG(old_error_flag); // turn errors back on 
        if(stat==CUBIT_FAILURE)
        {
            PRINT_WARNING("At least one of the fillets which met your requirements \n"
                "           can't be preview due to the curve's geometry\n");
        }

        //output the number of holes or slots which were found
        PRINT_INFO("Found %d fillets which met idealization parameters\n\n", master_curve_remove_list.size());
        return CUBIT_SUCCESS;
    }
    else
    {
        DLIList<BodySM*> new_bodysm_list;
        bool old_error_flag = GET_ERROR_FLAG();
        SET_ERROR_FLAG(false); // don't throw any tweak_remove errors

        //pass master_curve_remove_list to the tweak_remove command
        CubitStatus stat = gme_ptr->tweak_remove(master_curve_remove_list, new_bodysm_list,CUBIT_FALSE, CUBIT_FALSE );
        if(stat==CUBIT_FAILURE)
        {
            PRINT_WARNING("At least one of the fillets which met your requirements \n"
                "           can't be tweaked due to the curve's geometry\n");
        }
        SET_ERROR_FLAG(old_error_flag); // turn errors back on 

        //update DAG
        DLIList<Body*> new_body_list;
        stat = finish_sm_op( old_body_list, new_bodysm_list ,new_body_list );
        //output the number of holes or slots which were found
        PRINT_INFO("Found %d fillets which met idealization parameters\n\n", master_curve_remove_list.size());
        return CUBIT_SUCCESS;
    }
}

//=============================================================================
// Description: The user selects a surface they would like to idealize and also selects a radius
//				size for holes and or selects a radius and length for slots.  The program will identify
//              'holes' and 'slots' which meet the users criteria and tweak remove them automatically.
//              There is also a preview and exclude curve capability.
// Author     : Jonathan Bugman
// Date       : 10/23/2008
//=============================================================================
CubitStatus GeometryModifyTool::idealize_hole_slot_geometry(DLIList<RefEntity*> idealize_entity,
                                                            DLIList<RefEntity*> exclude_entity,
                                                            double arc_radius,
                                                            double slot_arc_radius,
                                                            double slot_length,
													        CubitBoolean preview)
{
    //cast the DLIList<RefEntity> to a DLIList<RefFace>
    DLIList<RefFace*> face_to_idealize;
    CAST_LIST(idealize_entity, face_to_idealize, RefFace);

    //grabbing geometry tolerance
    double geo_tol = GeometryQueryTool::instance()->get_sme_resabs_tolerance();

    //grabbing all the curve loops ONLY from surfaces which are a sheet body
    int y=0, i=0, j=0;
    DLIList<RefFace*> sheet_body_idealize_face;

    for(i=0; i<face_to_idealize.size(); i++)
    {
        RefFace* target_face = face_to_idealize.get_and_step();
        DLIList<Shell*> shell_list;
        target_face->shells(shell_list);
        for(j=0; j<shell_list.size(); j++)
        {
            Shell* target_shell = shell_list.get_and_step();
            if(target_face->is_nonmanifold( (GroupingEntity*)target_shell ) )
            {         
                sheet_body_idealize_face.append(face_to_idealize[i]);
            }
        }
    }

    //if no faces to idealize that pass requirements, error out a warning message
    if(sheet_body_idealize_face.size()==0)
    {
        //I'm returning success here even though no surfaces found that meet shell requirements set above
        {
            PRINT_INFO( "Failed to find any feature(s) which met user requirements\n\n" );
        }
        return CUBIT_SUCCESS;
    }

    //temp_body_ptr = face_to_idealize[y]->body();
    //if(temp_body_ptr->is_sheet_body())
    //{
    //    sheet_body_idealize_face.append(face_to_idealize[y]);
    //}

    face_to_idealize.clean_out();

    //switching the DLIList<RefFace> to DLIList<Surface>
    GeometryModifyEngine* gme_ptr;
    DLIList <Body*> old_body_list;
    DLIList<Surface*> sheet_body_idealize_surface;
    gme_ptr = tweak_setup(sheet_body_idealize_face,"idealize",old_body_list,sheet_body_idealize_surface);
    sheet_body_idealize_face.clean_out();

    //grab all the loops from each sheet body surface
    DLIList <LoopSM*> idealize_loopSM_list;
    for(y=0;y<sheet_body_idealize_surface.size();y++)
    {
            sheet_body_idealize_surface[y]->loopsms(idealize_loopSM_list);
    }
    sheet_body_idealize_surface.clean_out();

    //this section is going to remove all excluded curves loopsm from the 'master' loopsm list
    if(exclude_entity.size()>0)
    {
        //cast the exclude DLIList<RefEntity> to DLIList<RefEdge>
        DLIList<RefEdge*> exclude_edge;
        CAST_LIST(exclude_entity, exclude_edge, RefEdge);
        
        //switching the DLIList<RefEdge> to DLIList<Curve>
        DLIList<Curve*> exclude_cuves;
        GeometryModifyEngine* gme_ptr1;
        gme_ptr1 = tweak_setup(exclude_edge,"idealize",old_body_list,exclude_cuves);
        exclude_edge.clean_out();

        //grabbing all the curve loops from the given excluded curves
        DLIList <LoopSM*> exclude_loops;
        for(y=0;y<exclude_cuves.size();y++)
        {
            exclude_cuves[y]->loopsms(exclude_loops);
        }
        exclude_cuves.clean_out();

        //remove the excluded loops from the list of sheet body loopsms
        idealize_loopSM_list -= exclude_loops;
    }

    //removing all the external loops from the list as they will not be tweak removed
    DLIList <LoopSM*> possible_internal_LoopSM_list;
    for(y=0;y<idealize_loopSM_list.size();y++)
    {
        if(idealize_loopSM_list[y]->loop_type() == LOOP_TYPE_HOLE)
        {
            possible_internal_LoopSM_list.append(idealize_loopSM_list[y]);
        }
    }
    idealize_loopSM_list.clean_out();
    DLIList <Curve*> hole_curves_to_remove;          //hole_curves_to_remove is the curves selected for removal out of the 'hole' search algorithm
    DLIList <Curve*> master_curve_remove_list;
    DLIList <LoopSM*> arc_LoopSM_list;
    DLIList <Surface*> temp_list;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //this begins the hole search algorithm
    //if no arc radius given, skip hole search algorithm
    if(arc_radius!=CUBIT_DBL_MAX)
    {
        DLIList <LoopSM*> not_hole_loop;        //loops which don't meet the all curves are arc type filter
        DLIList <Curve*> possible_internal_arcs;

        //search through possible internal curves filtering only for curves of type arc
        //if one of the curves in a loop is not an arc, add that loop to the not_hole_loop list
        for(y=0;y<possible_internal_LoopSM_list.size();y++)
        {
            possible_internal_arcs.clean_out();
            possible_internal_LoopSM_list[y]->curves(possible_internal_arcs);
            for(i=0;i<possible_internal_arcs.size();i++)
            {
                temp_list.clean_out();
                possible_internal_arcs[i]->surfaces(temp_list);
                //check whether or not curve is of arc type and whether it is attached to more than one surface
                if( possible_internal_arcs[i]->geometry_type() != ARC_CURVE_TYPE || temp_list.size() != 1)
                {               
                    not_hole_loop.append(possible_internal_LoopSM_list[y]);
                    break;
                }
            }
        }

        //change name of possible_internal_LoopSM_list to arc_LoopSM_list
        arc_LoopSM_list = possible_internal_LoopSM_list;
        //subtract from the possible loops the loops which don't have curves which are all arcs or are attached to more than two surfaces
        arc_LoopSM_list-=not_hole_loop;
        not_hole_loop.clean_out();

        //this next filter checks to make sure that all arcs of the same loop share the same
        //radius center within the geometry tolerance
        CubitVector arc_center_point, arc_center_point1;
        DLIList <LoopSM*> not_center_arc_loop;
        double rad_distance, temp_arc_radius , temp_arc_radius1;

        //this for loop is going to check that each loops arc curves have the same center radius point
        //if not you can remove that loop as a possibility for being added to the tweak_remove command
        for(y=0;y<arc_LoopSM_list.size();y++)
        {
            //clean out curve list before grabbing a new loop
            hole_curves_to_remove.clean_out();
            arc_LoopSM_list[y]->curves(hole_curves_to_remove);
            //iterate across the hole_curves_to_remove size
            for (i=0;i<hole_curves_to_remove.size();i++)
            {
                //if you are on the first index, we need to set a baseline radius point
                if(i==0)
                {
                    hole_curves_to_remove[i]->get_center_radius(arc_center_point,temp_arc_radius);
                    //if this is the only arc in the loop go ahead and check if it meets specified arc parameter
                    //if it doesn't meet the users parameter add the loop to the not_center_arc_loop list
                    if(temp_arc_radius >= arc_radius && hole_curves_to_remove.size()==1)
                    {
                        not_center_arc_loop.append(arc_LoopSM_list[y]);
                        break;
                    }
                }
                //now compare the other arc center points to the baseline, if it ever fails the users parameter
                //add the loop to the not_center_arc_loop list
                else
                {
                    hole_curves_to_remove[i]->get_center_radius(arc_center_point1,temp_arc_radius1);
                    rad_distance = arc_center_point.distance_between_squared(arc_center_point1);
                    if(rad_distance > geo_tol || temp_arc_radius >= arc_radius)
                    {
                        not_center_arc_loop.append(arc_LoopSM_list[y]);
                        break;
                    }
                }
            }
        }

        //remove loops which didn't have perfect circular holes from the arc_loopsm_list
        arc_LoopSM_list -= not_center_arc_loop;
        for(y=0;y<arc_LoopSM_list.size();y++)
        {
            arc_LoopSM_list[y]->curves(hole_curves_to_remove);
        }
        master_curve_remove_list+=hole_curves_to_remove;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //this begins the slot search algorithm
    DLIList<LoopSM*> removable_slot_loop;
    if(slot_arc_radius!=CUBIT_DBL_MAX || slot_length!=CUBIT_DBL_MAX)
    {
        DLIList<LoopSM*> four_curve_possible_slot;
        DLIList<LoopSM*> possible_slot;
        DLIList<Curve*> internal_curves_in_loop;

        //checks to make sure the loop has only four curves - may want to expand this in the future
        for(y=0;y<possible_internal_LoopSM_list.size();y++)
        {
            possible_internal_LoopSM_list[y]->curves(internal_curves_in_loop);
            if(internal_curves_in_loop.size()==4)
            {
                four_curve_possible_slot.append(possible_internal_LoopSM_list[y]);
            }
            internal_curves_in_loop.clean_out();
        }


        //check to make sure it alternates straight line, arc, etc...
        for(y=0;y<four_curve_possible_slot.size();y++)
        {
            four_curve_possible_slot[y]->curves(internal_curves_in_loop);

            if(internal_curves_in_loop[0]->geometry_type() == ARC_CURVE_TYPE &&
                internal_curves_in_loop[1]->geometry_type() == STRAIGHT_CURVE_TYPE &&
                internal_curves_in_loop[2]->geometry_type() == ARC_CURVE_TYPE &&
                internal_curves_in_loop[3]->geometry_type() == STRAIGHT_CURVE_TYPE)
            {
                int num_of_surfs=0;
                for(i=0;i<internal_curves_in_loop.size();i++)
                {
                    temp_list.clean_out();
                    internal_curves_in_loop[i]->surfaces(temp_list);
                    num_of_surfs=num_of_surfs + temp_list.size();
                }
                if(num_of_surfs==4)
                {
                    possible_slot.append(four_curve_possible_slot[y]);
                }
            }
            else if(internal_curves_in_loop[0]->geometry_type() == STRAIGHT_CURVE_TYPE &&
                internal_curves_in_loop[1]->geometry_type() == ARC_CURVE_TYPE &&
                internal_curves_in_loop[2]->geometry_type() == STRAIGHT_CURVE_TYPE &&
                internal_curves_in_loop[3]->geometry_type() == ARC_CURVE_TYPE)
            {
                int num_of_surfs=0;
                for(i=0;i<internal_curves_in_loop.size();i++)
                {
                    temp_list.clean_out();
                    internal_curves_in_loop[i]->surfaces(temp_list);
                    num_of_surfs=num_of_surfs + temp_list.size();
                }
                if(num_of_surfs==4)
                {
                    possible_slot.append(four_curve_possible_slot[y]);
                }
            }
            internal_curves_in_loop.clean_out();
        }

        CubitVector arc_center_point;
        double temp_arc_radius = CUBIT_DBL_MAX, curve_length = CUBIT_DBL_MAX;

        //check to make sure that the rad and/or length meet users parameters
        for(y=0;y<possible_slot.size();y++)
        {
            possible_slot[y]->curves(internal_curves_in_loop);
            //if user specified rad, then passed rad_counter should = 2 after for loop completes
            //if length specified, length_counter should =2 after for loop completes
            int rad_counter = 0, length_counter = 0;
            for(i=0;i<internal_curves_in_loop.size();i++)
            {
                //if curve is an arc and user specified a radius enter if statement
                if( internal_curves_in_loop[i]->geometry_type() == ARC_CURVE_TYPE && slot_arc_radius!=CUBIT_DBL_MAX )
                {
                    //check the radius against the user inputed value
                    internal_curves_in_loop[i]->get_center_radius(arc_center_point,temp_arc_radius);
                    if(temp_arc_radius <= slot_arc_radius)
                    {
                        //if it passes rad test, add to rad_counter
                        rad_counter++;
                    }
                }
                else if(internal_curves_in_loop[i]->geometry_type() == STRAIGHT_CURVE_TYPE && slot_length!=CUBIT_DBL_MAX )
                {
                    //check the length against the user inputed value
                    curve_length = internal_curves_in_loop[i]->get_arc_length();
                    if(curve_length <= slot_length)
                    {
                        //if it passes rad test, add to length_counter
                        length_counter++;
                    }
                }
            }

            //checks that if user specified length and radius constraint that its parameter passes for all four curves
            if(slot_length!=CUBIT_DBL_MAX && slot_arc_radius!=CUBIT_DBL_MAX  && rad_counter==2 && length_counter==2)
            {
                removable_slot_loop.append(possible_slot[y]);
            }

            //if user only specified one length or arc parameter, it only needs to meet 2 criteria
            else if((slot_length!=CUBIT_DBL_MAX  && length_counter==2) || (slot_arc_radius!=CUBIT_DBL_MAX  && rad_counter==2))
            {
                removable_slot_loop.append(possible_slot[y]);
            }
            internal_curves_in_loop.clean_out();
        }
        //add removable loops curves to the master_curve_remove_list list
        for(y=0;y<removable_slot_loop.size();y++)
        {
            removable_slot_loop[y]->curves(master_curve_remove_list);
        }
    }

    //if no arcs are found to be removed, warn the user.
    if(master_curve_remove_list.size()==0)
    {
        //I'm returning success here even though no curves were found
        {
            PRINT_INFO( "Failed to find any feature(s) which met user requirements\n\n" );
        }
        return CUBIT_SUCCESS;
    }
    else if(preview == CUBIT_TRUE)
    {
        GfxPreview::clear();

        for(i=0; i<master_curve_remove_list.size();i++)
        {
            CubitStatus result;
            GMem g_mem;

            // get the graphics
            result = master_curve_remove_list[i]->get_geometry_query_engine()->
                get_graphics( master_curve_remove_list[i], &g_mem );

            if (result==CUBIT_FAILURE || g_mem.pointListCount == 0)
            {
                PRINT_WARNING("Unable to preview a curve\n" );
                double len = master_curve_remove_list[i]->
                    length_from_u(master_curve_remove_list[i]->start_param(),master_curve_remove_list[i]->end_param());

                PRINT_WARNING("Curve len: %f\n",len);
            }

            // Draw the polyline
            GfxPreview::draw_polyline( g_mem.point_list(), g_mem.pointListCount, CUBIT_BLUE );

        }

        //output the number of holes or slots which were found
        PRINT_INFO("Found %d holes and %d slots which met idealization parameters\n\n", arc_LoopSM_list.size(),removable_slot_loop.size());
        GfxPreview::flush();
        return CUBIT_SUCCESS;
    }
    else
    {
        DLIList<BodySM*> new_bodysm_list;
        //pass master_curve_remove_list to the tweak_remove command
        CubitStatus stat = gme_ptr->tweak_remove(master_curve_remove_list, new_bodysm_list,CUBIT_FALSE, CUBIT_FALSE );

        //update DAG
        DLIList<Body*> new_body_list;
        stat = finish_sm_op( old_body_list, new_bodysm_list ,new_body_list );
        //output the number of holes or slots which were found
        PRINT_INFO("Found %d holes, and %d slots which met idealization parameters\n\n", arc_LoopSM_list.size(), removable_slot_loop.size());
        return CUBIT_SUCCESS;
    }
}

//=============================================================================
// Description: Create drop down surfaces from the external edges of a surface to another surface with options
//	        Use it primarially when you have a doubler plate situation and you need
//		to connect them together with surfaces to represent a weld
// Author     : Jonathan Bugman
// Date       : 02/07/2008
//=============================================================================
CubitStatus GeometryModifyTool::create_surface_doubler(DLIList<RefEntity*> doubler_entity,
                                                       DLIList<RefEntity*> target_entity,
                                                       DLIList<Body*> &body_list_out,
                                                       CubitBoolean internal_flg,
                                                       CubitBoolean extend_flg,
                                                       CubitPlane *limit_plane,
                                                       CubitVector sweep_direction,
                                                       CubitBoolean preview)
{
	//need to switch the DLIList<RefEntity> to a DLIList<RefFace>
	DLIList<RefFace*> doubler_face;
	DLIList<RefFace*> target_face;
	CAST_LIST( doubler_entity, doubler_face, RefFace);
	CAST_LIST( target_entity, target_face, RefFace);

	DLIList <Surface*> doubler_surface;
	DLIList <Surface*> target_surface;
	DLIList <Body*> old_body_list;
	DLIList<RefFace*> tweak_face;
	DLIList<Surface*> tweak_surface;
	tweak_face+=doubler_face;
	tweak_face+=target_face;
	GeometryModifyEngine* gme_ptr;
    int ii=0, i=0;

	gme_ptr = tweak_setup(tweak_face,"doubler",old_body_list,tweak_surface);
    int z;
	for(z=0;z<doubler_face.size();z++)
	{
		doubler_surface.append(tweak_surface[z]);
	}
	for(z=0;z<tweak_face.size();z++)
    {
        target_surface.append(tweak_surface[z]);
    }

    DLIList<BodySM*> all_kept_bodies;
    DLIList<BodySM*> body_convert;
    DLIList<Surface*> copied_doubler_surface;
    DLIList <BodySM*> tweak_target_bodySM;

    for(z=0;z<doubler_surface.size();z++)
    {
        copied_doubler_surface.append(gme_ptr->make_Surface(doubler_surface[z]));
        body_convert.append(copied_doubler_surface[z]->bodysm());
    }

    //each workflow is dependent on whether or not a sweep_direction is specified
    if(sweep_direction == CubitVector(0,0,0))
    {
		DLIList<BodySM*> united_bodies;
		DLIList<BodySM*> separate_bodies;
		//this section takes all the doublers, unites them, and then splits them.  If only one body
		//then skip the aforementioned two steps.
		if(doubler_surface.size()==1)
		{
			separate_bodies=body_convert;
		}
		else
		{
			if(gme_ptr->unite(body_convert,united_bodies) == CUBIT_FAILURE || united_bodies.size()==0 )
			{
				PRINT_ERROR( "Command failed at unite command\n" );
				return CUBIT_FAILURE;
			}

			if(gme_ptr->split_body(united_bodies[0],separate_bodies) == CUBIT_FAILURE || separate_bodies.size()==0)
			{
				PRINT_ERROR( "Command failed at separate command\n" );
				return CUBIT_FAILURE;
			}
		}
		for(z=0; z<separate_bodies.size();z++)
		{
			DLIList<Surface*> body_surface;
			separate_bodies[z]->surfaces(body_surface);

			DLIList<CubitVector> doubler_surface_center_points;
			int d=0;
			for(d=0;d<body_surface.size();d++)
			{
				doubler_surface_center_points.append(body_surface[d]->bounding_box().center());
			}
			CubitVector doubler_normal,doubler_center_pt,doubler_to_target_vector,target_center_pt;
			double extrude_distance = 0.0;

			//make sure that the normal of the surface is pointing towards the target surface
			//the thicken command thickens in the direction of the surface normal, normals are checked using dot product check
			CubitVector center_pt = body_surface[0]->bounding_box().center();
			body_surface[0]->closest_point(center_pt,&doubler_center_pt,&doubler_normal);
            //adding this for loop because of l bracket doublers may have the first target surface perpendicular to an opposite side doubler surface
            //resulting in the code erroneously failing
            double dot=0.0;
            int mm;
            for(mm=0; mm<target_surface.size();mm++)
            {
                target_surface[mm]->closest_point_trimmed(doubler_center_pt, target_center_pt);
                doubler_to_target_vector = target_center_pt - doubler_center_pt;
                dot = doubler_to_target_vector.x()*doubler_normal.x()+doubler_to_target_vector.y()*doubler_normal.y()+doubler_to_target_vector.z()*doubler_normal.z();
                if(fabs(dot)>1E-6)
                {
                    mm=target_surface.size();
                }
            }
            if(fabs(dot)<1E-6)
            {
                PRINT_ERROR( "Doubler and target surface are touching or are perpendicular to each other\n" );
                for(ii =0;ii<separate_bodies.size();ii++)
                {
                    GeometryQueryEngine* gqe = separate_bodies[ii]->get_geometry_query_engine();
                    gqe->delete_solid_model_entities(separate_bodies[ii]);
                }
                return CUBIT_FAILURE;
            }
            else if(dot < 0)
            {
                if(gme_ptr->reverse_body(separate_bodies[z])==CUBIT_FAILURE)
                {
                    PRINT_ERROR( "Command failed at reverse body command.\n" );
                    for(ii =0;ii<separate_bodies.size();ii++)
                    {
                        GeometryQueryEngine* gqe = separate_bodies[ii]->get_geometry_query_engine();
                        gqe->delete_solid_model_entities(separate_bodies[ii]);
                    }
                    return CUBIT_FAILURE;
                }
                extrude_distance = 0.0001;
            }
            else
			{
				extrude_distance = 0.0001;
			}

			DLIList<BodySM*> thickened_doubler_bodySM;
			DLIList<BodySM*> current_body;
			current_body.append(separate_bodies[z]);

            if(gme_ptr->thicken(current_body,thickened_doubler_bodySM,extrude_distance,CUBIT_FALSE) == CUBIT_FAILURE || thickened_doubler_bodySM.size()==0)
            {
                PRINT_ERROR( "Command failed at thicken command, this may be due to using a non-ACIS geometry engine\n" );
                for(ii =0;ii<separate_bodies.size();ii++)
                {
                    GeometryQueryEngine* gqe = separate_bodies[ii]->get_geometry_query_engine();
                    gqe->delete_solid_model_entities(separate_bodies[ii]);
                }
                return CUBIT_FAILURE;
            }

			//need to grab the newly created surface opposite the user selected one from the thicken function to carry it through for the tweak target
			DLIList<Surface*> thicken_surfaces;
			thickened_doubler_bodySM[0]->surfaces(thicken_surfaces);
			DLIList <Surface*> post_thicken_doublers;

			int y=0;
			for(y=0;y<doubler_surface_center_points.size();y++)
			{
				doubler_center_pt = doubler_surface_center_points[y];
				int r=0;
				for(r=0;r<thicken_surfaces.size();r++)
				{
					CubitVector test_center_pt = thicken_surfaces[r]->bounding_box().center();
					if((test_center_pt-doubler_center_pt).length()<=.000001)
					{
						post_thicken_doublers.append(thicken_surfaces[r]);
					}
				}
			}

			DLIList <LoopSM*> doubler_loopSM_list;
			DLIList <Curve*> doubler_external_curves;
			Curve* test_external_curves1;
			Curve* test_external_curves2;

			//need to do this in order to grab all curves, not just external IMPORTANT:THIS HAS TO BE DONE BEFORE THE tweak? COMMAND!
			for(y=0;y<post_thicken_doublers.size();y++)
			{
				post_thicken_doublers[y]->loopsms(doubler_loopSM_list);
			}
			for(i=0;i<doubler_loopSM_list.size();i++)
			{
				doubler_loopSM_list[i]->curves(doubler_external_curves);
			}

			doubler_loopSM_list.clean_out();
            tweak_target_bodySM.clean_out();
			DLIList <LoopSM*> test_loopSM_list;
			DLIList <Curve*> thicken_external_curves;
			DLIList <Surface*> tweak_target_surface = thicken_surfaces;

            //stepping through the surfaces from the thicken body
            for(i=0; i < thicken_surfaces.size(); i++)
            {
                thicken_surfaces[i]->loopsms(test_loopSM_list);
                //grabbing the external curves from the current thicken_surface
                test_loopSM_list[0]->curves(thicken_external_curves);
                test_loopSM_list.clean_out();
                int j=0;
                for(j=0;j<thicken_external_curves.size();j++)
                {
                    //step through the first curve
                    test_external_curves1 = thicken_external_curves[j];
                    int k=0;
                    for(k=0; k<doubler_external_curves.size();k++)
                    {
                        //while stepping through the doubler plate curves, compare them to the test_test_surface curves
                        test_external_curves2 = doubler_external_curves[k];

						//if the two are equal, they are touching the doulber and therefore are either the side surfaces or the doubler
						if(test_external_curves2 == test_external_curves1)
						{
							//remove the surface from the tweak_target_surface list
							tweak_target_surface.remove_all_with_value(thicken_surfaces[i]);
							break;
						}
					}
					if(test_external_curves2 == test_external_curves1)
					{
						break;
					}

				}
				thicken_external_curves.clean_out();
            }

            //pass the found opposite surface into the tweak_target routine
            if(gme_ptr->tweak_target(tweak_target_surface,target_surface,tweak_target_bodySM,extend_flg,limit_plane) == CUBIT_FAILURE || tweak_target_bodySM.size()==0)
            {
                PRINT_ERROR( "Command failed at Tweak_Target routine\n" );
                for(ii =0;ii<thickened_doubler_bodySM.size();ii++)
                {
                    GeometryQueryEngine* gqe = thickened_doubler_bodySM[ii]->get_geometry_query_engine();
                    gqe->delete_solid_model_entities(thickened_doubler_bodySM[ii]);
                }
                return CUBIT_FAILURE;
            }

			//fill out a tweak_body_surface list from tweak_target routine
			DLIList<Surface*> tweak_body_surfaces;
			tweak_target_bodySM[0]->surfaces(tweak_body_surfaces);
			DLIList <Curve*> tweak_external_curves;
			doubler_external_curves.clean_out();
			
			//refilling DLIList's as needed based on internal_flg
			//if we are not keeping internal surfaces we do not want it's curves in the doubler_external_curves list
			//otherwise if we are, we do want the curves in the list for the following sections for loop
            if(internal_flg==CUBIT_FALSE)
            {
                int j=0;
                for(i=0;i<post_thicken_doublers.size();i++)
                {
                    post_thicken_doublers[i]->loopsms(doubler_loopSM_list);
                    for(j=0;j<doubler_loopSM_list.size();j++)
                    {
                      LoopType loop_type = doubler_loopSM_list[j]->loop_type();
                      if(loop_type == LOOP_TYPE_EXTERNAL ||
                         loop_type == LOOP_TYPE_U_PERIODIC ||
                         loop_type == LOOP_TYPE_V_PERIODIC)
                      {
                          doubler_loopSM_list[j]->curves(doubler_external_curves);
                          break;
                      }
                    }
                    doubler_loopSM_list.clean_out();
                }
            }
			else
			{
				for(i=0;i<post_thicken_doublers.size();i++)
				{
					post_thicken_doublers[i]->loopsms(doubler_loopSM_list);
				}
				for(i=0;i<doubler_loopSM_list.size();i++)
				{
					doubler_loopSM_list[i]->curves(doubler_external_curves);
				}

			}

            DLIList <Surface*> surfaces_to_keep;
            for(i=0;i<tweak_body_surfaces.size();i++)
            {
                tweak_body_surfaces[i]->loopsms(test_loopSM_list);

                if(test_loopSM_list.size()==0)
                {
                    PRINT_ERROR( "Command failed to find any doubler drop down curves\n" );
                    for(ii =0;ii<thickened_doubler_bodySM.size();ii++)
                    {
                        GeometryQueryEngine* gqe = thickened_doubler_bodySM[ii]->get_geometry_query_engine();
                        gqe->delete_solid_model_entities(thickened_doubler_bodySM[ii]);
                    }
                    return CUBIT_FAILURE;
                }

                test_loopSM_list[0]->curves(tweak_external_curves);
                test_loopSM_list.clean_out();

				int j=0;
				for(j=0;j<tweak_external_curves.size();j++)
				{
					test_external_curves1 = tweak_external_curves[j];

					int k=0;
					for(k=0; k<doubler_external_curves.size();k++)
					{
						//while stepping through the doubler plate curves, compare them to the test_loop_list
						test_external_curves2 = doubler_external_curves[k];

						if(test_external_curves2 == test_external_curves1)
						{
							surfaces_to_keep.append(tweak_body_surfaces[i]);
							break;
						}
					}
					if(test_external_curves2 == test_external_curves1)
					{
						break;
					}
				}
				tweak_external_curves.clean_out();
			}

            if(surfaces_to_keep.size()==0)
            {
                PRINT_ERROR( "Failed to find and keep surfaces\n" );
                for(ii =0;ii<tweak_target_bodySM.size();ii++)
                {
                    GeometryQueryEngine* gqe = tweak_target_bodySM[ii]->get_geometry_query_engine();
                    gqe->delete_solid_model_entities(tweak_target_bodySM[ii]);
                }
                return CUBIT_FAILURE;
            }

			//do this to remove the copied_doubler_surface since we no longer need the surface anymore
			int c=0;
			for(c=0;c<post_thicken_doublers.size();c++)
			{
				surfaces_to_keep.remove_all_with_value(post_thicken_doublers[c]);
			}

			DLIList <Surface*> surfaces_to_remove = tweak_body_surfaces;
			surfaces_to_remove -= surfaces_to_keep;
			DLIList<BodySM*> resulting_bodies;

			//remove all surfaces in the surfaces_to_remove list
            if(gme_ptr->tweak_remove(surfaces_to_remove,resulting_bodies,CUBIT_FALSE) == CUBIT_FAILURE || resulting_bodies.size()==0)
            {
                PRINT_ERROR( "Command failed at Tweak_Remove routine\n" );
                for(ii =0;ii<tweak_target_bodySM.size();ii++)
                {
                    GeometryQueryEngine* gqe = tweak_target_bodySM[ii]->get_geometry_query_engine();
                    gqe->delete_solid_model_entities(tweak_target_bodySM[ii]);
                }
                return CUBIT_FAILURE;
            }
            all_kept_bodies+=resulting_bodies;
        }
    }	
	else
	{
		DLIList<BodySM*> swept_bodies;
		DLIList<BodySM*> swept_doubler_bodySM;

		//take the copied_doubler_surface and extrude it along the sweep_direction to create a body
		for(z=0;z<copied_doubler_surface.size();z++)
		{
			DLIList<GeometryEntity*> DLIList_copied_doubler_surface;
            DLIList_copied_doubler_surface.append(copied_doubler_surface[z]);
            if(gme_ptr->sweep_translational(DLIList_copied_doubler_surface,swept_doubler_bodySM,sweep_direction*0.0001,0.0,0,CUBIT_FALSE,CUBIT_FALSE) == CUBIT_FAILURE || swept_doubler_bodySM.size()==0)
            {
                PRINT_ERROR( "Command failed at sweep->extrude command\n" );
                for(ii =0;ii<body_convert.size();ii++)
                {
                    GeometryQueryEngine* gqe = body_convert[ii]->get_geometry_query_engine();
                    gqe->delete_solid_model_entities(body_convert[ii]);
                }
                return CUBIT_FAILURE;
            }
			swept_bodies+=swept_doubler_bodySM;
			swept_doubler_bodySM.clean_out();
		}

		DLIList<BodySM*> united_bodies;
		DLIList<BodySM*> separate_bodies;
		//if more than one body, unite and split the newly created bodies, if only one body skip this step
		//as the unite will fail
		if(swept_bodies.size()==1)
		{
			separate_bodies=swept_bodies;
		}
		else
		{
			if(gme_ptr->unite(swept_bodies,united_bodies) == CUBIT_FAILURE || united_bodies.size()==0 )
			{
				PRINT_ERROR( "Command failed at unite command\n" );
				return CUBIT_FAILURE;
			}

			if(gme_ptr->split_body(united_bodies[0],separate_bodies) == CUBIT_FAILURE || separate_bodies.size()==0)
			{
				PRINT_ERROR( "Command failed at separate command\n" );
				return CUBIT_FAILURE;
			}
		}

		//create another copy of copied_doubler_surface since copied_doubler_surface is going to be manipulated
		DLIList<Surface*> temp_copied_doubler_surface=copied_doubler_surface;
		for(z=0;z<separate_bodies.size();z++)
		{
			//need to grab the newly created surface opposite the user selected one from the thicken function to carry it through for the tweak target
			//this will need to be changed to a for loop to account for multiple thickened bodies if we impliment multiple doubler surfaces
			DLIList<Surface*> thicken_surfaces;
			separate_bodies[z]->surfaces(thicken_surfaces);

			//initializing a lot of variables to be used in the next few steps
			DLIList<Surface*> master_surface_remove_list;
			DLIList <Curve*> thicken_external_curves;
			DLIList<Surface*> tweak_target_surface;
			DLIList<Surface*> surfaces_to_remove;

			//using a centerpoint of the surfaces, I want to find out which surface from the recently swept bodies corresponds to the surface of the body closest the target
			//this has to be done because sweep_translational moves the source surface.  Thicken on the otherhand does and doesn't based on number of surfaces being thickened.
			int y=0;
			for(y=0;y<temp_copied_doubler_surface.size();y++)
			{
				CubitVector doubler_center_pt = temp_copied_doubler_surface[y]->bounding_box().center();
				int r=0;
				for(r=0;r<thicken_surfaces.size();r++)
				{
					CubitVector test_center_pt = thicken_surfaces[r]->bounding_box().center();
					if((test_center_pt-doubler_center_pt).length()<=.000001)
					{
						tweak_target_surface.append(thicken_surfaces[r]);
						surfaces_to_remove.append(temp_copied_doubler_surface[y]);
					}
				}
			}
			//remove the manipulated surfaces from the temp_copied_doubler_surface list
			temp_copied_doubler_surface-=surfaces_to_remove;
			surfaces_to_remove.clean_out();

			//grab all the curves of the doubler
			DLIList <Curve*> doubler_external_curves;
			DLIList <LoopSM*> doubler_loopSM_list;
			for(y=0;y<tweak_target_surface.size();y++)
			{
				tweak_target_surface[y]->loopsms(doubler_loopSM_list);
			}

			for(i=0;i<doubler_loopSM_list.size();i++)
			{
				doubler_loopSM_list[i]->curves(doubler_external_curves);
			}

			Surface* temp_test_surface;
			DLIList <Surface*> doubler_surface_for_this_body = thicken_surfaces;
			Curve* test_external_curves1;
			Curve* test_external_curves2;
			DLIList <LoopSM*> test_loopSM_list;

			for(i=0; i < thicken_surfaces.size(); i++)
			{
				//step through the thickened bodies surfaces
				thicken_surfaces[i]->loopsms(test_loopSM_list);
				temp_test_surface = thicken_surfaces[i];
				//grabbing the external curve loop from the face and making it a DLIList <RefEdge*>
				test_loopSM_list[0]->curves(thicken_external_curves);

				int j=0;
				for(j=0;j<thicken_external_curves.size();j++)
				{
					//step through the loop list
					test_external_curves1 = thicken_external_curves[j];

					int k=0;
					for(k=0; k<doubler_external_curves.size();k++)
					{
						//while stepping through the doubler plate curves, compare them to the test_loop_list
						test_external_curves2 = doubler_external_curves[k];

						if(test_external_curves2 == test_external_curves1)
						{
							doubler_surface_for_this_body.remove_all_with_value(thicken_surfaces[i]);
							break;
						}
					}
					if(test_external_curves2 == test_external_curves1)
					{
						break;
					}
				}

				thicken_external_curves.clean_out();
				thicken_external_curves.reset();
                test_loopSM_list.clean_out();
            }

            //DLIList <BodySM*> tweak_target_bodySM
            tweak_target_bodySM.clean_out();
            if(gme_ptr->tweak_target(tweak_target_surface,target_surface,tweak_target_bodySM,extend_flg,limit_plane) == CUBIT_FAILURE || tweak_target_bodySM.size()==0)
            {
                PRINT_ERROR( "Command failed at Tweak_Target routine\n" );
                for(ii =0;ii<body_convert.size();ii++)
                {
                    GeometryQueryEngine* gqe = body_convert[ii]->get_geometry_query_engine();
                    gqe->delete_solid_model_entities(body_convert[ii]);
                }
                return CUBIT_FAILURE;
            }

			//make note of the doubler_surface and add it to the delete list
			master_surface_remove_list+=doubler_surface_for_this_body;

			//clean out these DLIList's as they will be used later on
			doubler_loopSM_list.clean_out();
			doubler_external_curves.clean_out();

			//refilling DLIList's as needed based on internal_flg
			//basically if surfaces share a curve, that surface will be kept so if you don't want the internal surfaces
			//create the list without the internal curves and they'll be removed
            if(internal_flg==CUBIT_FALSE)
            {
                int j=0;
                for(i=0;i<doubler_surface_for_this_body.size();i++)
                {
                    doubler_surface_for_this_body[i]->loopsms(doubler_loopSM_list);
                    for(j=0;j<doubler_loopSM_list.size();j++)
                    {
                      LoopType loop_type = doubler_loopSM_list[j]->loop_type();
                      if(loop_type == LOOP_TYPE_EXTERNAL ||
                         loop_type == LOOP_TYPE_U_PERIODIC ||
                         loop_type == LOOP_TYPE_V_PERIODIC)
                      {
                          doubler_loopSM_list[j]->curves(doubler_external_curves);
                          break;
                      }

                    }
                    doubler_loopSM_list.clean_out();
                }
            }
            else
			{
				for(i=0;i<doubler_surface_for_this_body.size();i++)
				{
					doubler_surface_for_this_body[i]->loopsms(doubler_loopSM_list);
				}
				for(i=0;i<doubler_loopSM_list.size();i++)
				{
					doubler_loopSM_list[i]->curves(doubler_external_curves);
				}
			}

			//recreate the thicken_surfaces list based now on the bodies after the tweak_target command
			thicken_surfaces.clean_out();
			tweak_target_bodySM[0]->surfaces(thicken_surfaces);

			DLIList <Surface*> surfaces_to_keep;
			surfaces_to_remove = thicken_surfaces;
			DLIList <Curve*> tweak_external_curves;
			
			for(i=0;i<thicken_surfaces.size();i++)
			{
				thicken_surfaces[i]->loopsms(test_loopSM_list);
				//grabs the external curves from face
				test_loopSM_list[0]->curves(tweak_external_curves);

				int j=0;
				for(j=0;j<tweak_external_curves.size();j++)
				{
					//step through the loop list
					test_external_curves1 = tweak_external_curves[j];

					int k=0;
					for(k=0; k<doubler_external_curves.size();k++)
					{
						//while stepping through the doubler plate curves, compare them to the test_loop_list
						test_external_curves2 = doubler_external_curves[k];

						if(test_external_curves2 == test_external_curves1)
						{
							surfaces_to_keep.append(thicken_surfaces[i]);
							break;
						}
					}
					if(test_external_curves2 == test_external_curves1)
					{
						break;
					}
				}
				test_loopSM_list.clean_out();
				tweak_external_curves.clean_out();
			}

			if(surfaces_to_keep.size()==0)
			{
				PRINT_ERROR( "Failed to find and keep tweak_target surfaces\n" );
                for(ii =0;ii<tweak_target_bodySM.size();ii++)
                {
                    GeometryQueryEngine* gqe = tweak_target_bodySM[ii]->get_geometry_query_engine();
                    gqe->delete_solid_model_entities(tweak_target_bodySM[ii]);
                }
				return CUBIT_FAILURE;
			}

			//subtract surfaces
			surfaces_to_remove -= surfaces_to_keep;
			master_surface_remove_list+=surfaces_to_remove;

			DLIList<BodySM*> resulting_bodies;

			if(gme_ptr->tweak_remove(master_surface_remove_list,resulting_bodies,CUBIT_FALSE) == CUBIT_FAILURE || resulting_bodies.size()==0)
			{
				PRINT_ERROR( "Command failed at Tweak_Remove routine\n" );
                for(ii =0;ii<tweak_target_bodySM.size();ii++)
                {
                    GeometryQueryEngine* gqe = tweak_target_bodySM[ii]->get_geometry_query_engine();
                    gqe->delete_solid_model_entities(tweak_target_bodySM[ii]);
                }
				return CUBIT_FAILURE;
			}

			//all_kept_bodies is a list of bodies that will eventually be passed into finish_sm_op at the end
			all_kept_bodies+=resulting_bodies;
		}
    }

    if(preview==CUBIT_FALSE)
    {
        DLIList<BodySM*> bodies_to_unite;
        //check to see if their is only one body.  If only one body skip over the unite and split because
        //the unite command will fail (there is a check at the beginning to return cubit_failure)
        //append the original doubler surfaces to the resulting body list

        for(i=0;i<doubler_surface.size();i++)
        {
            all_kept_bodies.insert_first(doubler_surface[i]->bodysm());
        }
        if(all_kept_bodies.size()!=1)
        {
            if(gme_ptr->unite(all_kept_bodies,bodies_to_unite) == CUBIT_FAILURE || bodies_to_unite.size()==0 )
            {
                PRINT_ERROR( "Command failed at unite command\n" );
                return CUBIT_FAILURE;
            }
            all_kept_bodies.clean_out();
            if(gme_ptr->split_body(bodies_to_unite[0],all_kept_bodies) == CUBIT_FAILURE || all_kept_bodies.size()==0)
            {
                PRINT_ERROR( "Command failed at separate command\n" );
                return CUBIT_FAILURE;
            }
        }
        else
        {
                PRINT_WARNING( "Command may have failed at finding doubler surface(s) and appending them to the drop-down surfaces\n" );
                return CUBIT_FAILURE;
        }
        
        //update DAG
		CubitStatus stat;
		stat = finish_sm_op( old_body_list, all_kept_bodies ,body_list_out );
		return CUBIT_SUCCESS;
	}
	else
	{
		DLIList<Curve*> kept_curves;
		for(i =0;i<all_kept_bodies.size();i++)
		{
			all_kept_bodies[i]->curves(kept_curves);
		}

        GfxPreview::clear();

        for(i=0; i<kept_curves.size();i++)
        {
            CubitStatus result;
            GMem g_mem;

            // get the graphics
            result = kept_curves[i]->get_geometry_query_engine()->
                get_graphics( kept_curves[i], &g_mem );

            if (result==CUBIT_FAILURE || g_mem.pointListCount == 0)
            {
                PRINT_WARNING("Unable to preview a curve\n" );;
                double len = kept_curves[i]->
                    length_from_u(kept_curves[i]->start_param(),kept_curves[i]->end_param());

                PRINT_WARNING("Curve len: %f\n",len);
            }

            // Draw the polyline
            GfxPreview::draw_polyline( g_mem.point_list(), g_mem.pointListCount, CUBIT_BLUE );
        }
        GfxPreview::flush();
        for(ii =0;ii<tweak_target_bodySM.size();ii++)
        {
            GeometryQueryEngine* gqe = tweak_target_bodySM[ii]->get_geometry_query_engine();
            gqe->delete_solid_model_entities(tweak_target_bodySM[ii]);
        }
        return CUBIT_SUCCESS;
    }
}

//=============================================================================
// Function   : tweak_bend
// Member Type: PUBLIC
// Description: Bend solid bodies based on a bend radius and angle
// Author     : Sam Showman
// Date       : 06/23/08
//=============================================================================
CubitStatus GeometryModifyTool::tweak_bend( DLIList<Body*> &bend_bodies,
                                            DLIList<Body*> &new_body_list,
                                            CubitVector& neutral_root,
                                            CubitVector& bend_axis,
                                            CubitVector& bend_direction,
                                            double radius,
                                            double angle,
                                            DLIList<CubitVector*> bend_regions,
                                            double width,
                                            CubitBoolean center_bend,
                                            int num_points,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview )
{
    if (CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE)
	{
		if (keep_old_body)
			CubitUndo::save_state();
		else
			CubitUndo::save_state_with_cubit_file( bend_bodies );
    }
    
    DLIList<BodySM*> new_body_sm_list;
    DLIList<BodySM*> bend_bodies_sm;
    GeometryModifyEngine* engine;
    engine = common_modify_engine(bend_bodies, bend_bodies_sm);
    
    if (!preview)
    {
        //do_attribute_setup();
        //push_vg_attributes_before_modify(bend_bodies_sm);
    }

    CubitStatus result = engine->tweak_bend(
            bend_bodies_sm,
            new_body_sm_list,
            neutral_root,
            bend_axis,
            bend_direction,
            radius,
            angle,
            bend_regions,
            width,
            center_bend,
            num_points,
            keep_old_body,
            preview);

    if (result == CUBIT_FAILURE)
    {
        if (!preview)
        {
            //remove_pushed_attributes(bend_bodies_sm, bend_bodies);
            //do_attribute_cleanup();
        }
        if (CubitUndo::get_undo_enabled())
            CubitUndo::remove_last_undo();
        return CUBIT_FAILURE;
    }

    if (preview == CUBIT_FALSE)
    {
        //restore_vg_after_modify(new_body_sm_list, bend_bodies, engine);
        //remove_pushed_attributes(new_body_sm_list, bend_bodies);

        // Update DAG
        CubitStatus stat = finish_sm_op( bend_bodies, new_body_sm_list, new_body_list );
        if (CubitUndo::get_undo_enabled())
        {
            if (stat == CUBIT_SUCCESS)
                CubitUndo::note_result_bodies( new_body_list );
            else
                CubitUndo::remove_last_undo();
        }

        //do_attribute_cleanup();

        // Update graphics
        /*
        bend_bodies.reset();
        int i = 0;
        for (i = 0; i < bend_bodies.size(); i++)
        {
            Body* body_ptr = bend_bodies.get_and_step();
            body_ptr->notify_all_observers( GEOMETRY_MODIFIED );
        }//*/
        
        // get list of entities to update
        DLIList<RefEntity*> entities_to_update;
        int i;
        for(i=0; i < new_body_sm_list.size(); i++)
        {
            BodySM* bodysm = new_body_sm_list.get_and_step();
            DLIList<TopologyBridge*> to_check;
            DLIList<TopologyBridge*> tmp;
            DLIList<Surface*> surfs;
            bodysm->surfaces(surfs);
            DLIList<Curve*> curves;
            bodysm->curves(curves);
            DLIList<Point*> points;
            bodysm->points(points);
            to_check.append(bodysm);
            to_check.append(bodysm->lump());
            CAST_LIST_TO_PARENT(surfs, tmp);
            to_check += tmp;
            CAST_LIST_TO_PARENT(curves, tmp);
            to_check += tmp;
            CAST_LIST_TO_PARENT(points, tmp);
            to_check += tmp;

            int k;
            for(k=0; k<to_check.size(); k++)
                if(BridgeManager* m = to_check.get_and_step()->bridge_manager())
                    if(TopologyEntity* t = m->topology_entity())
                        entities_to_update.append(CAST_TO(t, RefEntity));

        }
        // Update graphics
        while (entities_to_update.size())
            entities_to_update.pop()->notify_all_observers(GEOMETRY_MODIFIED);

        return stat;
    }

    return CUBIT_SUCCESS;
}



//=============================================================================
// Description: Chamfer curves on solid and sheet bodies.  The left and right
//              offsets are with respect to the curve direction.  If the given
//              right offset is negative, the left offset is used.  Users can
//              preview to clarify the meaning of left and right.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_chamfer( DLIList<RefEdge*> &ref_edge_list,
                                               double left_offset,
                                               DLIList<Body*> &new_body_list,
                                               double right_offset,
                                               CubitBoolean keep_old_body,
                                               CubitBoolean preview )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_edge_list, "Chamfering", old_body_list, curve_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_edge_list );
  }

  // Do chamfering
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_chamfer( curve_list, left_offset, new_bodysm_list,
    right_offset, keep_old_body, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    CubitStatus stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    if( CubitUndo::get_undo_enabled() )
    {
      if( stat == CUBIT_SUCCESS )
        CubitUndo::note_result_bodies( new_body_list );
      else
        CubitUndo::remove_last_undo();
    }

    return stat;
  }
  return CUBIT_SUCCESS;
}

void GeometryModifyTool::propagate_from_small_edge(RefEdge *edge,
                                                 DLIList<RefEdge*> &small_edges,
                                                 DLIList<RefFace*> &narrow_faces,
                                                 DLIList<RefFace*> &processed_faces,
                                                 double small_edge_length)
{
  int i, j;
  // First find any small edges connected to this edge.
  DLIList<RefVertex*> edge_verts;
  edge->ref_vertices(edge_verts);
  for(i=edge_verts.size(); i--;)
  {
    RefVertex *vert = edge_verts.get_and_step();
    DLIList<RefEdge*> vert_edges;
    vert->ref_edges(vert_edges);
    for(j=vert_edges.size(); j--;)
    {
      RefEdge *cur_edge = vert_edges.get_and_step();
      if(cur_edge != edge && !cur_edge->marked())
      {
        // Mark of > 0 means it has been processed.
        cur_edge->marked(1);
        if(cur_edge->get_arc_length() < small_edge_length)
        {
          small_edges.append(cur_edge);
          // Mark of 2 means it is a small edge.
          cur_edge->marked(2);
          propagate_from_small_edge(cur_edge, small_edges, narrow_faces,
            processed_faces, small_edge_length);
        }
      }
    }
  }
  // Now look at adjacent narrow faces and recursively process them.
  DLIList<RefFace*> edge_faces;
  edge->ref_faces(edge_faces);
  for(i=edge_faces.size(); i--;)
  {
    RefFace *cur_face = edge_faces.get_and_step();
    if(!cur_face->marked())
    {
      cur_face->marked(1);
      if(GeomMeasureTool::narrow_region_exists(cur_face, small_edge_length))
      {
        DLIList<CubitVector> split_pos1_list;
        DLIList<CubitVector> split_pos2_list;
        GeomMeasureTool::find_split_points_for_narrow_regions(cur_face,
          small_edge_length, split_pos1_list, split_pos2_list);
        if(split_pos1_list.size() == 0)
        {
          narrow_faces.append_unique(cur_face);
          propagate_over_narrow_face(cur_face, edge, processed_faces,
            small_edges, narrow_faces, small_edge_length);
        }
      }
    }
  }
}

void GeometryModifyTool::propagate_over_narrow_face(RefFace *narrow_face,
                                                  RefEdge *edge,
                                                  DLIList<RefFace*> &processed_faces,
                                                  DLIList<RefEdge*> &small_edges,
                                                  DLIList<RefFace*> &narrow_faces,
                                                  double small_edge_length)
{
  int i, j;
  processed_faces.append(narrow_face);
  DLIList<RefEdge*> face_edges;
  narrow_face->ref_edges(face_edges);
  for(i=face_edges.size(); i--;)
  {
    RefEdge *cur_edge = face_edges.get_and_step();
    if(cur_edge != edge && !cur_edge->marked())
    {
      cur_edge->marked(1);
      if(cur_edge->get_arc_length() < small_edge_length)
      {
        cur_edge->marked(2);
        small_edges.append(cur_edge);
        propagate_from_small_edge(cur_edge, small_edges,
                narrow_faces, processed_faces, small_edge_length);
      }
      else
      {
        DLIList<RefFace*> edge_faces;
        cur_edge->ref_faces(edge_faces);
        for(j=edge_faces.size(); j--;)
        {
          RefFace *cur_face = edge_faces.get_and_step();
          if(cur_face != narrow_face)
          {
            if(!cur_face->marked())
            {
              cur_face->marked(1);
              if(GeomMeasureTool::narrow_region_exists(cur_face, small_edge_length))
              {
                DLIList<CubitVector> split_pos1_list;
                DLIList<CubitVector> split_pos2_list;
                GeomMeasureTool::find_split_points_for_narrow_regions(cur_face,
                  small_edge_length, split_pos1_list, split_pos2_list);
                if(split_pos1_list.size() == 0)
                {
                  narrow_faces.append_unique(cur_face);
                  propagate_over_narrow_face(cur_face, edge, processed_faces,
                    small_edges, narrow_faces, small_edge_length);
                }
              }
            }
          }
        }
      }
    }
  }
}

CubitStatus GeometryModifyTool::remove_topology( DLIList<RefEdge*> &ref_edge_list,
                                                  DLIList<RefFace*> &ref_face_list,
                                                  double backoff_distance,
                                                  double small_curve_size,
                                                  DLIList<Body*> &new_body_list,
                                                  CubitBoolean propagate,
                                                  CubitBoolean preview)
{
  int i, j;
  CubitStatus ret = CUBIT_SUCCESS;
  DLIList<Surface*> surf_list;
  DLIList<Curve*> curve_list;
  DLIList<Body*> old_body_list;
  GeometryModifyEngine *gme_ptr1=NULL, *gme_ptr2=NULL, *gme_ptr=NULL;
  Body *b1=NULL, *b2=NULL, *b=NULL;

  if(ref_edge_list.size())
    b1 = ref_edge_list.get()->body();
  if(ref_face_list.size())
    b2 = ref_face_list.get()->body();

  if(b1 && b2)
  {
    if(b1 == b2)
      b = b1;
  }
  else if(b1)
    b = b1;
  else if(b2)
    b = b2;

  if(b)
    old_body_list.append(b);
  else
  {
    PRINT_ERROR("Failed to find an owning body for the topology being removed.\n");
    ret = CUBIT_FAILURE;
  }

  if(ret == CUBIT_SUCCESS)
  {
    if (!okay_to_modify( old_body_list, "REMOVE_TOPOLOGY" ))
      ret = CUBIT_FAILURE;
    else
    {
      // Remove any edges from the list that aren't small enough.
      DLIList<RefEdge*> edges_to_remove;
      for(i=ref_edge_list.size(); i--;)
      {
        RefEdge* re = ref_edge_list.get_and_step();
        if(re->get_arc_length() > small_curve_size)
        {
          edges_to_remove.append(re);
          PRINT_INFO("Ignoring curve %d as it is not a small curve based on the input. "
            "Try a larger small_curve_size value.\n", re->id());
        }
      }
      ref_edge_list -= edges_to_remove;

      // Remove any faces from the list that don't have at least one small edge.
      DLIList<RefFace*> faces_to_remove;
      for(i=ref_face_list.size(); i--;)
      {
        DLIList<RefEdge*> face_edges;
        RefFace* rf = ref_face_list.get_and_step();
        rf->ref_edges(face_edges);
        int face_ok = 0;
        for(j=face_edges.size(); j && !face_ok; j--)
        {
          RefEdge* cur_edge = face_edges.get_and_step();
          if(cur_edge->get_arc_length() <= small_curve_size)
            face_ok = 1;
        }
        if(!face_ok)
        {
          faces_to_remove.append(rf);
          PRINT_INFO("Ignoring surface %d as it does not have at least one small curve in it based on the input. "
            "Try a larger small_curve_size value.\n", rf->id());
        }
      }
      ref_face_list -= faces_to_remove;
    }

    if(ref_face_list.size() > 0 || ref_edge_list.size() > 0)
    {
      // If told to do so propagate the topology to be removed to include
      // other narrow surfaces and small edges.
      if(propagate)
      {
        // Get all of the small edges into a single list.
        DLIList<RefEdge*> small_edges = ref_edge_list;
        for(i=ref_face_list.size(); i--;)
        {
          RefFace *face = ref_face_list.get_and_step();
          DLIList<RefEdge*> edges;
          face->ref_edges(edges);
          for(j=edges.size(); j--;)
          {
            RefEdge *edge = edges.get_and_step();
            if(edge->get_arc_length() < small_curve_size)
              small_edges.append(edge);
          }
        }
        small_edges.uniquify_ordered();

        DLIList<RefFace*> processed_faces;
        DLIList<RefEdge*> copy_of_small_edges = small_edges;
        DLIList<RefFace*> narrow_faces;

        DLIList<RefFace*> all_faces;
        DLIList<RefEdge*> all_edges;

        // Set all of the marked flags to 0.
        b->ref_faces(all_faces);
        for(i=all_faces.size(); i>0; i--)
          all_faces.get_and_step()->marked(0);

        b->ref_edges(all_edges);
        for(i=all_edges.size(); i>0; i--)
          all_edges.get_and_step()->marked(0);

        // Mark of >0 means it has been processed.
        // Mark of 2 means it is a small edge.
        for(i=small_edges.size(); i>0; i--)
          small_edges.get_and_step()->marked(2);

        // First look at all of the edges connected to small edges
        // to see if there are other small edges.
        while(copy_of_small_edges.size())
        {
          RefEdge *edge = copy_of_small_edges.extract();
          propagate_from_small_edge(edge, small_edges,
                  narrow_faces, processed_faces, small_curve_size);
        }

        ref_face_list += narrow_faces;
        ref_face_list.uniquify_ordered();

        ref_edge_list = small_edges;
        ref_edge_list.uniquify_ordered();
        // Append to face list here so we don't lose the ones that were passed in.
      }
    }
    else
    {
      PRINT_WARNING("No entities to remove.\n");
      ret = CUBIT_FAILURE;
    }
  }

  if(ret == CUBIT_SUCCESS)
  {
    if(ref_edge_list.size())
    {
      for(i=ref_edge_list.size(); i--;)
      {
        RefEdge *re = ref_edge_list.get_and_step();
        Curve *cur = re->get_curve_ptr();
        DLIList<TopologyBridge*> tmp_curve_list;
        GeometryQueryEngine *gqe = cur->get_geometry_query_engine();
        gqe->get_underlying_curves(cur, tmp_curve_list);
        if(tmp_curve_list.size() == 0)
          tmp_curve_list.append(cur);
        for(int p=tmp_curve_list.size(); p--;)
        {
          Curve *crv = dynamic_cast<Curve*>(tmp_curve_list.get_and_step());
          if(!gme_ptr1)
            gme_ptr1 = get_engine(crv);
          curve_list.append(crv);
        }
      }
    }
    if(ref_face_list.size())
    {
      for(i=ref_face_list.size(); i--;)
      {
        RefFace *rf = ref_face_list.get_and_step();
        Surface *sur = rf->get_surface_ptr();
        DLIList<TopologyBridge*> tmp_surf_list;
        GeometryQueryEngine *gqe = sur->get_geometry_query_engine();
        gqe->get_underlying_surfaces(sur, tmp_surf_list);
        if(tmp_surf_list.size() == 0)
          tmp_surf_list.append(sur);
        for(int p=tmp_surf_list.size(); p--;)
        {
          Surface *srf = dynamic_cast<Surface*>(tmp_surf_list.get_and_step());
          if(!gme_ptr2)
            gme_ptr2 = get_engine(srf);
          surf_list.append(srf);
        }
      }
    }
    if(gme_ptr1 && gme_ptr2)
    {
      if(gme_ptr1 == gme_ptr2)
        gme_ptr = gme_ptr1;
    }
    else if(gme_ptr1)
      gme_ptr = gme_ptr1;
    else if(gme_ptr2)
      gme_ptr = gme_ptr2;

    if(!gme_ptr)
    {
      PRINT_ERROR("Failed to find a geometry modify engine.\n");
      ret = CUBIT_FAILURE;
    }
  }

  if(ret == CUBIT_SUCCESS)
  {
    if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
    {
      DLIList<Body*> bodies;
      int i;
      for( i=ref_face_list.size(); i--; )
      {
        RefFace* ref_face = ref_face_list.get_and_step();
        bodies.append( ref_face->body() );
      }
      for( i=ref_edge_list.size(); i--; )
      {
        RefEdge* ref_edge = ref_edge_list.get_and_step();
        bodies.append( ref_edge->body() );
      }
      bodies.uniquify_unordered();
      CubitUndo::save_state_with_cubit_file( bodies );
    }

    if(preview == CUBIT_FALSE)
      do_attribute_setup();

    DLIList<BodySM*> body_sm_list(old_body_list.size());
    GeometryModifyEngine* gme = common_modify_engine(old_body_list, body_sm_list);

    if(preview == CUBIT_FALSE)
      push_vg_attributes_before_modify(body_sm_list);

    DLIList<BodySM*> new_bodysm_list;
    if(gme_ptr->remove_topology(curve_list, surf_list, backoff_distance, small_curve_size,
        new_bodysm_list, preview) == CUBIT_FAILURE)
    {
      if( CubitUndo::get_undo_enabled() )
        CubitUndo::remove_last_undo();

      if(!preview)
        remove_pushed_attributes(new_bodysm_list, old_body_list);

      ret = CUBIT_FAILURE;
    }
    else
    {
      if( preview == CUBIT_FALSE )
      {
        restore_vg_after_modify(new_bodysm_list, old_body_list, gme);
        remove_pushed_attributes(new_bodysm_list, old_body_list);

        // Update DAG
        ret = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

        if( CubitUndo::get_undo_enabled() )
        {
          if( ret == CUBIT_FAILURE)
            CubitUndo::remove_last_undo();
          else
            CubitUndo::note_result_bodies( new_body_list );
        }
      }
    }

    if(preview == CUBIT_FALSE)
      do_attribute_cleanup();
  }

  return ret;
}

//=============================================================================
// Description: Chamfer vertices on solid or sheet bodies.  On a solid body
//              there can be up to 3 offsets; on a sheet body up to 2 offsets.
//              The offsets are in the direction of the supplied edges.  If
//              multiple vertices are supplied, only one offset value is
//              allowed and the edges are not used.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus
GeometryModifyTool::tweak_chamfer( DLIList<RefVertex*> &ref_vertex_list,
                                   double offset1,
                                   DLIList<Body*> &new_body_list,
                                   RefEdge *edge1,
                                   double offset2,
                                   RefEdge *edge2,
                                   double offset3,
                                   RefEdge *edge3,
                                   CubitBoolean keep_old_body,
                                   CubitBoolean preview )
{

  if( offset1 <= 0.0 )
  {
    PRINT_ERROR( "Chamfer radius not specified.\n" );
    return CUBIT_FAILURE;
  }

  DLIList<Point*> point_list(ref_vertex_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_vertex_list, "Chamfering", old_body_list, point_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  RefVertex *ref_vertex_ptr = ref_vertex_list.get();
  if( ref_vertex_list.size() == 1 && (edge1 || edge2 || edge3) )
  {
    // Make sure input edges are attached to the input vertex
    DLIList<RefEdge*> ref_edge_list;
    ref_vertex_ptr->ref_edges( ref_edge_list );

    if( edge1 && !ref_edge_list.is_in_list( edge1 ) )
    {
      PRINT_ERROR( "Specified curve %d not attached to vertex %d\n", edge1->id(),
        ref_vertex_ptr->id() );
      return CUBIT_FAILURE;
    }

    if( edge2 && !ref_edge_list.is_in_list( edge2 ) )
    {
      PRINT_ERROR( "Specified curve %d not attached to vertex %d\n", edge1->id(),
        ref_vertex_ptr->id() );
      return CUBIT_FAILURE;
    }

    if( edge3 && !ref_edge_list.is_in_list( edge3 ) )
    {
      PRINT_ERROR( "Specified curve %d not attached to vertex %d\n", edge1->id(),
        ref_vertex_ptr->id() );
      return CUBIT_FAILURE;
    }

    // Make sure offsets supplied
    if( edge1 && offset1 < 0.0 )
    {
      PRINT_ERROR( "Offset for curve %d specified incorrectly.\n", edge1->id() );
      return CUBIT_FAILURE;
    }

    if( edge2 && offset2 < 0.0 )
    {
      PRINT_ERROR( "Offset for curve %d specified incorrectly.\n", edge2->id() );
      return CUBIT_FAILURE;
    }

    if( edge3 && offset3 < 0.0 )
    {
      PRINT_ERROR( "Offset for curve %d specified incorrectly.\n", edge3->id() );
      return CUBIT_FAILURE;
    }
  }

  if( point_list.size() > 1 && offset2 > 0.0 )
  {
    PRINT_ERROR( "Cannot supply multiple radii when chamfering multiple vertices.\n" );
    return CUBIT_FAILURE;
  }

  Curve *curve1 = NULL;
  Curve *curve2 = NULL;
  Curve *curve3 = NULL;

  if( edge1 )
  {
    TopologyBridge* bridge = 0;
    GeometryModifyEngine* edge_gme = get_engine(edge1, &bridge);
    if( gme_ptr != edge_gme )
    {
      PRINT_ERROR( "Specified curve %d must belong to same geometry engine as vertex %d.\n",
        edge1->id(), ref_vertex_ptr->id()  );
      return CUBIT_FAILURE;
    }
    curve1 = dynamic_cast<Curve*>(bridge);
  }

  if( edge2 )
  {
    TopologyBridge* bridge = 0;
    GeometryModifyEngine* edge_gme = get_engine(edge2, &bridge);
    if( gme_ptr != edge_gme )
    {
      PRINT_ERROR( "Specified curve %d must belong to same geometry engine as vertex %d.\n",
        edge2->id(), ref_vertex_ptr->id()  );
      return CUBIT_FAILURE;
    }
    curve2 = dynamic_cast<Curve*>(bridge);
  }

  if( edge3 )
  {
    TopologyBridge* bridge = 0;
    GeometryModifyEngine* edge_gme = get_engine(edge3, &bridge);
    if( gme_ptr != edge_gme )
    {
      PRINT_ERROR( "Specified curve %d must belong to same geometry engine as vertex %d.\n",
        edge3->id(), ref_vertex_ptr->id()  );
      return CUBIT_FAILURE;
    }
    curve3 = dynamic_cast<Curve*>(bridge);
  }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_vertex_list );
  }

  // Do chamfering
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_chamfer( point_list, offset1, new_bodysm_list, curve1, offset2,
    curve2, offset3, curve3, keep_old_body, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    CubitStatus stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    if( CubitUndo::get_undo_enabled()  )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old_body )
        CubitUndo::note_result_bodies( new_body_list );
    }
    return stat;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Create a round fillet (or blend) at the given curves on solid
//              or sheet bodies.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_fillet( DLIList<RefEdge*> &ref_edge_list,
                                              double radius,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_edge_list, "Filleting", old_body_list, curve_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_edge_list );
  }

  // Do filleting
  DLIList<BodySM*> new_bodysm_list;
  if( gme_ptr->tweak_fillet(curve_list, radius, new_bodysm_list, keep_old_body,
    preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    CubitStatus stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    if( CubitUndo::get_undo_enabled() )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old_body )
        CubitUndo::note_result_bodies( new_body_list );
    }

    return stat;
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Create a round fillet (or blend) at the given curves on a solid
//              or sheet body.  The fillet has a variable radius from the start
//              to the end of the curve.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_fillet( RefEdge *ref_edge_ptr,
                                              double start_radius,
                                              double end_radius,
                                              Body *&new_body_ptr,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<Curve*> curve_list(1);
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  DLIList<RefEdge*> ref_edge_list(1);
  ref_edge_list.append( ref_edge_ptr );

  gme_ptr = tweak_setup( ref_edge_list, "Filleting", old_body_list, curve_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
    {
      DLIList<RefEdge*> edges(1);
      edges.append( ref_edge_ptr );
      CubitUndo::save_state_with_cubit_file( edges );
    }
  }

  // Do filleting
  BodySM *new_bodysm_ptr;
  Curve *curve_ptr = curve_list.get();
  CubitStatus stat = gme_ptr->tweak_fillet( curve_ptr, start_radius,
                                            end_radius, new_bodysm_ptr,
                                            keep_old_body, preview );

  if( CubitUndo::get_undo_enabled() )
    if( stat == CUBIT_FAILURE)
      CubitUndo::remove_last_undo();

  if( stat == CUBIT_FAILURE )
    return stat;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    DLIList<BodySM*> new_bodysm_list;
    new_bodysm_list.append( new_bodysm_ptr );
    DLIList<Body*> new_body_list;
    stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    if( CubitUndo::get_undo_enabled()  )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old_body )
        CubitUndo::note_result_bodies( new_body_list );
    }

    new_body_ptr = new_body_list.get();

    return stat;
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Create a round fillet (or blend) at the given vertices on sheet
//              bodies.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus
GeometryModifyTool::tweak_fillet( DLIList<RefVertex*> &ref_vertex_list,
                                  double radius,
                                  DLIList<Body*> &new_body_list,
                                  CubitBoolean keep_old_body,
                                  CubitBoolean preview )
{
  DLIList<Point*> point_list(ref_vertex_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_vertex_list, "Filleting", old_body_list, point_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_vertex_list );
  }

  // Do filleting
  DLIList<BodySM*> new_bodysm_list;
  CubitStatus stat = gme_ptr->tweak_fillet( point_list, radius,
                                            new_bodysm_list, keep_old_body,
                                            preview );

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
    if( stat == CUBIT_FAILURE)
      CubitUndo::remove_last_undo();

  if( stat == CUBIT_FAILURE )
    return stat;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    if( CubitUndo::get_undo_enabled()  )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old_body )
        CubitUndo::note_result_bodies( new_body_list );
    }

    return stat;
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified faces of a volume or volumes along a vector.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_move( DLIList<RefFace*>& ref_face_list,
                                            const CubitVector &delta,
                                            DLIList<Body*>& new_body_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview )
{
  DLIList<Surface*> surface_list(ref_face_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_face_list, "Moving", old_body_list, surface_list, CUBIT_TRUE );
  if (!gme_ptr)
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
    {
      int i;
      DLIList<RefEdge*> ref_edges;
      for( i=ref_face_list.size(); i--; )
        ref_face_list.get_and_step()->ref_edges( ref_edges );
      ref_edges.uniquify_unordered();
      CubitUndo::save_state_with_cubit_file( ref_edges );
    }
  }

  int i;
  DLIList<BodySM*> body_sms;
  for(i=old_body_list.size(); i--;)
  {
    BodySM* bsm = old_body_list.get_and_step()->get_body_sm_ptr();
    if(bsm)
      body_sms.append_unique(bsm);
  }

  if(!preview)
  {
    do_attribute_setup();
    push_vg_attributes_before_modify(body_sms);
  }

  // Do move
  DLIList<BodySM*> new_bodysm_list;
  CubitStatus stat = gme_ptr->tweak_move( surface_list, delta,
                                          new_bodysm_list, keep_old_body,
                                          preview );


  if( CubitUndo::get_undo_enabled() )
    if( stat == CUBIT_FAILURE)
      CubitUndo::remove_last_undo();

  if( stat == CUBIT_FAILURE )
  {
    if(!preview)
    {
      remove_pushed_attributes(new_bodysm_list, old_body_list);
      do_attribute_cleanup();
    }
    return stat;
  }
  else
  {
    if(!preview)
    {
      restore_vg_after_modify(new_bodysm_list, old_body_list, gme_ptr);
      remove_pushed_attributes(new_bodysm_list, old_body_list);
    }
  }

  // loop body sm list and find surfaces that need updating.
  // this is to account for some cases where the topology doesn't change, but the geometry does.
  DLIList<RefEntity*> entities_to_update;
  for(i=0; i<new_bodysm_list.size(); i++)
  {
    BodySM* bodysm = new_bodysm_list.get_and_step();
    DLIList<Surface*> surfs;
    bodysm->surfaces(surfs);
    int j;
    // find a surface that is also found in our input list
    for(j=0; j<surfs.size(); j++, surfs.step())
    {
      BridgeManager* man = surfs.get()->bridge_manager();
      if(man)
      {
        RefFace* ref_face = CAST_TO(man->topology_entity(), RefFace);
        if(ref_face && ref_face_list.is_in_list(ref_face))
        {
          // get neighbors
          DLIList<Point*> neighbor_points;
          surfs.get()->points(neighbor_points);
          DLIList<Curve*> neighbor_curves;
          DLIList<Surface*> neighbor_surfaces;
          DLIList<TopologyBridge*> neighbors;
          DLIList<TopologyBridge*> tmp;
          int k;
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->surfaces(neighbor_surfaces);
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->curves(neighbor_curves);

          CAST_LIST_TO_PARENT(neighbor_points, tmp);
          neighbors += tmp;
          neighbor_curves.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_curves, tmp);
          neighbors += tmp;
          neighbor_surfaces.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_surfaces, tmp);
          neighbors += tmp;
          neighbors.append(surfs.get()->lump());
          neighbors.append(surfs.get()->bodysm());

          for(k=0; k<neighbors.size(); k++)
            if(BridgeManager* m = neighbors.get_and_step()->bridge_manager())
              if(TopologyEntity* t = m->topology_entity())
                entities_to_update.append(CAST_TO(t, RefEntity));
        }
      }
    }
  }

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    if( CubitUndo::get_undo_enabled()  )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old_body )
        CubitUndo::note_result_bodies( new_body_list );
    }

    // Update graphics
    while (entities_to_update.size())
      entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );

    do_attribute_cleanup();
  }


  //collect all the new faces
  DLIList<RefFace*> new_faces;
  for( i=new_body_list.size(); i--; )
  {
    Body *new_body = new_body_list.get_and_step();
    DLIList<RefFace*> tmp_faces;
    new_body->ref_faces( tmp_faces );
    new_faces += tmp_faces;
  }

  //unmerge any merged adjacent surfaces or
  //merged curves in unmerged adjacent surfaces
  DLIList<RefFace*> adjacent_faces_to_unmerge;
  DLIList<RefEdge*> adjacent_edges_to_unmerge;
  for(i=ref_face_list.size(); i--;)
  {
    RefFace *tweaked_face = ref_face_list.get_and_step();
    if( !new_faces.move_to( tweaked_face ) )
      continue;

    //get all the edges of the face you tweaked
    DLIList<RefEdge*> tweaked_face_edges;
    tweaked_face->ref_edges( tweaked_face_edges );
    adjacent_edges_to_unmerge += tweaked_face_edges;

    //get all the adjacent faces to this edge
    int j;
    for( j=tweaked_face_edges.size(); j--; )
    {
      RefEdge *tmp_edge = tweaked_face_edges.get_and_step();
      DLIList<RefFace*> tmp_faces;
      tmp_edge->ref_faces( tmp_faces );
      tmp_faces.remove( tweaked_face );
      adjacent_faces_to_unmerge += tmp_faces;
    }

    //get all edges not in the surface,
    //sharing vertices with the surface
    DLIList<RefVertex*> ref_vertices;
    tweaked_face->ref_vertices( ref_vertices );
    for( j=ref_vertices.size(); j--; )
    {
      RefVertex *tmp_vert = ref_vertices.get_and_step();
      DLIList<RefEdge*>  ref_edges;
      tmp_vert->ref_edges( ref_edges );

      int k;
      for( k=ref_edges.size(); k--; )
      {
        RefEdge *tmp_edge = ref_edges.get_and_step();
        if( !tweaked_face_edges.move_to( tmp_edge ) )
          adjacent_edges_to_unmerge.append( tmp_edge );
      }
    }
  }

  //unmerge any adjacent faces
  adjacent_faces_to_unmerge.uniquify_unordered();
  for( i=adjacent_faces_to_unmerge.size(); i--; )
  {
    RefFace *ref_face = adjacent_faces_to_unmerge.get_and_step();

    DLIList<TopologyBridge*> bridge_list;
    ref_face->bridge_manager()->get_bridge_list(bridge_list);
    if (bridge_list.size() > 1)
    {
      if( MergeTool::instance()->unmerge( ref_face ) )
        PRINT_WARNING("Unmerging Surface %d\n", ref_face->id() );
    }
  }

  //unmerge any adjacent edges
  adjacent_edges_to_unmerge.uniquify_unordered();
  for( i=adjacent_edges_to_unmerge.size(); i--; )
  {
    RefEdge *ref_edge = adjacent_edges_to_unmerge.get_and_step();
    DLIList<TopologyBridge*> bridge_list;
    ref_edge->bridge_manager()->get_bridge_list(bridge_list);
    if (bridge_list.size() > 1)
    {
      if( MergeTool::instance()->unmerge( ref_edge) )
        PRINT_WARNING("Unmerging Curve %d\n", ref_edge->id() );
    }
  }


  return stat;
}

//=============================================================================
// Description: Tweak specified curves of a sheet body along a vector.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_move( DLIList<RefEdge*>& ref_edge_list,
                                            const CubitVector &delta,
                                            DLIList<Body*>& new_body_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_edge_list, "Moving", old_body_list, curve_list );
  if (!gme_ptr)
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_edge_list );
  }

  // Do move
  DLIList<BodySM*> new_bodysm_list;
  CubitStatus stat = gme_ptr->tweak_move( curve_list, delta,
                                          new_bodysm_list, keep_old_body,
                                          preview );

  if( CubitUndo::get_undo_enabled() )
    if( stat == CUBIT_FAILURE)
      CubitUndo::remove_last_undo();

  if( stat == CUBIT_FAILURE )
    return stat;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    if( CubitUndo::get_undo_enabled()  )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old_body )
        CubitUndo::note_result_bodies( new_body_list );
    }
  }

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return stat;
}

//=============================================================================
// Description: Tweak specified faces of a volume or volumes by offsetting
//              those faces by the offset distance.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_offset( DLIList<RefFace*> &ref_face_list,
                                              double offset_distance,
                                              DLIList<RefFace*> *add_ref_face_list_ptr,
                                              DLIList<double> *add_offset_list_ptr,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<RefFace*> all_ref_face_list(ref_face_list.size());
  all_ref_face_list = ref_face_list;
  if( add_ref_face_list_ptr->size() )
    all_ref_face_list += *add_ref_face_list_ptr;

  DLIList<Surface*> surface_list(ref_face_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_face_list, "Offsetting", old_body_list, surface_list, CUBIT_TRUE );
  if (!gme_ptr)
    return CUBIT_FAILURE;

  DLIList<Surface*> add_surface_list;
  if( add_ref_face_list_ptr && add_ref_face_list_ptr->size() )
  {
    DLIList<Body*> old_body_list2;
    GeometryModifyEngine* gme_ptr2 = tweak_setup( *add_ref_face_list_ptr, "Offsetting", 
      old_body_list2, add_surface_list, CUBIT_TRUE );
    if (!gme_ptr2)
      return CUBIT_FAILURE;
    if( gme_ptr != gme_ptr2 )
    {
      PRINT_ERROR("Offsetting surfaces on volumes containing surfaces from different\n"
        "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
    }
    old_body_list.merge_unique( old_body_list2 );
  }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
    {
      int i;
      DLIList<RefEdge*> ref_edges;
      for( i=ref_face_list.size(); i--; )
        all_ref_face_list.get_and_step()->ref_edges( ref_edges );
      ref_edges.uniquify_unordered();
      CubitUndo::save_state_with_cubit_file( ref_edges );
    }
  }

  int i;
  if(!preview)
  {
    DLIList<BodySM*> body_sms;
    for(i=old_body_list.size(); i--;)
    {
      BodySM* bsm = old_body_list.get_and_step()->get_body_sm_ptr();
      if(bsm)
        body_sms.append_unique(bsm);
    }

    do_attribute_setup();
    push_vg_attributes_before_modify(body_sms);
  }

  // Do offset
  DLIList<BodySM*> new_bodysm_list;
  CubitStatus stat;
  if( add_surface_list.size() )
    stat = gme_ptr->tweak_offset( surface_list, offset_distance, 
                                  &add_surface_list, add_offset_list_ptr,
                                  new_bodysm_list, keep_old_body, preview );
  else
    stat = gme_ptr->tweak_offset( surface_list, offset_distance, NULL, NULL,
                                  new_bodysm_list, keep_old_body, preview );

  if( CubitUndo::get_undo_enabled() )
    if( stat == CUBIT_FAILURE)
      CubitUndo::remove_last_undo();

  if( stat == CUBIT_FAILURE )
  {
    if(!preview)
    {
      remove_pushed_attributes(new_bodysm_list, old_body_list);
      do_attribute_cleanup();
    }
    return stat;
  }
  else
  {
    if(!preview)
    {
      restore_vg_after_modify(new_bodysm_list, old_body_list, gme_ptr);
      remove_pushed_attributes(new_bodysm_list, old_body_list);
    }
  }

  // Collect all of the old faces to be compared later with the new faces...DJQ
  DLIList<RefFace*> old_faces;
  for (i = 0; i < old_body_list.size(); i++)
  {
      Body *old_body = old_body_list.get_and_step();
      DLIList<RefFace*> tmp_faces;
      old_body->ref_faces(tmp_faces);
      old_faces +=tmp_faces;
  }

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    if( CubitUndo::get_undo_enabled()  )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old_body )
        CubitUndo::note_result_bodies( new_body_list );
    }
    do_attribute_cleanup();
  }

  //collect all the new faces
  DLIList<RefFace*> new_faces;
  for( i=new_body_list.size(); i--; )
  {
    Body *new_body = new_body_list.get_and_step();
    DLIList<RefFace*> tmp_faces;
    new_body->ref_faces( tmp_faces );
    new_faces += tmp_faces;
  }

  // Compare the new_faces list with the old_faces list to determine which faces are created
  // Add these faces to the all_ref_face_list to check for its neighbors...DJQ
  DLIList<RefFace*> difference = new_faces;
  difference -= old_faces;
  all_ref_face_list += difference;

  // loop body sm list and find surfaces that need updating.
  // this is to account for some cases where the topology
  //doesn't change, but the geometry does.

  DLIList<RefEntity*> entities_to_update;
  for(i=0; i<new_bodysm_list.size(); i++)
  {
    BodySM* bodysm = new_bodysm_list.get_and_step();
    DLIList<Surface*> surfs;
    bodysm->surfaces(surfs);
    int j;
    // find a surface that is also found in our input list
    for(j=0; j<surfs.size(); j++, surfs.step())
    {
      BridgeManager* man = surfs.get()->bridge_manager();
      if(man)
      {
        RefFace* ref_face = CAST_TO(man->topology_entity(), RefFace);
        if(ref_face && all_ref_face_list.is_in_list(ref_face))
        {
          // get neighbors
          DLIList<Point*> neighbor_points;
          surfs.get()->points(neighbor_points);
          DLIList<Curve*> neighbor_curves;
          DLIList<Surface*> neighbor_surfaces;
          DLIList<TopologyBridge*> neighbors;
          DLIList<TopologyBridge*> tmp;
          int k;
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->surfaces(neighbor_surfaces);
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->curves(neighbor_curves);

          CAST_LIST_TO_PARENT(neighbor_points, tmp);
          neighbors += tmp;
          neighbor_curves.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_curves, tmp);
          neighbors += tmp;
          neighbor_surfaces.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_surfaces, tmp);
          neighbors += tmp;
          neighbors.append(surfs.get()->lump());
          neighbors.append(surfs.get()->bodysm());

          for(k=0; k<neighbors.size(); k++)
            if(BridgeManager* m = neighbors.get_and_step()->bridge_manager())
              if(TopologyEntity* t = m->topology_entity())
                entities_to_update.append(CAST_TO(t, RefEntity));
        }
      }
    }
  }


  //unmerge any merged adjacent surfaces or
  //merged curves in unmerged adjacent surfaces
  DLIList<RefFace*> adjacent_faces_to_unmerge;
  DLIList<RefEdge*> adjacent_edges_to_unmerge;
  for(i=ref_face_list.size(); i--;)
  {
    RefFace *tweaked_face = ref_face_list.get_and_step();
    if( !new_faces.move_to( tweaked_face ) )
      continue;

    //get all the edges of the face you tweaked
    DLIList<RefEdge*> tweaked_face_edges;
    tweaked_face->ref_edges( tweaked_face_edges );
    adjacent_edges_to_unmerge += tweaked_face_edges;

    //get all the adjacent faces to this edge
    int j;
    for( j=tweaked_face_edges.size(); j--; )
    {
      RefEdge *tmp_edge = tweaked_face_edges.get_and_step();
      DLIList<RefFace*> tmp_faces;
      tmp_edge->ref_faces( tmp_faces );
      tmp_faces.remove( tweaked_face );
      adjacent_faces_to_unmerge += tmp_faces;
    }

    //get all edges not in the surface,
    //sharing vertices with the surface
    DLIList<RefVertex*> ref_vertices;
    tweaked_face->ref_vertices( ref_vertices );
    for( j=ref_vertices.size(); j--; )
    {
      RefVertex *tmp_vert = ref_vertices.get_and_step();
      DLIList<RefEdge*>  ref_edges;
      tmp_vert->ref_edges( ref_edges );

      int k;
      for( k=ref_edges.size(); k--; )
      {
        RefEdge *tmp_edge = ref_edges.get_and_step();
        if( !tweaked_face_edges.move_to( tmp_edge ) )
          adjacent_edges_to_unmerge.append( tmp_edge );
      }
    }
  }

  //unmerge any adjacent faces
  adjacent_faces_to_unmerge.uniquify_unordered();
  for( i=adjacent_faces_to_unmerge.size(); i--; )
  {
    RefFace *ref_face = adjacent_faces_to_unmerge.get_and_step();

    DLIList<TopologyBridge*> bridge_list;
    ref_face->bridge_manager()->get_bridge_list(bridge_list);
    if (bridge_list.size() > 1)
    {
      if( MergeTool::instance()->unmerge( ref_face ) )
        PRINT_WARNING("Unmerging Surface %d\n", ref_face->id() );
    }
  }

  //unmerge any adjacent edges
  adjacent_edges_to_unmerge.uniquify_unordered();
  for( i=adjacent_edges_to_unmerge.size(); i--; )
  {
    RefEdge *ref_edge = adjacent_edges_to_unmerge.get_and_step();
    DLIList<TopologyBridge*> bridge_list;
    ref_edge->bridge_manager()->get_bridge_list(bridge_list);
    if (bridge_list.size() > 1)
    {
      if( MergeTool::instance()->unmerge( ref_edge) )
        PRINT_WARNING("Unmerging Curve %d\n", ref_edge->id() );
    }
  }

  // Update graphics
  while (entities_to_update.size())
    entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );

  return stat;
}

//=============================================================================
// Description: Tweak specified curves of a sheet body or bodies by offsetting
//              those curves by the offset distance.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_offset( DLIList<RefEdge*> &ref_edge_list,
                                              double offset_distance,
                                              DLIList<RefEdge*> *add_ref_edge_list_ptr,
                                              DLIList<double> *add_offset_list_ptr,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<RefEdge*> all_ref_edge_list(ref_edge_list.size());
  all_ref_edge_list = ref_edge_list;
  if( add_ref_edge_list_ptr )
    all_ref_edge_list += *add_ref_edge_list_ptr;

  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_edge_list, "Offsetting", old_body_list, curve_list );
  if (!gme_ptr)
    return CUBIT_FAILURE;

  DLIList<Curve*> add_curve_list;
  if( add_ref_edge_list_ptr && add_ref_edge_list_ptr->size() )
  {
    DLIList<Body*> old_body_list2;
    GeometryModifyEngine* gme_ptr2 = tweak_setup( *add_ref_edge_list_ptr, "Offsetting", 
      old_body_list2, add_curve_list );
    if (!gme_ptr2)
      return CUBIT_FAILURE;
    if( gme_ptr != gme_ptr2 )
    {
      PRINT_ERROR("Offsetting curves on entities containing surfaces from different\n"
      "       geometry engines is not allowed.\n");
      return CUBIT_FAILURE;
    }
    old_body_list.merge_unique( old_body_list2 );
  }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( all_ref_edge_list );
  }

  // Do offset
  DLIList<BodySM*> new_bodysm_list;
  CubitStatus stat;
  if( add_curve_list.size() )
    stat = gme_ptr->tweak_offset( curve_list, offset_distance, &add_curve_list,
    add_offset_list_ptr, new_bodysm_list, keep_old_body, preview );
  else
    stat = gme_ptr->tweak_offset( curve_list, offset_distance, NULL, NULL,
    new_bodysm_list, keep_old_body, preview );

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
    if( stat == CUBIT_FAILURE)
      CubitUndo::remove_last_undo();

  if( stat == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Update DAG
  if( preview == CUBIT_FALSE )
  {
    stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    if( CubitUndo::get_undo_enabled()  )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old_body )
        CubitUndo::note_result_bodies( new_body_list );
    }
  }

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return stat;
}


CubitStatus GeometryModifyTool::tweak_remove_individually(
                                              DLIList<RefFace*> &ref_face_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean keep_surface,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_face_list );
  }

  // Split things up if individual
  //build a surfae to volume map
  std::map<RefFace*, RefVolume*> surface_to_volume_map;

  int i;
  for (i = ref_face_list.size(); i--; )
  {
    RefFace *tmp_face = ref_face_list.get_and_step();
    RefVolume *tmp_vol = tmp_face->ref_volume();
    surface_to_volume_map.insert( std::map<RefFace*, RefVolume*>::value_type( tmp_face, tmp_vol));
  }

  DLIList<RefFace*> one_ref_face;
  DLIList<Body*> tmp_new_body_list;
  CubitStatus total_rv = CUBIT_FAILURE;
  bool extend = true;

  // Succeed if any one surface succeeds.
  for (i = ref_face_list.size(); i--; )
  {
    //make sure that the surface to remove is still in the body...
    //that it hasn't been removed from a previous tweak operation
    RefFace *tmp_face = ref_face_list.get_and_step();
    std::map<RefFace*, RefVolume*>::iterator tmp_iter;
    tmp_iter = surface_to_volume_map.find( tmp_face );
    RefVolume *tmp_vol = tmp_iter->second;
    DLIList<RefFace*> ref_face_list;
    tmp_vol->ref_faces( ref_face_list );
    if( !ref_face_list.move_to( tmp_face ) )
      continue;

    one_ref_face.clean_out();
    one_ref_face.append( tmp_face );
    int id = one_ref_face.get()->id();

    //See if the owning body of the face is a multi-volume body
    Body *owning_body = one_ref_face.get()->body();
    int number_volumes_before = owning_body->num_ref_volumes();

    tmp_new_body_list.clean_out();

    CubitStatus rv = this->tweak_remove(one_ref_face, tmp_new_body_list,
      extend, keep_surface, keep_old_body, preview );
    if (rv)
    {
      total_rv = CUBIT_SUCCESS;
      if( !preview )
        PRINT_INFO("Successfully removed Surface %d\n\n", id);
      else
        PRINT_INFO("Successfully removed Surface %d in preview\n\n", id);

      //see if we have a multi-volume body or multiple bodies
      //if so, we know the original volume was destroyed, so we
      //cannot remove any more surfaces because the check above is
      //highly likely to crash
      bool volume_destroyed = false;
      int number_volumes_after = tmp_new_body_list.get()->num_ref_volumes();
      if( number_volumes_after > number_volumes_before  )
        volume_destroyed = true;

      new_body_list += tmp_new_body_list;

      if( volume_destroyed == true )
      {
        PRINT_WARNING("Unable to remove more surfaces because multiple bodies\n"
                      "       have been produced from removing surfaces individually\n" );
        if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
        {
          if( total_rv == CUBIT_FAILURE )  //didn't remove any surfaces
            CubitUndo::remove_last_undo();
          else
            CubitUndo::note_result_bodies( new_body_list );
        }

        return total_rv;
      }
    }
    else
    {
      if( !preview )
        PRINT_INFO("Unable to remove Surface %d\n\n", id);
      else
        PRINT_INFO("Unable to remove Surface %d in preview\n\n", id);
    }
  }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( total_rv == CUBIT_FAILURE)
      CubitUndo::remove_last_undo();
    else if( keep_old_body || keep_surface )
      CubitUndo::note_result_bodies( new_body_list );
  }

  return total_rv;
}


//=============================================================================
// Description: Function to remove surfaces from a body and then extend the
//              remaining surfaces to fill the gap or hole.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_remove_together(
                                              DLIList<RefFace*> &ref_face_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean extend_adjoining,
                                              CubitBoolean keep_surface,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_face_list );
  }

  CubitStatus stat = tweak_remove( ref_face_list, new_body_list,
                                   extend_adjoining, keep_surface,
                                   keep_old_body, preview );

  if( stat == CUBIT_FAILURE )
  {
    if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
      CubitUndo::remove_last_undo();

    return CUBIT_FAILURE;
  }


  if( preview == CUBIT_FALSE )
  {
    if( CubitUndo::get_undo_enabled()  )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old_body || keep_surface )
        CubitUndo::note_result_bodies( new_body_list );
    }
  }

  return CUBIT_SUCCESS;
}

//private funcion...should not be called from outside this class
CubitStatus GeometryModifyTool::tweak_remove( DLIList<RefFace*> &ref_face_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean extend_adjoining,
                                              CubitBoolean keep_surface,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
   DLIList<Surface*> surface_list(ref_face_list.size());
   DLIList<Body*> old_body_list;
   GeometryModifyEngine* gme_ptr;

   // clear any preview previews
   GfxPreview::clear();

   //collect all neighboring surfaces to those in the list
   int i,j;
   DLIList<RefFace*> neighboring_surfaces;
   for( i=ref_face_list.size(); i--; )
   {
     RefFace *tmp_face = ref_face_list.get_and_step();
     DLIList<RefEdge*> ref_edge_list;
     tmp_face->ref_edges( ref_edge_list );
     for( j=ref_edge_list.size(); j--; )
       ref_edge_list.get_and_step()->ref_faces( neighboring_surfaces );
   }

  //uniquify and add other surfaces
  neighboring_surfaces.uniquify_unordered();
  neighboring_surfaces += ref_face_list;

   gme_ptr = tweak_setup( ref_face_list, "Removing", old_body_list, surface_list, CUBIT_TRUE );
   if (!gme_ptr)
     return CUBIT_FAILURE;

   DLIList<Surface*> kept_surface_list;
   if( keep_surface )
   {
     int kk;
     for( kk=surface_list.size(); kk--; )
     {
       Surface *new_surf = gme_ptr->make_Surface( surface_list.get_and_step() );
       kept_surface_list.append( new_surf );
     }
   }

   DLIList<BodySM*> body_sms;
   for(i=old_body_list.size(); i--;)
   {
     BodySM* bsm = old_body_list.get_and_step()->get_body_sm_ptr();
     if(bsm)
       body_sms.append_unique(bsm);
   }

   if(!preview)
   {
    do_attribute_setup();
    push_vg_attributes_before_modify(body_sms);
   }

   // Do remove
   DLIList<BodySM*> new_bodysm_list;
   CubitStatus removal_status =
       gme_ptr->tweak_remove( surface_list,
                              new_bodysm_list,
                              extend_adjoining,
                              keep_old_body,
                              preview );


   if( removal_status == CUBIT_FAILURE )
   {
     if(!preview)
      remove_pushed_attributes(new_bodysm_list, old_body_list);
     if( keep_surface )
     {
       int kk;
       for( kk=kept_surface_list.size(); kk--; )
       {
         Surface *surf = kept_surface_list.get_and_step();
         gme_ptr->get_gqe()->delete_solid_model_entities( surf );
       }
     }
     if(!preview)
      do_attribute_cleanup();
     return CUBIT_FAILURE;
   }
   else
   {
     if(!preview)
     {
        restore_vg_after_modify(new_bodysm_list, old_body_list, gme_ptr);
        remove_pushed_attributes(new_bodysm_list, old_body_list);
     }
   }

  // loop body sm list and find surfaces that need updating.
  // this is to account for some cases where the topology doesn't change, but the geometry does.
  DLIList<RefEntity*> entities_to_update;
  for(i=0; i<new_bodysm_list.size(); i++)
  {
    BodySM* bodysm = new_bodysm_list.get_and_step();
    DLIList<Surface*> surfs;
    bodysm->surfaces(surfs);
    int j;
    // find a surface that is also found in our input list
    for(j=0; j<surfs.size(); j++, surfs.step())
    {
      BridgeManager* man = surfs.get()->bridge_manager();
      if(man)
      {
        RefFace* ref_face = CAST_TO(man->topology_entity(), RefFace);
        if( ref_face && neighboring_surfaces.is_in_list(ref_face) )
        {
          // get neighbors
          DLIList<Point*> neighbor_points;
          surfs.get()->points(neighbor_points);
          DLIList<Curve*> neighbor_curves;
          DLIList<Surface*> neighbor_surfaces;
          DLIList<TopologyBridge*> neighbors;
          DLIList<TopologyBridge*> tmp;
          int k;
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->surfaces(neighbor_surfaces);
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->curves(neighbor_curves);

          CAST_LIST_TO_PARENT(neighbor_points, tmp);
          neighbors += tmp;
          neighbor_curves.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_curves, tmp);
          neighbors += tmp;
          neighbor_surfaces.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_surfaces, tmp);
          neighbors += tmp;
          neighbors.append(surfs.get()->lump());
          neighbors.append(surfs.get()->bodysm());

          for(k=0; k<neighbors.size(); k++)
          {
            if(BridgeManager* m = neighbors.get_and_step()->bridge_manager())
            {
              if(TopologyEntity* t = m->topology_entity())
              {
                entities_to_update.append(CAST_TO(t, RefEntity));
                //RefEntity *ref_ent = CAST_TO(t, RefEntity );
              }
            }
          }
        }
      }
    }
  }


  DLIList<Body*> kept_surface_bodies;
  if( preview == CUBIT_FALSE && keep_surface )
  {
    int kk;
    for( kk=kept_surface_list.size(); kk--; )
    {
      Surface *surf = kept_surface_list.get_and_step();
      Body *new_body = make_Body( surf );
      kept_surface_bodies.append( new_body );
    }
  }

    // Update DAG
  CubitStatus stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

  if(!preview)
    do_attribute_cleanup();

  if( keep_surface )
    new_body_list += kept_surface_bodies;

  if( stat == CUBIT_FAILURE)
      return CUBIT_FAILURE;

    // Update graphics
  while (entities_to_update.size())
    entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}


//=============================================================================
// Description: Function to remove curves from a sheet body and then extend the
//              remaining curves or fill the gap or hole.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_remove( DLIList<RefEdge*> &ref_edge_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  // clear any preview previews
  GfxPreview::clear();

  gme_ptr = tweak_setup( ref_edge_list, "Removing", old_body_list, curve_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_edge_list );
  }

  // Do remove
  DLIList<BodySM*> new_bodysm_list;
  CubitStatus stat = gme_ptr->tweak_remove( curve_list, new_bodysm_list,
                                            keep_old_body, preview );
  //collect all neighboring curves to those in the list
   int i, j;
   DLIList<RefEdge*> neighboring_curves;
   for( i=ref_edge_list.size(); i--; )
   {
     RefEdge *tmp_edge = ref_edge_list.get_and_step();
     DLIList<RefVertex*> tmp_ref_vertex_list;
     tmp_edge->ref_vertices( tmp_ref_vertex_list );
     for( j=tmp_ref_vertex_list.size(); j--; )
       tmp_ref_vertex_list.get_and_step()->ref_edges( neighboring_curves );
   }

  //uniquify and add other curves
  neighboring_curves.uniquify_unordered();
  //neighboring_curves += ref_edge_list;

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
    if( stat == CUBIT_FAILURE )
      CubitUndo::remove_last_undo();

  if( stat == CUBIT_FAILURE)
    return stat;

  // loop body sm list and find curves that need updating.
  // this is to account for some cases where the topology doesn't change, but the geometry does.
  DLIList<RefEntity*> entities_to_update;
  for(i=0; i<new_bodysm_list.size(); i++)
  {
    BodySM* bodysm = new_bodysm_list.get_and_step();
    DLIList<Curve*> curves;
    bodysm->curves(curves);
    int j;
    // find a surface that is also found in our input list
    for(j=0; j<curves.size(); j++, curves.step())
    {
      BridgeManager* man = curves.get()->bridge_manager();
      if(man)
      {
        RefEdge* ref_edge = CAST_TO(man->topology_entity(), RefEdge);
        if( ref_edge && neighboring_curves.is_in_list(ref_edge) )
        {
          // get neighbors
          DLIList<Point*> neighbor_points;
          curves.get()->points(neighbor_points);
          DLIList<Curve*> neighbor_curves;
          DLIList<Surface*> neighbor_surfaces;
          DLIList<TopologyBridge*> neighbors;
          DLIList<TopologyBridge*> tmp;
          int k;
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->surfaces(neighbor_surfaces);
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->curves(neighbor_curves);

          CAST_LIST_TO_PARENT(neighbor_points, tmp);
          neighbors += tmp;
          neighbor_curves.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_curves, tmp);
          neighbors += tmp;
          neighbor_surfaces.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_surfaces, tmp);
          neighbors += tmp;
          //neighbors.append(curves.get()->lump());
          //neighbors.append(curves.get()->bodysm());

          for(k=0; k<neighbors.size(); k++)
          {
            if(BridgeManager* m = neighbors.get_and_step()->bridge_manager())
            {
              if(TopologyEntity* t = m->topology_entity())
              {
                entities_to_update.append(CAST_TO(t, RefEntity));
                //RefEntity *ref_ent = CAST_TO(t, RefEntity );
              }
            }
          }
        }
      }
    }
  }

  if(preview == CUBIT_FALSE )
  {
    // Update DAG
    stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    // Update graphics
    while (entities_to_update.size())
      entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );

    if( CubitUndo::get_undo_enabled()  )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old_body )
        CubitUndo::note_result_bodies( new_body_list );
    }
  }

  return stat;
}


//=============================================================================
// Description: Tweak specified faces of a volume or volumes up to a target
//              surface or set of connected target surfaces.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_target( DLIList<RefFace*> &ref_face_list,
                                              DLIList<RefFace*> &target_face_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean extend_flg,
                                              CubitPlane *limit_plane,
                                              CubitBoolean reverse_flg,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<Surface*> surface_list(ref_face_list.size());
  DLIList<Surface*> target_surf_list(target_face_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine *gme_ptr1, *gme_ptr2;

  gme_ptr1 = tweak_setup( ref_face_list, "Tweaking", old_body_list, surface_list, CUBIT_TRUE );
  if (!gme_ptr1)
    return CUBIT_FAILURE;

  DLIList<Body*> old_body_list2;
  gme_ptr2 = tweak_setup( target_face_list, "Tweaking", old_body_list2, target_surf_list, CUBIT_TRUE );
  if (!gme_ptr2)
    return CUBIT_FAILURE;

  if( gme_ptr1 != gme_ptr2 )
  {
    PRINT_ERROR( "Target surfaces must belong to same geometry engine as tweaked surfaces.\n" );
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old_body )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_face_list );
  }

  int i;
  DLIList<Body*> all_bodies = old_body_list;
  all_bodies += old_body_list2;
  DLIList<BodySM*> body_sms;
  for(i=all_bodies.size(); i--;)
  {
    BodySM* bsm = all_bodies.get_and_step()->get_body_sm_ptr();
    if(bsm)
      body_sms.append_unique(bsm);
  }

  if(!preview)
  {
    do_attribute_setup();
    push_vg_attributes_before_modify(body_sms);
  }

  // Do tweak to target
  DLIList<BodySM*> new_bodysm_list;
  CubitStatus stat = gme_ptr1->tweak_target( surface_list, target_surf_list,
                                             new_bodysm_list, extend_flg,
                                             limit_plane, reverse_flg, 
                                             keep_old_body, preview );


  if( stat == CUBIT_FAILURE )
  {
    if(!preview)
      remove_pushed_attributes(new_bodysm_list, all_bodies);
    if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
      CubitUndo::remove_last_undo();
    if(!preview)
      do_attribute_cleanup();
    return CUBIT_FAILURE;
  }
  else
  {
    if(!preview)
    {
      restore_vg_after_modify(new_bodysm_list, all_bodies, gme_ptr1);
      remove_pushed_attributes(new_bodysm_list, all_bodies);
    }
  }

  // Collect all the old_faces to be compared against new_faces later...DJQ
  DLIList<RefFace*> old_faces;
  for (i = 0; i < old_body_list.size(); i++)
  {
      Body *old_body = old_body_list.get_and_step();
      DLIList<RefFace*> tmp_faces;
      old_body->ref_faces(tmp_faces);
      old_faces +=tmp_faces;
  }
  
  // Update DAG
  stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

  if(!preview)
    do_attribute_cleanup();

  if( stat == CUBIT_FAILURE )
  {
    if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
      CubitUndo::remove_last_undo();
    return CUBIT_FAILURE;
  }

  //collect all the new faces
  DLIList<RefFace*> new_faces;
  for( i=new_body_list.size(); i--; )
  {
    Body *new_body = new_body_list.get_and_step();
    DLIList<RefFace*> tmp_faces;
    new_body->ref_faces( tmp_faces );
    new_faces += tmp_faces;
  }

  // Compare the new_faces list with the old_faces list to determine which faces are created
  // Add these faces to the all_ref_face_list to check for its neighbors...DJQ
  DLIList<RefFace*> difference = new_faces;
  difference -= old_faces;
  ref_face_list += difference;


  // loop body sm list and find surfaces that need updating.
  // this is to account for some cases where the topology doesn't change, but the geometry does.
  DLIList<RefEntity*> entities_to_update;
  for(i=0; i<new_bodysm_list.size(); i++)
  {
    BodySM* bodysm = new_bodysm_list.get_and_step();
    DLIList<Surface*> surfs;
    bodysm->surfaces(surfs);
    int j;
    // find a surface that is also found in our input list
    // some times the target surface gets modified
    for(j=0; j<surfs.size(); j++, surfs.step())
    {
      BridgeManager* man = surfs.get()->bridge_manager();
      if(man)
      {
        RefFace* ref_face = CAST_TO(man->topology_entity(), RefFace);
        if(ref_face && (ref_face_list.is_in_list(ref_face) ||
            target_face_list.is_in_list(ref_face)))
        {
          // get neighbors
          DLIList<Point*> neighbor_points;
          surfs.get()->points(neighbor_points);
          DLIList<Curve*> neighbor_curves;
          DLIList<Surface*> neighbor_surfaces;
          DLIList<TopologyBridge*> neighbors;
          DLIList<TopologyBridge*> tmp;
          int k;
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->surfaces(neighbor_surfaces);
          for(k=0; k<neighbor_points.size(); k++)
            neighbor_points.get_and_step()->curves(neighbor_curves);

          //for (int i = 0; i < neighbor_points.size(); i++)
          //{
          //    int id = neighbor_points.get_and_step()->get_saved_id();
          //    int temp = id;
          //}
          CAST_LIST_TO_PARENT(neighbor_points, tmp);
          neighbors += tmp;
          neighbor_curves.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_curves, tmp);
          neighbors += tmp;
          neighbor_surfaces.uniquify_unordered();
          CAST_LIST_TO_PARENT(neighbor_surfaces, tmp);
          neighbors += tmp;
          neighbors.append(surfs.get()->lump());
          neighbors.append(surfs.get()->bodysm());

          for(k=0; k<neighbors.size(); k++)
            if(BridgeManager* m = neighbors.get_and_step()->bridge_manager())
              if(TopologyEntity* t = m->topology_entity())
                entities_to_update.append(CAST_TO(t, RefEntity));
        }
      }
    }
  }

  //unmerge any merged adjacent surfaces or
  //merged curves in unmerged adjacent surfaces
  DLIList<RefFace*> adjacent_faces_to_unmerge;
  DLIList<RefEdge*> adjacent_edges_to_unmerge;
  for(i=ref_face_list.size(); i--;)
  {
    RefFace *tweaked_face = ref_face_list.get_and_step();
    if( !new_faces.move_to( tweaked_face ) )
      continue;

    //get all the edges of the face you tweaked
    DLIList<RefEdge*> tweaked_face_edges;
    tweaked_face->ref_edges( tweaked_face_edges );
    adjacent_edges_to_unmerge += tweaked_face_edges;

    //get all the adjacent faces to this edge
    int j;
    for( j=tweaked_face_edges.size(); j--; )
    {
      RefEdge *tmp_edge = tweaked_face_edges.get_and_step();
      DLIList<RefFace*> tmp_faces;
      tmp_edge->ref_faces( tmp_faces );
      tmp_faces.remove( tweaked_face );
      adjacent_faces_to_unmerge += tmp_faces;
    }

    //get all edges not in the surface,
    //sharing vertices with the surface
    DLIList<RefVertex*> ref_vertices;
    tweaked_face->ref_vertices( ref_vertices );
    for( j=ref_vertices.size(); j--; )
    {
      RefVertex *tmp_vert = ref_vertices.get_and_step();
      DLIList<RefEdge*>  ref_edges;
      tmp_vert->ref_edges( ref_edges );

      int k;
      for( k=ref_edges.size(); k--; )
      {
        RefEdge *tmp_edge = ref_edges.get_and_step();
        if( !tweaked_face_edges.move_to( tmp_edge ) )
          adjacent_edges_to_unmerge.append( tmp_edge );
      }
    }
  }

  //unmerge any adjacent faces
  adjacent_faces_to_unmerge.uniquify_unordered();
  for( i=adjacent_faces_to_unmerge.size(); i--; )
  {
    RefFace *ref_face = adjacent_faces_to_unmerge.get_and_step();

    DLIList<TopologyBridge*> bridge_list;
    ref_face->bridge_manager()->get_bridge_list(bridge_list);
    if (bridge_list.size() > 1)
    {
      if( MergeTool::instance()->unmerge( ref_face ) )
        PRINT_WARNING("Unmerging Surface %d\n", ref_face->id() );
    }
  }

  //unmerge any adjacent edges
  adjacent_edges_to_unmerge.uniquify_unordered();
  for( i=adjacent_edges_to_unmerge.size(); i--; )
  {
    RefEdge *ref_edge = adjacent_edges_to_unmerge.get_and_step();
    DLIList<TopologyBridge*> bridge_list;
    ref_edge->bridge_manager()->get_bridge_list(bridge_list);
    if (bridge_list.size() > 1)
    {
      if( MergeTool::instance()->unmerge( ref_edge) )
        PRINT_WARNING("Unmerging Curve %d\n", ref_edge->id() );
    }
  }

  if( CubitUndo::get_undo_enabled() && keep_old_body )
    CubitUndo::note_result_bodies( new_body_list );

  // Update graphics
  while (entities_to_update.size())
    entities_to_update.pop()->notify_all_observers( GEOMETRY_MODIFIED );


  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified faces of a volume or volumes up to a target
//              plane.
// Author     : Steve Storm
// Date       : 05/10/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_target( DLIList<RefFace*> &ref_face_list,
                                              CubitPlane &plane,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean reverse_flg,
                                              CubitBoolean keep_old_body,
                                              CubitBoolean preview )
{
  DLIList<Surface*> surface_list(ref_face_list.size()+1);
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_face_list, "Tweaking", old_body_list, surface_list, CUBIT_TRUE );
  if (!gme_ptr)
    return CUBIT_FAILURE;

  // Create a target Surface from the plane (arbitrarily 10.0 X 10.0 in size)

  // Get corners of the surface
  CubitVector normal = plane.normal();
  CubitVector x, y;
  normal.orthogonal_vectors( x, y );
  CubitVector p1 = plane.point_on_plane();
  CubitVector p2, p3, p4;
  p1.next_point( x, 5.0, p1 );
  p1.next_point( y, 5.0, p1 );
  p1.next_point( -x, 10.0, p2 );
  p2.next_point( -y, 10.0, p3 );
  p3.next_point( x, 10.0, p4 );

  BodySM* bodysm_ptr = gme_ptr->planar_sheet( p1, p2, p3, p4 );
  if( !bodysm_ptr )
  {
    PRINT_ERROR( "unable to create target surface from plane\n" );
    return CUBIT_FAILURE;
  }

  DLIList<Surface*> target_surf_list;
  bodysm_ptr->surfaces( target_surf_list );
  if( !target_surf_list.size() )
  {
    PRINT_ERROR( "unable to create target surface from plane\n" );
    return CUBIT_FAILURE;
  }

  int i;
  DLIList<BodySM*> body_sms;
  for(i=old_body_list.size(); i--;)
  {
    BodySM* bsm = old_body_list.get_and_step()->get_body_sm_ptr();
    if(bsm)
      body_sms.append_unique(bsm);
  }

  if(!preview)
  {
    do_attribute_setup();
    push_vg_attributes_before_modify(body_sms);
  }

  // Do tweak to target
  DLIList<BodySM*> new_bodysm_list;
  CubitStatus stat = gme_ptr->tweak_target( surface_list, target_surf_list,
                                            new_bodysm_list, CUBIT_TRUE,
                                            NULL, reverse_flg, keep_old_body,
                                            preview );


  if( stat == CUBIT_FAILURE )
  {
    if(!preview)
      remove_pushed_attributes(new_bodysm_list, old_body_list);
    if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
      CubitUndo::remove_last_undo();
    if(!preview)
      do_attribute_cleanup();
    return CUBIT_FAILURE;
  }
  else
  {
    if(!preview)
    {
      restore_vg_after_modify(new_bodysm_list, old_body_list, gme_ptr);
      remove_pushed_attributes(new_bodysm_list, old_body_list);
    }
  }

  // Delete temporary sheet body
  bodysm_ptr->get_geometry_query_engine()->delete_solid_model_entities( bodysm_ptr );

  // Update DAG
  stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

  if(!preview)
    do_attribute_cleanup();

  if( stat == CUBIT_FAILURE )
  {
    if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
      CubitUndo::remove_last_undo();
    return CUBIT_FAILURE;
  }

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a target surface or set of target surfaces.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_target( DLIList<RefEdge*> &ref_edge_list,
                                              DLIList<RefFace*> &target_face_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean extend_flg,
                                              CubitPlane *limit_plane,
                                              CubitBoolean reverse_flg,
                                              CubitBoolean keep_old,
                                              CubitBoolean preview,
                                              double max_area_increase /*= 0%*/ )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Surface*> target_surf_list(target_face_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine *gme_ptr1, *gme_ptr2;

  gme_ptr1 = tweak_setup( ref_edge_list, "Tweaking", old_body_list, curve_list );
  if( !gme_ptr1 )
    return CUBIT_FAILURE;

  DLIList<Body*> old_body_list2;
  gme_ptr2 = tweak_setup( target_face_list, "Tweaking", old_body_list2, target_surf_list );
  if( !gme_ptr2 )
    return CUBIT_FAILURE;

  if( gme_ptr1 != gme_ptr2 )
  {
    PRINT_ERROR( "Target surface(s) must belong to same geometry engine as tweaked curves.\n" );
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_edge_list );
  }

  // Do tweak to target
  DLIList<BodySM*> new_bodysm_list;
  CubitStatus stat = gme_ptr1->tweak_target( curve_list, target_surf_list,
                                             new_bodysm_list, extend_flg,
                                             limit_plane, reverse_flg,
                                             keep_old, preview, max_area_increase );

  if( stat == CUBIT_FAILURE )
  {
    if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
      CubitUndo::remove_last_undo();
    return CUBIT_FAILURE;
  }

  // Update DAG
  stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );


  if( stat == CUBIT_FAILURE )
  {
    if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
      CubitUndo::remove_last_undo();
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old )
      CubitUndo::note_result_bodies( new_body_list );
  }

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a target plane.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_target( DLIList<RefEdge*> &ref_edge_list,
                                              CubitPlane &plane,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean reverse_flg,
                                              CubitBoolean keep_old,
                                              CubitBoolean preview )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine* gme_ptr;

  gme_ptr = tweak_setup( ref_edge_list, "Tweaking", old_body_list, curve_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  // Create a target Surface from the plane (arbitrarily 10.0 X 10.0 in size)

  // Get corners of the surface
  CubitVector normal = plane.normal();
  CubitVector x, y;
  normal.orthogonal_vectors( x, y );
  CubitVector p1 = plane.point_on_plane();
  CubitVector p2, p3, p4;
  p1.next_point( x, 5.0, p1 );
  p1.next_point( y, 5.0, p1 );
  p1.next_point( -x, 10.0, p2 );
  p2.next_point( -y, 10.0, p3 );
  p3.next_point( x, 10.0, p4 );

  BodySM* bodysm_ptr = gme_ptr->planar_sheet( p1, p2, p3, p4 );
  if( !bodysm_ptr )
  {
    PRINT_ERROR( "unable to create target surface from plane\n" );
    return CUBIT_FAILURE;
  }

  DLIList<Surface*> target_surf_list;
  bodysm_ptr->surfaces( target_surf_list );
  if( !target_surf_list.size() )
  {
    PRINT_ERROR( "unable to create target surface from plane\n" );
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_edge_list );
  }

  // Do tweak to target
  DLIList<BodySM*> new_bodysm_list;
  CubitStatus stat = gme_ptr->tweak_target( curve_list, target_surf_list,
                                            new_bodysm_list, CUBIT_TRUE,
                                            NULL, reverse_flg, keep_old,
                                            preview );

  if( stat == CUBIT_FAILURE )
  {
    if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
      CubitUndo::remove_last_undo();
    return CUBIT_FAILURE;
  }

  // Delete temporary sheet body
  bodysm_ptr->get_geometry_query_engine()->delete_solid_model_entities( bodysm_ptr );

  // Update DAG
  stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

  if( stat == CUBIT_FAILURE )
  {
    if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
      CubitUndo::remove_last_undo();
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() && keep_old )
    CubitUndo::note_result_bodies( new_body_list );

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Description: Tweak specified edges of a sheet body or bodies up to target
//              curves that are part of a sheet body.  The target is a surface
//              created by thickening the owning surface of the target curves.
// Author     : Steve Storm
// Date       : 03/25/05
//=============================================================================
CubitStatus GeometryModifyTool::tweak_target( DLIList<RefEdge*> &ref_edge_list,
                                              DLIList<RefEdge*> &target_edge_list,
                                              DLIList<Body*> &new_body_list,
                                              CubitBoolean extend_flg,
                                              CubitPlane *limit_plane,
                                              CubitBoolean reverse_flg,
                                              CubitBoolean keep_old,
                                              CubitBoolean preview,
                                              double max_area_increase /*= 0*/ )
{
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<Curve*> target_curve_list(target_edge_list.size());
  DLIList<Body*> old_body_list;
  GeometryModifyEngine *gme_ptr1, *gme_ptr2;

  gme_ptr1 = tweak_setup( ref_edge_list, "Tweaking", old_body_list, curve_list );
  if( !gme_ptr1 )
    return CUBIT_FAILURE;

  DLIList<Body*> old_body_list2;
  gme_ptr2 = tweak_setup( target_edge_list, "Tweaking", old_body_list2, target_curve_list );
  if( !gme_ptr2 )
    return CUBIT_FAILURE;

  if( gme_ptr1 != gme_ptr2 )
  {
    PRINT_ERROR( "Target curves must belong to same geometry engine as tweaked curves.\n" );
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_edge_list );
  }

  // Do tweak to target
  DLIList<BodySM*> new_bodysm_list;
  CubitStatus stat = gme_ptr1->tweak_target( curve_list, target_curve_list,
                                             new_bodysm_list, extend_flg,
                                             limit_plane, reverse_flg,
                                             keep_old, preview, max_area_increase );

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
    if( stat == CUBIT_FAILURE )
      CubitUndo::remove_last_undo();

  if( stat == CUBIT_FAILURE)
    return stat;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    if( CubitUndo::get_undo_enabled()  )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old )
        CubitUndo::note_result_bodies( new_body_list );
    }
  }

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return stat;
}

//=============================================================================
// Description: Tweak specified vertex of a sheet body to a given location.
//              The given vertex must be part of a planar surface or surfaces
//              attached to linear curves only.  The given location will be
//              projected to be on the planar surface(s) before being used.
// Author     : Steve Storm
// Date       : 09/09/08
//=============================================================================
CubitStatus
GeometryModifyTool::tweak_target( RefVertex *ref_vertex_ptr,
                                  DLIList<RefFace*> &modify_ref_face_list,
                                  CubitVector &target_loc,
                                  Body *&new_Body_ptr,
                                  CubitBoolean keep_old,
                                  CubitBoolean preview )
{
  if( modify_ref_face_list.size() == 0 )
  {
    PRINT_ERROR( "No surfaces found to modify\n" );
    return CUBIT_FAILURE;
  }

  int i;
  RefFace *ref_face_ptr;

  // Check to make sure vertex is part of all given surfaces
  modify_ref_face_list.reset();
  for( i=modify_ref_face_list.size(); i--; )
  {
    ref_face_ptr = modify_ref_face_list.get_and_step();

    if( !ref_face_ptr->is_directly_related( ref_vertex_ptr ) )
    {
      PRINT_ERROR( "Vertex %d is not part of 'modify' Surface %d\n", 
        ref_vertex_ptr->id(), ref_face_ptr->id() );
      return CUBIT_FAILURE;
    }
  }
  
  GeometryModifyEngine *gme_ptr;
  DLIList<RefVertex*> ref_vertex_list(1);
  ref_vertex_list.append( ref_vertex_ptr );
  DLIList<Body*> old_body_list;
  DLIList<Point*> point_list(1);
  gme_ptr = tweak_setup( ref_vertex_list, "Tweaking", old_body_list, point_list );
  if( !gme_ptr )
    return CUBIT_FAILURE;

  // We already made sure the vertex is part of all of the modify faces, so
  // just use the common_modify_engine function to get the surface_list
  DLIList<Surface*> surface_list;
  if( !common_modify_engine( modify_ref_face_list, surface_list ) )
    return CUBIT_FAILURE;

  // Make sure part of a sheet body, not a solid body
  Body *body_ptr = ref_vertex_ptr->body();
  if( !body_ptr->is_sheet_body() )
  {
    PRINT_ERROR( "Vertex %d is not in a sheet body - Body %d is solid.\n"
      "       Tweaking a vertex to a target currently not possible on solid bodies.\n",
      ref_vertex_ptr->id(), body_ptr->id() );
    return CUBIT_FAILURE;
  }

  // Make sure all the given surfaces are planar
  modify_ref_face_list.reset();
  for( i=modify_ref_face_list.size(); i--; )
  {
    ref_face_ptr = modify_ref_face_list.get_and_step();

    if( !ref_face_ptr->is_planar() )
    {
      PRINT_ERROR( "Surfaces to modify must be planar - Surface %d is not planar\n", ref_face_ptr->id() );
      return CUBIT_FAILURE;
    }
  }

  // Make sure all curves (on modify surfaces) attached to the given vertex are linear
  // Get all attached curves
  int j;
  DLIList<RefEdge*> ref_edge_list;
  ref_vertex_ptr->ref_edges( ref_edge_list );
  RefEdge *ref_edge_ptr;
  for( i=ref_edge_list.size(); i--; )
  {
    ref_edge_ptr = ref_edge_list.get_and_step();

    // Check to see if this edge is linear, if it is in one of the modify surfaces
    modify_ref_face_list.reset();
    for( j=modify_ref_face_list.size(); j--; )
    {
      ref_face_ptr = modify_ref_face_list.get_and_step();
      if( ref_face_ptr->is_directly_related( ref_edge_ptr ) )
      {
        Curve *curve_ptr = ref_edge_ptr->get_curve_ptr();
        GeometryType curve_type = curve_ptr->geometry_type();
        if( curve_type != STRAIGHT_CURVE_TYPE )
        {
          PRINT_ERROR( "Curve %d is not linear. Curves that are in the 'modify'\n"
            "       surfaces attached to the tweak vertex must be linear.\n", 
            ref_edge_ptr->id() );
          return CUBIT_FAILURE;
        }
        else
          break;
      }
    }
  }

  // Project the location to the given surfaces and make sure all of these locations
  // are the same
  modify_ref_face_list.reset();

  CubitVector ref_loc( target_loc ); // Reference location
  ref_face_ptr = modify_ref_face_list.get_and_step();
  ref_face_ptr->move_to_surface( ref_loc );
  int ref_surf_id = ref_face_ptr->id();

  for( i=modify_ref_face_list.size()-1; i--; )
  {
    ref_face_ptr = modify_ref_face_list.get_and_step();

    CubitVector proj_loc( target_loc );
    ref_face_ptr->move_to_surface( proj_loc );

    if( !ref_loc.about_equal( proj_loc ) )
    {
      PRINT_ERROR( "Target location must project to all of the 'modify' surfaces at\n"
        "       exactly the same location - the tweaked Vertex %d will move to this\n"
        "       common (same) projected location, tweaking the modify surfaces with it.\n"
        "       Given target location = %f, %f, %f\n"
        "       Projected location on Surface %d = %f, %f, %f\n"
        "       Projected location on Surface %d = %f, %f, %f\n",
        ref_vertex_ptr->id(),
        target_loc.x(), target_loc.y(), target_loc.z(),
        ref_surf_id, ref_loc.x(), ref_loc.y(), ref_loc.z(),
        ref_face_ptr->id(), proj_loc.x(), proj_loc.y(), proj_loc.z() );

      return CUBIT_FAILURE;
    }
  }

  // It looks like the inputs are valid so get ready to do the tweak
  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
  {
    if( keep_old )
      CubitUndo::save_state();
    else
      CubitUndo::save_state_with_cubit_file( ref_vertex_list );
  }

  Point *point_ptr = point_list.get();

  // Do tweak to target
  BodySM *new_bodysm_ptr;
  CubitStatus stat = gme_ptr->tweak_target( point_ptr, surface_list, 
                                            ref_loc,
                                            new_bodysm_ptr, 
                                            keep_old, preview );

  if( CubitUndo::get_undo_enabled() && preview == CUBIT_FALSE )
    if( stat == CUBIT_FAILURE )
      CubitUndo::remove_last_undo();

  if( stat == CUBIT_FAILURE)
    return stat;

  DLIList<BodySM*> new_bodysm_list;
  new_bodysm_list.append( new_bodysm_ptr );
  DLIList<Body*> new_body_list;

  if( preview == CUBIT_FALSE )
  {
    // Update DAG
    stat = finish_sm_op( old_body_list, new_bodysm_list, new_body_list );

    if( CubitUndo::get_undo_enabled() )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
      else if( keep_old )
        CubitUndo::note_result_bodies( new_body_list );
    }
  }

  // Update graphics
  DLIList<Body*> moved_bodies(new_body_list);
  moved_bodies.intersect(old_body_list);
  while (moved_bodies.size())
    moved_bodies.pop()->notify_sub_all_observers( GEOMETRY_MODIFIED );

  return stat;
}

// KGM
#if 0
bool GeometryModifyTool::contains_intermediate_geometry(DLIList<RefFace*>& ref_face_list) const
{
  for (int i = 0; i < ref_face_list.size(); i++)
    if (GeometryQueryTool::instance()->
        contains_intermediate_geometry(ref_face_list.next(i)))
      return true;

  return false;
}
#endif

/*
//the following surface tool operations added by Tyronne Lim (CAT) ********************
CubitStatus GeometryModifyTool::create_net_surface( DLIList<Surface*>& ref_face_list, BodySM *& new_body,
                                                    DLIList<DLIList<CubitVector*>*> &vec_lists_u,
                                                    DLIList<DLIList<CubitVector*>*> &vec_lists_v,
                                                    double net_tol, CubitBoolean heal )
{ 
   GeometryModifyEngine* GMEPtr = get_engine(ref_face_list.get());
   return GMEPtr->create_net_surface( ref_face_list, new_body, vec_lists_u, vec_lists_v, net_tol, heal );
} */

CubitStatus GeometryModifyTool::create_net_surface( DLIList<RefEdge*>& u_curves, DLIList<RefEdge*>& v_curves,
                                                    double net_tol, CubitBoolean heal )
{
  DLIList<TopologyBridge*> bridge_list;
  DLIList<TopologyEntity*> entity_list;
  CAST_LIST_TO_PARENT( u_curves, entity_list );

  GeometryModifyEngine* GME_ptr =
    common_modify_engine( entity_list, bridge_list );
  if(! GME_ptr )
  {
     PRINT_ERROR("Cannot construct a Surface using entities that do "
                 "not share a common GeometryModifyEngine.\n");
     return CUBIT_FAILURE;
  }

  DLIList<Curve*> curves_in_u(bridge_list.size());
  CAST_LIST( bridge_list, curves_in_u, Curve );

  bridge_list.clean_out();
  entity_list.clean_out();
  CAST_LIST_TO_PARENT( v_curves, entity_list );

  GeometryModifyEngine* dummy_GME_ptr =
    common_modify_engine( entity_list, bridge_list );
  if(! dummy_GME_ptr || dummy_GME_ptr != GME_ptr )
  {
     PRINT_ERROR("Cannot construct a Surface using entities that do "
                 "not share a common GeometryModifyEngine.\n");
     return CUBIT_FAILURE;
  }

  DLIList<Curve*> curves_in_v(bridge_list.size());
  CAST_LIST( bridge_list, curves_in_v, Curve );

  BodySM *new_smbody = NULL;
  if( !GME_ptr->create_net_surface( curves_in_u, curves_in_v, new_smbody, net_tol, heal ) )
    return CUBIT_FAILURE;

  if (CubitUndo::get_undo_enabled())
    CubitUndo::save_state();

  Body* new_body = NULL;
  new_body = GeometryQueryTool::instance()->make_Body( new_smbody );
  if( new_body )
  {
    if (CubitUndo::get_undo_enabled())
      CubitUndo::note_result_body(new_body);

    return CUBIT_SUCCESS;
  }
  else
  {
    if (CubitUndo::get_undo_enabled())
      CubitUndo::remove_last_undo();

    return CUBIT_FAILURE;
  }
}

CubitStatus GeometryModifyTool::create_offset_surface( RefFace* ref_face_ptr, double offset_distance )
{
  TopologyBridge *bridge_ptr = NULL;
  GeometryModifyEngine* GMEPtr = get_engine(ref_face_ptr, &bridge_ptr );

  Surface *tmp_surf = NULL;
  tmp_surf = CAST_TO( bridge_ptr, Surface );

  BodySM *new_smbody;
  if( !GMEPtr->create_offset_surface( tmp_surf, new_smbody, offset_distance ) )
    return CUBIT_FAILURE;

  if (CubitUndo::get_undo_enabled())
    CubitUndo::save_state();

  Body* new_body = NULL;
  new_body = GeometryQueryTool::instance()->make_Body( new_smbody );
  if( new_body )
  {
    if (CubitUndo::get_undo_enabled())
      CubitUndo::note_result_body(new_body);

    return CUBIT_SUCCESS;
  }
  else
  {
    if (CubitUndo::get_undo_enabled())
      CubitUndo::remove_last_undo();

    return CUBIT_FAILURE;
  }
}

//-----------------------------------------------------------------------------
// Purpose       : This function creates a sheet body or bodies by offsetting
//                 the given faces. An optional additional face list and double
//                 list (must be same length) allow different offset distances
//                 for different faces. Adjoining faces are extended or trimmed
//                 to remain joined in the new sheet body.  Radial faces that
//                 cannot be so offset are removed and the resulting wound
//                 healed by the surrounding faces.
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 05/04/08
//-----------------------------------------------------------------------------
CubitStatus
GeometryModifyTool::create_offset_sheet( DLIList<RefFace*> &ref_face_list,
                                         double offset_distance,
                                         DLIList<RefFace*> *add_ref_face_list_ptr,
                                         DLIList<double> *add_offset_list_ptr,
                                         DLIList<Body*> &new_body_list,
                                         CubitBoolean preview )
{
  if( !ref_face_list.size() )
  {
    return CUBIT_SUCCESS;
  }
  else
  {
    for (int i = 0; i < ref_face_list.size(); i++)
    {
      RefFace *temp = ref_face_list.get_and_step();
      if (temp->get_surface_ptr()->is_cylindrical() == 16)
      {
        DLIList<Curve*> curves;
        CubitVector loc;
        double rad;
        temp->get_surface_ptr()->curves(curves);
        curves.reset();
        int j = 0;
        while (curves.get()->geometry_type() != ELLIPSE_CURVE_TYPE && curves.get()->geometry_type() != ARC_CURVE_TYPE)
        {
          curves.step();
          if (j == curves.size())
          {
              break;
          }   
          j++;
        }
        curves.get()->get_center_radius(loc, rad);
        if (offset_distance >= rad)
        {
          CubitVector norm, close, result;
          double angle;
          temp->get_surface_ptr()->closest_point(loc, &close, &norm);
          result.set(loc.x()-close.x(), loc.y()-close.y(), loc.z()-close.z());
          angle = result.interior_angle(norm);
          if (angle < GEOMETRY_RESABS && angle > -GEOMETRY_RESABS)
          {
            PRINT_ERROR("Offset is greater than the radius of curvature for surface %i.\n", temp->id());
            PRINT_WARNING("No body will be created for surface %i.\n", temp->id());
            if (ref_face_list.size() > 1)
            {
              ref_face_list.remove_all_with_value(temp);
            }
            else
            {
              return CUBIT_FAILURE;
            }
          }
        }
      }
    }
    for (int i = 0; i < add_ref_face_list_ptr->size(); i++)
    {
      RefFace *temp = ref_face_list.get_and_step();
      if (temp->get_surface_ptr()->is_cylindrical() == 16)
      {
        DLIList<Curve*> curves;
        CubitVector loc;
        double rad;
        temp->get_surface_ptr()->curves(curves);
        curves.reset();
        while (curves.get()->geometry_type() != ELLIPSE_CURVE_TYPE && curves.get()->geometry_type() != ARC_CURVE_TYPE)
        {
          curves.step();
        }
        curves.get()->get_center_radius(loc, rad);
        if (offset_distance >= rad)
        {
          CubitVector norm, close, result;
          double angle;
          temp->get_surface_ptr()->closest_point(loc, &close, &norm);
          result.set(loc.x()-close.x(), loc.y()-close.y(), loc.z()-close.z());
          angle = result.interior_angle(norm);
          if (angle < GEOMETRY_RESABS && angle > -GEOMETRY_RESABS)
          {
            PRINT_ERROR("Offset is greater than the radius of curvature for surface %i.\n", temp->id());
            PRINT_WARNING("No body will be created for surface %i.\n", temp->id());
            add_ref_face_list_ptr->remove_all_with_value(temp);
          }
        }
      }
    }
  }

  DLIList<RefFace*> all_ref_face_list(ref_face_list.size());
  all_ref_face_list = ref_face_list;
  if( add_ref_face_list_ptr->size() )
    all_ref_face_list += *add_ref_face_list_ptr;

  // Check for virtual geometry
  DLIList<RefEntity*> ref_ent_list;
  CAST_LIST_TO_PARENT(all_ref_face_list, ref_ent_list);
  if ( GeometryQueryTool::instance()->contains_intermediate_geometry(ref_ent_list) )
  {
    PRINT_ERROR("OFFSETTING surfaces containing virtual geometry is not\n"
      "       allowed. Delete virtual geometry on these surfaces\n"
      "       before operation.\n" );
    return CUBIT_FAILURE;
  }

  // Look for a common GeometryModifyEngine for all of the RefFaces
  int count = all_ref_face_list.size();
  DLIList<TopologyBridge*> bridge_list(count);
  DLIList<TopologyEntity*> entity_list(count);
  CAST_LIST_TO_PARENT( all_ref_face_list, entity_list );

  GeometryModifyEngine* GME_ptr =
    common_modify_engine( entity_list, bridge_list );
  if(! GME_ptr )
  {
     PRINT_ERROR("Cannot construct offset sheet(s) using surfaces that\n"
                 "       do not share a common geometry engine.\n");
     return CUBIT_FAILURE;
  }

  // Get Surfaces from the RefFaces
  DLIList<Surface*> surface_list(ref_face_list.size());
  int i;
  RefFace *ref_face_ptr;
  for( i=ref_face_list.size(); i--; )
  {
    ref_face_ptr = ref_face_list.get_and_step();
    surface_list.append( ref_face_ptr->get_surface_ptr() );
  }

  DLIList<Surface*> add_surf_list;
  if( add_ref_face_list_ptr->size() )
  {
    for( i=add_ref_face_list_ptr->size(); i--; )
    {
      ref_face_ptr = add_ref_face_list_ptr->get_and_step();
      add_surf_list.append( ref_face_ptr->get_surface_ptr() );
    }
  }

  DLIList<BodySM*> BodySM_list;
  if( add_surf_list.size() )
  {
    if( GME_ptr->create_offset_sheet( surface_list, offset_distance, &add_surf_list, 
      add_offset_list_ptr, BodySM_list, preview ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;
  }
  else
  {
    if( GME_ptr->create_offset_sheet( surface_list, offset_distance, NULL, 
      NULL, BodySM_list, preview ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;
  }

  if( !BodySM_list.size() )
    return CUBIT_FAILURE;

  BodySM *bodysm_ptr;
  Body* body_ptr;
  for( i=BodySM_list.size(); i--; )
  {
    DLIList<Surface*> surfs;
    bodysm_ptr = BodySM_list.get_and_step();
    bodysm_ptr->surfaces(surfs);
    if (!surfs.size())
    {
      PRINT_WARNING("Empty body created.  Body deleted.\n");
      PRINT_WARNING("Empty body likely due to an offset larger than the radius of curvature of a surface.\n");
      bodysm_ptr->get_geometry_query_engine()->delete_solid_model_entities(bodysm_ptr);
      break;
    }

    body_ptr = GeometryQueryTool::instance()->make_Body( bodysm_ptr );
    if( body_ptr )
      new_body_list.append( body_ptr );
  }

  if( new_body_list.size() )
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

CubitStatus
GeometryModifyTool::create_offset_body( Body *body_ptr, Body *&new_body,
                                        double offset_distance )
{
  GeometryModifyEngine* GMEPtr = get_engine(body_ptr);

  BodySM *body_sm = body_ptr->get_body_sm_ptr();
  if (!body_sm)
  {
    PRINT_ERROR("Body %d is invalid -- no attached BodySM.\n", body_ptr->id());
    return CUBIT_FAILURE;
  }

  BodySM *new_body_sm;
  if( !GMEPtr->create_offset_body( body_sm, new_body_sm, offset_distance ) )
    return CUBIT_FAILURE;

  if (CubitUndo::get_undo_enabled())
    CubitUndo::save_state();

  new_body = GeometryQueryTool::instance()->make_Body( new_body_sm );

  if( new_body )
  {
    if (CubitUndo::get_undo_enabled())
      CubitUndo::note_result_body(new_body);

    return CUBIT_SUCCESS;
  }
  else
  {
    if (CubitUndo::get_undo_enabled())
      CubitUndo::remove_last_undo();

    return CUBIT_FAILURE;
  }
}

CubitStatus GeometryModifyTool::create_skin_surface( DLIList<RefEdge*>& ref_edges, Body*& new_body,
                                                     DLIList<RefEdge*>& guides)
{
  DLIList<TopologyBridge*> bridge_list;
  DLIList<TopologyEntity*> entity_list;
  CAST_LIST_TO_PARENT( ref_edges, entity_list );

  DLIList<TopologyBridge*> guide_bridge_list;
  DLIList<TopologyEntity*> guide_entity_list;
  CAST_LIST_TO_PARENT(guides, guide_entity_list);

  if (ref_edges.size() < 2)
  {
     PRINT_ERROR("Must specify at least 2 curves to create a skinned surface.\n");
     return CUBIT_FAILURE;
  }

  GeometryModifyEngine* GME_ptr =
    common_modify_engine( entity_list, bridge_list );
  if(! GME_ptr )
  {
     PRINT_ERROR("Cannot construct a Surface using entities that do "
                 "not share a common GeometryModifyEngine.\n");
     return CUBIT_FAILURE;
  }

  if (guides.size() > 0)
  {
      GeometryModifyEngine* GME_ptr2 =
          common_modify_engine( guide_entity_list, guide_bridge_list );
      if (GME_ptr != GME_ptr2)
      {
          PRINT_ERROR("Performing create skin with geometry from\n"
			"different modeling engines is not allowed.\n");
		  return CUBIT_FAILURE;
      }
  }

  DLIList<Curve*> curves_to_skin(bridge_list.size());
  CAST_LIST( bridge_list, curves_to_skin, Curve );

  DLIList<Curve*> guide_curves(guide_bridge_list.size());
  CAST_LIST(guide_bridge_list, guide_curves, Curve);

  BodySM *new_body_sm = NULL;
  if( !GME_ptr->create_skin_surface( curves_to_skin, new_body_sm, guide_curves ) )
    return CUBIT_FAILURE;

  if (CubitUndo::get_undo_enabled())
    CubitUndo::save_state();

  new_body = NULL;
  new_body = GeometryQueryTool::instance()->make_Body( new_body_sm );

  if( new_body )
  {
    if (CubitUndo::get_undo_enabled())
      CubitUndo::note_result_body(new_body);

    return CUBIT_SUCCESS;
  }
  else
  {
    if (CubitUndo::get_undo_enabled())
      CubitUndo::remove_last_undo();

    return CUBIT_FAILURE;
  }
}


CubitStatus GeometryModifyTool::loft_surfaces_to_body( RefFace *face1, const double &takeoff1,
                                                      RefFace *face2, const double &takeoff2,
                                                      Body*& new_body,
                                                      CubitBoolean arc_length_option, CubitBoolean twist_option,
                                                      CubitBoolean align_direction, CubitBoolean perpendicular,
                                                      CubitBoolean simplify_option)
{

    DLIList<RefFace*> loft_faces;
    loft_faces.append(face1);
    loft_faces.append(face2);
    DLIList<Surface*> loft_surfaces;

    // Get engine and correspoding geom entities
    GeometryModifyEngine* result_ptr;
    result_ptr = common_modify_engine( loft_faces, loft_surfaces );
    if (!result_ptr)
    {
        PRINT_ERROR("Loft surfaces on volumes containing surfaces from different\n"
            "       geometry engines is not allowed.\n");
        return CUBIT_FAILURE;
    }

    if(2!=loft_surfaces.size())
        return CUBIT_FAILURE;

    loft_surfaces.reset();
    BodySM* new_body_sm = 0;

    if( CubitUndo::get_undo_enabled() )
      CubitUndo::save_state();

    CubitStatus result = result_ptr->loft_surfaces_to_body(
        loft_surfaces.get_and_step(),
        takeoff1,
        loft_surfaces.get_and_step(),
        takeoff2,
        new_body_sm,
        arc_length_option,
        twist_option,
        align_direction,
        perpendicular,
        simplify_option);

    if(result && new_body_sm)
    {
        new_body = GeometryQueryTool::instance()->make_Body(new_body_sm);

        if( CubitUndo::get_undo_enabled() )
          CubitUndo::note_result_body(new_body);
    }
    else
    {
      if( CubitUndo::get_undo_enabled() )
        CubitUndo::remove_last_undo();
    }
   return result;
}

CubitStatus GeometryModifyTool::create_surface( DLIList<RefVertex*> &vert_list, 
                                                Body *&new_body )
{
   //determine which vertices are free and which are not...
   //copy ones that are not
   vert_list.reset();
   DLIList<Point*> points;
   GeometryModifyEngine *GMEPtr = get_engine( vert_list.get()->get_point_ptr() );

   DLIList<RefVertex*> free_ref_vertices;
   int i;
   for( i=vert_list.size(); i--; )
   {
     RefVertex *tmp_vert = vert_list.get_and_step();

     if( tmp_vert->num_parent_ref_entities() == 0 )
       free_ref_vertices.append( tmp_vert );

     Point *tmp_point = tmp_vert->get_point_ptr();
     if( GMEPtr != get_engine( tmp_point ) )
     {
       PRINT_INFO("Vertices are not from same modeling engine.\n");
       return CUBIT_FAILURE;
     }

     if( tmp_vert->get_parents() == 0 )
     {
       points.append( tmp_point ); 
     }
     else
     {
       Point *new_point = GMEPtr->make_Point( tmp_vert->coordinates() );
       points.append( new_point );
     }
   }

   BodySM* body_sm = NULL;
   CubitStatus stat = GMEPtr->create_surface( points, body_sm );

   if( stat == CUBIT_FAILURE )
     return stat;

   if (CubitUndo::get_undo_enabled())
     CubitUndo::save_state();

   if( body_sm )
     new_body = GeometryQueryTool::instance()->make_Body(body_sm);

   if (new_body)
   {
     if (CubitUndo::get_undo_enabled())
       CubitUndo::note_result_body(new_body);

     stat = CUBIT_SUCCESS;
   }
   else
   {
     if (CubitUndo::get_undo_enabled())
       CubitUndo::remove_last_undo();

     stat = CUBIT_FAILURE;
   }

   for( i=free_ref_vertices.size(); i--; )
   {
     RefVertex *free_vertex = free_ref_vertices.get_and_step();
     CubitObserver::notify_static_observers( free_vertex, TOP_LEVEL_ENTITY_DESTRUCTED );
     CGMHistory::Event evt(CGMHistory::TOP_LEVEL_ENTITY_DELETED, free_vertex );
     GeometryQueryTool::instance()->history().add_event(evt);
   }

   return stat;
}

CubitStatus GeometryModifyTool::create_surface( DLIList<CubitVector*>& vec_list, 
                                                Body *&new_body,
                                                RefFace *ref_face_ptr,
                                                CubitBoolean project_points )
{
   GeometryModifyEngine* GMEPtr = gmeList.get();
   if( ref_face_ptr )
   {
     if( GMEPtr != get_engine(ref_face_ptr) )
     {
       PRINT_ERROR("Geometry engine of Surface %d is not the active geometry engine.\n", ref_face_ptr->id() );
       PRINT_INFO("      Use command \"Set Geometry Engine ...\" to set to correct engine.\n");
       return CUBIT_FAILURE;
     }
   }
   BodySM* body_sm = NULL;
   Surface *project_to_surface = NULL;
   if( ref_face_ptr )
     project_to_surface = ref_face_ptr->get_surface_ptr();
   CubitStatus stat = GMEPtr->create_surface( vec_list, body_sm, project_to_surface, project_points );

   if( stat == CUBIT_FAILURE )
     return stat;

   if (CubitUndo::get_undo_enabled())
     CubitUndo::save_state();

   if( body_sm )
     new_body = GeometryQueryTool::instance()->make_Body(body_sm);

   if (new_body)
   {
     if (CubitUndo::get_undo_enabled())
       CubitUndo::note_result_body(new_body);

     stat = CUBIT_SUCCESS;
   }
   else
   {
     if (CubitUndo::get_undo_enabled())
       CubitUndo::remove_last_undo();

     stat = CUBIT_FAILURE;
   }
   return stat;
}

CubitStatus GeometryModifyTool::create_weld_surface( CubitVector &root,
                                                    RefFace *ref_face1,
                                                    double leg1,
                                                    RefFace *ref_face2,
                                                    double leg2,
                                                    Body *&new_body )
{
    GeometryModifyEngine* GMEPtr;
    DLIList<RefFace*> ref_faces;
    ref_faces.append(ref_face1);
    ref_faces.append(ref_face2);
    DLIList<Surface*> surfaces;

    GMEPtr = common_modify_engine(ref_faces, surfaces);
    if (!GMEPtr)
    {
        PRINT_ERROR("Loft surfaces on volumes containing surfaces from different\n"
            "       geometry engines is not allowed.\n");
        return CUBIT_FAILURE;
    }

    if (CubitUndo::get_undo_enabled())
      CubitUndo::save_state();

    surfaces.reset();
    BodySM* new_body_sm = 0;
    CubitStatus result = GMEPtr->create_weld_surface(
        root,
        surfaces.get_and_step(),
        leg1,
        surfaces.get_and_step(),
        leg2,
        new_body_sm );

    if(result && new_body_sm)
    {
        new_body = GeometryQueryTool::instance()->make_Body(new_body_sm);
        if (new_body)
        {
          if (CubitUndo::get_undo_enabled())
            CubitUndo::note_result_body(new_body);
        }
        else
        {
          if (CubitUndo::get_undo_enabled())
            CubitUndo::remove_last_undo();
        }
    }
    return result;
}
//end of surface tool operations ******************************************************

//-------------------------------------------------------------------------
// Purpose       : Remove entity names from dead entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/08/03
//-------------------------------------------------------------------------
void GeometryModifyTool::remove_dead_entity_names( RefEntity* entity ) const
{
  TopologyEntity* topo_ent = dynamic_cast<TopologyEntity*>(entity);
  if (topo_ent->bridge_manager()->topology_bridge() == NULL)
    RefEntityName::instance()->remove_refentity_name( entity, CUBIT_TRUE );

  DLIList<RefEntity*> children;
  entity->get_child_ref_entities(children);
  children.last();
  for (int i = children.size(); i--; )
  {
    //PRINT_INFO("Removing dead entity on %s %d\n", children.get()->class_name(), children.get()->id() );
    remove_dead_entity_names( children.step_and_get() );
  }
}

//-------------------------------------------------------------------------
// Purpose       : Destroy or update modified body, as appropriate.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/08/03
//-------------------------------------------------------------------------
Body* GeometryModifyTool::update_body( Body* body ) const
{
  BodySM* body_sm = body->get_body_sm_ptr();
  if (body_sm)
    return GeometryQueryTool::instance()->make_Body(body_sm);

  GeometryQueryTool::instance()->destroy_dead_entity(body);
  return 0;
}

CubitStatus GeometryModifyTool::tolerant_imprint( DLIList<RefFace*> &ref_faces,
                                                  DLIList<RefEdge*> &ref_edge_list,
                                                  DLIList<Body*> &new_bodies,
                                                  bool merge )
{
  if( ref_faces.size() > 2 )
    return CUBIT_FAILURE;

  ref_faces.reset();
  RefFace *face1 = ref_faces.get_and_step();
  RefFace *face2 = ref_faces.get_and_step();

  if(ref_edge_list.size() > 0)
  {
    //collect all the bodies containing any edge on these 2 surfaces
    //so you can merge them afterward
    DLIList<Body*> bodies_to_merge;
    if( merge )
    {
      DLIList<RefEdge*> tmp_edges;
      face1->ref_edges( tmp_edges);

      int j;
      for( j=tmp_edges.size(); j--; )
      {
        RefEdge *tmp_edge = tmp_edges.get_and_step();
        DLIList<Body*> body_list;
        tmp_edge->bodies( body_list );
        bodies_to_merge += body_list;
      }

      for( j=ref_edge_list.size(); j--; )
      {
        RefEdge *tmp_edge = ref_edge_list.get_and_step();
        DLIList<Body*> body_list;
        tmp_edge->bodies( body_list );
        bodies_to_merge += body_list;
      }

      tmp_edges.clean_out();
      face2->ref_edges( tmp_edges );

      for( j=tmp_edges.size(); j--; )
      {
        RefEdge *tmp_edge = tmp_edges.get_and_step();
        DLIList<Body*> body_list;
        tmp_edge->bodies( body_list );
        bodies_to_merge += body_list;
      }
      bodies_to_merge.uniquify_unordered();
    }

    DLIList<RefEdge*> edges_to_imprint_onto_face1;
    DLIList<RefEdge*> edges_to_imprint_onto_face2;

    //sort edges...
    //edges on face1 and not on face2 will be imprinted on face2
    //edges on face2 and not on face1 will be imprinted on face1
    int i;
    for(i=ref_edge_list.size(); i--; )
    {
      RefEdge *tmp_edge = ref_edge_list.get_and_step();
      DLIList<RefFace*> tmp_faces;
      tmp_edge->ref_faces( tmp_faces );

      if( tmp_faces.move_to( face1 ) && !tmp_faces.move_to( face2 ) )
        edges_to_imprint_onto_face2.append( tmp_edge );
      else if( tmp_faces.move_to( face2 ) && !tmp_faces.move_to( face1 ) )
        edges_to_imprint_onto_face1.append( tmp_edge );
      else
        PRINT_ERROR("Will not imprint curve %d onto either surface.\n", tmp_edge->id() );
    }

    //if there are edges to impint onto both surfaces
    if( edges_to_imprint_onto_face1.size() &&
        edges_to_imprint_onto_face2.size() )
    {
      //get the modify engine
      DLIList<Surface*> surf_list( 1 );
      DLIList<Curve*> curve_list( edges_to_imprint_onto_face2.size() );
      DLIList<RefFace*> ref_face_list( 1 );
      ref_face_list.append( face2 );
      GeometryModifyEngine* gme = common_modify_engine( ref_face_list,
                                                        edges_to_imprint_onto_face2,
                                                        surf_list,
                                                        curve_list );
      if ( !gme )
      {
        PRINT_ERROR("Performing IMPRINT with entities containing geometry from\n"
                    "different modeling engines is not allowed.\n" );
        return CUBIT_FAILURE;
      }

      //copy the specified boundary curves of face1 to imprint onto face2...
      //these could be stale after we imprint face1
      DLIList<Curve*> copied_curves;
      for(i=curve_list.size(); i--; )
      {
        Curve *copied_curve = gme->make_Curve( curve_list.get_and_step() );
        copied_curves.append( copied_curve );
      }

    if( CubitUndo::get_undo_enabled() )
      CubitUndo::save_state_with_cubit_file( ref_faces );

      //do the imprint onto face1
      Body *new_body = NULL;
      CubitStatus status;
      status = tolerant_imprint( face1, edges_to_imprint_onto_face1, new_body );

      //if we failed...get out
      if( status == CUBIT_FAILURE )
      {
        //delete the copied curves
        for( i=copied_curves.size(); i--; )
          gme->get_gqe()->delete_solid_model_entities( copied_curves.get_and_step() );

        if( CubitUndo::get_undo_enabled() )
          CubitUndo::remove_last_undo();

        return CUBIT_FAILURE;
      }

      DLIList<Body*> original_body_list;
      face2->bodies( original_body_list );

      //get the Surface* ptr
      Surface *surface2 = face2->get_surface_ptr();

      DLIList<BodySM*> body_sm_list;
      for(i=original_body_list.size(); i--;)
        body_sm_list.append(original_body_list.get_and_step()->get_body_sm_ptr());

      int process_composites = 0;
      if(contains_composites(original_body_list))
        process_composites = 1;

      if(process_composites)
      {
        // Push virtual attributes down to solid model topology before
        // doing the imprint.
        do_attribute_setup();
        push_vg_attributes_before_modify(body_sm_list);
        // This must be done after pushing the vg atts because it uses them.
        push_imprint_attributes_before_modify(body_sm_list);
      }

      //do the imprint onto face2
      BodySM *new_bodysm = NULL;
      DLIList<TopologyBridge*> new_tbs, att_tbs, temporary_bridges;
      status = gme->tolerant_imprint_surface_with_curves( surface2, copied_curves, temporary_bridges, new_bodysm,
        &new_tbs, &att_tbs);
      temporary_bridges.uniquify_ordered();
      while(temporary_bridges.size())
        delete temporary_bridges.pop();

      DLIList<BodySM*> new_body_list;
      if(new_bodysm)
        new_body_list.append(new_bodysm);


      //delete the copied curves
      for( i=copied_curves.size(); i--; )
        gme->get_gqe()->delete_solid_model_entities( copied_curves.get_and_step() );

      if( status == CUBIT_FAILURE )
      {
        if( CubitUndo::get_undo_enabled() )
          CubitUndo::remove_last_undo();

        if(process_composites)
        {
          remove_pushed_attributes(new_body_list, original_body_list);
          do_attribute_cleanup();
        }

        return CUBIT_FAILURE;
      }
      else
      {
        if(process_composites)
        {
          // Analyze the results and adjust virtual attributes as necessary.
        DLIList<TopologyBridge*> tb_list;
        CAST_LIST(new_body_list, tb_list, TopologyBridge);
          GeometryQueryTool::instance()->ige_attribute_after_imprinting(new_tbs, att_tbs,
            tb_list, original_body_list);

          // Clean up attributes.
          remove_imprint_attributes_after_modify(body_sm_list, new_body_list);

          // Restore the virtual geometry.
          restore_vg_after_modify(new_body_list, original_body_list, gme);
          remove_pushed_attributes(new_body_list, original_body_list);
        }
      }

      DLIList<BodySM*> new_bodies;
      new_bodies.append( new_bodysm );
      DLIList<Body*> result_bodies;
      status = finish_sm_op( original_body_list, new_bodies, result_bodies );

      if(process_composites)
        do_attribute_cleanup();

      if( status == CUBIT_FAILURE )
      {
        if( CubitUndo::get_undo_enabled() )
          CubitUndo::remove_last_undo();
        return CUBIT_FAILURE;
      }
    }
    //user specified edges that will only imprint onto face1...do it anyway
    else if( edges_to_imprint_onto_face1.size() )
    {
      if( CubitUndo::get_undo_enabled() )
      {
        DLIList<RefFace*> tmp_faces(1);
        tmp_faces.append( face1 );
        CubitUndo::save_state_with_cubit_file( tmp_faces );
      }

      bool undo_enabled = CubitUndo::get_undo_enabled();
      CubitUndo::set_undo_enabled( false );

      Body *new_body = NULL;
      CubitStatus stat = tolerant_imprint( face1, edges_to_imprint_onto_face1,
                                          new_body, merge );

    if( undo_enabled )
      CubitUndo::set_undo_enabled( true );

    if( CubitUndo::get_undo_enabled() )
    {
      if( stat == CUBIT_FAILURE )
        CubitUndo::remove_last_undo();
    }

      return stat;
    }
    //user specified edges that will only imprint onto face2...do it anyway
    else if( edges_to_imprint_onto_face2.size() )
    {
      if( CubitUndo::get_undo_enabled() )
      {
        DLIList<RefFace*> tmp_faces(1);
        tmp_faces.append( face2 );
        CubitUndo::save_state_with_cubit_file( tmp_faces );
      }

      bool undo_enabled = CubitUndo::get_undo_enabled();
      CubitUndo::set_undo_enabled( false );

      Body *new_body = NULL;
      CubitStatus stat = tolerant_imprint( face2, edges_to_imprint_onto_face2,
                                          new_body, merge );

      if( undo_enabled )
        CubitUndo::set_undo_enabled( true );

      if( CubitUndo::get_undo_enabled() )
      {
        if( stat == CUBIT_FAILURE )
          CubitUndo::remove_last_undo();
      }

      return stat;
    }

    if( merge )
      MergeTool::instance()->merge_bodies( bodies_to_merge );
  }
  else
  {
    DLIList<RefFace*> faces_not_to_merge;
    Body *body1 = face1->body();
    Body *body2 = face2->body();
    /*
    if(merge)
    {
      DLIList<RefFace*> tmp_faces;
      body1->ref_faces(tmp_faces);
      if(tmp_faces.move_to(face1))
        tmp_faces.extract();
      faces_not_to_merge = tmp_faces;

      tmp_faces.clean_out();
      body2->ref_faces(tmp_faces);
      if(tmp_faces.move_to(face2))
        tmp_faces.extract();
      faces_not_to_merge += tmp_faces;
    }
    */

    //get the modify engine
    DLIList<Surface*> surf_list( 1 );
    DLIList<RefFace*> ref_face_list( 2 );
    ref_face_list.append( face1 );
    ref_face_list.append( face2 );
    GeometryModifyEngine* gme = common_modify_engine( ref_face_list, surf_list);

    if ( !gme )
    {
      PRINT_ERROR("Performing IMPRINT with entities containing geometry from\n"
                  "different modeling engines is not allowed.\n" );
      return CUBIT_FAILURE;
    }

    if( CubitUndo::get_undo_enabled() )
      CubitUndo::save_state_with_cubit_file( ref_faces );

    //do the imprint onto face1
    DLIList<BodySM*> new_bodysm_list;
    CubitStatus status = gme->tolerant_imprint(surf_list, new_bodysm_list);

    //if we failed...get out
    if( status == CUBIT_FAILURE )
    {
      if( CubitUndo::get_undo_enabled() )
        CubitUndo::remove_last_undo();

      return CUBIT_FAILURE;
    }

    DLIList<Body*> result_bodies;
    DLIList<Body*> original_body_list;
    original_body_list.append(body1);
    original_body_list.append(body2);
    status = finish_sm_op( original_body_list, new_bodysm_list, result_bodies );

    /*
    if( merge )
    {
      DLIList<RefFace*> faces_to_merge, tmp_faces;
      body1->ref_faces(faces_to_merge);
      body2->ref_faces(tmp_faces);
      faces_to_merge -= faces_not_to_merge;
      if(faces_to_merge.size() > 1)
      {
        MergeTool::instance()->merge_reffaces(faces_to_merge);
      }
    }
    */
  }

  return CUBIT_SUCCESS;
}


CubitStatus GeometryModifyTool::tolerant_imprint( RefFace *ref_face,
                                                  DLIList<RefEdge*> &ref_edge_list,
                                                  Body *&new_body,
                                                  bool merge )
{
  int i;
  DLIList<Body*> bodies_to_merge;
  if( merge )
  {
    DLIList<RefEdge*> tmp_edges;
    ref_face->ref_edges( tmp_edges );

    for( i=tmp_edges.size(); i--; )
    {
      RefEdge *tmp_edge = tmp_edges.get_and_step();
      DLIList<Body*> body_list;
      tmp_edge->bodies( body_list );
      bodies_to_merge += body_list;
    }

    for( i=ref_edge_list.size(); i--; )
    {
      RefEdge *tmp_edge = ref_edge_list.get_and_step();
      DLIList<Body*> body_list;
      tmp_edge->bodies( body_list );
      bodies_to_merge += body_list;
    }

    bodies_to_merge.uniquify_unordered();
  }


  DLIList<Body*> original_body_list;
  ref_face->bodies( original_body_list );

  DLIList<Surface*> surf_list( 1 );
  DLIList<Curve*> curve_list(ref_edge_list.size());
  DLIList<RefFace*> ref_face_list( 1 );
  ref_face_list.append( ref_face );
  GeometryModifyEngine* gme = common_modify_engine( ref_face_list,
                                                    ref_edge_list,
                                                    surf_list,
                                                    curve_list,
                                                    true);
  if ( !gme )
  {
    PRINT_ERROR("Performing IMPRINT with entities containing geometry from\n"
                "different modeling engines is not allowed.\n" );
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<RefFace*> ref_faces(1);
    ref_faces.append( ref_face );
    CubitUndo::save_state_with_cubit_file( ref_faces );
  }

  DLIList<BodySM*> body_sm_list;
  for(i=original_body_list.size(); i--;)
    body_sm_list.append(original_body_list.get_and_step()->get_body_sm_ptr());

  int process_composites = 0;
  if(contains_composites(original_body_list))
    process_composites = 1;

  DLIList<TopologyBridge*> tb_list;
  if(process_composites)
  {
    // Turn certain attributes on.
    do_attribute_setup();
    // Push virtual attributes down to solid model topology before
    // doing the imprint.
    push_vg_attributes_before_modify(body_sm_list);
    DLIList<TopologyBridge*> tmp_tb_list;
    CAST_LIST(surf_list, tmp_tb_list, TopologyBridge);
    // Put "ORIGINAL" attributes on the bodies being imprinted and
    // the curves as these originally existed.
    push_named_attributes_to_curves_and_points(tmp_tb_list, "ORIGINAL");
    CAST_LIST(curve_list, tb_list, TopologyBridge);
    push_named_attributes_to_curves_and_points(tb_list, "ORIGINAL");
    push_named_attributes_to_curves_and_points(tb_list, "IMPRINTER");
  }

  CubitStatus status = CUBIT_FAILURE;
  DLIList<BodySM*> new_body_list;
  DLIList<TopologyBridge*> new_tbs, att_tbs;
  // The bridges doing the imprinting often get split during the process but
  // because of the way we are making copies, the IMPRINTER attribute doesn't
  // get propagated to them.  temporary_bridges will be filled in with any
  // additional IMPRINTER bridges we need to consider below when deciding whether to
  // keep composite attributes.
  DLIList<TopologyBridge*> temporary_bridges;
  for(i=surf_list.size(); i>0; i--)
  {
    Surface *cur_surf = surf_list.get_and_step();
    BodySM *new_body_sm = NULL;
    CubitStatus tmp_status = gme->tolerant_imprint_surface_with_curves(
                                                cur_surf, curve_list,
                                                temporary_bridges,
                                                new_body_sm,
                                                &new_tbs, &att_tbs);
    if(new_body_sm)
      new_body_list.append(new_body_sm);
    if(tmp_status == CUBIT_SUCCESS)
      status = tmp_status;
  }

  temporary_bridges.uniquify_ordered();

  if( status == CUBIT_FAILURE )
  {
    if( CubitUndo::get_undo_enabled() )
      CubitUndo::remove_last_undo();

    if(process_composites)
    {
      remove_pushed_attributes(new_body_list, original_body_list);
      do_attribute_cleanup();
    }
    return CUBIT_FAILURE;
  }
  else
  {
    if(process_composites)
    {
      // Analyze the results and adjust virtual attributes as necessary.
      DLIList<TopologyBridge*> tmp_tb_list;
      CAST_LIST(new_body_list, tmp_tb_list, TopologyBridge);
      tb_list.merge_unique(tmp_tb_list);
      // The bridges coming back in temporary_bridges may not have IMPRINTER
      // attributes on them becuase of the way they were generated below.  Make
      // sure they get IMPRINTER attributes.
      push_named_attributes_to_curves_and_points(temporary_bridges, "IMPRINTER");
      tb_list += temporary_bridges;
      GeometryQueryTool::instance()->ige_attribute_after_imprinting(new_tbs, att_tbs,
        tb_list, original_body_list);

      while(temporary_bridges.size())
        delete temporary_bridges.pop();

      // Clean up attributes.
      remove_imprint_attributes_after_modify(body_sm_list, new_body_list);

      // Restore the virtual geometry.
      restore_vg_after_modify(new_body_list, original_body_list, gme);
      remove_pushed_attributes(new_body_list, original_body_list);
    }
  }

  DLIList<Body*> result_bodies;
  status = finish_sm_op( original_body_list, new_body_list, result_bodies );

  if(process_composites)
    do_attribute_cleanup();

  if( status == CUBIT_FAILURE )
  {
    if( CubitUndo::get_undo_enabled() )
      CubitUndo::remove_last_undo();
    return CUBIT_FAILURE;
  }

  if( merge )
    MergeTool::instance()->merge_bodies( bodies_to_merge );

  if( result_bodies.size() == 1 )
    new_body = result_bodies.get();
  else
    return CUBIT_FAILURE;

  return status;
}

CubitStatus GeometryModifyTool::tolerant_imprint( DLIList<Body*> &bodies,
                                      DLIList<Body*> &new_bodies, bool merge )
{
  //make sure all bodies are from the same modify engine
  DLIList<BodySM*> body_sm_list;
  GeometryModifyEngine* gme = common_modify_engine(bodies, body_sm_list);
  if ( !gme )
  {
    PRINT_ERROR("Performing IMPRINT with volumes containing geometry\n"
                "from different modeling engines is not allowed.\n"
                "Delete uncommon geometry on these volumes before operation.\n\n");
    return CUBIT_FAILURE;
  }


  //make sure that merge tolerance is not inapproiate for model
  int i;
  CubitBox bounding_box( CubitVector(0,0,0),
                         CubitVector(CUBIT_DBL_MAX,
                                     CUBIT_DBL_MAX,
                                     CUBIT_DBL_MAX ) );
  for( i=bodies.size(); i--; )
  {
    CubitBox tmp_box = bodies.get_and_step()->bounding_box();
    if(bounding_box.max_x() == CUBIT_DBL_MAX)
      bounding_box = tmp_box;
    else if( tmp_box.diagonal().length_squared() <
        bounding_box.diagonal().length_squared() )
      bounding_box = tmp_box;
  }

  //get the merge tolerance
  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;

  //if the merge tolerance is greater than 1/10th the length of the
  //diagonal of the bounding box of the smallest volume, fail!
  double tenth_smallest_bbox = 0.1*(bounding_box.diagonal().length());
  if( tolerance > tenth_smallest_bbox )
  {
    PRINT_ERROR("Merge tolerance is set excessively high.  Must be lower than %f\n",
                 tenth_smallest_bbox );
    PRINT_INFO("       (Merge tolerance must be less than than 1/10th of the diagonal\n"
                        "of the bounding box of the smallest volume)\n");
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state_with_cubit_file( bodies );

  body_sm_list.clean_out();
  DLIList<Body*> old_body_list;
  old_body_list += bodies;
  bodies.reset();
  for( i=bodies.size(); i--; )
    body_sm_list.append( bodies.get_and_step()->get_body_sm_ptr() );

  int process_composites = 0;
  if(contains_composites(bodies))
    process_composites = 1;

  if(process_composites)
  {
    // Push virtual attributes down to solid model topology before
    // doing the imprint.
    do_attribute_setup();
    push_vg_attributes_before_modify(body_sm_list);
    // This must be done after pushing the vg atts because it uses them.
    DLIList<TopologyBridge*> tb_list;
    CAST_LIST(body_sm_list, tb_list, TopologyBridge);
    push_named_attributes_to_curves_and_points(tb_list, "IMPRINTER");
    push_named_attributes_to_curves_and_points(tb_list, "ORIGINAL");
  }

  DLIList<BodySM*> new_body_list;
  DLIList<TopologyBridge*> new_tbs, att_tbs;
  CubitStatus result = gme->tolerant_imprint( body_sm_list, new_body_list, &new_tbs, &att_tbs );


  if(result == CUBIT_FAILURE)
  {
    if(process_composites)
    {
      remove_pushed_attributes(new_body_list, bodies);
      do_attribute_cleanup();
    }
    return result;
  }
  else
  {
    if(process_composites)
    {
      // Analyze the results and adjust virtual attributes as necessary.
       DLIList<TopologyBridge*> tb_list;
       CAST_LIST(new_body_list, tb_list, TopologyBridge);
      GeometryQueryTool::instance()->ige_attribute_after_imprinting(new_tbs, att_tbs,
        tb_list, bodies);

      // Clean up attributes.
      remove_imprint_attributes_after_modify(body_sm_list, new_body_list);

      // Restore the virtual geometry.
      restore_vg_after_modify(new_body_list, bodies, gme);
      remove_pushed_attributes(new_body_list, bodies);
    }
  }

  result = finish_sm_op( bodies, body_sm_list, new_bodies );

  if(process_composites)
    do_attribute_cleanup();

  if(result == CUBIT_FAILURE)
  {
    if( CubitUndo::get_undo_enabled() )
      CubitUndo::remove_last_undo();
    return CUBIT_FAILURE;
  }

  if( merge )
    MergeTool::instance()->merge_bodies( bodies );

  return CUBIT_SUCCESS;
}

CubitStatus GeometryModifyTool::remove_curve_slivers( DLIList<Body*> &bodies,
                                              double lengthlimit )
{
  DLIList<BodySM*> body_sm_list;
  GeometryModifyEngine* gme = common_modify_engine(bodies, body_sm_list);
  if ( !gme )
  {
    PRINT_ERROR("Curve sliver removal only supported on ACIS geometry\n");
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state_with_cubit_file( bodies );

  CubitStatus status = CUBIT_FAILURE;
  DLIList<BodySM*> modified_bodies;
  int i;
  for( i=body_sm_list.size(); i--; )
  {
    BodySM *tmp_body_sm = body_sm_list.get_and_step();
    if( gme->remove_curve_slivers( tmp_body_sm, lengthlimit ) == CUBIT_SUCCESS )
    {
      modified_bodies.append( tmp_body_sm );
      status = CUBIT_SUCCESS;
    }
  }

  if( status == CUBIT_FAILURE )
  {
    PRINT_WARNING("Did not remove any sliver curves\n");
    if( CubitUndo::get_undo_enabled() )
      CubitUndo::remove_last_undo();
    return CUBIT_FAILURE;
  }

  DLIList<Body*> dummy_list;
  status = finish_sm_op( bodies, modified_bodies, dummy_list );
  if( status == CUBIT_FAILURE )
  {
    PRINT_WARNING("Did not remove any sliver curves\n");
    if( CubitUndo::get_undo_enabled() )
      CubitUndo::remove_last_undo();
    return CUBIT_SUCCESS;
  }

  return status;
}

void GeometryModifyTool::split_surface_with_narrow_region(RefFace *face,
                                                      DLIList<CubitVector> &split_pos1_list,
                                                      DLIList<CubitVector> &split_pos2_list)
{
  int k;
  if(split_pos1_list.size() > 0)
  {
    DLIList<DLIList<CubitVector*>*> vec_lists;
    DLIList<CubitVector*> pt_list;
    for(k=split_pos1_list.size(); k--;)
    {
      CubitVector split_pos1 = split_pos1_list.get_and_step();
      CubitVector split_pos2 = split_pos2_list.get_and_step();
      face->move_to_surface(split_pos1);
      face->move_to_surface(split_pos2);
      DLIList<CubitVector*> *vec_list = new DLIList<CubitVector*>;
      vec_list->append( new CubitVector(split_pos1));
      vec_list->append( new CubitVector(split_pos2));
      vec_lists.append( vec_list );
    }

    GeometryModifyTool::instance()->split_surface(face,
      pt_list, vec_lists );

    while( vec_lists.size() )
    {
      DLIList<CubitVector*> *vec_list = vec_lists.pop();
      while( vec_list->size() ) delete vec_list->pop();
      delete vec_list;
    }
    /*
    while( pt_list.size() )
      delete( pt_list.pop() );
      */
  }
}

void GeometryModifyTool::fixup_merged_entities( DLIList<int> &merged_surface_ids,
                                                DLIList<int> &merged_curve_ids ) const
{
  //use ids to find surviving merged entities
  DLIList<RefFace*> ref_face_list;
  DLIList<RefEdge*> ref_edge_list;

  int i;
  //see what merged survived operation
  for( i=merged_surface_ids.size(); i--; )
  {
    int face_id = merged_surface_ids.get_and_step();
    RefFace *surviving_merged_face = RefEntityFactory::instance()->get_ref_face( face_id );
    if( surviving_merged_face )
      ref_face_list.append( surviving_merged_face );
  }

  //see what merged survived operation
  for( i=merged_curve_ids.size(); i--; )
  {
    int edge_id = merged_curve_ids.get_and_step();
    RefEdge *surviving_merged_edge = RefEntityFactory::instance()->get_ref_edge( edge_id );
    if( surviving_merged_edge )
      ref_edge_list.append( surviving_merged_edge );
  }

  //fix up merged faces -- some might need to be reversed
  for(i=ref_face_list.size(); i--; )
  {
    RefFace *merged_face = ref_face_list.get_and_step();
    BasicTopologyEntity *bte = static_cast<BasicTopologyEntity*>(merged_face);

    //get the first bridge of the entity
    DLIList<TopologyBridge*> face_bridge_list;
    bte->bridge_manager()->get_bridge_list( face_bridge_list );

    //if there are 2 bridges in the list, it's still merged...do nothing
    if( face_bridge_list.size() > 1 )
      continue;

    //get the center of the RefFace
    CubitVector center = merged_face->center_point();

    //get the normal according to the RefFace
    CubitVector ref_face_normal = merged_face->normal_at( center );

    //get the normal at the center from the underlying Surface
    Surface *surface_ptr = CAST_TO( face_bridge_list.get(), Surface );
    CubitVector surface_normal;
    surface_ptr->closest_point( center, NULL, &surface_normal );

    //if normals are opposite, flip sense of surface_ptr
    if( fabs(ref_face_normal.interior_angle( surface_normal ) - 180 ) < 0.1  )
      merged_face->reverse_normal();

    //One more thing.....if surface is a composite, update the graphics
    //on the hidden curve...could have been deleted.
    if ( GeometryQueryTool::instance()->ige_is_composite( surface_ptr ) )
      merged_face->notify_all_observers( TOPOLOGY_MODIFIED );
  }

  //fix up merged edges -- some might need to be reversed
  for(i=ref_edge_list.size(); i--; )
  {
    RefEdge *merged_edge = ref_edge_list.get_and_step();
    BasicTopologyEntity *bte = static_cast<BasicTopologyEntity*>(merged_edge);

    //get the first bridge of the entity
    DLIList<TopologyBridge*> edge_bridge_list;
    bte->bridge_manager()->get_bridge_list( edge_bridge_list );

    //get start/end points of the edge
    CubitVector edge_start_point = merged_edge->start_vertex()->coordinates();
    CubitVector edge_end_point = merged_edge->end_vertex()->coordinates();

    //get start/end point of the curve
    edge_bridge_list.reset();
    Curve *curve_ptr = CAST_TO( edge_bridge_list.get(), Curve);
    DLIList<Point*> tmp_points;
    curve_ptr->points( tmp_points );
    CubitVector curve_start_point = tmp_points.get_and_step()->coordinates();
    CubitVector curve_end_point = tmp_points.get_and_step()->coordinates();

    //check to see if curve sense needs to be reversed
    if( edge_start_point.distance_between( curve_start_point ) < GEOMETRY_RESABS &&
          edge_end_point.distance_between( curve_end_point ) < GEOMETRY_RESABS )
    {
      //do nothing...everything is fine
      continue;
    }
    else
    {
      if( edge_start_point.distance_between( curve_end_point ) < GEOMETRY_RESABS &&
            edge_end_point.distance_between( curve_start_point ) < GEOMETRY_RESABS )
      {
        //switch sense of ref entity
        merged_edge->reverse_tangent();
      }
    }
  }
}

void GeometryModifyTool::get_merged_curve_and_surface_ids(
                                          DLIList<Body*> &bodies,
                                          DLIList<int> &merged_surface_ids,
                                          DLIList<int> &merged_curve_ids ) const
{
  int i;
  for( i=bodies.size(); i--; )
  {
    DLIList<RefEntity*> merged_children;

    MergeTool::instance()->contains_merged_children( bodies.get_and_step(),
                                                     merged_children );

    int j;
    for( j=merged_children.size(); j--; )
    {
      RefEntity *ref_ent = merged_children.get_and_step();

      RefFace *ref_face = CAST_TO( ref_ent, RefFace );

      if( ref_face )
        merged_surface_ids.append( ref_face->id() );
      else
      {
        RefEdge *ref_edge = CAST_TO( ref_ent, RefEdge );

        if( ref_edge )
          merged_curve_ids.append( ref_edge->id() );
      }
    }
  }
}

void GeometryModifyTool::plane_preview(DLIList<Body*>& body_list,
                                       const CubitVector &pt1,
                                       const CubitVector &pt2,
                                       const CubitVector &pt3)
{
  CubitPlane plane;
  if( plane.mk_plane_with_points( pt1, pt2, pt3) == CUBIT_FAILURE)
  {
       PRINT_INFO( "Unable to create plane from given information.\n" );
       return ;
  }

  CubitBox bounding_box;
  Body* body_ptr = body_list.get_and_step();
  bounding_box = body_ptr->bounding_box();

  int i;
  for( i=1; i<body_list.size(); i++ )
  {
     body_ptr = body_list.get_and_step();
     bounding_box |= body_ptr->bounding_box();
  }

  int extension_type = 1;
  double extension = 10; //10%
  CubitVector p1, p2, p3, p4;

  if( AnalyticGeometryTool::instance()->
         min_pln_box_int_corners( plane, bounding_box, extension_type,
         extension, p1, p2, p3, p4 ) == CUBIT_FAILURE )
  {
     PRINT_INFO( "Unable to create plane from given information.\n" );
     return ;
  }

  GPoint gp[4];
  gp[0].x=p1.x(); gp[0].y=p1.y(); gp[0].z=p1.z();
  gp[1].x=p2.x(); gp[1].y=p2.y(); gp[1].z=p2.z();
  gp[2].x=p3.x(); gp[2].y=p3.y(); gp[2].z=p3.z();
  gp[3].x=p4.x(); gp[3].y=p4.y(); gp[3].z=p4.z();

  // clear previous previews
  GfxPreview::clear();

  // Get the color to draw in
  int color = CUBIT_BLUE;
  GfxPreview::draw_quad(gp, color);
  GfxPreview::flush();
  return;
}

void GeometryModifyTool::march_path(CubitVector &start_pos,
                                    RefFace *start_face,
                                    CubitVector &march_dir,  // should be normalized
                                    double &step_size)
{
  double cos_45 = 0.70710678118654752440084436210485;
  double geo_tol = GeometryQueryTool::instance()->get_sme_resabs_tolerance();
  CubitVector point_on_surf = start_pos;
  RefVolume *v = start_face->ref_volume();
  RefFace *cur_face = start_face;
  CubitVector norm = cur_face->normal_at(point_on_surf, v); 
  CubitVector sweep_dir = norm;

  RefEdge *snap_edge = NULL;
  bool snapped_to_edge_last_time = false;
  bool turned = false;
  CubitVector old_pos = point_on_surf;
  while(!turned)
  {
    CubitVector new_pos;
    if(snapped_to_edge_last_time)
    {
      // Just set the new position to the position on the
      // edge.  This will force us to jump out without doing
      // anything and then on the next loop we will start 
      // onto the new face.
      new_pos = old_pos;
      snapped_to_edge_last_time = false;
      continue;
    }
    else
    {
      // Calculate a new step along the vector.
      new_pos = old_pos + step_size * march_dir;
    }

//GfxDebug::draw_point(new_pos, 5);
//GfxDebug::flush();

    cur_face->get_surface_ptr()->closest_point_trimmed(new_pos, point_on_surf);

//GfxDebug::draw_point(point_on_surf, 6);
//GfxDebug::flush();

    norm = cur_face->normal_at(point_on_surf, v); 
    if(sweep_dir % norm < cos_45)
    {
      turned = true;
      /*
GfxDebug::draw_point(old_pos, 3);
GfxDebug::draw_point(point_on_surf, 3);
GfxDebug::draw_line(old_pos, point_on_surf, 3);
GfxDebug::flush();
*/
    }
    else
    {
      bool snapping_to_edge = true;
      CubitVector proj_dir = point_on_surf - new_pos;
      double proj_dist = proj_dir.length();
      if(proj_dist < geo_tol)
        snapping_to_edge = false;
      else
      {
        proj_dir /= proj_dist;
        double dot = proj_dir % norm;
        if(dot > .99 || dot < -.99)
          snapping_to_edge = false;
      }
      if(!snapping_to_edge)
      {
        snap_edge = NULL;
GfxDebug::draw_point(old_pos, 3);
GfxDebug::draw_point(point_on_surf, 3);
GfxDebug::draw_line(old_pos, point_on_surf, 3);
GfxDebug::flush();
        // didn't snap to boundary
        old_pos = point_on_surf;
      }
      else
      {
        // probably snapped to the boundary
        DLIList<RefEdge*> face_edges;
        RefEdge *best_edge = NULL;
        cur_face->ref_edges(face_edges);
        int i;
        DLIList<RefEdge*> possible_edges;
        CubitVector closest;
        for(i=face_edges.size(); i>0; i--)
        {
          RefEdge *e = face_edges.get_and_step();
//GfxDebug::draw_ref_edge(e, 8);
//GfxDebug::flush();
          e->closest_point_trimmed(point_on_surf, closest);
//GfxDebug::draw_point(point_on_surf, 66);
//GfxDebug::draw_point(closest, 77);
//GfxDebug::flush();
          double cur_dist = (closest - point_on_surf).length();
          if(cur_dist < geo_tol)
          {
            possible_edges.append(e);
          }
        }
        if(possible_edges.size() == 1)
          best_edge = possible_edges.get();
        else if(possible_edges.size() > 1)
        {
          int h;
          double smallest_dist = CUBIT_DBL_MAX;
          for(h=possible_edges.size(); h>0; h--)
          {
            RefEdge *ce = possible_edges.get_and_step();
            ce->closest_point_trimmed(old_pos, closest);
            double cur_dist = (old_pos-closest).length();
            if(cur_dist < smallest_dist)
            {
              smallest_dist = cur_dist;
              best_edge = ce;
            }
          }
        }
        if(best_edge)
        {
          if(snap_edge && snap_edge == best_edge)
          {
GfxDebug::draw_point(old_pos, 3);
GfxDebug::draw_point(point_on_surf, 3);
GfxDebug::draw_line(old_pos, point_on_surf, 3);
GfxDebug::flush();
            old_pos = point_on_surf;
            snapped_to_edge_last_time = true;
          }
          else
          {
            snap_edge = best_edge;
            cur_face = best_edge->other_face(cur_face, v);
            CubitVector old_pos_save = old_pos;
            old_pos = closest;
            i = 0;
            snapped_to_edge_last_time = true;

            GeometryModifyEngine *gme = get_engine((TopologyBridge*)best_edge->get_curve_ptr());
            if(gme)
            {
              Point *pt1 = gme->make_Point(old_pos_save);
              Point *pt2 = gme->make_Point(new_pos);
              CubitVector const* pt3 = NULL;
              Curve *crv = gme->make_Curve(STRAIGHT_CURVE_TYPE, pt1, pt2, pt3,CUBIT_FORWARD);
              if(crv)
              {
                CubitVector pos1, pos2;
                double dist;
                GeometryQueryTool::instance()->entity_entity_distance(crv, best_edge->get_curve_ptr(), pos1,
                  pos2, dist);
                GfxDebug::draw_point(pos2, 9);
                GfxDebug::flush();
                old_pos = pos2;
                delete crv;
                delete pt1;
                delete pt2;
GfxDebug::draw_point(old_pos_save, 3);
GfxDebug::draw_point(pos2, 3);
GfxDebug::draw_line(old_pos_save, pos2, 3);
GfxDebug::flush();
              }
            }
          }
        }
      }
    }
  }

 // GfxDebug::draw_point(point_on_surf, 4);
 // GfxDebug::flush();
}

CubitStatus GeometryModifyTool::stitch( DLIList<Body*> &bodies_to_stitch,
                                        DLIList<Body*> &result_list,
                                        bool tighten_gaps,
                                        double tolerance )
{
  if (!okay_to_modify( bodies_to_stitch, "STITCH" ))
    return CUBIT_FAILURE;

  //get all the BodySMs from 'bodies_to_stitch'
  int i;
  for( i=bodies_to_stitch.size(); i--; )
  {
    Body *tmp_body = bodies_to_stitch.get_and_step();
    if( !tmp_body->is_sheet_body() )
    {
      PRINT_ERROR("Can't stitch body %d.  It's a solid body\n", tmp_body->id() ); 
      return CUBIT_FAILURE;
    }
  }
  
  DLIList<TopologyEntity*> entity_list(bodies_to_stitch.size());
  DLIList<TopologyBridge*> bridge_list(bodies_to_stitch.size());
  DLIList<BodySM*>         bodysm_list(bodies_to_stitch.size());
  CAST_LIST_TO_PARENT(bodies_to_stitch, entity_list);
  GeometryModifyEngine *gme = common_modify_engine(entity_list, bridge_list);
  
  if( entity_list.size() != bridge_list.size() )
  {
    PRINT_ERROR("Cannot stitch entities of different geometry engines.\n");
    return CUBIT_FAILURE;
  }

  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state_with_cubit_file( bodies_to_stitch );

  CAST_LIST(bridge_list, bodysm_list, BodySM);
  DLIList<BodySM*> new_bodies;
  CubitStatus result = gme->stitch( bodysm_list, new_bodies, tighten_gaps, tolerance );

  if( result == CUBIT_FAILURE && CubitUndo::get_undo_enabled() )
    CubitUndo::remove_last_undo();

  if (!finish_sm_op(bodies_to_stitch, new_bodies, result_list))
  {
    result = CUBIT_FAILURE;
    CubitUndo::remove_last_undo();
  }

  if( CubitUndo::get_undo_enabled() )
  {
    if( result == CUBIT_SUCCESS )
      CubitUndo::note_result_bodies( result_list );
    else
      CubitUndo::remove_last_undo();
  }

  return result;
}



CubitStatus GeometryModifyTool::discover_topology(RefFace *surf, CubitVector &pos,
                                                  double &step_size,
                                                  int num_subdivisions)
{
  CubitStatus ret = CUBIT_SUCCESS;

  CubitVector point_on_surf;
  surf->get_surface_ptr()->closest_point_trimmed(pos, point_on_surf);
  RefVolume *v = surf->ref_volume();
  CubitVector norm = surf->normal_at(point_on_surf, v); 
  CubitVector march_dir = norm * CubitVector(1,0,0);
  if(march_dir.length() < .001)
  {
    march_dir = norm * CubitVector(0,1,0);
    if(march_dir.length() < .001)
    {
      march_dir = norm * CubitVector(0,0,1);
      if(march_dir.length() < .001)
      {
        PRINT_ERROR("Couldn't find a good march direction.\n");
        ret = CUBIT_FAILURE;
      }
    }
  }

  if(ret == CUBIT_SUCCESS)
  {
    // Get initial 4 directions.
    march_dir.normalize();
    DLIList<CubitVector> march_directions;
    CubitVector v1 = march_dir;
    CubitVector v3 = -march_dir;
    CubitVector v2 = norm*march_dir;
    CubitVector v4 = -v2;
    march_directions.append(v1);
    march_directions.append(v2);
    march_directions.append(v3);
    march_directions.append(v4);

    // Now subdivide directions further if requested.
    int i;
    // make sure we start at the end to process the all of the original directions correctly
    march_directions.last();  
    for(i=march_directions.size(); i>0; i--)
    {
      CubitVector cur_dir = march_directions.get_and_step();
      CubitVector next_dir = march_directions.get();
      subdivide_pie(cur_dir, next_dir, num_subdivisions, march_directions);
    }

    for(i=march_directions.size(); i>0; i--)
    {
      CubitVector cur_dir = march_directions.get_and_step();
      march_path(point_on_surf, surf, cur_dir, step_size);
    }
  }

  return ret;
}

void GeometryModifyTool::subdivide_pie(CubitVector &dir1, CubitVector &dir2, int num_subdivisions,
                                       DLIList<CubitVector> &all_directions)
{
  if(num_subdivisions > 0)
  {
    CubitVector mid = dir1 + dir2;
    mid.normalize();
    all_directions.append(mid);
    if(num_subdivisions > 1)
    {
      subdivide_pie(dir1, mid, num_subdivisions-1, all_directions);
      subdivide_pie(mid, dir2, num_subdivisions-1, all_directions);
    }
  }
}

//-------------------------------------------------------------------------
// traverse the body object and calls premodify function on each structure
//-------------------------------------------------------------------------
void GeometryModifyTool::body_premodify(Body* body) const
{
	// skip if the mesh autodelete is disabled
	if (!get_mesh_autodelete()) { return; }

	// volumes
	DLIList<RefVolume*> temp_vols;
	body->ref_volumes(temp_vols);
	for (int v = 0; v < temp_vols.size(); v++) {
		RefVolume* volume = temp_vols.get_and_step();
		volume->geometry_premodify();

		// faces
		DLIList<RefFace*> temp_faces;
		volume->ref_faces(temp_faces);
		for (int f = 0; f < temp_faces.size(); f++) {
			RefFace* face_ptr = temp_faces.get_and_step();
			face_ptr->geometry_premodify();

			//edges
			DLIList<RefEdge*> temp_edges;
			face_ptr->ref_edges(temp_edges);
			for (int e = 0; e < temp_edges.size(); e++) {
				RefEdge* edge_ptr = temp_edges.get_and_step();
				edge_ptr->geometry_premodify();

				// vertices
				DLIList<RefVertex*> temp_vertices;
				edge_ptr->ref_vertices(temp_vertices);
				for (int vertices = 0; vertices < temp_vertices.size(); vertices++) {
					RefVertex* vertex_ptr = temp_vertices.get_and_step();
					vertex_ptr->geometry_premodify();
				}
			}
		}
	}
}


CubitStatus GeometryModifyTool::prepare_surface_sweep(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Surface*> &surfaces,
                              const CubitVector& sweep_vector,
                              bool sweep_perp,
                              bool through_all,
                              bool outward,
                              bool up_to_next,
                              Surface *stop_surf,
                              Curve *curve_to_sweep_along,
                              BodySM* &cutting_tool_ptr ,
                              const CubitVector* point,
                              double *angle)
{
  GeometryModifyEngine* gme = get_engine(blank_bodies.get());

  if(surfaces.size() == 0 )
    return CUBIT_FAILURE;

  DLIList<GeometryEntity*> ref_ent_list;
  Surface * temp_face = NULL;
  for(int i = 0; i < surfaces.size(); i++)
    {
      //copy the faces before sweep
      temp_face = gme->make_Surface(surfaces.get_and_step());
      if (temp_face)
        ref_ent_list.append((GeometryEntity*)temp_face);
    }

  BodySM* to_body = NULL;
  CubitStatus stat = CUBIT_SUCCESS;
  if(up_to_next && blank_bodies.size() > 1) //unite all bland_bodies
    {
       DLIList<BodySM*> newBodies;
       DLIList<BodySM*> copied_bodies;
       for(int i = 0; i < blank_bodies.size(); i++)
         copied_bodies.append(gme->copy_body(blank_bodies.get_and_step()));

       stat = gme->unite(copied_bodies, newBodies);
       if(stat == CUBIT_FAILURE)
         {
           PRINT_ERROR("Cannot use 'up_to_next' option with specified geometry\n");
           PRINT_INFO("Try the 'stop surface <id>' option instead\n");
           return stat;
         }
       to_body = newBodies.get();
    }

  else if(up_to_next && blank_bodies.size() == 1)
    to_body = gme->copy_body(blank_bodies.get());

  DLIList<BodySM*> swept_bodies;
  if (point && angle) //sweep_surface_rotated
    stat = gme->sweep_rotational(ref_ent_list,swept_bodies,*point,
                           sweep_vector, *angle,0, 0.0,0,false,false,
                           false,stop_surf, to_body);  
	
  else
  {
    CubitVector tmp_sweep_vector = sweep_vector;

    //get model bbox info...will scale sweep vector by its diagonal
    //so that we go far enough
    if( through_all || stop_surf || up_to_next)
    {
      CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
      tmp_sweep_vector.normalize();
      tmp_sweep_vector*=(2*bounding_box.diagonal().length());
    }
    
    //see if we're sweeping along a specified curve
    if( curve_to_sweep_along )
      {
        DLIList<Curve*> curves_to_sweep_along;
        curves_to_sweep_along.append(curve_to_sweep_along);
        stat = gme->sweep_along_curve(ref_ent_list, swept_bodies,
                                 curves_to_sweep_along, 0.0,0,false,stop_surf,
                                 to_body);
      }

    else if (sweep_perp )
      stat = gme->sweep_perpendicular(ref_ent_list, swept_bodies,
                                 tmp_sweep_vector.length(),0.0,0,!outward,false,
                                 stop_surf, to_body);
    else
      stat = gme->sweep_translational(ref_ent_list, swept_bodies,
                                 tmp_sweep_vector,0.0,0, false, false, stop_surf,
                                 to_body);
  }

  if(stat == CUBIT_FAILURE || swept_bodies.size() == 0)
    {
       //delete copied faces
       GeometryEntity * temp_entity = NULL;
       for(int i = ref_ent_list.size();i--;)
         {
            temp_entity = ref_ent_list.get_and_step();
            if (temp_entity)
              gme->get_gqe()->delete_solid_model_entities( (Surface*)temp_entity);
         }

       return stat;
    }

  //if there are more than 1, unite them all
  DLIList<BodySM*>  newBodies;
  if (swept_bodies.size() > 1)
     stat = gme->unite(swept_bodies, newBodies);
  else
     newBodies = swept_bodies;

  if(stat == CUBIT_FAILURE || newBodies.size()!= 1)
    {
       PRINT_ERROR("webcut tool body is not created from acis.\n");
       //delete the swept_bodies
       BodySM* tmp_body = NULL;
       for (int i = swept_bodies.size(); i--;)
         {
           tmp_body= swept_bodies.get_and_step();
           if (tmp_body)
             gme->get_gqe()->delete_solid_model_entities(tmp_body);
         }

       //delete copied faces
       GeometryEntity * temp_entity = NULL;
       for(int i = ref_ent_list.size();i--;)
         {
            temp_entity = ref_ent_list.get_and_step();
            if (temp_entity)
              gme->get_gqe()->delete_solid_model_entities( (Surface*)temp_entity);
         }
       return CUBIT_FAILURE;
    }
  
    cutting_tool_ptr = newBodies.get();
    return stat;
}
