//-------------------------------------------------------------------------
// Filename      : AcisModifyEngine.cpp
//
// Purpose       : Performs all acis-based geometry modification or creation
//
// Special Notes :
//
// Creator       : Tim Tautges
//
// Creation Date :
//
// Owner         : Tim Tautges
//-------------------------------------------------------------------------
// ********** BEGIN STANDARD INCLUDES         **********
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
// ********** END STANDARD INCLUDES           **********

// ********** BEGIN ACIS INCLUDES             **********
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/acis.hxx"
#include "kernel/kernutil/version/version.hxx"
#include "operator/kernapi/api/operapi.hxx"
#include "boolean/kernapi/api/boolapi.hxx"
#include "constrct/kernapi/api/cstrapi.hxx"
#include "constrct/sg_husk/utils/utils.hxx"
#include "cover/kernapi/api/coverapi.hxx"
#include "euler/kernapi/api/eulerapi.hxx"
#include "intersct/kernapi/api/intrapi.hxx"
#include "sbool/kernapi/api/sboolapi.hxx"
#include "kernel/kernutil/law/generic_graph.hxx"
#include "ct_husk/api/ctapi.hxx"
#include "ct_husk/classes/cell.hxx"
#include "boolean/sg_husk/sanity/err_ent.hxx"
#else
#include "acis.hxx"
#include "version.hxx"
#include "operapi.hxx"
#include "boolapi.hxx"
#include "sboolapi.hxx"
#include "cstrapi.hxx"
#include "coverapi.hxx"
#include "eulerapi.hxx"
#include "intrapi.hxx"
#include "stchapi.hxx"
#include "utils.hxx"
#include "generic_graph.hxx"
#include "cell.hxx"
#include "ctapi.hxx"
#include "err_ent.hxx"
#include "acistype.hxx"
#include "entity_simplify.hxx"
#include "esplit.hxx"
#endif

#if CUBIT_ACIS_VERSION >= 800
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernapi/api/acis_options.hxx"
#else
#include "acis_options.hxx"
#endif
#endif

#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernapi/api/api.hxx"
#include "kernel/kernapi/api/kernapi.hxx"
#include "sweep/kernapi/api/sweepapi.hxx"
#include "intersct/sg_husk/query/ptentrel.hxx"
#include "kernel/kerndata/geometry/getbox.hxx"
#include "kernel/kerndata/data/debug.hxx"
#include "kernel/kerndata/data/entity.hxx"
#include "kernel/kerndata/data/datamsc.hxx"
#include "kernel/kernapi/api/api.err"
#include "kernel/kerndata/attrib/attrib.hxx"
#include "kernel/kerndata/geom/pcurve.hxx"
#include "kernel/kerndata/geom/point.hxx"
#include "kernel/kerndata/geom/allcurve.hxx"
#include "kernel/kerndata/geom/allsurf.hxx"
#include "kernel/kerndata/geom/plane.hxx"
#include "kernel/kerndata/geom/cone.hxx"
#include "kernel/kerndata/geom/sphere.hxx"
#include "kernel/kerndata/geom/torus.hxx"
#include "kernel/kerndata/geom/spline.hxx"
#include "kernel/kerndata/geom/transfrm.hxx"
#include "kernel/kerndata/geom/surface.hxx"
#include "kernel/kerndata/geometry/geometry.hxx"
#include "intersct/kerndata/geometry/geomutil.hxx"
#include "kernel/kerndata/lists/lists.hxx"
#include "kernel/kerngeom/surface/surdef.hxx"
#include "kernel/kerngeom/surface/allsfdef.hxx"
#include "kernel/kerngeom/curve/curdef.hxx"
#include "constrct/kernwire/sweeping/sweeping.hxx"
#include "shl_husk/api/shl_api.hxx"
#else
#include "api.hxx"
#include "kernapi.hxx"
#include "sweepapi.hxx"
#include "ptentrel.hxx"
#include "getbox.hxx"
#include "debug.hxx"
#include "entity.hxx"
#include "datamsc.hxx"
#include "api.err"
#include "attrib.hxx"
#include "pcurve.hxx"
#include "point.hxx"
#include "allcurve.hxx"
#include "allsurf.hxx"
#include "plane.hxx"
#include "cone.hxx"
#include "sphere.hxx"
#include "torus.hxx"
#include "spline.hxx"
#include "transfrm.hxx"
#include "surface.hxx"
#include "geometry.hxx"
#include "geomutil.hxx"
#include "lists.hxx"
#include "surdef.hxx"
#include "vector_utils.hxx"
#include "allsfdef.hxx"
#include "curdef.hxx"
#include "sweeping.hxx"
#include "shl_api.hxx"
#include "warp_api.hxx"
#endif

#if CUBIT_ACIS_VERSION < 800
#include "sweep/kernwire/sweeping/piershee.hxx"
#endif


#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernutil/tensor/tensor.hxx"
#include "kernel/kerndata/top/alltop.hxx"
#include "kernel/kerndata/top/wire.hxx"
#include "kernel/kerndata/transent/transent.hxx"
#include "faceter/api/af_api.hxx"
#include "kernel/sg_husk/query/q_wire.hxx"
#include "intersct/sg_husk/query/sgquery.hxx"
#include "kernel/kernint/d3_chk/chk_stat.hxx"
#include "intersct/kernint/d3_chk/chk.hxx"
#include "kernel/sg_husk/api/sgapi.err"
#include "faceter/meshmgr/ppm.hxx"
#include "faceter/meshmgr/ppmface.hxx"
#include "kernel/kerndata/savres/fileinfo.hxx"
#include "kernel/spline/bs3_crv/sp3crtn.hxx"
#include "kernel/spline/api/spl_api.hxx"
#include "intersct/spline/bs3_crv/bs3cutil.hxx"
#include "kernel/spline/bs2_crv/sp2crtn.hxx"
#include "ga_husk/api/ga_api.hxx"
#include "kernel/sg_husk/sweep/swp_spl.hxx"
#include "kernel/geomhusk/copyent.hxx"
#include "kernel/kerndata/geom/cnstruct.hxx"
#include "skin/kernapi/api/skinapi.hxx"
#include "offset/kernapi/api/ofstapi.hxx"
#include "kernel/kernutil/law/law.hxx"
#include "kernel/geomhusk/wire_qry.hxx"
#include "kernel/kernint/intcucu/intcucu.hxx"
#include "kernel/geomhusk/entwray.hxx"
#include "constrct/geomhusk/wire_utl.hxx"
#include "intersct/sg_husk/query/sgquertn.hxx"
#include "kernel/kerndata/top/body.hxx"
#else
#include "tensor.hxx"
#include "alltop.hxx"
#include "wire.hxx"
#include "transent.hxx"
#include "af_api.hxx"
#include "q_wire.hxx"
#include "sgquery.hxx"
#include "chk_stat.hxx"
#include "chk.hxx"
//#include "sgapi.err"
#include "ppm.hxx"
#include "ppmface.hxx"
#include "fileinfo.hxx"
#include "sp3crtn.hxx"
#include "spl_api.hxx"
#include "bs3cutil.hxx"
#include "sp2crtn.hxx"
#include "ga_api.hxx"
#include "swp_spl.hxx"
#include "copyent.hxx"
#include "cnstruct.hxx"
#include "skinapi.hxx"
#include "ofstapi.hxx"
#include "law.hxx"
#include "wire_qry.hxx"
#include "intcucu.hxx"
#include "entwray.hxx"
#include "wire_utl.hxx"
#include "sgquertn.hxx"
#include "body.hxx"
#endif

#ifdef ACIS_HEALER
#if CUBIT_ACIS_VERSION < 1100
#include "healhusk/heal_api/heal_api.hxx"
#else
#include "heal_api.hxx"
#endif
#endif

#if CUBIT_ACIS_VERSION < 1100
#include "kernel/geomhusk/getowner.hxx"
#include "baseutil/vector/vector.hxx"
#include "baseutil/vector/unitvec.hxx"
#include "baseutil/vector/position.hxx"
#include "baseutil/vector/transf.hxx"
#include "baseutil/vector/box.hxx"
#include "baseutil/errorsys/errmsg.hxx"
#include "baseutil/logical.h"
#include "baseutil/debug/module.hxx"
#include "baseutil/option/option.hxx"
#include "sweep/sg_husk/sweep/swp_opts.hxx"
#include "kernel/geomhusk/getowner.hxx"
#else
#include "getowner.hxx"
#include "vector.hxx"
#include "unitvec.hxx"
#include "position.hxx"
#include "transf.hxx"
#include "box.hxx"
#include "errmsg.hxx"
#include "logical.h"
#include "module.hxx"
#include "option.hxx"
#include "swp_opts.hxx"
#include "getowner.hxx"
#endif

#ifdef ACIS_IGES_TRANSLATOR
#if CUBIT_ACIS_VERSION < 1100
   #include "igeshusk/igs_api/routines.hxx"
#else
#include "acisiges_api.hxx"
#endif
#endif
#ifdef ACIS_PROE_TRANSLATOR
   #include "proe/proehusk/api/proeapi.hxx"
#endif

#ifdef ACIS_LOCAL_OPS
#if CUBIT_ACIS_VERSION < 1100
#include "rem_husk/api/rem_api.hxx"
#else
#include "rem_api.hxx"
#endif
#endif

#ifdef ACIS_STEP_TRANSLATOR
#if CUBIT_ACIS_VERSION < 1100
   #include "stephusk/include/acis_inc.hxx"
   #include "stephusk/include/stpapi.hxx"
   #include "stephusk/util/apiinit.hxx"
   #include "stephusk/stepwrit/api/apiwrite.hxx"
   #include "stephusk/stepread/api/apiread.hxx"
   #include "stephusk/util/stphead.hxx"
   #include "stephusk/util/stp_res.hxx"
#else
#include "acisstep_api.hxx"
#endif
#endif

// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********
#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "CubitUtil.hpp"
#include "GeometryDefines.h"
#include "CubitPlane.hpp"
#include "CpuTimer.hpp"
#include "ProgressTool.hpp"
#include "AppUtil.hpp"

#include "GeometryEntity.hpp"
#include "LumpACIS.hpp"
#include "SurfaceACIS.hpp"
#include "CurveACIS.hpp"
#include "PointACIS.hpp"

#include "BodyACIS.hpp"
#include "ShellACIS.hpp"
#include "LoopACIS.hpp"
#include "CoEdgeACIS.hpp"

#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"

#include "AcisQueryEngine.hpp"
#include "AcisModifyEngine.hpp"
#include "AcisHealerTool.hpp"
#include "AcisTweakTool.hpp"
#include "AcisSurfaceTool.hpp"
#include "AcisFacetManager.hpp"
#include "AcisEdgeTool.hpp"

#include "DoubleListItem.hpp"

#include "attrib_cubit_owner.hpp"

#include "CubitVector.hpp"

#include "AnalyticGeometryTool.hpp"
#include "MergeTool.hpp"
#include "SurfaceOverlapTool.hpp"
#include "GeometryUtil.hpp"
#include "RefEntityFactory.hpp"

#include "DLIList.hpp"

#include "CastTo.hpp"
#include "DAG.hpp"

#include "GMem.hpp"

#include "Body.hpp"
#include "attrib_snl_simple.hpp"

// still need this for vertex-id-rearranging junk
#include "RefVertex.hpp"

// include the DagDrawingTool header FOR DEBUGGING PURPOSES
#include "DagDrawingTool.hpp"

// including the AcisDrawTool for debugging
#include "AcisDrawTool.hpp"
#include "GfxDebug.hpp"

AcisModifyEngine* AcisModifyEngine::instance_ = 0;

AcisModifyEngine::~AcisModifyEngine()
{
   api_terminate_booleans();
   api_terminate_covering();
   api_terminate_sweeping();
   api_terminate_skinning();
   api_terminate_offsetting();
   api_terminate_shelling ();
#ifdef ACIS_LOCAL_OPS
   api_terminate_face_removal();
#endif

   instance_ = 0;
}

AcisModifyEngine::AcisModifyEngine()
{
  assert( !instance_ );

    // add this modify engine to geometrymodifytool
  GeometryModifyTool::instance()->add_gme(this);

   API_BEGIN;
     //initialize the correct modules
   api_initialize_booleans();
   api_initialize_covering();
   api_initialize_sweeping();
   api_initialize_skinning();
   api_initialize_offsetting();
   api_initialize_shelling ();
#ifdef ACIS_LOCAL_OPS
   api_initialize_face_removal();
#endif

   API_END;
}

CubitString AcisModifyEngine::identify()
{
   CubitString number_string = get_major_version();
   number_string += ".";
   number_string += CubitString(get_minor_version());
   number_string += ".";
   number_string += CubitString(get_point_version());
   CubitString version_string = "ACIS Version ";
   version_string += number_string;
   return version_string;
}

//-------------------------------------------------------------------------
// Purpose       : Return the type of this class.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/18/96
//-------------------------------------------------------------------------

const type_info& AcisModifyEngine::entity_type_info() const
{
   return typeid(AcisModifyEngine);
}

//-------------------------------------------------------------------------
// Purpose       : Create a sphere and return a pointer to the Body.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/18/96
//-------------------------------------------------------------------------

BodySM* AcisModifyEngine::sphere(double radius) const
{
   BODY* BODYPtr = NULL;

     // Create a BODY that represents the cuboid
   outcome result = api_make_sphere ( radius, BODYPtr );
   if (!result.ok())
   {
      PRINT_ERROR("In AcisModifyEngine::sphere, Line# %d\n"
                  "       Problems creating a sphere\n", __LINE__ );
      return NULL ;
   }

     // Build a Body from the BODY and return it.
   BodySM *this_bodysm = AcisQueryEngine::instance()->populate_topology_bridges(BODYPtr);
   return this_bodysm;
}

//-------------------------------------------------------------------------
// Purpose       : Create a brick and return a pointer to the Body.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/14/96
//-------------------------------------------------------------------------
BodySM* AcisModifyEngine::brick(double width, double depth,
                                double height) const
{
   BODY* BODYPtr = NULL;

     // Create a BODY that represents the cuboid
   outcome result = api_make_cuboid ( width, depth, height, BODYPtr );
   if (!result.ok())
   {
      PRINT_ERROR("In AcisModifyEngine::brick, Line# %d\n"
                  "       Problems creating a brick\n", __LINE__);
      return NULL ;
   }

   BodySM *this_bodysm = AcisQueryEngine::instance()->populate_topology_bridges(BODYPtr);
   return this_bodysm;
}

BodySM* AcisModifyEngine::brick( const CubitVector &center,
                                 const CubitVector axes[3],
                                 const CubitVector &extension ) const
{
   BODY* BODY_ptr;

   BODY_ptr = make_brick_BODY( center, axes, extension );

   BodySM *this_bodysm = AcisQueryEngine::instance()->populate_topology_bridges( BODY_ptr );
   return this_bodysm;
}

//-------------------------------------------------------------------------
// Purpose       : Create a prism and return a pointer to the Body.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/23/96
//-------------------------------------------------------------------------

BodySM* AcisModifyEngine::prism( double height, int sides,
                                 double major, double minor) const
{
   BODY* BODYPtr = NULL;

     // Create a BODY that represents the prism
   outcome result = api_make_prism ( height, major, minor, sides, BODYPtr );
   if (!result.ok())
   {
      PRINT_ERROR("In AcisModifyEngine::prism, Line# %d\n"
                  "       Problems creating a prism\n", __LINE__);
      return NULL ;
   }

     // Build a Body from the BODY and return it.
   BodySM *this_bodysm = AcisQueryEngine::instance()->populate_topology_bridges(BODYPtr);
   return this_bodysm;
}

//-------------------------------------------------------------------------
// Purpose       : Create a pyramid and return a pointer to the Body.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/23/96
//-------------------------------------------------------------------------

BodySM* AcisModifyEngine::pyramid( double height, int sides, double major,
                                   double minor, double top ) const
{
   BODY* BODYPtr = NULL;
     // Create a BODY that represents the pyramid
   outcome result = api_make_pyramid ( height, major, minor, top,
                                       sides, BODYPtr );
   if (!result.ok())
   {
      PRINT_ERROR("In AcisModifyEngine::pyramid, Line# %d\n"
                  "       Problems creating a pyramid\n", __LINE__);
      return NULL ;
   }
     // Build a Body from the BODY and return it.
   BodySM *this_bodysm = AcisQueryEngine::instance()->populate_topology_bridges(BODYPtr);
   return this_bodysm;
}

//-------------------------------------------------------------------------
// Purpose       : Create a cylinder and return a pointer to the Body.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/23/96
//-------------------------------------------------------------------------

BodySM* AcisModifyEngine::cylinder( double hi, double r1, double r2,
                                    double r3 ) const
{
   BODY* BODYPtr = NULL;
     // Create a BODY that represents the cylinder
   outcome result = api_make_frustum ( hi, r1, r2, r3, BODYPtr );
   if (!result.ok())
   {
      PRINT_ERROR("In AcisModifyEngine::cylinder, Line# %d\n"
                  "       Problems creating a cylinder\n", __LINE__);
      return NULL ;
   }

     // Build a Body from the BODY and return it.
   BodySM *this_bodysm = AcisQueryEngine::instance()->populate_topology_bridges(BODYPtr);
   return this_bodysm;
}

//-------------------------------------------------------------------------
// Purpose       : Create a torus and return a pointer to the Body.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/23/96
//-------------------------------------------------------------------------

BodySM* AcisModifyEngine::torus( double r1, double r2 ) const
{
   BODY* BODYPtr = NULL;
     // Create a BODY that represents the torus
   outcome result = api_make_torus ( r1, r2, BODYPtr );
   if (!result.ok())
   {
      PRINT_ERROR("In AcisModifyEngine::torus, Line# %d\n"
                  "       Problems creating a torus\n", __LINE__);
      return NULL ;
   }
     // Build a Body from the BODY and return it.
   BodySM *this_bodysm = AcisQueryEngine::instance()->populate_topology_bridges(BODYPtr);
   return this_bodysm;
}

//-------------------------------------------------------------------------
// Purpose       : Create a planar sheet body and return a pointer to the Body.
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 10/19/98
//-------------------------------------------------------------------------
BodySM* AcisModifyEngine::planar_sheet( const CubitVector& p1,
                                        const CubitVector& p2,
                                        const CubitVector& p3,
                                        const CubitVector& p4 ) const
{
   BODY* new_BODY = make_planar_quad_BODY( p1, p2, p3, p4 );

   // Build a VGI Body (i.e, the entire VGI structure) from this
   // ACIS BODY.
   BodySM *this_bodysm = AcisQueryEngine::instance()->populate_topology_bridges(new_BODY);

   return this_bodysm;
}

//-------------------------------------------------------------------------
// Purpose       : Return a copy of the input Body. The underlying
//                 ACIS BODY is also copied.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/24/96
//-------------------------------------------------------------------------

BodySM* AcisModifyEngine::copy_body ( BodySM* OSMEPtr) const
{
     // Make sure that we have a valid pointer
   assert(OSMEPtr != NULL) ;

   BodyACIS* bodyACISPtr = CAST_TO(OSMEPtr, BodyACIS) ;

     // Make sure that we have a valid BodyACIS
   assert(bodyACISPtr != NULL) ;

     // Get the ACIS BODY
   BODY* BODY_ptr = bodyACISPtr->get_BODY_ptr() ;

     // Construct a new BODY from what we have
   BODY* new_BODY_ptr = this->copy_BODY(BODY_ptr) ;

     // Build a VGI Body and return it.
   BodySM *this_bodysm = AcisQueryEngine::instance()->populate_topology_bridges(new_BODY_ptr);
   return this_bodysm;
}

CubitStatus AcisModifyEngine::imprint_BODYs (BODY* body1_ptr,
                                             BODY* body2_ptr) const
{
   CpuTimer imprint_BODYs_timer;

     // Imprints the 2 input ACIS BODYs. If something goes wrong, state is
     // restored and CUBIT_FALSE is returned.

   CubitStatus status = CUBIT_FAILURE;

     // Use the BEGIN/END pair so that a roll-back is done by ACIS if something
     // goes wrong in the api_imprint or api_check_entity operations
   if ( GeometryModifyTool::get_all_edges_imprint() )
   {
      api_set_int_option("all_free_edges", TRUE );
   }

   outcome result = api_imprint(body1_ptr, body2_ptr);

   //Gets rid of sliver curves/surfaces that could get produced by the imprint.
   AcisModifyEngine::instance()->cleanup_slivers( body1_ptr );
   AcisModifyEngine::instance()->cleanup_slivers( body2_ptr );

   if ( GeometryModifyTool::get_all_edges_imprint() )
   {
      api_set_int_option("all_free_edges", FALSE );
   }

   if (DEBUG_FLAG(18))
   {
      if (result.ok())
      {
         status = CUBIT_SUCCESS;

           // Now that the imprint operation has apparently succeeded (:-)
           // test the resulting BODYs for correctness. The output of the
           // api_check_entity procedure is sent to the file whose pointer
           // is debug_file_ptr (an ACIS global file pointer).

           // MJP Note:
           // I do this test here because there are bugs in the imprint
           // operation that result in invalid BODYs even though the api_imprint
           // call may return a successful status. I've been told that the imprint
           // operation will be a lot more stable in Version 2 -- and where have
           // we heard this kind of stuff before?? :-)

         FILE* save_dfp = NULL;

           // Save the original value of the ACIS Debug File Pointer
         save_dfp = debug_file_ptr;

           // Check the first entity and save the results of the check in the file,
           // api_check_entity_output.txt
         debug_file_ptr = fopen ("api_check_entity_output.txt", "w");
         result = api_check_entity (body1_ptr, (ENTITY_LIST *)NULL);
         if (debug_file_ptr) fclose (debug_file_ptr);

         if (result.ok())
         {
              // Check the second entity and save the results of the check in the
              // file, api_check_entity_output.txt
            debug_file_ptr = fopen("api_check_entity_output.txt", "w");
            result = api_check_entity (body2_ptr, (ENTITY_LIST *)NULL);
            if (debug_file_ptr) fclose (debug_file_ptr);

            if (result.ok())
            {
               status = CUBIT_SUCCESS;
            }

            else
            {
               status = CUBIT_FAILURE;
               AcisQueryEngine::instance()->ACIS_API_error (result);
            }
         }
         else
         {
            status = CUBIT_FAILURE;
            AcisQueryEngine::instance()->ACIS_API_error (result);
         }
           // Restore the original value of the ACIS Debug File Pointer
         debug_file_ptr = save_dfp;
      }
      else
      {
         status = CUBIT_FAILURE;
         AcisQueryEngine::instance()->ACIS_API_error (result);
      }
   }
   if (result.ok())
       status = CUBIT_SUCCESS;
   else
   {
       // mieht just be no overlap, check smaller boxes
     SPAbox box1 = get_body_box(body1_ptr);
     SPAbox box2 = get_body_box(body2_ptr);

     if (!body1_ptr || !body2_ptr || (box1 && box2))
     {
       if (DEBUG_FLAG(95))
       {
         PRINT_DEBUG_95( "Bodies were not "
                         "imprinted.\n");
         AcisQueryEngine::instance()->ACIS_API_error(result, "imprint Bodies");
       }
       status = CUBIT_FAILURE;
     }
   }
   PRINT_DEBUG_3( "CPU time taken to imprint two BODYs: "
               "%f secs\n", imprint_BODYs_timer.cpu_secs() ) ;

   return status;
}

//-------------------------------------------------------------------------
// Purpose       : Return a copy of the input ACIS BODY.  If so
//                 instructed, remove all the attributes attached to
//                 the original BODY in its copy.
//
// Special Notes :
//
// Creator       :
//
// Creation Date :
//-------------------------------------------------------------------------
BODY* AcisModifyEngine::copy_BODY (BODY* body_ptr,
                                   bool remove_attribs ) const
{
     // If the input pointer is NULL, return a NULL pointer
   if (body_ptr == NULL)
   {
      PRINT_ERROR("Cannot copy ACIS BODY - input pointer is NULL.\n");
      return (BODY *) NULL;
   }

   BODY* new_body_ptr = NULL;

   if( remove_attribs )
     ATTRIB_CUBIT_OWNER::set_copyable( false );

   outcome result = api_copy_body(body_ptr, new_body_ptr);

   if (!result.ok())
   {
      AcisQueryEngine::instance()->ACIS_API_error (result);
      new_body_ptr = NULL;
   }
     // Before returning this new BODY, first make sure that all the CUBIT
     // OWNER attributes are cleaned out of it
   if( remove_attribs )
     ATTRIB_CUBIT_OWNER::set_copyable( true );
/*
   if ( remove_attribs )
   {
     AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(new_body_ptr);
   }
*/
   // Also make sure that any healing attributes are removed from it
#ifdef ACIS_HEALER
   AcisHealerTool::instance()->end_BODY_for_healing( new_body_ptr );
#endif
   return new_body_ptr;
}

CubitStatus AcisModifyEngine::surface_intersection( Surface *surface1_ptr,
                                                    Surface *surface2_ptr,
                                                    DLIList<Curve*> &intersect_graph,
                                                    const double ) const
{
      //First make sure the entities are ACIS.
  SurfaceACIS *surf_ACIS_1= CAST_TO(surface1_ptr,SurfaceACIS);
  SurfaceACIS *surf_ACIS_2= CAST_TO(surface2_ptr,SurfaceACIS);
  if ( !surf_ACIS_1 || !surf_ACIS_2 )
  {
    PRINT_ERROR("In AcisModifyEngine::surface_intersection\n"
                "Not all entities are from the same engine.\n");
    return CUBIT_FAILURE;
  }
    //Get the original bodies.
  FACE* FACE1_ptr =  surf_ACIS_1->get_FACE_ptr();
  FACE* FACE2_ptr =  surf_ACIS_2->get_FACE_ptr();
  BODY* intersect_BODY = NULL;

    //Now get the intersection.
  CubitBoolean intersecting = CUBIT_FALSE;
  CubitStatus ret = FACE_intersection(FACE1_ptr, FACE2_ptr,
                                      intersect_BODY, intersecting);
  if ( ret != CUBIT_SUCCESS )
    return ret;
  if ( !intersecting )
  {
      //No intersection, we should have caught this before, but
      //just return with an empty edge list.
    return CUBIT_SUCCESS;
  }
  ret = convert_WIRE( intersect_BODY, intersect_graph );
  if ( ret != CUBIT_SUCCESS )
    return ret;

    //Okay we should have our fill of curves.
    //Delete the intersection graph wire, which should
    //be fine since we copied the curves.
  api_del_entity(intersect_BODY);

    //return sucess!
  return CUBIT_SUCCESS;
}

CubitStatus AcisModifyEngine::get_mid_plane( const CubitVector &point_1,
                                             const CubitVector &point_2,
                                             const CubitVector &point_3,
                                             BodySM *body_to_trim_to,
                                             BodySM *&midplane_body ) const
{

  DLIList<BodySM*> body_list;
  body_list.append(body_to_trim_to);
  SPAbox super_box;
  CubitStatus rv = AcisQueryEngine::instance()->create_super_acis_bounding_box( body_list, super_box );
  if (rv != CUBIT_SUCCESS)
    return CUBIT_FAILURE;

  BODY *sheet_BODY = this->create_infinite_plane_cutting_tool(
    point_1, point_2, point_3, super_box, CUBIT_TRUE);

  BODY* BODY_to_trim_to = AcisQueryEngine::get_BODY(body_to_trim_to);
  bool delete_attribs = CUBIT_TRUE;
  BODY* BODY_copy = this->copy_BODY(BODY_to_trim_to, delete_attribs);
  if (BODY_copy == NULL)
  {
    PRINT_ERROR("Cannot make a copy of the parent ACIS BODY of volume.\n"
                "       mid plane creation failed.\n" );
    return CUBIT_FAILURE;
  }

  outcome result = api_intersect( BODY_copy, sheet_BODY );

  if ( !result.ok() )
  {
    api_del_entity(BODY_copy);
    api_del_entity(sheet_BODY);
    PRINT_ERROR("Problems creating midplane.\n");
    return CUBIT_FAILURE;
  }

  //Create a new cubit body for the imprinted body.
  midplane_body = AcisQueryEngine::instance()->
    populate_topology_bridges(sheet_BODY);

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : get_spheric_mid_surface
// Member Type: PUBLIC
// Description: Calculates a mid-surface between 2 spheric surfaces.
// Author     : Philippe Pebay
// Date       : 03/07/06
//=============================================================================
CubitStatus AcisModifyEngine::get_spheric_mid_surface( Surface *surface_ptr1,
				  Surface *surface_ptr2,
				  BodySM *body_to_trim_to,
				  BodySM *&midsurface_body ) const
{
  SurfaceACIS *surf_ACIS1 = CAST_TO(surface_ptr1, SurfaceACIS );
  FACE* FACE1_ptr = make_FACE( surf_ACIS1->get_FACE_ptr() ); 
  int surface_type = (&(FACE1_ptr->geometry()->equation()))->type();
  if( surface_type != sphere_type )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_spheric_mid_surface\n"
		 "       The CGM sphere is not an ASIS sphere (%d %d).\n", surface_type, SPHERE_SURFACE_TYPE );
    return CUBIT_FAILURE;
  }
  SurfaceACIS *surf_ACIS2 = CAST_TO(surface_ptr2, SurfaceACIS );
  FACE* FACE2_ptr = make_FACE( surf_ACIS2->get_FACE_ptr() ); 
  surface_type = (&(FACE2_ptr->geometry()->equation()))->type();
  if( surface_type != sphere_type )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_spheric_mid_surface\n"
		 "       The CGM sphere is not an ASIS sphere (%d %d).\n", surface_type, SPHERE_SURFACE_TYPE );
    return CUBIT_FAILURE;
  }

//   surface const* acis_surface1 = &(FACE1_ptr->geometry()->equation());
//   sphere const* sphere_surface1 = (sphere *)acis_surface1;

  SURFACE const* acis_SURFACE1 = FACE1_ptr->geometry();
  SPHERE const* SPHERE_SURFACE1 = (SPHERE *)acis_SURFACE1;
  SPAposition center1 = SPHERE_SURFACE1->centre();

//   surface const* acis_surface2 = &(FACE2_ptr->geometry()->equation());
//   sphere const* sphere_surface2 = (sphere *)acis_surface2;
  SURFACE const* acis_SURFACE2 = FACE2_ptr->geometry();
  SPHERE const* SPHERE_SURFACE2 = (SPHERE *)acis_SURFACE2;
  SPAposition center2 = SPHERE_SURFACE2->centre();

  if ( center1 != center2 )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_spheric_mid_surface\n"
 		 "       Not supported yet: spheres do not have the same center.\n" );
    return CUBIT_FAILURE;
  }
  
  CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
  
  BODY *sheet_BODY;
  double radius = ( SPHERE_SURFACE1->radius() + SPHERE_SURFACE2->radius() ) / 2.;
  SPAunit_vector lat( 1., 0., 0. );
  SPAunit_vector lon( 0., 0., 1. );
  double const slat = - CUBIT_PI;
  double const elat = CUBIT_PI;
  double const slon = 0.;
  double const elon = 2. * CUBIT_PI;
  FACE *new_FACE_ptr = NULL;
  outcome result = api_make_spface( center1, 
				    radius, 
				    lat, 
				    lon,
				    slat, 
				    elat, 
				    slon, 
				    elon, 
				    new_FACE_ptr );
  if (!result.ok() || new_FACE_ptr == NULL ||
      new_FACE_ptr->geometry() == NULL )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_spheric_mid_surface\n"
		 "       ACIS api_make_spface function failed.\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    return CUBIT_FAILURE;
  }
  FACE *face_list[1];
  face_list[0] = new_FACE_ptr;
  result = api_sheet_from_ff( 1, face_list, sheet_BODY );
  if ( !result.ok() || sheet_BODY == NULL || sheet_BODY->lump() == NULL
      || sheet_BODY->lump()->shell() == NULL ||
      sheet_BODY->lump()->shell()->first_face() == NULL )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_spheric_mid_surface\n"
		 "       ACIS api_sheet_from_ff function failed.\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    return CUBIT_FAILURE;
  }
  result = api_body_to_2d( sheet_BODY );
  if ( !result.ok() )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_spheric_mid_surface\n"
		 "       ACIS api_body_to_2d function failed.\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    return CUBIT_FAILURE;
  } 

  BODY* BODY_to_trim_to = AcisQueryEngine::get_BODY(body_to_trim_to);
  bool delete_attribs = CUBIT_TRUE;
  BODY* BODY_copy = this->copy_BODY(BODY_to_trim_to, delete_attribs);
  if ( BODY_copy == NULL )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_spheric_mid_surface\n"
		 "       Cannot make a copy of the parent ACIS BODY of volume.\n" );
    return CUBIT_FAILURE;
  }

  result = api_intersect( BODY_copy, sheet_BODY );
  if ( !result.ok() )
  {
    api_del_entity( BODY_copy );    
    api_del_entity( sheet_BODY );    
    PRINT_ERROR( "In AcisModifyEngine::get_spheric_mid_surface\n"
		 "       ACIS api_intersect function failed.\n" );
    return CUBIT_FAILURE;
  }

  midsurface_body = AcisQueryEngine::instance()->populate_topology_bridges(sheet_BODY);

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : get_conic_mid_surface
// Member Type: PUBLIC
// Description: Calculates a mid-surface between 2 conic surfaces.
// Author     : Philippe Pebay
// Date       : 03/07/06
//=============================================================================
CubitStatus AcisModifyEngine::get_conic_mid_surface( Surface *surface_ptr1,
				  Surface *surface_ptr2,
				  BodySM *body_to_trim_to,
				  BodySM *&midsurface_body ) const
{
  SurfaceACIS *surf_ACIS1 = CAST_TO(surface_ptr1, SurfaceACIS );
  FACE* FACE1_ptr = make_FACE( surf_ACIS1->get_FACE_ptr() ); 
  int surface_type = (&(FACE1_ptr->geometry()->equation()))->type();
  if( surface_type != cone_type )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_conic_mid_surface\n"
		 "       The CGM cone is not an ASIS cone (%d %d).\n", surface_type, CONE_SURFACE_TYPE );
    return CUBIT_FAILURE;
  }
  SurfaceACIS *surf_ACIS2 = CAST_TO(surface_ptr2, SurfaceACIS );
  FACE* FACE2_ptr = make_FACE( surf_ACIS2->get_FACE_ptr() ); 
  surface_type = (&(FACE2_ptr->geometry()->equation()))->type();
  if( surface_type != cone_type )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_conic_mid_surface\n"
		 "       The CGM cone is not an ASIS cone (%d %d).\n", surface_type, CONE_SURFACE_TYPE );
    return CUBIT_FAILURE;
  }

  surface const* acis_surface1 = &(FACE1_ptr->geometry()->equation());
  cone const* cone_surface1 = (cone *)acis_surface1;
  double sine_angle1 = cone_surface1->sine_angle;
  double cosine_angle1 = cone_surface1->cosine_angle;

  surface const* acis_surface2 = &(FACE2_ptr->geometry()->equation());
  cone const* cone_surface2 = (cone *)acis_surface2;
  double sine_angle2 = cone_surface2->sine_angle;
  double cosine_angle2 = cone_surface2->cosine_angle;

  if ( ( sine_angle1 != 0.0 ) || ( fabs( cosine_angle1 ) != 1.0 ) 
       ||
       ( sine_angle2 != 0.0 ) || ( fabs( cosine_angle2 ) != 1.0 ) )
  {
    if( ( sine_angle1 != sine_angle2 ) || ( cosine_angle1 != cosine_angle2 ) )
    {
      PRINT_ERROR( "In AcisModifyEngine::get_conic_mid_surface\n"
		   "       Not supported yet: cones do not have the same major angle.\n" );
      return CUBIT_FAILURE;
    }
  }

  ellipse const* base1 = &cone_surface1->base;
  double ratio1 = base1->radius_ratio;
  ellipse const* base2 = &cone_surface2->base;
  double ratio2 = base2->radius_ratio;
  if( ratio1 != ratio2 )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_conic_mid_surface\n"
		 "       Not supported yet: cones do not have the same radius ratios.\n" );
    return CUBIT_FAILURE;
  }

  SPAunit_vector normal1 = base1->normal;
  SPAunit_vector normal2 = base2->normal;
  if ( ( normal1 != normal2 ) && ( - normal1 != normal2 ) )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_conic_mid_surface\n"
		 "       Not supported yet: cones do not have the same axis of symmetry.\n" );
    return CUBIT_FAILURE;
  }

  SPAposition const centre1 = base1->centre;
  SPAposition const centre2 = base2->centre;

  CubitVector null_vec( 0., 0., 0. );
  CubitVector c1_to_c2( centre2.x() - centre1.x(), 
			centre2.y() - centre1.y(), 
			centre2.z() - centre1.z() );
  if ( ! c1_to_c2.about_equal( null_vec ) )
  {
    CubitVector comm_dir( normal1.x(), normal1.y(), normal1.z() );
    c1_to_c2.normalize();
    if ( ( ! c1_to_c2.about_equal( comm_dir ) ) && ( ! (-c1_to_c2).about_equal( comm_dir ) ) )
    {
      CubitVector const point_1( ( centre1.x() + centre2.x() ) * .5, 
				 ( centre1.y() + centre2.y() ) * .5, 
				 ( centre1.z() + centre2.z() ) * .5 );
      CubitVector const point_2 = point_1 + comm_dir;
      CubitVector const point_3 = point_1 + ( comm_dir * c1_to_c2 );
      CubitStatus ret = this->get_mid_plane( point_1, point_2, point_3, 
					     body_to_trim_to, midsurface_body );

      return ret;
    }
  }

  CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
  
  BODY *sheet_BODY;
  CubitBoolean cylinder = CUBIT_FALSE;
  if ( ( sine_angle1 == 0.0 ) && ( fabs( cosine_angle1 ) == 1.0 ) )
  {
    cylinder = CUBIT_TRUE;
  }
  double const st_ang = 0.;
  double const end_ang = 2. * CUBIT_PI;
  double const height = 2. * (bounding_box.diagonal()).length();
  CubitVector center_vec;
  CubitVector temp_norm( normal1.x(), normal1.y(), normal1.z());
  if ( ! cylinder )
  {
    if ( sine_angle1 > 0.0 && cosine_angle1 > 0.0 )
    {
      sine_angle1 = - sine_angle1;
      temp_norm *= - 1.0;
    }
    else if ( sine_angle1 < 0.0 && cosine_angle1 < 0.0 )
    {
      cosine_angle1 = - cosine_angle1;
      temp_norm *= - 1.0;
    }
    //make center of surface along axis, 1/2 height of 
    //bounding box of entire model. 
    CubitVector root_vec( centre1.x(), centre1.y(), centre1.z() );
    center_vec = root_vec + ( -.5 * height) * temp_norm; 
   }
  else
  {
    CubitVector root_vec( centre1.x(), centre1.y(), centre1.z() );
    center_vec = .5 * height * temp_norm;
    center_vec = root_vec - center_vec;
  }
  SPAposition const center ( center_vec.x(), center_vec.y(), center_vec.z() );

  SPAvector const major_axis1 = base1->major_axis; 
  SPAvector const major_axis2 = base2->major_axis; 
  SPAvector newmajor_axis;
  if ( cosine_angle1 * cosine_angle2 < 0. )
   {
     newmajor_axis = ( major_axis1 - major_axis2 ) * .5;
   }
  else
   {
     newmajor_axis = ( major_axis1 + major_axis2 ) * .5;
   }

  FACE *new_FACE_ptr = NULL;
  outcome result = api_make_cnface( center, 
				    normal1, 
				    newmajor_axis,
				    ratio1, 
				    sine_angle1, 
				    cosine_angle1, 
				    st_ang,
				    end_ang, 
				    height, 
				    new_FACE_ptr );
  if (!result.ok() || new_FACE_ptr == NULL ||
      new_FACE_ptr->geometry() == NULL )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_conic_mid_surface\n"
		 "       ACIS api_make_cnface function failed.\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    return CUBIT_FAILURE;
  }
  FACE *face_list[1];
  face_list[0] = new_FACE_ptr;
  result = api_sheet_from_ff( 1, face_list, sheet_BODY );
  if ( !result.ok() || sheet_BODY == NULL || sheet_BODY->lump() == NULL
      || sheet_BODY->lump()->shell() == NULL ||
      sheet_BODY->lump()->shell()->first_face() == NULL )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_conic_mid_surface\n"
		 "       ACIS api_sheet_from_ff function failed.\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    return CUBIT_FAILURE;
  }
  result = api_body_to_2d( sheet_BODY );
  if ( !result.ok() )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_conic_mid_surface\n"
		 "       ACIS api_body_to_2d function failed.\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    return CUBIT_FAILURE;
  } 
  
  BODY* BODY_to_trim_to = AcisQueryEngine::get_BODY(body_to_trim_to);
  bool delete_attribs = CUBIT_TRUE;
  BODY* BODY_copy = this->copy_BODY(BODY_to_trim_to, delete_attribs);
  if ( BODY_copy == NULL )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_conic_mid_surface\n"
		 "       Cannot make a copy of the parent ACIS BODY of volume.\n" );
    return CUBIT_FAILURE;
  }

  result = api_intersect( BODY_copy, sheet_BODY );
  if ( !result.ok() )
  {
    api_del_entity( BODY_copy );    
    api_del_entity( sheet_BODY );    
    PRINT_ERROR( "In AcisModifyEngine::get_conic_mid_surface\n"
		 "       ACIS api_intersect function failed.\n" );
    return CUBIT_FAILURE;
  }

  midsurface_body = AcisQueryEngine::instance()->populate_topology_bridges(sheet_BODY);

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : get_toric_mid_surface
// Member Type: PUBLIC
// Description: Calculates a mid-surface between 2 toric surfaces.
// Author     : Philippe Pebay
// Date       : 03/07/06
//=============================================================================
CubitStatus AcisModifyEngine::get_toric_mid_surface( Surface *surface_ptr1,
				  Surface *surface_ptr2,
				  BodySM *body_to_trim_to,
				  BodySM *&midsurface_body ) const
{
  SurfaceACIS *surf_ACIS1 = CAST_TO(surface_ptr1, SurfaceACIS );
  FACE* FACE1_ptr = make_FACE( surf_ACIS1->get_FACE_ptr() ); 
  int surface_type = (&(FACE1_ptr->geometry()->equation()))->type();
  if( surface_type != torus_type )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_toric_mid_surface\n"
		 "       The CGM cone is not an ASIS torus (%d %d).\n", surface_type, CONE_SURFACE_TYPE );
    return CUBIT_FAILURE;
  }
  SurfaceACIS *surf_ACIS2 = CAST_TO(surface_ptr2, SurfaceACIS );
  FACE* FACE2_ptr = make_FACE( surf_ACIS2->get_FACE_ptr() ); 
  surface_type = (&(FACE2_ptr->geometry()->equation()))->type();
  if( surface_type != torus_type )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_toric_mid_surface\n"
		 "       The CGM cone is not an ASIS torus (%d %d).\n", surface_type, CONE_SURFACE_TYPE );
    return CUBIT_FAILURE;
  }

//   surface const* acis_surface1 = &(FACE1_ptr->geometry()->equation());
//   torus const* torus_surface1 = (torus *)acis_surface1;
  SURFACE const* acis_SURFACE1 = FACE1_ptr->geometry();
  TORUS const* TORUS_SURFACE1 = (TORUS *)acis_SURFACE1;
  double major_radius1 = TORUS_SURFACE1->major_radius();

//   surface const* acis_surface2 = &(FACE2_ptr->geometry()->equation());
//   torus const* torus_surface2 = (torus *)acis_surface2;
  SURFACE const* acis_SURFACE2 = FACE2_ptr->geometry();
  TORUS const* TORUS_SURFACE2 = (TORUS *)acis_SURFACE2;
  double major_radius2 = TORUS_SURFACE2->major_radius();
  if ( major_radius1 != major_radius2 )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_toric_mid_surface\n"
 		 "       Not supported yet: torii do not have the same radii.\n" );
    return CUBIT_FAILURE;
  }

  SPAposition center1 = TORUS_SURFACE1->centre();
  SPAposition center2 = TORUS_SURFACE2->centre();
  if ( center1 != center2 )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_toric_mid_surface\n"
 		 "       Not supported yet: torii do not have the same center.\n" );
    return CUBIT_FAILURE;
  }
  
  SPAunit_vector normal1 = TORUS_SURFACE1->normal();
  SPAunit_vector normal2 = TORUS_SURFACE2->normal();
  if ( ( normal1 != normal2 ) && ( - normal1 != normal2 ) )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_conic_mid_surface\n"
		 "       Not supported yet: torii do not have the same axis of symmetry.\n" );
    return CUBIT_FAILURE;
  }

  CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();

  BODY *sheet_BODY;
  double minor_radius = ( TORUS_SURFACE1->minor_radius() + TORUS_SURFACE2->minor_radius() ) / 2.;
  SPAposition pnt( major_radius1 + minor_radius, 0., 0. );
  double const uf = 0.;
  double const ut = 2. * CUBIT_PI;
  double const vf = 0.;
  double const vt = 2. * CUBIT_PI;
  FACE *new_FACE_ptr = NULL;
  outcome result = api_make_trface( center1, 
				    normal1, 
				    major_radius1, 
				    minor_radius, 
				    pnt,
				    uf, 
				    ut, 
				    vf, 
				    vt, 
				    new_FACE_ptr );
  if (!result.ok() || new_FACE_ptr == NULL ||
      new_FACE_ptr->geometry() == NULL )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_toric_mid_surface\n"
		 "       ACIS api_make_trface function failed.\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    return CUBIT_FAILURE;
  }
  FACE *face_list[1];
  face_list[0] = new_FACE_ptr;
  result = api_sheet_from_ff( 1, face_list, sheet_BODY );
  if ( !result.ok() || sheet_BODY == NULL || sheet_BODY->lump() == NULL
      || sheet_BODY->lump()->shell() == NULL ||
      sheet_BODY->lump()->shell()->first_face() == NULL )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_toric_mid_surface\n"
		 "       ACIS api_sheet_from_ff function failed.\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    return CUBIT_FAILURE;
  }
  result = api_body_to_2d( sheet_BODY );
  if ( !result.ok() )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_toric_mid_surface\n"
		 "       ACIS api_body_to_2d function failed.\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    return CUBIT_FAILURE;
  } 

  BODY* BODY_to_trim_to = AcisQueryEngine::get_BODY(body_to_trim_to);
  bool delete_attribs = CUBIT_TRUE;
  BODY* BODY_copy = this->copy_BODY(BODY_to_trim_to, delete_attribs);
  if ( BODY_copy == NULL )
  {
    PRINT_ERROR( "In AcisModifyEngine::get_toric_mid_surface\n"
		 "       Cannot make a copy of the parent ACIS BODY of volume.\n" );
    return CUBIT_FAILURE;
  }

  result = api_intersect( BODY_copy, sheet_BODY );
  if ( !result.ok() )
  {
    api_del_entity( BODY_copy );    
    api_del_entity( sheet_BODY );    
    PRINT_ERROR( "In AcisModifyEngine::get_toric_mid_surface\n"
		 "       ACIS api_intersect function failed.\n" );
    return CUBIT_FAILURE;
  }

  midsurface_body = AcisQueryEngine::instance()->populate_topology_bridges(sheet_BODY);

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer curves on solid bodies.  The left and right offsets are
//              with respect to the curve direction.  If the given right offset
//              is negative, the left offset is used.  Users can preview to
//              clarify the meaning of left and right.
// Author     : Steve Storm
// Date       : 03/26/05
//=============================================================================
CubitStatus AcisModifyEngine::tweak_chamfer( DLIList<Curve*> &curve_list,
                                             double left_offset,
                                             DLIList<BodySM*> &new_bodysm_list,
                                             double right_offset,
                                             CubitBoolean keep_old_body,
                                             CubitBoolean preview ) const
{
    return AcisTweakTool::instance()->tweak_chamfer( curve_list, left_offset,
      new_bodysm_list, right_offset, keep_old_body, preview );
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer vertices on solid or sheet bodies.  On a solid body
//              there can be up to 3 offsets; on a sheet body up to 2 offsets.
//              The offsets are in the direction of the supplied edges.  If
//              multiple vertices are supplied, only one offset value is
//              allowed and the edges are not used.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus
AcisModifyEngine::tweak_chamfer( DLIList<Point*> &ref_vertex_list,
                                 double offset1,
                                 DLIList<BodySM*> &new_bodysm_list,
                                 Curve *edge1,
                                 double offset2,
                                 Curve *edge2,
                                 double offset3,
                                 Curve *edge3,
                                 CubitBoolean keep_old_body,
                                 CubitBoolean preview ) const
{
  return AcisTweakTool::instance()->tweak_chamfer( ref_vertex_list, offset1,
    new_bodysm_list, edge1, offset2, edge2, offset3, edge3, keep_old_body,
    preview );
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on solid
//              bodies.
// Author     : Steve Storm
// Date       :
//=============================================================================
CubitStatus AcisModifyEngine::tweak_fillet( DLIList<Curve*> &curve_list,
                                            double radius,
                                            DLIList<BodySM*> &new_bodysm_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview ) const
{
    return AcisTweakTool::instance()->tweak_fillet( curve_list, radius,
      new_bodysm_list, keep_old_body, preview );
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on a solid
//              body.  The fillet has a variable radius from the start to the
//              end of the curve.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus AcisModifyEngine::tweak_fillet( Curve *curve_ptr,
                                            double start_radius,
                                            double end_radius,
                                            BodySM *&new_bodysm_ptr,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview ) const
{
  return AcisTweakTool::instance()->tweak_fillet( curve_ptr, start_radius,
    end_radius, new_bodysm_ptr, keep_old_body, preview );
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given vertices on sheet
//              bodies.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus
AcisModifyEngine::tweak_fillet( DLIList<Point*> &ref_vertex_list,
                                double radius,
                                DLIList<BodySM*> &new_bodysm_list,
                                CubitBoolean keep_old_body,
                                CubitBoolean preview ) const
{
  return AcisTweakTool::instance()->tweak_fillet( ref_vertex_list, radius,
    new_bodysm_list, keep_old_body, preview );
}

//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes along a vector.
// Author     : Steve Storm
// Date       :
//=============================================================================
CubitStatus AcisModifyEngine::tweak_move( DLIList<Surface*> &surface_list,
                                          const CubitVector &delta,
                                          DLIList<BodySM*> &new_bodysm_list,
                                          CubitBoolean keep_old_body,
                                          CubitBoolean preview ) const
{
   return AcisTweakTool::instance()->tweak_move( surface_list, delta,
     new_bodysm_list, keep_old_body, preview );
}

//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body along a vector.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus AcisModifyEngine::tweak_move( DLIList<Curve*> &curve_list,
                                          const CubitVector &delta,
                                          DLIList<BodySM*> &new_bodysm_list,
                                          CubitBoolean keep_old_body,
                                          CubitBoolean preview ) const
{
  return AcisTweakTool::instance()->tweak_move( curve_list, delta,
     new_bodysm_list, keep_old_body, preview );
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes by offsetting
//              those faces by the offset distance.
// Author     : Steve Storm
// Date       :
//=============================================================================
CubitStatus AcisModifyEngine::tweak_offset( DLIList<Surface*> &surface_list,
                                            double offset_distance,
                                            DLIList<BodySM*> &new_bodysm_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview ) const
{
   return AcisTweakTool::instance()->tweak_offset( surface_list,
     offset_distance, new_bodysm_list, keep_old_body, preview );
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body or bodies by offsetting
//              those curves by the offset distance.
// Author     :
// Date       :
//=============================================================================
CubitStatus AcisModifyEngine::tweak_offset( DLIList<Curve*> &curve_list,
                                            double offset_distance,
                                            DLIList<BodySM*> &new_bodysm_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview ) const
{
  return AcisTweakTool::instance()->tweak_offset( curve_list, offset_distance,
    new_bodysm_list, keep_old_body, preview );
}

//=============================================================================
// Function   : tweak_remove
// Member Type: PUBLIC
// Description: Function to remove surfaces from a body and then extend the
//              remaining surfaces to fill the gap or hole.
// Author     : Steve Storm
// Date       :
//=============================================================================
CubitStatus AcisModifyEngine::tweak_remove( DLIList<Surface*> &ref_face_list,
                                            DLIList<BodySM*> &new_bodysm_list,
                                            CubitBoolean extend_adjoining,
                                            CubitBoolean keep_surface,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview ) const
{
   return AcisTweakTool::instance()->tweak_remove( ref_face_list,
     new_bodysm_list, extend_adjoining, keep_surface, keep_old_body, preview );
}

//================================================================================
// Description: Remove curves from a sheet body or bodies and then extend the
//              remaining curves, remove the hole or remove the hardline.
// Author     : Steve Storm
// Date       : 03/16/05
//================================================================================
CubitStatus AcisModifyEngine::tweak_remove( DLIList<Curve*> &curve_list,
                                            DLIList<BodySM*> &new_bodysm_list,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview ) const
{
  return AcisTweakTool::instance()->tweak_remove( curve_list, new_bodysm_list,
    keep_old_body, preview );
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes up to a set target
//              surfaces.
// Author     : Steve Storm
// Date       :
//=============================================================================
CubitStatus AcisModifyEngine::tweak_target( DLIList<Surface*> &surface_list,
                                            DLIList<Surface*> &target_face_list,
                                            DLIList<BodySM*> &new_bodysm_list,
                                            CubitBoolean reverse_flg,
                                            CubitBoolean keep_old_body,
                                            CubitBoolean preview ) const
{
   return AcisTweakTool::instance()->tweak_target( surface_list, target_face_list,
     new_bodysm_list, reverse_flg, keep_old_body, preview );
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a set of target surfaces.
// Author     : Steve Storm
// Date       : 03/01/05
//=============================================================================
CubitStatus
AcisModifyEngine::tweak_target( DLIList<Curve*> &curve_list,
                                DLIList<Surface*> &target_surf_list,
                                DLIList<BodySM*> &new_body_list,
                                CubitBoolean reverse_flg,
                                CubitBoolean keep_old_bodies,
                                CubitBoolean preview ) const
{
  return AcisTweakTool::instance()->tweak_target( curve_list, target_surf_list,
    new_body_list, reverse_flg, keep_old_bodies, preview );
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a sheet body or bodies up to a set of
//              target curves that are part of a sheet body.  The target is a
//              surface created by thickening the owning surface of the target
//              curves.
// Author     : Steve Storm
// Date       : 03/01/05
//=============================================================================
CubitStatus
AcisModifyEngine::tweak_target( DLIList<Curve*> &curve_list,
                                DLIList<Curve*> &target_curve_list,
                                DLIList<BodySM*> &new_body_list,
                                CubitBoolean reverse_flg,
                                CubitBoolean keep_old_bodies,
                                CubitBoolean preview ) const
{
  return AcisTweakTool::instance()->tweak_target( curve_list, target_curve_list,
    new_body_list, reverse_flg, keep_old_bodies, preview );
}

//================================================================================
// Description: Creates a net surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus AcisModifyEngine::create_net_surface( DLIList<Surface*>& /*ref_face_list*/, BodySM *& new_body,
                                                  DLIList<DLIList<CubitVector*>*> &vec_lists_u,
                                                  DLIList<DLIList<CubitVector*>*> &vec_lists_v,
                                                  double net_tol, CubitBoolean heal ) const
{
   return AcisSurfaceTool::instance()->create_net_surface( new_body, vec_lists_u, vec_lists_v,
                                                           net_tol, heal );
}

//================================================================================
// Description: Creates a net surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus AcisModifyEngine::create_net_surface( DLIList<Curve*>& u_curves, DLIList<Curve*>& v_curves,
                                                  BodySM *& new_body, double net_tol, CubitBoolean heal ) const
{
   return AcisSurfaceTool::instance()->create_net_surface( u_curves, v_curves, new_body, net_tol, heal );
}

//================================================================================
// Description: Creates an offset surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus AcisModifyEngine::create_offset_surface( Surface* ref_face_ptr, BodySM*& new_body, double offset_distance ) const
{
   return AcisSurfaceTool::instance()->create_offset_surface( ref_face_ptr, new_body, offset_distance );
}

//================================================================================
// Description: Creates an offset body.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus AcisModifyEngine::create_offset_body( BodySM* body_ptr, BodySM*& new_body, double offset_distance ) const
{
   return AcisSurfaceTool::instance()->create_offset_body( body_ptr, new_body, offset_distance );
}

//================================================================================
// Description: Creates a skin surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus AcisModifyEngine::create_skin_surface( DLIList<Curve*>& curves, BodySM*& new_body ) const
{
   return AcisSurfaceTool::instance()->create_skin_surface( curves, new_body );
}

//================================================================================
// Description: Creates a body from lofting surfaces.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus AcisModifyEngine::loft_surfaces( Surface *face1, const double &takeoff1,
                                             Surface *face2, const double &takeoff2,
                                             BodySM*& new_body,
                                             CubitBoolean arc_length_option, CubitBoolean twist_option,
                                             CubitBoolean align_direction, CubitBoolean perpendicular,
                                             CubitBoolean simplify_option ) const
{
   return AcisSurfaceTool::instance()->loft_surfaces( face1, takeoff1, face2, takeoff2, new_body,
                                                      arc_length_option, twist_option, align_direction,
                                                      perpendicular, simplify_option );
}

//================================================================================
// Description: Creates a body by lofting surfaces between bodies.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus AcisModifyEngine::loft_surfaces_to_body( Surface *face1, const double &takeoff1,
                                                     Surface *face2, const double &takeoff2,
                                                     BodySM*& new_body,
                                                     CubitBoolean arc_length_option, CubitBoolean twist_option,
                                                     CubitBoolean align_direction, CubitBoolean perpendicular,
                                                     CubitBoolean simplify_option ) const
{
   return AcisSurfaceTool::instance()->loft_surfaces_to_body( face1, takeoff1, face2, takeoff2, new_body,
                                                              arc_length_option, twist_option, align_direction,
                                                              perpendicular, simplify_option );
}

//================================================================================
// Description: Creates a surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus AcisModifyEngine::create_surface( DLIList<CubitVector*>& vec_list, BodySM *&new_body, Surface *ref_face_ptr,
                                              CubitBoolean project_points ) const
{
   return AcisSurfaceTool::instance()->create_surface( vec_list, new_body, ref_face_ptr, project_points );
}

//================================================================================
// Description: Creates a weld surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus AcisModifyEngine::create_weld_surface( CubitVector &root,
                                                   Surface *ref_face1, double leg1, Surface *ref_face2, double leg2,
                                                   BodySM *&new_body ) const
{
   return AcisSurfaceTool::instance()->create_weld_surface( root, ref_face1, leg1, ref_face2, leg2, new_body );
}

//-------------------------------------------------------------------------
// Purpose       : Perform the acis slice routine.  Returns the acis
//                 intersection graph.
//
// Special Notes : The api_slice routine works a bit strangely.
//                 If no intersection is detected (i.e., the 2
//                 BODYs do not even touch at a single point, then
//                 it not only returns a NULL BODY, but also
//                 returns an outcome whose status is not ok.
//
//                 Note also that this is different from the next function
//                 in that the intersection graph is wanted so the
//                 roll back mechanism implored there to save having
//                 to make copies for safety can't be done here.
//
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::FACE_intersection( FACE *FACE1_ptr,
                                                 FACE *FACE2_ptr,
                                                 BODY *&intersect_graph,
                                                 CubitBoolean &interfering ) const
{
    // Make sure we have 2 non-NULL, valid BODYs first.
  if ( FACE1_ptr == NULL || !IS_ENTITY_TYPE( FACE1_ptr, FACE ) ||
       FACE2_ptr == NULL || !IS_ENTITY_TYPE( FACE2_ptr, FACE ) )
  {
    PRINT_ERROR("BUG: In AcisModifyEngine::BODY_intersection\n"
                "     One or both input FACE's are NULL or invalid.\n");
    assert(0);
    return CUBIT_FAILURE;
  }
    //make copies to make sure nothing bad happens to the FACES.
  FACE* copy_1 = make_FACE(FACE1_ptr);
  FACE* copy_2 = make_FACE(FACE2_ptr);

  intersect_graph = NULL;

    // Create the intersection graph between the given bodies
  outcome result = api_fafa_int(copy_1, copy_2,
                                intersect_graph);

  if (!result.ok() || intersect_graph == NULL )
  {
      // Find out if the returned "slice" is empty
    interfering = CUBIT_FALSE;
  }
  else
  {
    interfering = CUBIT_TRUE;
  }
    //now delete the copies.
  api_delent(copy_1);
  api_delent(copy_2);
    //return success.
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Converts the wire body into a list of Curves.
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::convert_WIRE( BODY *wire_BODY,
                                            DLIList <Curve*> &curve_list ) const
{
    //Now get the edges from the wire body (int_graph), and create
    //ref edges.
    DLIList <EDGE*> EDGE_list;
    WIRE *cur_WIRE = NULL;
    LUMP *cur_LUMP = NULL;
    // Make sure the wire is OK, and get the appropriate entity
    // Get a WIRE if it's there
    if ( wire_BODY->wire() == NULL ||
        wire_BODY->wire()->coedge() == NULL ||
        wire_BODY->wire()->coedge()->edge() == NULL )
    {
        // If there was no WIRE, see if there's a lump
        if(wire_BODY->lump() != NULL &&
            wire_BODY->lump()->shell() != NULL &&
            wire_BODY->lump()->shell()->wire() != NULL &&
            wire_BODY->lump()->shell()->wire()->coedge() != NULL &&
            wire_BODY->lump()->shell()->wire()->coedge()->edge() != NULL )
        {
            cur_LUMP = wire_BODY->lump();
        }
        else
        {
            PRINT_ERROR("Non-existent wire body in convert_WIRE function.\n");
            return CUBIT_FAILURE;
        }
    }
    else
    {
        cur_WIRE = wire_BODY->wire();
    }
    LUMP *start_LUMP = cur_LUMP;
    do{
        if((wire_BODY->wire() == NULL ||
            wire_BODY->wire()->coedge() == NULL ||
            wire_BODY->wire()->coedge()->edge() == NULL ) &&
            wire_BODY->lump() != NULL &&
            wire_BODY->lump()->shell() != NULL &&
            wire_BODY->lump()->shell()->wire() != NULL &&
            wire_BODY->lump()->shell()->wire()->coedge() != NULL &&
            wire_BODY->lump()->shell()->wire()->coedge()->edge() != NULL )
        {
            cur_WIRE = cur_LUMP->shell()->wire();
        }

        WIRE *start_WIRE = cur_WIRE;

        if (cur_WIRE != NULL){
            do {
                COEDGE *start_COEDGE = cur_WIRE->coedge();
                COEDGE *cur_COEDGE = start_COEDGE;
                do{
                    EDGE *cur_EDGE = cur_COEDGE->edge();
                    //Just get one of the edges.  Often you have two coedges
                    //for the same edge.
                    if ( EDGE_list.move_to(cur_EDGE) )
                    {
                        // additional check to avoid an infinite loop
                        if(cur_COEDGE == cur_COEDGE->next())
                            break;

                        cur_COEDGE = cur_COEDGE->next();
                        continue;
                    }
                    else
                    {
                        EDGE_list.append(cur_EDGE);
                    }
                    EDGE *new_EDGE_ptr = NULL;
                    outcome result = api_edge(cur_EDGE, new_EDGE_ptr);
                    if ( !result.ok() || new_EDGE_ptr == NULL )
                    {
                        PRINT_ERROR("In AcisModifyEngine::body_intersection\n"
                            "       Cannot copy curve in intersection graph.\n");
                        return CUBIT_FAILURE;
                    }
                    // Create a new CurveACIS object
                    // Use the new EDGE to create a new CurveACIS
                    Curve *new_curve =
                        AcisQueryEngine::instance()->populate_topology_bridges(new_EDGE_ptr);
                    curve_list.append(new_curve);
                    cur_COEDGE = cur_COEDGE->next();
                } while (cur_COEDGE && (cur_COEDGE != start_COEDGE));
                cur_WIRE = cur_WIRE->next();
            } while (cur_WIRE && (cur_WIRE != start_WIRE));
        }

        //we may not have a lump, ACIS sometimes will return a wire
        //as the intersection and other times a LUMP...
        if (cur_LUMP)
        {
            cur_LUMP = cur_LUMP->next();
        }
    }while (cur_LUMP && (start_LUMP != cur_LUMP));
    //At this point we are done!
    return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Determines whether the 2 input ACIS BODYs overlap
//                 in any way at all.
// Special Notes : The api_slice routine works a bit strangely.
//                 If no intersection is detected (i.e., the 2
//                 BODYs do not even touch at a single point, then
//                 it not only returns a NULL BODY, but also
//                 returns an outcome whose status is not ok.
//
//-------------------------------------------------------------------------
CubitBoolean AcisModifyEngine::BODYs_interfering (
    BODY* body1_ptr,
    BODY* body2_ptr) const
{

   CpuTimer timer;

     // Efficient algorithm to test interference between 2 ACIS BODYs. The
     // entire intersection boolean operation is *not* performed if api_slice
     // is used.

   CubitBoolean interfering = CUBIT_TRUE;

     // Use the NOP pair to ensure that the roll-back is done by ACIS after
     // the test, so that the input BODYs are unaffected by this operation.
   API_NOP_BEGIN;

   BODY* slice = NULL;

     // Make sure we have 2 non-NULL, valid BODYs first.
   if ( body1_ptr == NULL || !IS_ENTITY_TYPE( body1_ptr, BODY ) ||
        body2_ptr == NULL || !IS_ENTITY_TYPE( body2_ptr, BODY ) )
   {
      PRINT_ERROR("BUG: In AcisModifyEngine::BODYs_interfering\n"
                  "     One or both input BODYs are NULL or invalid.\n");
      assert(0);
      return CUBIT_FALSE;
   }

     // Create the intersection graph between the given bodies
   outcome status = api_slice(body1_ptr, body2_ptr,
                              *(SPAunit_vector *)NULL_REF, slice);

   if (!status.ok() || slice == NULL )
       interfering = CUBIT_FALSE;
   else
     interfering = CUBIT_TRUE;

     // Set the result variable to true so that the API_NOP_END
     // can roll back the bulletin board to ensure that the
     // ENTITYs that were modified within this NOP block will be
     // returned to their original state.
   result = outcome(0);

   API_NOP_END;


   PRINT_DEBUG_3( "CPU time taken to check if two BODYs are "
               "interfering: %f secs\n", timer.cpu_secs());

     // Return the interference status
   return interfering;
}

CubitStatus AcisModifyEngine::stitch_FACEs( DLIList<FACE*> &faces,
                                                    FACE *&result_FACE,
                                                    BODY *&stitched_BODY )
{
  //Make copies of each FACE
  ENTITY_LIST copied_FACES;
  BODY *sheet_BODY_ptr;
  int i;
  outcome result;

  BODY *master = NULL;
  for( i=faces.size(); i--;)
  {
    FACE *tmp_FACE = faces.get_and_step();
    ENTITY *new_FACE = NULL;
    result = api_copy_entity( (ENTITY*)tmp_FACE, new_FACE );

    FACE *face_list[1];
    face_list[0] = (FACE*)new_FACE;
    result = api_sheet_from_ff( 1, face_list, sheet_BODY_ptr );
    result = api_body_to_2d( sheet_BODY_ptr );

    if( !result.ok() )
    {
      PRINT_ERROR("Problems copying surface(s)\n");
      //delete all other surface
      for(i=0; i<copied_FACES.count(); i++)
        api_delent( copied_FACES[i]);
      return CUBIT_FAILURE;
    }
    else
    {
      if( master == NULL )
        master = sheet_BODY_ptr;
      else
        copied_FACES.add( sheet_BODY_ptr );
    }
  }

  //stitch together to form a sheet body.
#if CUBIT_ACIS_VERSION < 1100

  for( i=0; i<copied_FACES.count(); i++ )
  {
    result = api_unite((BODY*)copied_FACES[i], master );
    if( !result.ok() )
    {
      PRINT_ERROR("Problems stitching surfaces\n");
      //delete all other surface
      for(; i<copied_FACES.count(); i++)
        api_delent( copied_FACES[i]);
      return CUBIT_FAILURE;
    }
  }
#else
  exact_stitch_options stitch_opts;
  result = api_stitch( master, copied_FACES, &stitch_opts );
#endif

  if( !result.ok() )
  {
    api_delent( master );
    stitched_BODY = NULL;
    for(i=0; i<copied_FACES.count(); i++)
      api_delent( copied_FACES[i]);
    result_FACE = NULL;
    return CUBIT_FAILURE;
  }

  copied_FACES.clear();
  stitched_BODY = master;

  //should only have 1 FACE
  api_clean_entity( (ENTITY*)master );
  api_get_faces( master, copied_FACES );

  if( copied_FACES.count() != 1 )
  {
    result_FACE = NULL;
    return CUBIT_SUCCESS;
  }
  else
    result_FACE = (FACE*)copied_FACES[0];

  return CUBIT_SUCCESS;
}

CubitStatus AcisModifyEngine::webcut_with_sweep_curves_rotated(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Curve*> &curves,
                              const CubitVector& point,
                              const CubitVector& sweep_axis,
                              double angle,
                              Surface *stop_surf,
                              DLIList<BodySM*> &results_list,
                              CubitBoolean imprint)
{
  if(curves.size() == 0 )
    return CUBIT_FAILURE;

  //see if surface-where-sweep-is-to-stop has been specified
  EDGE *EDGE_to_sweep = NULL;
  FACE *stop_FACE= NULL;
  if( stop_surf )
  {
    SurfaceACIS *surf_ACIS = CAST_TO(stop_surf, SurfaceACIS );
    stop_FACE = surf_ACIS->get_FACE_ptr();
  }

  //sweep each curve about the specified axis
  int i,j;
  DLIList<FACE*> FACES_for_cutting;
  for(i=curves.size(); i--; )
  {
    Curve *curve = curves.get_and_step();

    CurveACIS *curve_ACIS = CAST_TO(curve, CurveACIS);
    EDGE *EDGE_ptr = curve_ACIS->get_EDGE_ptr();

    if (!EDGE_ptr)
    {
      PRINT_ERROR("Unable to get ACIS EDGE from Curve\n" );
      for( j=FACES_for_cutting.size(); j--; )
        api_delent( FACES_for_cutting.get_and_step() );
      return CUBIT_FAILURE;
    }

    outcome result = api_edge( EDGE_ptr, EDGE_to_sweep);
    if (!result.ok() || EDGE_to_sweep == NULL )
    {
      AcisQueryEngine::instance()->ACIS_API_error (result);
      for( j=FACES_for_cutting.size(); j--; )
        api_delent( FACES_for_cutting.get_and_step() );
      PRINT_ERROR( "Unable to copy ACIS curve\n" );
      return CUBIT_FAILURE;
    }

    BODY *BODY_of_swept_EDGE = NULL;
    CubitStatus status;

    status = sweep_EDGE_about_axis( EDGE_to_sweep, BODY_of_swept_EDGE,
                                    sweep_axis, point, angle, 0, 0.0, 0,
                                    false, false, stop_FACE );

    if( status == CUBIT_FAILURE )
    {
      api_delent( EDGE_to_sweep );
      int j;
      for( j=FACES_for_cutting.size(); j--; )
        api_delent( FACES_for_cutting.get_and_step() );
      return CUBIT_FAILURE;
    }

    FACE *FACE_ptr = BODY_of_swept_EDGE->lump()->shell()->first_face();
    FACES_for_cutting.append( FACE_ptr );
  }

  //stitch FACEs together into one BODY
  //they may not stitch into on single FACE... a sheet-body
  //with multiple FACEs will still webcut though
  BODY *tool_BODY = NULL;
  if( FACES_for_cutting.size() > 1 )
  {
    FACE *dummy_FACE = NULL;
    CubitStatus status;
    status = stitch_FACEs( FACES_for_cutting, dummy_FACE, tool_BODY );

    if( status == CUBIT_FAILURE )
    {
      PRINT_ERROR("Swept curves produce surfaces that cannot be stitched together.\n");
      for( i=FACES_for_cutting.size(); i--; )
        api_delent( FACES_for_cutting.get_and_step() );
      return CUBIT_FAILURE;
    }
  }
  else
    tool_BODY = AcisQueryEngine::instance()->get_BODY_of_ENTITY( FACES_for_cutting.get() );

  //get rid of any unnecessary edges
  api_clean_entity( (ENTITY*)tool_BODY);

  //webcut each blank BODY with the tool BODY
  int num_cut=0;
  DLIList<BodySM*> temp_new_bodies;

  for ( i = blank_bodies.size(); i > 0; i-- )
  {
     BodySM *blank_body = blank_bodies.get_and_step();
     BODY *blank_BODY = AcisQueryEngine::get_BODY(blank_body );
      //now webcut it.
     BODY *new_BODY_1, *new_BODY_2;

     CubitStatus status = webcut_with_sheet( blank_BODY,
                                             tool_BODY,
                                             new_BODY_1,
                                             new_BODY_2,
                                             imprint );
     if ( status == CUBIT_SUCCESS )
     {
       DLIList<BODY*> new_BODY_list;
       new_BODY_list.append(new_BODY_1);
       new_BODY_list.append(new_BODY_2);
       temp_new_bodies.clean_out();

       CubitStatus result = get_new_Body(blank_body, blank_BODY,
                                         new_BODY_list, temp_new_bodies,
                                         CUBIT_FALSE);
       if ( result != CUBIT_SUCCESS )
       {
         api_delent( tool_BODY );
         PRINT_ERROR("Problems with building volume.\n");
         num_cut = 0;
         return CUBIT_FAILURE;
       }
       else
       {
         results_list += temp_new_bodies;
         num_cut++;
       }
     }
   }

  return num_cut > 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

CubitStatus AcisModifyEngine::webcut_with_sweep_surfaces_rotated(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Surface*> &surfaces,
                              const CubitVector& point,
                              const CubitVector& sweep_vector,
                              double angle,
                              Surface *stop_surf,
                              bool up_to_next,
                              DLIList<BodySM*> &results_list,
                              CubitBoolean imprint )
{
  FACE *FACE_to_sweep = NULL;

  int i;
  CubitStatus status;
  DLIList<FACE*> FACES_to_sweep;
  if( surfaces.size() > 1 )
  {
    DLIList<FACE*> faces;
    for(i=surfaces.size(); i--;)
    {
      SurfaceACIS *surf_ACIS = CAST_TO(surfaces.get_and_step(), SurfaceACIS );
      faces.append( surf_ACIS->get_FACE_ptr() );
    }

    //attempt to stitch FACEs together, to get a single FACE
    BODY *stitched_BODY = NULL;
    status = stitch_FACEs( faces, FACE_to_sweep, stitched_BODY );

    //can't stitch surfaces together, so we'll sweep them separately
    if( status == CUBIT_FAILURE )
    {
      PRINT_ERROR("Input surfaces could not be stitched together to form a single surface.\n");
      PRINT_INFO("       Input surfaces must share coinsident curves to be stitchable.\n");
      return CUBIT_FAILURE;
    }
    else if( !FACE_to_sweep ) //couldn't stitch into one FACE....will sweep surfaces individually
    {
      if( stitched_BODY ) api_delent( stitched_BODY );
      for(i=surfaces.size(); i--;)
      {
        SurfaceACIS *surf_ACIS = CAST_TO(surfaces.get_and_step(), SurfaceACIS );
        FACE_to_sweep = make_FACE( surf_ACIS->get_FACE_ptr() );
        FACES_to_sweep.append( FACE_to_sweep );
      }
    }
    else
      FACES_to_sweep.append( FACE_to_sweep );
  }
  else
  {
    Surface *surface = surfaces.get();
    SurfaceACIS *surf_ACIS = CAST_TO(surface, SurfaceACIS );
    FACE_to_sweep = make_FACE( surf_ACIS->get_FACE_ptr() );
    FACES_to_sweep.append( FACE_to_sweep );
  }

  //see if surface-where-sweep-is-to-stop has been specified
  FACE *stop_FACE = NULL;
  if( stop_surf )
  {
    SurfaceACIS *surf_ACIS = CAST_TO(stop_surf, SurfaceACIS );
    stop_FACE = surf_ACIS->get_FACE_ptr();
  }

  //sweep each surface about the specified axis
  bool volume_is_negative;
  DLIList<BODY*> swept_BODYs;
  for( i=FACES_to_sweep.size(); i--; )
  {
    FACE_to_sweep = FACES_to_sweep.get_and_step();

    status = sweep_FACE_about_axis( FACE_to_sweep, sweep_vector, point, angle,
                                    volume_is_negative, true, 0, 0.0, 0, false, stop_FACE );

    if( status == CUBIT_FAILURE )
    {
      api_delent( FACE_to_sweep );
      for(i=swept_BODYs.size(); i--;)
        api_delent( swept_BODYs.get_and_step() );
      return CUBIT_FAILURE;
    }

    BODY *swept_BODY = AcisQueryEngine::instance()->get_BODY_of_ENTITY(FACE_to_sweep);
    swept_BODYs.append( swept_BODY );
  }

  //if there are more than 1, unite them all
  BODY *tool_BODY = swept_BODYs.get_and_step();
  if( swept_BODYs.size() > 1 )
  {
    for(i=swept_BODYs.size()-1; i--;)
    {
      outcome result = api_unite( swept_BODYs.get_and_step(), tool_BODY );
      if( !result.ok() )
      {
        api_delent( tool_BODY );
        for(; i--;  )
          api_delent( swept_BODYs.get_and_step() );
        return CUBIT_FAILURE;
      }
    }
  }

  //get rid of any unnecessary edges
  api_clean_entity( (ENTITY*)tool_BODY);

  //trim tool_body if 'up_to_next_surface' option is specified
  if( up_to_next)
  {
    //get a vertex on the original face that we're sweeping
    SurfaceACIS *surf_ACIS = CAST_TO(surfaces.get(), SurfaceACIS );

    ENTITY_LIST vertices;
    api_get_vertices( surf_ACIS->get_FACE_ptr(), vertices );

    VERTEX *surf_VERT = NULL;
    if( vertices.count() )
      surf_VERT = static_cast<VERTEX*>( vertices[0] );

    if( !trim_up_to_next_surface( tool_BODY, blank_bodies, surf_VERT ) )
    {
      api_delent( tool_BODY );
      PRINT_ERROR("Cannot use 'up_to_next' option with specified geometry\n");
      PRINT_INFO("Try the 'stop surface <id>' option instead\n");
      return CUBIT_FAILURE;
    }
  }

  return webcut( blank_bodies, tool_BODY, results_list, imprint );

}

CubitStatus AcisModifyEngine::webcut_with_sweep_curves(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Curve*> &curves,
                              const CubitVector& sweep_vector,
                              bool through_all,
                              Surface *stop_surf,
                              Curve *curve_to_sweep_along,
                              DLIList<BodySM*> &results_list,
                              CubitBoolean imprint )
{
  EDGE *EDGE_to_sweep = NULL;
  if(curves.size() == 0 )
    return CUBIT_FAILURE;

  CubitVector tmp_sweep_vector = sweep_vector;

  //get model bbox info...will scale sweep vector by its diagonal
  //so that we go far enough
  FACE *stop_FACE= NULL;
  if( through_all || stop_surf )
  {
    CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
    tmp_sweep_vector.normalize();
    tmp_sweep_vector*=(2*bounding_box.diagonal().length());
  }

  //see if surface-where-sweep-is-to-stop has been specified
  if( stop_surf )
  {
    SurfaceACIS *surf_ACIS = CAST_TO(stop_surf, SurfaceACIS );
    stop_FACE = surf_ACIS->get_FACE_ptr();
  }

  //see if we're sweeping along a specified curve
  BODY *WIRE_ptr = NULL;
  if( curve_to_sweep_along )
  {
    CurveACIS *curve_ACIS = CAST_TO(curve_to_sweep_along, CurveACIS);
    EDGE *EDGE_ptr = curve_ACIS->get_EDGE_ptr();

    if (!EDGE_ptr)
    {
      PRINT_ERROR("Unable to get ACIS EDGE from Curve\n" );
      return CUBIT_FAILURE;
    }

    EDGE *EDGE_to_sweep_along = NULL;
    outcome result = api_edge( EDGE_ptr, EDGE_to_sweep_along );
    if (!result.ok() || EDGE_to_sweep_along  == NULL )
    {
      AcisQueryEngine::instance()->ACIS_API_error (result);
      PRINT_ERROR( "Unable to copy ACIS curve\n" );
      return CUBIT_FAILURE;
    }

    EDGE* EDGES[1];
    EDGES[0] = EDGE_to_sweep_along;
    result = api_make_ewire( 1, EDGES, WIRE_ptr);
    if (!result.ok())
    {
      AcisQueryEngine::instance()->ACIS_API_error (result);
      PRINT_ERROR( "Unable to make ACIS wire body from curve to sweep along\n" );
    }
  }

  //sweep each curve
  int i,j;
  DLIList<FACE*> FACES_for_cutting;
  for(i=curves.size(); i--; )
  {
    Curve *curve = curves.get_and_step();

    CurveACIS *curve_ACIS = CAST_TO(curve, CurveACIS);
    EDGE *EDGE_ptr = curve_ACIS->get_EDGE_ptr();

    if (!EDGE_ptr)
    {
      PRINT_ERROR("Unable to get ACIS EDGE from Curve\n" );
      return CUBIT_FAILURE;
    }

    outcome result = api_edge( EDGE_ptr, EDGE_to_sweep);
    if (!result.ok() || EDGE_to_sweep == NULL )
    {
      AcisQueryEngine::instance()->ACIS_API_error (result);
      PRINT_ERROR( "Unable to copy ACIS curve\n" );
      for( j=FACES_for_cutting.size(); j--; )
        api_delent( FACES_for_cutting.get_and_step() );
      if( WIRE_ptr )
        api_delent( WIRE_ptr );
      return CUBIT_FAILURE;
    }

    BODY *BODY_of_swept_EDGE = NULL;
    CubitStatus status;

    if( WIRE_ptr )
      status = sweep_EDGE_along_WIRE( EDGE_to_sweep, WIRE_ptr, BODY_of_swept_EDGE,
                                      0.0, 0, false, stop_FACE );
    else
      status = sweep_EDGE_along_vector( EDGE_to_sweep, BODY_of_swept_EDGE,
                                        tmp_sweep_vector, 0.0, 0, false, stop_FACE );

    if( status == CUBIT_FAILURE )
    {
      api_delent( EDGE_to_sweep );
      for( j=FACES_for_cutting.size(); j--; )
        api_delent( FACES_for_cutting.get_and_step() );
      if( WIRE_ptr )
        api_delent( WIRE_ptr );
      return CUBIT_FAILURE;
    }

    FACE *FACE_ptr = BODY_of_swept_EDGE->lump()->shell()->first_face();
    FACES_for_cutting.append( FACE_ptr );
  }

  if( WIRE_ptr )
    api_delent( WIRE_ptr );

  //stitch FACEs together into one BODY
  //they may not stitch into on single FACE... a sheet-body
  //with multiple FACEs will still webcut though
  BODY *tool_BODY = NULL;
  if( FACES_for_cutting.size() > 1 )
  {
    FACE *dummy_FACE = NULL;
    CubitStatus status;
    status = stitch_FACEs( FACES_for_cutting, dummy_FACE, tool_BODY );

    if( status == CUBIT_FAILURE )
    {
      PRINT_ERROR("Swept curves produce surfaces that cannot be stitched together.\n");
      for( i=FACES_for_cutting.size(); i--; )
        api_delent( FACES_for_cutting.get_and_step() );
      return CUBIT_FAILURE;
    }
  }
  else
    tool_BODY = AcisQueryEngine::instance()->get_BODY_of_ENTITY( FACES_for_cutting.get() );

  //get rid of any unnecessary edges
  api_clean_entity( (ENTITY*)tool_BODY);

  //webcut each blank BODY with the tool BODY
  int num_cut=0;
  DLIList<BodySM*> temp_new_bodies;
  for ( i = blank_bodies.size(); i > 0; i-- )
  {
     BodySM *blank_body = blank_bodies.get_and_step();
     BODY *blank_BODY = AcisQueryEngine::get_BODY(blank_body );
     BODY *new_BODY_1, *new_BODY_2;

     CubitStatus status = webcut_with_sheet( blank_BODY,
                                             tool_BODY,
                                             new_BODY_1,
                                             new_BODY_2,
                                             imprint );
     if ( status == CUBIT_SUCCESS )
     {
       DLIList<BODY*> new_BODY_list;
       new_BODY_list.append(new_BODY_1);
       new_BODY_list.append(new_BODY_2);
       temp_new_bodies.clean_out();

       CubitStatus result = get_new_Body(blank_body, blank_BODY,
                                         new_BODY_list, temp_new_bodies,
                                         CUBIT_FALSE);
       if ( result != CUBIT_SUCCESS )
       {
         PRINT_ERROR("Problems with building volume.\n");
         num_cut = 0;
         return CUBIT_FAILURE;
       }
       else
       {
         results_list += temp_new_bodies;
         num_cut++;
       }
     }
   }

  return num_cut > 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}


CubitStatus AcisModifyEngine::webcut_with_sweep_surfaces(
                              DLIList<BodySM*> &blank_bodies,
                              DLIList<Surface*> &surfaces,
                              const CubitVector& sweep_vector,
                              bool sweep_perp,
                              bool through_all,
                              bool outward,
                              bool up_to_next,
                              Surface *stop_surf,
                              Curve *curve_to_sweep_along,
                              DLIList<BodySM*> &results_list,
                              CubitBoolean imprint )
{
  FACE *FACE_to_sweep = NULL;
  SurfaceACIS *surf_ACIS = NULL;

  CubitVector tmp_sweep_vector = sweep_vector;

  int i,j;
  DLIList<FACE*> FACES_to_sweep;

  //if multiple FACEs stitch, them together into one BODY,
  //they may not stitch into on single FACE... a sheet-body
  //with multiple FACEs will still webcut though
  if( surfaces.size() > 1 )
  {
    DLIList<FACE*> faces;
    for(i=surfaces.size(); i--;)
    {
      SurfaceACIS *surf_ACIS = CAST_TO(surfaces.get_and_step(), SurfaceACIS );
      faces.append( surf_ACIS->get_FACE_ptr() );
    }

    //attempt to stitch FACEs together, to get a single FACE
    BODY *dummy_BODY = NULL;
    CubitStatus status;
    status = stitch_FACEs( faces, FACE_to_sweep, dummy_BODY );

    if( status == CUBIT_FAILURE )
    {
      PRINT_ERROR("Input surfaces could not be stitched together.\n");
      PRINT_INFO("       Input surfaces must share coinsident curves to be stitchable.\n");
      return CUBIT_FAILURE;
    }
    else if( !FACE_to_sweep )
    {
      //can't stitch surfaces into single FACE, so we'll sweep them separately
      if( dummy_BODY ) api_delent( dummy_BODY );
      for(i=surfaces.size(); i--;)
      {
        surf_ACIS = CAST_TO(surfaces.get_and_step(), SurfaceACIS );
        FACE_to_sweep = make_FACE( surf_ACIS->get_FACE_ptr() );
        FACES_to_sweep.append( FACE_to_sweep );
      }
    }
    else
      FACES_to_sweep.append( FACE_to_sweep );
  }
  else //if only 1 surface
  {
    Surface *surface = surfaces.get();
    surf_ACIS = CAST_TO(surface, SurfaceACIS );
    FACE_to_sweep = make_FACE( surf_ACIS->get_FACE_ptr() );
    FACES_to_sweep.append( FACE_to_sweep );
  }

  FACE *stop_FACE= NULL;
  if( sweep_perp == true )
  {
    //if( !stop_surf )
    //  tmp_sweep_vector.set( sweep_vector.x(),0,0);
    if( through_all || stop_surf )
    {
      CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
      tmp_sweep_vector.set(1,0,0);
      tmp_sweep_vector = 2*(bounding_box.diagonal());
    }
  }
  else if( through_all || stop_surf || up_to_next  )
  {
    CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
    tmp_sweep_vector.normalize();
    tmp_sweep_vector*=(2*bounding_box.diagonal().length());
  }

  //see if surface-where-sweep-is-to-stop has been specified
  if( stop_surf )
  {
    //get the SURFACE
    SurfaceACIS *surf_ACIS = CAST_TO(stop_surf, SurfaceACIS );
    stop_FACE = surf_ACIS->get_FACE_ptr();
  }

  //see if we're sweeping along a specified curve
  BODY* WIRE_ptr = NULL;
  if( curve_to_sweep_along )
  {
    CurveACIS *curve_ACIS = CAST_TO(curve_to_sweep_along, CurveACIS);
    EDGE *EDGE_ptr = curve_ACIS->get_EDGE_ptr();

    if (!EDGE_ptr)
    {
      PRINT_ERROR("Unable to get ACIS EDGE from Curve\n" );
      for( j=FACES_to_sweep.size(); j--; )
        api_delent( FACES_to_sweep.get_and_step() );
      return CUBIT_FAILURE;
    }

    EDGE *EDGE_to_sweep_along = NULL;
    outcome result = api_edge( EDGE_ptr, EDGE_to_sweep_along );
    if (!result.ok() || EDGE_to_sweep_along  == NULL )
    {
      AcisQueryEngine::instance()->ACIS_API_error (result);
      PRINT_ERROR( "Unable to copy ACIS curve\n" );
      for( j=FACES_to_sweep.size(); j--; )
        api_delent( FACES_to_sweep.get_and_step() );
      return CUBIT_FAILURE;
    }

    EDGE* EDGES[1];
    EDGES[0] = EDGE_to_sweep_along;
    result = api_make_ewire( 1, EDGES, WIRE_ptr);
    if (!result.ok())
    {
      AcisQueryEngine::instance()->ACIS_API_error (result);
      PRINT_ERROR( "Unable to make ACIS wire body from curve to sweep along\n" );
      for( j=FACES_to_sweep.size(); j--; )
        api_delent( FACES_to_sweep.get_and_step() );
      return CUBIT_FAILURE;
    }
  }

  bool volume_neg;
  CubitStatus status;

  //sweep each surface
  DLIList<BODY*> swept_BODYs;
  for( i=FACES_to_sweep.size(); i--; )
  {
    FACE_to_sweep = FACES_to_sweep.get_and_step();

    if( sweep_perp )
      status = sweep_FACE_perpendicular( FACE_to_sweep, tmp_sweep_vector.length(), volume_neg,
                                         outward, 0.0, 0, false, stop_FACE );
    else if( WIRE_ptr )
      status = sweep_FACE_along_WIRE( FACE_to_sweep, WIRE_ptr, 0.0, 0, false, stop_FACE );
    else
      status = sweep_FACE_along_vector( FACE_to_sweep, tmp_sweep_vector, volume_neg, true, 0.0, 0,
                                        false, stop_FACE );

    if( status == CUBIT_FAILURE )
    {
      api_delent( FACE_to_sweep );
      for(i=swept_BODYs.size(); i--;)
        api_delent( swept_BODYs.get_and_step() );
      if( WIRE_ptr )
        api_delent( WIRE_ptr );
      return CUBIT_FAILURE;
    }
    BODY *swept_BODY = AcisQueryEngine::instance()->get_BODY_of_ENTITY(FACE_to_sweep);
    swept_BODYs.append( swept_BODY );
  }

  //clean up
  if( WIRE_ptr )
    api_delent( WIRE_ptr );

  //if there are more than 1, unite them all
  swept_BODYs.reset();
  BODY *tool_BODY = swept_BODYs.extract();
  if( swept_BODYs.size() > 0 )
  {
    for(i=swept_BODYs.size(); i--;)
    {
      outcome result = api_unite( swept_BODYs.get(), tool_BODY );
      if( !result.ok() )
      {
        PRINT_ERROR( "Problems uniting swept surfaces for webcut\n" );
        AcisQueryEngine::instance()->ACIS_API_error (result);
        api_delent( tool_BODY );
        for(i=swept_BODYs.size(); i--;  )
        {
          if( swept_BODYs.get() )
            api_delent( swept_BODYs.get_and_step() );
        }
        return CUBIT_FAILURE;
      }
      else
        swept_BODYs.change_to(NULL);
      swept_BODYs.step();
    }
  }

  //remove unnecessary edges
  api_clean_entity( (ENTITY*)tool_BODY);

  //trim tool_body if 'up_to_next_surface' option is specified
  if( up_to_next)
  {
    //get a vertex on the original face that we're sweeping
    ENTITY_LIST vertices;
    api_get_vertices( surf_ACIS->get_FACE_ptr(), vertices );

    VERTEX *surf_VERT = NULL;
    if( vertices.count() )
      surf_VERT = static_cast<VERTEX*>( vertices[0] );

    if( !trim_up_to_next_surface( tool_BODY, blank_bodies, surf_VERT ) )
    {
      api_delent( tool_BODY );
      PRINT_ERROR("Cannot use 'up_to_next' option with specified geometry\n");
      PRINT_INFO("Try the 'stop surface <id>' option instead\n");
      return CUBIT_FAILURE;
    }
  }

  return webcut( blank_bodies, tool_BODY, results_list, imprint);

}

CubitStatus AcisModifyEngine::trim_up_to_next_surface( BODY *&tool_BODY,
                                                       DLIList<BodySM*> &blank_bodies,
                                                       VERTEX *surf_VERT )
{
  outcome result;

  //make a copy of all the blank BODYs and unite them
  BODY *blank_BODY = AcisQueryEngine::get_BODY( blank_bodies.get_and_step() );
  BODY *master_blank_BODY= this->copy_BODY(blank_BODY, true );

  int i;
  for( i=blank_bodies.size()-1; i--; )
  {
    blank_BODY = AcisQueryEngine::get_BODY( blank_bodies.get_and_step() );
    BODY *copy_blank_BODY = this->copy_BODY(blank_BODY, true );
    result = api_unite( copy_blank_BODY, master_blank_BODY);

    if( !result.ok() )
    {
      PRINT_ERROR("Problems uniting blank volumes with 'up_to_next' option.\n");
      return CUBIT_FAILURE;
    }
  }

  //get all the cells
  generic_graph *boolean_graph;
  api_selective_boolean_stage1( master_blank_BODY, tool_BODY, boolean_graph );

  //find all cells that are common to tool and blank (with kind 0 AND kind 1)
  int num_verts = boolean_graph->number_of_vertices();
  gvertex **vertex_array = boolean_graph->get_vertices( num_verts );
  ENTITY_LIST common_vertices;
  for(i=num_verts; i--; )
  {
    if( vertex_array[i]->is_kind(0) && vertex_array[i]->is_kind(1) )
      common_vertices.add( vertex_array[i]->get_entity() );
  }

  ENTITY *vertex_to_keep = NULL;

  //Want the closest cell to point on swept surface....'surf_VERTEX'
  if( common_vertices.count() == 1 )
    vertex_to_keep = common_vertices[0];
  else
  {
    SPAposition point_on_surf = surf_VERT->geometry()->coords();
    SPAposition dummy_pos;
    double min_distance = -1, tmp_distance = 0;
    for( i=common_vertices.count(); i--; )
    {
      if( is_CELL( common_vertices[i] ) )
      {
        CELL *tmp_cell = NULL;
        tmp_cell = static_cast<CELL*>( common_vertices[i] );
        BODY *tmp_body = NULL;
        api_ct_copy_cell( tmp_cell, tmp_body );

        api_entity_point_distance( tmp_body, point_on_surf, dummy_pos, tmp_distance );
        api_delent( tmp_body );

        if( min_distance < 0 )
        {
          min_distance = tmp_distance;
          vertex_to_keep = common_vertices[i];
        }
        else if( tmp_distance < min_distance )
        {
          min_distance = tmp_distance;
          vertex_to_keep = common_vertices[i];
        }
      }
    }
  }

  ENTITY_LIST ent_list;
  if( vertex_to_keep )
    ent_list.add( vertex_to_keep );

  api_selective_boolean_stage2( master_blank_BODY, ent_list );

  tool_BODY = master_blank_BODY;

  if( !vertex_to_keep )
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Sweep the input FACE along the specified SPAvector
//                 a distance equal to the length of the SPAvector.
//
// Special Notes : volume_is_negative is set to CUBIT_TRUE if the volume
//                 of the resulting BODY becomes negative after the
//                 sweep operation.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date :
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::sweep_FACE_along_vector(
    FACE *& FACE_ptr,
    const CubitVector& sweep_vector,
    bool& volume_is_negative,
    bool primary_side,
    double draft_angle,
    int draft_type,
    bool rigid,
    FACE *to_FACE) const
{
     // Make sure the distance of the sweep is not zero
   if ( sweep_vector.length() < (1.0E5 * CUBIT_DBL_MIN) )
   {
      PRINT_WARNING("Distance of sweep is almost zero. No sweep was performed.\n");
      return CUBIT_FAILURE;
   }

     // Make sure the draft angle is reasonable
   if (draft_angle >= CUBIT_PI/2.0 || draft_angle <= -(CUBIT_PI/2.0))
   {
      PRINT_ERROR("Draft angle cannot be greater than %e radians\n",
                  CUBIT_PI/2.0);
      return CUBIT_FAILURE;
   }

     // Make sure the draft type is either 0 or 1 (this is an ACISism, of course!)
   if ( !(draft_type == 0 || draft_type == 1) )
   {
      PRINT_ERROR("Draft type must be either 0 or 1\n");
      return CUBIT_FAILURE;
   }

     // Sweep ho!!
   CubitStatus sweep_status = CUBIT_SUCCESS;
   SPAvector ACIS_sweep_vector;
   ACIS_sweep_vector.set_x(sweep_vector.x());
   ACIS_sweep_vector.set_y(sweep_vector.y());
   ACIS_sweep_vector.set_z(sweep_vector.z());

   logical ACIS_primary_side = CUBIT_TRUE;
   if (!primary_side)
      ACIS_primary_side = CUBIT_FALSE;

     // Turn off the merge option in ACIS so FACEs used to sweep won't be
     // "simplified" in their topology (i.e., 5-sided squared faces won't
     // be reduced to 4-sided square faces) -- in other words, turn on
     // non-regularized booleans.
   outcome result = api_set_int_option ( "merge",
                                 CUBIT_FALSE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_FACE_along_vector.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

   sweep_options ACIS_sweep_options;
   ACIS_sweep_options.set_draft_angle(draft_angle);
   ACIS_sweep_options.set_gap_type(draft_type);
   ACIS_sweep_options.set_which_side(ACIS_primary_side);
   ACIS_sweep_options.set_rigid(rigid);

   FACE *FACE_to_sweep_to = NULL;
   if( to_FACE )
   {
     //get face type
     SURFACE *SURFACE_ptr = to_FACE->geometry();
     int type = SURFACE_ptr->identity();

     //if face is planar...vector must not be parallel to surf normal
     if (type == PLANE_TYPE)
     {
        FACE_to_sweep_to = make_FACE( to_FACE, true );
     }
     //vector must intersect surf
     else
     {
        FACE_to_sweep_to = make_FACE( to_FACE );
     }

     SPAtransf trans = get_face_trans(FACE_to_sweep_to);
     surface *eq_surf = FACE_to_sweep_to->geometry()->trans_surface(trans );
     ACIS_sweep_options.set_to_face( eq_surf );
   }

   BODY* BODYptr = NULL;
   result = api_sweep_with_options (FACE_ptr, ACIS_sweep_vector,
                                    &ACIS_sweep_options, BODYptr);

   if (BODYptr)
   {
     //when setting the sweep to end at a face, result is BODYptr.
     if( to_FACE )
     {
       //swept body & to_FACE must intersect some
       BODY *tmp_body = AcisQueryEngine::instance()->get_BODY_of_ENTITY( FACE_to_sweep_to );
       if ( !BODYs_interfering( tmp_body, BODYptr ) )
       {
        PRINT_ERROR("Terminating surface must intersect swept tool volume.\n");
        api_delent( tmp_body );
        return CUBIT_FAILURE;
       }
       FACE_ptr = BODYptr->lump()->shell()->first_face();
     }
     else
     {
        PRINT_ERROR("Error during ACIS sweep operation, api_sweep_with_options\n");
        PRINT_ERROR("API wants to create a new volume.\n");
        api_delent( BODYptr );
        sweep_status = CUBIT_FAILURE;
     }
   }

   if (!result.ok())
   {
     if(rigid)
        PRINT_WARNING("Because surface is not planar, the 'rigid' keyword will be ignored.\n");

     PRINT_DEBUG_99("api_sweep_with_options failed, reverting to api_sw_face_vec.\n");

     result = api_sw_face_vec (FACE_ptr, ACIS_primary_side, ACIS_sweep_vector,
                               draft_angle, draft_type);

     if (!result.ok())
     {

       PRINT_ERROR("Error during ACIS sweep operation, api_sweep_with_options\n");
       PRINT_INFO("in sweep_FACE_along_vector\n");
       AcisQueryEngine::instance()->ACIS_API_error(result);
       sweep_status = CUBIT_FAILURE;
     }
   }

   if( FACE_to_sweep_to )
     api_delent( FACE_to_sweep_to );

   if (!result.ok())
   {
       PRINT_ERROR("Error during ACIS sweep operation, api_sweep_with_options\n");
       PRINT_INFO("in sweep_FACE_along_vector\n");
       AcisQueryEngine::instance()->ACIS_API_error(result);
       sweep_status = CUBIT_FAILURE;
   }

     // Now set the option back to its default
   result = api_set_int_option ( "merge",
                                 CUBIT_TRUE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_FACE_along_vector.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

     // Finally, check to make sure that the volume of the swept BODY is positive.
     // Note that ACIS allows the creation of swept solids that have
     // negative volume (all the surface normals point inward -- i.e.,
     // the swept solid is topologically identical to a void)
   if ( AcisQueryEngine::instance()->volume(AcisQueryEngine::instance()->get_BODY_of_ENTITY(FACE_ptr)) < 0.0)
   {
      sweep_status = CUBIT_FAILURE;
      volume_is_negative = true;

        // Set the value of result such that the API_END macro performs
        // a bulletin board rollback operation to the corresponding API_BEGIN
        //result = outcome(API_FAILED);
   }

   return sweep_status;
}

CubitStatus AcisModifyEngine::sweep_FACE_perpendicular(
    FACE *& FACE_ptr,
    double distance,
    bool& volume_is_negative,
    bool primary_side,
    double draft_angle,
    int draft_type,
    bool rigid,
    FACE *to_FACE) const
{
     // Make sure the distance of the sweep is not zero
   if ( distance < (1.0E5 * CUBIT_DBL_MIN) )
   {
      PRINT_WARNING("Distance of sweep is almost zero. No sweep was performed.\n");
      return CUBIT_FAILURE;
   }

     // Make sure the draft angle is reasonable
   if (draft_angle >= CUBIT_PI/2.0 || draft_angle <= -(CUBIT_PI/2.0))
   {
      PRINT_ERROR("Draft angle cannot be greater than %e radians\n",
                  CUBIT_PI/2.0);
      return CUBIT_FAILURE;
   }

     // Make sure the draft type is either 0 or 1 (this is an ACISism, of course!)
   if ( !(draft_type == 0 || draft_type == 1) )
   {
      PRINT_ERROR("Draft type must be either 0 or 1\n");
      return CUBIT_FAILURE;
   }

     // Sweep ho!!
   CubitStatus sweep_status = CUBIT_SUCCESS;

   logical ACIS_primary_side = CUBIT_FALSE;
   if (!primary_side)
      ACIS_primary_side = CUBIT_TRUE;

     // Turn off the merge option in ACIS so FACEs used to sweep won't be
     // "simplified" in their topology (i.e., 5-sided squared faces won't
     // be reduced to 4-sided square faces) -- in other words, turn on
     // non-regularized booleans.
   outcome result = api_set_int_option ( "merge",
                                 CUBIT_FALSE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_FACE_perpendicular.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

   sweep_options ACIS_sweep_options;
   ACIS_sweep_options.set_draft_angle(draft_angle);
   ACIS_sweep_options.set_gap_type(draft_type);
   ACIS_sweep_options.set_which_side(ACIS_primary_side);
   ACIS_sweep_options.set_rigid(rigid);

   FACE *FACE_to_sweep_to = NULL;
   if( to_FACE )
   {
     //get face type
     SURFACE *SURFACE_ptr = to_FACE->geometry();
     int type = SURFACE_ptr->identity();

     //if face is planar...vector must not be parallel to surf normal
     if (type == PLANE_TYPE)
     {
        FACE_to_sweep_to = make_FACE( to_FACE, true );
     }
     //vector must intersect surf
     else
     {
        FACE_to_sweep_to = make_FACE( to_FACE );
     }

     SPAtransf trans = get_face_trans(FACE_to_sweep_to);
     surface *eq_surf = FACE_to_sweep_to->geometry()->trans_surface(trans );
     ACIS_sweep_options.set_to_face( eq_surf );
   }

   BODY* BODYptr = NULL;
   result = api_sweep_with_options (FACE_ptr, distance,
                                    &ACIS_sweep_options, BODYptr);

   if (BODYptr)
   {
     //when setting the sweep to end at a face, result is BODYptr.
     if( to_FACE)
     {
       //swept body & to_FACE must intersect some
       BODY *tmp_body = AcisQueryEngine::instance()->get_BODY_of_ENTITY( FACE_to_sweep_to );
       if ( !BODYs_interfering( tmp_body, BODYptr ) )
       {
        PRINT_ERROR("Terminating surface must intersect swept tool volume.\n");
        api_delent( tmp_body );
        return CUBIT_FAILURE;
       }
       FACE_ptr = BODYptr->lump()->shell()->first_face();
     }
     else
     {
        PRINT_ERROR("Error during ACIS sweep operation, api_sweep_with_options\n");
        PRINT_ERROR("API wants to create a new volume \n");
        api_delent( BODYptr );
        sweep_status = CUBIT_FAILURE;
     }
   }

   if( FACE_to_sweep_to )
     api_delent( FACE_to_sweep_to );

   if (!result.ok())
   {
      PRINT_ERROR("Error during ACIS sweep operation,  api_sweep_with_options \n");
      PRINT_INFO("in sweep_FACE_perpendicular\n");
      PRINT_INFO("If face is not planar, perpendicular doesn't make sense.\n"
                 "Try sweeping along a vector instead.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);
      sweep_status = CUBIT_FAILURE;
   }

     // Now set the option back to its default
   result = api_set_int_option ( "merge",
                                 CUBIT_TRUE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_FACE_perpendicular.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

     // Finally, check to make sure that the volume of the swept BODY is positive.
     // Note that ACIS allows the creation of swept solids that have
     // negative volume (all the surface normals point inward -- i.e.,
     // the swept solid is topologically identical to a void)
   if ( AcisQueryEngine::instance()->volume(AcisQueryEngine::instance()->get_BODY_of_ENTITY(FACE_ptr)) < 0.0)
   {
      sweep_status = CUBIT_FAILURE;
      volume_is_negative = CUBIT_TRUE;

        // Set the value of result such that the API_END macro performs
        // a bulletin board rollback operation to the corresponding API_BEGIN
        //result = outcome(API_FAILED);
   }

   return sweep_status;
}

CubitStatus AcisModifyEngine::sweep_EDGE_along_vector(
    EDGE* EDGE_ptr,
    BODY *&BODYptr,
    const CubitVector& sweep_vector,
    double draft_angle,
    int draft_type,
    bool rigid,
    FACE *to_FACE) const
{
     // Make sure the distance of the sweep is not zero
   if ( sweep_vector.length() < (1.0E5 * CUBIT_DBL_MIN) )
   {
      PRINT_WARNING("Distance of sweep is almost zero. No sweep was performed.\n");
      return CUBIT_FAILURE;
   }

     // Make sure the draft angle is reasonable
   if (draft_angle >= CUBIT_PI/2.0 || draft_angle <= -(CUBIT_PI/2.0))
   {
      PRINT_ERROR("Draft angle cannot be greater than %e radians\n",
                  CUBIT_PI/2.0);
      return CUBIT_FAILURE;
   }

     // Make sure the draft type is either 0 or 1 (this is an ACISism, of course!)
   if ( !(draft_type == 0 || draft_type == 1) )
   {
      PRINT_ERROR("Draft type must be either 0 or 1\n");
      return CUBIT_FAILURE;
   }

     // Sweep ho!!
   CubitStatus sweep_status = CUBIT_SUCCESS;
   SPAvector ACIS_sweep_vector;
   ACIS_sweep_vector.set_x(sweep_vector.x());
   ACIS_sweep_vector.set_y(sweep_vector.y());
   ACIS_sweep_vector.set_z(sweep_vector.z());

     // Turn off the merge option in ACIS so EDGEs used to sweep won't be
     // "simplified" in their topology (i.e., 5-sided squared faces won't
     // be reduced to 4-sided square faces) -- in other words, turn on
     // non-regularized booleans.
   outcome result = api_set_int_option ( "merge",
                                 CUBIT_FALSE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_EDGE_along_vector.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

   sweep_options ACIS_sweep_options;
   ACIS_sweep_options.set_draft_angle(draft_angle);
   ACIS_sweep_options.set_gap_type(draft_type);
   ACIS_sweep_options.set_rigid(rigid);
   EDGE *copied_EDGE_ptr = NULL;
   result = api_edge( EDGE_ptr, copied_EDGE_ptr);

   FACE *FACE_to_sweep_to = NULL;
   if( to_FACE )
   {
     //get face type
     SURFACE *SURFACE_ptr = to_FACE->geometry();
     int type = SURFACE_ptr->identity();

     //if face is planar...vector must not be parallel to surf normal
     if (type == PLANE_TYPE)
     {
        FACE_to_sweep_to = make_FACE( to_FACE, true );
     }
     //vector must intersect surf
     else
     {
        FACE_to_sweep_to = make_FACE( to_FACE );
     }

     SPAtransf trans = get_face_trans(FACE_to_sweep_to);
     surface *eq_surf = FACE_to_sweep_to->geometry()->trans_surface(trans );
     ACIS_sweep_options.set_to_face( eq_surf );
   }

   if( !result.ok() )
   {
     PRINT_ERROR( "Unable to copy edge\n");
     return CUBIT_FAILURE;
   }
   ATTRIB_CUBIT_OWNER::remove_cubit_owner( (ENTITY *)copied_EDGE_ptr,
                                           CUBIT_TRUE );
   result = api_sweep_with_options (copied_EDGE_ptr, ACIS_sweep_vector,
                                    &ACIS_sweep_options, BODYptr);

   if (BODYptr)
   {
       sweep_status = CUBIT_SUCCESS;
   }
   if (!result.ok())
   {
      PRINT_ERROR("Error during ACIS sweep operation, api_sweep_with_options\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);
      if ( copied_EDGE_ptr != NULL )
        api_delent(copied_EDGE_ptr);
      sweep_status = CUBIT_FAILURE;
   }

     // Now set the option back to its default
   result = api_set_int_option ( "merge",
                                 CUBIT_TRUE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_EDGE_along_vector.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

   if( FACE_to_sweep_to )
     api_delent( FACE_to_sweep_to );

     // Finally, check to make sure that the volume of the swept BODY is positive.
     // Note that ACIS allows the creation of swept solids that have
     // negative volume (all the surface normals point inward -- i.e.,
     // the swept solid is topologically identical to a void)

   return sweep_status;
}


//-------------------------------------------------------------------------
// Purpose       : Sweep the input FACE along the specified EDGE
//
// Special Notes : Due to improved checking before sweeping, the
//                 volume_is_negative flag has been removed (not necessary)
//
// Creator       : Eric Nielsen
//
// Creation Date : 8/7/98
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::sweep_FACE_along_WIRE(
    FACE *& FACE_ptr,
    BODY* wire_BODY_ptr,
    double draft_angle,
    int draft_type,
    bool rigid,
    FACE *to_FACE ) const
{
     // Make sure the draft angle is reasonable
   if (draft_angle >= CUBIT_PI/2.0 || draft_angle <= -(CUBIT_PI/2.0))
   {
      PRINT_ERROR("Draft angle cannot be greater than %e radians\n",
                  CUBIT_PI/2.0);
      return CUBIT_FAILURE;
   }

   CubitStatus sweep_status = CUBIT_SUCCESS;

   sweep_options ACIS_sweep_options;

     // Turn off the merge option in ACIS so FACEs used to sweep won't be
     // "simplified" in their topology (i.e., 5-sided squared faces won't
     // be reduced to 4-sided square faces) -- in other words, turn on
     // non-regularized booleans.
   outcome result = api_set_int_option ( "merge", CUBIT_FALSE );

   if (!result.ok())
   {
     PRINT_WARNING("In AcisModifyEngine::sweep_FACE_along_WIRE.\n"
                   "         Error in api_set_int_option.\n");
     AcisQueryEngine::instance()->ACIS_API_error(result);

       // As the problem is not really serious, continue and make sure that
       // the value of result is set to "OK" so no rollback occurs
     result = outcome(0);
   }

   ACIS_sweep_options.set_draft_angle(draft_angle);
   ACIS_sweep_options.set_gap_type(draft_type);
   ACIS_sweep_options.set_rigid(rigid);

   FACE *FACE_to_sweep_to = NULL;
   if( to_FACE )
   {
     //get face type
     SURFACE *SURFACE_ptr = to_FACE->geometry();
     int type = SURFACE_ptr->identity();

     //if face is planar...vector must not be parallel to surf normal
     if (type == PLANE_TYPE)
     {
        FACE_to_sweep_to = make_FACE( to_FACE, true );
     }
     //vector must intersect surf
     else
     {
        FACE_to_sweep_to = make_FACE( to_FACE );
     }

     SPAtransf trans = get_face_trans(FACE_to_sweep_to);
     surface *eq_surf = FACE_to_sweep_to->geometry()->trans_surface(trans );
     ACIS_sweep_options.set_to_face( eq_surf );
   }

   // Do the real work!!
   BODY* BODYptr = NULL;
   result = api_sweep_with_options (FACE_ptr, (EDGE*)wire_BODY_ptr, &ACIS_sweep_options, BODYptr);

   if (BODYptr)
   {
     //when setting the sweep to end at a face, result is BODYptr.
     if( to_FACE )
     {
       //swept body & to_FACE must intersect some
       BODY *tmp_body = AcisQueryEngine::instance()->get_BODY_of_ENTITY( FACE_to_sweep_to );
       if ( !BODYs_interfering( tmp_body, BODYptr ) )
       {
        PRINT_ERROR("Terminating surface must intersect swept tool volume.\n");
        api_delent( tmp_body );
        return CUBIT_FAILURE;
       }
       FACE_ptr = BODYptr->lump()->shell()->first_face();
     }
     else
     {
        PRINT_ERROR("Error during ACIS sweep operation, api_sweep_with_options\n");
        PRINT_ERROR("API wants to create a new volume \n");
        api_delent( BODYptr );
        sweep_status = CUBIT_FAILURE;
     }
   }

   if( FACE_to_sweep_to )
     api_delent( FACE_to_sweep_to );

   if (!result.ok())
   {
      PRINT_ERROR("Error during ACIS sweep operation, api_sweep_with_options\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);
      sweep_status = CUBIT_FAILURE;
   }

     // Now set the option back to its default
   result = api_set_int_option ( "merge", CUBIT_TRUE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_FACE_along_WIRE.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

   return sweep_status;

}

CubitStatus AcisModifyEngine::sweep_EDGE_along_WIRE(
    EDGE* EDGE_ptr,
    BODY* wire_BODY_ptr,
    BODY *&BODYptr,
    double draft_angle,
    int draft_type,
    bool rigid,
    FACE *to_FACE) const
{

     // Make sure the draft angle is reasonable
   if (draft_angle >= CUBIT_PI/2.0 || draft_angle <= -(CUBIT_PI/2.0))
   {
      PRINT_ERROR("Draft angle cannot be greater than %e radians\n",
                  CUBIT_PI/2.0);
      return CUBIT_FAILURE;
   }

   CubitStatus sweep_status = CUBIT_SUCCESS;

   sweep_options ACIS_sweep_options;

     // Turn off the merge option in ACIS so FACEs used to sweep won't be
     // "simplified" in their topology (i.e., 5-sided squared faces won't
     // be reduced to 4-sided square faces) -- in other words, turn on
     // non-regularized booleans.
   outcome result = api_set_int_option ( "merge", CUBIT_FALSE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_EDGE_along_WIRE.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

   ACIS_sweep_options.set_draft_angle(draft_angle);
   ACIS_sweep_options.set_gap_type(draft_type);
   ACIS_sweep_options.set_rigid(rigid);
   EDGE *copied_EDGE_ptr = NULL;
   result = api_edge( EDGE_ptr, copied_EDGE_ptr );
   if( !result.ok() )
   {
     PRINT_ERROR( "Unable to copy edge\n");
     return CUBIT_FAILURE;
   }

   //if edge-to-sweep has start or end vertex coincident w/
   //end vertex on path wire, flip all edges in path wire
   SPAposition end_coords;
   get_wire_end_position(wire_BODY_ptr, end_coords);

   SPAposition tmp_coords = EDGE_ptr->start()->geometry()->coords();
   CubitVector start_coord( tmp_coords.x(), tmp_coords.y(), tmp_coords.z() );

   tmp_coords = EDGE_ptr->end()->geometry()->coords();
   CubitVector end_coord( tmp_coords.x(), tmp_coords.y(), tmp_coords.z() );

   CubitVector wire_end_coord( end_coords.x(), end_coords.y(), end_coords.z() );

   if ( (GeometryQueryTool::instance()->about_spatially_equal(start_coord, wire_end_coord,
                                       GeometryQueryTool::get_geometry_factor())) ||
        (GeometryQueryTool::instance()->about_spatially_equal(end_coord, wire_end_coord,
                                       GeometryQueryTool::get_geometry_factor())) )
   {
     //reverse all edges in wire body
     ENTITY_LIST EDGE_list;
     api_get_edges( wire_BODY_ptr, EDGE_list );

     EDGE_list.init();
     ENTITY *tmp_ent = NULL;
     while( (tmp_ent = EDGE_list.next()) != NULL )
     {
       EDGE *tmp_edge = (EDGE*)tmp_ent;
       VERTEX *old_start = tmp_edge->start();
       VERTEX *old_end = tmp_edge->end();
       tmp_edge->set_sense( tmp_edge->sense() == FORWARD ? REVERSED : FORWARD );
       tmp_edge->set_start( old_end );
       tmp_edge->set_end( old_start );
     }
   }

   ATTRIB_CUBIT_OWNER::remove_cubit_owner( (ENTITY *)copied_EDGE_ptr,
                                           CUBIT_TRUE );

   FACE *FACE_to_sweep_to = NULL;
   if( to_FACE )
   {
     //get face type
     SURFACE *SURFACE_ptr = to_FACE->geometry();
     int type = SURFACE_ptr->identity();

     //if face is planar...vector must not be parallel to surf normal
     if (type == PLANE_TYPE)
     {
        FACE_to_sweep_to = make_FACE( to_FACE, true );
     }
     //vector must intersect surf
     else
     {
        FACE_to_sweep_to = make_FACE( to_FACE );
     }

     SPAtransf trans = get_face_trans(FACE_to_sweep_to);
     surface *eq_surf = FACE_to_sweep_to->geometry()->trans_surface(trans );
     ACIS_sweep_options.set_to_face( eq_surf );
   }

     // Do the real work!!
   result = api_sweep_with_options (copied_EDGE_ptr, wire_BODY_ptr, &ACIS_sweep_options, BODYptr);

   if (BODYptr)
   {
      sweep_status = CUBIT_SUCCESS;
   }

   if (!result.ok())
   {
      PRINT_ERROR("Error during ACIS sweep operation, api_sweep_with_options\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);
      if ( copied_EDGE_ptr != NULL )
        api_delent(copied_EDGE_ptr);

      sweep_status = CUBIT_FAILURE;
   }

   if( FACE_to_sweep_to )
     api_delent( FACE_to_sweep_to );

     // Now set the option back to its default
   result = api_set_int_option ( "merge", CUBIT_TRUE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_EDGE_along_WIRE.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

   return sweep_status;
}

//-------------------------------------------------------------------------
// Purpose       : Sweep the input FACE about the specified axis,
//                 by the specified angle.
//
// Special Notes : volume_is_negative is set to CUBIT_TRUE if the volume
//                 of the resulting BODY becomes negative after the
//                 sweep operation.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date :
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::sweep_FACE_about_axis(
    FACE *&FACE_ptr,
    const CubitVector& axis_unit_vector,
    const CubitVector& point,
    double angle,
    bool &volume_is_negative,
    bool primary_side,
    int steps,
    double draft_angle,
    int draft_type,
    bool rigid,
    FACE *to_FACE ) const
{
     // Sweep ho!!
   CubitStatus sweep_status = CUBIT_SUCCESS;

     // Convert the Cubit point and SPAvector into an ACIS point and SPAvector
     // before calling the member function that calls the ACIS API
     // function.
   SPAunit_vector ACIS_axis_unit_vector(
       axis_unit_vector.x(),
       axis_unit_vector.y(),
       axis_unit_vector.z() );
   SPAposition ACIS_point( point.x(), point.y(), point.z() );
   logical ACIS_primary_side = CUBIT_TRUE;
   if (!primary_side) ACIS_primary_side = CUBIT_FALSE;

     // Make sure the distance of the sweep is not zero
   if (angle == 0.0)
   {
      PRINT_ERROR("The angle of the sweep must be > 0.0. Nothing done.\n");
      return CUBIT_FAILURE;
   }

     // Make sure the draft angle is reasonable
   if (draft_angle >= CUBIT_PI/2.0 || draft_angle <= -(CUBIT_PI/2.0))
   {
      PRINT_ERROR("Draft angle cannot be greater than %e radians\n",
                  CUBIT_PI/2.0);
      return CUBIT_FAILURE;
   }

     // Make sure the draft type is either 0 or 1 (this is an ACISism, of course!)
   if ( !(draft_type == 0 || draft_type == 1) )
   {
      PRINT_ERROR("Draft type must be either 0 or 1\n");
      return CUBIT_FAILURE;
   }

     // Make sure we have a reasonable number of steps
   if (steps < 0)
   {
      PRINT_WARNING("Number of steps must be a positive number.\n"
                    "         It is set to 0 and the sweep operation will proceed.\n");
      steps = 0;
   }

     // Turn off the merge option in ACIS so FACEs used to sweep won't be
     // "simplified" in their topology (i.e., 5-sided squared faces won't
     // be reduced to 4-sided square faces) -- in other words, turn on
     // non-regularized booleans
   outcome result = api_set_int_option ( "merge",
                                         CUBIT_FALSE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_FACE_about_axis.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

   sweep_options ACIS_sweep_options;
   ACIS_sweep_options.set_draft_angle(draft_angle);
   ACIS_sweep_options.set_gap_type(draft_type);
   ACIS_sweep_options.set_which_side(ACIS_primary_side);
   ACIS_sweep_options.set_steps(steps);
   ACIS_sweep_options.set_sweep_angle(angle);
   ACIS_sweep_options.set_rigid(rigid);

   FACE *FACE_to_sweep_to = NULL;
   if( to_FACE )
   {
     //get face type
     SURFACE *SURFACE_ptr = to_FACE->geometry();
     int type = SURFACE_ptr->identity();

     //if face is planar...vector must not be parallel to surf normal
     if (type == PLANE_TYPE)
     {
        FACE_to_sweep_to = make_FACE( to_FACE, true );
     }
     //vector must intersect surf
     else
     {
        FACE_to_sweep_to = make_FACE( to_FACE );
     }

     SPAtransf trans = get_face_trans(FACE_to_sweep_to);
     surface *eq_surf = FACE_to_sweep_to->geometry()->trans_surface(trans );
     ACIS_sweep_options.set_to_face( eq_surf );
   }

   BODY* face_body_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY(FACE_ptr);
   BODY* BODYptr = NULL;
   result = api_sweep_with_options(FACE_ptr,
                                  ACIS_point, ACIS_axis_unit_vector,
                                  &ACIS_sweep_options, BODYptr);

   if (BODYptr)
   {
     //when setting the sweep to end at a face, result is BODYptr.
     if( to_FACE )
     {
       //swept body & to_FACE must intersect some
       BODY *tmp_body = AcisQueryEngine::instance()->get_BODY_of_ENTITY( FACE_to_sweep_to );
       if ( !BODYs_interfering( tmp_body, BODYptr ) )
       {
        PRINT_ERROR("Terminating surface must intersect swept tool volume.\n");
        api_delent( tmp_body );
        return CUBIT_FAILURE;
       }
       FACE_ptr = BODYptr->lump()->shell()->first_face();
     }
     else
     {
        PRINT_ERROR("Error during ACIS sweep operation, api_sweep_with_options\n");
        PRINT_ERROR("API wants to create a new volume \n");
        api_delent( BODYptr );
        sweep_status = CUBIT_FAILURE;
     }
   }

   if( FACE_to_sweep_to )
     api_delent( FACE_to_sweep_to );

   if (!result.ok())
   {
      PRINT_ERROR("Error during ACIS sweep operation, api_sw_face_axis\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);
      sweep_status = CUBIT_FAILURE;
   }

     // Now set the option back to its default
   result = api_set_int_option ( "merge",
                                 CUBIT_TRUE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_FACE_about_axis.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

     // Finally, check to make sure that the volume of the swept BODY is positive.
     // Note that ACIS allows the creation of swept solids that have
     // negative volume (all the surface normals point inward -- i.e.,
     // the swept solid is topologically identical to a void)
   BODY *new_body = NULL;
   if(BODYptr)
   {
      new_body = BODYptr;
   }
   else if(face_body_ptr)
   {
      new_body = face_body_ptr;
   }

   if(new_body)
   {
      if ( AcisQueryEngine::instance()->volume(new_body) < 0.0)
      {
         sweep_status = CUBIT_FAILURE;
         volume_is_negative = CUBIT_SUCCESS;

         // Set the value of result such that the API_END macro performs
         // a bulletin board rollback operation to the corresponding API_BEGIN
         //result = outcome(API_FAILURE);
      }
   }

   return sweep_status;
}

CubitStatus AcisModifyEngine::sweep_EDGE_about_axis(
    EDGE* EDGE_ptr,
    BODY *&BODYptr,
    const CubitVector& axis_unit_vector,
    const CubitVector& point,
    double angle,
    int steps,
    double draft_angle,
    int draft_type,
    bool make_solid,
    bool rigid,
    FACE *to_FACE) const
{
     // Sweep ho!!
   CubitStatus sweep_status = CUBIT_SUCCESS;

     // Convert the Cubit point and SPAvector into an ACIS point and SPAvector
     // before calling the member function that calls the ACIS API
     // function.
   SPAunit_vector ACIS_axis_unit_vector(
       axis_unit_vector.x(),
       axis_unit_vector.y(),
       axis_unit_vector.z() );
   SPAposition ACIS_point( point.x(), point.y(),  point.z() );
   if (angle == 0.0)
   {
      PRINT_ERROR("The angle of the sweep must be > 0.0. Nothing done.\n");
      return CUBIT_FAILURE;
   }

     // Make sure the draft angle is reasonable
   if (draft_angle >= CUBIT_PI/2.0 || draft_angle <= -(CUBIT_PI/2.0))
   {
      PRINT_ERROR("Draft angle cannot be greater than %e radians\n",
                  CUBIT_PI/2.0);
      return CUBIT_FAILURE;
   }

     // Make sure the draft type is either 0 or 1 (this is an ACISism, of course!)
   if ( !(draft_type == 0 || draft_type == 1) )
   {
      PRINT_ERROR("Draft type must be either 0 or 1\n");
      return CUBIT_FAILURE;
   }

     // Make sure we have a reasonable number of steps
   if (steps < 0)
   {
      PRINT_WARNING("Number of steps must be a positive number.\n"
                    "         It is set to 0 and the sweep operation will proceed.\n");
      steps = 0;
   }

     // Turn off the merge option in ACIS so FACEs used to sweep won't be
     // "simplified" in their topology (i.e., 5-sided squared faces won't
     // be reduced to 4-sided square faces) -- in other words, turn on
     // non-regularized booleans
   outcome result = api_set_int_option ( "merge",
                                         CUBIT_FALSE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_EDGE_about_axis.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

   sweep_options ACIS_sweep_options;
   ACIS_sweep_options.set_draft_angle(draft_angle);
   ACIS_sweep_options.set_gap_type(draft_type);
   ACIS_sweep_options.set_steps(steps);
   ACIS_sweep_options.set_sweep_angle(angle);
   ACIS_sweep_options.set_rigid(rigid);
   EDGE *copied_EDGE_ptr = NULL;
   result = api_edge( EDGE_ptr, copied_EDGE_ptr );
   if( !result.ok() )
   {
     PRINT_ERROR( "Unable to copy edge\n");
     return CUBIT_FAILURE;
   }
   ATTRIB_CUBIT_OWNER::remove_cubit_owner( (ENTITY *)copied_EDGE_ptr,
                                           CUBIT_TRUE );

   FACE *FACE_to_sweep_to = NULL;
   if( to_FACE )
   {
     //get face type
     SURFACE *SURFACE_ptr = to_FACE->geometry();
     int type = SURFACE_ptr->identity();

     //if face is planar...vector must not be parallel to surf normal
     if (type == PLANE_TYPE)
     {
        FACE_to_sweep_to = make_FACE( to_FACE, true );
     }
     //vector must intersect surf
     else
     {
        FACE_to_sweep_to = make_FACE( to_FACE );
     }

     SPAtransf trans = get_face_trans(FACE_to_sweep_to);
     surface *eq_surf = FACE_to_sweep_to->geometry()->trans_surface(trans );
     ACIS_sweep_options.set_to_face( eq_surf );
   }

   if(make_solid)
   {
     ACIS_sweep_options.set_close_to_axis(TRUE);
     EDGE* EDGES[1];
     EDGES[0] = copied_EDGE_ptr;
     BODY* EDGE_BODY_ptr = NULL;

     result = api_make_ewire(1, EDGES, EDGE_BODY_ptr);

//     WIRE* WIRE_ptr = EDGE_BODY_ptr->lump()->shell()->wire();
     result = api_sweep_with_options(EDGE_BODY_ptr,
                                     ACIS_point, ACIS_axis_unit_vector,
                                     &ACIS_sweep_options, BODYptr);
     BODYptr = EDGE_BODY_ptr;

   }
   else {
     result = api_sweep_with_options(copied_EDGE_ptr,
                                     ACIS_point, ACIS_axis_unit_vector,
                                     &ACIS_sweep_options, BODYptr);
   }


   if (BODYptr)
   {
       sweep_status = CUBIT_SUCCESS;
   }
   if (!result.ok())
   {
      PRINT_ERROR("Error during ACIS sweep operation, api_sweep_with_options\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);
      if ( copied_EDGE_ptr != NULL )
        api_delent(copied_EDGE_ptr);

      sweep_status = CUBIT_FAILURE;
   }

   if( FACE_to_sweep_to )
     api_delent( FACE_to_sweep_to );

     // Now set the option back to its default
   result = api_set_int_option ( "merge",
                                 CUBIT_TRUE );
   if (!result.ok())
   {
      PRINT_WARNING("In AcisModifyEngine::sweep_EDGE_about_axis.\n"
                    "         Error in api_set_int_option.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);

        // As the problem is not really serious, continue and make sure that
        // the value of result is set to "OK" so no rollback occurs
      result = outcome(0);
   }

     // Finally, check to make sure that the volume of the swept BODY is positive.
     // Note that ACIS allows the creation of swept solids that have
     // negative volume (all the surface normals point inward -- i.e.,
     // the swept solid is topologically identical to a void)

   return sweep_status;
}

//-------------------------------------------------------------------------
// Purpose       : Make a (planar) BODY, given the 4 input points.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date :
//-------------------------------------------------------------------------
BODY* AcisModifyEngine::make_planar_quad_BODY ( const CubitVector& V1,
                                                const CubitVector& V2,
                                                const CubitVector& V3,
                                                const CubitVector& V4 ) const
{
   GeometryQueryTool* GT_ptr = GeometryQueryTool::instance();
   BODY* Body1 = NULL;

     // Debug info
   PRINT_DEBUG_18(
               "Debug output from AcisModifyEngine::make_planar_quad_BODY\n");

     // First make sure that none of the points have the same coordinates.
   if ( (GT_ptr->about_spatially_equal(V1, V2,
                                       GeometryQueryTool::get_geometry_factor())) ||
        (GT_ptr->about_spatially_equal(V1, V3,
                                       GeometryQueryTool::get_geometry_factor())) ||
        (GT_ptr->about_spatially_equal(V1, V4,
                                       GeometryQueryTool::get_geometry_factor())) ||
        (GT_ptr->about_spatially_equal(V2, V3,
                                       GeometryQueryTool::get_geometry_factor())) ||
        (GT_ptr->about_spatially_equal(V2, V4,
                                       GeometryQueryTool::get_geometry_factor())) ||
        (GT_ptr->about_spatially_equal(V3, V4,
                                       GeometryQueryTool::get_geometry_factor())) )
   {
      PRINT_ERROR("In AcisModifyEngine::make_planar_quad_BODY.\n"
                  "       Some of the input points are spatially equivalent.\n"
                  "       Cannot create a planar quadrilateral ACIS BODY.\n");
      return NULL;
   }

     // Make sure that the 4 input points are coplanar.
     // If the 4 points are coplanar, then the scalar triple product
     // of the 3 vectors joining points 1 and 2, 2 and 3, and 3 and 4, will
     // equal zero.
     //     For 4 coplanar points, Vector12 . (Vector23 * Vector34) = 0

   CubitVector vector12 = V2 - V1;
   CubitVector vector23 = V3 - V2;
   CubitVector vector34 = V4 - V3;

     // normalize the vectors
   vector12.normalize();
   vector23.normalize();
   vector34.normalize();

   double tripleProd = vector12 % (vector23 * vector34);

     // MJP Note:
     // The following tests to make sure that tripleProd is "0.0". But this
     // cannot be tested exactly, so we have to approximate the notion of zero.
     // I have run into a case where a not-so-large (size about 2000 units)
     // set of 4 points results in a tripleProd value of approximately
     // -3.5E-6 whose fabs is 3.5E-6 which is larger than SPAresabs (1E-6).
     // I have put in an arbitrary factor of 100 here!!!!!!! This raises
     // our notion of "zero" from 1E-6 to 1E-4......:-)
   if (fabs(tripleProd) > 100.0 * GEOMETRY_RESABS)
   {
      PRINT_ERROR("In AcisModifyEngine::make_planar_quad_BODY.\n"
                  "       The 4 input points are not coplanar.\n"
                  "       Cannot create a planar quadrilateral ACIS BODY.\n");
      return NULL;
   }

     // Initialize the bulletin board (this is done so that if there is an
     // error while building this body, ACIS will automatically roll its
     // state back to what it was before this next statement was executed.
   int success = CUBIT_TRUE;
   API_BEGIN;

         // Create an ACIS Body using the direct interface (it's tons of fun :-)

         // Create ACIS points
   APOINT* P1 = new APOINT (V1.x(), V1.y(), V1.z());
   APOINT* P2 = new APOINT (V2.x(), V2.y(), V2.z());
   APOINT* P3 = new APOINT (V3.x(), V3.y(), V3.z());
   APOINT* P4 = new APOINT (V4.x(), V4.y(), V4.z());

     // Create ACIS straight lines between the 4 input points
   STRAIGHT* S1 = new STRAIGHT(      P1->coords(),
                                     normalise(P2->coords() - P1->coords()) );
   STRAIGHT* S2 = new STRAIGHT(      P2->coords(),
                                     normalise(P3->coords() - P2->coords()) );
   STRAIGHT* S3 = new STRAIGHT(      P3->coords(),
                                     normalise(P4->coords() - P3->coords()) );
   STRAIGHT* S4 = new STRAIGHT(      P4->coords(),
                                     normalise(P1->coords() - P4->coords()) );

     // Compute the unit SPAvector of the plane containing the 4 points
   CubitVector planeNormal;
   planeNormal = (vector12 * vector23);
   planeNormal.normalize();

     // Create an ACIS surface
   PLANE* Plane1 = new PLANE( P1->coords(),
                              SPAunit_vector(planeNormal.x(),
                                          planeNormal.y(),
                                          planeNormal.z() ) );

     // Create ACIS vertices
   VERTEX* Vertex1 = new VERTEX (P1);
   VERTEX* Vertex2 = new VERTEX (P2);
   VERTEX* Vertex3 = new VERTEX (P3);
   VERTEX* Vertex4 = new VERTEX (P4);

     // Create ACIS edges
   EDGE* Edge1 = new EDGE (Vertex1, Vertex2, S1, FORWARD);
   EDGE* Edge2 = new EDGE (Vertex2, Vertex3, S2, FORWARD);
   EDGE* Edge3 = new EDGE (Vertex3, Vertex4, S3, FORWARD);
   EDGE* Edge4 = new EDGE (Vertex4, Vertex1, S4, FORWARD);

     // Create an ACIS loop of coedges
   COEDGE* Coedge1 = new COEDGE (Edge1, FORWARD, NULL, NULL);
   COEDGE* Coedge2 = new COEDGE (Edge2, FORWARD, Coedge1, NULL);
   COEDGE* Coedge3 = new COEDGE (Edge3, FORWARD, Coedge2, NULL);
 /*COEDGE* Coedge4 =*/ new COEDGE (Edge4, FORWARD, Coedge3, Coedge1);
   LOOP* Loop1 = new LOOP (Coedge1, NULL);

     // Create an ACIS face
   FACE* Face1 = new FACE (Loop1, NULL, Plane1, FORWARD);

     // Create an ACIS shell
   SHELL* Shell1 = new SHELL (Face1, NULL, NULL);

     // Create an ACIS lump
   LUMP* Lump1 = new LUMP (Shell1, NULL);

     // Create the ACIS body
   Body1 = new BODY (Lump1);

     // Make sure all's well with this BODY and then 2D-it!
   if (Body1 != NULL)
   {
        // Finally, make sure this BODY consists of a DOUBLE-SIDED,
        // BOTH_OUTSIDE FACE
      result = api_body_to_2d(Body1);
      if (!result.ok())
      {
         PRINT_ERROR("In AcisModifyEngine::make_planar_quad_BODY.\n"
                     "  Unable to convert FACE to a 2D BODY using api_body_to2d.\n");
         success = CUBIT_FALSE;
      }
   }

     // Scream bloody murder....only if necessary, of course :-)
   if (Body1 == NULL)
   {
      PRINT_ERROR("In AcisModifyEngine::make_planar_quad_BODY.\n"
                  "       Problem creating a Planar Quadrilateral ACIS BODY.\n");
      success = CUBIT_FALSE;
   }

     // Terminate the bulletin board
   API_END;

   if (success)
   {
       // Debug info
     PRINT_DEBUG_18("  Successfully created a planar quadrilateral ACIS BODY.\n");
     return Body1;
   }

   else
   {
     return NULL;
   }
}

//-------------------------------------------------------------------------
// Purpose       : Make a brick (cuboid), given the location, orientation
//                 and size.
//
// Special Notes : Extension is equivalent to 1/2 width, depth, height.
//
// Creator       : Steve Storm
//
// Creation Date : 10/09/00
//-------------------------------------------------------------------------
BODY*
AcisModifyEngine::make_brick_BODY( const CubitVector &center,
                                   const CubitVector axes[3],
                                   const CubitVector &extension ) const
{
   BODY* BODY_ptr = NULL;

   double width = 2.0*extension.x();
   double depth = 2.0*extension.y();
   double height = 2.0*extension.z();

     // Create a BODY that represents the cuboid
   outcome result = api_make_cuboid ( width, depth, height, BODY_ptr );
   if (!result.ok())
   {
      PRINT_ERROR("In AcisGeometryEngine::make_brick_BODY, Line# %d\n"
                  "     Problems creating a brick\n", __LINE__);
      return NULL ;
   }

   // Transform the brick into SPAposition

   SPAposition brick_center( center.x(), center.y(), center.z() );
   SPAunit_vector x_axis_this( axes[0].x(), axes[0].y(), axes[0].z() );
   SPAunit_vector y_axis_this( axes[1].x(), axes[1].y(), axes[1].z() );

   SPAtransf transformation = coordinate_transf( brick_center, x_axis_this, y_axis_this );

   // Concatenate the input transformation to the existing one
   result = api_apply_transf ( BODY_ptr, transformation );
   if (!result.ok())
   {
      AcisQueryEngine::instance()->ACIS_API_error ( result );
      PRINT_ERROR("In AcisGeometryEngine::make_brick_BODY, Line# %d\n"
                  "     Problems creating a brick\n", __LINE__);
      return NULL;
   }
   TRANSFORM* identity = new TRANSFORM(scale_transf(1.0));
   result = api_change_body_trans ( BODY_ptr, identity, FALSE );
   identity->lose();
   if (!result.ok())
   {
      AcisQueryEngine::instance()->ACIS_API_error ( result );
      PRINT_ERROR("In AcisGeometryEngine::make_brick_BODY, Line# %d\n"
                  "     Problems creating a brick\n", __LINE__);
      return NULL;

   }
     //The bounding SPAbox does not get re-computed, so we need
     //to do it.
   AcisQueryEngine::instance()->clear_bounding_box( BODY_ptr );
   AcisQueryEngine::instance()->bounding_box( BODY_ptr );

   return BODY_ptr;
}

//-------------------------------------------------------------------------
// Purpose       : "webcut" the input BODY (bodyPtr) using the
//                 input cutting tool BODY (CuttingToolPtr).
//                 The resulting BODYs are returned as webcutBody1
//                 and webcutBody2.
//
// Special Notes : The input BODY is not affected by the operation
//
// Creator       :
//
// Creation Date :
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::webcut_BODY( BODY* bodyPtr,
                                           BODY* CuttingToolPtr,
                                           BODY*& webcutBody1,
                                           BODY*& webcutBody2) const
{
     // Algorithm for webcut of the ACIS BODY, bodyPtr, using the solid
     // ACIS BODY, CuttingToolPtr, as the Cutting Tool:
     // 1) Save 2 copies of bodyPtr as body1Save and body2Save
     // 2) Save 2 copies of CuttingToolPtr as CuttingToolSave1 and CuttingToolSave2
     // 3) Subtract CuttingToolSave1 from bodySave1 => returns Body2
     // 4) Intersect CuttingToolSave2 with bodySave2 => returns Body3
     // 5) Clean up and return Body2 and Body3
     //

     // MJP Note:
     // By the end of this procedure, whether the operation is successful or not,
     // BODYs pointed to by CuttingToolPtr and bodyPtr are returned intact.

     // MJP Note:
     // The final 2 ACIS BODYs are *not* combined into a single BODY with
     // multiple LUMPs as Acis does not consider a BODY with multiple LUMPs some
     // of whose FACEs overlap in space, a valid BODY. LUMPs of a BODY need
     // to be disjoint.

   int webcut_debug_flag = 18;
   int webcut_debug_flag_on = DEBUG_FLAG(webcut_debug_flag);

     // Debug info
   if (webcut_debug_flag_on)
   {
      PRINT_DEBUG_18("\n...Performing the webcut operation on the Acis BODY\n");
   }

     // Set the output webcut BODY pointers to NULL
   webcutBody1 = NULL;
   webcutBody2 = NULL;

     // Make sure that none of the input BODY pointers are NULL.
     // Note: webcutBody1 and webcutBody2 are "empty" Acis BODYs that were
     //       created by the calling procedure and will probably be set to
     //       NULL
   if (bodyPtr == NULL)
   {
      PRINT_ERROR("The pointer to the input volume to be webcut is NULL.\n");
      return CUBIT_FAILURE;
   }
   if (CuttingToolPtr == NULL)
   {
      PRINT_ERROR("The pointer to the input CuttingTool volume is NULL.\n");
      return CUBIT_FAILURE;
   }

     // Save the BODY being cut, if debug is on
   if (webcut_debug_flag_on)
   {
      AcisQueryEngine::instance()->save_ENTITY_as_sat_file(bodyPtr, "OrigBody.sat", "w");
   }

     // Make a copy of the body to be decomposed as well as the
     // CuttingTool.
   BODY* body1Save = NULL;
   BODY* body2Save = NULL;
   BODY* CuttingToolSave = NULL;
   BODY* CuttingToolSave2 = NULL;

   CubitBoolean copy_attribs = GeometryModifyTool::instance()->get_new_ids();
   CubitBoolean delete_owners = CUBIT_TRUE;
   if (copy_attribs == CUBIT_TRUE) delete_owners = CUBIT_FALSE;

   body1Save = this->copy_BODY(bodyPtr, copy_attribs);
   if ( body1Save == NULL )
   {
        // If the original ACIS BODY was destroyed, then make sure that the
        // associated Body is deactivated and removed.
      if (bodyPtr == NULL)
      {
         PRINT_ERROR("During the webcut operation, the original\n"
                     "         ACIS volume was destroyed.\n"
                     "      ...Deleting the original volume...\n");
      }
      return CUBIT_FAILURE;
   }

   CuttingToolSave = this->copy_BODY(CuttingToolPtr);
   if ( CuttingToolSave == NULL )
   {
      AcisQueryEngine::instance()->delete_ACIS_BODY(body1Save, delete_owners);
      return CUBIT_FAILURE;
   }

     // Subtract CuttingToolSave from body1Save
     // NOTE: in the api_subtract, the resulting body is put into the
     //       body being subtracted from.
   outcome result = api_subtract(CuttingToolSave, body1Save);
   if ( !result.ok() )
   {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      if (CuttingToolSave != NULL)
          AcisQueryEngine::instance()->delete_ACIS_BODY(CuttingToolSave, CUBIT_TRUE);
      if (body1Save != NULL)
          AcisQueryEngine::instance()->delete_ACIS_BODY(body1Save, delete_owners);
      return CUBIT_FAILURE;
   }

     // As the subtract operation was successful, the blank has already been deleted.
     // So set the pointer to NULL so we don't try to delete it again.
   CuttingToolSave = NULL;

     // Before proceeding any further, make sure that result of the
     // subtraction is not a NULL BODY or a BODY that has no LUMPs.
   if ( (body1Save == NULL) || (body1Save->lump() == NULL) )
   {
      if (body1Save != NULL)
          AcisQueryEngine::instance()->delete_ACIS_BODY(body1Save, delete_owners);
      PRINT_INFO("INFO: Cutting Tool overlaps the original volume,\n"
                 "      Or cutting plane does not pass through volume.\n"
                 "         The original volume is unaffected.\n" );
      return CUBIT_FAILURE;
   }

     // Make another copy of the original body.
     // NOTE: I want to preserve the original body in case something
     //       goes haywire before we're done :-)
   body2Save = this->copy_BODY(bodyPtr, copy_attribs);
   if ( body2Save == NULL )
   {
      if (body1Save != NULL)
          AcisQueryEngine::instance()->delete_ACIS_BODY(body1Save, delete_owners);
      if (bodyPtr == NULL)
      {
         PRINT_ERROR("During the webcut operation, the original\n"
                     "         volume was destroyed.\n"
                     "      ...Deleting the original volume...\n");
      }
      return CUBIT_FAILURE;
   }

     // Make another copy of the CuttingTool body
   CuttingToolSave2 = this->copy_BODY(CuttingToolPtr);
   if ( CuttingToolSave2 == NULL )
   {
      if (body1Save != NULL)
          AcisQueryEngine::instance()->delete_ACIS_BODY(body1Save, delete_owners);
      if (body2Save != NULL)
          AcisQueryEngine::instance()->delete_ACIS_BODY(body2Save, delete_owners);
      if (bodyPtr == NULL)
      {
         PRINT_ERROR("During the webcut operation, the original\n"
                     "         volume was destroyed.\n"
                     "      ...Deleting the original volume...\n");
      }

      return CUBIT_FAILURE;
   }

     // Now intersect CuttingToolSave2 with body2Save
     // Note: If the outcome is successful, the result is the second BODY
     //       and the first BODY is deleted.
     // Debug output
   if (webcut_debug_flag_on)
   {
        //Save the result of the subtract operation
      AcisQueryEngine::instance()->save_ENTITY_as_sat_file(CuttingToolPtr, "Subtract.sat", "w");

        //Save the result of the intersect operation
      AcisQueryEngine::instance()->save_ENTITY_as_sat_file(bodyPtr, "Intersect.sat", "w");
   }
   result = api_intersect(CuttingToolSave2, body2Save);
   if ( !result.ok() )
   {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      if (CuttingToolSave2 != NULL)
          AcisQueryEngine::instance()->delete_ACIS_BODY(CuttingToolSave2, CUBIT_TRUE);
      if (body1Save != NULL)
          AcisQueryEngine::instance()->delete_ACIS_BODY(body1Save, delete_owners);
      if (body2Save != NULL)
          AcisQueryEngine::instance()->delete_ACIS_BODY(body2Save, delete_owners);
      return CUBIT_FAILURE;
   }

     // As the intersect operation was successful, the blank has already
     // been deleted.
     // So set the pointer to NULL so we don't try to delete it again.
   CuttingToolSave2 = NULL;

     // Debug output
   if (webcut_debug_flag_on)
   {
        //Save the result of the subtract operation
      AcisQueryEngine::instance()->save_ENTITY_as_sat_file(body1Save, "Subtract.sat", "w");

        //Save the result of the intersect operation
      AcisQueryEngine::instance()->save_ENTITY_as_sat_file(body2Save, "Intersect.sat", "w");
   }

     // Check to make sure that the BODY resulting from the intersection
     // operation is not NULL and it has at lease one valid LUMP. This
     // would occur if there is no intersection between the 2 bodies
     // being operated on.
   if ( (body2Save == NULL) || (body2Save->lump() == NULL) )
   {
      PRINT_WARNING("Cutting Tool does not intersect the original volume.\n"
                    "         The original volume is restored.\n" );

      AcisQueryEngine::instance()->delete_ACIS_BODY(body2Save, delete_owners);
      return CUBIT_FAILURE;
   }
     // Finally, set the BODY pointers, webcutBodyPtr1 and webcutBodyPtr2, to
     // the resulting Acis BODYs (after the imprint operation)
   webcutBody1 = body1Save;
   webcutBody2 = body2Save;
   AcisQueryEngine::instance()->clear_bounding_box( webcutBody1 );
   AcisQueryEngine::instance()->bounding_box( webcutBody1 );
   AcisQueryEngine::instance()->clear_bounding_box( webcutBody2 );
   AcisQueryEngine::instance()->bounding_box( webcutBody2 );

   return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Create a VERTEX given a location.
//
// Special Notes : At this time, the input VERTEX_type argument is not
//                 used.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 03-27-97
//-------------------------------------------------------------------------
VERTEX* AcisModifyEngine::make_VERTEX( CubitVector const& point ) const
{
  CubitStatus status = CUBIT_SUCCESS;
  VERTEX* VERTEX_ptr = NULL;
  API_BEGIN;

      // Create a new ACIS VERTEX -- it is not attached to a BODY
    APOINT* APOINT_ptr = new APOINT ( point.x(), point.y(), point.z() );

    if ( APOINT_ptr != NULL )
    {
      VERTEX_ptr = new VERTEX ( APOINT_ptr );
    }

    else
    {
      PRINT_ERROR("In AcisModifyEngine::make_VERTEX\n"
                  "       ACIS could not make an APOINT object.\n"
                  "       Cannot make VERTEX.\n");
      status = CUBIT_FAILURE;
    }

      // Complain if ACIS wasn't able to create the VERTEX
    if ( VERTEX_ptr == NULL && status == CUBIT_SUCCESS)
    {
      PRINT_ERROR("In AcisModifyEngine::make_VERTEX\n"
                  "       ACIS could not make a VERTEX object.\n");
      status = CUBIT_FAILURE;
    }

  API_END;
  if ( status == CUBIT_SUCCESS )
    return VERTEX_ptr;

  return (VERTEX*)NULL;
}

//-------------------------------------------------------------------------
// Purpose       : Create an EDGE of type EDGE_type.
//
// Special Notes :
//
// For STRAIGHT_CURVE_TYPE:
//    intermediate_point_ptr = NULL
//    sense is not used
//
// For PARABOLA_CURVE_TYPE
//    intermediate_point_ptr is the tip of the parabola
//    sense is not used
//
// For ELLIPSE_CURVE_TYPE
//    intermediate_point_ptr is the center of the ellipse
//    sense is used to determine which part of the ellipse is required
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 03-27-97
//-------------------------------------------------------------------------
EDGE* AcisModifyEngine::make_EDGE( GeometryType EDGE_type,
                                     VERTEX* start_VERTEX,
                                     VERTEX* end_VERTEX,
                                     SPAposition* intermediate_point_ptr,
                                   bool sense ) const
{
     // Make sure the input VERTEX'es are non-NULL -- these are needed by
     // all the EDGE generation routines.
   if ( start_VERTEX == NULL || end_VERTEX == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_EDGE\n"
                  "       An input end VERTEX is invalid (NULL).\n");
      assert (start_VERTEX != NULL && end_VERTEX != NULL);
      return (EDGE *)NULL;
   }

     // Make a straight EDGE
   if ( EDGE_type == STRAIGHT_CURVE_TYPE )
   {
      return make_straight_EDGE(start_VERTEX, end_VERTEX);
   }

     // Make a parabolic EDGE
   else if ( EDGE_type == PARABOLA_CURVE_TYPE )
   {
      return make_parabolic_EDGE( start_VERTEX,
                                  end_VERTEX,
                                  intermediate_point_ptr);
   }

     // Make an elliptical EDGE
   else if ( EDGE_type == ELLIPSE_CURVE_TYPE )
   {
      return make_elliptical_EDGE( start_VERTEX,
                                   end_VERTEX,
                                   intermediate_point_ptr,
                                   sense);
   }

   else if (EDGE_type == ARC_CURVE_TYPE)
   {
     return make_circular_EDGE(start_VERTEX,
                               end_VERTEX,
                               intermediate_point_ptr);
   }

     // OOps!!!!
   else
   {
      PRINT_ERROR("In AcisModifyEngine::make_EDGE\n"
                  "       Invalid curve type.\n");
      assert (0 == 1);
      return (EDGE *)NULL;
   }
}

// This function will create a supper cone using the info
// from the existing cone.  It is actually a face or sheet
// body that will result from this.  The height of the cone
// is determined by a factor of the model limits currently
// in cubit.
FACE* AcisModifyEngine::make_type_FACE( CONE *CONE_ptr ) const
{
     //Given this lets create another face that is extended from it.
   CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();

   BODY *sheet_body;
   CubitBoolean cylinder = CUBIT_FALSE;

   double const st_ang = 0.0;
   double const end_ang = 2*CUBIT_PI;
   double cost = CONE_ptr->cosine_angle();
   double sint = CONE_ptr->sine_angle();
   if ( sint == 0.0 && fabs(cost) == 1.0 )
   {
       //This is a cylinder.
     cylinder = CUBIT_TRUE;
   }
   SPAvector major_a = CONE_ptr->major_axis();
   SPAunit_vector normal_a = CONE_ptr->direction();
   double const radius_ratio = CONE_ptr->radius_ratio();
   double const height = 2*(bounding_box.diagonal()).length();
   SPAposition root = CONE_ptr->root_point();
   CubitVector elipse_radius, center_vec;
   CubitVector temp_norm( normal_a.x(), normal_a.y(), normal_a.z());
   if ( !cylinder )
   {
     if ( sint > 0.0 && cost > 0.0 )
     {
       sint = -sint;
       temp_norm *= -1.0;
     }
     else if ( sint < 0.0 && cost < 0.0 )
     {
       cost = -cost;
       temp_norm *= -1.0;
     }
       //First we need to calculate the new radius.
       //Also calc the old height to find the new center SPAposition.
     double theta_rad = acos(cost);
     elipse_radius.set(major_a.x(), major_a.y(), major_a.z() );
     double old_radius = elipse_radius.length();
     double old_height = old_radius/tan(theta_rad);

     //make center of surface along axis, 1/2 height of
     //bounding box of entire model.
     CubitVector root_vec( root.x(), root.y(), root.z() );
     center_vec = root_vec + (-0.5*height)*temp_norm;

     //calculate new radius
     double new_radius = (old_height+ 0.5*height) * tan (theta_rad );
     elipse_radius.normalize();
     elipse_radius *= new_radius;
   }
   else
   {
     CubitVector root_vec( root.x(), root.y(), root.z() );
     center_vec = .5*height*temp_norm;
     center_vec = root_vec - center_vec;
     elipse_radius.set(major_a.x(), major_a.y(), major_a.z() );
   }
   SPAvector const major_axis(elipse_radius.x(), elipse_radius.y(), elipse_radius.z() );

   SPAunit_vector const normal_axis(temp_norm.x(), temp_norm.y(), temp_norm.z());
   SPAposition const center ( center_vec.x(), center_vec.y(), center_vec.z() );


   FACE *new_FACE_ptr = NULL;
   outcome result = api_make_cnface( center, normal_axis, major_axis,
                             radius_ratio, sint, cost, st_ang,
                             end_ang, height, new_FACE_ptr );
   if (!result.ok() || new_FACE_ptr == NULL ||
       new_FACE_ptr->geometry() == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(cone)\n"
                  "       ACIS api_make_cnface function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (FACE *)NULL;
   }
   FACE *face_list[1];
   face_list[0] = new_FACE_ptr;
   result = api_sheet_from_ff( 1, face_list, sheet_body );
   if (!result.ok() || sheet_body == NULL || sheet_body->lump() == NULL
       || sheet_body->lump()->shell() == NULL ||
       sheet_body->lump()->shell()->first_face() == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(cone)\n"
                  "       ACIS api_sheet_from_ff function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (FACE *)NULL;
   }
   result =  api_body_to_2d( sheet_body );
   if (!result.ok())
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(cone)\n"
                  "       ACIS api_body_to_2d function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (FACE *)NULL;
   }
   int debug_flg = 0;
     // Save the BODY being cut, if debug is on
   if (debug_flg)
   {
      AcisQueryEngine::instance()->save_ENTITY_as_sat_file(sheet_body, "cone_sheet.sat", "w");
   }

   return sheet_body->lump()->shell()->first_face();
}
// This function will create a sphere face using the info
// from the existing sphere.  It is actually a face or sheet
// body that will result from this. The new spher face or sheet
// will be the grand sphere itself, not a section of it.
FACE* AcisModifyEngine::make_type_FACE( SPHERE *SPHERE_ptr ) const
{
   BODY *sheet_body;

   BODY *sphere_body;
   SPAposition const center = SPHERE_ptr->centre();
   double const radius = SPHERE_ptr->radius();

     //make a sphere body then move it to the give place.
   outcome result = api_make_sphere( radius, sphere_body );
   if (!result.ok() || sphere_body == NULL ||
       sphere_body->lump() == NULL ||
       sphere_body->lump()->shell() == NULL ||
       sphere_body->lump()->shell()->first_face() == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(sphere)\n"
                  "       ACIS api_make_sphere function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (FACE *)NULL;
   }
     //now move the body.
   SPAvector translation_vector( center.x(), center.y(), center.z() );
   result = api_apply_transf ( sphere_body,
                               translate_transf (translation_vector) );
   if ( !result.ok() )
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(sphere)\n"
                  "       ACIS api_aply_transf function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error ( result, "transforming body" );
      return (FACE*)NULL;
   }
   TRANSFORM* identity = new TRANSFORM(scale_transf(1.0));
   result = api_change_body_trans ( sphere_body, identity, FALSE );
   identity->lose();
   if (!result.ok())
   {
      AcisQueryEngine::instance()->ACIS_API_error ( result, "transforming body" );
      return (FACE*)NULL;

   }
     //The bounding SPAbox does not get re-computed, so we need
     //to do it.
   AcisQueryEngine::instance()->clear_bounding_box( sphere_body );
   AcisQueryEngine::instance()->bounding_box(sphere_body);


   FACE *face_list[1];
   face_list[0] = sphere_body->lump()->shell()->first_face();
   result = api_sheet_from_ff( 1, face_list, sheet_body );
   if (!result.ok() || sheet_body == NULL || sheet_body->lump() == NULL
       || sheet_body->lump()->shell() == NULL ||
       sheet_body->lump()->shell()->first_face() == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(sphere)\n"
                  "       ACIS api_sheet_from_ff function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (FACE *)NULL;
   }
   result =  api_body_to_2d( sheet_body );
   if (!result.ok())
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(sphere)\n"
                  "       ACIS api_body_to_2d function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (FACE *)NULL;
   }
   AcisQueryEngine::instance()->delete_ACIS_BODY(sphere_body, CUBIT_TRUE);

   return sheet_body->lump()->shell()->first_face();
}
// This function will create a torus face using the info
// from the existing torus.  It is actually a face or sheet
// body that will result from this. The new spher face or sheet
// will be the grand torus itself, not a section of it.
FACE* AcisModifyEngine::make_type_FACE( TORUS *TORUS_ptr ) const
{
  if ( TORUS_ptr == NULL )
  {
    PRINT_ERROR("NULL pointer passed to make face from. (torus).\n");
    return (FACE*)NULL;
  }

   BODY *sheet_body;
   SPAposition const center = TORUS_ptr->centre();
   SPAunit_vector const normal = TORUS_ptr->normal();
   double const uf = 0.0;
   double const ut = 2*CUBIT_PI;
   double const vf = 0.0;
   double const vt = 2*CUBIT_PI;
   double const major_radius = TORUS_ptr->major_radius();
   double const minor_radius = TORUS_ptr->minor_radius();
   CubitVector temp_center(center.x(), center.y(), center.z());
   CubitVector temp_normal(normal.x(), normal.y(), normal.z());
   temp_normal.normalize();
   CubitVector temp_point = temp_center + major_radius*temp_normal;
   SPAposition const point(temp_point.x(), temp_point.y(), temp_point.z());

     //make a torus body then move it to the give place.
   FACE *new_FACE_ptr = NULL;
   outcome result = api_make_trface( center, normal, major_radius,
                                     minor_radius, point,
                                     uf, ut, vf, vt, new_FACE_ptr );
   if (!result.ok() || new_FACE_ptr == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(torus)\n"
                  "       ACIS api_make_trface function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (FACE *)NULL;
   }
   FACE *face_list[1];
   face_list[0] = new_FACE_ptr;
   result = api_sheet_from_ff( 1, face_list, sheet_body );
   if (!result.ok() || sheet_body == NULL || sheet_body->lump() == NULL
       || sheet_body->lump()->shell() == NULL ||
       sheet_body->lump()->shell()->first_face() == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(torus)\n"
                  "       ACIS api_sheet_from_ff function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (FACE *)NULL;
   }
   result =  api_body_to_2d( sheet_body );
   if (!result.ok())
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(torus)\n"
                  "       ACIS api_body_to_2d function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (FACE *)NULL;
   }
   int debug_flg = 0;
     // Save the BODY being cut, if debug is on
   if (debug_flg)
   {
      AcisQueryEngine::instance()->save_ENTITY_as_sat_file(sheet_body, "torus_sheet.sat", "w");
   }

   return sheet_body->lump()->shell()->first_face();
}
// This function will create a spline face using the info
// from the existing spline.  It is actually a face or sheet
// body that will result from this. The new spher face or sheet
// will be the grand spline itself, not a section of it.
FACE* AcisModifyEngine::make_type_FACE( SPLINE *SPLINE_ptr ) const
{
   BODY *sheet_body;
   PRINT_WARNING("Extents of spline will be only as far as\n"
                 "the spline definition, this really isn't an infinite surface.\n");

     //use the geometry to make the spline body.
   surface const *this_spline = &(SPLINE_ptr->equation());

     //make a spline body then move it to the give place.
   outcome result = api_make_spline( *((spline*)this_spline), sheet_body );
   if (!result.ok() || sheet_body == NULL ||
       sheet_body->lump() == NULL ||
       sheet_body->lump()->shell() == NULL ||
       sheet_body->lump()->shell()->first_face() == NULL ||
       sheet_body->lump()->shell()->first_face()->geometry() == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(spline)\n"
                  "       ACIS api_make_spline function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (FACE *)NULL;
   }
   result =  api_body_to_2d( sheet_body );
   if (!result.ok())
   {
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(spline)\n"
                  "       ACIS api_body_to_2d function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (FACE *)NULL;
   }

   return sheet_body->lump()->shell()->first_face();
}
// This function will create a plane face using the info
// from the existing plane.  It is actually a face or sheet
// body that will result from this. The new spher face or sheet
// will be the grand plane itself, not a section of it.
FACE* AcisModifyEngine::make_type_FACE( PLANE *PLANE_ptr ) const
{
     //first lets find 3 points.
   SPAposition point_1, point_2, point_3;
   SPApar_pos const pos_1( 0, 1 );
   SPApar_pos const pos_2( 1, 0 );
   SPApar_pos const pos_3( 0, 0 );
   PLANE_ptr->equation().eval( pos_1, point_1 );
   PLANE_ptr->equation().eval( pos_2, point_2 );
   PLANE_ptr->equation().eval( pos_3, point_3 );

   CubitBox cubit_super_box = GeometryQueryTool::instance()->model_bounding_box();
   SPAbox super_box = AcisQueryEngine::bounding_box(cubit_super_box);

   CubitVector vec_1( point_1.x(), point_1.y(), point_1.z() );
   CubitVector vec_2( point_2.x(), point_2.y(), point_2.z() );
   CubitVector vec_3( point_3.x(), point_3.y(), point_3.z() );

   BODY *sheet_body = this->create_infinite_plane_cutting_tool(
       vec_1, vec_2, vec_3, super_box, CUBIT_TRUE );
   if (sheet_body == NULL ||
       sheet_body->lump() == NULL ||
       sheet_body->lump()->shell() == NULL ||
       sheet_body->lump()->shell()->first_face() == NULL ||
       sheet_body->lump()->shell()->first_face()->geometry() == NULL )
   {
        //webcut_failed(refbody_cleanup_list, BODY_cleanup_list);
      PRINT_ERROR("In AcisModifyEngine::make_type_FACE(plane)"
                  "       Cannot create infinite plane cutting tool using \n"
                  "       the give surface. \n");
      return (FACE*)NULL;
   }
   return sheet_body->lump()->shell()->first_face();
}
//-------------------------------------------------------------------------
// Purpose       : Create a face from the face pointer.  This will create
//                 a sheet body too.  It copies all the underlying data
//                 structure.
//
// Creator       : David White
//
// Creation Date : 10/6/97
//-------------------------------------------------------------------------
FACE* AcisModifyEngine::make_FACE(FACE* FACE_ptr,
                                    CubitBoolean extended_from) const
{
  if (FACE_ptr == NULL)
    return NULL;

  if (extended_from == CUBIT_TRUE)
  {
      // We need to get the type of surface.
      // For some reason I can't get it to cast on the fly.
    SURFACE *SURFACE_ptr = FACE_ptr->geometry();
    int type = SURFACE_ptr->identity();
    if (type == CONE_TYPE)
      return make_type_FACE((CONE*)SURFACE_ptr);
    else if (type == SPHERE_TYPE)
      return make_type_FACE((SPHERE*)SURFACE_ptr);
    else if (type == TORUS_TYPE)
      return make_type_FACE((TORUS*)SURFACE_ptr);
    else if (type == SPLINE_TYPE)
      return make_type_FACE((SPLINE*)SURFACE_ptr);
    else if (type == PLANE_TYPE)
      return make_type_FACE((PLANE*)SURFACE_ptr);
    else
      return NULL;
  }

  FACE *face_list[1];
  BODY *sheet_body;
  face_list[0] = FACE_ptr;

    // Make a FACE from copy_FACE.
  outcome rc = api_sheet_from_ff(1, face_list, sheet_body);

    // Now check to see if we created the sheet successfully.
  if (!rc.ok())
  {
    PRINT_ERROR("In AcisModifyEngine::make_FACE\n"
                "       ACIS api_sheet_from_ff function failed.\n\n");
    AcisQueryEngine::instance()->ACIS_API_error (rc);
    return (FACE *)NULL;
  }

    // Convert the FACE to a true sheet BODY.
  outcome rc_2 =  api_body_to_2d(sheet_body);
  if (!rc_2.ok())
  {
    PRINT_ERROR("In AcisModifyEngine::make_FACE\n"
                "       ACIS api_body_to_2d function failed.\n\n");
    AcisQueryEngine::instance()->ACIS_API_error (rc_2);
    return (FACE *)NULL;
  }

    // Make sure we were successful
  if ( sheet_body->lump() == NULL ||
       sheet_body->lump()->shell() == NULL ||
       sheet_body->lump()->shell()->first_face() == NULL )
  {
    PRINT_ERROR("In AcisModifyEngine::make_FACE\n"
                "       ACIS api_sheet_from_ff function failed.\n\n");
    return (FACE*)NULL;
  }

  return sheet_body->lump()->shell()->first_face();
}

//-------------------------------------------------------------------------
// Purpose       : Creates an ACIS FACE, given the type of surface geometry
//                 and an ordered list of bounding EDGEs.
//
// Special Notes : The list of EDGEs are assumed to form the outer LOOP
//                 of the FACE to be created.  The input EDGEs are not
//                 guaranteed to exist after this function returns.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------
FACE* AcisModifyEngine::make_FACE( GeometryType FACE_type,
                                   DLIList<EDGE*>& input_EDGE_list,
                                   Surface *old_surface_ptr) const
{
    // Make sure a supported type of surface is being requested.
  if ( FACE_type != PLANE_SURFACE_TYPE  &&
       FACE_type != BEST_FIT_SURFACE_TYPE)
  {
    PRINT_WARNING("In AcisGeometryEngine::make_FACE\n"
                  "         At this time, cannot make a FACE that isn't"
                  " planar or best fit.\n");
    return (FACE *)NULL;
  }

    // Set the FACE pointer, if requested.
  FACE *fit_FACE = NULL;
  if ( old_surface_ptr != NULL )
  {
    SurfaceACIS *surf_ACIS = CAST_TO(old_surface_ptr, SurfaceACIS );
    fit_FACE = surf_ACIS->get_FACE_ptr();
  }

    // Copy the original EDGEs and put them into an array.
    // Also initialize variables that will be used/returned by ACIS.
  int edge_count = input_EDGE_list.size();
  FACE* FACE_ptr = NULL;
  BODY* wire_BODY = NULL;
  EDGE** EDGEs = new EDGE* [edge_count];
  ENTITY* temp_ENTITY_ptr = NULL;
  input_EDGE_list.reset();
  int i;

  for (i = 0 ; i < edge_count ; i++)
  {
    copy_single_entity( (ENTITY *)input_EDGE_list.get_and_step(), temp_ENTITY_ptr);
    EDGEs[i] = (EDGE*)temp_ENTITY_ptr;
      // Deactivate the RefEntities attached to the EDGE and its VERTEXes.
    mark_owners_deactivated_flag(EDGEs[i], CUBIT_TRUE, CUBIT_TRUE);
  }

  // Make a wire BODY from the copied EDGEs.
  // Coincident VERTEXes will be deleted by ACIS.
  //  outcome result = api_make_ewire(edge_count, EDGEs, wire_BODY) ;
  int num_wires = 0;
  BODY **wireBODYs;
  ENTITY_LIST wire_BODYs;
  outcome result = api_make_ewires(edge_count, EDGEs, num_wires, wireBODYs);
  for(i=0; i<num_wires; i++ )
    wire_BODYs.add( wireBODYs[i] );

  if ( !result.ok() || wire_BODYs.count() == 0 )
  {
    PRINT_ERROR("In AcisGeometryEngine::make_FACE\n"
                "       ACIS api_make_ewire function failed.\n\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    if (wire_BODYs.count())
      for(i=0; i<wire_BODYs.count(); i++ )
        api_delent(wire_BODYs[i]);
    else
    {
      while (edge_count--)
        api_delent(EDGEs[edge_count]);
    }

    delete [] EDGEs;
    return (FACE *)NULL;
  }

    // Check to make sure the wire is closed
  for(i=0; i<wire_BODYs.count(); i++ )
  {
    result = api_closed_wire( (BODY*)wire_BODYs[i] );
    if( !result.ok() )
    {
      PRINT_ERROR( "ACIS reports that a closed loop was not found\n" );
        // Delete the wire BODY.  This will also delete
        // each copied EDGE
      for(i=0; i<wire_BODYs.count(); i++ )
        api_delent(wire_BODYs[i]);
        // Delete the array
      delete [] EDGEs;
      return (FACE *)NULL;
    }
  }


  //check planarity & co-planarity
  if( wire_BODYs.count() > 1 )
  {
    CubitVector master_centroid;
    CubitVector master_normal;
    for(i=0; i<wire_BODYs.count(); i++ )
    {
      SPAposition centroid;
      SPAunit_vector normal;

      if( !is_planar_wire( (WIRE*)wire_BODYs[i], centroid, normal ) )
      {
        PRINT_ERROR( "Multiple loops must be planar.\n" );
        for(i=0; i<wire_BODYs.count(); i++ )
          api_delent(wire_BODYs[i]);
        delete [] EDGEs;
        return (FACE *)NULL;
      }
      //make the plane
      if( i==0 )
      {
        master_centroid.set( centroid.x(), centroid.y(), centroid.z() );
        master_normal.set( normal.x(), normal.y(), normal.z() );
      }
      else
      {
        CubitVector tmp_centroid( centroid.x(), centroid.y(), centroid.z() );
        CubitVector tmp_normal( normal.x(), normal.y(), normal.z() );
        CubitVector tmp_vec( tmp_centroid - master_centroid );
        tmp_vec.normalize();
        double dot = master_normal % tmp_vec;
        tmp_normal = tmp_normal * master_normal;

        if( tmp_normal.length() > GEOMETRY_RESABS ||
            fabs(dot) > GEOMETRY_RESABS )
        {
          PRINT_ERROR( "Multiple loops must be co-planar.\n" );
          for(i=0; i<wire_BODYs.count(); i++ )
            api_delent(wire_BODYs[i]);
          delete [] EDGEs;
          return (FACE *)NULL;
        }
      }
    }
  }

  //check intersecting
  if( wire_BODYs.count() > 1 )
  {
    for(i=0; i<wire_BODYs.count(); i++ )
    {
      ENTITY *ent1 = wire_BODYs[i];

      int j;
      for(j=0; j<wire_BODYs.count(); j++ )
      {
        if( i==j )
          continue;

        ENTITY *ent2 = wire_BODYs[j];
        logical touching = FALSE;

        api_entity_entity_touch( ent1, ent2, touching );

        if( touching )
        {
          PRINT_ERROR( "Multiple loops cannot intersect.\n" );
          for(i=0; i<wire_BODYs.count(); i++ )
            api_delent(wire_BODYs[i]);
          delete [] EDGEs;
          return (FACE *)NULL;
        }
      }
    }
  }

  //sort wires so that the outermost is first
  wire_BODYs.init();
  ENTITY_LIST tmp_wire_BODYs( wire_BODYs ) ;
  ENTITY *outermost_wire= NULL;
  double longest_diagonal = 0.0;
  if( wire_BODYs.count() > 1 )
  {
    for(i=0; i<wire_BODYs.count(); i++ )
    {
      SPAbox tmp_box;
      CubitBox tmp_cbox;
      ENTITY *tmp_wire = wire_BODYs.next();
      if(i==0)
      {
        tmp_box = AcisQueryEngine::instance()->bounding_box( (BODY*)tmp_wire);
        outermost_wire = tmp_wire;
        tmp_cbox = AcisQueryEngine::bounding_box( tmp_box );
        longest_diagonal = tmp_cbox.diagonal().length();
      }
      else
      {
        tmp_box = AcisQueryEngine::instance()->bounding_box( (BODY*)tmp_wire);
        tmp_cbox = AcisQueryEngine::bounding_box( tmp_box );

        if( tmp_cbox.diagonal().length() > longest_diagonal )
        {
          longest_diagonal = tmp_cbox.diagonal().length();
          outermost_wire = tmp_wire;
        }
      }
    }

    //now that we found the outermost, put it first in the list
    wire_BODYs.clear();
    tmp_wire_BODYs.remove( outermost_wire );
    wire_BODYs.add( outermost_wire );

    tmp_wire_BODYs.init();
    for(i=0; i<tmp_wire_BODYs.count(); i++ )
      wire_BODYs.add( tmp_wire_BODYs.next() );
  }

    // Use the WIRE to make a FACE.
    // Note that the call to api_cover_wires creates not
    // only the FACE but the entire topological data
    // structure above FACE (till BODY).
  ENTITY_LIST acis_FACE_list;
  if ( fit_FACE != NULL )
  {
    wire_BODY = static_cast<BODY*>(wire_BODYs[0]);
    outcome result_2 = api_cover_wires ( wire_BODY,
                                         fit_FACE->geometry()->equation(),
                                         acis_FACE_list );

      // Now check to see if we covered the wire successfully
    if ( !result_2.ok() || acis_FACE_list.count() ==  0 ||
         wire_BODY == NULL )
    {
      PRINT_ERROR("Error in constructing surface on top of another"
                  "surface, the edges must fit the surface...\n");
      AcisQueryEngine::instance()->ACIS_API_error (result_2);
      acis_FACE_list.clear();
      if (wire_BODY)
        api_delent(wire_BODY);
      else
          // TO DO: Delete each EDGE in EDGEs
        void (0);
      delete [] EDGEs;
      return (FACE*)NULL;
    }
  }
  else
  {
    outcome surface_creation;
    bool valid_loops = true;
    outcome loop_check;

    //can create non-planar surfs with one wire
    if( wire_BODYs.count() == 1 )
    {
      wire_BODY = static_cast<BODY*>(wire_BODYs[0]);
      surface_creation = api_cover_wires ( wire_BODY,
                                           *(surface*)NULL_REF,
                                           acis_FACE_list );
    }
    else
    {
      wire_BODY = NULL;
      surface_creation = api_cover_wire_loops( wire_BODYs,
                                               wire_BODY);
      //check all nested in largest
      //check doubly nested
      api_get_faces( wire_BODY, acis_FACE_list);
      loop_check = api_check_face_loops( static_cast<FACE*>(acis_FACE_list[0]) );
      if( !loop_check.ok() )
        valid_loops = false;
    }


    if ( !surface_creation.ok() || wire_BODY == NULL ||
         !valid_loops || acis_FACE_list.count() != 1 )
    {
      PRINT_ERROR("In AcisGeometryEngine::make_FACE\n"
                  "       problem covering wires.\n\n");
      if( !surface_creation.ok() )
        AcisQueryEngine::instance()->ACIS_API_error (surface_creation);
      if( !valid_loops )
        AcisQueryEngine::instance()->ACIS_API_error( loop_check );

      acis_FACE_list.clear();
      if (wire_BODY)
        api_delent(wire_BODY);
      else
      {
          // TO DO: delete each EDGE in EDGEs
        void(0);

        wire_BODYs.init();
        ENTITY *ent_to_del;
        while( (ent_to_del = wire_BODYs.next())!=NULL )
          api_delent( ent_to_del );
      }
      delete [] EDGEs;
      return (FACE *)NULL;
    }
  }
  if (acis_FACE_list.count() != 1)
  {
    PRINT_ERROR("ACIS api_cover_wires function created "
                "more than one FACE.\n"
                "       At this time, we cannot deal with this.\n");
    acis_FACE_list.clear();
    if (wire_BODY)
      api_delent(wire_BODY);
    else
        // TO DO: delete each EDGE in EDGEs
      void(0);
    delete [] EDGEs;
    return (FACE *)NULL;
  }

    //make this a 2d body.
  outcome result_2 = api_body_to_2d( wire_BODY );
  if ( !result_2.ok() )
  {
    PRINT_ERROR("In AcisGeometryEngine::make_FACE\n"
                "       ACIS api_body_to_2d function failed.\n\n");
    AcisQueryEngine::instance()->ACIS_API_error (result_2);
    acis_FACE_list.clear();
    if (wire_BODY)
      api_delent(wire_BODY);
    else
        // TO DO: delete each EDGE in EDGEs
      void(0);
    delete [] EDGEs;
    return (FACE *)NULL;
  }

  FACE_ptr = (FACE*) acis_FACE_list[0];
  if ( FACE_ptr == NULL || FACE_ptr->geometry() == NULL )
  {
    PRINT_ERROR("In AcisGeometryEngine::make_FACE\n"
                "FACE or FACE's geometry was NULL.\n");
    acis_FACE_list.clear();
    if (wire_BODY)
      api_delent(wire_BODY);
    else
        // TO DO: delete each EDGE in EDGEs
      void(0);
    delete [] EDGEs;
    return (FACE*)NULL;
  }
  acis_FACE_list.clear();
  if ( !result.ok() )
  {
    if (wire_BODY)
      api_delent(wire_BODY);
    else
        // TO DO: delete each EDGE in EDGEs
      void(0);
    delete [] EDGEs;
    return (FACE*)NULL;
  }

    // This will reactivate all the RefEntities that need to be kept.
  mark_owners_deactivated_flag(FACE_ptr, CUBIT_FALSE, CUBIT_TRUE);

    // Where coincident vertices were merged by ACIS,
    // make sure the lower ID is used.
  input_EDGE_list.reset();

  for (i = 0; i < input_EDGE_list.size(); i++)
  {
      // get the copy and original EDGE
    EDGE* cur_EDGE_new = EDGEs[i];
    EDGE* cur_EDGE_old = input_EDGE_list.get_and_step();

      // If the original start RefVertex is deactivated,
      // see if we need to swap it with the replacement RefVertex
    VERTEX* cur_VERT_new = cur_EDGE_new->start();
    VERTEX* cur_VERT_old = cur_EDGE_old->start();

    for (int j = 2; j--; )
    {
      RefVertex* ref_vert_old =
        CAST_TO(ATTRIB_CUBIT_OWNER::get_topology_entity(cur_VERT_old),
                RefVertex);

      AcisBridge* bridge = ATTRIB_CUBIT_OWNER::cubit_owner(cur_VERT_old);

      if (ref_vert_old && bridge_deactivated(bridge))
      {
        RefVertex* ref_vert_new =
          CAST_TO(ATTRIB_CUBIT_OWNER::get_topology_entity(cur_VERT_new),
                  RefVertex);
        if (ref_vert_new && ref_vert_old->id() < ref_vert_new->id())
        {
          int id_old = ref_vert_old->id();
          int id_new = ref_vert_new->id();
          ref_vert_old->set_id(0);
          ref_vert_new->set_id(id_old);
          ref_vert_old->set_id(id_new);
        }
      }

/*
      PointACIS* ref_vert_old =
        CAST_TO(ATTRIB_CUBIT_OWNER::cubit_owner(cur_VERT_old),PointACIS);

      if (ref_vert_old && bridge_deactivated(ref_vert_old) && ref_vert_old->owner())
      {
        //RefVertex* ref_vert_new =
        //  CAST_TO(ATTRIB_CUBIT_OWNER::get_topology_entity(cur_VERT_new),
        //          RefVertex);
        PointACIS* ref_vert_new =
          CAST_TO(ATTRIB_CUBIT_OWNER::cubit_owner(cur_VERT_new),PointACIS);
        if (ref_vert_new)
        {
          ref_vert_old->owner()->swap_bridge(ref_vert_old,ref_vert_new,true);
        }
      }
*/
      cur_VERT_new = cur_EDGE_new->end();
      cur_VERT_old = cur_EDGE_old->end();
    }
  }

    // Get rid of the stuff we aren't using.  This will delete the RefVertices
    // that correspond to the VERTEXes deleted by ACIS in api_make_ewire.
    // We also delete the original ACIS EDGEs since they are no longer
    // being used.
  for (i = input_EDGE_list.size(); i--; )
  {
    AcisQueryEngine::instance()->unhook_ENTITY_from_VGI(input_EDGE_list.get());
    api_delent(input_EDGE_list.get_and_step());
  }
  //GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  cleanout_deactivated_geometry();
  delete [] EDGEs;

  return FACE_ptr ;
}

//-------------------------------------------------------------------------
// Purpose       : Create a straight EDGE given 2 VERTEX'es.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 03-27-97
//-------------------------------------------------------------------------
EDGE* AcisModifyEngine::make_straight_EDGE( VERTEX* start_VERTEX,
                                              VERTEX* end_VERTEX ) const
{
     // Make a straight EDGE
   EDGE* EDGE_ptr = NULL;
   outcome result = api_mk_ed_line ( start_VERTEX->geometry()->coords(),
                             end_VERTEX->geometry()->coords(),
                             EDGE_ptr );
   if (!result.ok() || EDGE_ptr == NULL)
   {
      PRINT_ERROR("In AcisModifyEngine::make_straight_EDGE\n"
                  "       ACIS api_mk_ed_line function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (EDGE *)NULL;
   }

   else
   {
        // Switch the end VERTEXes of the new EDGE to the input VERTEXes. Keep
        // in mind that the back links (from the VERTEX back to the EDGE)
        // need to be modified as well.
      EDGE_ptr->set_start(start_VERTEX);
      start_VERTEX->add_edge(EDGE_ptr);

      EDGE_ptr->set_end(end_VERTEX);
      end_VERTEX->add_edge(EDGE_ptr);

        // Return the newly created EDGE
      return EDGE_ptr;
   }
}

//-------------------------------------------------------------------------
// Purpose       : Create a parabolic EDGE given 2 (end) VERTEX'es and the
//                 tip of the parabola.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 03-27-97
//-------------------------------------------------------------------------
EDGE* AcisModifyEngine::make_parabolic_EDGE(
    VERTEX* from,
    VERTEX* to,
    SPAposition* peak_ptr ) const
{
   if ( peak_ptr == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_parabolic_EDGE\n"
                  "       Input parabola tip is invalid (NULL pointer).\n");
      assert ( peak_ptr != NULL );
      return (EDGE *)NULL;
   }

     // Convert positions to CubitVectors (legacy code!!!)
   SPAposition from_position = from->geometry()->coords();
   SPAposition to_position = to->geometry()->coords();

   CubitVector from_vec( from_position.x(),
                         from_position.y(),
                         from_position.z() );
   CubitVector to_vec( to_position.x(),
                       to_position.y(),
                       to_position.z() );
   CubitVector peak_vec( peak_ptr->x(),
                         peak_ptr->y(),
                         peak_ptr->z() );

   CubitVector width_vec, midpoint_vec, height_vec;
   width_vec = to_vec - from_vec;
   midpoint_vec = (from_vec + to_vec)/2.0;
   height_vec = peak_vec - midpoint_vec;

     // rho = 0.5 is a (spline) parabola -
     // see ACIS Geometric Modeler API Reference 3-138 ver. 1.7
   double rho = 0.5;
   EDGE* EDGE_ptr = NULL;
     // Find the focus of this parabola.
     // Since for a parabola with its peak at the origin, y = (1/(4*a))*x^2,
     // and since we have restricted this parabola to be symmetric (per the FastQ
     // method, see the FastQ file getwt.f), we can use the following relationship
     // to determine "a", the distance the focus lies from the peak on the line
     // formed by the peak and the midpoint of the start and end points
   //double a = width_vec.length_squared()/(16. * height_vec.length());
   CubitVector height_vec_unit = height_vec;
   height_vec_unit.normalize();

     // See Geometric Modeling by M. E. Mortenson, pg. 82, Fig 2.28 c, we need p2
     // as the center_vec from which to get the tangents
     // Our approach assumes that midpoint_vec, peak, and center_vec will all lie
     // on the same line, a valid assumption with the restrictions from FastQ
   CubitVector vec_E_C = peak_vec - midpoint_vec;
   CubitVector center_vec = ((2./rho)*vec_E_C + to_vec + from_vec)/2.;

   SPAvector start_tan = SPAvector ( center_vec.x() - from_vec.x(),
                               center_vec.y() - from_vec.y(),
                               center_vec.z() - from_vec.z() );

   SPAvector   end_tan = SPAvector ( center_vec.x() - to_vec.x(),
                               center_vec.y() - to_vec.y(),
                               center_vec.z() - to_vec.z() );

   SPAunit_vector start_tan_unit = normalise( start_tan );
   SPAunit_vector   end_tan_unit = normalise( end_tan );

   outcome result = api_mk_ed_conic ( from_position,
                                      start_tan_unit,
                                      to_position,
                                      end_tan_unit,
                                      rho,
                                      EDGE_ptr );
   if (!result.ok() || EDGE_ptr == NULL)
   {
      PRINT_ERROR("ACIS api_mk_ed_conic function failed\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return (EDGE *)NULL;
   }

   else
   {
        // Switch the end VERTEXes of the EDGE to the input VERTEXes. Keep
        // in mind that the back links (from the VERTEX back to the EDGE)
        // need to be modified as well.
      EDGE_ptr->set_start ( from );
      from->add_edge( EDGE_ptr );

      EDGE_ptr->set_end ( to );
      to->add_edge( EDGE_ptr );

      return EDGE_ptr;
   }
}

//-------------------------------------------------------------------------
// Purpose       : Create an elliptic EDGE given 2 (end) VERTEX'es, the
//                 center of the ellipse and a sense value.
//
// Special Notes : The sense value is used to determine which portion of
//                 the ellipse to keep (there are two ways to go from
//                 the first point to second one, following the curve of
//                 the ellipse). If sense is TRUE, then we keep the
//                 counterclockwise portion; otherwise, we keep the
//                 clockwise portion.
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------

EDGE* AcisModifyEngine::make_elliptical_EDGE( VERTEX* from,
                                                VERTEX* to,
                                                SPAposition* ctr_ptr,
                                                bool sense ) const
{
   if ( ctr_ptr == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_elliptical_EDGE\n"
                  "       Input ellipse center is invalid (NULL pointer).\n");
      return (EDGE *)NULL;
   }

     // Convert positions to CubitVectors (legacy code!!!)
   SPAposition from_position = from->geometry()->coords();
   SPAposition to_position = to->geometry()->coords();

   CubitVector from_vec( from_position.x(),
                         from_position.y(),
                         from_position.z() );
   CubitVector to_vec( to_position.x(),
                       to_position.y(),
                       to_position.z() );
   CubitVector ctr_vec( ctr_ptr->x(),
                        ctr_ptr->y(),
                        ctr_ptr->z() );

     // delta_vee_1 gives major axis of ellipse (circle)
   SPAvector delta_vee_1 = from_position - *ctr_ptr;
   SPAvector delta_vee_2 = to_position   - *ctr_ptr;
   if ( delta_vee_1.len() < CUBIT_RESABS*.000001 )
   {
     PRINT_ERROR("Zero length distance to center point.\n");
     return (EDGE*) NULL;
   }
   if ( delta_vee_2.len() < CUBIT_RESABS*.000001 )
   {
     PRINT_ERROR("Zero length distance to center point.\n");
     return (EDGE*) NULL;
   }
   SPAunit_vector normal = normalise( delta_vee_1 * delta_vee_2 );
     // Since this is a normalised SPAvector, its length should be 1.0.
     // If it is not, then the delta_vee_1 and delta_vee_2
     // vectors are collinear and we have to manually set the normal
     // SPAvector. Assume that this is being done in FASTQ and set the normal
     // to be the Z axis.
   if (normal.len() == 0.0)
   {
        // The start and end vertex are in the same SPAposition, or they are
        // 180 degrees opposed.
      normal.set_x(0.0);
      normal.set_y(0.0);
      normal.set_z(1.0);
   }
   double radius_ratio = delta_vee_2.len() / delta_vee_1.len();

   EDGE* EDGE_ptr = NULL;

     // if radius_ratio is not 1.0, then this code probably won't produce
     // the intended geometry since delta_vee_1 does not necessarily
     // represent the major axis of an ellipse unless it is circular
   outcome result;
   if ( radius_ratio <= 0.999 || radius_ratio > 1.001 )
   {
      EDGE_ptr = make_spiral_EDGE(from, to, ctr_ptr, sense);
   }
   else
   {
      PRINT_DEBUG_5( "radius ratio = %f\n", radius_ratio);

        // Calculate the included angle between from and to on the
        // circular ellipse
      CubitVector cubit_normal (normal.x(), normal.y(), normal.z());
      CubitVector cubit_from = from_vec - ctr_vec;
      CubitVector cubit_to = to_vec - ctr_vec;
      double incl_angle = cubit_normal.vector_angle ( cubit_from,
                                                      cubit_to );
      if (incl_angle == 0.0)
          incl_angle = 2.0 * CUBIT_PI;

      double start_angle = 0.0;
      double end_angle = incl_angle;

      if (radius_ratio > 1.0)
      {
         delta_vee_1 = delta_vee_2;
         radius_ratio = 1.0 / radius_ratio;
         start_angle = -incl_angle;
         end_angle = 0.0;
      }


      result = api_mk_ed_ellipse ( *ctr_ptr,
                                   sense ? normal : -normal,
                                   delta_vee_1,
                                   radius_ratio,
                                   start_angle,
                                   end_angle,
                                   EDGE_ptr );
   }
   if (!result.ok() || EDGE_ptr == NULL)
   {
      PRINT_ERROR("In AGE:make_elliptical_EDGE.\n"
                  "ACIS api_mk_ed_ellipse function failed\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);

      if ( from == NULL || to == NULL )
      {
         PRINT_ERROR("In AGE:make_elliptical_EDGE.\n"
                     "       Input VERTEX is NULL\n");
      }
      return (EDGE *)NULL;
   }

   else
   {
        // Switch the end VERTEXes of the EDGE to the input VERTEXes. Keep
        // in mind that the back links (from the VERTEX back to the EDGE)
        // need to be modified as well.
      EDGE_ptr->set_start ( from );
      from->add_edge( EDGE_ptr );

      EDGE_ptr->set_end ( to );
      to->add_edge( EDGE_ptr );

      return EDGE_ptr ;
   }
}

//-------------------------------------------------------------------------
// Purpose       : Create a circle which passes through the 3 given points
//
// Special Notes :
//
// Creator       : Darryl Melander
//
// Creation Date : 12/10/99
//-------------------------------------------------------------------------
EDGE* AcisModifyEngine::make_circular_EDGE( VERTEX* from,
                                              VERTEX* to,
                                              SPAposition* third_pt ) const
{
  assert (from != NULL);
  assert (to != NULL);
  assert (third_pt != NULL);

    // get the 3 points
  SPAposition start_pt = from->geometry()->coords();
  SPAposition end_pt = to->geometry()->coords();

    // create the arc
  EDGE *e_arc;
  outcome result = api_curve_arc_3pt(start_pt, *third_pt, end_pt,
                                     FALSE, e_arc);
  if (!result.ok())
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    return NULL;
  }

    // Switch the end VERTEXes of the EDGE to the input VERTEXes. Keep
    // in mind that the back links (from the VERTEX back to the EDGE)
    // need to be modified as well.
  e_arc->set_start(from);
  from->add_edge(e_arc);
  e_arc->set_end (to);
  to->add_edge(e_arc);

  return e_arc;
}


//-------------------------------------------------------------------------
// Purpose       : Create a logarithmic spiral given the start and end
//                 points and a point through which the axis of the spiral
//                 passes.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------

EDGE* AcisModifyEngine::make_spiral_EDGE( VERTEX* from,
                                            VERTEX* to,
                                            SPAposition* ctr_ptr,
                                            bool sense) const
{
   if ( ctr_ptr == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_spiral_EDGE\n"
                  "       Input spiral center is invalid (NULL pointer).\n");
      return (EDGE *)NULL;
   }

     // Convert positions to CubitVectors (legacy code!!!)
   SPAposition from_position = from->geometry()->coords();
   SPAposition to_position = to->geometry()->coords();

     // NOTE: This function assumes that vertices from, to, and ctr lie in the
     // X-Y plane.
   if (from_position.z() != to_position.z() ||
       from_position.z() != ctr_ptr->z()) {
      PRINT_ERROR("In AcisModifyEngine::make_spiral_EDGE\n"
                  "       Points must be in the constant-z plane.\n");
      return (EDGE *)NULL;
   }

   CubitVector from_vec( from_position.x(),
                         from_position.y(),
                         from_position.z() );
   CubitVector to_vec( to_position.x(),
                       to_position.y(),
                       to_position.z() );
   //CubitVector ctr_vec( ctr_ptr->x(),
   //                     ctr_ptr->y(),
   //                     ctr_ptr->z() );

   SPAvector delta_vee_1 = from_position - *ctr_ptr;
   SPAvector delta_vee_2 = to_position - *ctr_ptr;

   CubitVector side1( delta_vee_1.x(),
                      delta_vee_1.y(),
                      delta_vee_1.z()) ;

   CubitVector side2( delta_vee_2.x(),
                      delta_vee_2.y(),
                      delta_vee_2.z()) ;

   CubitVector center( ctr_ptr->x(),
                       ctr_ptr->y(),
                       ctr_ptr->z()) ;

     // Calcualate parameters for the spiral. See FASTQ source, arcpar.f
   double radius1 = side1.length();
   double radius2 = side2.length();

   double theta1 = atan2(side1.y(), side1.x());
   double theta2 = atan2(side2.y(), side2.x());
   if (sense && theta2 <= theta1)
   {
      theta2 += 2.0 * CUBIT_PI;
   }
   if (!sense && theta1 <= theta2)
   {
      theta1 += 2.0 * CUBIT_PI;
   }

   double total_angle = theta2 - theta1;

   double xk = log(radius2 / radius1) / (theta2 - theta1);
   double xa = radius2 / exp(xk * theta2);

     // Calculate points at 2.0 degree intervals
   int number_points = int(total_angle * 180.0 / CUBIT_PI / 2.0 + 0.5);
   number_points = number_points < 20 ? 20 : number_points;

   SPAposition *pos_array = new SPAposition [number_points];

   pos_array[0].set_x(from_vec.x());
   pos_array[0].set_y(from_vec.y());
   pos_array[0].set_z(from_vec.z());

   pos_array[number_points - 1].set_x(to_vec.x());
   pos_array[number_points - 1].set_y(to_vec.y());
   pos_array[number_points - 1].set_z(to_vec.z());

   for (int i=1; i < number_points-1; i++)
   {
      double theta = theta1 + total_angle * ((double)i/(double)number_points);
      double radius = xa * exp(xk * theta);
      double x = center.x() + radius * cos(theta);
      double y = center.y() + radius * sin(theta);

      pos_array[i].set_x(x);
      pos_array[i].set_y(y);
      pos_array[i].set_z(to_vec.z());
   }

     // calculate default tangent directions using up to three points
     // at the ends of the point array.
   SPAunit_vector uv1(0, 0, 0);
   SPAunit_vector uv2(0, 0, 0);

   EDGE* EDGE_ptr = NULL;
   outcome result;

     // fifth SPAparameter 0.0 -> interpolate the sequence of points.
   result = api_mk_ed_cubic(number_points, pos_array,
                            uv1, uv2, 0.0, EDGE_ptr);

     // Free memory
   delete [] pos_array;

   if ( !result.ok() || EDGE_ptr == NULL )
   {
      PRINT_ERROR("ACIS api_mk_ed_cubic function failed\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);

      if ( from == NULL || to == NULL )
      {
         PRINT_ERROR("In AGE:make_spiral_EDGE.\n"
                     "       Input VERTEX is NULL\n");
         return (EDGE *)NULL;
      }
   }

   else
   {
        // Switch the end VERTEXes of the EDGE to the input VERTEXes. Keep
        // in mind that the back links (from the VERTEX back to the EDGE)
        // need to be modified as well.
      EDGE_ptr->set_start ( from );
      from->add_edge( EDGE_ptr );

      EDGE_ptr->set_end ( to );
      to->add_edge( EDGE_ptr );
   }

   return EDGE_ptr;
}
//-------------------------------------------------------------------------
// Purpose       : Create a curve on a surface.
//
// Special Notes : If the input RefFace pointer is not NULL, then
//                 the points are first moved to the RefFace's surface
//                 before the interpolation is done.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 05/15/97
//-------------------------------------------------------------------------
EDGE* AcisModifyEngine::make_surface_EDGE( VERTEX* from_VERTEX,
                                             VERTEX* to_VERTEX,
                                             FACE *FACE_ptr,
                                             const CubitVector &plane_normal,
                                             const CubitVector *third_point) const
{
    //Now create a cutting plane using the normal passed in.  This
    //normal should point sort of tangential to the surface.  We
    //want the plane to intersect the surface.
  SPAposition const plane_point = from_VERTEX->geometry()->coords();
  SPAunit_vector const normal_vector( plane_normal.x(), plane_normal.y(),
                                   plane_normal.z() );
  plane const new_plane( plane_point, normal_vector );
  SPAposition point_1, point_2, point_3;
  SPApar_pos const pos_1( 0, 1 );
  SPApar_pos const pos_2( 1, 0 );
  SPApar_pos const pos_3( 0, 0 );
  new_plane.eval( pos_1, point_1 );
  new_plane.eval( pos_2, point_2 );
  new_plane.eval( pos_3, point_3 );

    // Get the size of the model so we can create an "infinite"
    // cutting plane
  CubitBox cubit_super_box = GeometryQueryTool::instance()->model_bounding_box();
  SPAbox super_box = AcisQueryEngine::bounding_box(cubit_super_box);
  CubitVector vec_1( point_1.x(), point_1.y(), point_1.z() );
  CubitVector vec_2( point_2.x(), point_2.y(), point_2.z() );
  CubitVector vec_3( point_3.x(), point_3.y(), point_3.z() );

    // Create the infinite plane, get a FACE
  BODY *sheet_body = this->create_infinite_plane_cutting_tool(
    vec_1, vec_2, vec_3, super_box, CUBIT_TRUE );
  if (sheet_body == NULL ||
      sheet_body->lump() == NULL ||
      sheet_body->lump()->shell() == NULL ||
      sheet_body->lump()->shell()->first_face() == NULL ||
      sheet_body->lump()->shell()->first_face()->geometry() == NULL )
  {
      //webcut_failed(refbody_cleanup_list, BODY_cleanup_list);
    PRINT_ERROR("In AcisModifyEngine::make_surface_EDGE"
                "       Cannot create infinite plane cutting tool using \n"
                "       the give surface. \n");
    return (EDGE*)NULL;
  }
  FACE *FACE_1 = sheet_body->lump()->shell()->first_face();
    // Copy the FACE
  FACE *FACE_2 = make_FACE( FACE_ptr );

    // Now we can create a body representing the intersection
    // of the infinite plane and the passed in surface
  EDGE* new_EDGE = NULL;
  BODY *intersect_BODY = NULL;
  LUMP *cur_LUMP = NULL;
  EDGE* cur_EDGE;
  COEDGE *cur_COEDGE;
  WIRE *cur_WIRE;

    // Get the intersection
  outcome result_1 = api_fafa_int ( FACE_1, FACE_2, intersect_BODY );

    // for some reason a NULL body still returns an ok result. - KGM
  if (intersect_BODY == NULL)
  {
    return (EDGE*)NULL;
  }

    // Make sure it was OK, and get the appropriate entity
    // Get a WIRE if it's there
  if ( !result_1.ok() || intersect_BODY->wire() == NULL ||
       intersect_BODY->wire()->coedge() == NULL ||
       intersect_BODY->wire()->coedge()->edge() == NULL )
  {
      // If there was no WIRE, see if there is a lump
    if ( result_1.ok() && intersect_BODY->lump() != NULL &&
         intersect_BODY->lump()->shell() != NULL &&
         intersect_BODY->lump()->shell()->wire() != NULL &&
         intersect_BODY->lump()->shell()->wire()->coedge() != NULL &&
         intersect_BODY->lump()->shell()->wire()->coedge()->edge() != NULL )
    {
      PRINT_DEBUG_18( "Using wire from lumps.\n");
      cur_LUMP = intersect_BODY->lump();
      cur_WIRE = cur_LUMP->shell()->wire();
    }
    else
    {
      PRINT_ERROR("ACIS api_fafa_int function failed\n");
      AcisQueryEngine::instance()->ACIS_API_error (result_1);
      return (EDGE*)NULL;
    }
  }
  else
  {
    PRINT_DEBUG_18( "Using wire straight from BODY\n");
    cur_WIRE = intersect_BODY->wire();
  }

    // Now find an EDGE that the points lie on.
    // Start by going through each LUMP, if LUMPs are being used.
  LUMP* start_LUMP = cur_LUMP;
  CURVE *geometry = NULL;
  do
  {
    PRINT_DEBUG_18( "**New Lump check**\n");

      // Now go through the wires, coedges, and edges
      // to find an edge that the points lie on
    WIRE* start_WIRE = cur_WIRE;
    do
    {
        // Get the first coedge in this wire
      COEDGE* start_COEDGE = cur_WIRE->coedge();
      cur_COEDGE = start_COEDGE;
        // Check it and loop through the rest
      do
      {
        PRINT_DEBUG_18( "  **New CoEdge check**\n");
        cur_EDGE =  cur_COEDGE->edge();
        geometry = cur_EDGE->geometry();
          // See if this is close enough to
          // the start and end points
        if (geometry != NULL)
        {
          SPAposition new_point;
            // Make sure start vertex is on curve
          geometry->equation().point_perp(from_VERTEX->geometry()->coords(),
                                          new_point);
          PRINT_DEBUG_18( "\n\tNew Point 1: (%g, %g, %g)\n"
                     "\t  From Vert: (%g, %g, %g)\n",
                     new_point.x(),
                     new_point.y(),
                     new_point.z(),
                     from_VERTEX->geometry()->coords().x(),
                     from_VERTEX->geometry()->coords().y(),
                     from_VERTEX->geometry()->coords().z());
          if (AcisQueryEngine::instance()->about_spatially_equal(from_VERTEX->geometry()->coords(),
                                          new_point, 1.0))
          {
              // Make sure end vertex is on curve
            geometry->equation().point_perp(to_VERTEX->geometry()->coords(),
                                            new_point);
            PRINT_DEBUG_18( "\tNew Point 2: (%g, %g, %g)\n"
                       "\t    To Vert: (%g, %g, %g) \n",
                       new_point.x(),
                       new_point.y(),
                       new_point.z(),
                       to_VERTEX->geometry()->coords().x(),
                       to_VERTEX->geometry()->coords().y(),
                       to_VERTEX->geometry()->coords().z());
            if (AcisQueryEngine::instance()->about_spatially_equal(to_VERTEX->geometry()->coords(),
                                            new_point, 1.0))
            {
                // You found it!!!
              cur_WIRE = NULL;
              start_LUMP = NULL;
              break;
            }
          }
        }
        cur_COEDGE = cur_COEDGE->next();
        geometry = NULL; // Make sure you indicate we haven't found geom
      }while (cur_COEDGE != NULL && cur_COEDGE != start_COEDGE);
        // You break out of this do-while loop when
        //   a) You find the correct EDGE
        //   b) There are no more COEDGEs in this WIRE
        //   c) You loop through all of the COEDGES in this WIRE

        // Print out some debug diagnostics
      if (DEBUG_FLAG(18))
      {
        if (cur_COEDGE == NULL)
          PRINT_INFO("Exit loop because CoEdge was NULL\n");
        else if (cur_COEDGE == start_COEDGE)
          PRINT_INFO("Exit loop because cur_COEDGE == start_COEDGE\n");
        else
          PRINT_INFO("Exit loop for Other Reason!!!\n");
      }

        // cur_WIRE is NULL if the correct EDGE has been found
      if (cur_WIRE)
        cur_WIRE = cur_WIRE->next();
    }while (cur_WIRE && cur_WIRE != start_WIRE);
      // You break out of this do-while loop when
      //   a) You find the correct EDGE
      //   b) There are no more WIREs in this WIRE list
      //   c) You loop through all the WIREs

      // Print out some debug diagnostics
    if (DEBUG_FLAG(18))
    {
      if (cur_WIRE == NULL)
        PRINT_DEBUG_18( "Exit loop because Wire was NULL\n");
      else if (cur_WIRE == start_WIRE)
        PRINT_DEBUG_18( "Exit loop because cur_WIRE == start_WIRE\n");
      else
        PRINT_DEBUG_18( "Exit loop for Other Reason!!!\n");
    }

      // If we are using lumps, go to the next lump
      // and get its first wire.
    if (start_LUMP)
    {
      cur_LUMP = cur_LUMP->next();
      if (cur_LUMP)
        cur_WIRE = cur_LUMP->shell()->wire();
    }

  }while (start_LUMP && cur_LUMP && start_LUMP != cur_LUMP);
    // You break out of this do-while loop when
    //   a) You find the correct EDGE
    //   b) There are no more LUMPs in the BODY
    //   c) You looped through all the LUMPs in the BODY
    //   d) We weren't using LUMPs to start with

    // Print out some debug diagnostics
  if (DEBUG_FLAG(18))
  {
    if (start_LUMP == NULL)
      PRINT_INFO("Exit loop because start_LUMP was NULL\n");
    else if (cur_LUMP == NULL)
      PRINT_INFO("Exit loop because cur_LUMP == NULL\n");
    else if (cur_LUMP == start_LUMP)
      PRINT_INFO("Exit loop because cur_WIRE == start_WIRE\n");
    else
      PRINT_INFO("Exit loop for Other Reason!!!\n");
  }

    // Make sure you found something valid
  if ( geometry == NULL )
  {
    PRINT_ERROR("ACIS failed to find geometry for new curve.\n");
    return (EDGE*)NULL;
  }

  PRINT_DEBUG_18( " The EDGE coords:\n"
             "   start: (%g, %g, %g)\n"
             "     end: (%g, %g, %g)\n",
             cur_EDGE->end_pos().x(),
             cur_EDGE->end_pos().y(),
             cur_EDGE->end_pos().z(),
             cur_EDGE->start_pos().x(),
             cur_EDGE->start_pos().y(),
             cur_EDGE->start_pos().z());

    // Handle special stuff for periodic curves
  int periodic = geometry->equation().periodic();
  double period = 0.0;
  if ( periodic )
  {
    period = geometry->equation().param_period();
  }

  double start_param, end_param;
  start_param = cur_EDGE->start_param();
  end_param = cur_EDGE->end_param();
  double curve_range = end_param-start_param;

  double v1_param, v2_param;
  SPAposition const v1_pos = from_VERTEX->geometry()->coords();
  SPAposition const v2_pos = to_VERTEX->geometry()->coords();
  v1_param = geometry->equation().param( v1_pos );
  v2_param = geometry->equation().param( v2_pos );

  // If periodic get the params in the range of the curve's param range.
  if(periodic)
  {
    // forward param range
    if(curve_range > 0.0)
    {
      while(v1_param < start_param - GEOMETRY_RESABS)
        v1_param += period;
      while(v1_param > end_param + GEOMETRY_RESABS)
        v1_param -= period;
      while(v2_param < start_param - GEOMETRY_RESABS)
        v2_param += period;
      while(v2_param > end_param + GEOMETRY_RESABS)
        v2_param -= period;
    }
    // reversed param range (will this ever happen?)
    else
    {
      while(v1_param < end_param - GEOMETRY_RESABS)
        v1_param += period;
      while(v1_param > start_param + GEOMETRY_RESABS)
        v1_param -= period;
      while(v2_param < end_param - GEOMETRY_RESABS)
        v2_param += period;
      while(v2_param > start_param + GEOMETRY_RESABS)
        v2_param -= period;
    }
  }

  double new_param_range = v2_param-v1_param;

  // I am using an integer to start with here just in 
  // case FORWARD and REVERSED are not defined to
  // be what we think.
  int sense_int = 1;
  if(curve_range*new_param_range < 0.0)
    sense_int = 0;  /* REVERSED */

    //if we have a third point and we are periodic... use
    //it to nail down the direction.
  if ( third_point != NULL && periodic )
  {
    SPAposition const some_point( third_point->x(),
                               third_point->y(),
                               third_point->z() );
    SPAposition closest_point;
    geometry->equation().point_perp( some_point, closest_point );
    double some_param = geometry->equation().param( closest_point );
    int val1 = (curve_range > 0.0);
    int val2 = (sense_int == 1);

    // val1 will equal val2 for two cases.  In either case v1_param will
    // be less than v2_param.
    if (val1 == val2)
    {
      if ( some_param < v1_param || some_param > v2_param )
        sense_int = !sense_int;
    }
    // val1 will NOT equal val2 for two cases.  In either case v1_param will
    // be greater than v2_param.
    else
    {
      if ( some_param < v2_param || some_param > v1_param )
        sense_int = !sense_int;
    }
  }

  REVBIT sense;
  if(sense_int == 0)
    sense = REVERSED;
  else
    sense = FORWARD;

    // Actually create the new EDGE
  API_BEGIN;
  new_EDGE = new EDGE ( from_VERTEX, to_VERTEX, geometry, sense );
  API_END;

    //Delete the body and the sheet_body and the other body associated with
    //this other face...
  AcisQueryEngine::instance()->delete_ACIS_BODY(intersect_BODY, CUBIT_TRUE);
  AcisQueryEngine::instance()->delete_ACIS_BODY(sheet_body, CUBIT_TRUE);
  AcisQueryEngine::instance()->delete_ACIS_BODY(AcisQueryEngine::instance()->get_BODY_of_ENTITY(FACE_2), CUBIT_TRUE);

  return new_EDGE;
}

//-------------------------------------------------------------------------
// Purpose       : Create a spline EDGE given a set of ACIS SPAposition objects.
//
// Special Notes : If the input RefFace pointer is not NULL, then
//                 the points are first moved to the RefFace's surface
//                 before the interpolation is done.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 05/15/97
//-------------------------------------------------------------------------
EDGE* AcisModifyEngine::make_spline_EDGE( VERTEX* from_VERTEX,
                                            VERTEX* to_VERTEX,
                                            SPAposition* pos_array,
                                            int number_points ) const
{
    // Create the spline curve and attach it to the EDGE
  EDGE* EDGE_ptr = NULL;
  outcome result = api_curve_spline(number_points, pos_array, 0, 0, EDGE_ptr,
                                    FALSE, from_VERTEX == to_VERTEX);
  if (result.ok() == FALSE || EDGE_ptr == NULL)
  {
    PRINT_ERROR("In AcisModifyEngine::make_spline_EDGE\n"
                "       ACIS function, api_mk_ed_cubic, failed\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    return (EDGE *)NULL;
  }
  else
  {
      // Switch the end VERTEXes of the EDGE to the input VERTEXes. Keep
      // in mind that the back links (from the VERTEX back to the EDGE)
      // need to be modified as well.
    EDGE_ptr->set_start ( from_VERTEX );
    from_VERTEX->add_edge( EDGE_ptr );

    EDGE_ptr->set_end ( to_VERTEX );
    to_VERTEX->add_edge( EDGE_ptr );

    return EDGE_ptr;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Boolean Operation of a list of Bodies: chop
//
// Special Notes :
//
// Creator       : Lingyun Pan, Cat
//
// Creation Date : 7/13/01
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::chop(DLIList<BodySM*> &bodysms,
                                   DLIList<BodySM*> &intersectBodies,
                                   DLIList<BodySM*> &outsideBodies,
                                   BodySM*& leftoverBody,
                                   bool keep_old ,
                                   bool nonreg) const
{
  leftoverBody = 0;

    // chops the blank with the  tool, returing the body formed by subtracting the tool from the blank,
	// and the body formed by intersecting the tool with the blank,
	// simultaneously.

  if (bodysms.size() <= 1)
  {
    PRINT_WARNING("There is only one volume in the list. Nothing modified\n");
    return CUBIT_FAILURE;
  }

    // Create new acis bodies that are copies of the original acis bodies.
    // If something goes wrong in the choping, we can then delete our copies
    // and the originals are intact.  If everything goes fine, we will delete
    // the originals at that time.

    // pass in keep_old to the delete_owner_attrib flag; if we're not keeping
    // old bodies, we want to keep the owner attrib, so we can pick up entities
    // that didn't change
  bool delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old);

  DLIList<BODY*> old_BODY_list;
  DLIList<BodySM*> old_BodySM_list;

  BodySM* BodySMPtr = bodysms.get_and_step();
  BODY* BODYPtr =  AcisQueryEngine::get_BODY(BodySMPtr);
  old_BodySM_list.append(BodySMPtr);
  old_BODY_list.append(BODYPtr);

  BODY *blank_BODY = this->copy_BODY(BODYPtr, delete_attribs);
  if (blank_BODY == NULL)
  {
    return CUBIT_FAILURE; // Nothing to clean up at this point.
  }
  
  BODY* outside_BODY=NULL;
  BODY* leftovers_BODY=NULL;
 // get the tool Body
  BodySMPtr = bodysms.get_and_step();
  BODYPtr =  AcisQueryEngine::get_BODY(BodySMPtr);

  BODY *tool_BODY = this->copy_BODY(BODYPtr, delete_attribs);
  if (tool_BODY == NULL)
  {
    AcisQueryEngine::instance()->delete_ACIS_BODY(blank_BODY, CUBIT_TRUE);
    return CUBIT_FAILURE;
  }

  // Chop the blank_BODY using the tool_BODY.
  // If this is successful, the result is blank_BODY and outside_BODY, leftovers_BODY
  outcome result = api_boolean_chop_body(tool_BODY, blank_BODY, nonreg,outside_BODY, leftovers_BODY);

  if (!result.ok() || blank_BODY == NULL ||  
       tool_BODY == NULL || blank_BODY->lump() == NULL )
  {


    if( !result.ok() )
    {
      PRINT_ERROR("In AcisModifyEngine::chop\n");
      AcisQueryEngine::instance()->ACIS_API_error(result, "chop Bodies");
      if (blank_BODY != NULL) 
        AcisQueryEngine::instance()->delete_ACIS_BODY(blank_BODY, CUBIT_TRUE);
      if (tool_BODY != NULL) 
        AcisQueryEngine::instance()->delete_ACIS_BODY(tool_BODY, CUBIT_TRUE);
    }

    if( blank_BODY->lump() == NULL )
    {
      PRINT_ERROR("Bodies do not intersect\n");
      if (blank_BODY != NULL) 
        AcisQueryEngine::instance()->delete_ACIS_BODY(blank_BODY, CUBIT_TRUE);
      if (leftovers_BODY != NULL) 
        AcisQueryEngine::instance()->delete_ACIS_BODY(leftovers_BODY, CUBIT_TRUE);
      if (outside_BODY != NULL) 
        AcisQueryEngine::instance()->delete_ACIS_BODY(outside_BODY, CUBIT_TRUE);
    }

    return CUBIT_FAILURE;
  }

  BodySM* intersectBody = get_new_Body(BodySMPtr, BODYPtr, blank_BODY, keep_old);
  intersectBodies.append( intersectBody );
  BodySM* outsideBody =  get_new_Body(old_BodySM_list, old_BODY_list, outside_BODY, keep_old);
  outsideBodies.append( outsideBody );

  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : To flip the ACIS normals of a face
//
// Special Notes :
//
// Creator       : Samuel Showman (CAT)
//
// Creation Date : 10/16/02
//-------------------------------------------------------------------------
CubitStatus
AcisModifyEngine::flip_normals( DLIList<Surface*>& ref_face_list ) const

{
  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() != 0);

  outcome result;
  int i;

  DLIList<SurfaceACIS*> copied_ref_face_list(ref_face_list.size());
  CAST_LIST(ref_face_list, copied_ref_face_list, SurfaceACIS);
  if (ref_face_list.size() != copied_ref_face_list.size())
  {
    PRINT_ERROR("Non-ACIS surface at %s:%d\n", __FILE__, __LINE__ );
    return CUBIT_FAILURE;
  }

  // Loop on FACEs.  We will work on reversing surfaces from one body at a time.
  copied_ref_face_list.reset();
  while( copied_ref_face_list.size() )
  {
    BODY *copied_BODY_ptr;
#ifdef BOYD07
    DLIList<SurfaceACIS*> reverse_face_list;
#endif
    DLIList<FACE*> reverse_FACE_list;
    if( get_copied_FACES_of_body( copied_ref_face_list, reverse_FACE_list,
      copied_BODY_ptr ) == CUBIT_FAILURE )
      break;

    // Get original Body and BODY
    BodySM *body_ptr = AcisQueryEngine::instance()->get_body_sm_of_ENTITY( copied_BODY_ptr );
    BODY *BODY_ptr = AcisQueryEngine::get_BODY(body_ptr);

    // Now cleanout the owner attributes from the copied BODY, if required
    if( delete_attribs )
      AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);

    reverse_FACE_list.reset();
    for( i=reverse_FACE_list.size(); i--; )
    {
      FACE *FACE_ptr = reverse_FACE_list.get_and_step();
      result = api_reverse_face( FACE_ptr );

      if (!result.ok())
      {
        //TODO: get ref_face from FACE
        //PRINT_INFO("Surface %d did not reverse\n", ref_face->id());
         PRINT_INFO("Surface did not reverse\n");
        return CUBIT_FAILURE;
      }
    }

    BodySM* new_body_ptr = get_new_Body( body_ptr, BODY_ptr, copied_BODY_ptr, CUBIT_FALSE );

    if( new_body_ptr!=NULL && new_body_ptr!=body_ptr )
    {
      PRINT_INFO( "Created new volume\n" );
    }
    else if( new_body_ptr!=NULL && new_body_ptr==body_ptr )
    {
      PRINT_INFO( "Modified volume\n" );
    }
    else
      PRINT_WARNING( "Volume was not modified\n" );
  }

  return CUBIT_SUCCESS;
}
//-------------------------------------------------------------------------
// Purpose       : To thicken sheet bodies
//
// Special Notes :
//
// Creator       : Samuel Showman (CAT)
//
// Creation Date : 10/14/02
//-------------------------------------------------------------------------
#if CUBIT_ACIS_VERSION < 800
CubitStatus AcisModifyEngine::thicken( DLIList<BodySM*>& ,
                                       DLIList<BodySM*>&,
                                       double ,
                                       bool ) const
{
  PRINT_ERROR( "This function is not available in ACIS versions below 800\n" );
  return CUBIT_FAILURE;
#else
CubitStatus AcisModifyEngine::thicken( DLIList<BodySM*>& bodies,
                                       DLIList<BodySM*>& new_body_list,
                                       double depth,
                                       bool both) const
{
   if (bodies.size() < 1)
   {
      PRINT_WARNING("There are no volumes in the list. Nothing modified\n");
      return CUBIT_FAILURE;
   }

   SPAposition box_l(0,0,0);
   SPAposition box_h(0,0,0);

   BODY* in_out_BODY;
   DLIList<BODY*> old_BODY_list;
   DLIList<BodySM*> old_Body_list;

   for(int i = bodies.size(); i >= 1 ; i-- ){

      BodySM* BodyPtr = bodies.get_and_step();
      BodySM* new_bodysm;
      BodyACIS* acis_body = dynamic_cast<BodyACIS*>(BodyPtr);
      if (!acis_body)
      {
        PRINT_ERROR("Non-ACIS BodySM passed to AcisModifyEngine::thicken.\n");
        continue;
      }

      BODY* BODYPtr =  AcisQueryEngine::get_BODY(BodyPtr);

      in_out_BODY = this->copy_BODY(BODYPtr, CUBIT_FALSE);

      if (in_out_BODY == NULL)
      {
         return CUBIT_FAILURE; // Nothing to clean up at this point.
      }

      if( !acis_body->is_sheet_body() )
      {
         PRINT_WARNING("Only surface volumes should be thickened\n"
            "       results may be different than expected\n");
      }

      outcome result = api_sheet_thicken(in_out_BODY, depth, both, box_l, box_h);

      if (!result.ok())//
      {
         if(result.error_number() == 41000)               // checking surface normals
         {
            PRINT_ERROR("In AcisModifyEngine::thicken\n"
                "       Surface normals may not be consistent\n"
                "       use: validate normal reverse\n");
            return CUBIT_FAILURE;
         }

         PRINT_ERROR("In AcisModifyEngine::thicken \n");
         AcisQueryEngine::instance()->ACIS_API_error (result);
         return CUBIT_FAILURE;
      }

      new_bodysm = get_new_Body(BodyPtr, BODYPtr, in_out_BODY, CUBIT_FALSE);
      if (new_bodysm)
        new_body_list.append(new_bodysm);
  }
  return CUBIT_SUCCESS;
#endif
}


//-------------------------------------------------------------------------
// Purpose       : Boolean Operation of a list of Bodies: unite
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/31/96
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::unite(DLIList<BodySM*> &bodies,
                                    DLIList<BodySM*> &newBodies,
                                      bool keep_old) const
{
    // Union all bodies in the bodies list into a single body.
  if (bodies.size() <= 1)
  {
    PRINT_WARNING("There is only one volume in the list. Nothing modified\n");
    return CUBIT_FAILURE;
  }

    // Create new acis bodies that are copies of the original acis bodies.
    // If something goes wrong in the unioning, we can then delete our copies
    // and the originals are intact.  If everything goes fine, we will delete
    // the originals at that time.

    // pass in keep_old to the delete_owner_attrib flag; if we're not keeping
    // old bodies, we want to keep the owner attrib, so we can pick up entities
    // that didn't change
   bool delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old);

   // check if the boolean operation is regularized or nonregularized
     BOOL_TYPE bool_type = UNION;
	CubitBoolean boolean_regularize = GeometryModifyTool::instance()->boolean_regularize();

	if (boolean_regularize == FALSE)
	{
		bool_type = NONREG_UNION;
	}

  DLIList<BODY*> old_BODY_list;
  DLIList<BodySM*> old_Body_list;

  BodySM* BodyPtr = bodies.get_and_step();
  BODY* BODYPtr =  AcisQueryEngine::get_BODY(BodyPtr);
  old_Body_list.append(BodyPtr);
  old_BODY_list.append(BODYPtr);

  BODY *master = this->copy_BODY(BODYPtr, delete_attribs);
  if (master == NULL)
  {
    return CUBIT_FAILURE; // Nothing to clean up at this point.
  }


    // Note: one iteration less because we already pulled out one body above
  for (int i = 1; i < bodies.size(); i++)
  {
    BodyPtr = bodies.get_and_step();
    BODYPtr =  AcisQueryEngine::get_BODY(BodyPtr);

    BODY *copy = this->copy_BODY(BODYPtr, delete_attribs);
    if (copy == NULL)
    {
      AcisQueryEngine::instance()->delete_ACIS_BODY(master, CUBIT_TRUE);
      return CUBIT_FAILURE;
    }


	// Do the union of the master and the copy.
    // If this is successful, the result is master and copy will be deleted
    // outcome result = api_unite(copy, master);
	// replaced with the api_boolean call in order to do nonregularize boolean operation

	 outcome result = api_boolean (copy,master,bool_type);



    if (!result.ok() || (master && master->lump() == NULL))
    {
      // If there are no LUMPs in the resulting BODY, return with a failure
      if (master && master->lump() == NULL)
        PRINT_ERROR("Unite produced volume with no volume - "
                    "check for (and remove) co-incident surfaces\n");
      else
        PRINT_ERROR("In AcisModifyEngine::unite\n");

      AcisQueryEngine::instance()->ACIS_API_error(result, "unite Bodies");
      if (master != NULL) AcisQueryEngine::instance()->delete_ACIS_BODY(master, CUBIT_TRUE);
      if (copy   != NULL) AcisQueryEngine::instance()->delete_ACIS_BODY(copy, CUBIT_TRUE);

      return CUBIT_FAILURE;
    }
    else {
      old_Body_list.append(BodyPtr);
      old_BODY_list.append(BODYPtr);
    }

		//For nonregularized boolean, find all internal surfaces and delete

	if ( bool_type == NONREG_UNION)
	{
		  outcome result1;
		  // ENTITY_LIST FACES;
          DLIList<FACE*>  FACE_list ;
		  FACE *this_FACE;
          AcisQueryEngine::instance()->get_FACEs(master,FACE_list);
		   // loop through all faces
		  for( int i = 0; i< FACE_list.size(); i++)
		  {
			  this_FACE = FACE_list.get_and_step();
			  assert( this_FACE != NULL );
			  // Make sure this is a DOUBLE_SIDED FACE
			  if (this_FACE->sides() == DOUBLE_SIDED)
			  {
				  if (this_FACE->cont()==BOTH_INSIDE)
				  {
					  // found internal faces
					  // Now we have the FACE - unhook it  from the BODY.  Keep track of new
					  // BODIES that are created as this is done.
					   PRINT_INFO( " Unhooking and deleting each internal surface...\n" );
					   BODY *new_BODY_ptr;
					   result1 = api_unhook_face( this_FACE, new_BODY_ptr );
					   if( !result1.ok() )
					   {
						   AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result1 );
						   PRINT_ERROR( " Face unhooking during rebuild of volume didn't work\n" );
						   return CUBIT_FAILURE;
					   }
					   AcisQueryEngine::instance()->delete_ACIS_BODY (new_BODY_ptr);
				  }
			  }
		  }

		 // heal the leftover body
//
//		   PRINT_INFO(" Healing the leftover volume...\n");
//		   if( AcisHealerTool::instance()->init_BODY_for_healing( master ) == CUBIT_SUCCESS )
//		   {
//			   int percent_before, percent_after, number_splines_simplified;
//			   if( AcisHealerTool::instance()->heal_BODY( master, percent_before,
//				   percent_after, number_splines_simplified ) == CUBIT_FAILURE )
//				   PRINT_ERROR( "Error healing the combined volume\n" );
//			   else
//				   PRINT_INFO( "Successfully healed the combined volume!\n" );
//			   AcisHealerTool::instance()->end_BODY_for_healing( master );
//		   }
	}

  }

  AcisQueryEngine::instance()->clear_bounding_box( master );
  AcisQueryEngine::instance()->bounding_box( master );

  BodySM *newBodyPtr = get_new_Body(old_Body_list, old_BODY_list, master, keep_old);
  newBodies.append( newBodyPtr );

  return CUBIT_SUCCESS;
}

CubitStatus AcisModifyEngine::subtract( DLIList<BodySM*> &tool_body_list,
                                        DLIList<BodySM*> &from_bodies,
                                        DLIList<BodySM*> &new_from_bodies,
                                        bool imprint,
                                        bool keep_old) const
{
    // pass in keep_old to the delete_owner_attrib flag; if we're not keeping
    // old bodies, we want to keep the owner attrib, so we can pick up entities
    // that didn't change
  bool delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old);

  // check if the boolean operation is regularized or nonregularized
     BOOL_TYPE bool_type = SUBTRACTION;
	CubitBoolean boolean_regularize = GeometryModifyTool::instance()->boolean_regularize();

	if (boolean_regularize == FALSE)
	{
		bool_type = NONREG_SUBTRACTION;
	}

  int ii;
  DLIList<CubitBox*> tool_boxes;

  CubitBox tool_box;


  DLIList<BODY*> tool_BODY_list;
  DLIList<BODY*> tool_BODY_list_copy;

  // get acis bodies, and copy them
  tool_body_list.reset();

  for (ii = tool_body_list.size(); ii > 0; ii--) {
    CubitBox *temp_box = new CubitBox(
	    AcisQueryEngine::instance()->bounding_box(tool_body_list.get()));
    tool_boxes.append(temp_box);

    BODY *BODYPtr1 = AcisQueryEngine::get_BODY(tool_body_list.get_and_step());
    tool_BODY_list.append(BODYPtr1);

    BODYPtr1 = this->copy_BODY(BODYPtr1, delete_attribs);
    tool_BODY_list_copy.append(BODYPtr1);
  }

  if (tool_BODY_list.size() != tool_body_list.size())
  {
    for(ii = tool_boxes.size(); ii>0; ii--)
       delete tool_boxes.get_and_step();
    return CUBIT_FAILURE;
  }


  DLIList<BODY*> from_BODY_list;
  DLIList<BODY*> from_BODY_list_copy;

  // get acis bodies, and copy them
  from_bodies.reset();

  for (ii = from_bodies.size(); ii > 0; ii--) {
    BODY *BODYPtr1 = AcisQueryEngine::get_BODY(from_bodies.get_and_step());
    from_BODY_list.append(BODYPtr1);

    BODYPtr1 = this->copy_BODY(BODYPtr1, delete_attribs);
    from_BODY_list_copy.append(BODYPtr1);
  }

  if (from_BODY_list.size() != from_bodies.size())
  {
    for(ii = tool_boxes.size(); ii>0; ii--)
       delete tool_boxes.get_and_step();
    return CUBIT_FAILURE;
  }


  // now, subtract the tool from the list of bodies
  tool_BODY_list.reset();
  tool_BODY_list_copy.reset();
  tool_body_list.reset();
  tool_boxes.reset();

  from_BODY_list.reset();
  from_BODY_list_copy.reset();
  from_bodies.reset();

  int fraction_remaining = 10;
  CubitStatus sub_stat = CUBIT_SUCCESS;

    // subtract the tool body from each body in the list
  CubitMessage* cmi = CubitMessage::instance();
  for (ii = 1; ii <= from_BODY_list.size(); ii++)
  {
    BODY *from_BODY = from_BODY_list.get();
    BODY *from_BODY_copy = from_BODY_list_copy.get();
    BodySM *from_Body = from_bodies.get();
    CubitBox box1 = AcisQueryEngine::instance()->bounding_box(from_Body);

    for(int jj = tool_BODY_list.size(); jj>0; jj--)
    {
      if (cmi->Interrupt())
      {
        PRINT_ERROR("Subtraction interrupted.  Aborting...\n");
        while (tool_boxes.size())
          delete tool_boxes.pop();
        while (tool_BODY_list_copy.size())
          AcisQueryEngine::instance()->delete_ACIS_BODY(tool_BODY_list_copy.pop(),CUBIT_TRUE);
        while (from_BODY_list_copy.size())
          AcisQueryEngine::instance()->delete_ACIS_BODY(from_BODY_list_copy.pop(),CUBIT_TRUE);
        return CUBIT_FAILURE;
      }

      tool_box = *tool_boxes.get();
      BODY* tool_BODY = tool_BODY_list.get();
      BODY *tool_BODY_copy = tool_BODY_list_copy.get();

        // first, check bounding SPAbox; if they don't intersect, don't do the subtract
      if ( tool_box.overlap( SPAresabs, box1 ) )
      {
          // the bodies overlap; proceed with the subtract

          // Subtract body1 from body2.
         // outcome result = api_subtract( tool_BODY_copy, from_BODY_copy );
        outcome result = api_boolean( tool_BODY_copy, from_BODY_copy, bool_type );
          // We may or may not get an error if the resulting BODY is
          // empty (body1 totally encloses body2).
          // Also check the result.
        if (from_BODY_copy == NULL ||
            (from_BODY_copy->lump() == NULL && from_BODY_copy->wire() == NULL) ||
            !result.ok())
        {
          if (!result.ok())
          {
            PRINT_ERROR("Subtraction operation failed.\n");
            AcisQueryEngine::instance()->ACIS_API_error(result, "subtract Bodies");
            if (tool_BODY_copy != NULL)
            {
              AcisQueryEngine::instance()->delete_ACIS_BODY(tool_BODY_copy, CUBIT_TRUE);
            }
          }
          else
             PRINT_ERROR("Subtraction operation failed.\n"
                         "       Empty volume created. Original volume is probably\n"
                         "       completely enclosed by the subtracting volume.\n");
          if (from_BODY_copy != NULL)
          {
            AcisQueryEngine::instance()->delete_ACIS_BODY(from_BODY_copy, CUBIT_TRUE);
          }

            // we had an error, so copy the original body into the copy list
          from_BODY_copy = this->copy_BODY(from_BODY);
          from_BODY_list_copy.change_to(from_BODY_copy);
          tool_BODY_copy = this->copy_BODY(tool_BODY);
          tool_BODY_list_copy.change_to(tool_BODY_copy);
          sub_stat = CUBIT_FAILURE;
        }
        else
        {
          // else subtract was successful; make another tool body
          tool_BODY_copy = this->copy_BODY(tool_BODY);
          if( imprint )
            imprint_BODYs( tool_BODY_copy, from_BODY_copy );
          tool_BODY_list_copy.change_to(tool_BODY_copy);
        }
      
				//For nonregularized boolean, find all internal surfaces and delete

	if ( bool_type == NONREG_SUBTRACTION)
	{
		  outcome result1;
		  // ENTITY_LIST FACES;
          DLIList<FACE*>  FACE_list ;
		  FACE *this_FACE;
          AcisQueryEngine::instance()->get_FACEs(from_BODY_copy,FACE_list);
		   // loop through all faces
		  for( int i = 0; i< FACE_list.size(); i++)
		  {
			  this_FACE = FACE_list.get_and_step();
			  assert( this_FACE != NULL );
			  // Make sure this is a DOUBLE_SIDED FACE
			  if (this_FACE->sides() == DOUBLE_SIDED)
			  {
				  if (this_FACE->cont()==BOTH_INSIDE)
				  {
					  // found internal faces
					  // Now we have the FACE - unhook it  from the BODY.  Keep track of new
					  // BODIES that are created as this is done.
					   PRINT_INFO( " Unhooking and deleting each internal surface...\n" );
					   BODY *new_BODY_ptr;
					   result1 = api_unhook_face( this_FACE, new_BODY_ptr );
					   if( !result1.ok() )
					   {
						   AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result1 );
						   PRINT_ERROR( " Face unhooking during rebuild of volume didn't work\n" );
						   return CUBIT_FAILURE;
					   }
					   AcisQueryEngine::instance()->delete_ACIS_BODY (new_BODY_ptr);
				  }
			  }
		  }

		 // heal the leftover body
//
//		   PRINT_INFO(" Healing the leftover volume...\n");
//		   if( AcisHealerTool::instance()->init_BODY_for_healing( master ) == CUBIT_SUCCESS )
//		   {
//			   int percent_before, percent_after, number_splines_simplified;
//			   if( AcisHealerTool::instance()->heal_BODY( master, percent_before,
//				   percent_after, number_splines_simplified ) == CUBIT_FAILURE )
//				   PRINT_ERROR( "Error healing the combined volume\n" );
//			   else
//				   PRINT_INFO( "Successfully healed the combined volume.\n" );
//			   AcisHealerTool::instance()->end_BODY_for_healing( master );
//		   }
	}

           // done with this j iteration; write out count, if necessary
        if (from_bodies.size() * tool_body_list.size() > 1)
        {
          int frac_done = (10 * ii) / (from_bodies.size()* tool_body_list.size());
          if ((10 - frac_done) < fraction_remaining)
          {
            PRINT_INFO("%d ", fraction_remaining);
            fraction_remaining--;
          }
        }
      }
      tool_boxes.step();
      tool_BODY_list.step();
      tool_BODY_list_copy.step();

    }

      // done with iteration over ii; step the lists
    from_BODY_list_copy.step();
    from_BODY_list.step();
    from_bodies.step();
  }

    // ok, we're done with all the imprints; construct new Body's for the new BODY's
  from_BODY_list_copy.reset();
  from_BODY_list.reset();
  from_bodies.reset();

  for (ii = 0; ii< from_BODY_list_copy.size(); ii++)
  {
      BODY *old_BODY = from_BODY_list.get();
      BODY *new_BODY = from_BODY_list_copy.get();

      BodySM *new_body = NULL;
      if (old_BODY != new_BODY)
         new_body = get_new_Body(from_bodies.get(), old_BODY, new_BODY, keep_old);

      if (new_body)
      {
        new_from_bodies.append(new_body);
        from_bodies.change_to(NULL);
      }

      // now step all the lists
    from_BODY_list.step();
    from_BODY_list_copy.step();
    from_bodies.step();
  }

  from_bodies.remove_all_with_value(NULL);
  from_bodies.reset();
  new_from_bodies.reset();

  for(ii = tool_BODY_list_copy.size(); ii>0; ii--)
     AcisQueryEngine::instance()->delete_ACIS_BODY(tool_BODY_list_copy.get_and_step(), CUBIT_TRUE);

  if (sub_stat == CUBIT_SUCCESS && !keep_old)
    AcisQueryEngine::instance()->delete_solid_model_entities(tool_body_list);

  for(ii = tool_boxes.size(); ii>0; ii--)
     delete tool_boxes.get_and_step();

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Boolean Operation of two Bodys: imprint
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/30/96
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::imprint(BodySM* BodyPtr1,
                                      BodySM* BodyPtr2,
                                      BodySM*& newBodyPtr1,
                                      BodySM*& newBodyPtr2,
                                      bool keep_old) const
{
    // There are two many places things can go wrong. When they do, in
    // addition to returning a CUBIT_FAILURE, we should return some sensible
    // values through the output arguments. Before any of the output Body's can
    // be created, the most sensible for these is NULL.
  newBodyPtr1 = NULL ;
  newBodyPtr2 = NULL ;

    // first, check bounding SPAbox; if they don't intersect, imprint will fail
  AcisQueryEngine* aqe = AcisQueryEngine::instance();
  BODY* BODYPtr1 =  aqe->get_BODY(BodyPtr1);
  BODY* BODYPtr2 =  aqe->get_BODY(BodyPtr2);

  if( (BODYPtr1 == NULL) || (BODYPtr2 == NULL) )
  {
    return CUBIT_FAILURE;
  }

  CubitBox box1 = aqe->bounding_box(aqe->bounding_box(BODYPtr1));
  CubitBox box2 = aqe->bounding_box(aqe->bounding_box(BODYPtr2));
  if (! box1.overlap(SPAresabs, box2))
    return CUBIT_FAILURE;

    // Make copies of the ACIS bodies. This is required because if the ultimate
    // result is a body with no LUMPs, then we need to be able to restore
    // the original BODY's and their Bodies.
  BODY* BODYPtr1Save;
  BODY* BODYPtr2Save;

    // pass in keep_old to the delete_owner_attrib flag; if we're not keeping
    // old bodies, we want to keep the owner attrib, so we can pick up entities
    // that didn't change
  bool delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old);
  BODYPtr1Save = this->copy_BODY(BODYPtr1, delete_attribs);
  if (BODYPtr1Save == NULL)
  {
    PRINT_ERROR("Copy of volume failed.\n       Imprint failed.\n");
    return CUBIT_FAILURE;
  }
  BODYPtr2Save = this->copy_BODY(BODYPtr2, delete_attribs);
  if (BODYPtr2Save == NULL)
  {
    AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr1Save, CUBIT_TRUE);
    PRINT_ERROR("Copy of volume failed.\n       Imprint failed.\n");
    return CUBIT_FAILURE;
  }

    // Keep track of the number of curves/surfaces/vertices in each body.  That
    // way we can tell if they really changed.
  int num_vol1_bef, num_face1_bef, num_edge1_bef, num_vertex1_bef;
  if( AcisQueryEngine::instance()->number_ENTITIES( BODYPtr1Save, num_vol1_bef, num_face1_bef, num_edge1_bef,
                       num_vertex1_bef ) == CUBIT_FAILURE )
  {
     PRINT_ERROR( "Unable to count number of entities in volume\n" );
     return CUBIT_FAILURE;
  }
  int num_vol2_bef, num_face2_bef, num_edge2_bef, num_vertex2_bef;
  if( AcisQueryEngine::instance()->number_ENTITIES( BODYPtr2Save, num_vol2_bef, num_face2_bef, num_edge2_bef,
                       num_vertex2_bef ) == CUBIT_FAILURE )
  {
     PRINT_ERROR( "Unable to count number of entities in volume\n" );
     return CUBIT_FAILURE;
  }

    // Now, imprint the 2 BODYs
  if ( GeometryModifyTool::get_all_edges_imprint() )
  {
    api_set_int_option("all_free_edges", TRUE );
  }

  outcome result = api_imprint( BODYPtr1Save, BODYPtr2Save );

   //Gets rid of sliver curves/surfaces that could get produced by the imprint.
   AcisModifyEngine::instance()->cleanup_slivers( BODYPtr1Save ); 
   AcisModifyEngine::instance()->cleanup_slivers( BODYPtr2Save ); 

  if ( GeometryModifyTool::get_all_edges_imprint() )
  {
    api_set_int_option("all_free_edges", FALSE );
  }

  if (!result.ok())
  {
      // might just be no overlap, check smaller boxes
    if (!BODYPtr1Save || !BODYPtr2Save ||
        box1.overlap(SPAresabs,box2) == CUBIT_TRUE)
    {
      if (DEBUG_FLAG(95))
      {
        AcisQueryEngine::instance()->ACIS_API_error(result, "imprint Bodies");
      }

      if (BODYPtr1Save != NULL)
      {
        AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr1Save, CUBIT_TRUE);
      }
      if (BODYPtr2Save != NULL)
      {
        AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr2Save, CUBIT_TRUE);
      }
      return CUBIT_FAILURE;
    }
  }


    // Compare body 1 before & after imprint
  int num_vol1_aft, num_face1_aft, num_edge1_aft, num_vertex1_aft;
  if( AcisQueryEngine::instance()->
      number_ENTITIES( BODYPtr1Save, num_vol1_aft, num_face1_aft, num_edge1_aft,
                       num_vertex1_aft ) == CUBIT_FAILURE )
  {
     PRINT_ERROR( "Unable to count number of entities in volume\n" );
    if (BODYPtr1Save != NULL)
    {
      AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr1Save, CUBIT_TRUE);
    }
    if (BODYPtr2Save != NULL)
    {
      AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr2Save, CUBIT_TRUE);
    }
     return CUBIT_FAILURE;
  }

  // Only create a new Body1 if it changed
  if( num_vol1_bef != num_vol1_aft ||
      num_face1_bef != num_face1_aft ||
      num_edge1_bef != num_edge1_aft ||
      num_vertex1_bef != num_vertex1_aft ||
      GeometryModifyTool::instance()->get_group_imprint() == CUBIT_FALSE)
  {
      DLIList<EDGE*> new_edges = AcisModifyEngine::instance()->find_new_EDGES(BODYPtr2Save);

      newBodyPtr1 = AcisQueryEngine::instance()->populate_topology_bridges(BODYPtr1Save);
      if (newBodyPtr1 == NULL )
      {
          PRINT_ERROR("Problems creating topology for imprinted volume.\n") ;
          if (BODYPtr1Save != NULL)
          {
              AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr1Save, CUBIT_TRUE);
          }
          if (BODYPtr2Save != NULL)
          {
              AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr2Save, CUBIT_TRUE);
          }
          return CUBIT_FAILURE;
      }

      // Add a imprint feature to the topo edges
      for (int edge_count = new_edges.size(); edge_count--; ) 
      {
          CubitSimpleAttrib *tmp_attrib = new CubitSimpleAttrib( "SOURCE_FEATURE", "IMPRINT" );
          new ATTRIB_SNL_SIMPLE( new_edges[edge_count], tmp_attrib );
          delete tmp_attrib;
      }

     //BodyPtr1->switch_entity_names(newBodyPtr1);
     if (GeometryModifyTool::get_old_names())
     {
       Body* model_ent = dynamic_cast<Body*>(BodyPtr1->topology_entity());
       if (model_ent)
       {
         Body* new_model_ent = dynamic_cast<Body*>(newBodyPtr1->topology_entity());
         if (!new_model_ent)
           new_model_ent = GeometryQueryTool::instance()->make_Body(newBodyPtr1);
         model_ent->switch_entity_names(new_model_ent);
       }
     }

     if (keep_old == CUBIT_FALSE)
       AcisQueryEngine::instance()->delete_solid_model_entities(BodyPtr1);
  }
  else
  {
     AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr1Save, CUBIT_TRUE);
  }

  // Compare body 2 before & after imprint
  int num_vol2_aft, num_face2_aft, num_edge2_aft, num_vertex2_aft;
  if( AcisQueryEngine::instance()->number_ENTITIES( BODYPtr2Save, num_vol2_aft, num_face2_aft, num_edge2_aft,
                       num_vertex2_aft ) == CUBIT_FAILURE )
  {
     PRINT_ERROR( "Unable to count number of entities in volume.\n");
    if (BODYPtr1Save != NULL)
    {
      AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr1Save, CUBIT_TRUE);
    }
    if (BODYPtr2Save != NULL)
    {
      AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr2Save, CUBIT_TRUE);
    }
     return CUBIT_FAILURE;
  }

  // Only create a new Body2 if it changed
  if( num_vol2_bef != num_vol2_aft ||
      num_face2_bef != num_face2_aft ||
      num_edge2_bef != num_edge2_aft ||
      num_vertex2_bef != num_vertex2_aft ||
      GeometryModifyTool::instance()->get_group_imprint() == CUBIT_FALSE)
  {
      DLIList<EDGE*> new_edges = AcisModifyEngine::instance()->find_new_EDGES(BODYPtr2Save);

      newBodyPtr2 = AcisQueryEngine::instance()->populate_topology_bridges(BODYPtr2Save);
      if (newBodyPtr2 == NULL )
      {
          PRINT_ERROR("Problems creating topology for imprinted volume.\n" ) ;
          if (BODYPtr2Save != NULL)
          {
              AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr2Save, CUBIT_TRUE);
          }
          return CUBIT_FAILURE;
      }

      // Add a imprint feature to the topo edges
      for (int edge_count = new_edges.size(); edge_count--; ) 
      {
          CubitSimpleAttrib *tmp_attrib = new CubitSimpleAttrib( "SOURCE_FEATURE", "IMPRINT" );
          new ATTRIB_SNL_SIMPLE( new_edges[edge_count], tmp_attrib );
          delete tmp_attrib;
      }

      if (GeometryModifyTool::get_old_names())
     {
       Body* model_ent = dynamic_cast<Body*>(BodyPtr2->topology_entity());
       if (model_ent)
       {
         Body* new_model_ent = dynamic_cast<Body*>(newBodyPtr2->topology_entity());
         if (!new_model_ent)
           new_model_ent = GeometryQueryTool::instance()->make_Body(newBodyPtr2);
         model_ent->switch_entity_names(new_model_ent);
       }
     }

     if (keep_old == CUBIT_FALSE)
       AcisQueryEngine::instance()->delete_solid_model_entities(BodyPtr2);
     //actuate_mesh_intervals( newBodyPtr2 );
  }
  else
  {
     AcisQueryEngine::instance()->delete_ACIS_BODY(BODYPtr2Save, CUBIT_TRUE);
  }

    // All's well with the world :-)
  return CUBIT_SUCCESS;

}

BodySM *AcisModifyEngine::get_new_Body(  BodySM *old_Body,
                                       BODY *old_BODY,
                                       BODY *new_BODY,
                                       const bool keep_old,
                                       const bool topology_check,
                                       const bool delete_old) const
{
  if (topology_check)
  {
      // must do topology check between old & new bodies.  In the case
      // of a topology check, we know there's only one old and new body
    int num_vol_bef, num_face_bef, num_edge_bef, num_vertex_bef;
    CubitStatus status1 = AcisQueryEngine::instance()->number_ENTITIES( old_BODY, num_vol_bef,
                                           num_face_bef, num_edge_bef,
                                           num_vertex_bef );

    int num_vol_aft, num_face_aft, num_edge_aft, num_vertex_aft;
    CubitStatus status2 = AcisQueryEngine::instance()->number_ENTITIES( new_BODY, num_vol_aft,
                                           num_face_aft, num_edge_aft,
                                           num_vertex_aft );

    if( status1 == CUBIT_FAILURE  ||
        status2 == CUBIT_FAILURE )
    {
      PRINT_ERROR( "Unable to count number of entities in volume\n" );

      return 0;
    }

    if( num_vol_bef == num_vol_aft &&
        num_face_bef == num_face_aft &&
        num_edge_bef == num_edge_aft &&
        num_vertex_bef == num_vertex_aft)
    {
      AcisQueryEngine::instance()->delete_ACIS_BODY(new_BODY, CUBIT_TRUE);
      return 0;
    }
  }

    // convert input to list-based version and pass to overloaded version of
    // get_new_body
  DLIList<BodySM*> Body_list(1);
  if(old_Body != NULL)
     Body_list.append(old_Body);
  DLIList<BODY*> BODY_list(1);
  if(old_BODY != NULL)
     BODY_list.append(old_BODY);

  return get_new_Body(Body_list, BODY_list, new_BODY, keep_old, delete_old);
}

BodySM *AcisModifyEngine::get_new_Body(DLIList<BodySM*> &old_Body_list,
                                       DLIList<BODY*> &old_BODY_list,
                                       BODY *new_BODY,
                                       const bool keep_old,
                                       const bool delete_old) const
{
  DLIList<BODY*> new_BODY_list;
  new_BODY_list.append(new_BODY);
  DLIList<BodySM*> new_Body_list;
  CubitStatus status = get_new_Body(old_Body_list, old_BODY_list,
                                    new_BODY_list, new_Body_list,
                                    keep_old, delete_old);

  if (status == CUBIT_FAILURE ||
      new_Body_list.size() == 0) return NULL;

  assert(new_Body_list.size() == 1);

  return new_Body_list.get();
}

CubitStatus AcisModifyEngine::get_new_Body( BodySM *old_Body,
                                             BODY *old_BODY,
                                             DLIList<BODY*> &new_BODY_list,
                                             DLIList<BodySM*> &new_Body_list,
                                             const bool keep_old,
                                             const bool delete_old) const
{
  DLIList<BODY*> old_BODY_list;
  old_BODY_list.append(old_BODY);
  DLIList<BodySM*> old_Body_list;
  old_Body_list.append(old_Body);
  return get_new_Body(old_Body_list, old_BODY_list,
                      new_BODY_list, new_Body_list,
                      keep_old, delete_old);
}

CubitStatus AcisModifyEngine::get_new_Body(DLIList<BodySM*> &old_Body_list,
                                             DLIList<BODY*> &old_BODY_list,
                                             DLIList<BODY*> &new_BODY_list,
                                             DLIList<BodySM*> &new_Body_list,
                                             const bool keep_old,
                                             const bool delete_old) const
{
  assert(old_Body_list.size() == old_BODY_list.size());

       //&& old_Body_list.size() > 0); It's possible now to have no
       //old bodies, so don't stress over this

    // Only create a new Body if it changed
  BodySM *new_Body = NULL;

    // body changed; need to build a new vgi body around
    // old_BODY, either by reusing old_Body or creating new one

  if (GeometryModifyTool::instance()->get_new_ids() == CUBIT_FALSE &&
      keep_old == CUBIT_FALSE)
  {
      // if we're reusing old refentities, prepare by marking
      // deactivated the ones that don't get reused, then delete
      // those
    int num_deactivd = 0;
    int i;
    for (i = old_BODY_list.size(); i > 0; i--)
      num_deactivd += mark_owners_deactivated_flag(old_BODY_list.get_and_step(),
                                                   CUBIT_TRUE);

    int num_reactivd =  0;
    for (i = new_BODY_list.size(); i > 0; i--)
      num_reactivd =  mark_owners_deactivated_flag(new_BODY_list.get_and_step(),
                                                   CUBIT_FALSE);
    if (DEBUG_FLAG(99))
    {
      PRINT_DEBUG_99( "Number of entities deactivated = %d,"
                      " reactivated = %d\n", num_deactivd, num_reactivd);
      print_deleted_reused(old_BODY_list, new_BODY_list);
    }

      // ok, now unhook the old body from the vgi data (only done if
      // we aren't keeping old geometry)
    if (keep_old == CUBIT_FALSE)
    {
      for (i = old_BODY_list.size(); i > 0; i--)
        AcisQueryEngine::instance()->delete_ACIS_BODY(old_BODY_list.get_and_step(), CUBIT_FALSE);
    }

      // now cleanout the deactivated geometry; automatically unhooks
      // from the dag
    cleanout_deactivated_geometry();
  }
    // now delete the old body, either obeying or not obeying the
    // deactivated flags (which is opposite of the
    // create_new_refents flag)
  if (GeometryModifyTool::instance()->get_new_ids() && !keep_old && delete_old)
  {
    for (int i = old_Body_list.size(); i > 0; i--)
    {
      BodySM *old_Body = old_Body_list.get_and_step();
      AcisQueryEngine::instance()->delete_solid_model_entities(old_Body);
    }
  }

  int i;
  for (i = new_BODY_list.size(); i > 0; i--)
  {
    new_Body = AcisQueryEngine::instance()->populate_topology_bridges(new_BODY_list.get_and_step());
    if (new_Body != NULL)
      new_Body_list.append(new_Body);
    else
    {
      PRINT_ERROR("Problems creating new Body from new BODY.\n"
                  "       get_new_Body failed.\n") ;
      return CUBIT_FAILURE;
    }
  }

  if( keep_old && GeometryModifyTool::get_old_names() )
  {
      // switch names with the first body only
    old_Body_list.reset();
    BodySM* old_bodysm = old_Body_list.get();
    Body* body = dynamic_cast<Body*>(old_bodysm->topology_entity());
    if (body)
    {
      Body* new_body = dynamic_cast<Body*>(new_Body->topology_entity());
      if (!new_body)
        new_body = GeometryQueryTool::instance()->make_Body(new_Body);
      body->switch_entity_names(new_body);
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus AcisModifyEngine::get_new_Body( DLIList<TopologyBridge*> &old_entity_list,
                                            DLIList<BODY*> &new_BODY_list,
                                            DLIList<BodySM*> &new_Body_list,
                                            const bool keep_old,
                                            const bool delete_old) const
{
  assert(old_entity_list.size() > 0 || keep_old );

    // don't do a topology check in this version of get_new_Body

    // need to build a new body or bodies

  if (GeometryModifyTool::instance()->get_new_ids() == CUBIT_FALSE &&
      keep_old == CUBIT_FALSE)
  {
      // if we're reusing old refentities, prepare by marking
      // deactivated the ones that don't get reused, then delete
      // those
    int i, j, num_deactivd = 0, num_reactivd = 0;

    ENTITY_LIST old_ENTITY_list, tmp_ENTITY_list;
    old_entity_list.reset();
    for (i = old_entity_list.size(); i--; )
    {
      tmp_ENTITY_list.clear();
      TopologyBridge* geom_ptr = old_entity_list.get_and_step();
      AcisBridge* bridge_ptr = dynamic_cast<AcisBridge*>(geom_ptr);
      if (!bridge_ptr)
      {
        PRINT_ERROR("Non-ACIS geometry at %s:%d\n",__FILE__,__LINE__);
        continue;
      }
      ENTITY* ENTITY_ptr = bridge_ptr->ENTITY_ptr();
      if (!ENTITY_ptr)
        continue;

      old_ENTITY_list.add(ENTITY_ptr);


      CubitStatus stat = AcisQueryEngine::instance()->
        get_child_ENTITYs( ENTITY_ptr, tmp_ENTITY_list, true);
      if (!stat)
      {
        PRINT_ERROR("AcisQueryEngine::get_child_ENTITYs failed at %s:%d\n",__FILE__,__LINE__);
        continue;
      }

      for (j = 0; j < tmp_ENTITY_list.count(); j++)
      {
        old_ENTITY_list.add( tmp_ENTITY_list[j] );
      }
    }

    for (i = 0; i < old_ENTITY_list.count(); i++)
    {
         // mark the owners deactivated, but don't recurse, since we've got all
         // the lower order topology already
       num_deactivd += mark_owners_deactivated_flag(old_ENTITY_list[i],
                                                    CUBIT_TRUE,
                                                    CUBIT_FALSE);
    }


    for (i = new_BODY_list.size(); i > 0; i--)
      num_reactivd += mark_owners_deactivated_flag(new_BODY_list.get_and_step(),
                                                   CUBIT_FALSE);
    if (DEBUG_FLAG(99))
    {
      PRINT_DEBUG_99( "Number of entities deactivated = %d,"
                      " reactivated = %d\n", num_deactivd, num_reactivd);
    }

      // ok, now unhook the old body from the vgi data (only done if
      // we aren't keeping old geometry)
    if (keep_old == CUBIT_FALSE)
    {
      ENTITY_LIST to_be_deleted;
      for (i = 0; i < old_ENTITY_list.count(); i++)
      {
          // need to check to make sure the old entity isn't in the new body
        ENTITY *this_ENTITY = old_ENTITY_list[i];
        BODY *this_BODY = AcisQueryEngine::instance()->get_BODY_of_ENTITY(this_ENTITY);
	        //if the entity is not in the new body, add it to a list of
          //entities to be deleted
        if (!new_BODY_list.move_to(this_BODY))
          to_be_deleted.add(this_ENTITY);
      }
      //delete the ENTITYs in the list
      AcisQueryEngine::instance()->delete_ACIS_ENTITY(to_be_deleted);
    }

      // now cleanout the deactivated geometry; automatically unhooks
      // from the dag
    cleanout_deactivated_geometry();
  }


    // now delete the old body, either obeying or not obeying the
    // deactivated flags (which is opposite of the
    // create_new_refents flag)
  int i;
  if (GeometryModifyTool::instance()->get_new_ids() && !keep_old && delete_old)
  {
    for (i = old_entity_list.size(); i > 0; i--)
    {
      TopologyBridge *old_entity = old_entity_list.get_and_step();
      if (BodySM* bod_ptr = dynamic_cast<BodySM*>(old_entity))
        AcisQueryEngine::instance()->delete_solid_model_entities(bod_ptr);
      else
        AcisQueryEngine::instance()->delete_solid_model_entities(
          dynamic_cast<GeometryEntity*>(old_entity),true);
    }
    cleanout_deactivated_geometry();
  }

  for (i = new_BODY_list.size(); i > 0; i--)
  {
    BodySM *new_Body = AcisQueryEngine::instance()->populate_topology_bridges(new_BODY_list.get_and_step());
    if (new_Body != NULL)
      new_Body_list.append(new_Body);
    else
    {
      PRINT_ERROR("Problems creating new Body from new BODY.\n"
                  "       get_new_Body failed.\n") ;
      return CUBIT_FAILURE;
    }
  }

  return CUBIT_SUCCESS;
}

int AcisModifyEngine::mark_owners_deactivated_flag(ENTITY *this_ENTITY,
                                                     bool flag,
                                                     bool recurse) const

{
    // marks as deactivated all the owners of the input BODY,
    // along with all its children
    //
    // order is important here; entities must be put in the
    // deactivated list from lower order to higher order, that saves
    // work when removing the DAG nodes from the DAG
    // (not sure it's absolutely necessary though, since the
    // remove_from_DAG() function calls that function recursively for
    // lower order entities)

  int num_deactivated = 0;

    // if we're not recursing, set the skip flag to true
  CubitBoolean skip = (recurse ? CUBIT_FALSE : CUBIT_TRUE);
  ENTITY_LIST ENTITIES;
  outcome result;

    // vertices
  if (IS_ENTITY_TYPE( this_ENTITY, VERTEX )) skip = CUBIT_TRUE;
  if (!skip) {
    result = api_get_vertices( this_ENTITY, ENTITIES);
    num_deactivated += mark_owners_deactivated_flag( ENTITIES, flag );
    ENTITIES.clear();
  }

    // edges
  if (IS_ENTITY_TYPE( this_ENTITY, EDGE )) skip = CUBIT_TRUE;
  if (!skip) {
    result = api_get_edges( this_ENTITY, ENTITIES);
    num_deactivated += mark_owners_deactivated_flag( ENTITIES, flag );
    ENTITIES.clear();
  }

    // coedges
  if (IS_ENTITY_TYPE( this_ENTITY, COEDGE )) skip = CUBIT_TRUE;
  if (!skip) {
    result = api_get_coedges( this_ENTITY, ENTITIES);
    num_deactivated += mark_owners_deactivated_flag( ENTITIES, flag );
    ENTITIES.clear();
  }

    // loops
  if (IS_ENTITY_TYPE( this_ENTITY, LOOP )) skip = CUBIT_TRUE;
  if (!skip) {
    result = api_get_loops( this_ENTITY, ENTITIES);
    num_deactivated += mark_owners_deactivated_flag( ENTITIES, flag );
    ENTITIES.clear();
  }

    // faces
  if (IS_ENTITY_TYPE( this_ENTITY, FACE )) skip = CUBIT_TRUE;
  if (!skip) {
    result = api_get_faces( this_ENTITY, ENTITIES);
    num_deactivated += mark_owners_deactivated_flag( ENTITIES, flag );
    ENTITIES.clear();
  }

    // shells
  if (IS_ENTITY_TYPE( this_ENTITY, SHELL )) skip = CUBIT_TRUE;
  if (!skip) {
    result = api_get_shells( this_ENTITY, ENTITIES);
    num_deactivated += mark_owners_deactivated_flag( ENTITIES, flag );
    ENTITIES.clear();
  }

    // volumes
  if (IS_ENTITY_TYPE( this_ENTITY, LUMP )) skip = CUBIT_TRUE;
  if (!skip) {
    result = api_get_lumps( this_ENTITY, ENTITIES);
    num_deactivated += mark_owners_deactivated_flag( ENTITIES, flag );
    ENTITIES.clear();
  }

    // now, finally, do the entity itself
  ENTITIES.add(this_ENTITY);
  num_deactivated += mark_owners_deactivated_flag( ENTITIES, flag );

  return num_deactivated;
}

int AcisModifyEngine::mark_owners_deactivated_flag( ENTITY_LIST& ENTITIES,
                                                    bool flag ) const

{
  int num_deactivated = 0;

  for (int i = 0; i < ENTITIES.count(); i++)
  {
    AcisBridge *acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(ENTITIES[i]);
    if (acis_bridge != NULL)
    {
      if (flag == false)
      {
          // if bridge was deactivated, re-activate it and...
        if (reactivate_bridge(acis_bridge))
        {
              // ...and get the pointers between
              // ACIS and the TB right.
          acis_bridge->ENTITY_ptr(ENTITIES[i]);
          num_deactivated++;
        }
        else // bridge was already reactivated
        {
            // if we're reactivating an entity and it has already been reactivated,
            // that means another acis entity is using it - make this acis entity
            // get a new one
          ATTRIB_CUBIT_OWNER::remove_cubit_owner(ENTITIES[i]);
        }
      }

        // if flag is true and acis_bridge is NOT deactivated,
        // deactivate it.
      else if (deactivate_bridge(acis_bridge))
      {
        num_deactivated++;
      }
    }
  }

  return num_deactivated;
}

void AcisModifyEngine::print_deleted_reused(DLIList<BODY*> &old_BODY_list,
                                              DLIList<BODY*> &new_BODY_list) const
{

    // print out what's saved and what's not
  DLIList<TopologyBridge*> deact_list, react_list;
  int i;
  for (i = old_BODY_list.size(); i > 0; i--)
    AcisQueryEngine::instance()->get_all_cubit_owners(old_BODY_list.get_and_step(), deact_list);
  for (i = new_BODY_list.size(); i > 0; i--)
    AcisQueryEngine::instance()->get_all_cubit_owners(new_BODY_list.get_and_step(), react_list);
  for (i = react_list.size(); i > 0; i--)
    deact_list.remove(react_list.get_and_step());
  PRINT_DEBUG_99( "Deleted entities:\n");
  for (i = deact_list.size(); i > 0; i--) {
    TopologyBridge* bridge = deact_list.get_and_step();
    TopologyEntity *entity = bridge->topology_entity();
    BasicTopologyEntity *bte = CAST_TO(entity, BasicTopologyEntity);
    GroupingEntity *grpe = CAST_TO(entity, GroupingEntity);
    SenseEntity *se = CAST_TO(entity, SenseEntity);
    if (bte != NULL) {
      PRINT_DEBUG_99( "%s %d\n",
                      bte->class_name(),
                      bte->id());
    }
    else if (se != NULL) {
      PRINT_DEBUG_99( "%s (lower BTE: %d)\n",
                      entity->class_name(),
                      (se->get_basic_topology_entity_ptr() ?
                      se->get_basic_topology_entity_ptr()->id() : 0));
    }
    else if (grpe != NULL) {
      PRINT_DEBUG_99( "%s (upper BTE: %d)\n",
                      entity->class_name(),
                      (grpe->get_basic_topology_entity_ptr() ?
                       grpe->get_basic_topology_entity_ptr()->id() : 0));
    }
  }
  PRINT_DEBUG_99( "ReUsed entities:\n");
  for (i = react_list.size(); i > 0; i--) {
    TopologyBridge* bridge = react_list.get_and_step();
    TopologyEntity *entity = bridge->topology_entity();
    BasicTopologyEntity *bte = CAST_TO(entity, BasicTopologyEntity);
    GroupingEntity *grpe = CAST_TO(entity, GroupingEntity);
    SenseEntity *se = CAST_TO(entity, SenseEntity);
    if (bte != NULL) {
      PRINT_DEBUG_99( "%s %d\n",
                      bte->class_name(),
                      bte->id());
    }
    else if (se != NULL) {
      PRINT_DEBUG_99( "%s (lower BTE: %d)\n",
                      entity->class_name(),
                      (se->get_basic_topology_entity_ptr() ?
                      se->get_basic_topology_entity_ptr()->id() : 0));
    }
    else if (grpe != NULL) {
      PRINT_DEBUG_99( "%s (upper BTE: %d)\n",
                      entity->class_name(),
                      (grpe->get_basic_topology_entity_ptr() ?
                       grpe->get_basic_topology_entity_ptr()->id() : 0));
    }
  }
}

//-------------------------------------------------------------------------
// Purpose       : Imprint multiple bodies at once
//
// Special Notes :
//
// Creator       : Tim Tautges
//
// Creation Date : 12/8/98
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::imprint(DLIList<BodySM*> &from_body_list,
                                      DLIList<BodySM*> &new_from_body_list,
                                      bool keep_old,
                                      DLIList<TopologyBridge*> *new_tbs,
                                      DLIList<TopologyBridge*> *att_tbs) const
{

  DLIList<BODY*> from_BODY_list;
  DLIList<BODY*> from_BODY_list_copy;
  CubitStatus success = CUBIT_SUCCESS;

  // get acis bodies, and copy them
  from_body_list.reset();
  int i, j;

    // pass in keep_old to the delete_owner_attrib flag; if we're not keeping
    // old bodies, we want to keep the owner attrib, so we can pick up entities
    // that didn't change
  bool delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old);

  const int num_bodies = from_body_list.size();
  from_body_list.reset();
  for (i = 0; i < num_bodies; i++) {
    BODY *BODYPtr1 = AcisQueryEngine::get_BODY(from_body_list.next(i));
    from_BODY_list.append(BODYPtr1);
    //BODY *BODYPtr1Save = this->copy_BODY(BODYPtr1, delete_attribs);
    //from_BODY_list_copy.append(BODYPtr1Save);
    from_BODY_list_copy.append(NULL);
  }

  if (from_BODY_list.is_in_list(NULL))
    return CUBIT_FAILURE;

  // set the free_edges acis option
  if ( GeometryModifyTool::get_all_edges_imprint() )
  {
    api_set_int_option("all_free_edges", TRUE );
  }

      // total number of imprints to be done
  int total_imprints = (from_body_list.size() *
                        (from_body_list.size()-1))/2;

  if( num_bodies > 2 )
  {
     char message[128];
     sprintf(message, "Imprinting %d ACIS Bodies", from_BODY_list.size() );
     AppUtil::instance()->progress_tool()->start(0, total_imprints, message);
  }

    // imprint each body with each one later in the list
  from_BODY_list.reset();
  from_BODY_list_copy.reset();
  for (i = 0; success && i < num_bodies - 1; i++) {
    //BODY *BODY_Ptr1 = from_BODY_list_copy.next(i);
    BodySM *Body_Ptr1 = from_body_list.next(i);
    CubitBox box1 = AcisQueryEngine::instance()->bounding_box(Body_Ptr1);

    for (j = i + 1; j < num_bodies; j++) {

      if (CubitMessage::instance()->Interrupt())
      {
        success = CUBIT_FAILURE;
        break;
      }

      //BODY *BODY_Ptr2 = from_BODY_list_copy.next(j);
      BodySM *Body_Ptr2 = from_body_list.next(j);
      CubitBox box2 = AcisQueryEngine::instance()->bounding_box(Body_Ptr2);

        // first, check bounding SPAbox; if they don't intersect, imprint will fail

      if (box1.overlap(SPAresabs, box2)) 
      {
        if (DEBUG_FLAG(95))
        {
          TopologyEntity *te1 = Body_Ptr1->topology_entity();
          TopologyEntity *te2 = Body_Ptr2->topology_entity();
          RefEntity *ref1 = CAST_TO(te1, RefEntity );
          RefEntity *ref2 = CAST_TO(te2, RefEntity );
          PRINT_INFO("Imprinting bodies %d and %d\n", ref1->id(), ref2->id() );
        }

        BODY* BODY_Ptr1 = from_BODY_list_copy.next(i);
        BODY* BODY_Ptr2 = from_BODY_list_copy.next(j);
        if (!BODY_Ptr1)
          BODY_Ptr1 = from_BODY_list.next(i);
        if (!BODY_Ptr2)
          BODY_Ptr2 = from_BODY_list.next(j);

        BODY_Ptr1 = copy_BODY(BODY_Ptr1, delete_attribs);
        BODY_Ptr2 = copy_BODY(BODY_Ptr2, delete_attribs);
        // the bodies overlap; proceed with the imprint
        outcome result = api_imprint(BODY_Ptr1, BODY_Ptr2);

        AcisModifyEngine::instance()->cleanup_slivers( BODY_Ptr1 );
        AcisModifyEngine::instance()->cleanup_slivers( BODY_Ptr2 );

        if (result.ok())
        {
          from_BODY_list_copy.step(i);
          if (from_BODY_list_copy.get())
            api_delent(from_BODY_list_copy.get());
          from_BODY_list_copy.change_to(BODY_Ptr1);
          from_BODY_list_copy.reset();

          from_BODY_list_copy.step(j);
          if (from_BODY_list_copy.get())
            api_delent(from_BODY_list_copy.get());
          from_BODY_list_copy.change_to(BODY_Ptr2);
          from_BODY_list_copy.reset();
        }
        else
        {
          if (DEBUG_FLAG(95)) {
            AcisQueryEngine::instance()->ACIS_API_error(result, "imprint Bodies");
          }

          if (BODY_Ptr1)
            api_delent(BODY_Ptr1);
          if (BODY_Ptr2)
            api_delent(BODY_Ptr2);
        }
      } // if overlap

      if( num_bodies > 2 )
        AppUtil::instance()->progress_tool()->step();
    }
  }

  if( num_bodies > 2 )
      AppUtil::instance()->progress_tool()->end();

  if( CubitMessage::instance()->Interrupt() )
    PRINT_INFO("Imprint aborted.\n");

  // new_ENTITIES will hold new entities that were created by 
  // the imprint operation.  att_ENTITIES will hold entities
  // from the original bodies that had composite attributes
  // on them.  At the end we will compare the new entities
  // with the ones that had attributes and if we find matches
  // we will put a composite attribute on the new entity.
  // This is how we are handling the imprinting of bodies
  // with composites.
  DLIList<ENTITY*> new_ENTITIES;
  DLIList<ENTITY*> att_ENTITIES;

  // Create/update ACIS Bridges.
  for ( i = 0; i < num_bodies; i++ )
  {
    BodySM* from_body = from_body_list.next(i);
    BODY*   from_BODY = from_BODY_list.next(i);
    BODY* new_BODY = from_BODY_list_copy.next(i);
    if (new_BODY)
    {
      // Get the new entities (just edges for now)
      // created by the imprint operation.
      DLIList<ENTITY*> cur_new_ENTITIES;
      if(new_tbs)
        get_new_ENTITIES(new_BODY, cur_new_ENTITIES);

      DLIList<EDGE*> new_edges = AcisModifyEngine::instance()->find_new_EDGES(new_BODY);

      BodySM* new_body = get_new_Body( from_body, from_BODY, new_BODY, keep_old, true );
      if (new_body)
      {
        // Add the new entities to the list we are
        // accumulating.
        if(new_tbs)
          new_ENTITIES += cur_new_ENTITIES;

        // Get entities with composite attributes.
        if(att_tbs)
          get_att_ENTITIES(new_BODY, att_ENTITIES, "COMPOSITE_GEOM");

        // Add a imprint feature to the topo edges
        for (int edge_count = new_edges.size(); edge_count--; ) 
        {
            CubitSimpleAttrib *tmp_attrib = new CubitSimpleAttrib( "SOURCE_FEATURE", "IMPRINT" );
            new ATTRIB_SNL_SIMPLE( new_edges[edge_count], tmp_attrib );
            delete tmp_attrib;
        }
        new_from_body_list.append(new_body);
      }
      else
      {
        // Get entities with composite attributes.
        if(att_tbs)
          get_att_ENTITIES(from_BODY, att_ENTITIES, "COMPOSITE_GEOM");
      }
    }
  }

  // Convert the new_ENTITIES/att_ENTITIES lists into
  // topology bridge lists.
  if(new_tbs)
  {
    for(i=new_ENTITIES.size(); i--;)
    {
      ENTITY *cur_ENT = new_ENTITIES.get_and_step();
      AcisBridge *acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(cur_ENT);
      if(dynamic_cast<TopologyBridge*>(acis_bridge))
        new_tbs->append_unique(dynamic_cast<TopologyBridge*>(acis_bridge));
    }
  }
  if(att_tbs)
  {
    for(i=att_ENTITIES.size(); i--;)
    {
      ENTITY *cur_ENT = att_ENTITIES.get_and_step();
      AcisBridge *acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(cur_ENT);
      if(dynamic_cast<TopologyBridge*>(acis_bridge))
        att_tbs->append_unique(dynamic_cast<TopologyBridge*>(acis_bridge));
    }
  }


  PRINT_INFO("Group imprint finished.\n");

  // reset the free_edges acis option
  api_set_int_option("all_free_edges", FALSE );

  return success;
}

//-------------------------------------------------------------------------
// Purpose       : Projects an EDGE to a FACE
//
// Special Notes :
//
//
// Creator       : Eric Nielsen
//
// Creation Date : 06/14/00
//-------------------------------------------------------------------------
EDGE* AcisModifyEngine::project_EDGE(EDGE* EDGE_in_ptr, FACE* FACE_ptr, bool print_error ) const
{
  curve* in_crv=NULL;
  surface *in_srf=NULL;
  curve * out_crv = NULL;
  EDGE * new_EDGE = NULL;

  if(FACE_ptr && FACE_ptr->geometry())
  {
    SPAtransf t= get_owner_transf(FACE_ptr);
    in_srf = FACE_ptr->geometry()->trans_surface(t);
  }

  if(EDGE_in_ptr && EDGE_in_ptr->geometry())
  {
    SPAtransf t= get_owner_transf(EDGE_in_ptr);
    in_crv = EDGE_in_ptr->geometry()->trans_curve(t);
  }

  if(in_crv && in_srf)
  {
    // create the curve
    SPAinterval in_range = EDGE_in_ptr->param_range();

    if( EDGE_in_ptr->sense() )
      in_range.negate();

    outcome result = api_project_curve_to_surface(*in_crv, in_range, *in_srf, out_crv);

    if (!result.ok() || out_crv == NULL) //out_crv check added 6/18/04 by aga@cat
    //api CAN return a null curve in some failures- if so we'll end the fcn here and return
    {
      if( print_error )
      {
        AcisQueryEngine::instance()->ACIS_API_error (result);
        PRINT_ERROR( "api_project_curve_to_surface failed ");
      }
      return (EDGE *)NULL;
    }

    // create the CURVE
    CURVE *the_curve=make_curve(*out_crv);

    // get ends of the edge
    SPAposition pt0=out_crv->eval_position(in_range.start_pt());
    SPAposition pt1=out_crv->eval_position(in_range.end_pt());

    ACIS_DELETE out_crv;

    // create the vertices
    VERTEX *svert= ACIS_NEW VERTEX( ACIS_NEW APOINT(pt0));

    // -- check for a closed edge
    VERTEX *evert= NULL;
    if ( pt0 != pt1 )
      evert = ACIS_NEW VERTEX( ACIS_NEW APOINT(pt1));
    else
      evert = svert;

    // make the edge
    new_EDGE= ACIS_NEW EDGE(svert,evert,the_curve,FORWARD,EDGE_cvty_unknown,in_range);
  }

  if(in_crv)
    ACIS_DELETE in_crv;
  if(in_srf)
    ACIS_DELETE in_srf;

  if(new_EDGE)
  {
    ENTITY *edge_entity = (ENTITY*)new_EDGE;
    api_simplify_entity( edge_entity );
    return new_EDGE;
  }
  else
    return (EDGE *)NULL;
}

//-------------------------------------------------------------------------
// Purpose       : Projects a list of RefEdges on to a list of RefFaces
//                 returning a list of new RefEdges
//
// Special Notes :
//
//
// Creator       : Eric Nielsen
//
// Creation Date : 06/14/00
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::project_edges( DLIList<Surface*> &ref_face_list,
                                             DLIList<Curve*> &ref_edge_list_in,
                                             DLIList<Curve*> &ref_edge_list_new,
                                             bool print_error ) const
{
  Curve* ref_edge = NULL;
  Surface* ref_face = NULL;
  EDGE *new_EDGE_ptr = NULL;

  for(int j=0; j< ref_face_list.size(); j++)
  {
    ref_face = ref_face_list.get_and_step();

    FACE *FACE_ptr = AcisQueryEngine::get_FACE( ref_face );

    for(int i=0; i< ref_edge_list_in.size(); i++)
    {
      ref_edge = ref_edge_list_in.get_and_step();

      EDGE *EDGE_ptr = AcisQueryEngine::get_EDGE(ref_edge);

      new_EDGE_ptr = project_EDGE(EDGE_ptr,FACE_ptr, print_error );

      if (!new_EDGE_ptr)
      {
        if( print_error )
          PRINT_ERROR( "Unable to project curve on surface\n" );
        continue;
      }
      else{
        Curve *curve_ptr = AcisQueryEngine::instance()->populate_topology_bridges(new_EDGE_ptr);
        ref_edge_list_new.append(curve_ptr);
      }

    }
  }

  if (ref_edge_list_new.size() == ref_edge_list_in.size())
    return CUBIT_SUCCESS;
  else
    return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : Projects a list of RefEdges on to a list of RefFaces
//                 and imprint the faces with the new RefEdges
//
// Special Notes :
//
//
// Creator       : Eric Nielsen
//
// Creation Date : 06/14/00
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::imprint_projected_edges(DLIList<Surface*> &ref_face_list,
                                                      DLIList<Curve*> &ref_edge_list_in,
                                                      DLIList<BodySM*>& new_body_list,
                                                      bool keep_old_body,
                                                      bool keep_free_edges) const
{
  CubitStatus status;

  DLIList<Curve*> ref_edge_list_new;

  project_edges(ref_face_list, ref_edge_list_in, ref_edge_list_new);

  // imprint Surface with curves
  status = imprint(ref_face_list, ref_edge_list_new, new_body_list, keep_old_body );

  if (keep_free_edges)
    return  status;

  PRINT_INFO( "Removing projected curves \n");
  for(int i=0; i< ref_edge_list_new.size();i++){
    // Now delete this RefEdge and its underlying solid model entities

    Curve* ref_edge = ref_edge_list_new.get_and_step();
    status = AcisQueryEngine::instance()->
      delete_solid_model_entities( ref_edge );
    if (status == CUBIT_FAILURE)
    {
      PRINT_ERROR("In GeometryQueryTool::delete_geometry\n"
        "       Could not delete RefEdge.\n"
        "       The Model database is likely corrupted "
        "due to\n       this unsuccessful deletion.\n" );
    }
  }
  return status;
}

//-------------------------------------------------------------------------
// Purpose       : Projects a list of RefEdges on to a list of RefFaces
//                 and imprint the bodies with the new RefEdges
//
// Special Notes :
//
//
// Creator       : Eric Nielsen
//
// Creation Date : 06/14/00
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::imprint_projected_edges(DLIList<Surface*> &ref_face_list,
                                                      DLIList<BodySM*> &body_list,
                                                      DLIList<Curve*> &ref_edge_list_in,
                                                      DLIList<BodySM*>& new_body_list,
                                                      bool keep_old_body,
                                                      bool keep_free_edges) const
{
  CubitStatus status;

  DLIList<Curve*> ref_edge_list_new;

  project_edges(ref_face_list, ref_edge_list_in, ref_edge_list_new);

  // imprint bodies with curves
  status = imprint(body_list,ref_edge_list_new,new_body_list, keep_old_body );

  if (keep_free_edges)
    return  status;

  PRINT_INFO( "Removing projected curves \n");
  for(int i=0; i< ref_edge_list_new.size();i++){
    // Now delete this RefEdge and its underlying solid model entities

    Curve* ref_edge = ref_edge_list_new.get_and_step();
    status = AcisQueryEngine::instance()->
      delete_solid_model_entities( ref_edge );
    if (status == CUBIT_FAILURE)
    {
      PRINT_ERROR("In GeometryModifyTool::delete_geometry\n"
        "       Could not delete RefEdge.\n"
        "       The Model database is likely corrupted "
        "due to\n       this unsuccessful deletion.\n" );
    }
  }
  return status;
}

//-------------------------------------------------------------------------
// Purpose       : Boolean Operation of body with list of bodies.
//
// Special Notes : INTERSECTION
//
// Creator       : David White
//
// Creation Date : 11/24/97
//------------------------------------------------------------------------
CubitStatus AcisModifyEngine::intersect( BodySM* tool_body,
                                           DLIList<BodySM*> &from_bodies,
                                           DLIList<BodySM*> &new_from_bodies,
                                           bool keep_old) const
{
    // pass in keep_old to the delete_owner_attrib flag; if we're not keeping
    // old bodies, we want to keep the owner attrib, so we can pick up entities
    // that didn't change
  bool delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old);


  int ii;
  CubitBox tool_box = AcisQueryEngine::instance()->bounding_box(tool_body);

  DLIList<BODY*> from_BODY_list;
  DLIList<BODY*> from_BODY_list_copy;

  // get acis bodies, and copy them
  from_bodies.reset();

  for (ii = from_bodies.size(); ii > 0; ii--) {
    BODY *BODYPtr1 = AcisQueryEngine::get_BODY(from_bodies.get_and_step());
    from_BODY_list.append(BODYPtr1);
    BODYPtr1 = this->copy_BODY(BODYPtr1, delete_attribs);
    from_BODY_list_copy.append(BODYPtr1);
  }

  if (from_BODY_list.size() != from_bodies.size())
    return CUBIT_FAILURE;

    // check if the boolean operation is regularized or nonregularized
  BOOL_TYPE bool_type = INTERSECTION;
  CubitBoolean boolean_regularize = GeometryModifyTool::instance()->boolean_regularize();

  if (boolean_regularize == FALSE)
  {
    bool_type = NONREG_INTERSECTION;
  }

  // now, intersect the tool with the list of bodies
  from_BODY_list.reset();
  from_BODY_list_copy.reset();
  from_bodies.reset();

  CubitStatus int_stat = CUBIT_SUCCESS;

    // intersect the tool body with each body in the list
  for (ii = 1; ii <= from_BODY_list.size(); ii++) {
    BODY *from_BODY = from_BODY_list.get();
    BODY *from_BODY_copy = from_BODY_list_copy.get();
    BodySM *from_Body = from_bodies.get();
    CubitBox box1 = AcisQueryEngine::instance()->bounding_box(from_Body);

      // first, check bounding SPAbox; if they don't intersect, don't do the intersection

    if (from_Body == tool_body) {
      PRINT_ERROR("Can't intersect volume with itself.\n");
    }

    else if (tool_box.overlap(SPAresabs, box1)) {

        // the bodies overlap; proceed with the intersect
      BODY* tool_BODY = AcisQueryEngine::get_BODY(tool_body);
      BODY *tool_BODY_copy = this->copy_BODY(tool_BODY);


        // Intersect body1 with body2.
        //outcome result = api_intersect( tool_BODY_copy, from_BODY_copy );
       outcome result = api_boolean(tool_BODY_copy, from_BODY_copy, bool_type );

        // We may or may not get an error if the resulting BODY is
        // empty
        // Also check the result.
      if (from_BODY_copy == NULL ||
          (from_BODY_copy->lump() == NULL && from_BODY_copy->wire() == NULL) ||
          !result.ok())
      {
        if (!result.ok())
        {
          PRINT_ERROR("Intersection operation failed.\n");
          AcisQueryEngine::instance()->ACIS_API_error(result, "intersect Bodies");
          int_stat = CUBIT_FAILURE;

          //I might need to take this out.  Not sure. CDE 10/17/2003
          if (tool_BODY_copy != NULL)
          {
            AcisQueryEngine::instance()->delete_ACIS_BODY(tool_BODY_copy, CUBIT_TRUE);
          }
        }
        else
          PRINT_ERROR("Intersection operation failed.\n"
                      "       Empty volume created and will be deleted.\n");
        if (from_BODY_copy != NULL)
        {
          AcisQueryEngine::instance()->delete_ACIS_BODY(from_BODY_copy, CUBIT_TRUE);
        }

          // we had an error, so put the original body into the copy list
        from_BODY_list_copy.change_to(from_BODY);
      }

      if ( bool_type == NONREG_INTERSECTION)
      {
		    outcome result1;
		    // ENTITY_LIST FACES;
            DLIList<FACE*>  FACE_list ;
		    FACE *this_FACE;
            AcisQueryEngine::instance()->get_FACEs(from_BODY_copy,FACE_list);
		     // loop through all faces
		    for( int i = 0; i< FACE_list.size(); i++)
		    {
			    this_FACE = FACE_list.get_and_step();
			    assert( this_FACE != NULL );
			    // Make sure this is a DOUBLE_SIDED FACE
			    if (this_FACE->sides() == DOUBLE_SIDED)
			    {
				    if (this_FACE->cont()==BOTH_INSIDE)
				    {
					    // found internal faces
					    // Now we have the FACE - unhook it  from the BODY.  Keep track of new
					    // BODIES that are created as this is done.
					     PRINT_INFO( " Unhooking and deleting each internal surface...\n" );
					     BODY *new_BODY_ptr;
					     result1 = api_unhook_face( this_FACE, new_BODY_ptr );
					     if( !result1.ok() )
					     {
						     AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result1 );
						     PRINT_ERROR( " Face unhooking during rebuild of volume didn't work\n" );
						     return CUBIT_FAILURE;
					     }
					     AcisQueryEngine::instance()->delete_ACIS_BODY (new_BODY_ptr);
				    }
			    }
		    }

		 // heal the leftover body
//
//		   PRINT_INFO(" Healing the leftover volume...\n");
//		   if( AcisHealerTool::instance()->init_BODY_for_healing( master ) == CUBIT_SUCCESS )
//		   {
//			   int percent_before, percent_after, number_splines_simplified;
//			   if( AcisHealerTool::instance()->heal_BODY( master, percent_before,
//				   percent_after, number_splines_simplified ) == CUBIT_FAILURE )
//				   PRINT_ERROR( "Error healing the combined volume\n" );
//			   else
//				   PRINT_INFO( "Successfully healed the combined volume.\n" );
//			   AcisHealerTool::instance()->end_BODY_for_healing( master );
//		   }
	    }

      // done with this j iteration; write out count, if necessary
//      if (from_bodies.size() > 1) {
//        int frac_done = (10 * ii) / from_bodies.size();
//        if ((10 - frac_done) < fraction_remaining) {
//          PRINT_INFO("%d ", fraction_remaining);
//          fraction_remaining--;
//        }
//      }
    }

    else {
      AcisQueryEngine::instance()->delete_ACIS_BODY(from_BODY_copy, CUBIT_TRUE);
      from_BODY_list_copy.change_to(from_BODY);
    }

      // done with iteration over ii; step the lists
    from_BODY_list_copy.step();
    from_BODY_list.step();
    from_bodies.step();
  }

//  PRINT_INFO("\n");

    // ok, we're done; construct new Body's for the new BODY's
  from_BODY_list_copy.reset();
  from_BODY_list.reset();
  from_bodies.reset();

  for (ii = from_BODY_list_copy.size(); ii > 0; ii--) {

    BODY *old_BODY = from_BODY_list.get();
    BODY *new_BODY = from_BODY_list_copy.get();

    BodySM *new_body = NULL;
    if (old_BODY != new_BODY)
      new_body = get_new_Body(from_bodies.get(), old_BODY, new_BODY, keep_old);
    else if (!keep_old && int_stat == CUBIT_SUCCESS) {
      //GeometryQueryTool::instance()->delete_Body(from_bodies.get());
      AcisQueryEngine::instance()->delete_solid_model_entities(from_bodies.get());
      from_bodies.change_to(NULL);
    }

    if (new_body) {
      new_from_bodies.append(new_body);
      from_bodies.change_to(NULL);
    }

      // now step all the lists
    from_BODY_list.step();
    from_BODY_list_copy.step();
    from_bodies.step();
  }

  from_bodies.remove_all_with_value(NULL);
  from_bodies.reset();
  new_from_bodies.reset();

  if (int_stat == CUBIT_SUCCESS && !keep_old )
    //GeometryQueryTool::instance()->delete_Body(tool_body);
    AcisQueryEngine::instance()->delete_solid_model_entities(tool_body);
/*
    // write the message about bodies kept and not kept
//NOTE TO JS:
  PRINT_INFO("Intersect finished; \n"
             "Old bodies retained: ");
  if (from_bodies.size() > 0) {
    for (ii = from_bodies.size(); ii > 0; ii--)
      PRINT_INFO("%d ", from_bodies.get_and_step()->id());
  }
  else {
    PRINT_INFO(" (none)");
  }

  PRINT_INFO("\n");

//NOTE TO JS:
  PRINT_INFO("New bodies created: ");
  if (new_from_bodies.size() > 0) {
    for (ii = new_from_bodies.size(); ii > 0; ii--)
      PRINT_INFO("%d ", new_from_bodies.get_and_step()->id());
  }
  else {
    PRINT_INFO(" (none)");
  }
  PRINT_INFO("\n");
*/
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Perform a rotational sweep
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/31/96
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::sweep_rotational(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  const CubitVector& point,
  const CubitVector& sweep_axis,
  double angle,
  int steps,
  double draft_angle,
  int draft_type,
  bool switchside,
  bool make_solid,
  bool rigid) const
{

  DLIList<FACE*> FACEs_to_be_swept_list;
  DLIList<EDGE*> EDGEs_to_be_swept_list;
  DLIList<Surface*> ref_faces_list;
  DLIList<Curve*> ref_edges_list;

    // First get the ACIS FACEs associated with the input RefFaces
  //DLIList<FACE*> temp_faces;
  int i;
  for (i = ref_ent_list.size(); i > 0; i--)
  {
    GeometryEntity *ref_ent = ref_ent_list.get_and_step();
    Surface *ref_face = CAST_TO(ref_ent, Surface);
    Curve* ref_edge = CAST_TO(ref_ent, Curve);

    if(ref_face != NULL)
    {
      FACE* FACE_ptr = AcisQueryEngine::get_FACE(ref_face);
      if (!FACE_ptr)
        return CUBIT_FAILURE;
      //SurfaceACIS::get_FACEs_of_RefFace( ref_face, temp_faces);
      //if (temp_faces.size() > 1) {
      //  PRINT_ERROR("Surface %d cannot be swept as the number of ACIS FACEs associated \n"
      //              "       with it is not 1.\n",
      //              ref_face->id());
      //  return CUBIT_FAILURE;
      //}
      //
      //if (temp_faces.size() > 0)
      //   FACEs_to_be_swept_list += temp_faces;
      //temp_faces.clean_out();
      FACEs_to_be_swept_list.append(FACE_ptr);
      ref_faces_list.append(ref_face);
    }
    else if (ref_edge != NULL)
    {
      //EDGEs_to_be_swept_list.append(CurveACIS::get_first_EDGE(ref_edge));
      EDGEs_to_be_swept_list.append(AcisQueryEngine::get_EDGE(ref_edge));
      ref_edges_list.append(ref_edge);
    }
  }

  ref_edges_list.reset();
  ref_faces_list.reset();

      // Next, call ACIS to do the real work for each of the FACEs in the list!!
  FACEs_to_be_swept_list.reset();
  EDGEs_to_be_swept_list.reset();
  FACE* FACE_ptr = NULL;
  FACE* copied_FACE_ptr = NULL;
  EDGE* copied_EDGE_ptr = NULL;
  Surface* refface_ptr = NULL;
  BodySM* body_ptr = NULL;
  BODY* BODY_ptr = NULL;
  BODY* copy_of_original_BODY_ptr = NULL;
  bool primary_side = true;
  bool skip_loop = false;
  CubitStatus sweep_status = CUBIT_FAILURE;
  for (i = 0; i < FACEs_to_be_swept_list.size(); i++)
  {
    skip_loop = CUBIT_FALSE;
    copied_FACE_ptr = NULL;

      // Get the RefFace we're going to sweep
    GeometryEntity* ref_ent = ref_ent_list.get_and_step();
    refface_ptr = CAST_TO(ref_ent, Surface);
    if(make_solid)
    {
      PRINT_WARNING("The make_solid option is ignored for surface\n");
    }

      // Get the FACE we're going to sweep
    FACE_ptr = FACEs_to_be_swept_list.get_and_step();

      // If the FACE is DOUBLE-SIDED; deal with the switchside option
    if (FACE_ptr->sides() && switchside)
    {
      primary_side = CUBIT_FALSE;
    }

      // Make sure the FACE has a plane surface as it's underlying geometry
    if (FACE_ptr->geometry()->equation().type() != plane_type)
    {
      PRINT_ERROR("One of the surfaces being swept is not planar.\n"
                  "       Skipping to the next surface.\n");
      skip_loop = CUBIT_TRUE;
    }

      // If this FACE is planar, continue
    if( !skip_loop )
    {
        // Get the associated Body and it's ACIS BODY
      body_ptr = AcisQueryEngine::instance()->get_body_sm_of_ENTITY(FACE_ptr);
      if (body_ptr != NULL)
      {
        BODY_ptr = AcisQueryEngine::get_BODY(body_ptr);

          // Make sure we have a valid ACIS BODY associated with this Body
          // and FACE
        if (BODY_ptr == NULL)
        {
          PRINT_ERROR("No valid volume associated with the RefFace.\n");
          skip_loop = CUBIT_TRUE;
        }
      }

      else
      {
        PRINT_WARNING("No volume associated with this surface.\n"
                      "         Sweeping will continue with the ACIS BODY.\n");
        BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY(FACE_ptr);
        if (BODY_ptr == NULL)
        {
          PRINT_ERROR("No valid volume associated with the RefFace.\n");
          skip_loop = CUBIT_TRUE;
        }
      }

        // If there were no problems finding an associated ACIS BODY, proceed
      if( !skip_loop )
      {
          // Make a copy of the original BODY before doing the sweep.
          // We will perform the sweep on the copy, not on the original body.
        copy_of_original_BODY_ptr = this->copy_BODY(BODY_ptr,
                                                    CUBIT_FALSE);

          // Find the FACE in the copied BODY that corresponds to the FACE
          // being swept in the original BODY.  Note that when we made a copy
          // of the original BODY, all the "parent" attributes were copied
          // over to the new BODY.
        if (copy_of_original_BODY_ptr)
        {
          copied_FACE_ptr = (FACE *) (AcisQueryEngine::instance()->
            get_ENTITY_of_entity( copy_of_original_BODY_ptr, refface_ptr) );
        }

        if (copied_FACE_ptr)
        {
            // Remove all the CUBIT_OWNER attributes from the copied BODY
            // as we're going to sweep it.
          if (GeometryModifyTool::instance()->get_new_ids() == CUBIT_TRUE)
            AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copy_of_original_BODY_ptr);

            // Rotate the FACE about the axis, sweep_vector.
          int reverse_failed = CUBIT_FALSE;
          bool volume_is_negative = false;
          if ( !this->sweep_FACE_about_axis( copied_FACE_ptr,
                                             sweep_axis,
                                             point,
                                             angle,
                                             volume_is_negative,
                                             primary_side,
                                             steps,
                                             draft_angle,
                                             draft_type,
                                             rigid) )
          {
              // Volume was negative, so we're going to try to reverse the
              // BODY
            if (volume_is_negative)
            {
                // Delete the BODY that was just created (the one with the
                // negative volume).
                //AcisQueryEngine::instance()->delete_ACIS_BODY(copy_of_original_BODY_ptr);

                // Make another copy of the original BODY before attempting
                // the sweep again. We will perform the sweep on this copy,
                // not on the original body.
              copy_of_original_BODY_ptr = this->copy_BODY(BODY_ptr, CUBIT_FALSE);

                // Find the FACE in the copied BODY that corresponds to the FACE
                // being swept in the original BODY.  Note that when we made a copy
                // of the original BODY, all the "parent" attributes were copied
                // over to the new BODY.
              if (copy_of_original_BODY_ptr)
              {
                copied_FACE_ptr = (FACE *) (AcisQueryEngine::instance()->
                  get_ENTITY_of_entity( copy_of_original_BODY_ptr, refface_ptr ));
              }

              if (copied_FACE_ptr)
              {
                  // Remove all the CUBIT_OWNER attributes from the copied BODY
                  // as we're going to sweep it.
                if (GeometryModifyTool::instance()->get_new_ids() == CUBIT_TRUE)
                  AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copy_of_original_BODY_ptr);

                  // Find out how many FACEs exist in the BODY that owns
                  // copied_FACE_ptr
                DLIList<FACE*> FACE_list;
                if (AcisQueryEngine::instance()->get_FACEs( (ENTITY *)copy_of_original_BODY_ptr, FACE_list ))
                {
                    // If copied_FACE_ptr is the only FACE in this BODY, then we can
                    // try reversing the BODY before attempting the sweep again
                  if (FACE_list.size() == 1)
                  {
                      // Make sure this is a DOUBLE_SIDED FACE
                    if (copied_FACE_ptr->sides() == DOUBLE_SIDED)
                    {
                        // Reverse the BODY that this FACE belongs to
                      if (!api_reverse_body(copy_of_original_BODY_ptr).ok())
                      {
                        reverse_failed = CUBIT_TRUE;
                      }
                    }

                      // The single FACE in this BODY is *not* DOUBLE_SIDED.
                      // Cannot do the reversal and sweep.
                    else
                    {
                      reverse_failed = CUBIT_TRUE;
                    }
                  }

                    // More than 1 FACE in this BODY. Cannot do the reversal and sweep.
                  else
                  {
                    reverse_failed = CUBIT_TRUE;
                  }
                }

                  // Couldn't get the FACES of the BODY being swept. Cannot do
                  // the reversal and sweep.
                else
                {
                  reverse_failed = CUBIT_TRUE;
                }
              }

                // Couldn't find the FACE in the copied BODY that corresponds to
                // the FACE being swept in the original BODY.  Cannot do the
                // reversal and sweep.
              else
              {
                reverse_failed = CUBIT_TRUE;
              }
            }

              // The sweep failed for reasons other than the creation of a negative
              // volume.  Doesn't make sense to do the reversal and sweep.
            else
            {
              reverse_failed = CUBIT_TRUE;
            }

              // If the reverse was successful, then try the sweep again
            if (reverse_failed == CUBIT_FALSE)
            {
              if ( !this->sweep_FACE_about_axis( copied_FACE_ptr,
                                                 sweep_axis,
                                                 point,
                                                 angle,
                                                 volume_is_negative,
                                                 primary_side,
                                                 steps,
                                                 draft_angle,
                                                 draft_type,
                                                 rigid) )
              {
                  // The attempt to reverse the BODY and sweep it again failed.
                PRINT_ERROR("Problem sweeping one of the surfaces.\n");

                  // Delete the BODY created by this sweep
                AcisQueryEngine::instance()->delete_ACIS_BODY( copy_of_original_BODY_ptr, CUBIT_TRUE);
                copy_of_original_BODY_ptr = NULL;

                skip_loop = CUBIT_TRUE;
              }
            }

            else
            {
              skip_loop = CUBIT_TRUE;
            }
          } // End of if branch for initial sweep attemp
        }

          // There were problems getting the FACE in the copied BODY, corresponding
          // to the FACE being swept in the original BODY. We can't proceed with the
          // sweep of this RefFace.
        else
        {
          skip_loop = CUBIT_TRUE;
        } // End of if branch if copied_FACE_ptr exists

          // As all went well (apparently :-), create a new Body using the new
          // swept solid and get rid of the original Body and its associated
          // ACIS BODY.
        if( !skip_loop )
        {
          BodySM* new_body_ptr = get_new_Body(body_ptr, BODY_ptr,
                                            copy_of_original_BODY_ptr,
                                            CUBIT_FALSE, CUBIT_FALSE);
          if (new_body_ptr)
          {
            sweep_status = CUBIT_SUCCESS;
            result_body_list.append(new_body_ptr);
          }
        } // End of if branch if skip_loop != TRUE of final new_body_ptr construction
        else
        {
            // Get rid of the copied body, if it exists
          if( copy_of_original_BODY_ptr )
          {
            AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copy_of_original_BODY_ptr);
            AcisQueryEngine::instance()->delete_ACIS_BODY( copy_of_original_BODY_ptr);
            copy_of_original_BODY_ptr = NULL;
          }
        }
      } // End of if branch if skip_loop != TRUE for problems finding assoc. ACIS BODY
    } // End of if branch if skip_loop != TRUE for non-planar FACE_ptr
  }


  for (i = 0; i < EDGEs_to_be_swept_list.size(); i++)
  {
      // Get the EDGE we're going to sweep
    copied_EDGE_ptr = EDGEs_to_be_swept_list.get_and_step();

    BODY* new_BODY_ptr = NULL;
    sweep_status = sweep_EDGE_about_axis( copied_EDGE_ptr,
                                          new_BODY_ptr,
                                          sweep_axis,
                                          point,
                                          angle,
                                          steps,
                                          draft_angle,
                                          draft_type,
                                          make_solid,
                                          rigid);

    if(sweep_status == CUBIT_FAILURE)
    {
      PRINT_ERROR("Problem sweeping curve.\n");
      continue;
    }

    BodySM* new_body_ptr = get_new_Body(body_ptr, BODY_ptr,
                                      new_BODY_ptr,
                                      CUBIT_FALSE, CUBIT_FALSE);
    if (new_body_ptr)
    {
      result_body_list.append(new_body_ptr);
      sweep_status = CUBIT_SUCCESS;
    }
  } // End of i loop, iterating over EDGEs_to_be_swept_list

  return sweep_status;
}

//-------------------------------------------------------------------------
// Purpose       : Perform a translational sweep operation
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/31/96
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::sweep_translational(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  const CubitVector& sweep_vector,
  double draft_angle,
  int draft_type,
  bool switchside,
  bool rigid) const
{
  DLIList<FACE*> FACEs_to_be_swept_list;
  DLIList<EDGE*> EDGEs_to_be_swept_list;
  DLIList<Surface*> ref_faces_list;
  DLIList<Curve*> ref_edges_list;

    // First get the ACIS FACEs associated with the input RefFaces
  //DLIList<FACE*> temp_faces;
  int i;
  for (i = ref_ent_list.size(); i > 0; i--)
  {
    GeometryEntity *ref_ent = ref_ent_list.get_and_step();
    Surface *ref_face = CAST_TO(ref_ent, Surface);
    Curve* ref_edge = CAST_TO(ref_ent, Curve);

    if(ref_face != NULL)
    {
      //SurfaceACIS::get_FACEs_of_RefFace( ref_face, temp_faces);
      //if (temp_faces.size() > 1) {
      //  PRINT_ERROR("Surface %d cannot be swept as the number of ACIS FACEs associated \n"
      //              "       with it is not 1.\n",
      //              ref_face->id());
      //  return CUBIT_FAILURE;
      //}
      //if (temp_faces.size() > 0)
      //   FACEs_to_be_swept_list += temp_faces;
      //temp_faces.clean_out();
      FACE* FACE_ptr = AcisQueryEngine::get_FACE(ref_face);
      if (!FACE_ptr)
        return CUBIT_FAILURE;
      FACEs_to_be_swept_list.append(FACE_ptr);
      ref_faces_list.append(ref_face);
    }
    else if (ref_edge != NULL)
    {
      //EDGEs_to_be_swept_list.append(CurveACIS::get_first_EDGE(ref_edge));
      EDGEs_to_be_swept_list.append(AcisQueryEngine::get_EDGE(ref_edge));
      ref_edges_list.append(ref_edge);
    }

  }

  ref_edges_list.reset();
  ref_faces_list.reset();

    // Now sweep them along the sweep_vector, a distance equal to the length of
    // the sweep SPAvector, with the input draft angle and draft type.

  FACEs_to_be_swept_list.reset();
  EDGEs_to_be_swept_list.reset();
  FACE* FACE_ptr = NULL;
  FACE* copied_FACE_ptr = NULL;
  EDGE* copied_EDGE_ptr = NULL;
  Surface* refface_ptr = NULL;
  BodySM* body_ptr = NULL;
  BODY* BODY_ptr = NULL;
  BODY* copy_of_original_BODY_ptr = NULL;
  bool primary_side = true;
  bool skip_loop = false;
  CubitStatus sweep_status = CUBIT_FAILURE;
  for (i = 0; i < FACEs_to_be_swept_list.size(); i++)
  {
    skip_loop = CUBIT_FALSE;
    copied_FACE_ptr = NULL;

      // Get the RefFace we're going to sweep
    refface_ptr = ref_faces_list.get_and_step();

      // Get the FACE we're going to sweep
    FACE_ptr = FACEs_to_be_swept_list.get_and_step();

      // If the FACE is DOUBLE-SIDED and the switchside option is TRUE
    if (FACE_ptr->sides() && switchside)
    {
      primary_side = CUBIT_FALSE;
    }

      // Get the associated Body and its ACIS BODY
    body_ptr = AcisQueryEngine::instance()->get_body_sm_of_ENTITY(FACE_ptr);
    if (body_ptr != NULL)
    {
      BODY_ptr = AcisQueryEngine::get_BODY(body_ptr);

        // Make sure we have a valid ACIS BODY associated with this Body
        // and FACE
      if (BODY_ptr == NULL)
      {
        PRINT_ERROR("No valid volume associated with the RefFace.\n");
        skip_loop = CUBIT_TRUE;
      }
    }

    else
    {
      PRINT_WARNING("No volume associated with this surface.\n"
                    "         Sweeping will continue with the ACIS BODY.\n");
      BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY(FACE_ptr);
      if (BODY_ptr == NULL)
      {
        PRINT_ERROR("No valid volume associated with the RefFace.\n");
        skip_loop = CUBIT_TRUE;
      }
    }

      // If there were no problems finding an associated ACIS BODY, proceed
    if( !skip_loop )
    {
        // Make a copy of the original BODY before doing the sweep.
        // We will perform the sweep on the copy, not on the original body.
      copy_of_original_BODY_ptr = this->copy_BODY(BODY_ptr, CUBIT_FALSE);

        // Find the FACE in the copied BODY that corresponds to the FACE
        // being swept in the original BODY.  Note that when we made a copy
        // of the original BODY, all the "parent" attributes were copied
        // over to the new BODY.
      if (copy_of_original_BODY_ptr)
      {
        copied_FACE_ptr = (FACE *) (AcisQueryEngine::instance()->
          get_ENTITY_of_entity( copy_of_original_BODY_ptr, refface_ptr ));
      }

      if (copied_FACE_ptr)
      {
          // Remove all the CUBIT_OWNER attributes from the copied BODY
          // as we're going to sweep it.
        if (GeometryModifyTool::instance()->get_new_ids() == CUBIT_TRUE)
          AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copy_of_original_BODY_ptr);

          // Sweep the FACE along the SPAvector, sweep_vector
        bool reverse_failed = false;
        bool volume_is_negative = false;

        if ( this->sweep_FACE_along_vector( copied_FACE_ptr,
                                            sweep_vector,
                                            volume_is_negative,
                                            primary_side,
                                            draft_angle,
                                            draft_type,
                                            rigid) == CUBIT_FAILURE)
        {
            // Volume was negative, so we're going to try to reverse the
            // BODY
          if (volume_is_negative)
          {
              // Delete the BODY that was just created (the one with the
              // negative volume).
            AcisQueryEngine::instance()->delete_ACIS_BODY(copy_of_original_BODY_ptr, CUBIT_TRUE);

              // Make another copy of the original BODY before attempting
              // the sweep again. We will perform the sweep on this copy,
              // not on the original body.
            copy_of_original_BODY_ptr = this->copy_BODY(BODY_ptr, CUBIT_FALSE);

              // Find the FACE in the copied BODY that corresponds to the FACE
              // being swept in the original BODY.  Note that when we made a copy
              // of the original BODY, all the "parent" attributes were copied
              // over to the new BODY.
            if (copy_of_original_BODY_ptr)
            {
              copied_FACE_ptr = (FACE *) (AcisQueryEngine::instance()->
                get_ENTITY_of_entity( copy_of_original_BODY_ptr, refface_ptr ));
            }

            if (copied_FACE_ptr)
            {
                // Remove all the CUBIT_OWNER attributes from the copied BODY
                // as we're going to sweep it.
              if (GeometryModifyTool::instance()->get_new_ids() == CUBIT_TRUE)
                AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copy_of_original_BODY_ptr);

                // Find out how many FACEs exist in the BODY that owns
                // copied_FACE_ptr
              DLIList<FACE*> FACE_list;
              if (AcisQueryEngine::instance()->get_FACEs( (ENTITY *)copy_of_original_BODY_ptr, FACE_list ))
              {
                  // If copied_FACE_ptr is the only FACE in this BODY, then we can
                  // try reversing the BODY before attempting the sweep again
                if (FACE_list.size() == 1)
                {
                    // Make sure this is a DOUBLE_SIDED FACE
                  if (copied_FACE_ptr->sides() == DOUBLE_SIDED)
                  {
                      // Reverse the BODY that this FACE belongs to
                    if (!api_reverse_body(copy_of_original_BODY_ptr).ok())
                    {
                      reverse_failed = CUBIT_TRUE;
                    }
                  }

                    // The single FACE in this BODY is *not* DOUBLE_SIDED.
                    // Cannot do the reversal and sweep.
                  else
                  {
                    reverse_failed = CUBIT_TRUE;
                  }
                }

                  // More than 1 FACE in this BODY. Cannot do the reversal and sweep.
                else
                {
                  reverse_failed = CUBIT_TRUE;
                }
              }

                // Couldn't get the FACES of the BODY being swept. Cannot do
                // the reversal and sweep.
              else
              {
                reverse_failed = CUBIT_TRUE;
              }
            }

              // Couldn't find the FACE in the copied BODY that corresponds to
              // the FACE being swept in the original BODY.  Cannot do the
              // reversal and sweep.
            else
            {
              reverse_failed = CUBIT_TRUE;
            }
          }

            // The sweep failed for reasons other than the creation of a negative
            // volume.  Doesn't make sense to do the reversal and sweep.
          else
          {
            reverse_failed = CUBIT_TRUE;
          }

            // If reverse was successful, then try the sweep again
          if ( reverse_failed == CUBIT_FALSE )
          {
            if ( this->sweep_FACE_along_vector( copied_FACE_ptr,
                                                sweep_vector,
                                                volume_is_negative,
                                                primary_side,
                                                draft_angle,
                                                draft_type,
                                                rigid) == CUBIT_FAILURE )
            {
                // The attempt to reverse the BODY and sweep it again failed.
              PRINT_ERROR("Problem sweeping one of the surfaces.\n");

                // Delete the BODY created by this sweep
              AcisQueryEngine::instance()->delete_ACIS_BODY( copy_of_original_BODY_ptr, CUBIT_TRUE );

              skip_loop = CUBIT_TRUE;
            }
          }
          else
          {
            skip_loop = CUBIT_TRUE;
          }
        } // End of if branch for initial sweep attemp
      }

        // There were problems getting the FACE in the copied BODY, corresponding
        // to the FACE being swept in the original BODY. We can't proceed with the
        // sweep of this RefFace.
      else
      {
        skip_loop = CUBIT_TRUE;
      } // End of if branch if copied_FACE_ptr exists

        // As all went well (apparently :-), create a new Body using the new
        // swept solid and get rid of the original Body and its associated
        // ACIS BODY.
      if( !skip_loop )
      {
        BodySM* new_body_ptr = get_new_Body(body_ptr, BODY_ptr,
                                          copy_of_original_BODY_ptr,
                                          CUBIT_FALSE, CUBIT_FALSE);
        if (new_body_ptr)
        {
          result_body_list.append(new_body_ptr);
          sweep_status = CUBIT_SUCCESS;
        }

      } // End of if branch if skip_loop != TRUE of final new_body_ptr construction
      else
      {
          // Get rid of the copied body, if it exists
        if( copy_of_original_BODY_ptr )
        {
          AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copy_of_original_BODY_ptr);
          AcisQueryEngine::instance()->delete_ACIS_BODY( copy_of_original_BODY_ptr);

          copy_of_original_BODY_ptr = NULL;
        }
      }
    } // End of if branch if skip_loop != TRUE for problems finding assoc. ACIS BODY
  } // End of i loop, iterating over FACEs_to_be_swept_list


    //Now we iterate over the EDGEs
  for (i = 0; i < EDGEs_to_be_swept_list.size(); i++)
  {
      // Get the EDGE we're going to sweep
    copied_EDGE_ptr = EDGEs_to_be_swept_list.get_and_step();

    BODY* new_BODY_ptr = NULL;
    sweep_status = sweep_EDGE_along_vector( copied_EDGE_ptr,
                                            new_BODY_ptr,
                                            sweep_vector,
                                            draft_angle,
                                            draft_type,
                                            rigid);
    if(sweep_status == CUBIT_FAILURE)
    {
      PRINT_ERROR("Problem sweeping curve.\n");
      continue;
    }

    BodySM* new_body_ptr = get_new_Body(body_ptr, BODY_ptr,
                                      new_BODY_ptr,
                                      CUBIT_FALSE, CUBIT_FALSE);
    if (new_body_ptr)
    {
      sweep_status = CUBIT_SUCCESS;
      result_body_list.append(new_body_ptr);
    }

  } // End of i loop, iterating over EDGEs_to_be_swept_list

  return sweep_status;
}

CubitStatus AcisModifyEngine::sweep_perpendicular(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  double distance,
  double draft_angle,
  int draft_type,
  bool switchside,
  bool rigid) const
{
  DLIList<FACE*> FACEs_to_be_swept_list;

  DLIList<Surface*> ref_faces_list;
  DLIList<Curve*> ref_edges_list;

    // First get the ACIS FACEs associated with the input RefFaces
  //DLIList<FACE*> temp_faces;
  int i;
  for (i = ref_ent_list.size(); i > 0; i--)
  {
    GeometryEntity *ref_ent = ref_ent_list.get_and_step();
    Surface *ref_face = CAST_TO(ref_ent, Surface);
    Curve* ref_edge = CAST_TO(ref_ent, Curve);

    if(ref_face != NULL)
    {
      //SurfaceACIS::get_FACEs_of_RefFace( ref_face, temp_faces);
      //if (temp_faces.size() > 1) {
      //  PRINT_ERROR("Surface %d cannot be swept as the number of ACIS FACEs associated \n"
      //              "       with it is not 1.\n",
      //              ref_face->id());
      //  return CUBIT_FAILURE;
      //}
      //if (temp_faces.size() > 0)
      //   FACEs_to_be_swept_list += temp_faces;
      //temp_faces.clean_out();
      FACE* FACE_ptr = AcisQueryEngine::get_FACE(ref_face);
      if (!FACE_ptr)
        return CUBIT_FAILURE;
      FACEs_to_be_swept_list.append(FACE_ptr);
      ref_faces_list.append(ref_face);
    }
    else if (ref_edge != NULL)
    {
      ref_edges_list.append(ref_edge);
    }
  }

  if(ref_edges_list.size())
     PRINT_ERROR("Curves cannot be swept perpendicularly, please use the vector sweep.\n");

  ref_faces_list.reset();

    // Now sweep them perpendicularly, a distance equal to the
    // distance given with the input draft angle and
    // draft type.

  FACEs_to_be_swept_list.reset();
  FACE* FACE_ptr = NULL;
  FACE* copied_FACE_ptr = NULL;

  Surface* refface_ptr = NULL;
  BodySM* body_ptr = NULL;
  BODY* BODY_ptr = NULL;
  BODY* copy_of_original_BODY_ptr = NULL;
  bool primary_side = CUBIT_TRUE;
  int skip_loop = CUBIT_FALSE;
  CubitStatus sweep_status = CUBIT_FAILURE;
  for (i = 0; i < FACEs_to_be_swept_list.size(); i++)
  {
    skip_loop = CUBIT_FALSE;
    copied_FACE_ptr = NULL;

      // Get the RefFace we're going to sweep
    refface_ptr = ref_faces_list.get_and_step();

      // Get the FACE we're going to sweep
    FACE_ptr = FACEs_to_be_swept_list.get_and_step();

      // If the FACE is DOUBLE-SIDED and the switchside option is TRUE
    if (FACE_ptr->sides() && switchside)
    {
      primary_side = CUBIT_FALSE;
    }

      // Get the associated Body and its ACIS BODY
    body_ptr = AcisQueryEngine::instance()->get_body_sm_of_ENTITY(FACE_ptr);
    if (body_ptr != NULL)
    {
      BODY_ptr = AcisQueryEngine::get_BODY(body_ptr);

        // Make sure we have a valid ACIS BODY associated with this Body
        // and FACE
      if (BODY_ptr == NULL)
      {
        PRINT_ERROR("No valid volume associated with the RefFace.\n");
        skip_loop = CUBIT_TRUE;
      }
    }

    else
    {
      PRINT_WARNING("No volume associated with this surface.\n"
                    "         Sweeping will continue with the ACIS BODY.\n");
      BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY(FACE_ptr);
      if (BODY_ptr == NULL)
      {
        PRINT_ERROR("No valid volume associated with the RefFace.\n");
        skip_loop = CUBIT_TRUE;
      }
    }

      // If there were no problems finding an associated ACIS BODY, proceed
    if( !skip_loop )
    {
        // Make a copy of the original BODY before doing the sweep.
        // We will perform the sweep on the copy, not on the original body.
      copy_of_original_BODY_ptr = this->copy_BODY(BODY_ptr, CUBIT_FALSE);

        // Find the FACE in the copied BODY that corresponds to the FACE
        // being swept in the original BODY.  Note that when we made a copy
        // of the original BODY, all the "parent" attributes were copied
        // over to the new BODY.
      if (copy_of_original_BODY_ptr)
      {
        copied_FACE_ptr = (FACE *) (AcisQueryEngine::instance()->
          get_ENTITY_of_entity( copy_of_original_BODY_ptr, refface_ptr ));
      }

      if (copied_FACE_ptr)
      {
          // Remove all the CUBIT_OWNER attributes from the copied BODY
          // as we're going to sweep it.
        if (GeometryModifyTool::instance()->get_new_ids() == CUBIT_TRUE)
          AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copy_of_original_BODY_ptr);

          // Sweep the FACE along the SPAvector, sweep_vector
        bool reverse_failed = CUBIT_FALSE;
        bool volume_is_negative = CUBIT_FALSE;

        if ( this->sweep_FACE_perpendicular( copied_FACE_ptr,
                                             distance,
                                             volume_is_negative,
                                             primary_side,
                                             draft_angle,
                                             draft_type,
                                             rigid) == CUBIT_FAILURE)
        {
            // Volume was negative, so we're going to try to reverse the
            // BODY
          if (volume_is_negative)
          {
              // Delete the BODY that was just created (the one with the
              // negative volume).
            AcisQueryEngine::instance()->delete_ACIS_BODY(copy_of_original_BODY_ptr, CUBIT_TRUE);

              // Make another copy of the original BODY before attempting
              // the sweep again. We will perform the sweep on this copy,
              // not on the original body.
            copy_of_original_BODY_ptr = this->copy_BODY(BODY_ptr, CUBIT_FALSE);

              // Find the FACE in the copied BODY that corresponds to the FACE
              // being swept in the original BODY.  Note that when we made a copy
              // of the original BODY, all the "parent" attributes were copied
              // over to the new BODY.
            if (copy_of_original_BODY_ptr)
            {
              copied_FACE_ptr = (FACE *) (AcisQueryEngine::instance()->
                get_ENTITY_of_entity( copy_of_original_BODY_ptr, refface_ptr ));
            }

            if (copied_FACE_ptr)
            {
                // Remove all the CUBIT_OWNER attributes from the copied BODY
                // as we're going to sweep it.
              if (GeometryModifyTool::instance()->get_new_ids() == CUBIT_TRUE)
                AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copy_of_original_BODY_ptr);

                // Find out how many FACEs exist in the BODY that owns
                // copied_FACE_ptr
              DLIList<FACE*> FACE_list;
              if (AcisQueryEngine::instance()->get_FACEs( (ENTITY *)copy_of_original_BODY_ptr, FACE_list ))
              {
                  // If copied_FACE_ptr is the only FACE in this BODY, then we can
                  // try reversing the BODY before attempting the sweep again
                if (FACE_list.size() == 1)
                {
                    // Make sure this is a DOUBLE_SIDED FACE
                  if (copied_FACE_ptr->sides() == DOUBLE_SIDED)
                  {
                      // Reverse the BODY that this FACE belongs to
                    if (!api_reverse_body(copy_of_original_BODY_ptr).ok())
                    {
                      reverse_failed = CUBIT_TRUE;
                    }
                  }

                    // The single FACE in this BODY is *not* DOUBLE_SIDED.
                    // Cannot do the reversal and sweep.
                  else
                  {
                    reverse_failed = CUBIT_TRUE;
                  }
                }

                  // More than 1 FACE in this BODY. Cannot do the reversal and sweep.
                else
                {
                  reverse_failed = CUBIT_TRUE;
                }
              }

                // Couldn't get the FACES of the BODY being swept. Cannot do
                // the reversal and sweep.
              else
              {
                reverse_failed = CUBIT_TRUE;
              }
            }

              // Couldn't find the FACE in the copied BODY that corresponds to
              // the FACE being swept in the original BODY.  Cannot do the
              // reversal and sweep.
            else
            {
              reverse_failed = CUBIT_TRUE;
            }
          }

            // The sweep failed for reasons other than the creation of a negative
            // volume.  Doesn't make sense to do the reversal and sweep.
          else
          {
            reverse_failed = CUBIT_TRUE;
          }

            // If reverse was successful, then try the sweep again
          if ( reverse_failed == CUBIT_FALSE )
          {
            if ( this->sweep_FACE_perpendicular( copied_FACE_ptr,
                                                 distance,
                                                 volume_is_negative,
                                                 primary_side,
                                                 draft_angle,
                                                 draft_type,
                                                 rigid) == CUBIT_FAILURE )
            {
                // The attempt to reverse the BODY and sweep it again failed.
              PRINT_ERROR("Problem sweeping one of the surfaces.\n");

                // Delete the BODY created by this sweep
              AcisQueryEngine::instance()->delete_ACIS_BODY( copy_of_original_BODY_ptr, CUBIT_TRUE );

              skip_loop = CUBIT_TRUE;
            }
          }
          else
          {
            skip_loop = CUBIT_TRUE;
          }
        } // End of if branch for initial sweep attemp
      }

        // There were problems getting the FACE in the copied BODY, corresponding
        // to the FACE being swept in the original BODY. We can't proceed with the
        // sweep of this RefFace.
      else
      {
        skip_loop = CUBIT_TRUE;
      } // End of if branch if copied_FACE_ptr exists

        // As all went well (apparently :-), create a new Body using the new
        // swept solid and get rid of the original Body and its associated
        // ACIS BODY.
      if( !skip_loop )
      {
        BodySM* new_body_ptr = get_new_Body(body_ptr, BODY_ptr,
                                          copy_of_original_BODY_ptr,
                                          CUBIT_FALSE, CUBIT_FALSE);
        if (new_body_ptr)
        {
          result_body_list.append(new_body_ptr);
          sweep_status = CUBIT_SUCCESS;
        }

      } // End of if branch if skip_loop != TRUE of final new_body_ptr construction
      else
      {
          // Get rid of the copied body, if it exists
        if( copy_of_original_BODY_ptr )
        {
          AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copy_of_original_BODY_ptr);
          AcisQueryEngine::instance()->delete_ACIS_BODY( copy_of_original_BODY_ptr);

          copy_of_original_BODY_ptr = NULL;
        }
      }
    } // End of if branch if skip_loop != TRUE for problems finding assoc. ACIS BODY
  } // End of i loop, iterating over FACEs_to_be_swept_list
  return sweep_status;
}

//-------------------------------------------------------------------------
// Purpose       : Perform a face sweep operation along a curve
//
// Special Notes : In ACIS v4.2, I was unable to get it to sweep along
//                 multiple curves, so for now it gives an error message
//                 if this is attempted.
//
//                 The function will not allow you to sweep back into
//                 a solid.
//
// Creator       : Steve Storm
//
// Creation Date : 8/7/98
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::sweep_along_curve(
  DLIList<GeometryEntity*>& ref_ent_list,
  DLIList<BodySM*>& result_body_list,
  DLIList<Curve*>& ref_edge_list,
  double draft_angle,
  int draft_type,
  bool rigid) const
{
   outcome result;
   DLIList<FACE*> FACEs_to_be_swept_list;
   DLIList<EDGE*> EDGEs_to_be_swept_list;
   DLIList<Surface*> ref_faces_list;
   DLIList<Curve*> ref_edges_list;


   // First get the ACIS FACEs associated with the input RefFaces
//   DLIList<FACE*> temp_faces;
   int i;
   for (i = ref_ent_list.size(); i > 0; i--)
   {
     GeometryEntity* ref_ent = ref_ent_list.get_and_step();
     Surface *ref_face = CAST_TO(ref_ent, Surface);
     Curve* ref_edge = CAST_TO(ref_ent, Curve);
     if(ref_face != NULL)
     {
//       SurfaceACIS::get_FACEs_of_RefFace( ref_face, temp_faces);
//       if (temp_faces.size() > 1) {
//         PRINT_ERROR("Surface %d cannot be swept as the number of ACIS FACEs associated \n"
//                     "       with it is not 1.\n",
//                     ref_face->id());
 //        return CUBIT_FAILURE;
//       }
//       if (temp_faces.size() > 0)
//          FACEs_to_be_swept_list += temp_faces;
//       temp_faces.clean_out();
       FACE* FACE_ptr = AcisQueryEngine::get_FACE(ref_face);
       if (!FACE_ptr)
        return CUBIT_FAILURE;
       FACEs_to_be_swept_list.append(FACE_ptr);
       ref_faces_list.append(ref_face);
     }
     else if (ref_edge != NULL)
     {
       //EDGEs_to_be_swept_list.append(CurveACIS::get_first_EDGE(ref_edge));
       EDGEs_to_be_swept_list.append(AcisQueryEngine::get_EDGE(ref_edge));
       ref_edges_list.append(ref_edge);
     }

   }

   // Now sweep them along the EDGE.

   BODY* EDGE_BODY_ptr = NULL;

   FACEs_to_be_swept_list.reset();
   EDGEs_to_be_swept_list.reset();
   ref_edges_list.reset();

   FACE* FACE_ptr = NULL;
   FACE* copied_FACE_ptr = NULL;
   Surface* ref_face_ptr = NULL;
   BodySM* body_ptr = NULL;
   BODY* BODY_ptr = NULL;
   BODY* copied_BODY_ptr = NULL;

   // Get the EDGE to sweep along.  No need to copy it as it is going into
   // a WIRE body (it's copied into the WIRE body).
   //
   Curve* ref_edge_ptr = NULL;

     // Sweep the FACE along the chain of EDGES.
   int num_edges = ref_edge_list.size();
   EDGE** EDGES = new EDGE*[num_edges];
   ref_edge_list.reset();
   for( int j=0; j<num_edges; j++ )
   {
     ref_edge_ptr = ref_edge_list.get_and_step();
     //Curve* curve_ptr = ref_edge_ptr->get_curve_ptr();
     //CurveACIS *curve_ACIS = CAST_TO(curve_ptr, CurveACIS);
     CurveACIS *curve_ACIS = CAST_TO(ref_edge_ptr, CurveACIS);
     EDGE *EDGE_ptr = curve_ACIS->get_EDGE_ptr();

     if (!EDGE_ptr)
     {
       PRINT_ERROR("Unable to get ACIS EDGE from Curve\n" );
       delete [] EDGES;
       return CUBIT_FAILURE;
     }

     EDGE *copied_EDGE_ptr = NULL;
     result = api_edge( EDGE_ptr, copied_EDGE_ptr );
     if (!result.ok() || copied_EDGE_ptr == NULL )
     {
       AcisQueryEngine::instance()->ACIS_API_error (result);
       PRINT_ERROR( "Unable to copy ACIS curve\n" );
       delete [] EDGES;
       return CUBIT_FAILURE;
     }
     EDGES[j] = copied_EDGE_ptr;
   }

     // Make a WIRE body to sweep along.  If we have to reverse the EDGE
     // (since ACIS only allows you to sweep from start to end), the wire
     // will be needed.
   result = api_make_ewire( num_edges, EDGES, EDGE_BODY_ptr );
   WIRE* WIRE_ptr = NULL;
   if (!result.ok())
   {
     AcisQueryEngine::instance()->ACIS_API_error (result);
     PRINT_ERROR( "Unable to make ACIS wire body from curve to sweep along\n" );
   }
   else
      WIRE_ptr = EDGE_BODY_ptr->lump()->shell()->wire();

   CubitStatus sweep_status = CUBIT_FAILURE;
   for (i = 0; i < FACEs_to_be_swept_list.size(); i++)
   {
      copied_FACE_ptr = NULL;

      // Get the RefFace we're going to sweep
      ref_face_ptr = ref_faces_list.get_and_step();

      // Get the FACE we're going to sweep
      FACE_ptr = FACEs_to_be_swept_list.get_and_step();

      // Get the associated Body and its ACIS BODY
      body_ptr = AcisQueryEngine::instance()->get_body_sm_of_ENTITY(FACE_ptr);
      if (body_ptr != NULL)
      {
         BODY_ptr = AcisQueryEngine::get_BODY(body_ptr);

         // Make sure we have a valid ACIS BODY associated with this Body
         // and FACE
         if (BODY_ptr == NULL)
         {
            PRINT_ERROR("No valid volume associated with Surface.\n" );
            continue;
         }
      }

      else
      {
          // This can happen when a Surface entity is created at the SM level for temp use
          // This function should only care about the ACIS BODY anyway
         BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY(FACE_ptr);
         if (BODY_ptr == NULL)
         {
            PRINT_ERROR("No valid volume associated with Surface\n" );
            continue;
         }
      }

      // Make a copy of the original BODY before doing the sweep.
      // We will perform the sweep on the copy, not on the original body.
      copied_BODY_ptr = this->copy_BODY(BODY_ptr, CUBIT_FALSE);

      if( copied_BODY_ptr == NULL )
      {
         PRINT_ERROR( "Unable to copy volume before sweeping\n" );
         continue;
      }

      // Find the FACE in the copied BODY that corresponds to the FACE
      // being swept in the original BODY.  Note that when we made a copy
      // of the original BODY, all the "parent" attributes were copied
      // over to the new BODY.
      copied_FACE_ptr = (FACE *) (AcisQueryEngine::instance()->get_ENTITY_of_entity(
         copied_BODY_ptr, ref_face_ptr ) );

      if( copied_FACE_ptr == NULL )
      {
         PRINT_ERROR( "Unable to find FACE to sweep in copied volume\n" );
         continue;
      }

        // Remove all the CUBIT_OWNER attributes from the copied BODY
        // as we're going to sweep it.
      if (GeometryModifyTool::instance()->get_new_ids() == CUBIT_TRUE)
        AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);

      // Check to see if the start of the edge is in the profile of
      // the face.  If not, check the end - if this is, we'll need to
      // reverse the WIRE body.
      //
      // When we support a chain of EDGES, here we will need to get the start
      // and end vertices of the WIRE.  This should be the only change needed.
      //VERTEX* start_VERTEX_ptr;// copied_EDGE_ptr->start();
      //VERTEX* end_VERTEX_ptr;// copied_EDGE_ptr->end();
      //if( start_VERTEX_ptr == end_VERTEX_ptr )
      //{
       //  PRINT_ERROR( "Edge start and end are the same - sweeping not possible\n" );
       //  return CUBIT_FAILURE;
      //}
      ENTITY *ent = WIRE_ptr;

      SPAposition acis_start_coords;
      get_wire_start_position(ent, acis_start_coords);
      CubitVector start_vec( acis_start_coords.x(), acis_start_coords.y(),
        acis_start_coords.z() );

      bool start_flg;
      Curve* start_curve = find_curve_by_end_coord( start_vec,
                                 ref_edge_list, start_flg );
      if( !start_curve )
      {
        PRINT_ERROR( "unable to find start vertex of sweep chain; aborting\n" );
        delete []EDGES;
        return CUBIT_FAILURE;
      }
      CubitVector start_position;
      start_curve->position_from_fraction( start_flg ? 0.0 : 1.0, start_position);

      SPAposition acis_end_coords;
      get_wire_end_position(ent, acis_end_coords);
      CubitVector end_vec( acis_end_coords.x(), acis_end_coords.y(),
        acis_end_coords.z() );

      Curve* end_curve = find_curve_by_end_coord( end_vec,
                              ref_edge_list, start_flg );
      if( !end_curve )
      {
        PRINT_ERROR( "unable to find end vertex of sweep chain; aborting\n" );
        delete []EDGES;
        return CUBIT_FAILURE;
      }
      CubitVector end_position;
      end_curve->position_from_fraction( start_flg ? 0.0 : 1.0, end_position );

      //SPAposition pos1, pos2;
      double distance1, distance2;
      SPAposition acis_closest_point;
      (copied_FACE_ptr->geometry()->equation()).point_perp ( acis_start_coords,
                                                             acis_closest_point );
      CubitVector closest_point_start( acis_closest_point.x(),
                                       acis_closest_point.y(),
                                       acis_closest_point.z() );
      distance1 = closest_point_start.distance_between( start_vec );

      //VERTEX* close_VERTEX_ptr = NULL;
      CubitVector closest_position;
      EDGE* start_EDGE_ptr = NULL;

      if( distance1 > SPAresabs )
      {
         // Try the other end
         (copied_FACE_ptr->geometry()->equation()).point_perp ( acis_end_coords,
                                                                acis_closest_point );
         CubitVector closest_point_end( acis_closest_point.x(), acis_closest_point.y(),
                                        acis_closest_point.z() );
         distance2 = closest_point_end.distance_between( end_vec );

         if( distance2 > SPAresabs )
         {
            PRINT_ERROR( "Neither chain end lies in the plane of surface\n" );
            continue;
         }
         else
         {
            //close_VERTEX_ptr = end_VERTEX_ptr;
           closest_position = end_position;
           start_curve = end_curve;
           start_EDGE_ptr = EDGES[ref_edge_list.size()-1];
         }

         // Reverse the WIRE in the BODY
         //api_reverse_body( EDGE_BODY_ptr );
         api_reverse_wire( WIRE_ptr );

      }
      else
      {
         //close_VERTEX_ptr = start_VERTEX_ptr;
        closest_position = start_position;
        start_EDGE_ptr = EDGES[0];
      }

      SPAposition acis_closest_position;
      SPAunit_vector acis_surface_normal;

      //APOINT* acis_point = close_VERTEX_ptr->geometry() ;
      //SPAposition acis_curve_position = acis_point->coords() ;
      SPAposition acis_curve_position( closest_position.x(),
        closest_position.y(), closest_position.z() );

      // Need to compute the FACE normal at the start of the sweep
      (copied_FACE_ptr->geometry()->equation()).point_perp ( acis_curve_position,
         acis_closest_position, acis_surface_normal );

      CubitVector surface_normal( acis_surface_normal.x(), acis_surface_normal.y(),
                                  acis_surface_normal.z() );

      // Check to see if SENSE of FACE is REVERSED.  The direction of the
      // ACIS face is the direction of sweep.
      if( copied_FACE_ptr->sense() == REVERSED )
      {
         if( copied_FACE_ptr->geometry()->equation().left_handed_uv() == CUBIT_FALSE )
            surface_normal *= -1.0;
      }
      else
      {
         if( copied_FACE_ptr->geometry()->equation().left_handed_uv() )
            surface_normal *= -1.0;
      }

      // Get curve tangent at start of the sweep
      SPAunit_vector acis_curve_tangent;

      (start_EDGE_ptr->geometry()->equation()).point_perp( acis_curve_position,
         acis_closest_position,
         acis_curve_tangent );

      CubitVector curve_tangent( acis_curve_tangent.x(), acis_curve_tangent.y(),
         acis_curve_tangent.z() );

      // Determine if the WIRE has been reversed - ??? there's a better way to
      // do this but I don't know it.  I would think you could just go by the
      // underlying EDGE but apparently when the WIRE get's reversed, the EDGE
      // doesn't change - the sweep uses the WIRE, not the EDGE.
      if( closest_position == end_position )
         curve_tangent *= -1.0;

      // Make sure the sense is correct based on the underlying EDGE sense
      if( start_EDGE_ptr->sense() == REVERSED )
         curve_tangent *= -1.0;


      // Properties of dot product
      // angle acute if dot product > 0
      // angle obtuse if dot product < 0
      // angle 90 deg if dot product = 0

      if( surface_normal%curve_tangent < 0 )
      {
         if( FACE_ptr->sides() != DOUBLE_SIDED )
         {
            PRINT_ERROR( "Attempting to sweep Surface back into the solid\n"
               "       Sweeping of Surface aborted\n" );
            api_delent(copied_FACE_ptr);
            continue;
         }
         api_reverse_face( copied_FACE_ptr );
      }

      sweep_status = sweep_FACE_along_WIRE( copied_FACE_ptr,
                                            EDGE_BODY_ptr,
                                            draft_angle,
                                            draft_type,
                                            rigid);
      if( sweep_status == CUBIT_FAILURE )
      {
         // Get rid of the copied body, if it exists
         if( copied_BODY_ptr )
         {
            AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);
            AcisQueryEngine::instance()->delete_ACIS_BODY( copied_BODY_ptr );
            copied_BODY_ptr = NULL;
            continue;
         }
      }

      // Create a new Body using the new swept solid and delete the
      // original Body and its associated ACIS BODY.
      BodySM* new_body_ptr = 0;
      if(body_ptr)
      new_body_ptr = get_new_Body(body_ptr, BODY_ptr, copied_BODY_ptr,
                                        CUBIT_FALSE, CUBIT_FALSE);
      else // this handles the case of a NULL body sweep surface (typically the result of a SM only entity)
           // in this case the caller must properly handle the sweep surface deletion
          new_body_ptr = AcisQueryEngine::instance()->populate_topology_bridges(copied_BODY_ptr);

      if (new_body_ptr)
      {
        sweep_status = CUBIT_SUCCESS;
        result_body_list.append(new_body_ptr);
      }
   } // End of i loop, iterating over FACEs_to_be_swept_list
     //Now let's iterate over the EDGEs_to_be_swept_list
   for (i = 0; i < EDGEs_to_be_swept_list.size(); i++)
   {
     EDGE *copied_EDGE_ptr = EDGEs_to_be_swept_list.get_and_step();
     BODY* new_BODY_ptr = NULL;

      sweep_status = sweep_EDGE_along_WIRE( copied_EDGE_ptr,
                                            EDGE_BODY_ptr,
                                            new_BODY_ptr,
                                            draft_angle,
                                            draft_type,
                                            rigid);
      if( sweep_status == CUBIT_FAILURE )
      {
        PRINT_ERROR("Problem sweeping curve.\n");
        continue;

      }

      // Create a new Body using the new swept solid

      BodySM* new_body_ptr = get_new_Body(body_ptr, BODY_ptr, new_BODY_ptr,
                                        CUBIT_FALSE, CUBIT_FALSE);
      if (new_body_ptr)
      {
         sweep_status = CUBIT_SUCCESS;
         result_body_list.append(new_body_ptr);
      }

   } // End of i loop, iterating over FACEs_to_be_swept_list
   delete []EDGES;

     // Remove the created WIRE body
   api_delent( EDGE_BODY_ptr );


   return sweep_status;
}

CubitStatus AcisModifyEngine::webcut_with_cylinder(
                                              DLIList<BodySM*> &webcut_body_list,
                                              double radius,
                                              const CubitVector &axis,
                                              const CubitVector &center,
                                              DLIList<BodySM*> &results_list,
                                              bool imprint )
{
     //create the cylinder body.
   SPAbox super_box;
   CubitStatus rv = AcisQueryEngine::instance()->create_super_acis_bounding_box( webcut_body_list, super_box );
   if (rv != CUBIT_SUCCESS)
      return CUBIT_FAILURE;
   double max_size =  AcisQueryEngine::instance()->get_max_size_of_acis_box( super_box );
     //lets find the distance to the center for each body and take
     //the max.
   double curr;
   CubitVector cent_bod;
   BodySM *body_ptr;
   for ( int ii = webcut_body_list.size(); ii > 0; ii-- )
   {
      body_ptr = webcut_body_list.get_and_step();
      cent_bod = AcisQueryEngine::instance()->bounding_box(body_ptr).center();
      cent_bod = cent_bod - center;
      curr = cent_bod.length();
      if ( curr > max_size )
          max_size = curr;
   }
   double height;
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

   BODY* cutting_tool_ptr;


   outcome result = api_make_frustum( height, radius, radius, radius, cutting_tool_ptr );
   if ( !result.ok() || cutting_tool_ptr == NULL ||
        cutting_tool_ptr->lump() == NULL )
   {

      PRINT_ERROR("In AcisModifyEngine::webcut_with_cylinder\n"
                  "       ACIS api_make_frustum function failed.\n\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return CUBIT_FAILURE;
   }
     // The current frustum is centered on the z axis.
   CubitVector axis2(0., 0., 1.0 );
     //now find the normal to the current axis and axis we want to be
     //at. This normal is where we will rotate about.
   CubitVector normal_axis = axis2 * axis;
   if ( normal_axis.length() > CUBIT_RESABS )
   {
        //angle in degrees.
      double angle = normal_axis.vector_angle( axis2, axis );
      SPAvector rotation_axis( normal_axis.x(), normal_axis.y(),
                            normal_axis.z() );
      SPAtransf transformation = rotate_transf( angle, rotation_axis );

      result = api_apply_transf( cutting_tool_ptr, transformation );
      if (!result.ok())
      {
         AcisQueryEngine::instance()->ACIS_API_error ( result, "rotating body" );
         return CUBIT_FAILURE;
      }
   }

   if ( center.length() > CUBIT_RESABS )
   {
        //move the body to the center SPAposition.
      SPAvector translation_vector( center.x(), center.y(), center.z() );
      SPAtransf transformation = translate_transf( translation_vector );
      result = api_apply_transf( cutting_tool_ptr, transformation );
      if (!result.ok())
      {
         AcisQueryEngine::instance()->ACIS_API_error ( result, "transforming body" );
         return CUBIT_FAILURE;
      }
   }
     //okay we should be ready now.

     // Use the BODY to perform webcut
   CubitStatus stat =
      this->webcut(webcut_body_list, cutting_tool_ptr, results_list, imprint ) ;

     // Delete the BODY that was created to be used as a tool
   AcisQueryEngine::instance()->delete_ACIS_BODY(cutting_tool_ptr, CUBIT_TRUE) ;
   cutting_tool_ptr = NULL;

   return stat;
}

CubitStatus
AcisModifyEngine::webcut_with_brick( DLIList<BodySM*>& webcut_body_list,
                                     const CubitVector &center,
                                     const CubitVector axes[3],
                                     const CubitVector &extension,
                                     DLIList<BodySM*> &results_list,
                                     bool imprint )
{
   // Create the brick to cut with
   BODY *cutting_tool_ptr = make_brick_BODY( center, axes, extension );
   if( cutting_tool_ptr == NULL )
      return CUBIT_FAILURE;

   CubitStatus stat =
       this->webcut(webcut_body_list, cutting_tool_ptr, results_list, imprint) ;

     // Delete the BODY that was created to be used as a tool
   AcisQueryEngine::instance()->delete_ACIS_BODY(cutting_tool_ptr, CUBIT_TRUE) ;

   return stat;
}

CubitStatus
AcisModifyEngine::webcut_with_planar_sheet( DLIList<BodySM*>& webcut_body_list,
                                            const CubitVector &center,
                                            const CubitVector axes[2],
                                            double width,
                                            double height,
                                            DLIList<BodySM*> &results_list,
                                            bool imprint )
{
   // Create the planar sheet to cut with
   CubitVector p1, p2, p3, p4;

   // Get the corners of the sheet
   center.next_point( axes[0], width/2.0, p1 );
   p1.next_point( axes[1], -height/2.0, p1 );
   p1.next_point( axes[1], height, p2 );
   p2.next_point( axes[0], -width, p3 );
   p3.next_point( axes[1], -height, p4 );

   BODY* cutting_tool_ptr = make_planar_quad_BODY( p1, p2, p3, p4 );
   if( cutting_tool_ptr == NULL )
      return CUBIT_FAILURE;

   BodySM *webcut_body;
   DLIList<BodySM*> temp_new_bodies;
   int num_cut = 0;
   for ( int ii = webcut_body_list.size(); ii > 0; ii-- )
   {
      webcut_body = webcut_body_list.remove();
      BODY *webcut_BODY = AcisQueryEngine::get_BODY(webcut_body );
      //now webcut with it.
      BODY *new_BODY_1, *new_BODY_2;
      CubitStatus status = webcut_with_sheet( webcut_BODY,
         cutting_tool_ptr,
         new_BODY_1,
         new_BODY_2,
         imprint );
      if ( status == CUBIT_SUCCESS )
      {
         DLIList<BODY*> new_BODY_list;
         new_BODY_list.append(new_BODY_1);
         new_BODY_list.append(new_BODY_2);
         temp_new_bodies.clean_out();

         CubitStatus result = get_new_Body(webcut_body, webcut_BODY,
            new_BODY_list, temp_new_bodies, CUBIT_FALSE);
         if ( result != CUBIT_SUCCESS )
         {
            PRINT_ERROR("Problems with building volume.\n");
            num_cut = 0;
            return CUBIT_FAILURE;
         }
         else
         {
            results_list += temp_new_bodies;
            num_cut++;
         }
      }
   }

     // Delete the BODY that was created to be used as a tool
   AcisQueryEngine::instance()->delete_ACIS_BODY(cutting_tool_ptr, CUBIT_TRUE) ;

   // Return CUBIT_SUCCESS as long as one body is webcut, earlier failures should
   // return CUBIT_FAILURE independently of these return statments
   if (num_cut>0)
   {
     return CUBIT_SUCCESS;
   }
   else
   {
     return CUBIT_FAILURE;
   }
}

//-------------------------------------------------------------------------
// Purpose       : creates a sheet body with the given curve loop
//                 uses the new sheet body to cut the body list.
//
// Special Notes :
//
// Creator       : Eric W. Nielsen
//
// Creation Date : 12/20/99
//-------------------------------------------------------------------------
CubitStatus
AcisModifyEngine::webcut_with_curve_loop( DLIList<BodySM*> &webcut_body_list,
                                          DLIList<Curve*> &refedge_list,
                                          DLIList<BodySM*>& results_list,
                                          bool imprint)
{
  Curve* ref_edge_ptr;
  int edge_count = refedge_list.size() ;
  EDGE** EDGEs = new EDGE* [edge_count] ;
  EDGE* EDGE_ptr = NULL;

  ENTITY_LIST acis_FACE_list ;
  int i;

  for( i=0; i<edge_count; i++ )
  {
     ref_edge_ptr = refedge_list.get_and_step();

     // Get the EDGE ptr from the RefEdge
     EDGE_ptr = AcisQueryEngine::get_EDGE( ref_edge_ptr );
     if( EDGE_ptr == NULL )
     {
        for( i=0; i<edge_count; i++ )
           api_delent( EDGEs[i] );
        delete [] EDGEs;
        PRINT_ERROR( "Curve is not an ACIS curve\n" );
        acis_FACE_list.clear();
        return CUBIT_FAILURE;
     }

     // copy each EDGE
     outcome result = api_edge( EDGE_ptr, EDGEs[i] );
     if (!result.ok())
     {
        for( i=0; i<edge_count; i++ )
           api_delent( EDGEs[i] );
        delete [] EDGEs;
        PRINT_ERROR( "Unable to copy curve; aborting\n" );
        AcisQueryEngine::instance()->ACIS_API_error (result);
        return CUBIT_FAILURE;
     }

  }

  BODY* wire_BODY = NULL ;

  // make a wire Body
  outcome result = api_make_ewire(edge_count, EDGEs, wire_BODY) ;
  if ( !result.ok() || wire_BODY == NULL )
  {
     PRINT_ERROR("In AcisModifyEngine::webcut_with_curve_loop\n"
        "       ACIS api_make_ewire function failed.\n\n");
     AcisQueryEngine::instance()->ACIS_API_error (result);
     for( i=0; i<edge_count; i++ )
        api_delent( EDGEs[i] );
     delete [] EDGEs;
     acis_FACE_list.clear();
     return CUBIT_FAILURE;
  }

  delete [] EDGEs;

  // Check to make sure the wire is closed
  result = api_closed_wire( wire_BODY );
  if( !result.ok() )
  {
    PRINT_ERROR( "ACIS reports that a closed loop was not found\n" );
    if( wire_BODY )
        api_delent( wire_BODY );
    acis_FACE_list.clear();
    return CUBIT_FAILURE;
  }

  // Cover the wire Body
  outcome result_1 = api_cover_wires ( wire_BODY,
     *(surface*)NULL_REF,
     acis_FACE_list );

  if ( !result_1.ok() || acis_FACE_list.count() == 0 ||
     wire_BODY == NULL )
  {
     PRINT_ERROR("In AcisModifyEngine::webcut_with_curve_loop\n"
        "       ACIS api_cover_wires function failed.\n\n");
     AcisQueryEngine::instance()->ACIS_API_error (result_1);
     acis_FACE_list.clear();
     if( wire_BODY )
        api_delent( wire_BODY );
     return CUBIT_FAILURE;
  }

  if (acis_FACE_list.count() != 1)
  {
     PRINT_ERROR("ACIS api_cover_wires function created "
        "more than one FACE.\n"
        "       At this time, we cannot deal with this.\n");
     acis_FACE_list.clear();
     if( wire_BODY )
        api_delent( wire_BODY );
     return CUBIT_FAILURE;
  }

  else
  {
     //make this a 2d body.
     outcome result_2 = api_body_to_2d( wire_BODY );
     if ( !result_2.ok() )
     {
        PRINT_ERROR("In AcisModifyEngine::webcut_with_curve_loop\n"
           "       ACIS api_body_to_2d function failed.\n\n");
        AcisQueryEngine::instance()->ACIS_API_error (result_2);
        acis_FACE_list.clear();
        if( wire_BODY )
           api_delent( wire_BODY );
        return CUBIT_FAILURE;
     }
     FACE *FACE_ptr = (FACE*) acis_FACE_list.next();
     if ( FACE_ptr == NULL || FACE_ptr->geometry() == NULL )
     {
        PRINT_ERROR("In AcisModifyEngine::webcut_with_curve_loop\n"
           "FACE or FACE's geometry was NULL.\n");
        if( wire_BODY )
           api_delent( wire_BODY );
        acis_FACE_list.clear();
        return CUBIT_FAILURE;
     }
  }

  BodySM* webcut_body;
  DLIList<BodySM*> temp_new_bodies;
  int num_cut = 0;
  for ( int ii = webcut_body_list.size(); ii > 0; ii-- )
  {
    webcut_body = webcut_body_list.remove();
    BODY *webcut_BODY = AcisQueryEngine::get_BODY(webcut_body );
      //now webcut with it.
    BODY *new_BODY_1, *new_BODY_2;
    CubitStatus status = webcut_with_sheet( webcut_BODY,
                                            wire_BODY,
                                            new_BODY_1,
                                            new_BODY_2,
                                            imprint );
    if ( status == CUBIT_SUCCESS )
    {
      DLIList<BODY*> new_BODY_list;
      new_BODY_list.append(new_BODY_1);
      new_BODY_list.append(new_BODY_2);
      temp_new_bodies.clean_out();

      CubitStatus cub_result = get_new_Body(webcut_body, webcut_BODY,
                                            new_BODY_list, temp_new_bodies,
                                            CUBIT_FALSE);
      if ( cub_result != CUBIT_SUCCESS )
      {
        PRINT_ERROR("Problems with building volume.\n");
        num_cut = 0;
        if( wire_BODY )
           api_delent( wire_BODY );
        return CUBIT_FAILURE;
      }
      else
      {
        results_list += temp_new_bodies;
        num_cut++;
      }
    }
  }

  if( wire_BODY )
     api_delent( wire_BODY );
  acis_FACE_list.clear();

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Sections the list of bodies (cuts in half).  This is
//                 when trying to see a cross section of the model.
//
// Special Notes : Rewritten - J.Kraftcheck - 2004-4-5
//
// Creator       : David White
//
// Creation Date : 12/10/97
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::section(  DLIList<BodySM*>& section_body_list,
                                        const CubitVector &point_1,
                                        const CubitVector &point_2,
                                        const CubitVector &point_3,
                                        DLIList<BodySM*>& new_body_list,
                                        bool keep_normal_side,
                                        bool keep_old,
                                        bool keep_both_sides )
{
  int i;
  bool delete_attribs =
     (GeometryModifyTool::instance()->get_new_ids() || keep_old);

   if (keep_both_sides == CUBIT_TRUE )
   {
      PRINT_ERROR("keeping both sides of section is not implemented.");
      return CUBIT_FAILURE;
   }

     // Compute the super bounding SPAbox
   SPAbox super_bounding_box;
   CubitStatus rv =
       AcisQueryEngine::instance()->create_super_acis_bounding_box (section_body_list,
                                                        super_bounding_box);
   if (rv != CUBIT_SUCCESS)
      return rv;

     // Construct a BODY from the given CubitVectors and the size of
     // the super bounding SPAbox.
   BODY *cutting_tool_ptr;
   cutting_tool_ptr = this->create_infinite_plane_cutting_tool(
       point_1,
       point_2,
       point_3,
       super_bounding_box );
   if (cutting_tool_ptr == NULL || cutting_tool_ptr->lump() == NULL )
   {
      PRINT_ERROR("failed in sectioning volumes.\n");
      return CUBIT_FAILURE;
   }

    // Copy BODYs
  i = section_body_list.size();
  DLIList<BODY*>  tool_BODY_list(i);
  DLIList<BODY*> blank_BODY_list(i);
  DLIList<BODY*> input_BODY_list(i);
  section_body_list.reset();
  tool_BODY_list.append( cutting_tool_ptr );
  for (i = section_body_list.size(); i > 0; i-- )
  {
    BodySM* body_sm = section_body_list.get_and_step();

    BODY* input_BODY = AcisQueryEngine::get_BODY( body_sm );
    if (!input_BODY) continue;

    //also make sure that the cutting plane and blank body intersect
    if ( !BODYs_interfering( input_BODY, cutting_tool_ptr) )
    {
      TopologyEntity *topo_entity = body_sm->topology_entity();
      RefEntity *ref_entity = CAST_TO(topo_entity, RefEntity);
      assert(ref_entity != 0);
      PRINT_WARNING("Volume %d does not interect sectioning plane\n", ref_entity->id() );
      section_body_list.remove( body_sm );
      continue;
    }

    input_BODY_list.append( input_BODY );

    BODY* blank_BODY = copy_BODY( input_BODY, delete_attribs );
    if (!blank_BODY) continue;
    blank_BODY_list.append( blank_BODY  );

    if (i < section_body_list.size()) // not first iteration
    {
      BODY* tool_BODY = copy_BODY( cutting_tool_ptr );
      if (tool_BODY)
        tool_BODY_list.append( tool_BODY );
    }
  }

  if( input_BODY_list.size() == 0 )
  {
    PRINT_WARNING("No geometry interects sectioning plane\n");
    return CUBIT_FAILURE;
  }

    // If there were any errors in the above loop, clean up and exit
  if ( input_BODY_list.size() < section_body_list.size() ||
       blank_BODY_list.size() < section_body_list.size() ||
        tool_BODY_list.size() < section_body_list.size() )
  {
     if (input_BODY_list.size() < section_body_list.size())
       PRINT_ERROR("Internal error: Invalid BodySM at %s:%d\n", __FILE__, __LINE__ );
     else
       PRINT_ERROR("Internal error: BODY copy failed at %s:%d\n", __FILE__, __LINE__ );
     while (tool_BODY_list.size())
       api_delent( input_BODY_list.pop() );
     while (blank_BODY_list.size())
       api_delent( input_BODY_list.pop() );
     return CUBIT_FAILURE;
  }

    // Deterine which type of Boolean to use.
  BOOL_TYPE bool_type;
  CubitBoolean nonreg = !GeometryModifyTool::instance()->boolean_regularize();
  if (keep_normal_side && !nonreg)
    bool_type = SUBTRACTION;
  else if (keep_normal_side && nonreg)
    bool_type = NONREG_SUBTRACTION;
  else if (!keep_normal_side && !nonreg)
    bool_type = INTERSECTION;
  else
    bool_type = NONREG_INTERSECTION;

    // Do Boolean for each blank body
  input_BODY_list.reset();
  blank_BODY_list.reset();
  tool_BODY_list.reset();
  for ( i = input_BODY_list.size(); i--; )
  {
    input_BODY_list.step();
    BODY* blank_BODY = blank_BODY_list.get_and_step();
    BODY*  tool_BODY =  tool_BODY_list.get_and_step();
    outcome result = api_boolean( tool_BODY, blank_BODY, bool_type );
    if (!result.ok())
    {
      PRINT_ERROR("Volume section failed.\n");
      blank_BODY_list.back();
      BODY *blank = blank_BODY_list.remove();
      AcisBridge *tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(blank);
      if( tb_ptr && keep_old == false )
      {
        BodyACIS *body_acis = CAST_TO(tb_ptr, BodyACIS );
        if( body_acis )
        {
          BodySM *body_sm = CAST_TO(body_acis, BodySM );
          TopologyEntity *topo_entity = body_sm->topology_entity();
          RefEntity *ref_entity = CAST_TO(topo_entity, RefEntity);
          section_body_list.remove( body_sm );
        }
      }

      api_delent( blank );
      tool_BODY_list.back();
      api_delent( tool_BODY_list.remove() );
      input_BODY_list.back();
      input_BODY_list.remove();
    }
  }

  if( input_BODY_list.size() == 0 )
    return CUBIT_FAILURE;

    // If a non-regularized operation, clean up internal surfaces
  if (nonreg)
  {
    DLIList<FACE*> FACE_list;
    input_BODY_list.reset();
    blank_BODY_list.reset();

    for (i = input_BODY_list.size(); i--; )
    {
      input_BODY_list.step();
      BODY* blank_BODY = blank_BODY_list.get_and_step();

      FACE_list.clean_out();
      AcisQueryEngine::instance()->get_FACEs( blank_BODY, FACE_list );

      FACE_list.reset();
      for (int j = FACE_list.size(); j--; )
      {
        FACE* FACE_ptr = FACE_list.get_and_step();
        assert( FACE_ptr != NULL );

        if (FACE_ptr->sides() == DOUBLE_SIDED &&
            FACE_ptr->cont() == BOTH_INSIDE)
        {
          PRINT_INFO("Unhooking and deleting internal surface.\n");
          BODY* FACE_BODY_ptr;
          outcome result = api_unhook_face( FACE_ptr, FACE_BODY_ptr );
          if (!result.ok())
          {
            AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
            PRINT_ERROR( " Face unhooking during rebuild of volume didn't work\n" );
          }
          else
          {
            api_delent( FACE_BODY_ptr );
          }
        }
      }
    }
  }

    // Create/update AcisBridges
  if( keep_old == true )
  {
    int k;
    for( k=blank_BODY_list.size(); k--; )
    {
      BodySM *this_bodysm = AcisQueryEngine::instance()->
            populate_topology_bridges( blank_BODY_list.get_and_step() );
      new_body_list.append( this_bodysm );
    }
  }
  else
    get_new_Body( section_body_list,
                  input_BODY_list,
                  blank_BODY_list,
                  new_body_list,
                  keep_old );

  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : webcuts a list of bodies using a plane.
//                 The newly created bodies are merged and imprinted
//                 depending on the respective flags.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::webcut( DLIList<BodySM*>& webcut_body_list,
                                      const CubitVector &vecVertex1,
                                      const CubitVector &vecVertex2,
                                      const CubitVector &vecVertex3,
                                      DLIList<BodySM*>& results_list,
                                      bool imprint ) const
{
   //int webcut_debug_flag = 18;
   //int webcut_debug_flag_on = DEBUG_FLAG(webcut_debug_flag);

   BODY* cutting_tool_ptr = NULL;

     // To be able to create a Cutting Tool of the appropriate ("infinite")
     // size, first generate the bounding SPAbox that surrounds all the BODYs
     // being webcut. The procedure, GeometryModifyTool::create_infinite_plane_cutting
     // _tool, will use the size of this "super" bounding SPAbox to determine an
     // appropriate "infinite-enough" size for the Cutting Tool.

     // Compute the super bounding SPAbox
   SPAbox super_bounding_box;
   CubitStatus rv =
       AcisQueryEngine::instance()->create_super_acis_bounding_box ( webcut_body_list, super_bounding_box );
   if (rv != CUBIT_SUCCESS)
      return CUBIT_FAILURE;

     // Construct a BODY from the given CubitVectors and the size of
     // the super bounding SPAbox.
   cutting_tool_ptr = this-> create_infinite_plane_cutting_tool(
       vecVertex1,
       vecVertex2,
       vecVertex3,
       super_bounding_box );

   if (cutting_tool_ptr == NULL)
   {
        //webcut_failed(refbody_cleanup_list, BODY_cleanup_list);
      PRINT_ERROR("In AcisModifyEngine::webcut(,Vector,,,),%d"
                  "       Cannot create infinite plane cutting tool using \n"
                  "       the vectors. \n", __LINE__ );
      return CUBIT_FAILURE;
   }

   else
   {
        // Use the BODY to perform webcut
      CubitStatus stat =
          this->webcut(webcut_body_list, cutting_tool_ptr, results_list, imprint) ;

        // Delete the BODY that was created to be used as a tool
      AcisQueryEngine::instance()->delete_ACIS_BODY(cutting_tool_ptr, CUBIT_TRUE) ;
      cutting_tool_ptr = NULL;

      return stat ;
   }
}

//-------------------------------------------------------------------------
// Purpose       : webcuts a list of bodies using another Body as the tool.
//                 The newly created bodies are merged and imprinted
//                 depending on the respective flags.
//
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 12/23/96
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::webcut( DLIList<BodySM*>& webcut_body_list,
                                      BodySM const* tool_body,
                                      DLIList<BodySM*>& results_list,
                                      bool imprint ) const
{
     // check the tool body if it is made of AcisModifyEngine
   if( !is_modify_engine(tool_body))
   {
      PRINT_ERROR("In AcisModifyEngine::webcut() %d, \n"
                  "   The tool volume is not made of AGE itself.\n", __LINE__);
      return CUBIT_FAILURE;
   }


   BODY* tool_BODY = NULL;

     //get the ACIS BODY of tool_body
   tool_BODY = AcisQueryEngine::get_BODY( const_cast<BodySM*>(tool_body) );
   if( tool_BODY == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::webcut() %d, \n"
                  "       There is no underlying tool volume.\n",__LINE__);
      return CUBIT_FAILURE;
   }

   return this->webcut( webcut_body_list, tool_BODY, results_list, imprint );
}

//-------------------------------------------------------------------------
// Purpose       : webcuts a list of bodies using another Body as the tool.
//                 The newly created bodies are merged and imprinted
//                 depending on the respective flags.
// Special Notes :
//
// Creator       : Jihong Ma
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::webcut( DLIList<BodySM*>& webcut_body_list,
                                      BODY* tool_BODY,
                                      DLIList<BodySM*>& results_list,
                                      bool imprint ) const
{
  CubitBoolean delete_bodies = (GeometryModifyTool::instance()->get_new_ids() ?
                                CUBIT_FALSE : CUBIT_TRUE);

    // Number of bodies that were webcut
   int count = 0 ;

   BodySM* Body_to_be_webcut = NULL;
   BODY* BODY_to_be_webcut = NULL;

     // List to maintain thew newly created BODYs
   DLIList<BODY*> new_webcut_BODY_list;

     // List of old bodies that did get webcut
   DLIList<BodySM*> old_body_list;
   DLIList<BODY*> old_BODY_list;

     // Variables to get the results of webcut from another routine.
   BODY* BODY_result1 = NULL;
   BODY* BODY_result2 = NULL;

     // For each Body in the input list, make sure that it is associated
     // with AGE before the webcutting is done.
   CubitStatus status = CUBIT_FAILURE;
   webcut_body_list.reset() ;
   int i;

     // step 1: bodies taken from webcut_body_list and webcut; if
     //    successful, the old Body and BODY are put on old_Body(BODY)_list
     //    and pieces put on new_webcut_BODY_list
     //   => 2 important lists:
     //      a) new_webcut_BODY_list (new pieces from bodies in old_...)
     //      b) results_list - bodies that weren't cut by the tool
   for( i = 0 ; i < webcut_body_list.size() ; i++)
   {
      Body_to_be_webcut = webcut_body_list.get_and_step();
      if( !is_modify_engine(Body_to_be_webcut))
      {
         PRINT_ERROR("In AcisModifyEngine::webcut(), \n"
                     "       Volumes are associated with a different \n"
                     "       geometric modeling engine.\n" );
      }

        //otherwise, get the BODY of the Body_to_be_webcut
      BODY_to_be_webcut = AcisQueryEngine::get_BODY( Body_to_be_webcut );
      if( BODY_to_be_webcut == NULL )
      {
         PRINT_ERROR("In AcisModifyEngine::webcut() %d, \n"
                     "       There is no underlying volume.\n",
                     __LINE__);
         assert(BODY_to_be_webcut != NULL) ;
      }

        //now can do a webcut, using a BODY tool, resulting two BODYs
      status = webcut_BODY( BODY_to_be_webcut, tool_BODY,
                            BODY_result1, BODY_result2);

      if (status == CUBIT_SUCCESS )
      {
           // Add the 2 new BODYs to the list of webcut BODYs
         new_webcut_BODY_list.append(BODY_result1);
         new_webcut_BODY_list.append(BODY_result2);
         old_body_list.append(Body_to_be_webcut);
         old_BODY_list.append(BODY_to_be_webcut);

           // Since we have successfully webcut this BODY, delete it.
           // Delete the corresponding Body too. (only do this if we're not
           // using persistent objects)

         if (delete_bodies == CUBIT_FALSE) {
           AcisQueryEngine::instance()->delete_ACIS_BODY(BODY_to_be_webcut) ;
           //AcisQueryEngine::instance()->delete_solid_model_entities(Body_to_be_webcut);
         }

           // Increment the counter
         count++ ;
      }

      else
      {
           //Since the body wasn't cut, it is a result...
         results_list.append_unique( Body_to_be_webcut );
      }
   }

     // A list to hold the newly created Bodies due to imprint operations
     // on the remaining Bodys in the Model. The new Bodys replace the old
     // Bodys if any of the new BODYs created by webcutting imprint
     // on the old Bodys.
   DLIList<BodySM*>  new_webcut_body_list;

   DLIList<BODY*> just_webcut_list;
     // A list to get the BODY's created that were webcut from the imprint
     // operation.


     // Perform the imprinting operations, if required. Note that this
     // operation needs to be performed *after* the model has been cleared of
     // Bodies that were already webcut. Also, note that it returns a list
     // of all the new Bodies created during the imprint operation. This
     // list is then enhanced by adding any new Bodies created due to the
     // actual cutting of an old Body, before it is sent into the webcut_merge
     // procedure. This is the list of Bodies that is "merged".
   if (new_webcut_BODY_list.size() > 0)
   {
     DLIList<BODY*> old_model_BODY_list;
     DLIList<BodySM*> old_model_Body_list;

      if ( imprint == CUBIT_TRUE )
      {
           // Work-in-progress message
         PRINT_INFO ("      Performing the imprint operations...\n");

           // Imprint away :)
         webcut_imprint(tool_BODY, old_body_list, new_webcut_BODY_list,
                        just_webcut_list, results_list,
                        old_model_Body_list, old_model_BODY_list) ;
      }
      if ( just_webcut_list.size() == 0 )
          just_webcut_list = new_webcut_BODY_list;

        // Work-in-progress message
      PRINT_INFO ("      Creating the new (webcut) volumes...\n");

        // Now create the new Bodies from the ACIS BODYs that were created
        // by webcutting the original BODYs using the Cutting Tool.
      BodySM* new_body_ptr = NULL;
      BODY* webcut_BODY_ptr = NULL;
      new_webcut_BODY_list.reset();
      old_body_list.reset();
      old_BODY_list.reset();
      old_model_Body_list.reset();
      old_model_BODY_list.reset();

      for ( i = 0 ; i < new_webcut_BODY_list.size(); i++ )
      {
        DLIList<BODY*> temp_BODY_list;
        DLIList<BodySM*> temp_body_list;
        BodySM *old_Body;
        BODY *old_BODY;

        webcut_BODY_ptr = new_webcut_BODY_list.get_and_step();
        temp_BODY_list.append(webcut_BODY_ptr);
          // in call to get_new_Body, don't delete old bodies if
          // they were successfully webcut (need to do it that way to
          // preserve name propagation as before); otherwise, delete the
          // imprinted bodies
        if (imprint == CUBIT_FALSE || just_webcut_list.move_to(webcut_BODY_ptr)) {
          webcut_BODY_ptr = new_webcut_BODY_list.get_and_step();
          temp_BODY_list.append(webcut_BODY_ptr);
          i++;
          old_Body = old_body_list.get_and_step();
          old_BODY = old_BODY_list.get_and_step();
          status = get_new_Body( old_Body, old_BODY, temp_BODY_list, temp_body_list,
                                 false, delete_bodies );
        }
        else {
          old_Body = old_model_Body_list.get_and_step();
          old_BODY = old_model_BODY_list.get_and_step();
          BodySM* new_Body = get_new_Body( old_Body, old_BODY, webcut_BODY_ptr,
                                           false, true, true );
          status = CUBIT_SUCCESS;
          if (new_Body)
            temp_body_list.append(new_Body);
        }

        if (status == CUBIT_FAILURE) {
          PRINT_ERROR("Failed to build new volume in AGE::webcut.\n");
        }

          // if we didn't get any new bodies back from get_new_Body,
          // that means that an imprinted body wasn't really imprinted
        if (temp_body_list.size() > 0) {
          temp_BODY_list.reset();
          temp_body_list.reset();
          for (int j = 0; j < temp_BODY_list.size(); j++) {
            webcut_BODY_ptr = temp_BODY_list.get_and_step();
            new_body_ptr = temp_body_list.get_and_step();

            new_webcut_body_list.append(new_body_ptr);
            results_list.append_unique(new_body_ptr);
          }
        }
        else
          new_webcut_body_list.append(old_Body);

      }
   }
   // Return CUBIT_SUCCESS as long as one body is webcut, earlier failures should
   // return CUBIT_FAILURE independently of these return statments
   if (count>0)
   {
     return CUBIT_SUCCESS;
   }
   else
   {
     return CUBIT_FAILURE;
   }
}


//-------------------------------------------------------------------------
// Purpose       : This function creates a Point, given coordinates.  The
//                 particular type of Point object that is created depends on
//                 the specific GeometryQueryEngine.  For example, this
//                 engine, AcisModifyEngine, creates a PointACIS and
//                 returns it.
//
// Special Notes :
// In this routine, a floating ACIS VERTEX is created. This VERTEX
// is not a part of any BODY and, hence, does not take part in
// any boolean operations.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 03/05/97
//-------------------------------------------------------------------------
Point* AcisModifyEngine::make_Point(CubitVector const& point) const
{
     // Create a new ACIS VERTEX -- it is not attached to a BODY
   VERTEX* VERTEX_ptr = this->make_VERTEX( point );
   if ( VERTEX_ptr == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_Point\n"
                  "       Cannot make Point object.\n");
      return (Point *)NULL;
   }

     // Create a new PointACIS object
   return AcisQueryEngine::instance()->populate_topology_bridges( VERTEX_ptr );

}



// This function creates a curve given an existing curve.
Curve* AcisModifyEngine::make_Curve(Curve *curve_ptr) const
{
  EDGE *EDGE_ptr, *new_EDGE_ptr = NULL;
  CurveACIS *curve_ACIS = CAST_TO(curve_ptr, CurveACIS);
  if (!curve_ACIS)
  {
     PRINT_ERROR("Cannot create an ACIS curve from the given curve.\n"
                 "Possible incompatible geometry engines.\n");
     return (Curve *)NULL;
  }
  EDGE_ptr = curve_ACIS->get_EDGE_ptr();

  outcome result = api_edge( EDGE_ptr, new_EDGE_ptr );

  if ( !result.ok() || new_EDGE_ptr == NULL )
  {
    PRINT_ERROR("In AcisModifyEngine::make_Curve\n"
                "       Cannot make Curve object.\n");
    return (Curve *)NULL;
  }

    // Create a new CurveACIS object
    // Use the new EDGE to create a new CurveACIS
  return AcisQueryEngine::instance()->populate_topology_bridges(new_EDGE_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : This function creates a Curve of type, curve_type,
//                 given the end points of the curve and an intermediate
//                 point.  The particular type of Curve object that is
//                 created depends on the specific modeling engine.  For
//                 example, if the engine is AcisModifyEngine, then a
//                 CurveACIS is created and returned.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------

Curve* AcisModifyEngine::make_Curve( GeometryType curve_type,
                                       Point const* point1_ptr,
                                       Point const* point2_ptr,
                                       CubitVector const* intermediate_point_ptr,
                                       CubitSense sense) const
{
     // Extract the VERTEXes from the Points and construct a SPAposition
     // from the cubit SPAvector.

     // Cast the Points to  PointACIS
   PointACIS* pointACIS_ptr1 = CAST_TO(const_cast<Point*>(point1_ptr), PointACIS) ;
   PointACIS* pointACIS_ptr2 = CAST_TO(const_cast<Point*>(point2_ptr), PointACIS) ;

     // Make sure that we don't get NULL pointers after the cast down.
   assert ( pointACIS_ptr1 != NULL && pointACIS_ptr2 != NULL ) ;

     // Get a reference to the VERTEX list of the PointACIS
   VERTEX* VERTEX_ptr1 = pointACIS_ptr1->get_VERTEX_ptr();
   VERTEX* VERTEX_ptr2 = pointACIS_ptr2->get_VERTEX_ptr();

   SPAposition* pos_ptr = NULL ;
   SPAposition pos ;

   if ( curve_type != STRAIGHT_CURVE_TYPE )
   {
      assert ( intermediate_point_ptr != NULL) ;

        // Construct a SPAposition from the given intermediate point.
      pos = SPAposition (intermediate_point_ptr->x(),
                      intermediate_point_ptr->y(),
                      intermediate_point_ptr->z()) ;
      pos_ptr = &pos ;
   }

   logical acis_sense = TRUE;
   if ( sense == CUBIT_REVERSED )
   {
      acis_sense = FALSE;
   }

     // Use the VERTEXes and the SPAposition to make a new EDGE
   EDGE* EDGE_ptr = this->make_EDGE(curve_type,
                                    VERTEX_ptr1,
                                    VERTEX_ptr2,
                                    pos_ptr,
                                    acis_sense == TRUE) ;

   if ( EDGE_ptr == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_Curve\n"
                  "       Cannot make Curve object.\n");
      return (Curve *)NULL;
   }

     // Use the new EDGE to create a new CurveACIS
   return AcisQueryEngine::instance()->populate_topology_bridges(EDGE_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : This function creates a Curve of type, curve_type,
//                 given the end points of the curve and an intermediate
//                 point.  The particular type of Curve object that is
//                 created depends on the specific modeling engine.  For
//                 example, if the engine is AcisModifyEngine, then a
//                 CurveACIS is created and returned.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 05/16/97
//-------------------------------------------------------------------------
Curve* AcisModifyEngine::make_Curve( GeometryType curve_type,
                                       Point const* point1_ptr,
                                       Point const* point2_ptr,
                                       DLIList<CubitVector*>& vector_list,
                                       Surface* ref_face_ptr ) const
{
     // We make splines here :-)
   if ( curve_type != SPLINE_CURVE_TYPE )
   {
      PRINT_ERROR("In AcisModifyEngine::make_Curve\n"
                  "       Can only make spline surves in this routine.\n");
      return (Curve *)NULL;
   }

     // Extract the VERTEXes from the Points.

     // Cast the Points to  PointACIS
   PointACIS* pointACIS_ptr1 = CAST_TO(const_cast<Point*>(point1_ptr), PointACIS) ;
   PointACIS* pointACIS_ptr2 = CAST_TO(const_cast<Point*>(point2_ptr), PointACIS) ;

     // Make sure that we don't get NULL pointers after the cast down.
   assert ( pointACIS_ptr1 != NULL && pointACIS_ptr2 != NULL ) ;

     // Get a reference to the VERTEX list of the PointACIS
   VERTEX* VERTEX_ptr1 = pointACIS_ptr1->get_VERTEX_ptr();
   VERTEX* VERTEX_ptr2 = pointACIS_ptr2->get_VERTEX_ptr();

     // Convert the list of CubitVectors into an array of ACIS position's.
     // MJP NOTE:
     // Add the positions of the 2 end VERTEX'es into this list.
   int number_points = vector_list.size();
   vector_list.reset();
   CubitVector* vector_ptr = NULL, on_surf;
   SPAposition* pos_array = new SPAposition [number_points];

     // Do the vector_list of points
   for (int i = 0; i < number_points; i++)
   {
      vector_ptr = vector_list.get_and_step();

        // If the input RefFace pointer is not NULL, then move this point
        // to the surface of the RefFace
      if (ref_face_ptr != NULL)
      {
         ref_face_ptr->closest_point(*vector_ptr, &on_surf);
         *vector_ptr = on_surf;
         //ref_face_ptr->move_to_surface(*vector_ptr);
      }

        // Set the coordinates of the SPAposition object
      pos_array[i].set_x(vector_ptr->x());
      pos_array[i].set_y(vector_ptr->y());
      pos_array[i].set_z(vector_ptr->z());
   }

     // Use the VERTEXes and the SPAposition array to make a new EDGE
   EDGE* EDGE_ptr = this->make_spline_EDGE( VERTEX_ptr1,
                                            VERTEX_ptr2,
                                            pos_array,
                                            number_points );
     // Free memory
   delete [] pos_array;

   if ( EDGE_ptr == NULL )
   {
      PRINT_ERROR("In AcisModifyEngine::make_Curve\n"
                  "       Cannot make Curve object.\n");
      return (Curve *)NULL;
   }

   return AcisQueryEngine::instance()->populate_topology_bridges(EDGE_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : This function creats a curve that is on a
//                 given surface.  The third point is used
//                 for curves that could be periodic to dertermine
//                 the correct direction.
//
//
// Creator       : David White
//
// Creation Date : 10/21/97
//-------------------------------------------------------------------------
Curve* AcisModifyEngine::make_Curve( Point const* point1_ptr,
                                     Point const* point2_ptr,
                                     Surface* ref_face_ptr,
                                     const CubitVector *third_point) const

{
   CubitVector point_1 = point1_ptr->coordinates();
   CubitVector point_2 = point2_ptr->coordinates();
   CubitVector surf_norm;
   ref_face_ptr->closest_point( point_1, NULL, &surf_norm );

   CubitVector plane_norm;

   CubitVector edge0 = point_2 - point_1;
   CubitVector edge1 = surf_norm;
   plane_norm = edge0*edge1;
   if ( plane_norm.length() < CUBIT_RESABS )
   {
        //try using the normal at point 2.
      ref_face_ptr->closest_point(point_2, NULL, &surf_norm);
      edge1 = surf_norm;
      edge0 = point_1 - point_2;
      plane_norm = edge0*edge1;
      if ( plane_norm.length() < CUBIT_RESABS )
      {
         if ( third_point == NULL )
         {
            PRINT_ERROR("Can't create curve on surface.\n"
                        "\tVertices are on opposite sides of surface.\n"
                        "\tUnable to infer a plane that will cut the surface.\n"
                        "\tTry specifying a 3rd vertex in the new curve.\n");
            return (Curve*)NULL;
         }
         else
         {
             // Try getting a normal from the third point
           edge0 = point_1 - *third_point;
           plane_norm = edge0*edge1;
           if ( plane_norm.length() < CUBIT_RESABS )
           {
             PRINT_ERROR("Can't create curve on surface.\n"
                         "\tVertices are on opposite sides of surface\n"
                         "\tand all three points are colinear.\n"
                         "\tTry specifying a non-colinear 3rd vertex.\n");
             return (Curve*)NULL;
           }
         }
      }
   }
   plane_norm.normalize();
     // Extract the VERTEXes from the Points.

     // Cast the Points to  PointACIS
   PointACIS* pointACIS_ptr1 = CAST_TO(const_cast<Point*>(point1_ptr), PointACIS) ;
   PointACIS* pointACIS_ptr2 = CAST_TO(const_cast<Point*>(point2_ptr), PointACIS) ;

     // Make sure that we don't get NULL pointers after the cast down.
   assert ( pointACIS_ptr1 != NULL && pointACIS_ptr2 != NULL ) ;

     // Get a reference to the VERTEX list of the PointACIS
   VERTEX* VERTEX_ptr1 = pointACIS_ptr1->get_VERTEX_ptr();
   VERTEX* VERTEX_ptr2 = pointACIS_ptr2->get_VERTEX_ptr();

   FACE *FACE_ptr = CAST_TO(ref_face_ptr, SurfaceACIS )->get_FACE_ptr();
     //now using the points, vertices, FACE and plane normal,
     //creat an edge.
   EDGE* EDGE_ptr = this->make_surface_EDGE( VERTEX_ptr1,
                                             VERTEX_ptr2,
                                             FACE_ptr,
                                             plane_norm,
                                             third_point);
     // Free memory
   if ( EDGE_ptr == NULL )
   {
     PRINT_ERROR("In AcisModifyEngine::make_Curve\n"
                 "       Cannot make Curve object.\n");
     return (Curve *)NULL;
   }

   return AcisQueryEngine::instance()->populate_topology_bridges(EDGE_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : Create a surface given a surface. From another surface
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 10/6/97
//-------------------------------------------------------------------------

Surface* AcisModifyEngine::make_Surface( Surface *old_surface_ptr,
                                           bool extended_from) const
{
     //Get the surface ACIS.
   SurfaceACIS *surf_ACIS = CAST_TO(old_surface_ptr, SurfaceACIS );
   if (surf_ACIS == NULL)
   {
      PRINT_ERROR("surface is not created from acis.\n");
      return (Surface*) NULL;
   }

     //Get the acis face.
   FACE *FACE_ptr = surf_ACIS->get_FACE_ptr();

     //use this to make a FACE.
   if (FACE_ptr == NULL )
   {
     PRINT_ERROR("surface is not created from acis.\n");
     return (Surface*) NULL;
   }
   FACE *new_FACE_ptr = NULL;

   new_FACE_ptr = make_FACE( FACE_ptr, extended_from );
   if (new_FACE_ptr == NULL)
   {
     PRINT_ERROR("In AcisModifyEngine::make_Surface\n"
                 "       Cannot make Surface object.\n");
     return (Surface *)NULL;
   }

   //populate the bridges from the body
   BODY *body_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY(new_FACE_ptr);
   BodySM *body_sm = AcisQueryEngine::instance()->populate_topology_bridges( body_ptr );
   DLIList<Surface*> surfs;
   body_sm->surfaces( surfs );
   Surface *surface = surfs.get();

   return surface;

}

//-------------------------------------------------------------------------
// Purpose       : This function creates a Surface of type,
//                 surface_type, given the list of curves.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------
Surface* AcisModifyEngine::make_Surface( GeometryType surface_type,
                                           DLIList<Curve*>& curve_list,
                                           Surface *old_surface_ptr,
                                           bool check_edges ) const
{
  DLIList<RefEdge*> copied_ref_edges;
  int i;
  if( check_edges )
  {
    DLIList<RefEdge*> edge_list;
    for( i=curve_list.size(); i--; )
    {
      RefEdge *tmp_ref_edge = NULL;
      tmp_ref_edge = CAST_TO( curve_list.get_and_step()->topology_entity(), RefEdge );
      if( tmp_ref_edge )
        edge_list.append( tmp_ref_edge );
    }

    //duplicate curves possibly

    //Modified by J.K. to avoid core dump if vertices are merged.
    //If vertices belong to edges other than those we are creating
    //the surface from, copy the corresponding edges.  This will
    //also pick up the case where edges already belong to a face.
    DLIList<RefEdge*> vtx_edges;
    DLIList<RefVertex*> vtx_list;
    for ( i = edge_list.size(); i > 0; i-- )
    {
      RefEdge *ref_edge = edge_list.get();

      CubitBoolean other_edge = CUBIT_FALSE;
      ref_edge->ref_vertices(vtx_list);
      while( vtx_list.size() )
      {
        vtx_edges.clean_out();
        vtx_list.pop()->ref_edges( vtx_edges );
        while( vtx_edges.size() )
        {
          if( ! edge_list.is_in_list( vtx_edges.pop() ) )
          {
            other_edge = CUBIT_TRUE;
            vtx_list.clean_out();
            break;
          }
        }
      }

      if ( other_edge || (ref_edge->get_parents() > 0) )
      {
        RefEdge *replacement_edge = GeometryModifyTool::instance()->make_RefEdge( ref_edge );
        if (!replacement_edge)
        {
          PRINT_WARNING("Creation of Surface Unsuccessful\n");
          return (Surface *)NULL;
        }
        edge_list.change_to( replacement_edge );
        copied_ref_edges.append( replacement_edge );
      }
      edge_list.step();
    }

    curve_list.clean_out();
    for( i=edge_list.size(); i--; )
    {
      Curve *tmp_curve = NULL;
      tmp_curve =
      CAST_TO( edge_list.get_and_step()->bridge_manager()->topology_bridge(), Curve);
      if( tmp_curve )
        curve_list.append( tmp_curve );
    }
  }

   DLIList<EDGE*> EDGE_list ;
     // List of EDGEs to use to create the FACE from which we can make
     // a SurfaceACIS.

     // Extract the EDGEs from the given Curves
   Curve const* curve_ptr = NULL ;
   CurveACIS* curveACIS_ptr = NULL ;
   EDGE* EDGE_ptr = NULL ;

   curve_list.reset() ;


   for ( i = 0 ; i < curve_list.size() ; i++ )
   {
      curve_ptr = curve_list.get_and_step() ;
      curveACIS_ptr = CAST_TO(const_cast<Curve*>(curve_ptr), CurveACIS) ;

      // Make sure that we don't get NULL pointers after the cast down.
      if ( curveACIS_ptr == NULL )
      {
         PRINT_ERROR("In AcisModifyEngine::make_Surface\n"
                     "       Got a NULL pointer to CurveACIS\n") ;
         assert ( curveACIS_ptr != NULL ) ;
      }

        // Get the EDGE
      EDGE_ptr = curveACIS_ptr->get_EDGE_ptr() ;

        // Add the EDGE to the list of EDGEs to be used to make a FACE
      EDGE_list.append(EDGE_ptr) ;
   }

     // Use the EDGEs to make a FACE
   FACE* FACE_ptr = this->make_FACE(surface_type, EDGE_list,
                                    old_surface_ptr) ;

   if (FACE_ptr == NULL)
   {
     //delete all the copied ref edges...copied just for creating this surface
     for(i=copied_ref_edges.size(); i--; )
       GeometryQueryTool::instance()->delete_RefEdge( copied_ref_edges.get_and_step() );

     PRINT_ERROR("In AcisModifyEngine::make_Surface\n"
                  "       Cannot make Surface object.\n");
     return (Surface *)NULL;
   }

     // make the topology bridges for the face
   Surface *surface = AcisQueryEngine::instance()->populate_topology_bridges(FACE_ptr);
   return surface;

/*

     // now generate the body
   DLIList<BODY*> new_BODIES;
   DLIList<Body*> new_bodies;
   new_BODIES.append(AcisQueryEngine::instance()->get_BODY_of_ENTITY(FACE_ptr));
   CubitStatus success = get_new_Body(old_entities,
                                      new_BODIES, new_bodies, CUBIT_FALSE);

   AcisBridge *tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(FACE_ptr);
   return CAST_TO(tb_ptr, Surface);
*/
}

//-------------------------------------------------------------------------
// Purpose       : This function creates a Lump of type, lump_type,
//                 given the list of surfaces.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------

Lump* AcisModifyEngine::make_Lump( DLIList<Surface*>& surface_list ) const
{
     // Extract the FACEs of the surface list. Use the FACEs to
     // create a SHELL.
   Surface const* surface_ptr = NULL ;
   SurfaceACIS* surfaceACIS_ptr = NULL ;
   FACE* first_FACE_ptr = NULL ;
   FACE* current_FACE_ptr = NULL ;
   FACE* previous_FACE_ptr = NULL ;

   surface_list.reset() ;
   for ( int i = 0 ; i < surface_list.size() ; i++ )
   {
      surface_ptr = surface_list.get_and_step() ;
      surfaceACIS_ptr = CAST_TO(const_cast<Surface*>(surface_ptr), SurfaceACIS) ;

        // Make sure that we have a valid SurfaceACIS pointer
      if ( surfaceACIS_ptr == NULL )
      {
         PRINT_ERROR("In AcisModifyEngine::make_Lump\n"
                     "       Got a NULL pointer to SurfaceACIS\n") ;
         assert ( surfaceACIS_ptr != NULL ) ;
      }

      current_FACE_ptr = surfaceACIS_ptr->get_FACE_ptr();

        // Connect up the previous FACE and the current FACE if
        // the current FACE is not the first FACE.
      if ( previous_FACE_ptr != NULL )
      {
         previous_FACE_ptr->set_next(current_FACE_ptr) ;
      }
      else
      {
         first_FACE_ptr = current_FACE_ptr ;
      }

      previous_FACE_ptr = current_FACE_ptr ;
   }

     // Create a new SHELL if the first FACE is not attached to
     // a SHELL.
   SHELL* SHELL_ptr = first_FACE_ptr->shell() ;
   if ( SHELL_ptr == NULL )
   {
      SHELL_ptr = new SHELL(first_FACE_ptr, NULL, NULL) ;
      if ( SHELL_ptr == NULL )
      {
         PRINT_ERROR("In AcisModifyEngine::make_Lump\n"
                     "       ACIS cannot make SHELL object.\n"
                     "       Cannot make Lump object.\n");
         return (Lump *)NULL;
      }
   }

     // Create a new LUMP and attach the SHELL to it.
   LUMP* LUMP_ptr = SHELL_ptr->lump();
   if (LUMP_ptr == NULL)
   {
      LUMP_ptr = new LUMP(SHELL_ptr, NULL) ;
      if ( LUMP_ptr == NULL )
      {
         PRINT_ERROR("In AcisModifyEngine::make_Lump\n"
                     "       ACIS cannot make LUMP object.\n"
                     "       Cannot make Lump object.\n");
         return (Lump *)NULL;
      }
   }

     // Create a LumpACIS using the LUMP_ptr
   LumpACIS* lumpACIS_ptr = new LumpACIS(LUMP_ptr) ;
   assert( lumpACIS_ptr != NULL );

   AcisQueryEngine::instance()->populate_topology_bridges(LUMP_ptr);
/*
   AcisBridge *tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(LUMP_ptr);

   Lump *this_lump = CAST_TO(tb_ptr, Lump);
   ShellSM *shellsm = this_lump->shellsm();

     // Use the SHELL to build VGI Shell and CoFace, using 0 for color
   Shell* Shell_ptr = GeometryQueryTool::instance()->make_Shell(shellsm) ;

   if (Shell_ptr == NULL)
   {
      PRINT_ERROR("In AcisModifyEngine::make_Lump\n"
                  "       Could not build a Shell for unknown reasons.\n");
      return (Lump *)NULL;
   }
*/
     // Return the new LumpACIS
   return lumpACIS_ptr ;
}

BodySM* AcisModifyEngine::make_BodySM( Surface *surface_ptr) const
{
    //- given a Surface, make a BodySM
    // make sure this Surface is the right type
  SurfaceACIS *sa_ptr = CAST_TO(surface_ptr, SurfaceACIS);
  assert(sa_ptr != 0);

    // now, get the acis BODY and populate the topology bridges
  FACE *face = sa_ptr->get_FACE_ptr();
  BODY *body_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY(face);

  return AcisQueryEngine::instance()->populate_topology_bridges(body_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : This function creates a BodySM given a list of Lumps.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 03/29/97
//-------------------------------------------------------------------------

BodySM* AcisModifyEngine::make_BodySM( DLIList<Lump*>& lump_list ) const
{
     // Extract the LUMPs of the lump list. Use the LUMPs to
     // create a BODY.
   Lump const* lump_ptr = NULL ;
   LumpACIS* lumpACIS_ptr = NULL ;
   LUMP* first_LUMP_ptr = NULL ;
   LUMP* current_LUMP_ptr = NULL ;
   LUMP* previous_LUMP_ptr = NULL ;

   lump_list.reset() ;
   for ( int i = 0 ; i < lump_list.size() ; i++ )
   {
      lump_ptr = lump_list.get_and_step() ;
      lumpACIS_ptr = CAST_TO(const_cast<Lump*>(lump_ptr), LumpACIS) ;

        // Make sure that we have a valid LumpACIS pointer
      if ( lumpACIS_ptr == NULL )
      {
         PRINT_ERROR("In AcisModifyEngine::make_Lump\n"
                     "       Input Lump is not a LumpACIS\n") ;
         assert ( lumpACIS_ptr != NULL ) ;
      }

      current_LUMP_ptr = lumpACIS_ptr->get_LUMP_ptr();

        // Connect up the previous LUMP and the current LUMP if
        // the current LUMP is not the first LUMP.
      if ( previous_LUMP_ptr != NULL )
      {
         previous_LUMP_ptr->set_next(current_LUMP_ptr) ;
      }
      else
      {
         first_LUMP_ptr = current_LUMP_ptr ;
      }

      previous_LUMP_ptr = current_LUMP_ptr ;
   }

     // Create a new BODY if the first LUMP is not attached to
     // a BODY.
   BODY* BODY_ptr = first_LUMP_ptr->body() ;
   if ( BODY_ptr == NULL )
   {
      BODY_ptr = new BODY(first_LUMP_ptr) ;
      if ( BODY_ptr == NULL )
      {
         PRINT_ERROR("In AcisModifyEngine::make_BodySM\n"
                     "       ACIS cannot make BODY object.\n"
                     "       Cannot make BodySM object.\n");
         return (BodySM *)NULL;
      }
   }

     // Use the BODY to build VGI Body and CoVolume(s)
   BodySM *this_bodysm = AcisQueryEngine::instance()->populate_topology_bridges(BODY_ptr);

     /* Commenting out following section - we shouldn't be making
        Body's yet!!!

   Body* body_ptr = GeometryModifyTool::instance()->make_Body(this_bodysm);
   assert( body_ptr != NULL );

     // First get the OSME associated with this Body
   BodySM* OSME_ptr =
       body_ptr->get_body_sm_ptr() ;

     // Cast the OSME to BodyACIS
   BodyACIS* bodyACIS_ptr = CAST_TO(OSME_ptr, BodyACIS) ;

     // Make sure that we have a valid BodyACIS
   assert ( bodyACIS_ptr != NULL ) ;
     */
     // Return the BodyACIS
   return this_bodysm;
}

CubitStatus
AcisModifyEngine::webcut_with_sheet( DLIList<BodySM*>& webcut_body_list,
                                     BodySM *sheet_body,
                                     DLIList<BodySM*> &new_bodies,
                                     bool imprint )
{
   int num_cut=0;
   DLIList<BodySM*> temp_new_bodies;

   for ( int ii = webcut_body_list.size(); ii > 0; ii-- )
   {
     BodySM *webcut_body = webcut_body_list.get_and_step();
     BODY *sheet_BODY = AcisQueryEngine::get_BODY(sheet_body );
     BODY *webcut_BODY = AcisQueryEngine::get_BODY(webcut_body );
      //now webcut it.
     BODY *new_BODY_1, *new_BODY_2;

     CubitStatus status = webcut_with_sheet( webcut_BODY,
                                             sheet_BODY,
                                             new_BODY_1,
                                             new_BODY_2,
                                             imprint );
     if ( status == CUBIT_SUCCESS )
     {
       DLIList<BODY*> new_BODY_list;

       new_BODY_list.append(new_BODY_1);
       new_BODY_list.append(new_BODY_2);
       temp_new_bodies.clean_out();

       CubitStatus result = get_new_Body(webcut_body, webcut_BODY,
                                         new_BODY_list, temp_new_bodies,
                                         CUBIT_FALSE);
       if ( result != CUBIT_SUCCESS )
       {
         PRINT_ERROR("Problems with building volume.\n");
         num_cut = 0;
         return CUBIT_FAILURE;
       }
       else
       {
         new_bodies += temp_new_bodies;
         num_cut++;
       }
     }
   }

  // Return CUBIT_SUCCESS as long as one body is webcut, earlier failures should
  // return CUBIT_FAILURE independently of these return statments
  return num_cut > 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}


// This function does the actuall webcuting with two acis bodies...
CubitStatus AcisModifyEngine::webcut_with_sheet( BODY *webcut_body,
                                                 BODY *sheet_body,
                                                 BODY *&webcut_body_1,
                                                 BODY *&webcut_body_2,
                                                 bool imprint )
{
  //don't try to webcut a sheet body with itself
  if( webcut_body == sheet_body )
  {
    return CUBIT_FAILURE;
  }

  CubitBoolean new_ids = GeometryModifyTool::instance()->get_new_ids();
  CubitBoolean not_new_ids = CUBIT_TRUE;
  if (new_ids) not_new_ids = CUBIT_FALSE;

  BODY *sheet_BODY = sheet_body;
  BODY *copy_sheet_BODY = this->copy_BODY(sheet_BODY);
  BODY *webcut_BODY = webcut_body;
  BODY *copy_webcut_BODY1 = copy_BODY(webcut_BODY,
                                      new_ids);
  BODY *copy_webcut_BODY2 = copy_BODY(webcut_BODY,
                                      new_ids);
//first test to see if these intersect.  If they don't
//then continueing on will not show the error and duplicate
//bodies will be created.
  if ( !BODYs_interfering( webcut_BODY, sheet_BODY ) )
  {
    PRINT_WARNING("Cutting Sheet does not intersect the original volume.\n"
                  "         The original volume is restored.\n" );
      //These don't intersect.  Clean up and exit.
    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_sheet_BODY, CUBIT_TRUE);
    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_webcut_BODY1, not_new_ids);
    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_webcut_BODY2, not_new_ids);

    return CUBIT_FAILURE;
  }

  //set all faces to one sided
  API_BEGIN;
  sg_body_to_1d( copy_sheet_BODY );
  API_END;

    //Now we have to do two subtracts. Lets subtract this face from the body
    //to be webcut.  Then reverse it and do it again.
  BODY *copy_2_sheet_BODY = copy_BODY( copy_sheet_BODY );

  outcome rc_3 = api_reverse_body( copy_2_sheet_BODY );
  if ( !rc_3.ok() || copy_2_sheet_BODY == NULL )
  {
    PRINT_ERROR("Error in webcutting volume with sheet.\n");
    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_sheet_BODY, CUBIT_TRUE);
    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_2_sheet_BODY, CUBIT_TRUE);
    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_webcut_BODY1, not_new_ids);
    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_webcut_BODY2, not_new_ids);
    AcisQueryEngine::instance()->ACIS_API_error (rc_3);
    return CUBIT_FAILURE;
  }

  outcome rc_4 = api_subtract( copy_sheet_BODY, copy_webcut_BODY1 );

  if ( !rc_4.ok() || copy_sheet_BODY == NULL ||
       copy_webcut_BODY1 == NULL ||
       copy_webcut_BODY1->lump() == NULL ||
       (!BodyACIS::is_sheet_body( copy_webcut_BODY1 ) &&
       !is_closed_solid_body( copy_webcut_BODY1 ) ))
  {
    CubitBoolean is_error = CUBIT_TRUE;
    if ( rc_4.ok() && copy_webcut_BODY1 && copy_sheet_BODY &&
        ( copy_webcut_BODY1->lump() == NULL ||
          !is_closed_solid_body( copy_webcut_BODY1 )))
    {
        //This is where the sheet just intersects with the shell
        //of the outside body.  should do nothing here.
      is_error = CUBIT_FALSE;
    }
    if( !rc_4.ok() )  //was already deleted if api_subtract was successfull
      AcisQueryEngine::instance()->delete_ACIS_BODY(copy_sheet_BODY, CUBIT_TRUE);

    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_2_sheet_BODY, CUBIT_TRUE);
    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_webcut_BODY1, not_new_ids);
    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_webcut_BODY2, not_new_ids);
    if ( is_error )
    {
      PRINT_ERROR("Error in webcutting volume with sheet.\n");
      AcisQueryEngine::instance()->ACIS_API_error (rc_4);
    }
    else
      PRINT_WARNING("Sheet grazes shell of volume.  Nothing was cut.\n");
    return CUBIT_FAILURE;
  }
  outcome rc_5 = api_subtract( copy_2_sheet_BODY, copy_webcut_BODY2 );

  if ( !rc_5.ok() || copy_2_sheet_BODY == NULL ||
       copy_webcut_BODY2 == NULL ||
       copy_webcut_BODY2->lump() == NULL ||
       (!BodyACIS::is_sheet_body( copy_webcut_BODY2 ) &&
       !is_closed_solid_body( copy_webcut_BODY2 ) ))
  {
    CubitBoolean is_error = CUBIT_TRUE;
    if ( rc_5.ok() && copy_2_sheet_BODY && copy_webcut_BODY2
         && (copy_webcut_BODY2->lump() == NULL ||
         !is_closed_solid_body( copy_webcut_BODY2 )))
      is_error = CUBIT_FALSE;

    if( !rc_5.ok() )  //was already deleted if api_subtract was successfull
      AcisQueryEngine::instance()->delete_ACIS_BODY(copy_2_sheet_BODY, CUBIT_TRUE);

    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_webcut_BODY1, not_new_ids);
    AcisQueryEngine::instance()->delete_ACIS_BODY(copy_webcut_BODY2, not_new_ids);
    if ( is_error )
    {
      PRINT_ERROR("Error in webcutting volume with sheet.\n");
      AcisQueryEngine::instance()->ACIS_API_error (rc_5);
    }
    else
      PRINT_WARNING("Sheet grazes shell of volume.  Nothing was cut.\n");
    return CUBIT_FAILURE;
  }

  //it is possible that the sheet body tool and sheet blank overlap 
  //use different webcut function then
  
  //if copy_webcut_BODY1 is a sheet body, see if resultant bodies overlap 
  if( BodyACIS::is_sheet_body( copy_webcut_BODY1 ) )
  {
    bool overlap = false;
    
    API_BEGIN;
    API_NOP_BEGIN;
      outcome my_result = api_intersect( copy_webcut_BODY1, copy_webcut_BODY2 ); 
      if( my_result.ok() && copy_webcut_BODY2->lump() )
        overlap = true;
    API_NOP_END;
    API_END;

    //if resultant bodies overlap...use appropriate webcut function 
    if( overlap ) 
    {
      AcisQueryEngine::instance()->delete_ACIS_BODY(copy_webcut_BODY1, not_new_ids);
      AcisQueryEngine::instance()->delete_ACIS_BODY(copy_webcut_BODY2, not_new_ids);

      CubitStatus stat = AcisModifyEngine::webcut_BODY( webcut_body,
                                                        sheet_body,
                                                        copy_webcut_BODY1,
                                                        copy_webcut_BODY2); 
      if( stat == CUBIT_FAILURE )
        return CUBIT_FAILURE;
    }
  }

    //okay everything should be okay now.
  if ( imprint == CUBIT_TRUE )
  {
    CubitStatus stat = imprint_BODYs( copy_webcut_BODY1,
                                      copy_webcut_BODY2 );
    if ( stat != CUBIT_SUCCESS ||
         copy_webcut_BODY1 == NULL ||
         copy_webcut_BODY2 == NULL )
    {
      PRINT_ERROR("Error in webcutting volume with sheet.\n");
      return CUBIT_FAILURE;
    }
  }
  webcut_body_1 =  copy_webcut_BODY1;
  webcut_body_2 =  copy_webcut_BODY2;
  AcisQueryEngine::instance()->clear_bounding_box( webcut_body_1 );
  AcisQueryEngine::instance()->bounding_box( webcut_body_1 );
  AcisQueryEngine::instance()->clear_bounding_box( webcut_body_2 );
  AcisQueryEngine::instance()->bounding_box( webcut_body_2 );

  return CUBIT_SUCCESS;
}

CubitStatus
AcisModifyEngine::webcut_with_extended_surf( DLIList<BodySM*> &webcut_body_list,
                                             Surface *extend_from,
                                             DLIList<BodySM*> &new_bodies,
                                             int &num_cut,
                                             bool imprint )
{
    //now lets create a new extended face from the surface.
  SurfaceACIS *surf_ACIS = CAST_TO(extend_from, SurfaceACIS );
  if(!surf_ACIS)
  {
     PRINT_ERROR("Couldn't extend face to use for webcutting.\n"
                 "Possible incompatible geometry engine for this operation.\n");
     return CUBIT_FAILURE;
  }
  FACE *extend_from_FACE = surf_ACIS->get_FACE_ptr();

  FACE *extended_FACE = make_FACE( extend_from_FACE, CUBIT_TRUE );
  if ( extended_FACE == NULL )
  {
    PRINT_ERROR("Couldn't extend face to use for webcutting.\n");
    return CUBIT_FAILURE;
  }
  BodySM *webcut_body;
  DLIList<BodySM*> temp_new_bodies;
  num_cut = 0;
  for ( int ii = webcut_body_list.size(); ii > 0; ii-- )
  {
    webcut_body = webcut_body_list.remove();
    BODY *tool_BODY = AcisQueryEngine::instance()->get_BODY_of_ENTITY(extended_FACE);
    BODY *webcut_BODY = AcisQueryEngine::get_BODY(webcut_body );
      //now webcut with it.
    BODY *new_BODY_1, *new_BODY_2;
    CubitStatus status = webcut_with_sheet( webcut_BODY,
                                            tool_BODY,
                                            new_BODY_1,
                                            new_BODY_2,
                                            imprint );
    if ( status == CUBIT_SUCCESS )
    {
      DLIList<BODY*> new_BODY_list;
      new_BODY_list.append(new_BODY_1);
      new_BODY_list.append(new_BODY_2);
      temp_new_bodies.clean_out();

      CubitStatus result = get_new_Body(webcut_body, webcut_BODY,
                                        new_BODY_list, temp_new_bodies,
                                        CUBIT_FALSE);
      if ( result != CUBIT_SUCCESS )
      {
        PRINT_ERROR("Problems with building volume.\n");
        num_cut = 0;
        return CUBIT_FAILURE;
      }
      else
      {
        new_bodies += temp_new_bodies;
        num_cut++;
      }
    }
  }
  return CUBIT_SUCCESS;
}

CubitStatus AcisModifyEngine::regularize_body( BodySM *body_ptr,
                                               BodySM *&new_body_ptr)
{
  CubitBoolean new_ids = GeometryModifyTool::instance()->get_new_ids();
  CubitBoolean not_new_ids = (new_ids ? CUBIT_FALSE : CUBIT_TRUE);

   BODY *old_BODY = AcisQueryEngine::get_BODY(body_ptr );
   if (!old_BODY)
      return CUBIT_FAILURE;
   BODY *acis_BODY = this->copy_BODY(old_BODY, new_ids);

   outcome result = api_clean_entity( (ENTITY*) acis_BODY);
   //outcome result = api_regularise_entity( (ENTITY*) acis_BODY);
   if (!result.ok() || acis_BODY == NULL ||
       acis_BODY->lump() == NULL )
   {
      PRINT_ERROR("found in AcisModifyEngine::regularize_body.\n"
                  "api_regularise_entity failed.\n");
      AcisQueryEngine::instance()->ACIS_API_error (result);
      AcisQueryEngine::instance()->delete_ACIS_BODY(acis_BODY, not_new_ids);
      return CUBIT_FAILURE;
   }

   new_body_ptr = get_new_Body(body_ptr, old_BODY, acis_BODY, false);

   if (new_body_ptr) return CUBIT_SUCCESS;
   else return CUBIT_FAILURE;

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

CubitStatus AcisModifyEngine::regularize_entity( GeometryEntity *old_refentity_ptr,
                                                    BodySM *&new_body_ptr)
{
  CubitBoolean new_ids = GeometryModifyTool::instance()->get_new_ids();
  CubitBoolean not_new_ids = (new_ids ? CUBIT_FALSE : CUBIT_TRUE);

  BODY *copied_BODY_ptr;
  SurfaceACIS *ref_face = CAST_TO(old_refentity_ptr, SurfaceACIS);
  CurveACIS *ref_edge= CAST_TO(old_refentity_ptr, CurveACIS);
  PointACIS *ref_vertex= CAST_TO(old_refentity_ptr, PointACIS);

  DLIList<SurfaceACIS*> reg_face_list;
  DLIList<FACE*> FACE_list;

  DLIList<CurveACIS*> reg_edge_list;
  DLIList<EDGE*> EDGE_list;

  DLIList<PointACIS*> reg_vertex_list;
  DLIList<VERTEX*> VERTEX_list;

  // get the copied_BODY_ptr and FACE_list to clean
  if (ref_face)
  {
	  reg_face_list.append(ref_face);
      if (get_copied_FACES_of_body( reg_face_list, FACE_list, copied_BODY_ptr )== CUBIT_FAILURE )
	  {
         PRINT_ERROR("Surface has no volume associated with it.\n");
          return CUBIT_FAILURE;
	  }
  }


// get the copied_BODY_ptr and EDGE_list to clean
  if (ref_edge)
  {
	  reg_edge_list.append(ref_edge);
      if(get_copied_EDGES_of_body( reg_edge_list, EDGE_list, copied_BODY_ptr )== CUBIT_FAILURE )
	  {
         PRINT_ERROR("Edge has no volume associated with it.\n");
         return CUBIT_FAILURE;
	  }
  }

 if (ref_vertex)
  {
	  reg_vertex_list.append(ref_vertex);
      if (get_copied_VERTICES_of_body( reg_vertex_list, VERTEX_list, copied_BODY_ptr )== CUBIT_FAILURE )
	  {
         PRINT_ERROR("Vertex has no volume associated with it.\n");
         return CUBIT_FAILURE;
	  }
  }

  if (!ref_face && !ref_edge && !ref_vertex)
  {
	  //complain
      PRINT_ERROR("RefEntity is neither a RefEdge, RefFace, nor RefVertex.\n");
	  PRINT_ERROR("found in AcisModifyEngine::regularize_refentity.\n"
					  "api_regularise_refentity aborted.\n");
	  AcisQueryEngine::instance()->delete_ACIS_BODY(copied_BODY_ptr, not_new_ids);
	   return CUBIT_FAILURE;
  }

  // Get original Body and BODY
  BodySM *body_ptr = AcisQueryEngine::instance()->get_body_sm_of_ENTITY( copied_BODY_ptr );

  BODY *BODY_ptr =  AcisQueryEngine::get_BODY(body_ptr);

  int i;
  FACE *acis_FACE;
  EDGE *acis_EDGE;
  VERTEX *acis_VERTEX;

 //clean FACE
  for (i=0;i<FACE_list.size();i++)
  {
	  acis_FACE = FACE_list.get_and_step();
	  outcome result = api_clean_entity( (ENTITY*) acis_FACE);
	  if (!result.ok() || acis_FACE == NULL )
	  {
		  PRINT_ERROR("found in AcisModifyEngine::regularize_refentity.\n"
					  "api_regularise_refentity failed.\n");
		 AcisQueryEngine::instance()->ACIS_API_error (result);
		 AcisQueryEngine::instance()->delete_ACIS_BODY(copied_BODY_ptr, not_new_ids);
		 return CUBIT_FAILURE;
	  }
  }

  // clean EDGE
  for (i=0;i<EDGE_list.size();i++)
  {
	  acis_EDGE = EDGE_list.get_and_step();
	  outcome result = api_clean_entity( (ENTITY*) acis_EDGE);
	  if (!result.ok() || acis_EDGE == NULL )
	  {
		  PRINT_ERROR("found in AcisModifyEngine::regularize_refentity.\n"
					  "api_regularise_refentity failed.\n");
		  AcisQueryEngine::instance()->ACIS_API_error (result);
		  AcisQueryEngine::instance()->delete_ACIS_BODY(copied_BODY_ptr, not_new_ids);
		  return CUBIT_FAILURE;
	  }
  }

  // clean VERTEX
  for (i=0;i<VERTEX_list.size();i++)
  {
	  acis_VERTEX = VERTEX_list.get_and_step();
	  outcome result = api_clean_entity( (ENTITY*) acis_VERTEX);
	  if (!result.ok() || acis_VERTEX == NULL )
	  {
		  PRINT_ERROR("found in AcisModifyEngine::regularize_refentity.\n"
					  "api_regularise_refentity failed.\n");
		  AcisQueryEngine::instance()->ACIS_API_error (result);
		  AcisQueryEngine::instance()->delete_ACIS_BODY(copied_BODY_ptr, not_new_ids);
		  return CUBIT_FAILURE;
	  }
  }
  //rebuild the cubit geometry
  new_body_ptr = get_new_Body( body_ptr,  BODY_ptr , copied_BODY_ptr , false);
  if (new_body_ptr) return CUBIT_SUCCESS;
	   else return CUBIT_FAILURE;
}

CubitStatus AcisModifyEngine::split_body( BodySM *body_ptr,
                                          DLIList<BodySM*> &new_bodies )
{
  bool delete_attribs = GeometryModifyTool::instance()->get_new_ids() != 0;
     //get the acis body.
   BODY* acis_BODY = NULL;
   BODY **new_BODY_list;
   int n_body;

   BODY *old_BODY = AcisQueryEngine::get_BODY(body_ptr );
   acis_BODY = this->copy_BODY(old_BODY, delete_attribs);

   outcome result = api_separate_body( acis_BODY, n_body, new_BODY_list );

   if (!result.ok())
   {
      PRINT_ERROR("Error in seperating volumes.\n");
      AcisQueryEngine::instance()->delete_ACIS_BODY(acis_BODY, CUBIT_TRUE);
      return CUBIT_FAILURE;
   }

   DLIList<BODY*> new_BODYs;
   for (int i = 0; i < n_body; i++)
     new_BODYs.append(new_BODY_list[i]);

   CubitStatus status = get_new_Body(body_ptr, old_BODY,
                                     new_BODYs, new_bodies,
                                     CUBIT_FALSE);
   return status;
}

CubitStatus AcisModifyEngine::split_periodic( BodySM *body_ptr,
                                              BodySM *&new_body)
{

   BODY* acis_BODY = NULL;

     //get the ACIS BODY of tool_body
   BODY *old_BODY = AcisQueryEngine::get_BODY( body_ptr );
   acis_BODY = this->copy_BODY(old_BODY,
                               GeometryModifyTool::instance()->get_new_ids());

   outcome some_result = api_split_periodic_faces( acis_BODY );

   if ( acis_BODY == NULL || !some_result.ok() )
   {
      if (acis_BODY != NULL )
          AcisQueryEngine::instance()->delete_ACIS_BODY( acis_BODY );

      PRINT_ERROR("Unable to split periodic surface.\n");
      return CUBIT_FAILURE;
   }

   new_body = get_new_Body(body_ptr, old_BODY, acis_BODY, false);

   if ( new_body == NULL )
   {
      PRINT_ERROR("Can not make a volume out of the split volumes.\n");
      return CUBIT_FAILURE;
   }

   else return CUBIT_SUCCESS;
}

EDGE *
AcisModifyEngine::make_spline_EDGE( DLIList<CubitVector*> &vec_list )
{
  int i;
  EDGE *new_EDGE_ptr = NULL;

  // Create a spline
  int num_pnts = vec_list.size();
  SPAposition* pos_array = new SPAposition[num_pnts];

  // Fill the vector_list of points
  vec_list.reset();
  for( i=0; i<num_pnts; i++)
  {
    // Set the coordinates of the SPAposition object
    pos_array[i].set_x(vec_list.get()->x());
    pos_array[i].set_y(vec_list.get()->y());
    pos_array[i].set_z(vec_list.get()->z());
    vec_list.step();
  }

  outcome result = api_curve_spline(num_pnts, pos_array, NULL, NULL, new_EDGE_ptr);
  delete [] pos_array;
  if( !result.ok() || new_EDGE_ptr == NULL )
  {
    get_acis_query_engine()->ACIS_API_error(result);
    PRINT_ERROR( "unable to make spline from list of positions\n" );
    return NULL;
  }

  return new_EDGE_ptr;
}

//-------------------------------------------------------------------------
// Purpose       : Reverse the ACIS BODY
//
// Special Notes : Moved from BodyACIS - j.k, 2004-5-25
//
// Creator       : Tim Tautges
//
// Creation Date : 3/18/97
//-------------------------------------------------------------------------
CubitStatus AcisModifyEngine::reverse_body( BodySM* body_to_reverse )
{
  BODY* body_ptr = AcisQueryEngine::get_BODY( body_to_reverse );
  if (!body_ptr)
    return CUBIT_FAILURE;

  outcome result = api_reverse_body( body_ptr );
  if (!result.ok())
  {
    AcisQueryEngine::instance()->ACIS_API_error( result, "reversing BODY." );
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

void AcisModifyEngine::webcut_cleanup( DLIList<BodySM*>& body_list,
                                         ENTITY_LIST& BODY_list) const
{
     // This procedure is called by webcut_handler and gets rid of
     // all the Body's and Acis BODY's created during the execution
     // of the webcut command. The command was unsuccessful and hence,
     // these entities need to be removed so the system can be returned
     // to its state before the command was invoked.

     // Delete the Bodies
   AcisQueryEngine::instance()->delete_solid_model_entities(body_list);
   //GeometryQueryTool::instance()->delete_Body(body_list, CUBIT_TRUE);

     // Delete the Acis BODYs
   AcisQueryEngine::instance()->delete_ACIS_BODY(BODY_list);

   return;

}

void AcisModifyEngine::webcut_imprint(
    BODY* cutting_tool_BODY_ptr,
    DLIList<BodySM*> &old_body_list,
    DLIList<BODY*>& new_webcut_BODY_list,
    DLIList<BODY*>& just_webcut_BODY_list,
    DLIList<BodySM*>& results_list,
    DLIList<BodySM*> &imprinted_Body_from_Model_list,
    DLIList<BODY*> &imprinted_BODY_from_Model_list) const
{
     // This procedure takes care of the imprint operations that are required
     // by the webcut operation. Each of the BODYs in the new_webcut_BODY_list
     // is tested against every other BODY in the model as well as every other
     // BODY in the same list for interference. If interference exists, then
     // the BODYs are imprinted with each other and new (corresponding) Bodies
     // are created. Pointers to these new Bodies are inserted into the list,
     // new_webcut_body_list.

     // MJP Note:
     // Note that this procedure takes care of imprinting the existing Bodies
     // (i.e., the ones in the model) and creating new Bodies when necessary
     // to replace these. However, new Body's are *not* created from the list
     // of ACIS BODYs in the input list, new_webcut_BODY_list. These are simply
     // returned to the calling routine, where it is expected that new Bodies
     // will be created, as required.

   int webcut_debug_flag = 18;
   int webcut_debug_flag_on = DEBUG_FLAG(webcut_debug_flag);

   int i = 0;
   int j = 0;

     // Make copies of the input BODYs
   DLIList<BODY*> copy_new_webcut_BODY_list ;
   DLIList<BODY*> new_results_BODYS;
   new_webcut_BODY_list.reset() ;
   for (i = 0; i < new_webcut_BODY_list.size() ; i++)
   {
      BODY* old_BODY = new_webcut_BODY_list.get_and_step() ;
      BODY* new_BODY = this->copy_BODY(old_BODY,
                                       GeometryModifyTool::instance()->get_new_ids()) ;

        // Make sure the the new BODY is a valid BODY
      assert( new_BODY != NULL) ;

        // Add the new BODY to the backup list
      copy_new_webcut_BODY_list.append(new_BODY) ;
   }

     // Loop over each of the new BODYs and imprint it with any other new BODY
     // that touches it.
   BODY* first_BODY_ptr = NULL;
   BODY* second_BODY_ptr = NULL;
   copy_new_webcut_BODY_list.reset() ;
   for (i = 0; i < (copy_new_webcut_BODY_list.size() - 1); i++)
   {
        // Get the first BODY
      first_BODY_ptr = copy_new_webcut_BODY_list.get_and_step();

      for (j = (i + 1); j < copy_new_webcut_BODY_list.size(); j++)
      {
           // Get the second BODY
         second_BODY_ptr = copy_new_webcut_BODY_list.get_and_step();

           // If these BODYs touch, then imprint them
         CubitBoolean result_interfer =
             BODYs_interfering( first_BODY_ptr,second_BODY_ptr );

         if ( result_interfer == CUBIT_TRUE )
         {
            CubitStatus resultImprint =
                imprint_BODYs( first_BODY_ptr, second_BODY_ptr );

            if( resultImprint == CUBIT_FAILURE )
            {
               PRINT_WARNING ("Imprint operation failed.\n"
                              "         Going back to the previous state");

                 // Delete the new BODYs and return.
               copy_new_webcut_BODY_list.reset() ;
               for (int k = 0 ; k < copy_new_webcut_BODY_list.size() ; k++)
               {
                  first_BODY_ptr = copy_new_webcut_BODY_list.get_and_step() ;
                  AcisQueryEngine::instance()->delete_ACIS_BODY(first_BODY_ptr) ;
               }
               return ;
            }
         }
      }

        // Go to the (i+1)th element of the list
      copy_new_webcut_BODY_list.reset() ;
      copy_new_webcut_BODY_list.step(i+1) ;

   }

/*
  CDE-- commented out this code.  For now, we've decided that when webcutting 
  with the 'imprint' option, the imprint should occur between bodies involved
  in the webcut.  This code was imprinting all resulting pieces of the webcut 
  with all the other Bodies in the model.

     // Get the Bodys in the Model
   DLIList<Body*> Body_in_Model_list_tmp;
   GeometryQueryTool::instance()->bodies(Body_in_Model_list_tmp);
   DLIList<BodySM*> Body_in_Model_list;
   Body_in_Model_list_tmp.reset();
   for (i = Body_in_Model_list_tmp.size(); i--; )
    Body_in_Model_list.append( Body_in_Model_list_tmp.get_and_step()->get_body_sm_ptr());

     // Loop over each of the new BODYs and imprint it with any of the
     // BODYs already existing in the model, if they touch.

     // List of new imprinted BODYs -- Bodies need to be created from these
   DLIList<BODY*> copy_imprinted_BODY_from_Model_list;

   BodySM* model_Body_ptr = NULL;
   BODY* model_BODY_ptr = NULL;
   BODY* new_model_BODY_ptr = NULL;
   BODY* new_BODY_ptr = NULL ;
   int copy_made = CUBIT_FALSE;

   Body_in_Model_list.reset();
     // Loop over all the existing Bodies
   for (i = 0; i < Body_in_Model_list.size(); i++)
   {
        // Get this existing Body's parent ACIS BODY
      model_Body_ptr = Body_in_Model_list.get_and_step() ;
      model_BODY_ptr = AcisQueryEngine::get_BODY(model_Body_ptr) ;

        // Make sure that model_BODY_ptr is not the Cutting Tool as the
        // user is allowed to use an existing Body as a Cutting Tool.
        // If this is the case, then don't involve it in the imprint
        // operations.
      if (model_BODY_ptr == cutting_tool_BODY_ptr || old_body_list.move_to(model_Body_ptr))
      {
           // Don't do anything
         continue;
      }

        // Loop over all the newly created (by webcutting) ACIS BODYs and test
        // the existing Body against them (for imprinting purposes, if there
        // is interference)
      copy_made = CUBIT_FALSE;
      new_model_BODY_ptr = NULL;
      copy_new_webcut_BODY_list.reset();

      for (j = 0; j < copy_new_webcut_BODY_list.size(); j++)
      {
           // Get the BODY that we're going to test the existing Body against
         new_BODY_ptr = copy_new_webcut_BODY_list.get_and_step() ;

           // DEBUG stuff
         if (webcut_debug_flag_on)
         {
              //Save the two BODYs
            AcisQueryEngine::instance()->save_ENTITY_as_sat_file(new_BODY_ptr, "BODY1.sat", "w");
            AcisQueryEngine::instance()->save_ENTITY_as_sat_file(model_BODY_ptr, "BODY2.sat", "w");
         }

           // If these BODYs touch, then ...
         if ( BODYs_interfering(new_BODY_ptr, model_BODY_ptr) == CUBIT_TRUE )
         {
              // DEBUG stuff
            if (webcut_debug_flag_on)
            {
               AcisQueryEngine::instance()->save_ENTITY_as_sat_file(new_BODY_ptr,
                                             "BODY1_after_test.sat",
                                             "w");
               AcisQueryEngine::instance()->save_ENTITY_as_sat_file(new_model_BODY_ptr,
                                             "BODY2_after_test.sat",
                                             "w");
            }

              // Because of the fact that the RefEntity datastructure is
              // nonmanifold and hence shares underlying ACIS ENTITYs, I have to:
              //   a) make a copy of the parent ACIS BODY of the Body being
              //      deactivated (this is done only once for the inner loop)
              //   b) then all imprinting is done with this copy instead of
              //     with the original parent ACIS BODY.
              // deactivateRefStructure takes care of detaching the RefEntities
              // and their associated ACIS ENTITYs (removing the double links
              // between them).
            if (!copy_made) // Make this copy at most once per outer loop iteration
            {
                 // Make a copy of the parent ACIS BODY.
               new_model_BODY_ptr = this->copy_BODY(model_BODY_ptr,
                                                    GeometryModifyTool::instance()->get_new_ids());

                 // Make sure a valid copy was made. If not, then we need to exit
                 // this loop (i.e., stop testing this Body for imprinting
                 // and go on to the next one).
               if (new_model_BODY_ptr == NULL)
               {
                  PRINT_ERROR("Cannot test imprinting for volume.\n" );
                  break;
               }
               else
               {
                  copy_made = CUBIT_TRUE;
                  if ( results_list.move_to(model_Body_ptr) )
                    new_results_BODYS.append( new_model_BODY_ptr );
                  imprinted_Body_from_Model_list.append(model_Body_ptr) ;
                  imprinted_BODY_from_Model_list.append(model_BODY_ptr) ;
                  copy_imprinted_BODY_from_Model_list.append(new_model_BODY_ptr) ;
               }
            }

              // Imprint the BODYs (use the new copied BODY). If the imprint
              // operation fails, the state of these BODYs is left unchanged.
            if (imprint_BODYs(new_BODY_ptr, new_model_BODY_ptr) ==  CUBIT_SUCCESS )
            {
                 // DEBUG stuff, if the imprint operation succeeded.
               if (webcut_debug_flag_on)
               {
                  AcisQueryEngine::instance()->save_ENTITY_as_sat_file(new_BODY_ptr,
                                                "BODY1_after_imprint.sat",
                                                "w");
                  AcisQueryEngine::instance()->save_ENTITY_as_sat_file(new_model_BODY_ptr,
                                                "BODY2_after_imprint.sat",
                                                "w");
               }
            }

            else
            {
               PRINT_WARNING ("Imprint operation failed.\n"
                              "         Going back to the previous state");

                 // Delete the new BODYs and return.
               copy_new_webcut_BODY_list.reset() ;
               int k;
               for (k = 0 ; k < copy_new_webcut_BODY_list.size() ; k++)
               {
                  new_BODY_ptr = copy_new_webcut_BODY_list.get_and_step() ;
                  AcisQueryEngine::instance()->delete_ACIS_BODY(first_BODY_ptr) ;
               }

               copy_imprinted_BODY_from_Model_list.reset() ;
               for (k = 0 ; k < copy_imprinted_BODY_from_Model_list.size() ; k++)
               {
                  new_model_BODY_ptr = copy_imprinted_BODY_from_Model_list.get_and_step() ;
                  AcisQueryEngine::instance()->delete_ACIS_BODY(new_model_BODY_ptr) ;
               }

               return ;
            }
         }
      } // Loop over the new (webcut) BODYs

   } // Loop over the existing Bodies in Model
*/
     // Now that we have come so far successfully, create the final list
     // of new BODYs that will be used to create VGI Bodys.
   new_webcut_BODY_list.clean_out() ;

     // First the BODYs that were copies of the webcut BODYs
   BODY* new_BODY_ptr = NULL ;
   copy_new_webcut_BODY_list.reset() ;
   for (i = 0 ; i < copy_new_webcut_BODY_list.size() ; i++)
   {
      new_BODY_ptr = copy_new_webcut_BODY_list.get_and_step() ;
      new_webcut_BODY_list.append(new_BODY_ptr) ;
   }
   just_webcut_BODY_list = new_webcut_BODY_list;

/*
     // Then the BODYs that were copies of Model BODYs
   copy_imprinted_BODY_from_Model_list.reset() ;
   for (i = 0 ; i < copy_imprinted_BODY_from_Model_list.size() ; i++)
   {
      new_BODY_ptr = copy_imprinted_BODY_from_Model_list.get_and_step() ;
      new_webcut_BODY_list.append(new_BODY_ptr) ;
   }
*/
     // now, two important lists:
     // a) just_webcut_BODY_list: bodies in original bodies-to-be-webcut list
     //        that were modified, either by webcut or by imprint
     // b) new_webcut_BODY_list: all bodies modified by webcut or imprint,
     //        including model bodies affected by imprint

   return;
}

BODY* AcisModifyEngine::create_infinite_plane_cutting_tool(
    const CubitVector &input_p1,
    const CubitVector &input_p2,
    const CubitVector &input_p3,
    const SPAbox& box_to_cut,
    bool just_face ) const
{
    // Create a plane from the vertices
    // Note:  The old infinite plane code always swept in the opposite direction
    //        from the normal of the passed in points, so we create the plane with
    //        the points in opposite order.
  CubitPlane p;
  p.mk_plane_with_points(input_p3, input_p2, input_p1);

    // Convert it to an ACIS plane
  SPAunit_vector plane_normal(p.normal().x(),
                           p.normal().y(),
                           p.normal().z());

  plane the_plane(SPAposition(input_p1.x(),
                           input_p1.y(),
                           input_p1.z()),
                  plane_normal);

  if (DEBUG_FLAG(18))
    PRINT_DEBUG_18(
                "plane created with root point (%g %g %g) and\n"
                "  normal (%g %g %g).\n",
                input_p1.x(),
                input_p1.y(),
                input_p1.z(),
                plane_normal.x(),
                plane_normal.y(),
                plane_normal.z());

    // Need to declare this up here so that it is outside of the API block
  BODY* infinite_plane_BODY = NULL;
  API_BEGIN;

    // Make a PLANE from it
  PLANE* the_PLANE = new PLANE(the_plane);

    // Make a FACE out of it
  FACE* planar_FACE = new FACE(NULL, NULL, the_PLANE, FORWARD);

    // Get the midpoint and diagonal length of the
    // SPAbox this needs to fit around.  Make diag length a little
    // bigger than necessary as a buffer.
  SPAposition midpoint = box_to_cut.mid();
  double diag_length = 1.01 * (box_to_cut.high() - midpoint).len();

    // Get the point on the plane that the SPAbox midpoint projects to.
  ENTITY_LIST temp_list;
  double distance;
  double *dist_ptr = &distance;
  ENTITY* entity_ptr = planar_FACE;
    // Test points behind plane
  result = api_ray_test_ents(midpoint, plane_normal, 1, 1, 1, &entity_ptr,
                             temp_list, dist_ptr);
  if (!result.ok())
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
  }
  if (temp_list.count() < 1)
  {
      // Test points in front of plane, if needed
    result = api_ray_test_ents(midpoint, -plane_normal, 1, 1, 1, &entity_ptr,
                               temp_list, dist_ptr);
    if (!result.ok())
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
    }
      // Indicate that the distance is in front of the plane
    distance = -distance;
  }

  if (DEBUG_FLAG(18))
  {
    PRINT_DEBUG_18(
                "Midpoint before moving it: %g, %g, %g.\n",
                midpoint.x(), midpoint.y(), midpoint.z());
  }

    // Move midpoint from SPAbox to plane
  midpoint += distance * plane_normal;

  if (DEBUG_FLAG(18))
  {
    PRINT_DEBUG_18(
                "just_face = %d.\n"
                "The center of the box was %s the plane.\n",
                just_face,
                distance >= 0 ? "behind" : "in front of");
    PRINT_DEBUG_18(
                "Distance to center of box was determined to be %g.\n",
                distance);
    PRINT_DEBUG_18( "Mapped point is %g, %g, %g.\n",
                midpoint.x(), midpoint.y(), midpoint.z());
  }

    // Get 4 points on the plane that form a square with sides of
    // length 2*diag_length.
  SPAvector in_plane_1 = diag_length * normalise(the_plane.u_deriv);
  SPAvector in_plane_2 = plane_normal * in_plane_1;
  SPAposition corners[4];
  corners[0] = midpoint - in_plane_1 - in_plane_2;
  corners[1] = midpoint + in_plane_1 - in_plane_2;
  corners[2] = midpoint + in_plane_1 + in_plane_2;
  corners[3] = midpoint - in_plane_1 + in_plane_2;
  if (DEBUG_FLAG(18))
  {
      // Uncomment this if you want to see the plane before it's swept
//     GfxDebug::draw_line(
//       corners[0].x(), corners[0].y(), corners[0].z(),
//       corners[1].x(), corners[1].y(), corners[1].z(), CUBIT_RED);
//     GfxDebug::draw_line(
//       corners[1].x(), corners[1].y(), corners[1].z(),
//       corners[2].x(), corners[2].y(), corners[2].z(), CUBIT_GREEN);
//     GfxDebug::draw_line(
//       corners[2].x(), corners[2].y(), corners[2].z(),
//       corners[3].x(), corners[3].y(), corners[3].z(), CUBIT_BLUE);
//     GfxDebug::draw_line(
//       corners[3].x(), corners[3].y(), corners[3].z(),
//       corners[0].x(), corners[0].y(), corners[0].z(), CUBIT_CYAN);
//     GfxDebug::draw_line(
//       midpoint.x(), midpoint.y(), midpoint.z(),
//       midpoint.x()+plane_normal.x(),
//       midpoint.y()+plane_normal.y(),
//       midpoint.z()+plane_normal.z(), CUBIT_BLUE);
//     GfxDebug::flush();

    PRINT_DEBUG_18( "Diag length is %g.\n"
                "Got four corners: (%g %g %g)\n"
                "                  (%g %g %g).\n"
                "                  (%g %g %g).\n"
                "                  (%g %g %g).\n",
                diag_length,
                corners[0].x(), corners[0].y(), corners[0].z(),
                corners[1].x(), corners[1].y(), corners[1].z(),
                corners[2].x(), corners[2].y(), corners[2].z(),
                corners[3].x(), corners[3].y(), corners[3].z());
  }
    // Make some geometry and topology from the corners
  STRAIGHT* S[4];
  S[0] = new STRAIGHT(corners[0], normalise(in_plane_1));
  S[1] = new STRAIGHT(corners[1], normalise(in_plane_2));
  S[2] = new STRAIGHT(corners[2], normalise(-in_plane_1));
  S[3] = new STRAIGHT(corners[3], normalise(-in_plane_2));
  APOINT* P[4];
  VERTEX* V[4];
  EDGE*   E[4];
  COEDGE* CE[4];
  int i;
  for ( i = 0; i < 4; i++)
  {
    P[i] = new APOINT (corners[i]);
    V[i] = new VERTEX (P[i]);
  }
  CE[0] = NULL;
  for (i = 0; i < 4; i++)
  {
    E[i] = new EDGE(V[i], V[(i+ 1)%4], S[i], FORWARD);
    CE[i] = new COEDGE(E[i], FORWARD,
                       i == 0 ? NULL : CE[i-1],
                       i == 3 ? CE[0] : NULL);
  }
  LOOP* LOOP1 = new LOOP (CE[0], NULL);
  FACE* FACE1 = new FACE (LOOP1, NULL, the_PLANE, FORWARD);
  SHELL* SHELL1 = new SHELL(FACE1, NULL, NULL);
  LUMP* LUMP1 = new LUMP(SHELL1, NULL);
  infinite_plane_BODY = new BODY(LUMP1);
    // At this point we are done with planar_FACE.  It has been replaced by
    // a FACE that is bounded by a LOOPm FACE1.
  planar_FACE->lose();
    // Let's make sure everything went as planned.
  if (infinite_plane_BODY != NULL)
  {
    result = api_body_to_2d(infinite_plane_BODY);
    if (!result.ok() || infinite_plane_BODY == NULL)
    {
      PRINT_ERROR("Unable to create infinite plane.\n");
      AcisQueryEngine::instance()->ACIS_API_error(result);
      infinite_plane_BODY = NULL;
    }
      // Just leave it as a 2D body if we are creating a sheet, or
      // if the plane never intersects the input box_to_cut.
    else if (!just_face && distance < diag_length && distance > -diag_length)
    {
        // Now sweep the FACE to make a 3D body
      if (DEBUG_FLAG(18))
        PRINT_DEBUG_18(
                    "Sweeping (%g %g %g).\n",
                    ((diag_length-distance)*plane_normal).x(),
                    ((diag_length-distance)*plane_normal).y(),
                    ((diag_length-distance)*plane_normal).z());

        // A negative distance means that it was in front of the plane, so
        // we need to sweep extra distance.  A positive distance means behind the
        // plane, so we don't need to sweep as far.
      sweep_options sweep_ops;
      logical ACIS_primary_side = true;
      sweep_ops.set_which_side(ACIS_primary_side);
      result = api_sweep_with_options(FACE1,
                                      (diag_length-distance)*plane_normal,
                                      &sweep_ops,
                                      infinite_plane_BODY);
      if (!result.ok())
      {
        PRINT_ERROR("In GeometryModifyTool::create_infinite_plane_cutting_tool.\n"
                    "       Problem sweeping the infinite plane.\n");
        AcisQueryEngine::instance()->ACIS_API_error(result);
        infinite_plane_BODY = NULL;
      }
    }
    else if (DEBUG_FLAG(18))
      PRINT_DEBUG_18( "Not sweeping the infinite plane.\n");
  }

  API_END;

  return infinite_plane_BODY;
}

//***********************************************************
// Function: create_solid_bodies_from_surfs
// Description: This function assumes that the reffaces sent into
// this function are either sheet bodies, or free surfaces.  This
// Will have been taken care of in the calling function.  GT?
// All the surfaces, in ACIS, will be turned into sheet bodies.  They
// will be united, then the void enclosed.  The result will be healed
// if the option is specified (default is to heal, from experience).
// This function will not create new reffaces.
//***********************************************************
CubitStatus AcisModifyEngine::create_solid_bodies_from_surfs( DLIList<Surface*> &ref_face_list,
                                                      DLIList<BodySM*> &new_bodies,
                                                      bool keep_old,
                                                      bool heal ) const
{
  int ii;
  Surface* ref_face;
  outcome rc;
    //First sort the ref_face_list to sheet_bodies and free faces.
  DLIList<BodySM*> body_list;
  DLIList<Surface*> free_faces;

  AcisQueryEngine* aqe = AcisQueryEngine::instance();

  ref_face_list.reset();
  for( ii = ref_face_list.size(); ii > 0; ii-- )
  {
    ref_face = ref_face_list.get_and_step();
    FACE* FACE_ptr = aqe->get_FACE(ref_face);
    BodySM* body_ptr = aqe->get_body_sm_of_ENTITY(FACE_ptr);

    if (body_ptr)
      body_list.append_unique(body_ptr);
    else
      free_faces.append_unique(ref_face);
  }

    //Now for all of the free_faces, create BODYS for them.
  DLIList<BODY*> BODY_list, old_BODY_list;
  for ( ii = free_faces.size(); ii > 0; ii-- )
  {
    ref_face = free_faces.get_and_step();
    FACE *acis_FACE = AcisQueryEngine::get_FACE(ref_face);
    if (!acis_FACE)
      return CUBIT_FAILURE;

    FACE *face_list[1];
    BODY *sheet_body;
    face_list[0] = acis_FACE;
    rc = api_sheet_from_ff(1, face_list, sheet_body);
    if (!rc.ok())
    {
      PRINT_ERROR("Couldn't build a volume with Surface.\n") ;
      for (int jj = BODY_list.size(); jj > 0; jj-- )
        AcisQueryEngine::instance()->delete_ACIS_BODY(BODY_list.remove());
      AcisQueryEngine::instance()->ACIS_API_error (rc);
      return CUBIT_FAILURE;
    }
    rc = api_body_to_2d(sheet_body);
    if (!rc.ok())
    {
      PRINT_ERROR("Couldn't build a volume with Surface.\n" );
      for (int jj = BODY_list.size(); jj > 0; jj-- )
        AcisQueryEngine::instance()->delete_ACIS_BODY(BODY_list.remove());
      AcisQueryEngine::instance()->ACIS_API_error (rc);
      return CUBIT_FAILURE;
    }
      // Make sure we were successful
    if ( sheet_body->lump() == NULL ||
         sheet_body->lump()->shell() == NULL ||
         sheet_body->lump()->shell()->first_face() == NULL )
    {
      PRINT_ERROR("Couldn't build a volume with Surface.\n" );
      for (int jj = BODY_list.size(); jj > 0; jj-- )
        AcisQueryEngine::instance()->delete_ACIS_BODY(BODY_list.remove());
      return CUBIT_FAILURE;
    }
    BODY_list.append(sheet_body);
  }
  BODY *first_BODY = NULL;

  for ( ii = body_list.size(); ii > 0; ii-- )
  {
    BodySM *body_ptr = body_list.get_and_step();
    BODY *BODY_ptr = AcisQueryEngine::get_BODY(body_ptr);
    if (!BODY_ptr)
      return CUBIT_FAILURE;

    if (first_BODY == NULL )
      first_BODY = BODY_ptr;
    BODY_list.append(BODY_ptr);
    old_BODY_list.append(BODY_ptr);
  }
    //Okay unite all of these bodies into one body.
  BODY_list.reset();
  int new_ids = GeometryModifyTool::instance()->get_new_ids();
  bool delete_attribs = (new_ids || keep_old);
  BODY *master;

  //if we're not healing, we're assuming the surfaces have gaps between each other
  //less than 1e-6, so we'll use stitching, otherwise we assume the gaps are bigger than
  //1e-6 and we'll us api_combine_body, and subsequently we'll heal
  ENTITY_LIST bodies_to_be_stitched;

  if ( first_BODY )
  {
    master = this->copy_BODY(first_BODY, delete_attribs);
    BODY_list.remove(first_BODY);
  }
  else
    master = this->copy_BODY(BODY_list.remove(), delete_attribs);

  if( !heal )
    bodies_to_be_stitched.add( master );

  for ( ii = BODY_list.size(); ii > 0; ii-- )
  {
    BODY *BODY_ptr = BODY_list.get_and_step();
    BODY *copy = this->copy_BODY(BODY_ptr, delete_attribs);
    if (copy == NULL)
    {
      AcisQueryEngine::instance()->delete_ACIS_BODY(master);
      return CUBIT_FAILURE;
    }
    if( !heal )
      bodies_to_be_stitched.add( copy );
    else
    {
      outcome result = api_combine_body(copy, master );

      if (!result.ok() || (!master ))
      {
        // If there are no LUMPs in the resulting BODY, return with a failure
        if (master == NULL)
          PRINT_ERROR("In AcisModifyEngine::create_solid_bodies_from_surfs\n");

        AcisQueryEngine::instance()->ACIS_API_error(result, "create body from surfs");
        if (master != NULL) AcisQueryEngine::instance()->delete_ACIS_BODY(master);
        if (copy   != NULL) AcisQueryEngine::instance()->delete_ACIS_BODY(copy);

        return CUBIT_FAILURE;
      }
    }
  }

  if( !heal )
  {
    ENTITY_LIST dummy_list, resultant_bodies;
    tolerant_stitch_options stitch_opts;
    //surfaces better be between GEOMETRY_RESABS or else fail
    stitch_opts.set_max_stitch_tol( GEOMETRY_RESABS );
    outcome result = api_stitch( bodies_to_be_stitched, dummy_list, resultant_bodies, &stitch_opts );

    master = (BODY*)dummy_list.first();

    if (!result.ok() || (!master ))
    {
      // If there are no LUMPs in the resulting BODY, return with a failure
      if (master == NULL)
        PRINT_ERROR("In AcisModifyEngine::create_solid_bodies_from_surfs\n");

      AcisQueryEngine::instance()->ACIS_API_error(result, "create body from surfs");
      if (master != NULL) AcisQueryEngine::instance()->delete_ACIS_BODY(master);

      return CUBIT_FAILURE;
    }
  }

    // If all went well, delete the old Bodies
  DLIList<TopologyBridge*> old_entities;
  body_list.reset() ;

  CAST_LIST_TO_PARENT(body_list, old_entities);

  for ( ii = free_faces.size(); ii > 0; ii-- )
    old_entities.append(free_faces.get_and_step());

/*
  if (!keep_old )
  {
    for( ii = body_list.size(); ii > 0; ii--)
    {
      Body* body_ptr = body_list.get_and_step() ;
      GeometryModifyTool::instance()->delete_Body(body_ptr) ;
    }
    for ( ii = free_faces.size(); ii > 0; ii-- )
    {
      RefFace *ref_ptr = free_faces.get_and_step();
      RefEntity *ref_entity_ptr = CAST_TO(ref_ptr, RefEntity);
      GeometryModifyTool::instance()->delete_RefEntity(ref_entity_ptr);
    }
  }
*/

  BODY *new_master = master;
#ifdef ACIS_HEALER
//  BODY *new_master = master;
  if ( heal )
  {

      // Now heal the combined body
    PRINT_INFO(" Healing the combined volume...\n");
    if( AcisHealerTool::instance()->init_BODY_for_healing( master ) == CUBIT_SUCCESS )
    {
      int percent_before, percent_after, number_splines_simplified;
      if( AcisHealerTool::instance()->heal_BODY( master, percent_before,
                                                 percent_after, number_splines_simplified ) == CUBIT_FAILURE )
        PRINT_ERROR( "Error healing the combined volume\n" );
      else
        PRINT_INFO( "Healed the combined volume to %d%% good geometry.\n", percent_after );
      AcisHealerTool::instance()->end_BODY_for_healing( master );
    }
  }
#else
  PRINT_WARNING("Healer was not included with this version.\n"
                " Volume healing is usually required for creating\n"
                " volumes out of surfaces.\n");
#endif

  //Gather all the FACEs from resultant BODY
  ENTITY_LIST tmp_FACEs;
  api_get_faces( new_master, tmp_FACEs );
  FACE **tmp_FACE_array = new FACE*[tmp_FACEs.count()];
  tmp_FACEs.init();
  ENTITY *tmp_ent = NULL;
  while( (tmp_ent = tmp_FACEs.next()) != NULL)
    tmp_FACE_array[ii++] = static_cast<FACE*>(tmp_ent);

  //lets make sure that there are not inappropriate intersections between FACEs
  ENTITY_LIST insane_ents;
  logical bad_ints = true;
  api_check_list_ff_ints( tmp_FACEs.count(), tmp_FACE_array, &insane_ents, bad_ints );

  delete [] tmp_FACE_array;

  if( insane_ents.count() || bad_ints )
  {
    PRINT_INFO("Bad surface-surface intersections\n");
    insane_ents.init();
    while( (tmp_ent=insane_ents.next()) != NULL)
    {
      ERROR_ENTITY *error_ent = static_cast<ERROR_ENTITY*>(tmp_ent);
      PRINT_ERROR("ACIS error number %d %s\n", error_ent->get_error_id(),
                              find_err_mess (error_ent->get_error_id()));
    }
    api_delent( new_master );
    return CUBIT_FAILURE;
  }

  //separate the bodies so that you have no multi-volume bodies
  ENTITY_LIST separated_bodies;
  BODY **new_BODY_list;
  int n_bodies;
  api_separate_body( new_master, n_bodies, new_BODY_list );

  DLIList<BODY*> new_BODIES;
  int i;
  for( i=0; i<n_bodies; i++ )
    new_BODIES.append( new_BODY_list[i] );

  CubitStatus success =
    get_new_Body(old_entities, new_BODIES, new_bodies, keep_old);

  if ( new_bodies.size() == 0 )
  {
    PRINT_ERROR("Problems creating new volume.\n" );
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}

// Currently no command line to this
CubitStatus AcisModifyEngine::webcut_across_translate(
                                                 DLIList<BodySM*>& body_list,
                                                 Surface* top_surf_ptr,
                                                 Surface* bottom_surf_ptr,
                                                 DLIList<BodySM*>& ,
                                                 bool  ) const
{
   CubitVector top_orig, top_norm, bottom_orig, bottom_norm;
   double top_orig_pnt[3], top_norm_vec[3], bottom_orig_pnt[3], bottom_norm_vec[3];
   //CubitVector top_sweep_vector, bottom_sweep_vector;
   double distance;
   int i, ii;
   BodySM* body_ptr;
   Surface* ref_face_ptr;
   FACE* FACE_ptr;
   FACE* copied_FACE_ptr;
   ENTITY *entity_ptr;
   DLIList<BODY*> tool_BODY_list;
   //int volume_is_negative;
   AgtSide side;
   AgtOrientation orien;
   outcome result;

   AnalyticGeometryTool* agt = AnalyticGeometryTool::instance();

   // Find top plane
   if( top_surf_ptr->get_point_normal( top_orig, top_norm ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Surface is not planar\n" );
      return CUBIT_FAILURE;
   }

   if( bottom_surf_ptr->get_point_normal( bottom_orig, bottom_norm ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Surface is not planar\n" );
      return CUBIT_FAILURE;
   }

   top_orig.get_xyz( top_orig_pnt ); top_norm.get_xyz( top_norm_vec );
   bottom_orig.get_xyz( bottom_orig_pnt ); bottom_norm.get_xyz( bottom_norm_vec );

   // Check distance between planes
   distance = agt->dist_pln_pln( top_orig_pnt, top_norm_vec, bottom_orig_pnt, bottom_norm_vec,
                                 &side, &orien );

   if( agt->dbl_eq( distance, 0.0 ) )
   {
      PRINT_ERROR( "Distance between planes is 0.0\n" );
      return CUBIT_FAILURE;
   }

   double top_bott[3];
   top_bott[0] = top_norm.x();
   top_bott[1] = top_norm.y();
   top_bott[2] = top_norm.z();
   if( side == AGT_NEG_SIDE ) {
      agt->reverse_vec( top_bott, top_bott );
   }

   // Top plane
   for( i=0; i<body_list.size(); i++ )
   {
      body_ptr = body_list.get_and_step();

      DLIList<Surface*> ref_face_list;
      AcisQueryEngine::instance()->get_surfs_on_plane( body_ptr, top_orig_pnt, top_norm_vec, ref_face_list );

      if( ref_face_list.size() == 0 )
         continue;

      for( ii=0; ii<ref_face_list.size(); ii++ )
      {
         ref_face_ptr = ref_face_list.get_and_step();

         // Copy each surface, reverse and extrude
         FACE_ptr = AcisQueryEngine::get_FACE( ref_face_ptr );

         if( FACE_ptr == NULL )
         {
            return CUBIT_FAILURE;
         }

         result = api_copy_entity_contents( (ENTITY *)FACE_ptr, entity_ptr);
         if( !result.ok() )
         {
            PRINT_ERROR( "Unable to copy face off volume\n" );
            continue;
         }

         copied_FACE_ptr = (FACE *)entity_ptr;
         ATTRIB_CUBIT_OWNER::remove_cubit_owner( (ENTITY *)copied_FACE_ptr,
                                                 CUBIT_TRUE );

         // Create a sheet body from this free-face
         BODY* sheet_BODY_ptr;
         result = api_sheet_from_ff( 1, &copied_FACE_ptr, sheet_BODY_ptr );
         if (!result.ok())
         {
            PRINT_ERROR("Error in api_sheet_from_ff\n");
            AcisQueryEngine::instance()->ACIS_API_error(result);
            PRINT_ERROR( "Unable to create sheet volume from face\n" );
            continue;
         }
         if( sheet_BODY_ptr == NULL )
         {
            PRINT_ERROR( "Unable to create sheet volume from face\n" );
            continue;
         }

         DLIList<FACE*> sheet_FACE_list;
         AcisQueryEngine::instance()->get_FACEs( (ENTITY *)sheet_BODY_ptr, sheet_FACE_list );
         if( sheet_FACE_list.size() != 1 )
         {
            PRINT_ERROR( "Sheet doesn't have one face\n" );
            continue;
         }
         FACE* sheet_FACE_ptr = sheet_FACE_list.get();

         // Get underlying plane normal.  The FACE normal might be reversed
         // of this.
         surface const* acis_surface = &(sheet_FACE_ptr->geometry()->equation());
         SPAposition sheet_surf_orig;
         SPAunit_vector sheet_surf_norm;
         acis_surface->planar( sheet_surf_orig, sheet_surf_norm );
         double sheet_norm[3];
         sheet_norm[0] = sheet_surf_norm.x();
         sheet_norm[1] = sheet_surf_norm.y();
         sheet_norm[2] = sheet_surf_norm.z();

         // Get sheet normal in same direction as ACIS FACE.  The direction of the
         // ACIS face is the direction of sweep.
         if( sheet_FACE_ptr->sense() == REVERSED )
            agt->reverse_vec( sheet_norm, sheet_norm );

         // Properties of dot product
         // angle acute if dot product > 0
         // angle obtuse if dot product < 0
         // angle 90 deg if dot product = 0

         if( agt->dot_vec( sheet_norm, top_bott ) < 0.0 )
         {
            //PRINT_INFO( "Reversing volume...\n" );
            api_reverse_body( sheet_BODY_ptr );
         }

         result = api_sw_face_norm( sheet_FACE_ptr, FORWARD, distance+distance*(.50), 0, 0 );

         if (!result.ok())
         {
            PRINT_ERROR("Error during ACIS sweep operation, api_sw_face_norm\n");
            AcisQueryEngine::instance()->ACIS_API_error(result);
            continue;
         }

         BODY* tool_BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY(sheet_FACE_ptr);

         if (is_BODY(tool_BODY_ptr))
         {
            result = api_body_to_2d(tool_BODY_ptr);
            if (!result.ok())
              PRINT_WARNING("Couldn't turn a tool into a 3D sheet\n" );

            tool_BODY_list.append( tool_BODY_ptr );

            // Move the body slightly so it overlaps
            // Create a SPAvector representing the offsets
            double fact = 0.35;
            SPAvector translation_vector( -top_bott[0]*distance*fact, -top_bott[1]*distance*fact, -top_bott[2]*distance*fact );
            SPAtransf tran_vec = translate_transf (translation_vector);

            // Concatenate the input transformation to the existing one
            result = api_apply_transf ( tool_BODY_ptr, tran_vec );
            if (!result.ok())
            {
               AcisQueryEngine::instance()->ACIS_API_error ( result, "transforming body" );
               continue;
            }
            TRANSFORM* identity = new TRANSFORM(scale_transf(1.0));
            result = api_change_body_trans ( tool_BODY_ptr, identity, FALSE );
            identity->lose();
            if (!result.ok())
            {
               AcisQueryEngine::instance()->ACIS_API_error ( result, "transforming body" );
               continue;

            }

         }

         else
         {
            PRINT_ERROR( "Unable to create tool from surface\n" );
         }
      }

      PRINT_INFO( "Got %d tools\n from the top\n", tool_BODY_list.size() );

      // For now just create tools
      for( ii=0; ii<tool_BODY_list.size(); ii++ )
      {
         BODY *temp_tool_BODY_ptr = tool_BODY_list.get_and_step();

         /*BodySM *new_body = */
         AcisQueryEngine::instance()->populate_topology_bridges(temp_tool_BODY_ptr);
      }

   }

   return CUBIT_SUCCESS;
}

// Imprints curves to bodies.
CubitStatus AcisModifyEngine::imprint( DLIList<BodySM*> &body_list,
                                       DLIList<Curve*> &ref_edge_list,
                                       DLIList<BodySM*>& new_body_list,
                                       bool keep_old_body,
                                       bool show_messages) const
{
   int i, j;
   outcome result;
   CubitStatus status = CUBIT_SUCCESS;
   CubitBoolean imprint_worked = CUBIT_FALSE;

   Curve* ref_edge_ptr;
   BodySM *body_ptr;

   BODY *BODY_ptr;
   BODY *copied_BODY_ptr;
   BODY *wire_BODY_ptr;

   EDGE *EDGE_ptr;
   EDGE *copied_EDGE_ptr;

   DLIList<BODY*> ewire_list;

   bool delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

   // Make wires which can be used to imprint onto the bodies
   ref_edge_list.reset();
   for( i=0; i<ref_edge_list.size(); i++ )
   {
      ref_edge_ptr = ref_edge_list.get_and_step();

      // Get the EDGE pointer
      EDGE_ptr = AcisQueryEngine::get_EDGE( ref_edge_ptr );

      if( EDGE_ptr == NULL )
      {
         return CUBIT_FAILURE;
      }

      // Copy the EDGE.  I've found that the ewire will reference
      // the original EDGE and that EDGE can get split later when
      // doing imprints.  Also note that the copied EDGE will get
      // deleted when deleting the ewire.
      result = api_edge( EDGE_ptr, copied_EDGE_ptr );
      if (!result.ok())
      {
         PRINT_ERROR( "Unable to copy curve; aborting\n" );
         AcisQueryEngine::instance()->ACIS_API_error (result);
         return CUBIT_FAILURE;
      }

      // Remove the CUBIT owner attribute from the EDGE
      ATTRIB_CUBIT_OWNER::remove_cubit_owner( (ENTITY *)copied_EDGE_ptr,
                                              CUBIT_TRUE );

      // Make an ewire from this EDGE
      result = api_make_ewire( 1, &copied_EDGE_ptr, wire_BODY_ptr );
      if (!result.ok())
      {
         AcisQueryEngine::instance()->ACIS_API_error (result);
         PRINT_ERROR( "Unable to make ACIS WIRE from curve; aborting\n" );
         for( j=0; j<ewire_list.size(); j++ )
            api_delent( ewire_list.get_and_step() );
         return CUBIT_FAILURE;
      }

      ewire_list.append( wire_BODY_ptr );
   }

   // Loop on bodies, imprinting each with all the ewires
   for( i=0; i<body_list.size(); i++ )
   {
      body_ptr = body_list.get_and_step();

      BODY_ptr = AcisQueryEngine::get_BODY(body_ptr);
      if( BODY_ptr == NULL )
      {
         PRINT_ERROR( "Unable to find ACIS BODY from volume; aborting\n" );
         for( j=0; j<ewire_list.size(); j++ )
            api_delent( ewire_list.get_and_step() );
         return CUBIT_FAILURE;
      }

      // Copy the body before working with it
      copied_BODY_ptr = copy_BODY( BODY_ptr, delete_attribs );

      // Imprint the wire-bodies onto the BODY
      api_set_int_option("all_free_edges", TRUE );
      ewire_list.reset();
      ref_edge_list.reset();
      for( j=0; j<ewire_list.size(); j++ )
      {
         wire_BODY_ptr = ewire_list.get_and_step();
         ref_edge_ptr = ref_edge_list.get_and_step();
         BODY* temp_BODY_ptr = copy_BODY( copied_BODY_ptr, delete_attribs );

         result = api_imprint( wire_BODY_ptr, copied_BODY_ptr );

         //Gets rid of sliver curves/surfaces that could get produced by the imprint.
         AcisModifyEngine::instance()->cleanup_slivers( copied_BODY_ptr );

         if( !result.ok() && result.error_number()!=200 )
         {
           if ( show_messages )
           {
             PRINT_ERROR( "problem imprinting curve onto volume\n" );
             AcisQueryEngine::instance()->ACIS_API_error(result);
           }
             //revert back to the old body.
           api_delent(copied_BODY_ptr);
           copied_BODY_ptr = temp_BODY_ptr;
           status = CUBIT_FAILURE;
         }
         else
           imprint_worked = CUBIT_TRUE;
      }
      if( imprint_worked )
      {
          api_set_int_option("all_free_edges", FALSE );

          DLIList<EDGE*> new_edges = AcisModifyEngine::instance()->find_new_EDGES(copied_BODY_ptr);

          BodySM* new_body_ptr = get_new_Body( body_ptr, BODY_ptr, copied_BODY_ptr,
              keep_old_body, CUBIT_TRUE );

          if( new_body_ptr )
          {
              // Add a imprint feature to the topo edges
              for (int edge_count = new_edges.size(); edge_count--; ) 
              {
                  CubitSimpleAttrib *tmp_attrib = new CubitSimpleAttrib( "SOURCE_FEATURE", "IMPRINT" );
                  new ATTRIB_SNL_SIMPLE( new_edges[edge_count], tmp_attrib );
                  delete tmp_attrib;
              }
              new_body_list.append(new_body_ptr);
          }
      }
   }

   // Free memory
   for( i=0; i<ewire_list.size(); i++ )
      api_delent( ewire_list.get_and_step() );
   if ( imprint_worked && status != CUBIT_SUCCESS )
     status = CUBIT_SUCCESS;
   return status;
}

// This is a special case - needed because if a curve hangs over
// several surfaces on a body and the user wants to split only one
// of the surfaces.  This method sometimes works when function to
// imprint a curve to a body fails.  The surfaces must be part of
// a body, but the curves just have to be valid ACIS EDGEs.
CubitStatus AcisModifyEngine::imprint( DLIList<Surface*> &surface_list,
                                       DLIList<Curve*> &curve_list,
                                       DLIList<BodySM*>& new_body_list,
                                       bool keep_old_body ) const
{
  int i;

  outcome result;

  // Copy the incoming surface_list since we will be pulling
  // surfaces out of it.
  DLIList<SurfaceACIS*> copied_surface_list(surface_list.size());
  Curve* curve_ptr;

  BodySM *body_ptr;
  BODY *BODY_ptr;

  EDGE *EDGE_ptr;
  EDGE *copied_EDGE_ptr;

  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  CAST_LIST( surface_list, copied_surface_list, SurfaceACIS );
  if (surface_list.size() != copied_surface_list.size())
  {
    PRINT_ERROR("Non-ACIS Surface at %s:%d\n", __FILE__, __LINE__ );
    return CUBIT_FAILURE;
  }

  // Need to copy all of the EDGE's before doing anything, if they are
  // part of a body being imprinted on they could be invalidated before
  // subsequent imprints.
  DLIList<EDGE*> EDGE_list;
  curve_list.reset();
  for( i=curve_list.size(); i--; )
  {
    curve_ptr = curve_list.get_and_step();

    // Get the EDGE pointer
    EDGE_ptr = AcisQueryEngine::get_EDGE( curve_ptr );

    if( EDGE_ptr == NULL )
    {
      PRINT_ERROR("Non-ACIS Curve at %s:%d\n", __FILE__, __LINE__ );
      return CUBIT_FAILURE;
    }

    // Copy the EDGE.
    result = api_edge( EDGE_ptr, copied_EDGE_ptr );
    if (!result.ok())
    {
      AcisQueryEngine::instance()->ACIS_API_error (result);
      PRINT_ERROR( "Unable to copy curve; aborting.\n" );
      return CUBIT_FAILURE;
    }

    // Remove the CUBIT owner attribute from the EDGE
    ATTRIB_CUBIT_OWNER::remove_cubit_owner( (ENTITY *)copied_EDGE_ptr, CUBIT_TRUE );

    // Add the EDGE to the list
    EDGE_list.append( copied_EDGE_ptr );
  }

  // Loop on FACEs.  We will work on imprinting surfaces from one body at a time.
  copied_surface_list.reset();
  while( copied_surface_list.size() )
  {
    DLIList<FACE*> imprint_FACE_list;
    BODY *copied_BODY_ptr;
    if( get_copied_FACES_of_body( copied_surface_list, imprint_FACE_list,
                                  copied_BODY_ptr ) == CUBIT_FAILURE )
      continue;

    // Get original Body and BODY
    body_ptr = AcisQueryEngine::instance()->get_body_sm_of_ENTITY( copied_BODY_ptr );
    BODY_ptr = AcisQueryEngine::get_BODY(body_ptr);

    // Now cleanout the owner attributes from the copied BODY, if required
    if( delete_attribs )
      AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);

    if( imprint( copied_BODY_ptr, imprint_FACE_list, EDGE_list ) == CUBIT_FAILURE )
    {
      // This only fails if there is a serious error - delete the copied Body
        api_delent( copied_BODY_ptr );
        continue;
    }

    DLIList<EDGE*> new_edges = AcisModifyEngine::instance()->find_new_EDGES(copied_BODY_ptr);

    BodySM* new_body_ptr = get_new_Body( body_ptr, BODY_ptr, copied_BODY_ptr,
        keep_old_body, CUBIT_TRUE );

    if( new_body_ptr )
    {
        // Add a imprint feature to the topo edges
        for (int edge_count = new_edges.size(); edge_count--; ) 
        {
            CubitSimpleAttrib *tmp_attrib = new CubitSimpleAttrib( "SOURCE_FEATURE", "IMPRINT" );
            new ATTRIB_SNL_SIMPLE( new_edges[edge_count], tmp_attrib );
            delete tmp_attrib;
        }

        new_body_list.append( new_body_ptr );
    }
  }

  // Free the memory allocated when copying edges
  for( i=0; i<EDGE_list.size(); i++ )
  {
    copied_EDGE_ptr = EDGE_list.get_and_step();
    api_delent( copied_EDGE_ptr );
  }

  return CUBIT_SUCCESS;
}

// Imprint EDGEs to select FACEs.  FACEs must be from same BODY.  Input FACEs
// are modified but EDGEs are not.
CubitStatus AcisModifyEngine::imprint( BODY *BODY_ptr,
                                       DLIList<FACE*> &FACE_list,
                                       DLIList<EDGE*> &EDGE_list ) const
{
  // Imprint the edges to the surfaces on this body
  outcome result;
  EDGE *EDGE_ptr;
  FACE *FACE_ptr;
  EDGE_list.reset();
  FACE_list.reset();
  int i, j;
  for( i=EDGE_list.size(); i--; )
  {
    // Imprint this edge selectively to the given faces on this body
    EDGE_ptr = EDGE_list.get_and_step();

    // Track the existing FACEs before imprinting this BODY with the edge.
    // Needed so we know which new surfaces result from the imprint.  These
    // new surfaces will need to be imprinted with subsequent edges.
    DLIList<FACE*> FACES_BEFORE;
    DLIList<FACE*> FACES_AFTER;
    AcisQueryEngine::instance()->get_FACEs( (ENTITY *)BODY_ptr, FACES_BEFORE );

    FACE_list.reset();
    for( j=FACE_list.size(); j--; )
    {
      // Get the FACE to imprint to
      FACE_ptr = FACE_list.get_and_step();

      if( imprint( BODY_ptr, FACE_ptr, EDGE_ptr ) == CUBIT_FAILURE )
      {
        PRINT_ERROR( "Unable to imprint EDGE to FACE\n" );
        continue;
      }

      // Check against an assumption that if invalid, will
      // cause a crash.  The assumption is that the imprinted
      // FACE_ptr will remain valid after each imprint.
      FACES_AFTER.clean_out();
      AcisQueryEngine::instance()->get_FACEs( (ENTITY *)BODY_ptr, FACES_AFTER );
      if( !FACES_AFTER.is_in_list( FACE_ptr ) )
      {
        PRINT_ERROR( "serious problem in imprint - an imprinted surface disappeared!\n"
          "       Model restored - operation aborted\n" );
        return CUBIT_FAILURE;
      }

     }
     // End of looping over FACES

     // Find which new surfaces were created - add these to the
     // FACE_list so that additional edges will imprint
     // on these as well.
     DLIList<FACE*> new_FACES;
     FACE* after_FACE_ptr;
     for( j=FACES_AFTER.size(); j--; )
     {
       after_FACE_ptr = FACES_AFTER.get_and_step();

       if( !FACES_BEFORE.is_in_list( after_FACE_ptr ) )
         new_FACES.append(after_FACE_ptr);
     }
     FACE_list += new_FACES;
  }

  return CUBIT_SUCCESS;
}

CubitStatus AcisModifyEngine::imprint( BODY *BODY_ptr,
                                       FACE *FACE_ptr,
                                       EDGE *EDGE_ptr ) const
{
  outcome result;

  BODY* tbody = NULL;
  BODY* bbody = NULL;    // Blank body
  FACE* bface = NULL;    // Blank face
  surface* tsurf = NULL; // Tool surface
  EDGE* edge = NULL;

  //EXCEPTION_TRY

  bbody = BODY_ptr;
  bface = FACE_ptr;
  edge = EDGE_ptr;

  // Get the geometry of the blank face.
  surface const& bsurf = bface->geometry()->equation();

  // Create an empty tool body
  tbody = ACIS_NEW BODY( (LUMP*) NULL );

  // Create the tool body which consists of a single face
  // with a "cross surface" geometric support and no loops.
  const curve &temp_cur = edge->geometry()->equation();
  SPAinterval cur_range = temp_cur.param_range();
  SPAinterval edge_range = edge->param_range();
  if (edge->sense() == REVERSED)
    edge_range = -edge_range;
  if (edge_range == cur_range)
    tsurf = cross_surface( edge->geometry()->equation(),
    bsurf, *(pcurve*) NULL_REF);
  else
  {
    // Create using a subsetted copy of the curve
    curve* subset_cur = temp_cur.subset(edge_range&cur_range);
    tsurf = cross_surface( *subset_cur, bsurf, *(pcurve*) NULL_REF);
    ACIS_DELETE subset_cur;
  }
  if( !tsurf )
  {
    PRINT_ERROR( "Unable to make imprint surface on body\n" );
    return CUBIT_FAILURE;
  }
  FACE* tface = ACIS_NEW FACE( NULL, NULL, make_surface( *tsurf ), FORWARD );
  ACIS_DELETE tsurf; // SRS Added 10-29-99
  tbody->set_lump( ACIS_NEW LUMP( ACIS_NEW SHELL( tface, NULL, NULL ), NULL ) );

  // Now, imprint the two bodies and throw away the tool body. The
  // resultant blank body will have the edge imprinted on the
  // selected face.
  EDGE* ssi_edges[1];
  FACE* tfaces[1];
  FACE* bfaces[1];

  ssi_edges[0] = edge;
  tfaces[0] = tface;
  bfaces[0] = bface;

  result = api_boolean_start( tbody, bbody );

  if( result.ok() )
  {
    // If the edge is known to lie on the face, some expensive
    // checking can be ommited by changing the last argument from
    // TRUE to FALSE.
    result = api_update_intersection( tface, *(SPAtransf *)NULL_REF, bface,
      *(SPAtransf *)NULL_REF, 1, ssi_edges, TRUE);

    if( result.ok() )
    {
      result = api_selectively_intersect( 1, tfaces, bfaces );
      if( result.ok() )
      {
        result = api_imprint_complete( tbody, bbody );
        if( !result.ok() )
        {
          PRINT_ERROR( "problem in completing imprint for curve on surface\n" );
          AcisQueryEngine::instance()->ACIS_API_error(result);
          api_delent( tbody );
          return CUBIT_FAILURE;
        }
      }
      else
      {
        PRINT_ERROR( "problem in selective intersect for curve on surface\n" );
        AcisQueryEngine::instance()->ACIS_API_error(result);
        api_delent( tbody );
        return CUBIT_FAILURE;
      }
    }
    else if( result.error_number()!=200 && result.error_number()!=615
      && result.error_number()!=715 )
    {
      // Error 200 can occur when the imprint does nothing (curve that's
      // a surface boundary, curve away from surface).  Errors 615/715 can occur
      // when imprinting a remnant of the same FACE again with the same curve.
      // The orignal FACE_ptr remains after an imprint - it lives as one
      // of the remnants.  Then, when we attempt to imprint the other pieces
      // with the same curves, we are getting error 715.  The same error
      // can occur when trying to use one of the surfaces own edges to
      // imprint with.  This has been reported to Spatial and logged as r6294.
      // It's not fatal so for now we ignore it.

      // We get these all the time, and they don't seem to hurt us.
      api_delent( tbody );
      return CUBIT_SUCCESS;
    }
  }
  else
  {
    PRINT_ERROR( "problem starting boolean for curve on surface\n" );
    AcisQueryEngine::instance()->ACIS_API_error(result);
    api_delent( tbody );
    return CUBIT_FAILURE;
  }

  //EXCEPTION_CATCH(TRUE)

  // Remove temporary objects
  api_delent( tbody );

  //EXCEPTION_END

  return CUBIT_SUCCESS;
}

//- Imprints a list of Surfaces with list of Curves, sorted per
//- Surface (ie., curve_lists_list is same length as surface_list).
//- This version is more efficient than the general-purpose one
//- above, as we know which curves to imprint with which surfaces.
//- All input surfaces must be from the same body.
CubitStatus AcisModifyEngine::imprint( DLIList<Surface*> &surface_list,
                                       DLIList<DLIList<Curve*>*> &curve_lists_list,
                                       BodySM*& new_body,
                                       bool keep_old_body ) const
{
  // Calling code has already ensured that all the Surfaces are from the
  // same body

  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  SurfaceACIS* surf_ptr = dynamic_cast<SurfaceACIS*>(surface_list.get());
  if( surf_ptr == NULL )
    return CUBIT_FAILURE;
  FACE *FACE_ptr = surf_ptr->get_FACE_ptr();
  if (!FACE_ptr)
  {
    PRINT_ERROR("Invalid SurfaceACIS (no FACE ptr) at %s:%d\n",__FILE__,__LINE__);
    return CUBIT_FAILURE;
  }
  BODY *BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY( FACE_ptr );
  if( BODY_ptr == NULL )
    return CUBIT_FAILURE;
  BodySM *body_ptr = AcisQueryEngine::instance()->get_body_sm_of_ENTITY( FACE_ptr );
  if( body_ptr == NULL )
    return CUBIT_FAILURE;

  int i, j;
  EDGE *EDGE_ptr;
  DLIList<DLIList<EDGE*>*> EDGE_lists_list;
  DLIList<EDGE*> *EDGE_list_ptr;
  DLIList<SurfaceACIS*> surf_ACIS_list;

  surface_list.reset();
  curve_lists_list.reset();
  for( i=surface_list.size(); i--; )
  {
    surf_ptr = dynamic_cast<SurfaceACIS*>(surface_list.get_and_step());
    if (!surf_ptr)
    {
      PRINT_ERROR("Non-ACIS surface at %s: %d\n",__FILE__,__LINE__);
      return CUBIT_FAILURE;
    }
    surf_ACIS_list.append( surf_ptr );

    // Make sure there is a FACE associated with this SurfaceACIS
    FACE_ptr = surf_ptr->get_FACE_ptr();
    if (!FACE_ptr)
    {
      PRINT_ERROR("Invalid SurfaceACIS (no FACE ptr) at %s:%d\n",__FILE__,__LINE__);
      return CUBIT_FAILURE;
    }

    EDGE_list_ptr = new DLIList<EDGE*>;
    EDGE_lists_list.append( EDGE_list_ptr );

    DLIList<Curve*> *curve_list_ptr = curve_lists_list.get_and_step();
    curve_list_ptr->reset();
    for( j=curve_list_ptr->size(); j--; )
    {
      Curve *curve_ptr = curve_list_ptr->get_and_step();
      CurveACIS* curve_acis_ptr = dynamic_cast<CurveACIS*>(curve_ptr);
      if (!curve_acis_ptr)
      {
        PRINT_ERROR("Non-ACIS curve at %s: %d\n",__FILE__,__LINE__);
        return CUBIT_FAILURE;
      }
      EDGE_ptr = curve_acis_ptr->get_EDGE_ptr();
      if( EDGE_ptr == NULL )
      {
        PRINT_ERROR("Invalid CurveACIS (no EDGE ptr) at %s:%d\n",__FILE__,__LINE__);
        return CUBIT_FAILURE;
      }
      EDGE_list_ptr->append( EDGE_ptr );
    }
  }

  // Copy the BODY and get the corresponding FACEs on the copied BODY that need
  // to be imprinted
  DLIList<FACE*> copied_FACE_list;
  BODY *copied_BODY_ptr;
  if( get_copied_FACES_of_body( surf_ACIS_list, copied_FACE_list, copied_BODY_ptr )
    == CUBIT_FAILURE )
  {
    return CUBIT_FAILURE;
  }

  // Now cleanout the owner attributes from the copied BODY, if required
  if( delete_attribs )
    AcisQueryEngine::instance()->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);

  // Imprint the temporary EDGEs to the FACEs
  if( imprint( copied_BODY_ptr, copied_FACE_list, EDGE_lists_list ) == CUBIT_FAILURE )
  {
    api_delent( copied_BODY_ptr );
    while( EDGE_lists_list.size() ) delete EDGE_lists_list.pop();
    return CUBIT_FAILURE;
  }
  DLIList<EDGE*> new_edges = AcisModifyEngine::instance()->find_new_EDGES(copied_BODY_ptr);
  new_body = get_new_Body( body_ptr, BODY_ptr, copied_BODY_ptr, false);

  // Add a imprint feature to the topo edges
  for (int edge_count = new_edges.size(); edge_count-- && new_body; ) 
  {
      CubitSimpleAttrib *tmp_attrib = new CubitSimpleAttrib( "SOURCE_FEATURE", "IMPRINT" );
      new ATTRIB_SNL_SIMPLE( new_edges[edge_count], tmp_attrib );
      delete tmp_attrib;
  }

  // Free memory allocated in EDGE_lists_list
  while( EDGE_lists_list.size() ) delete EDGE_lists_list.pop();

  return new_body ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

DLIList<EDGE*> AcisModifyEngine::find_new_EDGES(BODY *copied_BODY_ptr)
{
    DLIList<EDGE*> return_list;
    ENTITY_LIST edge_list;
    api_get_edges(copied_BODY_ptr, edge_list);
    AcisBridge *ab_edge_ptr = NULL;
    for(int i=0; i<edge_list.count(); i++)
    {
        ab_edge_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( edge_list[i] );
        if( !ab_edge_ptr )
        {
            return_list.append((EDGE*)edge_list[i]);
            //AcisDrawTool::instance()->draw_EDGE( (EDGE*)edge_list[i], CUBIT_BLUE );
        }
    }
    // GfxDebug::flush();
    return return_list;
}

// This version is more efficient if it is known exactly which EDGEs
// need to be imprinted with which FACEs (ie., it does not attempt to
// imprint all input EDGEs to all input FACEs - only those in corresponding
// list positions)
CubitStatus AcisModifyEngine::imprint( BODY *BODY_ptr,
                                       DLIList<FACE*> &FACE_list,
                                       DLIList<DLIList<EDGE*>*> &EDGE_lists_list ) const
{
  // Imprint the edges to the surfaces on this body
  outcome result;
  FACE *FACE_ptr = 0;
  EDGE *EDGE_ptr = 0;

  int i, j, k, m;

  DLIList<EDGE*> *EDGE_list_ptr;
  EDGE_lists_list.reset();
  FACE_list.reset();
  for( i=EDGE_lists_list.size(); i--; )
  {
    EDGE_list_ptr = EDGE_lists_list.get_and_step();
    FACE_ptr = FACE_list.get_and_step();

    DLIList<FACE*> temp_FACE_list;
    temp_FACE_list.append( FACE_ptr );

    EDGE_list_ptr->reset();
    for( j=EDGE_list_ptr->size(); j--; )
    {
      // Imprint this edge list selectively to the current FACE
      EDGE_ptr = EDGE_list_ptr->get_and_step();

      // Track the existing FACEs before imprinting this BODY with the edge.
      // Needed so we know which new surfaces result from the imprint.  These
      // new surfaces will need to be imprinted with subsequent edges.
      DLIList<FACE*> FACES_BEFORE;
      DLIList<FACE*> FACES_AFTER;

      temp_FACE_list.reset();
      for( k=0; k<temp_FACE_list.size(); k++ )
      {
        FACE_ptr = temp_FACE_list.get_and_step();

        FACES_BEFORE.clean_out();
        AcisQueryEngine::instance()->get_FACEs( (ENTITY *)BODY_ptr, FACES_BEFORE );

        if( imprint( BODY_ptr, FACE_ptr, EDGE_ptr ) == CUBIT_FAILURE )
        {
          PRINT_ERROR( "Unable to imprint EDGE to FACE\n" );
          continue;
        }

        // Check against an assumption that if invalid, will
        // cause a crash.  The assumption is that the imprinted
        // FACE_ptr will remain valid after each imprint.
        FACES_AFTER.clean_out();
        AcisQueryEngine::instance()->get_FACEs( (ENTITY *)BODY_ptr, FACES_AFTER );
        if( !FACES_AFTER.is_in_list( FACE_ptr ) )
        {
          PRINT_ERROR( "serious problem in imprint - an imprinted surface disappeared!\n"
            "       Model restored - operation aborted\n" );
          return CUBIT_FAILURE;
        }

        // Find which new surfaces were created - add these to the
        // temp_FACE_list so that additional edges will imprint
        // on these as well.
        DLIList<FACE*> new_FACES;
        FACE* after_FACE_ptr;
        for( m=FACES_AFTER.size(); m--; )
        {
          after_FACE_ptr = FACES_AFTER.get_and_step();

          if( !FACES_BEFORE.is_in_list( after_FACE_ptr ) )
            new_FACES.append(after_FACE_ptr);
        }
        temp_FACE_list += new_FACES;

        // Set list back to correct location
        temp_FACE_list.move_to( FACE_ptr );
        temp_FACE_list.step();
     } // End of looping over temp_FACE_list
    }
  }
  return CUBIT_SUCCESS;
}

// Imprints locations to bodies (for splitting curves or putting hardpoints
// on surfaces).
CubitStatus AcisModifyEngine::imprint( DLIList<BodySM*> &body_list,
                                       DLIList<CubitVector*> &vector_list,
                                       DLIList<BodySM*>& new_body_list,
                                       bool keep_old_body,
                                       DLIList<TopologyBridge*> *new_tbs,
                                       DLIList<TopologyBridge*> *att_tbs ) const
{
   // Until Spatial fixes imprinting with vertices, implement a hack
   // (Steve Storm, Caterpillar Inc.)
   int i, j, k;
   outcome result;
   double tol = 1e-3;
   //CubitStatus status = CUBIT_SUCCESS;

   CubitVector* vector_ptr;
   BodySM *body_ptr;

   VERTEX *VERTEX_ptr;
   EDGE *EDGE_ptr;
   FACE *FACE_ptr;

   BODY *BODY_ptr;
   BODY *copied_BODY_ptr;

   //CubitVector vertex_coords;

   bool delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

   // Set an option for imprinting
   api_set_int_option("all_free_edges", TRUE );

  // new_ENTITIES will hold new entities that were created by 
  // the imprint operation.  att_ENTITIES will hold entities
  // from the original bodies that had composite attributes
  // on them.  At the end we will compare the new entities
  // with the ones that had attributes and if we find matches
  // we will put a composite attribute on the new entity.
  // This is how we are handling the imprinting of bodies
  // with composites.
  DLIList<ENTITY*> new_ENTITIES;
  DLIList<ENTITY*> att_ENTITIES;

   body_list.reset();
   for( i=0; i<body_list.size(); i++ )
   {
      body_ptr = body_list.get_and_step();

      BODY_ptr = AcisQueryEngine::get_BODY( body_ptr );
      if( BODY_ptr == NULL )
         return CUBIT_FAILURE;

      // Copy the body before working with it
      copied_BODY_ptr = copy_BODY( BODY_ptr, delete_attribs );

      // Get all the vertices of the BODY
      DLIList<VERTEX*> VERTEX_list;
      AcisQueryEngine::instance()->get_VERTICEs( (ENTITY*)BODY_ptr, VERTEX_list );

      // Get all the FACES of the *original* BODY
      DLIList<FACE*> FACE_list;
      AcisQueryEngine::instance()->get_FACEs( (ENTITY*)BODY_ptr, FACE_list );

      // Loop through the passed-in locations, first checking for coincidence
      // w/vertices, then if on curves, then on faces.
      vector_list.reset();
      for( j=0; j<vector_list.size(); j++ )
      {
         vector_ptr = vector_list.get_and_step();
         SPAposition test_point( vector_ptr->x(), vector_ptr->y(), vector_ptr->z() );

         CubitBoolean on_vertex = CUBIT_FALSE;
         VERTEX_list.reset();
         for( k=0; k<VERTEX_list.size(); k++ )
         {
            VERTEX_ptr = VERTEX_list.get_and_step();

            SPAposition VERTEX_coords = VERTEX_ptr->geometry()->coords();

            if (( fabs(VERTEX_coords.x() - vector_ptr->x()) < tol) &&
               ( fabs(VERTEX_coords.y() - vector_ptr->y()) < tol) &&
               ( fabs(VERTEX_coords.z() - vector_ptr->z()) < tol))
            {
               on_vertex = CUBIT_TRUE;
               break; // Skip - coincident with a vertex
            }
         }

         if( on_vertex )
            continue; // No imprinting necessary

         // Get all the EDGES of the copied BODY
         DLIList<EDGE*> EDGE_list;
         AcisQueryEngine::instance()->get_EDGEs( (ENTITY*)copied_BODY_ptr, EDGE_list );

         // Check the EDGES
         CubitBoolean on_edge = CUBIT_FALSE;
         EDGE_list.reset();
         for( k=0; k<EDGE_list.size(); k++ )
         {
            EDGE_ptr = EDGE_list.get_and_step();

            if (EDGE_ptr->geometry() == NULL)
               continue; // Check next edge

            //const curve* acis_curve = &EDGE_ptr->geometry()->equation();

            SPAtransf ftrans;

            point_edge_containment pe_rel;

            ENTITY *ENTITY_ptr = NULL;
            SPAparameter param;

            pe_rel = sg_point_in_edge( test_point, EDGE_ptr, ftrans, ENTITY_ptr, param );

            if( pe_rel == point_off_edge )
               continue; // Check next edge
            else if( pe_rel == point_on_edge )
            {
              on_edge = CUBIT_TRUE;
              VERTEX *new_vertex = ACIS_NEW VERTEX( ACIS_NEW APOINT( test_point) );

              API_BEGIN;
              sg_split_edge_at_vertex( EDGE_ptr, new_vertex );  
              API_END;

              if( !result.ok() && result.error_number()!=200 )
              {
                AcisQueryEngine::instance()->ACIS_API_error(result);
                PRINT_ERROR( "problem imprinting curve onto volume\n");
              }

              break;
            }
         }

         if( on_edge )
            continue; // Already imprinted onto the body - check next vertex

         // Check the FACES
         FACE_list.reset();
         for( k=0; k<FACE_list.size(); k++ )
         {
            FACE_ptr = FACE_list.get_and_step();

            if (FACE_ptr->geometry() == NULL)
               continue; // Check next face

            SPAtransf ftrans;

            point_face_containment pf_rel = point_in_face( test_point, FACE_ptr, ftrans );
            if( pf_rel == point_outside_face )
               continue; // Check next face
            else if( pf_rel == point_inside_face )
            {
               // Create a short curve normal to the surface at the vertex and imprint it
               SPAunit_vector norm;
               norm = FACE_ptr->geometry()->equation().point_normal( test_point );

               if ( FACE_ptr->sense() == REVERSED )
                  norm = -norm;

               // Find the end point
               SPAposition end_point;
               end_point.set_x( test_point.x() + (tol * norm.x()) );
               end_point.set_y( test_point.y() + (tol * norm.y()) );
               end_point.set_z( test_point.z() + (tol * norm.z()) );

               // Create the curve
               EDGE* new_EDGE_ptr = NULL;
               result = api_mk_ed_line ( test_point, end_point, new_EDGE_ptr );
               if (!result.ok() || new_EDGE_ptr == NULL)
               {
                  AcisQueryEngine::instance()->ACIS_API_error (result);
                  PRINT_WARNING("Unable to make short normal curve for imprint operation\n" );
                  continue;
               }

               // Create a wire from the short curve
               BODY *wire_BODY_ptr = NULL;
               result = api_make_ewire( 1, &new_EDGE_ptr, wire_BODY_ptr );
               if (!result.ok())
               {
                  AcisQueryEngine::instance()->ACIS_API_error (result);
                  PRINT_WARNING( "Unable to make ACIS WIRE from short imprint curve\n" );

                  if( wire_BODY_ptr )
                     api_delent( wire_BODY_ptr );
                  else
                     api_delent( new_EDGE_ptr );
                  continue;
               }

               // Imprint the wire-body onto the BODY
               result = api_imprint( wire_BODY_ptr, copied_BODY_ptr );

               if( !result.ok() && result.error_number()!=200 )
               {
                  AcisQueryEngine::instance()->ACIS_API_error(result);
                  PRINT_ERROR( "problem imprinting short curve onto volume\n");
               }

               api_delent( wire_BODY_ptr );

               break; // Done imprinting with this vertex

            }
            else if( pf_rel == point_boundary_face )
            {
               PRINT_WARNING( "Algorithm error - point found on boundary of surface - please report this\n" );
               continue;
            }
            else if( pf_rel == point_unknown_face )
               continue; // Check next face
         }
      }

      AcisModifyEngine::instance()->cleanup_slivers( copied_BODY_ptr );

      // Get the new entities (just edges for now)
      // created by the imprint operation.
      DLIList<ENTITY*> cur_new_ENTITIES;
      if(copied_BODY_ptr && new_tbs)
        get_new_ENTITIES(copied_BODY_ptr, cur_new_ENTITIES);

      BodySM* new_body_ptr = get_new_Body( body_ptr, BODY_ptr,
                               copied_BODY_ptr, keep_old_body, CUBIT_TRUE );
      if (new_body_ptr)
      {
        new_body_list.append(new_body_ptr);

        // Add the new entities to the list we are
        // accumulating.
        if(new_tbs)
          new_ENTITIES += cur_new_ENTITIES;

        // Get entities with composite attributes.
        if(copied_BODY_ptr && att_tbs)
          get_att_ENTITIES(copied_BODY_ptr, att_ENTITIES, "COMPOSITE_GEOM");

      }
      else
      {
        // Get entities with composite attributes.
        if(att_tbs)
          get_att_ENTITIES(BODY_ptr, att_ENTITIES, "COMPOSITE_GEOM");
      }

   }

    // Convert the new_ENTITIES/att_ENTITIES lists into
    // topology bridge lists.
    if(new_tbs)
    {
      for(i=new_ENTITIES.size(); i--;)
      {
        ENTITY *cur_ENT = new_ENTITIES.get_and_step();
        AcisBridge *acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(cur_ENT);
        if(dynamic_cast<TopologyBridge*>(acis_bridge))
          new_tbs->append_unique(dynamic_cast<TopologyBridge*>(acis_bridge));
      }
    }
    if(att_tbs)
    {
      for(i=att_ENTITIES.size(); i--;)
      {
        ENTITY *cur_ENT = att_ENTITIES.get_and_step();
        AcisBridge *acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(cur_ENT);
        if(dynamic_cast<TopologyBridge*>(acis_bridge))
          att_tbs->append_unique(dynamic_cast<TopologyBridge*>(acis_bridge));
      }
    }

   api_set_int_option("all_free_edges", FALSE );

   return CUBIT_SUCCESS;

// When Spatial gets vertex body to body imprinting working, you should be
// able to replace all of the above code with this code.
#if 0
   int i, j;
   outcome result;
   CubitStatus status = CUBIT_SUCCESS;

   RefVertex* ref_vertex_ptr;
   Body *body_ptr;

   BODY *BODY_ptr;
   BODY *copied_BODY_ptr;

   CubitVector vertex_coords;
   BODY* vertex_BODY_ptr;
   DLIList<BODY*> vertex_BODY_list;

   bool delete_attribs =
      (GeometryModifyTool::instance()->new_ids() || keep_old_body);

   // Create a list of BODYs that consist of single points.  These
   // will be used for imprinting onto the BODYs.
   vector_list.reset();
   SPAposition acis_position[1];
   for( i=0; i<vector_list.size(); i++ )
   {
      vector_ptr = vector_list.get_and_step();

      acis_position[0].set_x( vector_ptr->x() );
      acis_position[0].set_y( vector_ptr->y() );
      acis_position[0].set_z( vector_ptr->z() );

      result = api_make_wire( NULL, 1, acis_position, vertex_BODY_ptr );
      if( !result.ok() )
      {
         AcisQueryEngine::instance()->ACIS_API_error (result);
         PRINT_ERROR( "unable to make degenerate wire body from location %f, %f, %f; aborting\n",
                      vector_ptr->x(), vector_ptr->y(), vector_ptr->z() );
         for( j=0; j<vertex_BODY_list.size(); j++ )
            api_delent( vertex_BODY_list.get_and_step() );
         return CUBIT_FAILURE;
      }

      vertex_BODY_list.append( vertex_BODY_ptr );
   }

   // Loop on bodies, imprinting each with all the vertex BODYs.
   for( i=0; i<body_list.size(); i++ )
   {
      body_ptr = body_list.get_and_step();

      BODY_ptr = AcisQueryEngine::get_BODY( body_ptr );
      if( BODY_ptr == NULL )
      {
//NOTE TO JS
         PRINT_ERROR( "Unable to get ACIS BODY from Body %d; aborting\n",
            body_ptr->id() );
         for( j=0; j<vertex_BODY_list.size(); j++ )
            api_delent( vertex_BODY_list.get_and_step() );
         return CUBIT_FAILURE;
      }

      // Copy the body before working with it
      copied_BODY_ptr = copy_BODY( BODY_ptr, delete_attribs );

      // Imprint the wire-bodies onto the BODY
      api_set_int_option("all_free_edges", TRUE );
      vertex_BODY_list.reset();
      vector_list.reset();
      for( i=0; i<vertex_BODY_list.size(); i++ )
      {
         vertex_BODY_ptr = vertex_BODY_list.get_and_step();
         vector_ptr = vector_list.get_and_step();

         result = api_imprint( vertex_BODY_ptr, copied_BODY_ptr );

         if( !result.ok() && result.error_number()!=200 )
         {
//NOTE to JS
            PRINT_ERROR( "problem imprinting location %f, %f, %f to body %d\n",
                         vector_ptr->x(), vector_ptr->y(), vector_ptr->z(),
                         body_ptr->id() );
            AcisQueryEngine::instance()->ACIS_API_error(result);
            status = CUBIT_FAILURE;
         }
      }
      api_set_int_option("all_free_edges", FALSE );

      Body* new_body_ptr = get_new_Body( body_ptr, BODY_ptr, copied_BODY_ptr,
                                         keep_old_body, CUBIT_TRUE );

      if( new_body_ptr!=NULL && new_body_ptr!=body_ptr )
      {
//NOTE TO JS
         PRINT_INFO( "Created new body %d\n", new_body_ptr->id() );
         new_body_list.append( new_body_ptr );
      }
      else if( new_body_ptr!=NULL && new_body_ptr==body_ptr )
      {
//NOTE TO JS
         PRINT_INFO( "Modified body %d\n", body_ptr->id() );
      }
//NOTE TO JS
      else
         PRINT_WARNING( "Body %d was not modified\n", body_ptr->id() );
   }

   // Free memory
   for( i=0; i<vertex_BODY_list.size(); i++ )
      api_delent( vertex_BODY_list.get_and_step() );

   return status;
#endif
}

CubitStatus
AcisModifyEngine::offset_curves( DLIList<Curve*>& ref_edge_list,
                                 DLIList<Curve*>& new_curves,
                                 double offset_distance,
                                 const CubitVector& offset_direction,
                                 int gap_type )
{
   int i;

   Curve* ref_edge_ptr;
   EDGE* EDGE_ptr;
   outcome result;

   // If there's just one curve, and it's linear, all we can do is copy it and move it.
   // SRS Note (1-12-01) - this code should be moved up into GeometryTool...but, we
   // would need to be able to move a curve.  When that capability exists move the
   // code.

   //CubitBoolean move_curve = CUBIT_FALSE;

   if( ref_edge_list.size() == 1 && (offset_direction.x() || offset_direction.y() ||
      offset_direction.z()) )
   {
      ref_edge_ptr = ref_edge_list.get();
      if( AcisQueryEngine::instance()->is_curve_app_straight( ref_edge_ptr ) == CUBIT_TRUE )
      {

         EDGE_ptr = AcisQueryEngine::get_EDGE(ref_edge_ptr);

         if( EDGE_ptr == NULL )
            return CUBIT_FAILURE;

         // Copy and move the EDGE
         EDGE* new_EDGE_ptr2 = copy_and_move_EDGE( EDGE_ptr, offset_direction,
                                                   offset_distance );

         Curve *curve_ptr = AcisQueryEngine::instance()->populate_topology_bridges( new_EDGE_ptr2 );
         new_curves.append(curve_ptr);
         return CUBIT_SUCCESS;
      }
      else
      {
         PRINT_WARNING( "Direction qualifier ignored - only valid for one straight curve\n" );
      }
   }
   else if( offset_direction.x() || offset_direction.y() || offset_direction.z() )
   {
      PRINT_WARNING( "Direction qualifier ignored - only valid for one straight curve\n" );
   }

   // Create a wire from the given EDGEs

   // Offset the RefEdges
   DLIList<EDGE*> old_EDGE_list;
   if( AcisQueryEngine::instance()->get_EDGEs_of_Curves( ref_edge_list, old_EDGE_list )
       == CUBIT_FAILURE )
      return CUBIT_FAILURE;

   DLIList<EDGE*> new_EDGE_list;
   if( offset_EDGES( old_EDGE_list, offset_distance, gap_type, new_EDGE_list ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Unable to offset the curves\n" );
      return CUBIT_FAILURE;
   }

   // Make free curves from the offset EDGEs
   new_EDGE_list.reset();
   for( i=0; i<new_EDGE_list.size(); i++ )
   {
      EDGE_ptr = new_EDGE_list.get_and_step();

      Curve *curve_ptr = AcisQueryEngine::instance()->populate_topology_bridges( EDGE_ptr );
      new_curves.append(curve_ptr);
   }

  return CUBIT_SUCCESS;
}

Curve*
AcisModifyEngine::trim_curve( Curve* trim_curve,
                              const CubitVector& trim_vector,
                              const CubitVector& keep_vector,
                              bool keep_old )
{
  //api_trim_curve(entity_with_ray &eray1, SPAposition* trim_pt, entity_with_ray *peray2);
   //The curve is trimmed to either a SPAposition or at its intersection with another curve. To trim
   //to a SPAposition, specify a SPAposition in trim_pt and NULL for eray2. The curve is trimmed to
   //the normal projection of the SPAposition onto the curve. To trim to the intersection with
   //another curve, specify NULL for trim_pt and an entity_with_ray for eray2. The curve
   //is trimmed to the intersection with the curve given in eray2 that is closest to the ray given
   //in eray2.


   EDGE* EDGE_ptr = AcisQueryEngine::get_EDGE( trim_curve );
   if( EDGE_ptr == NULL )
      return 0;

   ENTITY_LIST edge_VERTS;
   api_get_vertices( EDGE_ptr, edge_VERTS );
   VERTEX *start_vertex = static_cast<VERTEX*>(edge_VERTS[0]);
   VERTEX *end_vertex = static_cast<VERTEX*>(edge_VERTS[1]);
   SPAposition s_vec = start_vertex->geometry()->coords();
   SPAposition e_vec = end_vertex->geometry()->coords();
   CubitVector s_vector( s_vec.x(), s_vec.y(), s_vec.z() );
   CubitVector e_vector( e_vec.x(), e_vec.y(), e_vec.z() );

   int start = 0, end = 0;
   if( keep_vector.distance_between( s_vector ) < 0.001 )
     start = 1;
   else if( keep_vector.distance_between( e_vector ) < 0.001 )
     end = 1;

   if( start || end )
   {
     double trim_param = trim_curve->u_from_position( trim_vector );
     double start_param = trim_curve->start_param();
     double end_param = trim_curve->end_param();

     //if trim_param is not within start and end params, we are extending
     //and we don't need keep_vector
     bool within_u_param_range = true;
     if( start_param < end_param )
     {
       if( trim_param > end_param || trim_param < start_param )
         within_u_param_range = false;
     }
     else
     {
       if( trim_param > start_param || trim_param < end_param )
         within_u_param_range = false;
     }

     if( within_u_param_range )
     {
       double keep_param;
       if( start )
         keep_param = start_param + .95*(trim_param-start_param);
       else
         keep_param = trim_param + .05*(end_param-trim_param);

       CubitVector tmp_keep_vec = keep_vector;
       trim_curve->position_from_u( keep_param, tmp_keep_vec );
     }
   }

   // Copy the edge to trim first
   EDGE* copied_EDGE_ptr;
   outcome result = api_edge( EDGE_ptr, copied_EDGE_ptr );

   SPAposition trim_pos( trim_vector.x(),trim_vector.y(),trim_vector.z() );
   SPAposition* ptrim_pos = &trim_pos;
   entity_with_ray *peray2 = NULL;

   SPAposition keep_pos( keep_vector.x(),keep_vector.y(),keep_vector.z() );
   SPAunit_vector keep_vec( 0.0, 0.0, 0.0 );
   entity_with_ray eray1( copied_EDGE_ptr, keep_pos, keep_vec );

   // Finally, do the real work
   result = api_trim_curve(eray1, ptrim_pos, peray2);

   if( !result.ok() )
   {
     AcisQueryEngine::instance()->ACIS_API_error (result);
     PRINT_ERROR( "Unable to trim curve\n" );
     api_delent( copied_EDGE_ptr );
     return 0;
   }

   // Make a free curve from the trimmed EDGE
   Curve *curve_ptr = AcisQueryEngine::instance()->populate_topology_bridges( copied_EDGE_ptr );
   if (keep_old)
     return curve_ptr;

   // Delete the old curve (trim_curve), if it is a free curve
   DLIList<Surface*> surfaces;
   trim_curve->surfaces( surfaces );
   if ( surfaces.size() == 0 )
   {
     CubitStatus status =
      AcisQueryEngine::instance()->delete_solid_model_entities(trim_curve);

       if (status == CUBIT_FAILURE)
       {
         PRINT_ERROR("In AcisModifyEngine::trim_curve\n"
           "       Could not delete Curve.\n"
           " The Model database is likely corrupted "
           "due to\nthis unsuccessful deletion.\n"  );
       }
   }

   return curve_ptr;
}


Curve*
AcisModifyEngine::create_arc_three( Point* vertex1, Point* vertex2,
                                      Point *vertex3, bool full )
{
  CubitVector vec1 = vertex1->coordinates();
  CubitVector vec2 = vertex2->coordinates();
  CubitVector vec3 = vertex3->coordinates();

  SPAposition pos1( vec1.x(), vec1.y(), vec1.z() );
  SPAposition pos2( vec2.x(), vec2.y(), vec2.z() );
  SPAposition pos3( vec3.x(), vec3.y(), vec3.z() );

  EDGE *EDGE_ptr;

  outcome result = api_curve_arc_3pt( pos1, pos2, pos3, full, EDGE_ptr );
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    PRINT_ERROR( "Unable to create arc from given vertices\n" );
    return NULL;
  }
 
  //reuse vertices
  if( full )
  {
    PointACIS *acis_point1 = CAST_TO( vertex1, PointACIS );
    EDGE_ptr->set_start(acis_point1->get_VERTEX_ptr() );
    EDGE_ptr->set_end(acis_point1->get_VERTEX_ptr() );
    acis_point1->get_VERTEX_ptr()->add_edge(EDGE_ptr);
  }
  else
  {

    PointACIS *acis_point1 = CAST_TO( vertex1, PointACIS );
    PointACIS *acis_point3 = CAST_TO( vertex3, PointACIS );

    EDGE_ptr->set_start(acis_point1->get_VERTEX_ptr() );
    acis_point1->get_VERTEX_ptr()->add_edge(EDGE_ptr);

    EDGE_ptr->set_end(acis_point3->get_VERTEX_ptr() );
    acis_point3->get_VERTEX_ptr()->add_edge(EDGE_ptr);
  }

  return AcisQueryEngine::instance()->populate_topology_bridges(EDGE_ptr);
}

Curve*
AcisModifyEngine::create_arc_three( Curve* ref_edge1, Curve* ref_edge2,
                                      Curve *ref_edge3, bool full )
{
  EDGE* EDGE_ptr1 = AcisQueryEngine::get_EDGE( ref_edge1 );
  EDGE* EDGE_ptr2 = AcisQueryEngine::get_EDGE( ref_edge2 );
  EDGE* EDGE_ptr3 = AcisQueryEngine::get_EDGE( ref_edge3 );
  if( !EDGE_ptr1 || !EDGE_ptr2 || !EDGE_ptr3)
    return NULL;

  SPAposition tang_pos( 0.0, 0.0, 0.0 );
  SPAunit_vector tang_vec( 0.0, 0.0, 0.0 );
  entity_with_ray eray1( EDGE_ptr1, tang_pos, tang_vec );
  entity_with_ray eray2( EDGE_ptr2, tang_pos, tang_vec );
  entity_with_ray eray3( EDGE_ptr3, tang_pos, tang_vec );

  EDGE *EDGE_ptr;

  outcome result = api_curve_arc_3curve( eray1, eray2, eray3, full, EDGE_ptr );
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    PRINT_ERROR( "Unable to create arc from given curves\n" );
    return NULL;
  }

  return AcisQueryEngine::instance()->populate_topology_bridges(EDGE_ptr);
}

Curve*
AcisModifyEngine::create_arc_center_edge( Point* vertex1, Point* vertex2,
                                          Point *vertex3, const CubitVector &normal,
                                          double radius, bool full )
{
  CubitVector vec1 = vertex1->coordinates(); // Center of arc
  CubitVector vec2 = vertex2->coordinates(); // Position on arc
  CubitVector vec3 = vertex3->coordinates(); // Position on arc

  // Re-adjust vec2 if radius was given
  if( radius != CUBIT_DBL_MAX )
  {
    CubitVector dir1( vec1, vec2 );
    vec1.next_point( dir1, radius, vec2 );
  }

  SPAposition pos1( vec1.x(), vec1.y(), vec1.z() );
  SPAposition pos2( vec2.x(), vec2.y(), vec2.z() );
  SPAposition pos3( vec3.x(), vec3.y(), vec3.z() );

  EDGE *EDGE_ptr;
  SPAunit_vector arc_normal;
  if( normal.length() != 0 )
  {
    arc_normal.set_x( normal.x() );
    arc_normal.set_y( normal.y() );
    arc_normal.set_z( normal.z() );
  }

  outcome result;
  if( full == CUBIT_FALSE )
  {
    if( normal.length() != 0.0 )
      result = api_curve_arc_center_edge( pos1, pos2, pos3, &arc_normal, EDGE_ptr );
    else
      result = api_curve_arc_center_edge( pos1, pos2, pos3, NULL, EDGE_ptr );
  }
  else
  {
    if( normal.length() == 0.0 )
    {
      CubitVector dir1( vec1, vec2 );
      CubitVector dir2( vec1, vec3 );
      CubitVector tmp_normal = dir1 * dir2; // Cross-product
      tmp_normal.normalize();
      arc_normal.set_x( tmp_normal.x() );
      arc_normal.set_y( tmp_normal.y() );
      arc_normal.set_z( tmp_normal.z() );
    }
    result = api_curve_arc_center_edge( pos1, pos2, pos2, &arc_normal, EDGE_ptr );
  }

  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    PRINT_ERROR( "Unable to create arc from given vertices\n" );
    return NULL;
  }

  //reuse vertices...if possible
  if( full )
  {
    PointACIS *acis_point2 = CAST_TO( vertex2, PointACIS );
    EDGE_ptr->set_start(acis_point2->get_VERTEX_ptr() );
    EDGE_ptr->set_end(acis_point2->get_VERTEX_ptr() );
    acis_point2->get_VERTEX_ptr()->add_edge(EDGE_ptr);
  }
  else
  {
    //reuse vertices only if it is coincident with 
    //new vertex on new curve

    //compare original start vs. new start vertex on EDGE 
    VERTEX *start_vert = EDGE_ptr->start();
    SPAposition temp_pos = start_vert->geometry()->coords();
    CubitVector tmp_pos1( temp_pos.x(), temp_pos.y(), temp_pos.z() ); 
    vec2 = vertex2->coordinates(); 
 
    if(tmp_pos1.distance_between( vec2 ) < GEOMETRY_RESABS )
    {
      PointACIS *acis_point2 = CAST_TO( vertex2, PointACIS );
      EDGE_ptr->set_start(acis_point2->get_VERTEX_ptr() );
      acis_point2->get_VERTEX_ptr()->add_edge(EDGE_ptr);
    }

    //compare original end vs. new end vertex on EDGE 
    VERTEX *end_vert = EDGE_ptr->end();
    temp_pos = end_vert->geometry()->coords();
    tmp_pos1.set( temp_pos.x(), temp_pos.y(), temp_pos.z() ); 
    vec3 = vertex3->coordinates(); 

    if(tmp_pos1.distance_between( vec3 ) < GEOMETRY_RESABS )
    {
      PointACIS *acis_point3 = CAST_TO( vertex3, PointACIS );
      EDGE_ptr->set_end(acis_point3->get_VERTEX_ptr() );
      acis_point3->get_VERTEX_ptr()->add_edge(EDGE_ptr);
    }
    else
    {
      //make sure this point isn't a free vertex
      if( !vertex3->topology_entity() )
        AcisQueryEngine::instance()->delete_solid_model_entities( vertex3 );
    }
  }


  return AcisQueryEngine::instance()->populate_topology_bridges(EDGE_ptr);
}

CubitStatus
AcisModifyEngine::create_curve_combine( DLIList<Curve*>& curve_list,
                                    Curve *&new_curve_ptr )
{
  CubitStatus result;

  result = AcisEdgeTool::instance()->create_curve_combine(curve_list,
                                                          new_curve_ptr);

  if (result)
  {
    RefEdge* new_ref_edge_ptr = GeometryQueryTool::instance()->make_free_RefEdge( new_curve_ptr );

    if (new_ref_edge_ptr)
      PRINT_INFO( "Created new \"combined\" curve %d\n", new_ref_edge_ptr->id() );
    else
      result = CUBIT_FAILURE;
  }

  return result;

}

CubitStatus
AcisModifyEngine::get_copied_FACES_of_body( DLIList<SurfaceACIS*>& ref_face_list,
                                            DLIList<FACE*>& FACE_list,
                                            DLIList<SurfaceACIS*>& removed_ref_faces,
                                            BODY*& copied_BODY_ptr ) const
{
  int i, num_faces;

  // Note: we will be pulling surfaces out of incoming ref_face_list.
  SurfaceACIS* ref_face_ptr;
  BodySM *body_ptr;

  BODY *BODY_ptr;
  BODY *BODY_ptr2;

  FACE *FACE_ptr;

  copied_BODY_ptr = NULL;

  outcome result;

  DLIList<SurfaceACIS*> body_face_list;

  ref_face_list.reset();
  ref_face_ptr = ref_face_list.remove();
  removed_ref_faces.append( ref_face_ptr );

  FACE_ptr = ref_face_ptr->get_FACE_ptr();
  if (!FACE_ptr)
  {
    PRINT_ERROR( "Unable to get Acis FACE from Surface\n"  );
     return CUBIT_FAILURE;
  }

  body_face_list.append_unique( ref_face_ptr );

  body_ptr = AcisQueryEngine::instance()->get_body_sm_of_ENTITY( FACE_ptr );
  if( body_ptr == NULL )
  {
    PRINT_ERROR( "Unable to get volume from FACE\n"  );
    return CUBIT_FAILURE;
  }

  BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY( FACE_ptr );
  if( BODY_ptr == NULL )
  {
    PRINT_ERROR( "Unable to get volume from FACE\n" );
    return CUBIT_FAILURE;
  }

  // Add all remaining faces from this same body to face list
  num_faces = ref_face_list.size();
  for( i=0; i<num_faces; i++ )
  {
    ref_face_ptr = ref_face_list.get();

    BODY_ptr2 = AcisQueryEngine::instance()->get_BODY_of_entity( ref_face_ptr );

    if( BODY_ptr == BODY_ptr2 )
    {
      body_face_list.append_unique( ref_face_ptr );
      ref_face_list.remove();
      removed_ref_faces.append( ref_face_ptr );
    }
    else
      ref_face_list.step();
  }

  // Now we have all the faces from this particular BODY.

  // Copy the body before working with it.  Here we keep the
  // attributes on the BODY so we can get back to the TE's
  // in CUBIT from the copy (to find which FACE's to correlate
  // to the RefFace's).
  copied_BODY_ptr = copy_BODY(BODY_ptr, CUBIT_FALSE);

  // Loop through the copied BODY's FACES to find which one's
  // match to the original BODY's RefFaces.
  ENTITY_LIST ENTITIES;
  result = api_get_faces( copied_BODY_ptr, ENTITIES );
  if( !result.ok() )
  {
    PRINT_ERROR( "Unable to get ACIS FACES from volume\n" );
    api_delent( copied_BODY_ptr );
    copied_BODY_ptr = NULL;
    return CUBIT_FAILURE;
  }

  ENTITIES.init();
  ENTITY* ENTITY_ptr;
  int index = -1;
  DLIList<int> index_list;
  DLIList<FACE*> temp_FACE_list;
  while( (ENTITY_ptr = ENTITIES.next() ) != NULL )
  {
    AcisBridge* ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(ENTITY_ptr);
    ref_face_ptr = CAST_TO( ab_ptr, SurfaceACIS );
    if( ref_face_ptr == NULL )
    {
      PRINT_ERROR( "Unable to find RefFace from ACIS FACE!\n" );
      api_delent( copied_BODY_ptr );
      copied_BODY_ptr = NULL;
      return CUBIT_FAILURE;
    }

    index = body_face_list.where_is_item( ref_face_ptr );
    if( index != -1 )
    {
      index_list.append( index );
      temp_FACE_list.append( (FACE *)ENTITY_ptr );
    }
    if( temp_FACE_list.size() == body_face_list.size() )
      break;
  }
  ENTITIES.clear();

  if( temp_FACE_list.size() != body_face_list.size() )
  {
    PRINT_ERROR( "Internal error correlating FACE lists\n" );
    api_delent( copied_BODY_ptr );
    copied_BODY_ptr = NULL;
    return CUBIT_FAILURE;
  }

  // Re-sort the list into the proper order
  temp_FACE_list.reset();
  index_list.reset();
  for( i=0; i<index_list.size(); i++ )
  {
    // Find the FACE with the corresponding index
    index = index_list.where_is_item( i );
    temp_FACE_list.reset();
    temp_FACE_list.step( index );

    FACE_list.append( temp_FACE_list.get() );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisModifyEngine::get_copied_FACES_of_body( DLIList<SurfaceACIS*>& ref_face_list,
                                              DLIList<FACE*>& FACE_list,
                                              BODY*& copied_BODY_ptr ) const
{
  DLIList<SurfaceACIS*> removed_faces;

  return get_copied_FACES_of_body( ref_face_list, FACE_list, removed_faces,
                                   copied_BODY_ptr );
}

CubitStatus
AcisModifyEngine::get_copied_EDGES_of_body( DLIList<CurveACIS*>& ref_edge_list,
                                            DLIList<EDGE*>& EDGE_list,
                                            DLIList<CurveACIS*>& removed_ref_edges,
                                            BODY*& copied_BODY_ptr ) const
{
  int i, num_edges;

  // Note: we will be pulling surfaces out of incoming ref_edge_list.
  CurveACIS* ref_edge_ptr;
  BodySM *body_ptr;

  BODY *BODY_ptr;
  BODY *BODY_ptr2;

  EDGE *EDGE_ptr;

  copied_BODY_ptr = NULL;

  outcome result;

  DLIList<CurveACIS*> body_edge_list;

  ref_edge_ptr = ref_edge_list.remove();
  removed_ref_edges.append( ref_edge_ptr );

  EDGE_ptr = ref_edge_ptr->get_EDGE_ptr();

  body_edge_list.append_unique( ref_edge_ptr );

  body_ptr = AcisQueryEngine::instance()->get_body_sm_of_ENTITY( EDGE_ptr );
  if( body_ptr == NULL )
  {
    PRINT_ERROR( "Unable to get volume from EDGE\n" );
    return CUBIT_FAILURE;
  }

  BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY( EDGE_ptr );
  if( BODY_ptr == NULL )
  {
    PRINT_ERROR( "Unable to get volume from EDGE\n" );
    return CUBIT_FAILURE;
  }

  // Add all remaining edges from this same body to edge list
  num_edges = ref_edge_list.size();
  for( i=0; i<num_edges; i++ )
  {
    ref_edge_ptr = ref_edge_list.get();

    BODY_ptr2 = AcisQueryEngine::instance()->get_BODY_of_entity( ref_edge_ptr );

    if( BODY_ptr == BODY_ptr2 )
    {
      body_edge_list.append_unique( ref_edge_ptr );
      ref_edge_list.remove();
      removed_ref_edges.append( ref_edge_ptr );
    }
    else
      ref_edge_list.step();
  }

  // Now we have all the edges from this particular BODY.

  // Copy the body before working with it.  Here we keep the
  // attributes on the BODY so we can get back to the TE's
  // in CUBIT from the copy (to find which EDGE's to correlate
  // to the RefEdge's).
  copied_BODY_ptr = copy_BODY(BODY_ptr, CUBIT_FALSE);

  // Loop through the copied BODY's EDGES to find which one's
  // match to the original BODY's RefEdges.
  ENTITY_LIST ENTITIES;
  result = api_get_edges(copied_BODY_ptr, ENTITIES);
  if( !result.ok() )
  {
    PRINT_ERROR( "Unable to get ACIS EDGES from volume\n" );
    api_delent( copied_BODY_ptr );
    copied_BODY_ptr = NULL;
    return CUBIT_FAILURE;
  }

  ENTITIES.init();
  ENTITY* ENTITY_ptr;
  int index = -1;
  DLIList<int> index_list;
  DLIList<EDGE*> temp_EDGE_list;
  while( (ENTITY_ptr = ENTITIES.next() ) != NULL )
  {
    AcisBridge* ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(ENTITY_ptr);
    ref_edge_ptr = CAST_TO( ab_ptr, CurveACIS );
    if( ref_edge_ptr == NULL )
    {
      PRINT_ERROR( "Unable to find Cubit curve from ACIS EDGE!\n" );
      api_delent( copied_BODY_ptr );
      copied_BODY_ptr = NULL;
      return CUBIT_FAILURE;
    }
    index = body_edge_list.where_is_item( ref_edge_ptr );
    if( index != -1 )
    {
      index_list.append( index );
      temp_EDGE_list.append( (EDGE *)ENTITY_ptr );
    }
    if( temp_EDGE_list.size() == body_edge_list.size() )
      break;
  }
  ENTITIES.clear();

  if( temp_EDGE_list.size() != body_edge_list.size() )
  {
    PRINT_ERROR( "Internal error correlating EDGE lists\n" );
    api_delent( copied_BODY_ptr );
    copied_BODY_ptr = NULL;
    return CUBIT_FAILURE;
  }

  // Re-sort the list into the proper order
  temp_EDGE_list.reset();
  index_list.reset();
  for( i=0; i<index_list.size(); i++ )
  {
    // Find the EDGE with the corresponding index
    index = index_list.where_is_item( i );
    temp_EDGE_list.reset();
    temp_EDGE_list.step( index );

    EDGE_list.append( temp_EDGE_list.get() );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisModifyEngine::get_copied_EDGES_of_body( DLIList<CurveACIS*>& ref_edge_list,
                                            DLIList<EDGE*>& EDGE_list,
                                            BODY*& copied_BODY_ptr ) const
{
  DLIList<CurveACIS*> removed_edges;

  return get_copied_EDGES_of_body( ref_edge_list, EDGE_list, removed_edges,
                                   copied_BODY_ptr );
}

//-------------------------------------------------------------------------
// Purpose       :  To get a pointer to a copy of the body along with a list of VERTEXs
//
// Special Notes :
//
// Creator       : Lingyun Pan (CAT)
//
// Creation Date : 07/19/01
//-------------------------------------------------------------------------

CubitStatus
AcisModifyEngine::get_copied_VERTICES_of_body(
                                     DLIList<PointACIS*>& ref_vertex_list,
                                     DLIList<VERTEX*>& VERTEX_list,
                                     DLIList<PointACIS*>& removed_ref_vertices,
                                     BODY*& copied_BODY_ptr ) const
{
  int i, num_vertices;

  // Note: we will be pulling vertices out of incoming ref_vertex_list.
  PointACIS* ref_vertex_ptr;
  BodySM *body_ptr;

  BODY *BODY_ptr;
  BODY *BODY_ptr2;

  VERTEX *VERTEX_ptr;

  copied_BODY_ptr = NULL;

  outcome result;

  DLIList<PointACIS*> body_vertex_list;

  ref_vertex_ptr =ref_vertex_list .remove();
  removed_ref_vertices.append( ref_vertex_ptr);

  VERTEX_ptr = ref_vertex_ptr->get_VERTEX_ptr();

  body_vertex_list.append_unique( ref_vertex_ptr );

  body_ptr = AcisQueryEngine::instance()->get_body_sm_of_ENTITY(VERTEX_ptr );
  if( body_ptr == NULL )
  {
    PRINT_ERROR( "Unable to get volume from VERTEX\n" );
    return CUBIT_FAILURE;
  }

  BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY( VERTEX_ptr );
  if( BODY_ptr == NULL )
  {
    PRINT_ERROR( "Unable to get volume from VERTEX\n" );
    return CUBIT_FAILURE;
  }

  // Add all remaining vertex from this same body to vertex list
  num_vertices = ref_vertex_list.size();
  for( i=0; i<num_vertices; i++ )
  {
    ref_vertex_ptr = ref_vertex_list.get();

    BODY_ptr2 = AcisQueryEngine::instance()->get_BODY_of_entity( ref_vertex_ptr );

    if( BODY_ptr == BODY_ptr2 )
    {
      body_vertex_list.append_unique( ref_vertex_ptr );
      ref_vertex_list.remove();
      removed_ref_vertices.append( ref_vertex_ptr );
    }
    else
      ref_vertex_list.step();
  }

  // Now we have all the vertices from this particular BODY.

  // Copy the body before working with it.  Here we keep the
  // attributes on the BODY so we can get back to the TE's
  // in CUBIT from the copy (to find which VERTICE's to correlate
  // to the RefVertice's).
  copied_BODY_ptr = copy_BODY(BODY_ptr, CUBIT_FALSE);

  // Loop through the copied BODY's VERTICEs to find which one's
  // match to the original BODY's RefVertices.
  ENTITY_LIST ENTITIES;
  result = api_get_vertices( copied_BODY_ptr, ENTITIES);
  if( !result.ok() )
  {
    PRINT_ERROR( "Unable to get ACIS VERTICES from a volume\n" );
    api_delent( copied_BODY_ptr );
    copied_BODY_ptr = NULL;
    return CUBIT_FAILURE;
  }

  ENTITIES.init();
  ENTITY* ENTITY_ptr;
  int index = -1;
  DLIList<int> index_list;
  DLIList<VERTEX*> temp_VERTEX_list;
  while( (ENTITY_ptr = ENTITIES.next() ) != NULL )
  {
    AcisBridge* ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(ENTITY_ptr);
    ref_vertex_ptr = CAST_TO(ab_ptr, PointACIS);
    if( ref_vertex_ptr == NULL )
    {
      PRINT_ERROR( "Unable to find RefVertex from ACIS VERTEX!\n" );
      VERTEX_list.clean_out();
      api_delent( copied_BODY_ptr );
      copied_BODY_ptr = NULL;
      VERTEX_list.clean_out();
      return CUBIT_FAILURE;
    }

    index = body_vertex_list.where_is_item( ref_vertex_ptr );
    if( index != -1 )
    {
      index_list.append( index );
      temp_VERTEX_list.append( (VERTEX *)ENTITY_ptr );
    }
    if( temp_VERTEX_list.size() == body_vertex_list.size() )
      break;
  }
  ENTITIES.clear();

  if( temp_VERTEX_list.size() != body_vertex_list.size() )
  {
    PRINT_ERROR( "Internal error correlating VERTEX lists\n" );
    api_delent( copied_BODY_ptr );
    copied_BODY_ptr = NULL;
    return CUBIT_FAILURE;
  }

  // Re-sort the list into the proper order
  temp_VERTEX_list.reset();
  index_list.reset();
  for( i=0; i<index_list.size(); i++ )
  {
    // Find the VERTEX with the corresponding index
    index = index_list.where_is_item( i );
    temp_VERTEX_list.reset();
    temp_VERTEX_list.step( index );

    VERTEX_list.append( temp_VERTEX_list.get() );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisModifyEngine::get_copied_VERTICES_of_body(
                                        DLIList<PointACIS*>& ref_vertex_list,
                                        DLIList<VERTEX*>& VERTEX_list,
                                        BODY*& copied_BODY_ptr ) const
{
  DLIList<PointACIS*> removed_vertices;

  return get_copied_VERTICES_of_body( ref_vertex_list, VERTEX_list,
                                      removed_vertices, copied_BODY_ptr );
}


GeometryQueryEngine *AcisModifyEngine::get_gqe()
{
  return AcisQueryEngine::instance();
}

CubitBoolean AcisModifyEngine::is_modify_engine(const TopologyBridge *tb_ptr) const
{

  if (CAST_TO(const_cast<TopologyBridge*>(tb_ptr), AcisBridge)) return CUBIT_TRUE;
  else return CUBIT_FALSE;
}

CubitStatus
AcisModifyEngine::get_offset_intersections( Curve* ref_edge1,
                                            Curve* ref_edge2,
                                            DLIList<CubitVector*>& intersection_list,
                                            double offset_distance,
                                            bool ext_first )
{
   int i;
   //int done = 0;
   bool bounded = !ext_first;

   // Get the ACIS EDGES
   EDGE* EDGE_ptr1 = AcisQueryEngine::get_EDGE(ref_edge1);
   EDGE* EDGE_ptr2 = AcisQueryEngine::get_EDGE(ref_edge2);
   if (!EDGE_ptr1 || !EDGE_ptr2)
    return CUBIT_FAILURE;

   // Determine which curves are straight.  Use this because some curves,
   // esp. from IGES, can be "straight" splines.
   bool straight1, straight2;
   straight1 = AcisQueryEngine::instance()->is_curve_app_straight( ref_edge1 );
   straight2 = AcisQueryEngine::instance()->is_curve_app_straight( ref_edge2 );

   // If first is straight we need to find the plane to offset it in.
   CubitVector pln_norm;
   EDGE* offset_EDGE_ptr = NULL;
   if( straight1 )
   {
      // Need to find the plane normal to do offsets in

      // Get point normal form
      CubitVector orig1, dir1;
      if( ref_edge1->get_point_direction( orig1, dir1 ) == CUBIT_FAILURE )
      {
        // This must be a straight spline.  Get it from endpoints.
        ref_edge1->position_from_fraction( 0.0, orig1 );
        CubitVector end_coords;
        ref_edge1->position_from_fraction( 1.0, end_coords );
        dir1 = end_coords - orig1;
        dir1.normalize();
      }
      if( straight2 )
      {
         // Get point normal form
         CubitVector orig2, dir2;
         if( ref_edge2->get_point_direction( orig2, dir2 ) == CUBIT_FAILURE )
         {
           // This must be a straight spline.  Get it from endpoints.
           ref_edge2->position_from_fraction( 0.0, orig2 );
           CubitVector end_coords;
           ref_edge2->position_from_fraction( 1.0, end_coords );
           dir2 = end_coords - orig2;
           dir2.normalize();
         }

         // Find SPAvector from orig1 to orig2
         CubitVector vec1 = orig2 - orig1;
         vec1.normalize();

         // Cross dir1 with vec1 to get plane normal.
         pln_norm = dir1 * vec1;

         // Make sure this isn't a zero SPAvector.
         if( pln_norm.length() == 0.0 )
         {
           // Cross dir1 with dir2 to get plane normal
           pln_norm = dir1 * dir2;
         }

         if( pln_norm.length() == 0.0 )
         {
            PRINT_ERROR( "Unable to calculate plane that the curves lie in\n" );
            return CUBIT_FAILURE;
         }
      }
      else // EDGE2 is not straight
      {
         // Proper normal should be normal to curve, but check to see if EDGE 1 is on
         // this plane.
         CubitVector orig2;
         AcisQueryEngine::instance()->get_EDGE_normal( EDGE_ptr2, pln_norm );

         double aln_orig[3], aln_vec[3], apln_orig[3], apln_norm[3];

         AnalyticGeometryTool *agt = AnalyticGeometryTool::instance();

         agt->copy_pnt( orig1, aln_orig ); agt->copy_pnt( dir1, aln_vec );
         ref_edge2->position_from_fraction( 0.0, orig2 ); //
         ref_edge2->position_from_u(ref_edge2->start_param(),orig2);
         agt->copy_pnt( orig2, apln_orig );
         agt->copy_pnt( pln_norm, apln_norm );

         if( !agt->is_ln_on_pln( aln_orig,aln_vec, apln_orig, apln_norm ) )
         {
            PRINT_ERROR( "Curve does not lie in the plane defined by Curve\n" );
            return CUBIT_FAILURE;
         }
      }

      // Get offset direction cross of plane normal and edge 1 direction
      CubitVector offset_dir = pln_norm*dir1;

      // Do first side
      offset_EDGE_ptr = copy_and_move_EDGE( EDGE_ptr1, offset_dir, offset_distance );
      if( offset_EDGE_ptr == NULL )
      {
         PRINT_ERROR( "unable to offset Curve for intersection calculation\n" );
         return CUBIT_FAILURE;
      }

      if( AcisQueryEngine::instance()->get_intersections( offset_EDGE_ptr, EDGE_ptr2,
                                              intersection_list, bounded ) ==
          CUBIT_FAILURE )
      {
         PRINT_ERROR( "unable to calculate intersections of offset Curves\n" );
         api_delent( offset_EDGE_ptr );
         return CUBIT_FAILURE;
      }

      api_delent( offset_EDGE_ptr );

      // Do second side
      if( offset_distance != 0.0 )
      {
         offset_dir = -offset_dir;
         offset_EDGE_ptr = copy_and_move_EDGE( EDGE_ptr1, offset_dir, offset_distance );
         if( offset_EDGE_ptr == NULL )
         {
            PRINT_ERROR( "unable to offset Curve for intersection calculation\n" );
            return CUBIT_FAILURE;
         }

         if( AcisQueryEngine::instance()->get_intersections( offset_EDGE_ptr, EDGE_ptr2, intersection_list, bounded ) == CUBIT_FAILURE )
         {
            api_delent( offset_EDGE_ptr );
            PRINT_ERROR( "unable to calculate intersections of offset Curves\n" );
            return CUBIT_FAILURE;
         }

         api_delent( offset_EDGE_ptr );
      }

  }
  else if( offset_distance != 0.0 )
  {
     // If there is an offset, offset the first curve to each "side" of the original
     // and do the intersections - i.e., try each side.  Note we already handled the
     // special case if the first curve is straight.

     // Do first side
     DLIList<EDGE*> EDGE_list;
     DLIList<EDGE*> new_EDGE_list;

     EDGE_list.append( EDGE_ptr1 );
     if( offset_EDGES( EDGE_list, offset_distance, 1, new_EDGE_list ) == CUBIT_FAILURE )
     {
        PRINT_ERROR( "Unable to offset the curves\n" );
        return CUBIT_FAILURE;
     }

     // Note: the new_EDGE_list can have more curves in it than EDGE_list
     for( i=0; i<new_EDGE_list.size(); i++ )
     {
        offset_EDGE_ptr = new_EDGE_list.get_and_step();
        if( AcisQueryEngine::instance()->get_intersections( offset_EDGE_ptr,
                  EDGE_ptr2, intersection_list, bounded ) == CUBIT_FAILURE )
        {
           api_delent( offset_EDGE_ptr );
           PRINT_ERROR( "unable to calculate intersections of offset Curves\n" );
           return CUBIT_FAILURE;
        }
        api_delent( offset_EDGE_ptr );
     }

     new_EDGE_list.clean_out();

     // Do other side
     if( offset_EDGES( EDGE_list, -offset_distance, 1, new_EDGE_list ) == CUBIT_FAILURE )
     {
        PRINT_ERROR( "Unable to offset the curves\n" );
        return CUBIT_FAILURE;
     }

     for( i=0; i<new_EDGE_list.size(); i++ )
     {
        offset_EDGE_ptr = new_EDGE_list.get_and_step();
        if( AcisQueryEngine::instance()->get_intersections( offset_EDGE_ptr, EDGE_ptr2, intersection_list, bounded ) == CUBIT_FAILURE )
        {
           api_delent( offset_EDGE_ptr );
           PRINT_ERROR( "unable to calculate intersections of offset Curves\n" );
           return CUBIT_FAILURE;
        }
        api_delent( offset_EDGE_ptr );
     }

     new_EDGE_list.clean_out();

  }
  else
  {
     // Just get the intersections
     if( AcisQueryEngine::instance()->get_intersections( EDGE_ptr1, EDGE_ptr2, intersection_list, bounded ) == CUBIT_FAILURE )
     {
        PRINT_ERROR( "unable to calculate intersections of Curves\n" );
        return CUBIT_FAILURE;
     }
  }

  // Cull out any points that are not on bounds of RefEdge 2
  if( bounded == CUBIT_FALSE )
  {
     CubitVector *vec_ptr;
     for( i=0; i<intersection_list.size(); i++ )
     {
        vec_ptr = intersection_list.get();

        if( ref_edge2->point_containment( *vec_ptr ) == CUBIT_PNT_OFF )
        {
           delete vec_ptr;
           intersection_list.remove();
        }
        else
           intersection_list.step();
     }
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisModifyEngine::get_offset_intersections( Curve* ref_edge_ptr,
                                            Surface* ref_face_ptr,
                                            DLIList<CubitVector*> &intersection_list,
                                            double offset_distance,
                                            bool ext_surf )
{
   outcome result;

   EDGE* EDGE_ptr = AcisQueryEngine::get_EDGE(ref_edge_ptr);
   FACE* FACE_ptr = AcisQueryEngine::get_FACE(ref_face_ptr);
   if (!EDGE_ptr || !FACE_ptr)
    return CUBIT_FAILURE;

   // If there is an offset, offset the surface to each "side" of the original
   // and do the intersections - i.e., try each side.
   if( offset_distance != 0.0 )
   {
      FACE* offset_FACE_ptr = NULL;

      result = api_offset_face ( FACE_ptr, offset_distance, offset_FACE_ptr );

      if( offset_FACE_ptr == NULL )
      {
         AcisQueryEngine::instance()->ACIS_API_error(result);
         PRINT_ERROR( "Unable to offset surface\n" );
         return CUBIT_FAILURE;
      }

      // This creates a FACE from the given FACE, only it can be extended out.
      FACE* ext_offset_FACE_ptr = make_FACE( offset_FACE_ptr, ext_surf );
      BODY *ext_offset_BODY_ptr =
        AcisQueryEngine::instance()->get_BODY_of_ENTITY( ext_offset_FACE_ptr );

      if( AcisQueryEngine::instance()->get_intersections( EDGE_ptr, ext_offset_FACE_ptr, intersection_list ) == CUBIT_FAILURE )
      {
         PRINT_ERROR("Problems finding intersections of curve and surface.\n");
         api_delent( offset_FACE_ptr );
         api_delent( ext_offset_BODY_ptr );
         return CUBIT_FAILURE;
      }

      api_delent( offset_FACE_ptr ); offset_FACE_ptr = NULL;
      api_delent( ext_offset_BODY_ptr ); ext_offset_BODY_ptr = NULL;
      ext_offset_FACE_ptr = NULL;

      // Try other side
      result = api_offset_face ( FACE_ptr, -offset_distance, offset_FACE_ptr );

      if( offset_FACE_ptr == NULL )
      {
         AcisQueryEngine::instance()->ACIS_API_error(result);
         PRINT_ERROR( "Unable to offset surface\n" );
         return CUBIT_FAILURE;
      }

      // This creates a FACE from the given FACE, only it is extended out.
      ext_offset_FACE_ptr = make_FACE( offset_FACE_ptr, CUBIT_TRUE );
      ext_offset_BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY( ext_offset_FACE_ptr );

      if( AcisQueryEngine::instance()->get_intersections( EDGE_ptr, ext_offset_FACE_ptr, intersection_list ) == CUBIT_FAILURE )
      {
         PRINT_ERROR("Problems finding intersections of curve and surface.\n" );
         api_delent( offset_FACE_ptr );
         api_delent( ext_offset_BODY_ptr );
         return CUBIT_FAILURE;
      }

      api_delent( offset_FACE_ptr );
      api_delent( ext_offset_BODY_ptr );
   }
   else
   {
      FACE* ext_FACE_ptr = make_FACE( FACE_ptr, ext_surf );
      BODY *ext_BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY( ext_FACE_ptr );
      if( AcisQueryEngine::instance()->get_intersections( EDGE_ptr, ext_FACE_ptr, intersection_list ) == CUBIT_FAILURE )
      {
         PRINT_ERROR("Problems finding intersections of curve and surface.\n" );
         api_delent( ext_BODY_ptr );
         return CUBIT_FAILURE;
      }
      api_delent( ext_BODY_ptr );
   }

   return CUBIT_SUCCESS;
}

EDGE *
AcisModifyEngine::copy_and_move_EDGE( EDGE *EDGE_ptr,
                                      const CubitVector &input_dir,
                                      double dist )
{
   // Setup translation transform
   SPAtransf tr;
   CubitVector dir(input_dir);
   dir.normalize(); // For now, allow this SPAvector to be modified
   CubitVector delta = dist*dir;
   SPAvector vec(delta.x(), delta.y(), delta.z() );
   tr = translate_transf( vec );

   EDGE* new_EDGE_ptr;
   outcome result = api_trans_edge( EDGE_ptr, tr, new_EDGE_ptr );
   if (!result.ok())
   {
     AcisQueryEngine::instance()->ACIS_API_error( result );
     return NULL;
   }

   // Remove the CUBIT owner attributes from the new edge
   ATTRIB_CUBIT_OWNER::remove_cubit_owner((ENTITY*)new_EDGE_ptr, CUBIT_TRUE);

   return new_EDGE_ptr;
}

CubitStatus
AcisModifyEngine::offset_EDGES( DLIList<EDGE*> &EDGE_list,
                                double offset_distance,
                                int gap_type,
                                DLIList<EDGE*> &new_EDGE_list )
{
   int i;
   EDGE *EDGE_ptr = NULL;
   //ENTITY *ENTITY_ptr = NULL;
   EDGE *tmp_EDGE_ptr = NULL;

   // Need to create a wire from the EDGEs

   // Copy them first
   int edge_count = EDGE_list.size();
   EDGE** EDGEs = new EDGE*[edge_count];

   EDGE_list.reset();
   for( i=0; i<EDGE_list.size(); i++ )
   {
      EDGE_ptr = EDGE_list.get_and_step();

      api_edge( EDGE_ptr, tmp_EDGE_ptr );

      EDGEs[i] = tmp_EDGE_ptr;

      // Remove the owner attribute from the copied edge & children
      ATTRIB_CUBIT_OWNER::remove_cubit_owner(EDGEs[i], CUBIT_TRUE);
   }

   // The EDGEs are put into the wire so don't delete them.
   BODY* wire_BODY = NULL;
   outcome result = api_make_ewire( edge_count, EDGEs, wire_BODY );

   if( !result.ok() || wire_BODY==NULL )
   {
     AcisQueryEngine::instance()->ACIS_API_error(result);
     PRINT_ERROR( "unable to make ACIS wire body from curves\n" );
     if( wire_BODY )
       api_delent( (ENTITY*)wire_BODY );
     else
     {
       for( i=0; i<edge_count; i++)
         api_delent( (ENTITY*)EDGEs[i] );
     }
     return CUBIT_FAILURE;
   }

   //The gap type is as follows; 0 = rounded like arcs, 1 = extended like lines, and 2 =
   //natural like curve extensions. Lines are extend with lines and circles are extended with
   //circles.

   law* offset_law = new constant_law( offset_distance );
   law *twist_law=NULL;
   BODY* offset_WIRE;
   logical trim_flg = TRUE;

   // need to check if wire is planar and calculate the
   // normal to the plane
   WIRE* this_wire = wire_BODY->wire() ?
      wire_BODY->wire() : wire_BODY->lump()->shell()->wire();

   SPAposition centroid;
   SPAunit_vector normal;
   if(!is_planar_wire(this_wire, centroid, normal))
   {
      PRINT_ERROR("selected curves do not form a planar chain\n");
      api_delent( (ENTITY*)wire_BODY );
      return CUBIT_FAILURE;
   }

   // Offset the wire to get the new curves
   result = api_offset_planar_wire( wire_BODY, offset_law, twist_law,
                                    normal, offset_WIRE, gap_type, trim_flg );

   if( !result.ok() )
   {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "unable to offset wire volume\n" );
      api_delent( (ENTITY*)wire_BODY );
      return CUBIT_FAILURE;
   }

   // Copy out the EDGES from the WIRE, then delete the WIRE
   DLIList<EDGE*> wire_edge_list;
   AcisQueryEngine::instance()->get_EDGEs( (ENTITY *)offset_WIRE, wire_edge_list );

   // Sometimes the function works but we don't get any edges
   if( !wire_edge_list.size() )
   {
     api_delent( (ENTITY*)wire_BODY );
     return CUBIT_FAILURE;
   }

   wire_edge_list.reset();
   EDGE *new_EDGE_ptr;
   for( i=0; i<wire_edge_list.size(); i++ )
   {
      EDGE_ptr = wire_edge_list.get_and_step();

      result = api_edge( EDGE_ptr, new_EDGE_ptr );
      if( !result.ok() )
      {
         AcisQueryEngine::instance()->ACIS_API_error(result);
         PRINT_INFO( "Unable to copy curve from new WIRE\n" );
         api_delent( (ENTITY*)wire_BODY );
         return CUBIT_FAILURE;
      }

      new_EDGE_list.append( new_EDGE_ptr );
   }

   // Delete the WIRE
   api_delent( (ENTITY*)wire_BODY );

   return CUBIT_SUCCESS;
}

Curve* AcisModifyEngine::find_curve_by_end_coord( const CubitVector& coords,
                                                  DLIList<Curve*>& curves,
                                                  bool& start_flag ) const
{
  Curve* curve = 0;
  CubitVector curve_pt;
  const double tol = AcisQueryEngine::instance()->get_sme_resabs_tolerance();
  const double tolsqr = tol*tol;

  DLIList<Curve*> tmp_list(curves);
  tmp_list.reset();
  for ( int i = curves.size(); i--; )
  {
    curve = tmp_list.get_and_step();

    curve->position_from_fraction( 0.0, curve_pt );
    if ( (curve_pt - coords).length_squared() < tolsqr )
    {
      start_flag = true;
      return curve;
    }

    curve->position_from_fraction( 1.0, curve_pt );
    if ( (curve_pt - coords).length_squared() < tolsqr )
    {
      start_flag = false;
      return curve;
    }
  }

  return 0;
}

bool AcisModifyEngine::bridge_deactivated(AcisBridge* bridge) const
{
  return deactivatedSet.count(bridge) > 0;
}

CubitStatus AcisModifyEngine::deactivate_bridge( AcisBridge* bridge ) const
{
  AcisModifyEngine* nonconst = AcisModifyEngine::instance();
  bool unique = nonconst->deactivatedSet.insert(bridge).second;
  return (CubitStatus)unique;
}

CubitStatus AcisModifyEngine::reactivate_bridge( AcisBridge* bridge ) const
{
  AcisModifyEngine* nonconst = AcisModifyEngine::instance();
  int count = nonconst->deactivatedSet.erase(bridge);
  return count > 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

void AcisModifyEngine::cleanout_deactivated_geometry() const
{
  AcisModifyEngine* nonconst = AcisModifyEngine::instance();
//  std::set<AcisBridge*>::iterator itor = nonconst->deactivatedSet.begin(),
//                                  end  = nonconst->deactivatedSet.end();
//  for ( ; itor != end; ++itor )
//    delete *itor;

  nonconst->deactivatedSet.clear();
}

CubitStatus AcisModifyEngine::scale( BodySM *&body, const CubitVector& f )
{
#if CUBIT_ACIS_VERSION < 1100
  api_initialize_operators();
#elif defined(WIN32) || defined(MACOSX)
  api_initialize_warp();
#endif
   BODY *new_BODY = NULL;

  if( f.x() != f.y() ||
      f.x() != f.z() ||
      f.y() != f.z() )
  {
    BODY *tmp_BODY = AcisQueryEngine::get_BODY( body );
    new_BODY = AcisModifyEngine::instance()->copy_BODY(tmp_BODY,
                               GeometryModifyTool::instance()->get_new_ids());
  }

  CubitStatus result = AcisQueryEngine::instance()->transform( new_BODY, scale_transf( f.x(), f.y(), f.z() ) );

#if CUBIT_ACIS_VERSION < 1100
  api_terminate_operators();
#elif defined(WIN32) || defined(MACOSX)
  api_terminate_warp();
#endif

  if( new_BODY )
  {
    BODY *old_BODY = AcisQueryEngine::get_BODY( body );
    BodySM *new_body = get_new_Body( body, old_BODY, new_BODY, false);
    body = new_body;
  }

  return result;
}

CubitStatus AcisModifyEngine::tolerant_imprint( DLIList<BodySM*> &bodies_in,
                                                 DLIList<BodySM*> &new_bodies,
                                                  DLIList<TopologyBridge*> *new_tbs,
                                                  DLIList<TopologyBridge*> *att_tbs ) const
{
  //make sure all bodies are from the same modify engine
  DLIList<BodySM*> new_body_list;

  ProgressTool *progress_tool = NULL;
  if( bodies_in.size() > 2 )
  {
     progress_tool = AppUtil::instance()->progress_tool();
     progress_tool->start(0, 100, "Tolerant Imprinting" );
  }

  CubitStatus status = imprint_overlapping_curves( bodies_in, new_body_list, progress_tool );

  status = imprint_overlapping_surfaces( bodies_in, new_body_list, progress_tool );

  status = imprint_overlapping_curves( bodies_in, new_body_list, progress_tool, new_tbs, att_tbs );

  if( progress_tool )
    progress_tool->end();

  return CUBIT_SUCCESS;
}

CubitStatus AcisModifyEngine::imprint_overlapping_curves( DLIList<BodySM*> &body_sms,
                                                          DLIList<BodySM*> &new_body_sms,
                                                          ProgressTool *progress_tool,
                                                          DLIList<TopologyBridge*> *new_tbs,
                                                          DLIList<TopologyBridge*> *att_tbs ) const
{
  //find all mergeable curves between the bodies
  CubitStatus status;
#ifdef BOYD16
  DLIList<RefEntity*> entities;
#endif
  DLIList< DLIList<Curve*>* > lists_of_mergeable_curves;

  status = MergeTool::instance()->find_only_mergeable_curves( body_sms, lists_of_mergeable_curves );
  
  int i,j;
  //increment 4%
  if( progress_tool )
  {
    for( i=4; i--; )
      progress_tool->step();
  }

  //find all overlapping curves between the bodies
  std::map<Curve*, DLIList<Curve*>* > curve_to_list_map;
  DLIList< DLIList<Curve*>* > lists_of_overlapping_curves;
  //for each list, curves 2-n overlap the first curve
  
  SurfaceOverlapTool::instance()->find_overlapping_curves( body_sms, lists_of_overlapping_curves,
                                                           curve_to_list_map );

  //increment 8%
  if( progress_tool )
  {
    for( i=9; i--; )
      progress_tool->step();
  }

  //remove idential groups and sub-groups of mergeable curves from
  //overlapping groups.  For example:  If curves A,B,C,D are in a mergeable group
  //the following groups will have no need for imprinting, becasue the curves are
  //mergeable all
  //A,B,C,D
  //B,C,D  is a sub-group
  //C,D  is a sub-group
  //These groups need not be considered for imprinting

  std::map<Curve*, DLIList<Curve*>*>::iterator curve_iter;
  lists_of_mergeable_curves.reset();

  DLIList< DLIList<Curve*>* > lists_to_delete;

  for( i=lists_of_mergeable_curves.size(); i--; )
  {
    DLIList<Curve*> *mergeable_curve_list = lists_of_mergeable_curves.get_and_step();

    for( j=mergeable_curve_list->size(); j--; )
    {
      Curve *tmp_curve = mergeable_curve_list->get_and_step();
      curve_iter = curve_to_list_map.find( tmp_curve );

      if( curve_iter == curve_to_list_map.end() )
        continue;

      DLIList<Curve*>* overlapping_curves = curve_iter->second;

      //make sure that list hasn't already been removed
      if( lists_to_delete.move_to( overlapping_curves ) )
        continue;

      //don't want to perturb the list...make a copy
      DLIList<Curve*> copy_list = *mergeable_curve_list;

      //if the list of overlapping curve is the same as the mergeable curve list,
      //remove the list list of overlapping curves

      if( copy_list == *overlapping_curves )
      {
        //remove overlapping_curve_list from the list
        lists_of_overlapping_curves.move_to( overlapping_curves );
        lists_of_overlapping_curves.change_to( NULL );

        //put in list for deletion
        lists_to_delete.append( overlapping_curves );

        //remove overlapping_curve_list from the map
        curve_to_list_map.erase( curve_iter );
        continue;
      }
      else if (overlapping_curves->size() < copy_list.size() )
      {
        //If every curve in overlapping list is in the merge list,
        //remove that overlapping list.
        int list_size_before = copy_list.size();
        copy_list.merge_unique( *overlapping_curves );
        if( list_size_before == copy_list.size() )
        {
          //remove overlapping_curve_list from the list
          lists_of_overlapping_curves.move_to( overlapping_curves );
          lists_of_overlapping_curves.change_to( NULL );

          //put in list for deletion
          lists_to_delete.append( overlapping_curves );
        }
      }
    }
  }

  //clean up lists
  for( i=lists_of_mergeable_curves.size(); i--; )
    delete lists_of_mergeable_curves.get_and_step();
  lists_of_overlapping_curves.remove_all_with_value( NULL );

  for( i=lists_to_delete.size(); i--; )
    delete lists_to_delete.get_and_step();

  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;

  DLIList<BodySM*> imprint_body_order;
  lists_of_overlapping_curves.reset();
  //Now setup multi-map for imprinting vertices onto bodies
  //This is how the curve - vertex imprinting should be done:
  // -imprint vertices of curve 1 in list onto curves 2 - n
  // -imprint vertcies of curve 2 - n onto curve 1
  std::multimap<BodySM*, CubitVector* > body_point_imprint_map;
  std::multimap<BodySM*, CubitVector* >::iterator tmp_iter, upper_iter;

  for( i=lists_of_overlapping_curves.size(); i--; )
  {
    DLIList<Curve*> *curve_list = lists_of_overlapping_curves.get_and_step();
    curve_list->reset();
    Curve *first_curve = curve_list->get_and_step();

    DLIList<Point*> curve_1_points;
    first_curve->points( curve_1_points );
    Point *s_point1 = curve_1_points.get_and_step();
    Point *e_point1 = curve_1_points.get_and_step();

    BodySM *first_curve_body = first_curve->bodysm();

    while( curve_list->size() > 1 )
    {
      Curve *other_curve = curve_list->pop();
      DLIList<Point*> curve_2_points;
      other_curve->points( curve_2_points );
      Point *s_point2 = curve_2_points.get_and_step();
      Point *e_point2 = curve_2_points.get_and_step();

      BodySM *other_curve_body = other_curve->bodysm();

      CubitVector close_pt;
      first_curve->closest_point_trimmed( s_point2->coordinates(), close_pt );

      if( close_pt.distance_between( s_point1->coordinates() ) > tolerance &&
          close_pt.distance_between( e_point1->coordinates() ) > tolerance )
      {
        //are there any other points within tolerance to this point
        //to be imprinted onto this body?
        tmp_iter = body_point_imprint_map.find( first_curve_body );
        bool consider_point = true;
        if( tmp_iter != body_point_imprint_map.end() )
        {
          CubitVector *tmp_vec = tmp_iter->second;
          upper_iter = body_point_imprint_map.upper_bound( tmp_iter->first );
          
          for(; tmp_iter!=upper_iter; ++tmp_iter)
          {
            tmp_vec = tmp_iter->second;
            if( close_pt.distance_between( *tmp_vec ) < tolerance )
            {
              consider_point = false;
              break;
            }
          }
        }

        if( consider_point ) 
        {
          //imprint body of curve1 with vertex at this location
          CubitVector *tmp_vec = new CubitVector( close_pt );
          body_point_imprint_map.insert( std::multimap<BodySM*,
              CubitVector*>::value_type( first_curve_body, tmp_vec ));
          imprint_body_order.append_unique( first_curve_body );
        }
      }

      first_curve->closest_point_trimmed( e_point2->coordinates(), close_pt );

      if( close_pt.distance_between( s_point1->coordinates() ) > tolerance &&
          close_pt.distance_between( e_point1->coordinates() ) > tolerance )
      {
        //are there any other points within tolerance to this point
        //to be imprinted onto this body?
        tmp_iter = body_point_imprint_map.find( first_curve_body );
        bool consider_point = true;
        if( tmp_iter != body_point_imprint_map.end() )
        {
          CubitVector *tmp_vec = tmp_iter->second;
          upper_iter = body_point_imprint_map.upper_bound( tmp_iter->first );
          
          for(; tmp_iter!=upper_iter; ++tmp_iter)
          {
            tmp_vec = tmp_iter->second;
            if( close_pt.distance_between( *tmp_vec ) < tolerance )
            {
              consider_point = false;
              break;
            }
          }
        }

        if( consider_point ) 
        {
          //imprint body of curve1 with vertex at this location
          CubitVector *tmp_vec = new CubitVector( close_pt );
          body_point_imprint_map.insert( std::multimap<BodySM*,
              CubitVector*>::value_type( first_curve_body, tmp_vec ));
          imprint_body_order.append_unique( first_curve_body );
        }
      }

      other_curve->closest_point_trimmed( s_point1->coordinates(), close_pt );

      if( close_pt.distance_between( s_point2->coordinates() ) > tolerance &&
          close_pt.distance_between( e_point2->coordinates() ) > tolerance )
      {
        //are there any other points within tolerance to this point
        //to be imprinted onto this body?
        tmp_iter = body_point_imprint_map.find( other_curve_body );
        bool consider_point = true;
        if( tmp_iter != body_point_imprint_map.end() )
        {
          CubitVector *tmp_vec = tmp_iter->second;
          upper_iter = body_point_imprint_map.upper_bound( tmp_iter->first );
          
          for(; tmp_iter!=upper_iter; ++tmp_iter)
          {
            tmp_vec = tmp_iter->second;
            if( close_pt.distance_between( *tmp_vec ) < tolerance )
            {
              consider_point = false;
              break;
            }
          }
        }

        if( consider_point ) 
        {
          //imprint body of curve1 with vertex at this location
          CubitVector *tmp_vec = new CubitVector( close_pt );
          body_point_imprint_map.insert( std::multimap<BodySM*,
              CubitVector*>::value_type( other_curve_body, tmp_vec ));
          imprint_body_order.append_unique( other_curve_body );
        }
      }

      other_curve->closest_point_trimmed( e_point1->coordinates(), close_pt );

      if( close_pt.distance_between( s_point2->coordinates() ) > tolerance &&
          close_pt.distance_between( e_point2->coordinates() ) > tolerance )
      {
        //Are there already any other points within tolerance to this point
        //to be imprinted onto this body?
        //If so, don't imprint this point.
        tmp_iter = body_point_imprint_map.find( other_curve_body );
        bool consider_point = true;
        if( tmp_iter != body_point_imprint_map.end() )
        {
          CubitVector *tmp_vec = tmp_iter->second;
          upper_iter = body_point_imprint_map.upper_bound( tmp_iter->first );
          
          for(; tmp_iter!=upper_iter; ++tmp_iter )
          {
            tmp_vec = tmp_iter->second;
            if( close_pt.distance_between( *tmp_vec ) < tolerance )
            {
              consider_point = false;
              break;
            }
          }
        }

        if( consider_point ) 
        {
          //imprint body of curve1 with vertex at this location
          CubitVector *tmp_vec = new CubitVector( close_pt );
          body_point_imprint_map.insert( std::multimap<BodySM*,
              CubitVector*>::value_type( other_curve_body, tmp_vec ));
          imprint_body_order.append_unique( other_curve_body );
        }
      }
    }
  }

  //Imprint Vertices-onto-Body
  imprint_body_order.reset();
  std::multimap<BodySM*, CubitVector*>::iterator point_iter, last_point_iter;
  for( i=imprint_body_order.size(); i--; )
  {
    BodySM *tmp_body_sm = imprint_body_order.get_and_step();
    DLIList<BodySM*> body_to_imprint;
    body_to_imprint.append( tmp_body_sm );

    DLIList<CubitVector*> points_to_imprint;

    point_iter = body_point_imprint_map.find( tmp_body_sm );
    last_point_iter = body_point_imprint_map.upper_bound( point_iter->first );

    for(; point_iter!=last_point_iter; ++point_iter )
      points_to_imprint.append( point_iter->second );

    DLIList<BodySM*> new_bodies;
    CubitStatus status = imprint( body_to_imprint, points_to_imprint, new_bodies, false, new_tbs, att_tbs );
    if( new_bodies.size() )
      new_body_sms.append( new_bodies.get() );

    for( j=points_to_imprint.size(); j--; )
      delete points_to_imprint.get_and_step();
  }

  //increment 4%
  if( progress_tool )
  {
    for( i=4; i--; )
      progress_tool->step();
  }

  return CUBIT_SUCCESS;
}

CubitStatus AcisModifyEngine::imprint_overlapping_surfaces( DLIList<BodySM*> &body_sms,
                                                            DLIList<BodySM*> &new_body_sms,
                                                            ProgressTool *progress_tool) const
{
  //Find all mergeable surfaces between the bodies
  DLIList< DLIList<Surface*>* > lists_of_mergeable_surfaces;
  CubitStatus status = MergeTool::instance()->
      find_only_mergeable_surfaces( body_sms, lists_of_mergeable_surfaces );

  if( progress_tool )
    progress_tool->percent( 0.24 );

  //Find all overlapping surfaces between the bodies
  DLIList<Surface*> overlapping_surfaces1;
  DLIList<Surface*> overlapping_surfaces2;
  SurfaceOverlapTool::instance()->
    find_overlapping_surfaces( body_sms, overlapping_surfaces1, overlapping_surfaces2 );

  if( progress_tool )
    progress_tool->percent( 0.57 );

  //We only want the non-mergeable overlapping surfaces.
  //Remove mergable surfaces from overlapping pairs.
  int i, j, k;
  for( i=lists_of_mergeable_surfaces.size(); i--; )
  {
    DLIList<Surface*> *tmp_surface_list = lists_of_mergeable_surfaces.get_and_step();
    for( j=tmp_surface_list->size(); j--; )
    {
      Surface *tmp_surface = tmp_surface_list->get_and_step();
      while( overlapping_surfaces1.move_to( tmp_surface ) )
        overlapping_surfaces1.change_to(NULL);
      while( overlapping_surfaces2.move_to( tmp_surface ) )
        overlapping_surfaces2.change_to(NULL);
    }
    delete tmp_surface_list;
  }

  if( progress_tool )
    progress_tool->percent( 0.59 );

  //Surface A might be mergeable/overlapping with B and just ever-so-slightly
  //overlapping with C. Get rid of C or other overlapping surfaces
  //that might be in the list.
  overlapping_surfaces1.reset();
  overlapping_surfaces2.reset();
  for( i=overlapping_surfaces1.size(); i--; )
  {
    Surface *tmp_surface1 = overlapping_surfaces1.get();
    Surface *tmp_surface2 = overlapping_surfaces2.get();
    
    if( tmp_surface1 != NULL )
    {
      if( tmp_surface2 == NULL )
        overlapping_surfaces1.change_to(NULL);
    }
    else
    {
      if( tmp_surface2 != NULL )
        overlapping_surfaces2.change_to(NULL);
    }

    overlapping_surfaces1.step();
    overlapping_surfaces2.step();
  }

  overlapping_surfaces1.remove_all_with_value(NULL);
  overlapping_surfaces2.remove_all_with_value(NULL);

  if ( overlapping_surfaces1.size() != overlapping_surfaces2.size() )
  {
    PRINT_ERROR("Number of overlapping pairs of surfaces is wrong.\n");
    return CUBIT_FAILURE;
  }

  if( progress_tool )
    progress_tool->percent( 0.62 );

  // Make sure we are working with ACIS surfaces and not
  // something at some higher level.
  DLIList<Surface*> surfs_to_append1;
  DLIList<Surface*> surfs_to_append2;
  for(j=overlapping_surfaces1.size(); j--;)
  {
    Surface *s1 = overlapping_surfaces1.get();
    Surface *s2 = overlapping_surfaces2.get();

    GeometryQueryEngine *gqe1 = s1->get_geometry_query_engine();
    GeometryQueryEngine *gqe2 = s2->get_geometry_query_engine();

    DLIList<TopologyBridge*> tbs1, tbs2;
    gqe1->get_underlying_surfaces(s1, tbs1);
    gqe2->get_underlying_surfaces(s2, tbs2);
    int num_surfs1 = tbs1.size();
    int num_surfs2 = tbs2.size();
    if(num_surfs1 > 0 || num_surfs2 > 0)
    {
      if(num_surfs1 == 0)
      {
        for(k=tbs2.size(); k--;)
        {
          surfs_to_append1.append(s1);
          surfs_to_append2.append(dynamic_cast<Surface*>(tbs2.get_and_step()));
        }
      }
      else if(num_surfs2 == 0)
      {
        for(k=tbs1.size(); k--;)
        {
          surfs_to_append2.append(s2);
          surfs_to_append1.append(dynamic_cast<Surface*>(tbs1.get_and_step()));
        }
      }
      else
      {
        for(k=tbs1.size(); k--;)
        {
          int u;
          Surface *tmp1 = dynamic_cast<Surface*>(tbs1.get_and_step());
          for(u=num_surfs2; u--;)
          {
            surfs_to_append1.append(tmp1);
            surfs_to_append2.append(dynamic_cast<Surface*>(tbs2.get_and_step()));
          }
        }
      }
      overlapping_surfaces1.change_to(NULL);
      overlapping_surfaces2.change_to(NULL);
    }

    overlapping_surfaces1.step();
    overlapping_surfaces2.step();
  }
  overlapping_surfaces1.remove_all_with_value(NULL);
  overlapping_surfaces2.remove_all_with_value(NULL);
  for(k=surfs_to_append1.size(); k--;)
  {
    overlapping_surfaces1.append(surfs_to_append1.get_and_step());
    overlapping_surfaces2.append(surfs_to_append2.get_and_step());
  }


  //Now that we have non-mergeable overlapping surface pairs, determine which
  //curves of surface1 to imprint onto surface2 and vice versa.
  //We eliminate curves between the 2 surfaces that overlap.
  DLIList<Surface*> surface_ordering_list;
  std::multimap<Surface*, Curve*> surface_curve_imprint_map;
  overlapping_surfaces1.reset();
  overlapping_surfaces2.reset();

  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;

  for( i=overlapping_surfaces1.size(); i--; )
  {
    Surface *overlapping_surface1 = overlapping_surfaces1.get_and_step();
    Surface *overlapping_surface2 = overlapping_surfaces2.get_and_step();

    //get all the curves of the surfaces
    DLIList<Curve*> curves2, curves1;
    overlapping_surface2->curves( curves2 );
    overlapping_surface1->curves( curves1 );

    // Make sure we are working with ACIS curves 
    // and not something at some higher level.
    DLIList<Curve*> curves_to_remove;
    DLIList<Curve*> curves_to_append;
    curves2.reset();
    for(j=curves2.size(); j--;)
    {
      Curve *c = curves2.get_and_step();
      GeometryQueryEngine *gqe = c->get_geometry_query_engine();
      DLIList<TopologyBridge*> tbs;
      gqe->get_underlying_curves(c, tbs);
      if(tbs.size() > 0)
      {
        curves_to_remove.append(c);
        for(k=tbs.size(); k--;)
          curves_to_append.append(dynamic_cast<Curve*>(tbs.get_and_step()));
      }
    }
    curves2 -= curves_to_remove;
    curves2 += curves_to_append;
    curves_to_remove.clean_out();
    curves_to_append.clean_out();
    curves1.reset();
    for(j=curves1.size(); j--;)
    {
      Curve *c = curves1.get_and_step();
      GeometryQueryEngine *gqe = c->get_geometry_query_engine();
      DLIList<TopologyBridge*> tbs;
      gqe->get_underlying_curves(c, tbs);
      if(tbs.size() > 0)
      {
        curves_to_remove.append(c);
        for(k=tbs.size(); k--;)
          curves_to_append.append(dynamic_cast<Curve*>(tbs.get_and_step()));
      }
    }
    curves1 -= curves_to_remove;
    curves1 += curves_to_append;

    //Get bounding box for later
    CubitBox bb_surf = overlapping_surface1->bounding_box();

    //remove overlapping curves between the 2 surfaces
    for( j=curves2.size(); j--; )
    {
      Curve *tmp_curve = curves2.get();

      int k;
      bool overlap = false;
      for( k=curves1.size(); k--; )
      {
        if( curves1.get() )
        {
          Curve *tmp_curve_2 = curves1.get();

          if( SurfaceOverlapTool::instance()->
            check_overlap( tmp_curve, tmp_curve_2, NULL))
          {
            //They should be about spatially equal too!
//            CubitSense dummy_sense;
//            if( MergeTool::instance()->about_spatially_equal( curves1.get(), curves2.get(),
//                      dummy_sense, GeometryQueryTool::get_geometry_factor() ) )
//            {
              curves1.change_to(NULL);
              curves2.change_to(NULL);
              overlap = true;
              break;
//            }
          }
        }
        curves1.step();
      }
      if( overlap )
      {
        curves2.step();
        continue;
      }

      double distance;
      CubitVector v1, v2;

      //bounding box check between curve and surface
      CubitBox bb_curve = tmp_curve->bounding_box();
      if( bb_surf.overlap( tolerance, bb_curve ) )
      {
        //measure distance between curve and surface
        overlapping_surface1->get_geometry_query_engine()->entity_entity_distance(
            overlapping_surface1, tmp_curve, v1, v2, distance );

        if( distance < tolerance )
        {
          //check multimap to make sure that this curve has not already
          //been added (from another surface that shares it)
          //or that is overlaps another curve in the multimap
          std::pair< std::multimap<Surface*, Curve* >::iterator,
                     std::multimap<Surface*, Curve* >::iterator > surface_curve_pair;
          surface_curve_pair = surface_curve_imprint_map.equal_range( overlapping_surface1 );

          std::multimap<Surface*, Curve*>::iterator iter = surface_curve_pair.first;
          std::multimap<Surface*, Curve*>::iterator end = surface_curve_pair.second;
          bool match_found = false;
          while( iter != end )
          {
            if( iter->second == tmp_curve )
            {
              match_found = true;
              break;
            }

            if( SurfaceOverlapTool::instance()->check_overlap( iter->second, tmp_curve, NULL ) )
            {
              match_found = true;
              break;
            }

            iter++;
          }

          if( match_found == false )
          {
            surface_curve_imprint_map.insert( std::multimap<Surface*, Curve*>::value_type(
                        overlapping_surface1, tmp_curve ));
            surface_ordering_list.append_unique( overlapping_surface1 );
          }
        }
      }
      curves2.step();
    }

    bb_surf = overlapping_surface2->bounding_box();
    curves1.remove_all_with_value( NULL );
    curves2.remove_all_with_value( NULL );
    //insert body and copied curves into list
    for( j=curves1.size(); j--; )
    {
      Curve *tmp_curve = curves1.get();

      double distance;
      CubitVector v1, v2;
      CubitBox bb_curve = tmp_curve->bounding_box();
      //bounding box check between curve and surface
      if( bb_surf.overlap( tolerance, bb_curve ) )
      {
        //measure distance between curve and surface
        overlapping_surface2->get_geometry_query_engine()->entity_entity_distance(
          overlapping_surface2, tmp_curve, v1, v2, distance );

        if( distance < tolerance )
        {
          //check multimap to make sure that this curve has not already
          //been added (from another surface that shares it)
          //or that is overlaps another curve in the multimap
          std::pair< std::multimap<Surface*, Curve* >::iterator,
                     std::multimap<Surface*, Curve* >::iterator > surface_curve_pair;
          surface_curve_pair = surface_curve_imprint_map.equal_range( overlapping_surface2 );

          std::multimap<Surface*, Curve*>::iterator end = surface_curve_pair.second;
          std::multimap<Surface*, Curve*>::iterator iter = surface_curve_pair.first;
          bool match_found = false;
          while( iter != end )
          {
            if( iter->second == tmp_curve )
            {
              match_found = true;
              break;
            }

            if( SurfaceOverlapTool::instance()->check_overlap( iter->second, tmp_curve, NULL ) )
            {
              match_found = true;
              break;
            }

            iter++;
          }

          if( match_found == false )
          {
            surface_curve_imprint_map.insert( std::multimap<Surface*, Curve*>::value_type(
                                overlapping_surface2, tmp_curve ));
            surface_ordering_list.append_unique( overlapping_surface2 );
          }

        }
      }
      curves1.step();
    }
  }

  if( progress_tool )
    progress_tool->percent( 0.64 );

  //We cannot use the original curves for imprinting....they'll be
  //getting modified.  Make copies of them here and build a new multimap
  std::multimap<Surface*, Curve*>::iterator surf_curve_iter, surf_curve_last;
  std::multimap<Surface*, Curve*> final_surface_curve_imprint_map;
  for( surf_curve_iter = surface_curve_imprint_map.begin();
       surf_curve_iter != surface_curve_imprint_map.end(); )
  {
    Surface *tmp_surf = surf_curve_iter->first;
    surf_curve_last = surface_curve_imprint_map.upper_bound( tmp_surf );

    for(; surf_curve_iter!=surf_curve_last; ++surf_curve_iter )
    {
      Curve *copied_curve = make_Curve( surf_curve_iter->second );
      final_surface_curve_imprint_map.insert( std::multimap<Surface*,
                  Curve*>::value_type( tmp_surf, copied_curve ));
    }
  }

  if( progress_tool )
    progress_tool->percent( 0.67 );

  double percentage_per_surface = 0;
  double delta = 0;

  if( progress_tool )
  {
    if( surface_ordering_list.size() )
      percentage_per_surface = 0.16 / surface_ordering_list.size();
    else
      progress_tool->percent( 0.83 );

  }

  //Now iterate of each surface that will be imprinted with curves.
  surface_ordering_list.reset();
  bool print_errors = false;
  for( i=surface_ordering_list.size(); i--; )
  {
    if( progress_tool )
    {
      delta += percentage_per_surface;
      progress_tool->percent( 0.67 + delta );
    }

    Surface *tmp_surface = surface_ordering_list.get_and_step();

    surf_curve_iter = final_surface_curve_imprint_map.find( tmp_surface );
    tmp_surface = surf_curve_iter->first;

    DLIList<Surface*> surface_list;
    surface_list.append( tmp_surface );

    DLIList<Curve*> curves_to_imprint;
    surf_curve_last = final_surface_curve_imprint_map.upper_bound( surf_curve_iter->first );
    //For each curve that might be imprinted on the surface...

     for(; surf_curve_iter!=surf_curve_last; ++surf_curve_iter )
     {
       DLIList<Curve*> tmp_curve_list1, tmp_curve_list2; 
       Curve *original_curve = surf_curve_iter->second;
       tmp_curve_list1.append( original_curve ); 

       //project
       CubitStatus projected = project_edges( surface_list, tmp_curve_list1, tmp_curve_list2, print_errors );

       if( projected == CUBIT_FAILURE )
       {
         AcisQueryEngine::instance()->delete_solid_model_entities( original_curve );
         continue;
       }
 
       Curve *projected_curve = tmp_curve_list2.get(); 
       
       //if midpoint of projected is far from original curve, continue 
       CubitVector original_curve_mid_point = original_curve->center_point();
       if( original_curve_mid_point.distance_between( projected_curve->center_point() ) > tolerance ) 
       {
         AcisQueryEngine::instance()->delete_solid_model_entities( original_curve );
         AcisQueryEngine::instance()->delete_solid_model_entities( projected_curve );
         continue;
       }
        
       //now intersect
       tmp_curve_list1.clean_out();
       CubitStatus intersected =  curve_surface_intersection( tmp_surface, projected_curve, tmp_curve_list1 ); 
      
       //delete the curves
       AcisQueryEngine::instance()->delete_solid_model_entities( original_curve );
       AcisQueryEngine::instance()->delete_solid_model_entities( projected_curve );

       if( intersected  == CUBIT_FAILURE) 
         continue;
       
       curves_to_imprint += tmp_curve_list1;
     }

     //filter out small curves
     int j;
     for( j=curves_to_imprint.size(); j--; )
     {
       Curve *tmp_curve = curves_to_imprint.get();
       double curve_length = tmp_curve->measure();
       if( curve_length < tolerance )
       {
         curves_to_imprint.change_to( NULL );
         AcisQueryEngine::instance()->delete_solid_model_entities( tmp_curve ); 
       }
       curves_to_imprint.step();
     }
     curves_to_imprint.remove_all_with_value( NULL );

    BodySM *new_body_sm = NULL;

    //Now imprint the curves onto the surface
    CubitStatus status;
    if( curves_to_imprint.size() )
    {
/*
      GfxDebug::clear();
      AcisDrawTool::instance()->draw_surface( tmp_surface, 4, true);
      GfxDebug::mouse_xforms();
      int kkk;
      PRINT_INFO("number of curve = %d\n", curves_to_imprint.size() );
      for( kkk=curves_to_imprint.size(); kkk--; )
      {
        AcisDrawTool::instance()->draw_curve( curves_to_imprint.get_and_step(), kkk+5 );
        GfxDebug::mouse_xforms();
      }
*/
      status = embed_curves_into_surface( tmp_surface, curves_to_imprint, new_body_sm );
    }

    if( new_body_sm == NULL )
      continue;
    new_body_sms.append_unique( new_body_sm );

    //Delete all the curves that were imprinted.  We're responsible for deleting them.
    for( j=curves_to_imprint.size(); j--; )
      AcisQueryEngine::instance()->delete_solid_model_entities( curves_to_imprint.get_and_step() );

  }
  return CUBIT_SUCCESS;
}


//WARNING!!! Only use this function if you KNOW that the curves you want
//to imprint line on the surface.....otherwise you'll get problems.
//api_embed_wire_in_faces returns unexpected results if the curve
//does not lie on the surface.
CubitStatus AcisModifyEngine::embed_curves_into_surface( Surface *surface,
                                                           DLIList<Curve*> &curves_to_imprint,
                                                           BodySM *&new_body ) const
{
  if( curves_to_imprint.size() == 0 )
    return CUBIT_FAILURE;

  BODY *copied_BODY = NULL;
  DLIList<FACE*> imprint_FACE_list;
  DLIList<SurfaceACIS*> acis_surfs;
  SurfaceACIS *tmp_acis_surf = CAST_TO( surface, SurfaceACIS );

  if( tmp_acis_surf == NULL )
    return CUBIT_FAILURE;

  FACE *FACE_ptr = tmp_acis_surf->get_FACE_ptr();
  BODY *BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY( FACE_ptr );
  BodySM *old_body_sm = AcisQueryEngine::instance()->get_body_sm_of_ENTITY( FACE_ptr );

  acis_surfs.append( tmp_acis_surf );
  if( get_copied_FACES_of_body( acis_surfs, imprint_FACE_list, copied_BODY ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( imprint_FACE_list.size() == 0 )
    return CUBIT_FAILURE;

  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;
  if( tolerance < 0.001 )
    tolerance = 0.001;

  int i;
  EDGE *EDGEs[ 1 ];
  ENTITY_LIST wire_BODYs;
#ifdef BOYD16
  outcome tmp_result;
#endif
  for( i=0; i<curves_to_imprint.size(); i++ )
  {
    //copy
    EDGE *EDGE_ptr = AcisQueryEngine::get_EDGE( curves_to_imprint.get_and_step() );

    EDGE *copied_EDGE_ptr = NULL;
    outcome tmp_result = api_edge( EDGE_ptr, copied_EDGE_ptr );
    ATTRIB_CUBIT_OWNER::remove_cubit_owner( (ENTITY *)copied_EDGE_ptr,
                                              CUBIT_TRUE );
    EDGEs[0] = copied_EDGE_ptr;
    BODY *tmp_BODY = NULL;
    tmp_result = api_make_ewire( 1, EDGEs, tmp_BODY );
    wire_BODYs.add( (ENTITY*)tmp_BODY );
  }

  int number_wire_bodies = wire_BODYs.count();

  if( number_wire_bodies == 0 )
    return CUBIT_FAILURE;

  BODY *wire_BODY_ptr = (BODY*)wire_BODYs[0];

  ENTITY_LIST *face_list = new ENTITY_LIST;
  FACE* imprint_FACE = imprint_FACE_list.get();
  face_list->add( (ENTITY*)imprint_FACE );

  ENTITY_LIST temp_face_list;
  api_get_faces( copied_BODY, temp_face_list );
  int number_faces_before = temp_face_list.count();

  API_BEGIN
    result = api_embed_wire_in_faces( wire_BODY_ptr, copied_BODY, face_list, tolerance );
  API_END

  for(i=1; i<number_wire_bodies; i++ )
  {
#ifdef BOYD16
    bool embedded = false;
#endif
    wire_BODY_ptr = (BODY*)wire_BODYs[i];

    //get all the faces on the body
    ENTITY_LIST FACEs;
    api_get_faces( copied_BODY, FACEs );

    //if number of faces on body did not change, no face got split, and we can
    //split the same FACE again.
    if( number_faces_before == FACEs.count()  )
    {
      result = api_embed_wire_in_faces( wire_BODY_ptr, copied_BODY, face_list, tolerance );
      continue;
    }

    //need to determine which face to embed curve onto
    //may be able to remove this later on and just embed into body
    //just embed the remining wires into the body
    //result = api_embed_wire_in_faces( wire_BODY_ptr, copied_BODY, NULL, tolerance );
    //if( !result.ok() )
    //  AcisQueryEngine::instance()->ACIS_API_error (result);

    //get midpoint of curve to embed
    ENTITY_LIST wire_edge_list;
    api_get_edges( wire_BODY_ptr, wire_edge_list );
    EDGE *wire_edge = (EDGE*)wire_edge_list[0];
    SPAposition mid_point = wire_edge->mid_pos();

    //get position at fraction 0.79 along curve
#ifdef BOYD16
    SPAposition start_point = wire_edge->start()->geometry()->coords();
#endif
    SPAparameter s_param = wire_edge->start_param();
    SPAparameter e_param = wire_edge->end_param();

    //SPAinterval param_interval = wire_edge->param_range();
    SPAparameter new_parameter = s_param + ((e_param - s_param ) * 0.79);

    if (wire_edge->sense() == REVERSED )
      new_parameter = -(new_parameter);

    SPAposition point_on_curve;
    wire_edge->geometry()->equation().eval(new_parameter, point_on_curve );

    //find surface which contains point at parameter 0.79 and midpoint
    //loop to find appropriate face
    ENTITY *temp_ENT;
    FACEs.init();
    FACE *embed_FACE = NULL;
    while( (temp_ENT = FACEs.next()) != NULL )
    {
      FACE *tmp_FACE = (FACE*)temp_ENT;
      SPAposition point1;
      api_find_cls_ptto_face(mid_point, tmp_FACE, point1 );
      if ( distance_to_point( point1, mid_point ) < 0.0005 )
      {
        api_find_cls_ptto_face(point_on_curve, tmp_FACE, point1 );

        if ( distance_to_point( point1, point_on_curve ) < 0.0005 )
        {
          embed_FACE = tmp_FACE;
          break;
        }
      }
    }
    if( embed_FACE )
    {
      face_list->clear();
      face_list->add( (ENTITY*)embed_FACE );
      result = api_embed_wire_in_faces( wire_BODY_ptr, copied_BODY, face_list, tolerance );
    }
  }
  
  AcisModifyEngine::instance()->cleanup_slivers( copied_BODY );

  //check to see if anything actually happened
  int num_vol_bef, num_face_bef, num_edge_bef, num_vertex_bef;
  CubitStatus status1 = AcisQueryEngine::instance()->number_ENTITIES( BODY_ptr, num_vol_bef,
                                           num_face_bef, num_edge_bef,
                                           num_vertex_bef );

  int num_vol_aft, num_face_aft, num_edge_aft, num_vertex_aft;
  CubitStatus status2 = AcisQueryEngine::instance()->number_ENTITIES( copied_BODY, num_vol_aft,
                                           num_face_aft, num_edge_aft,
                                           num_vertex_aft );

  if( status1 == CUBIT_FAILURE  ||
      status2 == CUBIT_FAILURE )
  {
    return CUBIT_FAILURE;
  }

  if( num_vol_bef == num_vol_aft &&
      num_face_bef == num_face_aft &&
      num_edge_bef == num_edge_aft &&
      num_vertex_bef == num_vertex_aft)
  {
    AcisQueryEngine::instance()->delete_ACIS_BODY(copied_BODY, CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }

  BodySM* new_body_ptr = get_new_Body( old_body_sm, BODY_ptr, copied_BODY,
                                       false, CUBIT_TRUE );

  new_body = new_body_ptr;
  return CUBIT_SUCCESS;
}

CubitStatus AcisModifyEngine::curve_surface_intersection( Surface *surface,
                                            Curve* curve,
                                            DLIList<Curve*> &new_curves ) const
{
  FACE *FACE_ptr = AcisQueryEngine::get_FACE( surface );
  EDGE *EDGE_ptr = AcisQueryEngine::get_EDGE( curve );

  ENTITY_LIST *ent_list = new ENTITY_LIST;

  outcome result = api_edfa_int( EDGE_ptr, FACE_ptr, ent_list );

  if( !result.ok() || ent_list->count() == 0 )
  {
    delete ent_list;
    return CUBIT_FAILURE;
  }

  ENTITY *tmp_ent;
  ent_list->init();
  while( (tmp_ent=ent_list->next()) != NULL )
  {
    if( is_EDGE( tmp_ent ) )
    {
      EDGE_ptr = static_cast<EDGE*>(tmp_ent);
      Curve *curve_ptr = AcisQueryEngine::instance()->populate_topology_bridges(EDGE_ptr);
      new_curves.append( curve_ptr );
    }
  }
  delete ent_list;

  if( new_curves.size() == 0 )
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

void AcisModifyEngine::get_new_ENTITIES(ENTITY *top_ENTITY, DLIList<ENTITY*> &new_ENTITIES) const
{
  ENTITY_LIST ENTITIES;
  int i;
  
  // vertices
  api_get_vertices( top_ENTITY, ENTITIES);
  for (i = 0; i < ENTITIES.count(); i++)
  {
    if (!ATTRIB_CUBIT_OWNER::cubit_owner(ENTITIES[i]))
    {
      new_ENTITIES.append_unique(ENTITIES[i]);
    }
  }
  ENTITIES.clear();

  // edges
  api_get_edges( top_ENTITY, ENTITIES);
  for (i = 0; i < ENTITIES.count(); i++)
  {
    if (!ATTRIB_CUBIT_OWNER::cubit_owner(ENTITIES[i]))
    {
      new_ENTITIES.append_unique(ENTITIES[i]);
    }
  }
  ENTITIES.clear();
/*
  // faces
  api_get_faces( top_ENTITY, ENTITIES);
  for (i = 0; i < ENTITIES.count(); i++)
  {
    if (!ATTRIB_CUBIT_OWNER::cubit_owner(ENTITIES[i]))
    {
      new_ENTITIES.append_unique(ENTITIES[i]);
    }
  }
  ENTITIES.clear();
  */
}

void AcisModifyEngine::get_att_ENTITIES(ENTITY *top_ENTITY, DLIList<ENTITY*> &att_ENTITIES, char *att_name) const
{
  ENTITY_LIST ENTITIES;
  int i;
  // vertices
  api_get_vertices( top_ENTITY, ENTITIES);
  for (i = 0; i < ENTITIES.count(); i++)
  {
    ATTRIB_SNL_SIMPLE *attribute =
      (ATTRIB_SNL_SIMPLE *) find_attrib(ENTITIES[i], 
                                        ATTRIB_SNL_TYPE,
                                        ATTRIB_SNL_SIMPLE_TYPE);
    for(;attribute != NULL;attribute = (ATTRIB_SNL_SIMPLE *)
          find_next_attrib(attribute,
                          ATTRIB_SNL_TYPE,
                          ATTRIB_SNL_SIMPLE_TYPE))
    {
      if( attribute->attribute_name() == att_name )
        att_ENTITIES.append_unique(ENTITIES[i]);
    }
  }
  ENTITIES.clear();

  // edges
  api_get_edges( top_ENTITY, ENTITIES);
  for (i = 0; i < ENTITIES.count(); i++)
  {
    ATTRIB_SNL_SIMPLE *attribute =
      (ATTRIB_SNL_SIMPLE *) find_attrib(ENTITIES[i], 
                                        ATTRIB_SNL_TYPE,
                                        ATTRIB_SNL_SIMPLE_TYPE);
    for(;attribute != NULL;attribute = (ATTRIB_SNL_SIMPLE *)
          find_next_attrib(attribute,
                          ATTRIB_SNL_TYPE,
                          ATTRIB_SNL_SIMPLE_TYPE))
    {
      if( attribute->attribute_name() == att_name )
        att_ENTITIES.append_unique(ENTITIES[i]);
    }
  }
  ENTITIES.clear();
/*
  // faces
  api_get_faces( top_ENTITY, ENTITIES);
  for (i = 0; i < ENTITIES.count(); i++)
  {
    ATTRIB_SNL_SIMPLE *attribute =
      (ATTRIB_SNL_SIMPLE *) find_attrib(ENTITIES[i], 
                                        ATTRIB_SNL_TYPE,
                                        ATTRIB_SNL_SIMPLE_TYPE);
    for(;attribute != NULL;attribute = (ATTRIB_SNL_SIMPLE *)
          find_next_attrib(attribute,
                          ATTRIB_SNL_TYPE,
                          ATTRIB_SNL_SIMPLE_TYPE))
    {
      if( attribute->attribute_name() == att_name )
        att_ENTITIES.append_unique(ENTITIES[i]);
    }
  }
  ENTITIES.clear();
  */
}

CubitStatus AcisModifyEngine::cleanup_slivers( BODY *body_to_cleanup )
{
   double sliver_curve_tol = GeometryQueryTool::get_sliver_curve_cleanup_tolerance();
   double sliver_surface_tol = GeometryQueryTool::get_sliver_surface_cleanup_tolerance();

   ENTITY_LIST dummy_list;
   api_detect_short_edges( body_to_cleanup, dummy_list, sliver_curve_tol, true );
   api_detect_sliver_faces( body_to_cleanup, dummy_list, sliver_surface_tol, true );

   return CUBIT_SUCCESS;
}

CubitStatus AcisModifyEngine::remove_curve_slivers( BodySM *body_sm, 
                                                    double lengthlimit ) const
{
  //get the BODY
  BodyACIS *tmp_body_acis = CAST_TO(body_sm, BodyACIS );
  BODY *tmp_BODY = tmp_body_acis->get_BODY_ptr();

  //copy the BODY 
  BODY *copy_BODY = AcisModifyEngine::instance()->copy_BODY( tmp_BODY, false );

  //remove the slivers
  ENTITY_LIST removed_curves;
  outcome result = api_detect_short_edges( copy_BODY, removed_curves, lengthlimit, true );

  if( !result.ok() )
  {
    PRINT_ERROR("Problem remove sliver curves\n");
    AcisQueryEngine::instance()->ACIS_API_error (result);
    api_delent( (ENTITY*)copy_BODY );
    return CUBIT_FAILURE;
  }

  TopologyEntity *te = body_sm->topology_entity();
  RefEntity *ref_ent = CAST_TO( te, RefEntity ); 
  PRINT_INFO("Removed %d sliver curves on body %d\n",
            removed_curves.count(), ref_ent->id() );

  if( removed_curves.count() == 0 )
  {
    api_delent( (ENTITY*)copy_BODY );
    return CUBIT_SUCCESS;
  }
  
  get_new_Body(body_sm, tmp_BODY, copy_BODY, CUBIT_FALSE);

  return CUBIT_SUCCESS;
}
