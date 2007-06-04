#ifdef ACIS_HEALER

//- Class:       AcisHealerTool
//- Description: Functions to heal bad ACIS geometry
// 
//- Owner: Steve Storm      
//- Checked by:
//- Version: $Id:

// ********** BEGIN ACIS INCLUDES             **********
#ifdef UNIX_PORT
#include "undefwin_.h"
#endif

#if CUBIT_ACIS_VERSION < 1100
#include "kernel/acis.hxx"
#include "kernel/kernapi/api/api.hxx"
#include "kernel/kernapi/api/kernapi.hxx"
#include "kernel/kerndata/top/body.hxx"
#include "kernel/kerndata/top/face.hxx"

#include "boolean/kernapi/api/boolapi.hxx"
#include "euler/kernapi/api/eulerapi.hxx"
#include "healhusk/heal_api/heal_api.hxx"

#include "kernel/kerndata/data/entity.hxx"
#include "kernel/kerndata/lists/lists.hxx"

#include "kernel/kerndata/top/alltop.hxx" // Required for is_TEDGE

#include "healhusk/util/bhl_opt.hxx"
#include "healhusk/util/tols.hxx"
#include "healhusk/tgtspl/net_heal.hxx"
#include "healhusk/attrib/at_aggr.hxx"
#include "healhusk/attrib/aggrsimg.hxx"
#include "healhusk/attrib/aggrstch.hxx"
#include "healhusk/attrib/hmaster.hxx"
#include "healhusk/attrib/entsimg.hxx"
#include "healhusk/attrib/hanalsol.hxx"
#include "healhusk/attrib/huvsolv.hxx"
#include "healhusk/simgeom/simgeom.hxx"
# include "healhusk/stitch/stitch.hxx"
# include "healhusk/stitch/stch_utl.hxx"

// Needed for advanced feedback
#include "healhusk/attrib/aggrgbld.hxx"
#include "healhusk/attrib/edgmbld.hxx"
#include "healhusk/attrib/vegmbld.hxx"
#include "healhusk/attrib/cegmbld.hxx"
#include "healhusk/attrib/fagmbld.hxx"
#include "healhusk/attrib/entstch.hxx"
#include "healhusk/attrib/curgmbld.hxx"
#include "healhusk/attrib/surgmbld.hxx"
#include "healhusk/attrib/pcgmbld.hxx"
#include "healhusk/attrib/shlgmbld.hxx"
#include "healhusk/attrib/lpgmbld.hxx"
#include "healhusk/attrib/lmpgmbld.hxx"
#include "healhusk/util/traverse.hxx"
#include "healhusk/attrib/hsharped.hxx"
#include "healhusk/attrib/hadvspl.hxx"
#include "healhusk/attrib/hreblend.hxx"
#else
#include "acis.hxx"
#include "api.hxx"
#include "kernapi.hxx" 
#include "body.hxx"
#include "face.hxx" 

#include "boolapi.hxx" 
#include "eulerapi.hxx" 
#include "heal_api.hxx"

#include "entity.hxx" 
#include "lists.hxx"

#include "alltop.hxx" // Required for is_TEDGE

# include "bhl_opt.hxx"
# include "tols.hxx"
# include "net_heal.hxx"
# include "at_aggr.hxx"
# include "aggrsimg.hxx"
# include "aggrstch.hxx"
# include "hmaster.hxx"
# include "entsimg.hxx"
# include "hanalsol.hxx"
# include "huvsolv.hxx"
# include "simgeom.hxx"
# include "stitch_heal.hxx"
# include "stch_utl_heal.hxx"


// Needed for advanced feedback
# include "aggrgbld.hxx"
# include "edgmbld.hxx"
# include "vegmbld.hxx"
# include "cegmbld.hxx"
# include "fagmbld.hxx"
# include "entstch.hxx"
# include "curgmbld.hxx"
# include "surgmbld.hxx"
# include "pcgmbld.hxx"
# include "shlgmbld.hxx"
# include "lpgmbld.hxx"
# include "lmpgmbld.hxx"
# include "traverse_heal.hxx"
# include "hsharped.hxx"
# include "hadvspl.hxx"
#endif

#ifdef UNIX_PORT
# include "defwin_.h"
#endif
// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********

#include "AcisHealerTool.hpp"
#include "AcisToolUtil.hpp"
#include "AcisQueryEngine.hpp"
#include "AcisModifyEngine.hpp"
#include "GeometryHealerTool.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "CastTo.hpp"
#include "CubitMessage.hpp"
#include "CubitUtil.hpp"
#include "RefEntityFactory.hpp"

#include "attrib_cubit_owner.hpp"
#include "BodyACIS.hpp"
#include "SurfaceACIS.hpp"
#include "RefEntity.hpp"
#include "Body.hpp"
#include "RefEdge.hpp"
#include "CoEdge.hpp"
#include "Loop.hpp"
#include "RefFace.hpp"
#include "Shell.hpp"
#include "RefVolume.hpp"
#include "RefVertex.hpp"

#include "DLIList.hpp"

// ********** END CUBIT INCLUDES              **********

AcisHealerTool* AcisHealerTool::instance_ = 0;

static logical doPreprocess = TRUE;    // This is normally always done
static logical userControlPreprocess = FALSE;

//Keep track of user defined simplify tols
static double simplifyTol = 0;           // User simplify tolerance
static double defaultSimplifyTol = CUBIT_DBL_MAX;
static logical simplifyTolSw = FALSE;    // Is user setting the simplify tolerance?
static logical doSimplify = TRUE;        // User setting for simplify (do it or not)
static logical userControlSimplify = FALSE;   // Is user controlling whether to simplify?
void set_user_defined_simplify_tols(BODY *body);

//Keep track of user defined stitch tols
static double stitchMinTol = 0;
static double defaultStitchMinTol = CUBIT_DBL_MAX;
static logical stitchMinTolSw = FALSE;
static double stitchMaxTol = 0;
static double defaultStitchMaxTol = CUBIT_DBL_MAX;
static logical stitchMaxTolSw = FALSE;
static logical doStitch = TRUE;
static logical userControlStitch = FALSE;
void set_user_defined_stitch_tols(BODY *body);

//logical hh_stitching_type = FALSE;

//Keep track of user defined geombuild tols
static double geomBuildTol = 0;
static double defaultGeomBuildTol = CUBIT_DBL_MAX;
static logical geomBuildTolSw = FALSE;
static logical doGeomBuild = TRUE;
static logical userControlGeomBuild = FALSE;

static double analyticTol = 0;
static double defaultAnalyticTol = CUBIT_DBL_MAX;
static logical analyticTolSw = FALSE;
static logical doAnalytic = TRUE;
static logical userControlAnalytic = FALSE;

static double isosplineTol = 0;
static double defaultIsosplineTol = CUBIT_DBL_MAX;
static logical isosplineTolSw = FALSE;
static logical doIsospline = TRUE;
static logical userControlIsospline = FALSE;

static double reblendTol = 0;
static double defaultReblendTol = CUBIT_DBL_MAX;
static logical reblendTolSw = FALSE;
static double reblendClassifyTol = 0;
static double defaultReblendClassifyTol = CUBIT_DBL_MAX;
static logical reblendClassifyTolSw = FALSE;
static logical doReblend = TRUE;
static logical userControlReblend = FALSE;

static logical doSharpEdge = TRUE;
static logical userControlSharpEdge = FALSE;

static logical doGenericSpline = TRUE;
static logical userControlGenericSpline = FALSE;

static logical doWrapup = TRUE;             // This is normally always done
static logical userControlWrapup = FALSE;

void set_user_defined_geombuild_tols(BODY *body);
//End user defined geombuild

static logical doPostprocess = TRUE;             // This is normally always done
static logical userControlPostprocess = FALSE;

static int hdebug = 0;
// Set to 1 for debugging

// Method: instance
// provides access to the unique model for this execution.
// sets up this instance on first access
AcisHealerTool* AcisHealerTool::instance()
{
  if( instance_ == 0 )
  {
    new AcisHealerTool();
    api_initialize_healing();
  }
  return instance_;
}

AcisHealerTool::AcisHealerTool() 
{
   instance_ = this;

   cleanAtt = CUBIT_FALSE;

   showMethod = 1; // 0-none, 1-highlight, 2-draw

   showSummary = CUBIT_TRUE;
   //- Flag to determine whether list a summary when showing bad geometry

   showDetails = CUBIT_FALSE;
   //- Flag to determine whether list details when showing bad geometry

   showOnHeal = CUBIT_TRUE;
   //- Flag to determine whether to show bad geometry after healing

   checkVolOnHeal = CUBIT_FALSE;
   checkSurfOnHeal = CUBIT_FALSE;
   checkCurveOnHeal = CUBIT_FALSE;
   volLimit = 0.1;
   surfLimit = 0.1;
   curveLimit = 0.1;
   //- Allow the healer to automatically check for small geometry after
   //- healing.

   showBadVertices = CUBIT_TRUE;
   showBadCurves = CUBIT_TRUE;
   showBadCoEdges = CUBIT_FALSE;
   showBadLoops = CUBIT_FALSE;
   showBadSurfaces = CUBIT_FALSE;
   showBadShells = CUBIT_FALSE;
   showBadVolumes = CUBIT_FALSE;
   showBadBodies = CUBIT_FALSE;
   //- What to show during analysis and show

   incPreprocess = CUBIT_TRUE; // This is normally always done
   incSimplify = CUBIT_TRUE;                        
   incStitch = CUBIT_TRUE;
   incGeombuild = CUBIT_TRUE;
   incAnalytic = CUBIT_TRUE; // The next 6 are ignored if Geombuild is off
   incIsospline = CUBIT_TRUE;
   incReblend = CUBIT_FALSE; // Until implemented by Spatial
   incSharpedge = CUBIT_TRUE;
   incGenericspline = CUBIT_TRUE;
   incWrapup = CUBIT_TRUE;      // Normally always done
   incPostprocess = CUBIT_TRUE; // Normally always done
   //- Controls for incremental healing steps

   GeometryHealerTool::instance()->add_ghe(this);
}

AcisHealerTool::~AcisHealerTool()
{
   api_terminate_healing();
   instance_ = NULL;
}

CubitStatus AcisHealerTool::init_BODY_for_healing( BODY *BODY_ptr )
{
   outcome result;
   
   result = api_hh_init_body_for_healing( BODY_ptr );
   if( !result.ok() )
   {
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
      return CUBIT_FAILURE;
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::analyze_BODY_for_healing( BODY *BODY_ptr )
{
   outcome result;

   result = api_hh_simplify_analyze( BODY_ptr );
   if( !result.ok())
   {
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
      return CUBIT_FAILURE;
   }

   result = api_hh_stitch_analyze( BODY_ptr );
   if( !result.ok())
   {
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
      return CUBIT_FAILURE;
   }
   
   //result = api_hh_geombuild_analyze( BODY_ptr );
   result = api_hh_analyze_body( BODY_ptr );
   if( !result.ok())
   {
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
      return CUBIT_FAILURE;
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::end_BODY_for_healing( BODY *BODY_ptr )
{
  outcome result;
  
  result = api_hh_end_body_for_healing( BODY_ptr );
  if( !result.ok() )
  {
    AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::heal_BODY( BODY *BODY_ptr, int& percent_before,
                                       int& percent_after, 
                                       int& num_splines_simplified,
                                       CubitBoolean make_tolerant, 
                                       FILE* logfile_ptr )
{
   CubitStatus status = CUBIT_SUCCESS;
   CubitBoolean did_geombuild = CUBIT_FALSE;
   outcome result;
   
   // Setup logfile ptr
   if( logfile_ptr != NULL )
      hh_set_bhl_log_file( logfile_ptr );

   // set up callback functions (prototypes in heal_api.hxx)
   hh_set_simplify_opt_callback (&set_user_defined_simplify_tols);
   hh_set_stitch_opt_callback (&set_user_defined_stitch_tols);
   hh_set_geombuild_opt_callback (&set_user_defined_geombuild_tols);

   if( !userControlPreprocess || (userControlPreprocess && doPreprocess) )
   {
     result = api_hh_preprocess( BODY_ptr );
     if( !result.ok() )
     {
       PRINT_ERROR( "error in preprocess\n" );
       AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
       return CUBIT_FAILURE;
     }
   }
   
   if( !userControlSimplify || (userControlSimplify && doSimplify) )
   {
#if CUBIT_ACIS_VERSION < 1100
      result = api_hh_simplify_auto( BODY_ptr);
#else
      hh_simplify_options simp_opts;
      simp_opts.set_do_curve_simplification(0);
      result = api_hh_simplify_auto( BODY_ptr, &simp_opts );
#endif 
      if( !result.ok() )
      {
         PRINT_ERROR( "error in simplify\n" );
         AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
         return CUBIT_FAILURE;
      }

		// Give some info as to simplification results
      ATTRIB_HH_AGGR_SIMPLIFY *att = find_aggr_simplify ( BODY_ptr );
      if( att )
      {
        num_splines_simplified = att->num_input_splines() - att->num_final_splines();
      }
   }

   if( !userControlStitch || (userControlStitch && doStitch) )
   {
      result = api_hh_stitch_auto( BODY_ptr );
      if( !result.ok() )
      {
         PRINT_ERROR( "error in stitch\n" );
         AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
         return CUBIT_FAILURE;
      }
   }
   
   if( !userControlGeomBuild || (userControlGeomBuild && doGeomBuild ) )
   {
      result = api_hh_geombuild_auto( BODY_ptr );
      if( !result.ok() )
      {
         PRINT_ERROR( "error in geombuild\n" );
         AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
         return CUBIT_FAILURE;
      }
      else
        did_geombuild = CUBIT_TRUE;
   }
   
   if( !userControlPostprocess || (userControlPostprocess && doPostprocess ) )
   {
     result = api_hh_postprocess( BODY_ptr );
     if( !result.ok() )
     {
       PRINT_ERROR( "error in postprocess\n" );
       AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
       return CUBIT_FAILURE;
     }
   }

   // reset log file
   hh_reset_bhl_log_file();

   // Retrieve healing results
   
   if( percentage_goodgeom_after( BODY_ptr, percent_after ) == CUBIT_FAILURE )
   {
     // No results available
     percent_before = -1;
     percent_after = -1;
   }
   else
   {
     if( did_geombuild )
       percentage_goodgeom_before( BODY_ptr, percent_before );
     else
       percent_before = percent_after; // Assume they are the same
   }

   // Make tolerant edges if requested.  If this is done before
   // retrieving results (above), the percent_before and
   // percent_after come out the same - at percent_after amount.
   if( make_tolerant == CUBIT_TRUE )
   {
      result = api_hh_make_tolerant( BODY_ptr );
      if( !result.ok())
      {
         AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
         status = CUBIT_FAILURE;
      }
      else
         status = CUBIT_SUCCESS;
   }

   result = api_reset_boxes( BODY_ptr );
   if( !result.ok())
   {
     AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
     status = CUBIT_FAILURE;
   }
   else
     status = CUBIT_SUCCESS;

   return status;
}

CubitBoolean AcisHealerTool::is_healer_engine( const TopologyBridge* tb_ptr )
{
   if (CAST_TO(const_cast<TopologyBridge*>(tb_ptr), AcisBridge)) return CUBIT_TRUE;
   else return CUBIT_FALSE;
}

CubitBoolean AcisHealerTool::is_healer_engine( const TopologyEntity* te_ptr )
{
   TopologyBridge *tb_ptr = te_ptr->bridge_manager()->topology_bridge();
   return is_healer_engine(tb_ptr);
}

CubitStatus AcisHealerTool::auto_heal_bodies( DLIList<Body*> &body_list,
                                              DLIList<Body*> &new_body_list,
                                              DLIList<TopologyEntity*> &bad_geometry,
                                              CubitBoolean rebuild,
                                              CubitBoolean keep_old,
                                              CubitBoolean make_tolerant,
                                              FILE* logfile_ptr )
{

  return heal_bodies( body_list, new_body_list, bad_geometry, rebuild, keep_old, 
                      make_tolerant, logfile_ptr );
}

BODY* AcisHealerTool::auto_rebuild_BODY( BODY *BODY_ptr, int& percent_before,
                                         int& percent_after, 
                                         CubitBoolean make_tolerant,
                                         FILE* logfile_ptr )
{
   // Unhook all the faces, heal each one individually, then stitch them back 
   // together, heal the final body

   int j, k, num_faces = 0;
   DLIList<FACE*> FACE_list;
   FACE *FACE_ptr;
   outcome result;

   int local_percent_before,
      local_percent_after,
      num_splines_simplified;
   
   // This is just to analyze the BODY to get the percentage of goodgeom
   // before healing, and to populate the logfile.  This slows down
   // healing, since this analysis is only used for this printout.
   // However, since the rebuild option is rarely used, this is okay.

   if( logfile_ptr ) 
      hh_set_bhl_log_file (logfile_ptr);
   
   if( init_BODY_for_healing( BODY_ptr ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "initializing Body for healing\n" );
      return NULL;
   }

   if( analyze_BODY_for_healing( BODY_ptr ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Unable to analyze BODY for healing\n" );
      end_BODY_for_healing( BODY_ptr );
      return NULL;
   }
   
   // Find percentage of goodgeometry before healing
   percentage_goodgeom_before( BODY_ptr, percent_before );
 
   // Clean out the attributes since not needed anymore
   if( end_BODY_for_healing( BODY_ptr ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "ending Body for healing\n" );
      return NULL;
   }

   // Turn off logfile, so that we don't get information about each face
   if( logfile_ptr ) 
      hh_set_bhl_log_file (NULL);

   // Get a list of all this BODY's FACE's
   AcisModifyEngine::instance()->get_acis_query_engine()->get_FACEs( (ENTITY *)BODY_ptr, FACE_list );
   num_faces = FACE_list.size();
   
   // Now we have the FACES - unhook them from the BODY.  Keep track of new
   // BODIES that are created as this is done.

   PRINT_INFO( " Unhooking and healing each individual surface...\n" );

   BODY** BODIES = new BODY*[num_faces];
   for( j=0; j<num_faces; j++ )
      BODIES[j] = NULL;
   for( j=0; j<num_faces; j++ )
   {
      FACE_ptr = FACE_list.get_and_step();
      BODY *new_BODY_ptr;
      result = api_unhook_face( FACE_ptr, new_BODY_ptr );
      if( !result.ok() )
      {
        AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
        for( k=0; k<=j; k++ )
        {
          if( BODIES[k] )
            AcisModifyEngine::instance()->get_acis_query_engine()->delete_ACIS_BODY(BODIES[k]);
        }
        delete [] BODIES;
        PRINT_ERROR( "   Face unhooking during rebuild of Body didn't work\n" );
        return NULL;
      }
      
      BODIES[j] = new_BODY_ptr;
      
      // Heal this sheet body by itself.  This sometimes helps final results.
      if( init_BODY_for_healing( new_BODY_ptr ) == CUBIT_FAILURE )
      {
         PRINT_WARNING( "   Unable to heal one of the detached surfaces\n" );
         continue;
      }
      if( heal_BODY( new_BODY_ptr, local_percent_before, local_percent_after,
         num_splines_simplified ) 
         == CUBIT_FAILURE )
      {
         PRINT_WARNING( "   Unable to heal one of the detached surfaces\n" );
      }
      if( end_BODY_for_healing( new_BODY_ptr ) == CUBIT_FAILURE )
      {
         PRINT_WARNING( "   Unable to cleanup healing attributes on one of the detached surfaces\n" );
      }
   }
   
   BODY *master = BODIES[0];
   
   PRINT_INFO( " Combining the unhooked surfaces and healing the new Body...\n" );

   // Combine the new BODIES
   for( j=1; j<num_faces; j++ )
   {
      // This deletes the first body & puts combined result into second
      result = api_combine_body(BODIES[j], master);
      
      if (!result.ok()) 
      {
         //PRINT_ERROR("In AcisHealerTool::auto_rebuild_body\n");
        AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error(result, "combine Bodies");
        AcisModifyEngine::instance()->get_acis_query_engine()->delete_ACIS_BODY(master);
         for( k=j; k<num_faces; k++ )
         {
           AcisModifyEngine::instance()->get_acis_query_engine()->delete_ACIS_BODY(BODIES[k]);
         }
         return NULL;
      }
   }
   
   delete [] BODIES;

   // Turn back on logfile to get information about healed body
   if( logfile_ptr ) 
   {
      hh_set_bhl_log_file( logfile_ptr );
      fprintf( logfile_ptr, "\nNOW STITCHING THE BODY BACK TOGETHER :\n" );
      fprintf( logfile_ptr, "**************************************\n\n" );
   }

   if( init_BODY_for_healing( master ) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Unable to initialize combined BODY for healing\n" );
      return NULL;
   }
   
   if( heal_BODY( master, local_percent_before, percent_after,
                  num_splines_simplified, make_tolerant) == CUBIT_FAILURE )
   {
      PRINT_ERROR( "Unable to heal combined BODY\n" );
      return NULL;
   }

   return master;
}

CubitStatus AcisHealerTool::overall_printout( BODY* BODY_ptr, FILE* file_ptr )
{
   
   ATTRIB_HH_AGGR_GEOMBUILD* att = find_aggr_geombuild( BODY_ptr );
   
   if (NULL == att)
      return CUBIT_FAILURE;
   else
   {
      bhl_anal_geometry_results anal_results;
      att->check(&anal_results);

      att->print( file_ptr );
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::analysis_printout( BODY* BODY_ptr, FILE* file_ptr)
{
   
   ATTRIB_HH_AGGR_GEOMBUILD* att = find_aggr_geombuild( BODY_ptr );
   
   if (NULL == att)
      return CUBIT_FAILURE;
   else
   {
      bhl_anal_geometry_results anal_results;
      att->check(&anal_results);

      att->print_analyze( file_ptr );
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::calculate_printout( BODY* BODY_ptr, FILE* file_ptr)
{
   
   ATTRIB_HH_AGGR_GEOMBUILD* att = find_aggr_geombuild( BODY_ptr );
   
   if (NULL == att)
      return CUBIT_FAILURE;
   else
   {
      bhl_anal_geometry_results anal_results;
      att->check(&anal_results);

      att->print_calculate( file_ptr );
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::fix_printout( BODY* BODY_ptr, FILE* file_ptr)
{
   
   ATTRIB_HH_AGGR_GEOMBUILD* att = find_aggr_geombuild( BODY_ptr );
   
   if (NULL == att)
      return CUBIT_FAILURE;
   else
   {
      bhl_anal_geometry_results anal_results;
      att->check(&anal_results);

      att->print_fix( file_ptr );
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::percentage_goodgeom_before( BODY* BODY_ptr, 
                                                        int& percent_goodgeom )
{
   ATTRIB_HH_AGGR_GEOMBUILD* att = find_aggr_geombuild( BODY_ptr );
   
   if (att == NULL)
   {
     percent_goodgeom = -1;
     return CUBIT_FAILURE;
   }
   else
   {
     bhl_anal_geometry_results in_anal_results;
     
     att->check(&in_anal_results);
     
     in_anal_results = ((ATTRIB_HH_AGGR_GEOMBUILD*)att)->analysis_results();
     
     percent_goodgeom = in_anal_results.healed_percentage;
   }
   
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::percentage_goodgeom_after( BODY* BODY_ptr, 
                                                       int& percent_goodgeom )
{
  ATTRIB_HH_AGGR_GEOMBUILD* att = find_aggr_geombuild( BODY_ptr );
  
  if (att == NULL)
  {
    percent_goodgeom = -1;
    return CUBIT_FAILURE;
  }
  else
  {
    bhl_anal_geometry_results out_anal_results;
   
    API_BEGIN
    att->check(&out_anal_results);
    API_END
    
    out_anal_results = ((ATTRIB_HH_AGGR_GEOMBUILD*)att)->output_analysis_results();
    
    percent_goodgeom = out_anal_results.healed_percentage;
  }
  
  return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::analyze_badgeom( DLIList<Body*> &body_list,
                                             DLIList<TopologyEntity*> &bad_geometry,
                                             FILE* logfile_ptr )
{
  CubitStatus status;
  status = get_badgeom( CUBIT_TRUE, body_list, CUBIT_FALSE, logfile_ptr, bad_geometry );
  if( cleanAtt )
    clean_attributes( body_list );
  return status;
}

CubitStatus AcisHealerTool::get_badgeom( DLIList<Body*> &body_list, 
                                         DLIList<TopologyEntity*> &bad_geometry )
{
  DLIList<Body*> analyzed;
  DLIList<Body*> not_analyzed;
  
  Body* body_ptr;
  BODY* BODY_ptr;

   // Loop through & find which bodies have been analyzed
  for( int i=0; i<body_list.size(); i++ )
  {
    body_ptr = body_list.get_and_step();
    BODY_ptr =  BodyACIS::get_BODY_ptr(body_ptr);

    if( BODY_ptr == NULL )
    {
      PRINT_ERROR( "Couldn't find ACIS BODY from CUBIT Body %d\n", body_ptr->id() );
      not_analyzed.append_unique( body_ptr );
      continue;
    }

    if( is_analyzed( BODY_ptr ) == CUBIT_TRUE )
      analyzed.append_unique( body_ptr );
    else
      not_analyzed.append_unique( body_ptr );
  }

  // Now do the work
  if( not_analyzed.size() == body_list.size() )
  {
    PRINT_ERROR( "none of the input bodies have been analyzed\n"
      "     Hint: try 'healer analyze body <id_list>' first\n" );
    return CUBIT_FAILURE;
  }

  if( not_analyzed.size() > 0 )
  {
    char pre[64];
    sprintf( pre, "WARNING: %d bodies haven't been analyzed: ", 
      not_analyzed.size() );
    DLIList<CubitEntity*> temp_list;
    CAST_LIST_TO_PARENT(not_analyzed, temp_list);    
    CubitUtil::list_entity_ids( pre, temp_list );
  }

  // Finally, show the bad geometry.
  return get_badgeom( CUBIT_FALSE, analyzed, CUBIT_FALSE, NULL, bad_geometry );
}

CubitStatus AcisHealerTool::get_tcurves( DLIList<Body*> &body_list,
                                         DLIList<RefEdge*> &t_curves  )
{

  if( body_list.size() == 0 )
    return CUBIT_SUCCESS;
  
  for( int i=0; i<body_list.size(); i++ )
  {
    Body* body_ptr = body_list.get_and_step();
    DLIList<RefEdge*> ref_edge_list;
    body_ptr->ref_edges( ref_edge_list );
    for( int j=0; j<ref_edge_list.size(); j++ )
    {
      RefEdge* ref_edge_ptr = ref_edge_list.get_and_step();
      if( ref_edge_ptr->is_tolerant() == CUBIT_TRUE )
      {
        t_curves.append( ref_edge_ptr );
      }
    }
  }
  if( showMethod == 1 || showMethod == 2 )
    PRINT_INFO( "%d tolerant curves found\n", t_curves.size() );
  else
  {
    DLIList<CubitEntity*> temp_list;
    CAST_LIST_TO_PARENT(t_curves, temp_list);    
    CubitUtil::list_entity_ids( 
      "Tolerant curves in requested bodies: ", temp_list );
  }
  
  return CUBIT_FAILURE;
}

CubitStatus AcisHealerTool::heal_incremental( DLIList<Body*> &body_list, 
                                              DLIList<Body*> &new_body_list,
                                              DLIList<TopologyEntity*> &bad_geometry,
                                              double simplify_tol,
                                              double stitch_min_tol, 
                                              double stitch_max_tol,
                                              double geom_build_tol,
                                              double analytic_tol,
                                              double isospline_tol,
                                              double reblend_classify_tol,
                                              double reblend_tol,
                                              CubitBoolean keep_old, 
                                              CubitBoolean make_tolerant,
                                              FILE* logfile_ptr )
{
  // Set-up the processes and tolerances.  The tolerances sent to this
  // function can be used to override the defaults.  If these aren't
  // supplied, and the defaults aren't set, the healer uses intelligent
  // tolerances.

  userControlPreprocess = TRUE;
  if( incPreprocess )
    doPreprocess = TRUE;
  else
    doPreprocess = FALSE;

  userControlSimplify = TRUE; // User is controlling whether to simplify
  if( incSimplify )
  {
    doSimplify = TRUE;
    set_simplify_tol( simplify_tol );
  }
  else
    doSimplify = FALSE;
  
  userControlStitch = TRUE;
  if( incStitch )
  {
    doStitch = TRUE;
    set_stitch_min_tol( stitch_min_tol );
    set_stitch_max_tol( stitch_max_tol );
  }
  else
    doStitch = FALSE;
  
  userControlGeomBuild = TRUE;
  if( incGeombuild && (incAnalytic || incIsospline || incReblend || incSharpedge || incGenericspline) )
  {
    doGeomBuild = TRUE;
    set_geombuild_tol( geom_build_tol );
    set_analytic_tol( analytic_tol );
    set_isospline_tol( isospline_tol );
    set_reblend_classify_tol( reblend_classify_tol );
    set_reblend_tol( reblend_tol );
  }
  else
    doGeomBuild = FALSE;
  
  // If geombuild, check each step of geombuild
  if( doGeomBuild==TRUE )
  {
    userControlAnalytic = TRUE;
    if( incAnalytic )
    {
      doAnalytic = TRUE;
      set_analytic_tol( analytic_tol );
    }
    else
      doAnalytic = FALSE;
    
    userControlIsospline = TRUE;
    if( incIsospline )
    {
      doIsospline = TRUE;
      set_isospline_tol( isospline_tol );
    }
    else
      doIsospline = FALSE;
    
    userControlReblend = TRUE;
    if( incReblend )
    {
      doReblend = TRUE;
      set_reblend_classify_tol( reblend_classify_tol );
      set_reblend_tol( reblend_tol );
    }
    else
      doReblend= FALSE;
    
    userControlSharpEdge = TRUE;
    if( incSharpedge )
      doSharpEdge = TRUE;
    else
      doSharpEdge= FALSE;
    
    userControlGenericSpline = TRUE;
    if( incGenericspline )
      doGenericSpline = TRUE;
    else
      doGenericSpline= FALSE;
    
    userControlWrapup = TRUE;
    if( incWrapup )
      doWrapup = TRUE;
    else
      doWrapup = FALSE;
  }
  
  userControlPostprocess = TRUE;
  if( incPostprocess )
    doPostprocess = TRUE;
  else
    doPostprocess = FALSE;

  // Give some information on what steps are being done
  incremental_presummary( body_list, logfile_ptr, keep_old );
  
  CubitStatus status = heal_bodies( body_list, new_body_list, 
     bad_geometry, CUBIT_FALSE, keep_old, make_tolerant, logfile_ptr );
  
  // Reset local switches to default (as autoheal would do it)
  reset_switches();
  
  return status;
}

// ********************************************************************
//
//  callback routine to set user defined simplification tolerance
//
// ********************************************************************
void set_user_defined_simplify_tols(BODY *body)
{
   if ( simplifyTolSw == TRUE )
   {
      ATTRIB_HH_AGGR_SIMPLIFY * att = (ATTRIB_HH_AGGR_SIMPLIFY *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_SIMPLIFY_TYPE);
      if (NULL != att){
         att->set_tol(simplifyTol);
      }
   }
   else if( defaultSimplifyTol != CUBIT_DBL_MAX )
   {
      ATTRIB_HH_AGGR_SIMPLIFY * att = (ATTRIB_HH_AGGR_SIMPLIFY *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_SIMPLIFY_TYPE);
      if (NULL != att){
         att->set_tol(defaultSimplifyTol);
      }
   }
   
   if ( userControlSimplify == TRUE )
   {
      ATTRIB_HH_AGGR_SIMPLIFY * att = (ATTRIB_HH_AGGR_SIMPLIFY *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_SIMPLIFY_TYPE);
      if (NULL != att){
         att->set_do_simplify(doSimplify);
      }
   }
}

// ********************************************************************
//
//  callback routine to set user defined stitch tolerances
//
// ********************************************************************
void set_user_defined_stitch_tols(BODY *body)
{
   if ( stitchMinTolSw == TRUE )
   {
      ATTRIB_HH_AGGR_STITCH * att = (ATTRIB_HH_AGGR_STITCH *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_STITCH_TYPE);
      if (NULL != att){
         att->set_min_tol(stitchMinTol);
      }
   }
   else if( defaultStitchMinTol != CUBIT_DBL_MAX )
   {
      ATTRIB_HH_AGGR_STITCH * att = (ATTRIB_HH_AGGR_STITCH *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_STITCH_TYPE);
      if (NULL != att){
         att->set_min_tol(defaultStitchMinTol);
      }
   }
   
   if ( stitchMaxTolSw == TRUE )
   {
      ATTRIB_HH_AGGR_STITCH * att = (ATTRIB_HH_AGGR_STITCH *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_STITCH_TYPE);
      if (NULL != att){
         att->set_max_tol(stitchMaxTol);
      }
   }
   else if( defaultStitchMaxTol != CUBIT_DBL_MAX )
   {
      ATTRIB_HH_AGGR_STITCH * att = (ATTRIB_HH_AGGR_STITCH *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_STITCH_TYPE);
      if (NULL != att){
         att->set_max_tol(defaultStitchMaxTol);
      }
   }
   
   if ( userControlStitch == TRUE )
   {
      ATTRIB_HH_AGGR_STITCH * att = (ATTRIB_HH_AGGR_STITCH *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_STITCH_TYPE);
      if (NULL != att){
         att->set_do_stitch(doStitch);
      }
   }
}

// ********************************************************************
//
//  callback routine to set user defined geombuild tolerances
//
// ********************************************************************
void set_user_defined_geombuild_tols(BODY *body)
{
   if ( userControlGeomBuild == TRUE )
   {
      ATTRIB_HH_AGGR_GEOMBUILD * att = (ATTRIB_HH_AGGR_GEOMBUILD *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_GEOMBUILD_TYPE);
      if (NULL != att){
         att->set_do_geombuild(doGeomBuild);
      }
   }
   if ( geomBuildTolSw == TRUE )
   {
      ATTRIB_HH_AGGR_GEOMBUILD * att = (ATTRIB_HH_AGGR_GEOMBUILD *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_GEOMBUILD_TYPE);
      if (NULL != att){
         att->set_tol(geomBuildTol);
      }
   }
   else if( defaultGeomBuildTol != CUBIT_DBL_MAX )
   {
      ATTRIB_HH_AGGR_GEOMBUILD * att = (ATTRIB_HH_AGGR_GEOMBUILD *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_GEOMBUILD_TYPE);
      if (NULL != att){
         att->set_tol(defaultGeomBuildTol);
      }
   }
   
   if ( analyticTolSw == TRUE )
   {
      ATTRIB_HH_AGGR_ANALYTIC * att = (ATTRIB_HH_AGGR_ANALYTIC *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_ANALYTIC_TYPE);
      if (NULL != att){
         att->set_tol(analyticTol);
      }
   }
   else if ( defaultAnalyticTol != CUBIT_DBL_MAX )
   {
      ATTRIB_HH_AGGR_ANALYTIC * att = (ATTRIB_HH_AGGR_ANALYTIC *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_ANALYTIC_TYPE);
      if (NULL != att){
         att->set_tol(defaultAnalyticTol);
      }
   }
   
   if ( isosplineTolSw == TRUE )
   {
      ATTRIB_HH_AGGR_ISOSPLINE * att = (ATTRIB_HH_AGGR_ISOSPLINE *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_ISOSPLINE_TYPE);
   }
   else if ( defaultIsosplineTol != CUBIT_DBL_MAX )
   {
      ATTRIB_HH_AGGR_ISOSPLINE * att = (ATTRIB_HH_AGGR_ISOSPLINE *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_ISOSPLINE_TYPE);
   }
   
   if ( userControlAnalytic == TRUE )
   {
      ATTRIB_HH_AGGR_ANALYTIC * att = (ATTRIB_HH_AGGR_ANALYTIC *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_ANALYTIC_TYPE);
      if (NULL != att){
         att->set_do_analytic(doAnalytic);
      }
   }
   if ( userControlIsospline == TRUE )
   {
      ATTRIB_HH_AGGR_ISOSPLINE * att = (ATTRIB_HH_AGGR_ISOSPLINE *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_ISOSPLINE_TYPE);
      if (NULL != att){
         att->set_do_isospline(doIsospline);
      }
   }
/*
   if( userControlReblend == TRUE )
   {
      ATTRIB_HH_AGGR_REBLEND * att = (ATTRIB_HH_AGGR_REBLEND *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_REBLEND_TYPE);
      if (NULL != att)
         att->set_do_reblend(doReblend);
   }
   if ( reblendClassifyTolSw == TRUE )
   {
      ATTRIB_HH_AGGR_REBLEND * att = (ATTRIB_HH_AGGR_REBLEND *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_REBLEND_TYPE);
      if (NULL != att){
         att->set_classify_tol(reblendClassifyTol);
      }
   }
   else if( defaultReblendClassifyTol != CUBIT_DBL_MAX )
   {
      ATTRIB_HH_AGGR_REBLEND * att = (ATTRIB_HH_AGGR_REBLEND *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_REBLEND_TYPE);
      if (NULL != att){
         att->set_classify_tol(defaultReblendClassifyTol);
      }
   }

   if ( reblendTolSw == TRUE )
   {
      ATTRIB_HH_AGGR_REBLEND * att = (ATTRIB_HH_AGGR_REBLEND *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_REBLEND_TYPE);
      if (NULL != att){
         att->set_tol(reblendClassifyTol);
      }
   }
   else if( defaultReblendTol != CUBIT_DBL_MAX )
   {
      ATTRIB_HH_AGGR_REBLEND * att = (ATTRIB_HH_AGGR_REBLEND *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_REBLEND_TYPE);
      if (NULL != att){
         att->set_tol(defaultReblendTol);
      }
   }

*/
   if( userControlSharpEdge )
   {
      ATTRIB_HH_AGGR_SHARP_EDGE * att = (ATTRIB_HH_AGGR_SHARP_EDGE *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_SHARP_EDGE_TYPE);
      if (NULL != att)
         att->set_do_sharp_edge(doSharpEdge);
   }

   if( userControlGenericSpline )
   {
      ATTRIB_HH_AGGR_GEN_SPLINE * att = (ATTRIB_HH_AGGR_GEN_SPLINE *)
         find_leaf_attrib(body,ATTRIB_HH_AGGR_GEN_SPLINE_TYPE);
      if (NULL != att)
         att->set_do_gen_spline(doGenericSpline);
   }
}

double AcisHealerTool::get_simplify_tol(ENTITY* ent)
{
   if (simplifyTolSw == TRUE)
      return simplifyTol;
   else if( defaultSimplifyTol != CUBIT_DBL_MAX )
      return defaultSimplifyTol;
   else if( ent != NULL )
   {
      if( ent->identity(BODY_LEVEL)==BODY_TYPE) {
         ATTRIB_HH_AGGR_SIMPLIFY * att = (ATTRIB_HH_AGGR_SIMPLIFY *)
            find_leaf_attrib(ent,ATTRIB_HH_AGGR_SIMPLIFY_TYPE);
         if (NULL != att)
            return att->tol();
      }
      else if( ent->identity(FACE_LEVEL)==FACE_TYPE) {
         ATTRIB_HH_ENT_SIMPLIFY_FACE * att = (ATTRIB_HH_ENT_SIMPLIFY_FACE*)
            find_leaf_attrib(ent,ATTRIB_HH_ENT_SIMPLIFY_FACE_TYPE);
         if (NULL != att)
            return att->tol();
      }
   }
   return -1.0;
}
void AcisHealerTool::set_simplify_tol( double tol )
{
   if( tol != CUBIT_DBL_MAX )
   {
      simplifyTol = tol;
      simplifyTolSw = TRUE;
   }
   else
      simplifyTolSw = FALSE;
}
double AcisHealerTool::get_default_simplify_tol()
{ 
  return defaultSimplifyTol; 
}
void AcisHealerTool::set_default_simplify_tol( double tol )
{
   defaultSimplifyTol = tol;
}

double AcisHealerTool::get_stitch_min_tol(ENTITY* ent)
{
   if (stitchMinTolSw == TRUE)
      return stitchMinTol;
   else if( defaultStitchMinTol != CUBIT_DBL_MAX )
      return defaultStitchMinTol;
   else if( ent != NULL )
   {
      if (stitchMinTolSw == FALSE){
         if( ent->identity(BODY_LEVEL)==BODY_TYPE) {
            ATTRIB_HH_AGGR_STITCH * att = (ATTRIB_HH_AGGR_STITCH *)
               find_leaf_attrib(ent,ATTRIB_HH_AGGR_STITCH_TYPE);
            if (NULL != att)
               return att->min_tol();
         }
      }
   }
   return -1.0;
}
void AcisHealerTool::set_stitch_min_tol( double tol )
{
   if( tol != CUBIT_DBL_MAX )
   {
      stitchMinTol = tol;
      stitchMinTolSw = TRUE;
   }
   else
      stitchMinTolSw = FALSE;
}
double AcisHealerTool::get_default_stitch_min_tol()
{
   return defaultStitchMinTol;
}
void AcisHealerTool::set_default_stitch_min_tol( double tol )
{
   defaultStitchMinTol = tol;
}

double AcisHealerTool::get_stitch_max_tol(ENTITY* ent)
{
   if (stitchMaxTolSw == TRUE)
      return stitchMaxTol;
   else if( defaultStitchMaxTol != CUBIT_DBL_MAX )
      return defaultStitchMaxTol;
   else if( ent != NULL )
   {
      if( ent->identity(BODY_LEVEL)==BODY_TYPE) {
         ATTRIB_HH_AGGR_STITCH * att = (ATTRIB_HH_AGGR_STITCH *)
            find_leaf_attrib(ent,ATTRIB_HH_AGGR_STITCH_TYPE);
         if (NULL != att)
            return att->max_tol();
      }
   }
   return -1.0;
}
void AcisHealerTool::set_stitch_max_tol( double tol )
{
   if( tol != CUBIT_DBL_MAX )
   {
      stitchMaxTol = tol;
      stitchMaxTolSw = TRUE;
   }
   else
      stitchMaxTolSw = FALSE;
}
double AcisHealerTool::get_default_stitch_max_tol()
{
   return defaultStitchMaxTol;
}
void AcisHealerTool::set_default_stitch_max_tol( double tol )
{
   defaultStitchMaxTol = tol;
}

double AcisHealerTool::get_geombuild_tol(ENTITY* ent)
{
   if (geomBuildTolSw == TRUE)
      return geomBuildTol;
   else if( defaultGeomBuildTol != CUBIT_DBL_MAX )
      return defaultGeomBuildTol;
   else if( ent != NULL )
   {
      if( ent->identity(BODY_LEVEL)==BODY_TYPE) {
         ATTRIB_HH_AGGR_GEOMBUILD * att = (ATTRIB_HH_AGGR_GEOMBUILD *)
            find_leaf_attrib(ent,ATTRIB_HH_AGGR_GEOMBUILD_TYPE);
         if (NULL != att)
            return att->tol();
      }
   }
   return -1.0;
}
void AcisHealerTool::set_geombuild_tol( double tol )
{
   if( tol != CUBIT_DBL_MAX )
   {
      geomBuildTol = tol;
      geomBuildTolSw = TRUE;
   }
   else
      geomBuildTolSw = FALSE;
}
double AcisHealerTool::get_default_geombuild_tol()
{
   return defaultGeomBuildTol;
}
void AcisHealerTool::set_default_geombuild_tol( double tol )
{
   defaultGeomBuildTol = tol;
}

double AcisHealerTool::get_analytic_tol(ENTITY* ent)
{
   if (analyticTolSw == TRUE)
      return analyticTol;
   else if( defaultAnalyticTol != CUBIT_DBL_MAX )
      return defaultAnalyticTol;
   else if( ent != NULL )
   {
      if( ent->identity(BODY_LEVEL)==BODY_TYPE) {
         ATTRIB_HH_AGGR_ANALYTIC * att = (ATTRIB_HH_AGGR_ANALYTIC *)
            find_leaf_attrib(ent,ATTRIB_HH_AGGR_ANALYTIC_TYPE);
         if (NULL != att)
            return att->tol();
      }
   }
   return -1.0;
}
void AcisHealerTool::set_analytic_tol( double tol )
{
   if( tol != CUBIT_DBL_MAX )
   {
      analyticTol = tol;
      analyticTolSw = TRUE;
   }
   else
      analyticTolSw = FALSE;
}
double AcisHealerTool::get_default_analytic_tol()
{
   return defaultAnalyticTol;
}
void AcisHealerTool::set_default_analytic_tol( double tol )
{
   defaultAnalyticTol = tol;
}

double AcisHealerTool::get_isospline_tol(ENTITY* ent)
{
   if (isosplineTolSw == TRUE)
      return isosplineTol;
   else if( defaultIsosplineTol != CUBIT_DBL_MAX )
      return defaultIsosplineTol;
   else if( ent != NULL )
   {
      if( ent->identity(BODY_LEVEL)==BODY_TYPE) {
         ATTRIB_HH_AGGR_ISOSPLINE * att = (ATTRIB_HH_AGGR_ISOSPLINE *)
            find_leaf_attrib(ent,ATTRIB_HH_AGGR_ISOSPLINE_TYPE);
      }
   }
   return -1.0;
}
void AcisHealerTool::set_isospline_tol( double tol )
{
   if( tol != CUBIT_DBL_MAX )
   {
      isosplineTol = tol;
      isosplineTolSw = TRUE;
   }
   else
      isosplineTolSw = FALSE;
}
double AcisHealerTool::get_default_isospline_tol()
{
   return defaultIsosplineTol;
}
void AcisHealerTool::set_default_isospline_tol( double tol )
{
   defaultIsosplineTol = tol;
}

double AcisHealerTool::get_reblend_classify_tol(ENTITY* ent)
{ 
   if (reblendClassifyTolSw == TRUE)
      return reblendClassifyTol;
   else if( defaultReblendClassifyTol != CUBIT_DBL_MAX )
      return defaultReblendClassifyTol;
   else if( ent != NULL )
   {
     /*
      if( ent->identity(BODY_LEVEL)==BODY_TYPE) {
         ATTRIB_HH_AGGR_REBLEND * att = (ATTRIB_HH_AGGR_REBLEND *)
            find_leaf_attrib(ent,ATTRIB_HH_AGGR_REBLEND_TYPE);
         if (NULL != att)
            return att->classify_tol();
      }
      */
   }
   return -1.0;  
}
void AcisHealerTool::set_reblend_classify_tol( double tol )
{
   if( tol != CUBIT_DBL_MAX )
   {
      reblendClassifyTol = tol;
      reblendClassifyTolSw = TRUE;
   }
   else
      reblendClassifyTolSw = FALSE;
}
double AcisHealerTool::get_default_reblend_classify_tol()
{
   return defaultReblendClassifyTol;
}
void AcisHealerTool::set_default_reblend_classify_tol( double tol )
{
   defaultReblendClassifyTol = tol;
}

double AcisHealerTool::get_reblend_tol(ENTITY* ent)
{
   if (reblendTolSw == TRUE)
      return reblendTol;
   else if( defaultReblendTol != CUBIT_DBL_MAX )
      return defaultReblendTol;
   else if( ent != NULL )
   {
     /*
      if( ent->identity(BODY_LEVEL)==BODY_TYPE) {
         ATTRIB_HH_AGGR_REBLEND * att = (ATTRIB_HH_AGGR_REBLEND *)
            find_leaf_attrib(ent,ATTRIB_HH_AGGR_REBLEND_TYPE);
         if (NULL != att)
            return att->tol();
      }
      */
   }
   return -1.0;
}
void AcisHealerTool::set_reblend_tol( double tol )
{
   if( tol != CUBIT_DBL_MAX )
   {
      reblendTol = tol;
      reblendTolSw = TRUE;
   }
   else
      reblendTolSw = FALSE;
}
double AcisHealerTool::get_default_reblend_tol()
{
   return defaultReblendTol;
}
void AcisHealerTool::set_default_reblend_tol( double tol )
{
   defaultReblendTol = tol;
}

void AcisHealerTool::reset_switches()
{
   // reset local switches

  userControlPreprocess = FALSE;
  doPreprocess = TRUE;

  simplifyTolSw = FALSE;        // User is not specifying simplify tolerance
  userControlSimplify = FALSE;  // User is not controlling whether or not to simplify
  doSimplify = TRUE;
  
  stitchMinTolSw = FALSE;
  stitchMaxTolSw = FALSE;
  userControlStitch = FALSE;
  doStitch = TRUE;
  
  geomBuildTolSw = FALSE;
  userControlGeomBuild = FALSE;
  doGeomBuild = TRUE;
  
  analyticTolSw = FALSE;
  userControlAnalytic = FALSE;
  doAnalytic = TRUE;
  
  isosplineTolSw = FALSE;
  userControlIsospline = FALSE;
  doIsospline = TRUE;
  
  reblendClassifyTolSw = FALSE;
  reblendTolSw = FALSE;
  userControlReblend = FALSE;
  doReblend = TRUE;
  
  userControlSharpEdge = FALSE;
  doSharpEdge = TRUE;
  
  userControlGenericSpline = FALSE;
  doGenericSpline = TRUE;

  userControlWrapup = FALSE;
  doWrapup = TRUE;

  userControlPostprocess = FALSE;
  doPostprocess = TRUE;
}

void AcisHealerTool::reset_default_tolerances()
{
   defaultSimplifyTol = CUBIT_DBL_MAX;
   defaultStitchMinTol = CUBIT_DBL_MAX;
   defaultStitchMaxTol = CUBIT_DBL_MAX;
   defaultGeomBuildTol = CUBIT_DBL_MAX;
   defaultAnalyticTol = CUBIT_DBL_MAX;
   defaultIsosplineTol = CUBIT_DBL_MAX;
   defaultReblendClassifyTol = CUBIT_DBL_MAX;
   defaultReblendTol = CUBIT_DBL_MAX;
}
void AcisHealerTool::list_default_tolerances()
{
   if( defaultSimplifyTol != CUBIT_DBL_MAX )
      PRINT_INFO( "Healer Default Simplify Tol = %f\n", defaultSimplifyTol );
   else
      PRINT_INFO( "Healer Default Simplify Tol = calculated (normally ~0.0001)\n" );

   if( defaultStitchMinTol != CUBIT_DBL_MAX )
      PRINT_INFO( "Healer Default Stitch Min Tol = %f\n", defaultStitchMinTol );
   else
      PRINT_INFO( "Healer Default Stitch Min Tol = calculated (normally ~10e-5)\n" );

   if( defaultStitchMaxTol != CUBIT_DBL_MAX )
      PRINT_INFO( "Healer Default Stitch Max Tol = %f\n", defaultStitchMaxTol );
   else
      PRINT_INFO( "Healer Default Stitch Max Tol = calculated (normally ~1.0)\n" );

   if( defaultGeomBuildTol != CUBIT_DBL_MAX )
      PRINT_INFO( "Healer Default Geombuild Tol = %f\n", defaultGeomBuildTol );
   else
      PRINT_INFO( "Healer Default Geombuild Tol = calculated (normally ~0.01)\n" );

   if( defaultAnalyticTol != CUBIT_DBL_MAX )
      PRINT_INFO( "Healer Default Analytic Tol = %f\n", defaultAnalyticTol );
   else
      PRINT_INFO( "Healer Default Analytic Tol = calculated (normally ~0.01)\n" );

   if( defaultIsosplineTol != CUBIT_DBL_MAX )
      PRINT_INFO( "Healer Default Isospline Tol = %f\n", defaultIsosplineTol );
   else
      PRINT_INFO( "Healer Default Isospline Tol = calculated (normally ~0.01)\n" );

// Future
//   if( defaultReblendClassifyTol != CUBIT_DBL_MAX )
//      PRINT_INFO( "Healer Default Reblend Classify Tol = %f\n", defaultReblendClassifyTol );
//   else
//      PRINT_INFO( "Healer Default Reblend Classify Tol = calculated\n" );
//
//   if( defaultReblendTol != CUBIT_DBL_MAX )
//      PRINT_INFO( "Healer Default Reblend Tol = %f\n", defaultReblendTol );
//   else
//      PRINT_INFO( "Healer Default Reblend Tol = calculated\n" );
}

CubitBoolean AcisHealerTool::is_simplify_initialized(ENTITY* ent)
{
   ATTRIB_HH_AGGR_SIMPLIFY * att = (ATTRIB_HH_AGGR_SIMPLIFY *)
      find_leaf_attrib(ent,ATTRIB_HH_AGGR_SIMPLIFY_TYPE);
   if (att!=NULL)
      return CUBIT_TRUE;
   else
      return CUBIT_FALSE;
}

CubitBoolean AcisHealerTool::is_stitch_initialized(ENTITY* ent)
{
   ATTRIB_HH_AGGR_STITCH * att = (ATTRIB_HH_AGGR_STITCH *)
      find_leaf_attrib(ent,ATTRIB_HH_AGGR_STITCH_TYPE);
   if (att!=NULL)
      return CUBIT_TRUE;
   else
      return CUBIT_FALSE;
}

CubitBoolean AcisHealerTool::is_geombuild_initialized(ENTITY* ent)
{
   ATTRIB_HH_AGGR_GEOMBUILD * att = (ATTRIB_HH_AGGR_GEOMBUILD *)
      find_leaf_attrib(ent,ATTRIB_HH_AGGR_GEOMBUILD_TYPE);
   if (att!=NULL)
      return CUBIT_TRUE;
   else
      return CUBIT_FALSE;
}

CubitBoolean AcisHealerTool::is_initialized(ENTITY* ent)
{
   if( hdebug )
   {
      if( is_simplify_initialized(ent) )
         PRINT_INFO("Simplify initialized\n");
      else
         PRINT_INFO("Not simplify initialized\n");
      if( is_stitch_initialized(ent) )
         PRINT_INFO("Stitch initialized\n");
      else
         PRINT_INFO("Not stitch initialized\n");
      if( is_geombuild_initialized(ent) )
         PRINT_INFO("Geombuild initialized\n");
      else
         PRINT_INFO("Not geombuild initialized\n");
   }

   if( is_simplify_initialized(ent) && is_stitch_initialized(ent) &&
      is_geombuild_initialized(ent) )
      return CUBIT_TRUE;
   else
      return CUBIT_FALSE;
}

CubitStatus AcisHealerTool::get_bad_vertices( BODY* BODY_ptr, DLIList<RefVertex*>& vertex_list)
{
   ENTITY_LIST ent_list;
   outcome result;
      
   result = api_hh_get_bad_vertices( BODY_ptr, ent_list );
   if( !result.ok() )
   {
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
      PRINT_ERROR( "unable to get bad vertices\n" );
      return CUBIT_FAILURE;
   }

   for (int i=0;i<ent_list.count();i++)
   {
      TopologyEntity* TE_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
      RefVertex* ref_vertex = CAST_TO( TE_ptr, RefVertex );
      if( ref_vertex == NULL )
      {
         PRINT_ERROR( "Couldn't find CUBIT vertex that corresponds to ACIS vertex!\n" );
         return CUBIT_FAILURE;
      }
      else
         vertex_list.append_unique(ref_vertex);
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::get_bad_edges( BODY* BODY_ptr, DLIList<RefEdge*>& edge_list)
{
   ENTITY_LIST ent_list;
   outcome result;
   
   result = api_hh_get_bad_edges( BODY_ptr, ent_list);
   if( !result.ok() )
   {
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
      PRINT_ERROR( "unable to get bad curves\n" );
      return CUBIT_FAILURE;
   }
   
   for (int i=0;i<ent_list.count();i++)
   {
      TopologyEntity* TE_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
      RefEdge* ref_edge = CAST_TO( TE_ptr, RefEdge );
      if( ref_edge == NULL )
      {
         PRINT_ERROR( "Couldn't find CUBIT curve that corresponds to ACIS edge!\n" );
         return CUBIT_FAILURE;
      }
      else
         edge_list.append_unique(ref_edge);
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::get_bad_coedges( BODY* BODY_ptr, DLIList<CoEdge*>& coedge_list)
{
   ENTITY_LIST ent_list;
   outcome result;
   
   result = api_hh_get_bad_coedges( BODY_ptr, ent_list);
   if( !result.ok() )
   {
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
      PRINT_ERROR( "unable to get bad coedges\n" );
      return CUBIT_FAILURE;
   }
   
   for (int i=0;i<ent_list.count();i++)
   {
      TopologyEntity* TE_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
      CoEdge* coedge = CAST_TO( TE_ptr, CoEdge );
      if( coedge == NULL )
      {
         PRINT_ERROR( "Couldn't find CUBIT coedge that corresponds to ACIS coedge!\n" );
         return CUBIT_FAILURE;
      }
      else
         coedge_list.append_unique(coedge);
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::get_bad_loops( BODY* BODY_ptr, DLIList<Loop*>& loop_list)
{
   ENTITY_LIST ent_list;
   outcome result;
   
   result = api_hh_get_bad_loops( BODY_ptr, ent_list );
   if( !result.ok() )
   {
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
      PRINT_ERROR( "unable to get bad loops\n" );
      return CUBIT_FAILURE;
   }
   
   for (int i=0;i<ent_list.count();i++)
   {
      TopologyEntity* TE_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
      Loop* loop = CAST_TO( TE_ptr, Loop );
      if( loop == NULL )
      {
         PRINT_ERROR( "Couldn't find CUBIT loop that corresponds to ACIS LOOP!\n" );
         return CUBIT_FAILURE;
      }
      else
         loop_list.append_unique(loop);
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::get_bad_faces( BODY* BODY_ptr, DLIList<RefFace*>& face_list)
{
   ENTITY_LIST ent_list;
   outcome result;
   
   result = api_hh_get_bad_faces( BODY_ptr, ent_list);
   if( !result.ok() )
   {
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
      PRINT_ERROR( "unable to get bad faces\n" );
      return CUBIT_FAILURE;
   }
   
   for (int i=0;i<ent_list.count();i++)
   {
      TopologyEntity* TE_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
      RefFace* ref_face = CAST_TO( TE_ptr, RefFace );
      if( ref_face == NULL )
      {
         PRINT_ERROR( "Couldn't find CUBIT surface that corresponds to ACIS face!\n" );
         return CUBIT_FAILURE;
      }
      else
         face_list.append_unique(ref_face);
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::get_bad_shells( BODY* BODY_ptr, DLIList<Shell*>& shell_list)
{
   ENTITY_LIST ent_list;
   outcome result;
   
   result = api_hh_get_bad_shells( BODY_ptr, ent_list);
   if( !result.ok() )
   {
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
      PRINT_ERROR( "unable to get bad shells\n" );
      return CUBIT_FAILURE;
   }
   
   for (int i=0;i<ent_list.count();i++)
   {
      TopologyEntity* TE_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
      Shell* shell = CAST_TO( TE_ptr, Shell );
      if( shell == NULL )
      {
         PRINT_ERROR( "Couldn't find CUBIT shell that corresponds to ACIS shell!\n" );
         return CUBIT_FAILURE;
      }
      else
         shell_list.append_unique(shell);
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::get_bad_volumes( BODY* BODY_ptr, DLIList<RefVolume*>& volume_list)
{
   ENTITY_LIST ent_list;
   outcome result;
   
   result = api_hh_get_bad_lumps( BODY_ptr, ent_list);
   if( !result.ok() )
   {
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error( result );
      PRINT_ERROR( "unable to get bad volumes\n" );
      return CUBIT_FAILURE;
   }
   
   for (int i=0;i<ent_list.count();i++)
   {
      TopologyEntity* TE_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
      RefVolume* volume = CAST_TO( TE_ptr, RefVolume );
      if( volume == NULL )
      {
         PRINT_ERROR( "Couldn't find CUBIT volume that corresponds to ACIS lump!\n" );
         return CUBIT_FAILURE;
      }
      else
         volume_list.append_unique(volume);
   }
   return CUBIT_SUCCESS;
}


// Gets coedges that do not lie on faces
CubitStatus AcisHealerTool::get_coedges_not_on_faces(BODY* body, DLIList<CoEdge*> coedge_list)
{
   CoEdge* coedge_ptr;
   ENTITY_LIST ent_list;
   api_get_coedges( body, ent_list);
   for (int i=0;i<ent_list.count();i++)
   {
      ATTRIB_HH_ENT_GEOMBUILD_COEDGE *att = (ATTRIB_HH_ENT_GEOMBUILD_COEDGE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_COEDGE_TYPE);
      if (NULL != att)
      {
         if(att->get_on_face()==1)
         {
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(((COEDGE *)ent_list[i]));
            coedge_ptr = CAST_TO( te_ptr, CoEdge );
            if( coedge_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find CoEdge from ACIS COEDGE!\n" );
               return CUBIT_FAILURE;
            }
            coedge_list.append_unique( coedge_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Gets coedges without partners
CubitStatus AcisHealerTool::get_coedges_no_partner(BODY *body, DLIList<CoEdge*> coedge_list)
{
   CoEdge* coedge_ptr;
   ENTITY_LIST ent_list;
   api_get_coedges(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      COEDGE* coedge = (COEDGE*)ent_list[i];
      if(coedge->partner() == NULL)
      {
         TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(((COEDGE *)ent_list[i]));
         coedge_ptr = CAST_TO( te_ptr, CoEdge );
         if( coedge_ptr == NULL )
         {
            PRINT_ERROR( "Unable to find CoEdge from ACIS COEDGE!\n" );
            return CUBIT_FAILURE;
         }
         coedge_list.append_unique( coedge_ptr );
      }
   }
   return CUBIT_SUCCESS;
}

// Get vertices that do not lie on associated faces
CubitStatus AcisHealerTool::get_vertices_not_on_faces(BODY *body, DLIList<RefVertex*> ref_vertex_list)			
{
   RefVertex* ref_vertex_ptr;
   ENTITY_LIST ent_list;
   api_get_vertices(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_VERTEX *att = (ATTRIB_HH_ENT_GEOMBUILD_VERTEX *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_VERTEX_TYPE);
      if (NULL != att){
         if(att->get_on_faces()==1){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_vertex_ptr = CAST_TO( te_ptr, RefVertex );
            if( ref_vertex_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefVertex from ACIS VERTEX!\n" );
               return CUBIT_FAILURE;
            }
            ref_vertex_list.append_unique( ref_vertex_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get vertices that do not lie on the associated edges
CubitStatus AcisHealerTool::get_vertices_not_on_edges(BODY *body, DLIList<RefVertex*> ref_vertex_list)
{
   RefVertex* ref_vertex_ptr;
   ENTITY_LIST ent_list;
   api_get_vertices(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_VERTEX *att = (ATTRIB_HH_ENT_GEOMBUILD_VERTEX *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_VERTEX_TYPE);
      if (NULL != att){
         if(att->get_on_edges()==1){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_vertex_ptr = CAST_TO( te_ptr, RefVertex );
            if( ref_vertex_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefVertex from ACIS VERTEX!\n" );
               return CUBIT_FAILURE;
            }
            ref_vertex_list.append_unique( ref_vertex_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get vertices where edges do not meet
CubitStatus AcisHealerTool::get_vertices_edges_dont_meet(BODY *body, DLIList<RefVertex*> ref_vertex_list)
{
   RefVertex* ref_vertex_ptr;
   ENTITY_LIST ent_list;
   api_get_vertices(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_VERTEX *att = (ATTRIB_HH_ENT_GEOMBUILD_VERTEX *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_VERTEX_TYPE);
      if (NULL != att){
         if(att->get_edges_meet()==0){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_vertex_ptr = CAST_TO( te_ptr, RefVertex );
            if( ref_vertex_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefVertex from ACIS VERTEX!\n" );
               return CUBIT_FAILURE;
            }
            ref_vertex_list.append_unique( ref_vertex_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get discontinuous curves
CubitStatus AcisHealerTool::get_discontinuous_curves(BODY *body, DLIList<RefEdge*> ref_edge_list)
{
   RefEdge* ref_edge_ptr;
   ENTITY_LIST ent_list;
   api_get_edges(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_EDGE *att = (ATTRIB_HH_ENT_GEOMBUILD_EDGE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_EDGE_TYPE);
      if (NULL != att){
         ATTRIB_HH_ENT_GEOMBUILD_CURVE *attc = (ATTRIB_HH_ENT_GEOMBUILD_CURVE *)
            find_leaf_attrib((ENTITY *)(att->new_geometry()),ATTRIB_HH_ENT_GEOMBUILD_CURVE_TYPE);
         if(attc && attc->get_continuity()==1){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_edge_ptr = CAST_TO( te_ptr, RefEdge );
            if( ref_edge_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefEdge from ACIS EDGE!\n" );
               return CUBIT_FAILURE;
            }
            ref_edge_list.append_unique( ref_edge_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get degenerate curves
CubitStatus AcisHealerTool::get_degenerate_curves(BODY *body, DLIList<RefEdge*> ref_edge_list)
{
   RefEdge* ref_edge_ptr;
   ENTITY_LIST ent_list;
   api_get_edges(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_EDGE *att = (ATTRIB_HH_ENT_GEOMBUILD_EDGE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_EDGE_TYPE);
      if (NULL != att){
         ATTRIB_HH_ENT_GEOMBUILD_CURVE *attc = (ATTRIB_HH_ENT_GEOMBUILD_CURVE *)
            find_leaf_attrib((ENTITY *)(att->new_geometry()),ATTRIB_HH_ENT_GEOMBUILD_CURVE_TYPE);
         if(attc && attc->get_degeneracy()==1){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_edge_ptr = CAST_TO( te_ptr, RefEdge );
            if( ref_edge_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefEdge from ACIS EDGE!\n" );
               return CUBIT_FAILURE;
            }
            ref_edge_list.append_unique( ref_edge_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get self intersecting curves
CubitStatus AcisHealerTool::get_self_intersecting_curves(BODY *body, DLIList<RefEdge*> ref_edge_list)
{
   RefEdge* ref_edge_ptr;
   ENTITY_LIST ent_list;
   api_get_edges(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_EDGE *att = (ATTRIB_HH_ENT_GEOMBUILD_EDGE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_EDGE_TYPE);
      if (NULL != att){
         ATTRIB_HH_ENT_GEOMBUILD_CURVE *attc = (ATTRIB_HH_ENT_GEOMBUILD_CURVE *)
            find_leaf_attrib((ENTITY *)(att->new_geometry()),ATTRIB_HH_ENT_GEOMBUILD_CURVE_TYPE);
         if(attc && attc->get_selfint()==1){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_edge_ptr = CAST_TO( te_ptr, RefEdge );
            if( ref_edge_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefEdge from ACIS EDGE!\n" );
               return CUBIT_FAILURE;
            }
            ref_edge_list.append_unique( ref_edge_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get periodic curves
CubitStatus AcisHealerTool::get_periodic_curves(BODY *body, DLIList<RefEdge*> ref_edge_list)
{
   RefEdge* ref_edge_ptr;
   ENTITY_LIST ent_list;
   api_get_edges(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_EDGE *att = (ATTRIB_HH_ENT_GEOMBUILD_EDGE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_EDGE_TYPE);
      if (NULL != att){
         ATTRIB_HH_ENT_GEOMBUILD_CURVE *attc = (ATTRIB_HH_ENT_GEOMBUILD_CURVE *)
            find_leaf_attrib((ENTITY *)(att->new_geometry()),ATTRIB_HH_ENT_GEOMBUILD_CURVE_TYPE);
         if(attc && attc->get_closure()==2){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_edge_ptr = CAST_TO( te_ptr, RefEdge );
            if( ref_edge_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefEdge from ACIS EDGE!\n" );
               return CUBIT_FAILURE;
            }
            ref_edge_list.append_unique( ref_edge_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get closed curves
CubitStatus AcisHealerTool::get_closed_curves(BODY *body, DLIList<RefEdge*> ref_edge_list)
{
   RefEdge* ref_edge_ptr;
   ENTITY_LIST ent_list;
   api_get_edges(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_EDGE *att = (ATTRIB_HH_ENT_GEOMBUILD_EDGE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_EDGE_TYPE);
      if (NULL != att){
         ATTRIB_HH_ENT_GEOMBUILD_CURVE *attc = (ATTRIB_HH_ENT_GEOMBUILD_CURVE *)
            find_leaf_attrib((ENTITY *)(att->new_geometry()),ATTRIB_HH_ENT_GEOMBUILD_CURVE_TYPE);
         if(attc && attc->get_closure()){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_edge_ptr = CAST_TO( te_ptr, RefEdge );
            if( ref_edge_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefEdge from ACIS EDGE!\n" );
               return CUBIT_FAILURE;
            }
            ref_edge_list.append_unique( ref_edge_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get short edges (length less than geometry tolerance)
CubitStatus AcisHealerTool::get_short_edges(BODY *body, DLIList<RefEdge*> ref_edge_list)
{
   RefEdge* ref_edge_ptr;
   ENTITY_LIST ent_list;
   api_get_edges(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_EDGE *att = (ATTRIB_HH_ENT_GEOMBUILD_EDGE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_EDGE_TYPE);
      if (NULL != att)
      {
          // MIN_EDGE_LENGTH was renamed in Acis 5.2, so that's why
          // we have the #ifdef's here.
         if(att->get_length() < BHL_MIN_EDGE_LENGTH )
         {
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_edge_ptr = CAST_TO( te_ptr, RefEdge );
            if( ref_edge_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefEdge from ACIS EDGE!\n" );
               return CUBIT_FAILURE;
            }
            ref_edge_list.append_unique( ref_edge_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get tangent edges (edge joins two tangent faces)
CubitStatus AcisHealerTool::get_tangent_edges(BODY *body, DLIList<RefEdge*> ref_edge_list)
{
   RefEdge* ref_edge_ptr;
   ENTITY_LIST ent_list;
   api_get_edges(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_EDGE *att = (ATTRIB_HH_ENT_GEOMBUILD_EDGE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_EDGE_TYPE);
      if (NULL != att){
         if(att->get_vexity () == TANGENT){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_edge_ptr = CAST_TO( te_ptr, RefEdge );
            if( ref_edge_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefEdge from ACIS EDGE!\n" );
               return CUBIT_FAILURE;
            }
            ref_edge_list.append_unique( ref_edge_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get convex edges 
CubitStatus AcisHealerTool::get_convex_edges(BODY *body, DLIList<RefEdge*> ref_edge_list)
{
   RefEdge* ref_edge_ptr;
   ENTITY_LIST ent_list;
   api_get_edges(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_EDGE *att = (ATTRIB_HH_ENT_GEOMBUILD_EDGE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_EDGE_TYPE);
      if (NULL != att){
         if(att->get_vexity () ==CONVEX){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_edge_ptr = CAST_TO( te_ptr, RefEdge );
            if( ref_edge_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefEdge from ACIS EDGE!\n" );
               return CUBIT_FAILURE;
            }
            ref_edge_list.append_unique( ref_edge_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get concave edges
CubitStatus AcisHealerTool::get_concave_edges(BODY *body, DLIList<RefEdge*> ref_edge_list)
{
   RefEdge* ref_edge_ptr;
   ENTITY_LIST ent_list;
   api_get_edges(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_EDGE *att = (ATTRIB_HH_ENT_GEOMBUILD_EDGE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_EDGE_TYPE);
      if (NULL != att){
         if(att->get_vexity () ==CONCAVE){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_edge_ptr = CAST_TO( te_ptr, RefEdge );
            if( ref_edge_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefEdge from ACIS EDGE!\n" );
               return CUBIT_FAILURE;
            }
            ref_edge_list.append_unique( ref_edge_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get loops that do not lie on their associated faces
CubitStatus AcisHealerTool::get_loops_not_on_faces(BODY *body, DLIList<Loop*> loop_list)
{
   Loop* loop_ptr;
   ENTITY_LIST ent_list;
   api_get_loops(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_LOOP *att = (ATTRIB_HH_ENT_GEOMBUILD_LOOP *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_LOOP_TYPE);
      if (NULL != att){
         if(att->get_on_face()==0){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            loop_ptr = CAST_TO( te_ptr, Loop );
            if( loop_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find Loop from ACIS LOOP!\n" );
               return CUBIT_FAILURE;
            }
            loop_list.append_unique( loop_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get loops with incorrect orientation
CubitStatus AcisHealerTool::get_loops_disoriented(BODY *body, DLIList<Loop*> loop_list)
{
   Loop* loop_ptr;
   ENTITY_LIST ent_list;
   api_get_loops(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_LOOP *att = (ATTRIB_HH_ENT_GEOMBUILD_LOOP *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_LOOP_TYPE);
      if (NULL != att){
         if(att->get_oriented()==1){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            loop_ptr = CAST_TO( te_ptr, Loop );
            if( loop_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find Loop from ACIS LOOP!\n" );
               return CUBIT_FAILURE;
            }
            loop_list.append_unique( loop_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get loops that have coedge gaps
CubitStatus AcisHealerTool::get_loops_gaps(BODY *body, DLIList<Loop*> loop_list)
{
   Loop* loop_ptr;
   ENTITY_LIST ent_list;
   api_get_loops(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_LOOP *att = (ATTRIB_HH_ENT_GEOMBUILD_LOOP *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_LOOP_TYPE);
      if (NULL != att){
         if(att->get_gaps()==1){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            loop_ptr = CAST_TO( te_ptr, Loop );
            if( loop_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find Loop from ACIS LOOP!\n" );
               return CUBIT_FAILURE;
            }
            loop_list.append_unique( loop_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get open loops
CubitStatus AcisHealerTool::get_loops_open (BODY *body, DLIList<Loop*> loop_list)
{
   Loop* loop_ptr;
   ENTITY_LIST ent_list;
   api_get_loops(body, ent_list);
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_LOOP *att = (ATTRIB_HH_ENT_GEOMBUILD_LOOP *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_LOOP_TYPE);
      if (NULL != att){
         if(att->get_closure()==1){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            loop_ptr = CAST_TO( te_ptr, Loop );
            if( loop_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find Loop from ACIS LOOP!\n" );
               return CUBIT_FAILURE;
            }
            loop_list.append_unique( loop_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get discontinuous surfaces
CubitStatus AcisHealerTool::get_discontinuous_surfaces(BODY *body, DLIList<RefFace*> ref_face_list)
{
   RefFace* ref_face_ptr;
   ENTITY_LIST ent_list;
   api_get_faces( body, ent_list );
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_FACE *att = (ATTRIB_HH_ENT_GEOMBUILD_FACE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_FACE_TYPE);
      if (NULL != att){
         ATTRIB_HH_ENT_GEOMBUILD_SURFACE *attc = (ATTRIB_HH_ENT_GEOMBUILD_SURFACE *)
            find_leaf_attrib((ENTITY *)(att->new_geometry()),ATTRIB_HH_ENT_GEOMBUILD_SURFACE_TYPE);
         if(attc && attc->get_continuity() != 0 || att->is_discontinuous()){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_face_ptr = CAST_TO( te_ptr, RefFace );
            if( ref_face_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefFace from ACIS FACE!\n" );
               return CUBIT_FAILURE;
            }
            ref_face_list.append_unique( ref_face_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get degenerate surfaces
CubitStatus AcisHealerTool::get_degenerate_surfaces(BODY *body, DLIList<RefFace*> ref_face_list)
{
   RefFace* ref_face_ptr;
   ENTITY_LIST ent_list;
   api_get_faces( body, ent_list );
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_FACE *att = (ATTRIB_HH_ENT_GEOMBUILD_FACE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_FACE_TYPE);
      if (NULL != att){
         ATTRIB_HH_ENT_GEOMBUILD_SURFACE *attc = (ATTRIB_HH_ENT_GEOMBUILD_SURFACE *)
            find_leaf_attrib((ENTITY *)(att->new_geometry()),ATTRIB_HH_ENT_GEOMBUILD_SURFACE_TYPE);
         if(attc && attc->get_degeneracy() != 0){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_face_ptr = CAST_TO( te_ptr, RefFace );
            if( ref_face_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefFace from ACIS FACE!\n" );
               return CUBIT_FAILURE;
            }
            ref_face_list.append_unique( ref_face_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get self intersecting surfaces
CubitStatus AcisHealerTool::get_self_intersecting_surfaces(BODY *body, DLIList<RefFace*> ref_face_list)
{
   RefFace* ref_face_ptr;
   ENTITY_LIST ent_list;
   api_get_faces( body, ent_list );
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_FACE *att = (ATTRIB_HH_ENT_GEOMBUILD_FACE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_FACE_TYPE);
      if (NULL != att){
         ATTRIB_HH_ENT_GEOMBUILD_SURFACE *attc = (ATTRIB_HH_ENT_GEOMBUILD_SURFACE *)
            find_leaf_attrib((ENTITY *)(att->new_geometry()),ATTRIB_HH_ENT_GEOMBUILD_SURFACE_TYPE);
         if(attc && attc->get_selfint()==1){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_face_ptr = CAST_TO( te_ptr, RefFace );
            if( ref_face_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefFace from ACIS FACE!\n" );
               return CUBIT_FAILURE;
            }
            ref_face_list.append_unique( ref_face_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get periodic surfaces
CubitStatus AcisHealerTool::get_periodic_surfaces(BODY *body, DLIList<RefFace*> ref_face_list)
{
   RefFace* ref_face_ptr;
   ENTITY_LIST ent_list;
   api_get_faces( body, ent_list );
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_FACE *att = (ATTRIB_HH_ENT_GEOMBUILD_FACE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_FACE_TYPE);
      if (NULL != att){
         ATTRIB_HH_ENT_GEOMBUILD_SURFACE *attc = (ATTRIB_HH_ENT_GEOMBUILD_SURFACE *)
            find_leaf_attrib((ENTITY *)(att->new_geometry()),ATTRIB_HH_ENT_GEOMBUILD_SURFACE_TYPE);
         if(attc && attc->get_closure()==2){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_face_ptr = CAST_TO( te_ptr, RefFace );
            if( ref_face_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefFace from ACIS FACE!\n" );
               return CUBIT_FAILURE;
            }
            ref_face_list.append_unique( ref_face_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

// Get closed surfaces
CubitStatus AcisHealerTool::get_closed_surfaces (BODY *body, DLIList<RefFace*> ref_face_list)
{
   RefFace* ref_face_ptr;
   ENTITY_LIST ent_list;
   api_get_faces( body, ent_list );
   for (int i=0;i<ent_list.count();i++){
      ATTRIB_HH_ENT_GEOMBUILD_FACE *att = (ATTRIB_HH_ENT_GEOMBUILD_FACE *)
         find_leaf_attrib(ent_list[i],ATTRIB_HH_ENT_GEOMBUILD_FACE_TYPE);
      if (NULL != att){
         ATTRIB_HH_ENT_GEOMBUILD_SURFACE *attc = (ATTRIB_HH_ENT_GEOMBUILD_SURFACE *)
            find_leaf_attrib((ENTITY *)(att->new_geometry()),ATTRIB_HH_ENT_GEOMBUILD_SURFACE_TYPE);
         if(attc && attc->get_closure()){
            TopologyEntity* te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
            ref_face_ptr = CAST_TO( te_ptr, RefFace );
            if( ref_face_ptr == NULL )
            {
               PRINT_ERROR( "Unable to find RefFace from ACIS FACE!\n" );
               return CUBIT_FAILURE;
            }
            ref_face_list.append_unique( ref_face_ptr );
         }
      }
   }
   return CUBIT_SUCCESS;
}

CubitStatus AcisHealerTool::heal_bodies( DLIList<Body*> &body_list, 
                                         DLIList<Body*> &new_body_list,
                                         DLIList<TopologyEntity*> &bad_geometry,
                                         CubitBoolean rebuild,
                                         CubitBoolean keep_old,
                                         CubitBoolean make_tolerant,
                                         FILE* logfile_ptr )
{
   DLIList<Body*> show_body_list;
   DLIList<Body*> modified_body_list;
   DLIList<Body*> untouched_body_list;

   CubitBoolean delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old || rebuild);
   
   // Now do the healing
   Body* body_ptr;
   BODY* BODY_ptr;
   int percent_before, percent_after;
   CubitStatus status = CUBIT_SUCCESS;
   int old_body_id;
   for( int i=0; i<body_list.size(); i++ )
   {
      if (cubit_intr)
      {
        PRINT_WARNING("Healing interrupted.\n");
        break;
      }
      
      body_ptr = body_list.get_and_step();
      old_body_id = body_ptr->id();

      //get the ids of the volumes in this body
      DLIList<RefVolume*> volumes_of_body;
      body_ptr->ref_volumes( volumes_of_body );
      DLIList<int> volume_ids;
      int kk;
      for(kk=volumes_of_body.size(); kk--;) 
        volume_ids.append( volumes_of_body.get_and_step()->id() );
      
      if( rebuild )
      {
          PRINT_INFO( "Rebuilding Copy of Volume %d", volume_ids.get_and_step()  );
          if( volume_ids.size() > 1 )
            for( kk=volume_ids.size()-1; kk--;)
            PRINT_INFO( ", %d", volume_ids.get_and_step() );
          PRINT_INFO( "...\n" ); 
      }
      else
      {
        if( delete_attribs )
          PRINT_INFO( "Healing Copy of Volume %d", volume_ids.get_and_step() );
        else
          PRINT_INFO( "Healing Volume %d", volume_ids.get_and_step() );

          if( volume_ids.size() > 1 )
            for( kk=volume_ids.size()-1; kk--;)
            PRINT_INFO( ", %d", volume_ids.get_and_step() );
          PRINT_INFO( "...\n" ); 
      }
      
      // Copy this BODY
      BODY_ptr =  BodyACIS::get_BODY_ptr( body_ptr );
      if( BODY_ptr == NULL )
      {
        PRINT_ERROR( "Couldn't find ACIS BODY from CUBIT body %d\n", body_ptr->id() );
        continue;
      }
      BODY* copied_BODY_ptr = AcisModifyEngine::instance()->copy_BODY(BODY_ptr, delete_attribs);
      
      if( logfile_ptr )
      {
         fprintf( logfile_ptr, "***************************************************************\n" );
         if( delete_attribs ) // Somewhat dangerous, but probably okay
            fprintf( logfile_ptr, "Body %d (created from Body %d):\n", 
                     RefEntityFactory::instance()->maximum_id( "body")+1,
            old_body_id );
         else
            fprintf( logfile_ptr, "Body %d:\n", body_ptr->id() );
         fprintf( logfile_ptr, "***************************************************************\n" );
      }
      
      // Heal the BODY
      BODY *new_BODY_ptr = NULL;
      if( rebuild )
      {
         new_BODY_ptr = auto_rebuild_BODY( copied_BODY_ptr, percent_before, 
                                           percent_after, make_tolerant,
                                           logfile_ptr );
         
         if( new_BODY_ptr == NULL )
         {
            PRINT_ERROR( "Error rebuilding Volume %d", volume_ids.get_and_step() );

            if( volume_ids.size() > 1 )
              for( kk=volume_ids.size()-1; kk--;)
              PRINT_INFO( ", %d", volume_ids.get_and_step()  );
            PRINT_INFO( "\n" ); 

            AcisModifyEngine::instance()->get_acis_query_engine()->delete_ACIS_BODY( copied_BODY_ptr );
            status = CUBIT_FAILURE;
            if( logfile_ptr )
               fprintf( logfile_ptr, " ERROR: new Volume(s) not created\n" );
            continue;
         }
         else
         {
            // Success
           if( percent_before == -1 )
             PRINT_INFO( "  Percentage good geometry before healing: not available\n" );
           else
             PRINT_INFO( "  Percentage good geometry before healing: %d%%\n", percent_before );
           if( percent_after == -1 )
             PRINT_INFO( "  Percentage good geometry after healing:  not available\n" );
           else
             PRINT_INFO( "  Percentage good geometry after healing:  %d%%\n", percent_after );

           PRINT_INFO( " Successfully rebuilt and healed Copy of Volume %d", volume_ids.get_and_step() );
           if( volume_ids.size() > 1 )
             for( kk=volume_ids.size()-1; kk--;)
             PRINT_INFO( ", %d", volume_ids.get_and_step() );
           PRINT_INFO( "!\n" ); 
            
           //The copied BODY is no longer needed, so delete it
           AcisModifyEngine::instance()->get_acis_query_engine()->delete_ACIS_BODY( copied_BODY_ptr );
         }
      }
      else
      {
         int num_splines_simplified = 0;
         if( heal_BODY( copied_BODY_ptr, percent_before, percent_after, num_splines_simplified,
                        make_tolerant, logfile_ptr ) == CUBIT_SUCCESS )
         {
           if( num_splines_simplified > 0 )
              PRINT_INFO( "  Simplified %d spline surfaces\n", num_splines_simplified );

           if( percent_before == -1 )
             PRINT_INFO( "  Percentage good geometry before healing: not available\n" );
           else
             PRINT_INFO( "  Percentage good geometry before healing: %d%%\n", percent_before );
           if( percent_after == -1 )
             PRINT_INFO( "  Percentage good geometry after healing:  not available\n" );
           else
             PRINT_INFO( "  Percentage good geometry after healing:  %d%%\n", percent_after );

           if( delete_attribs )
             PRINT_INFO( " Successfully healed Copy of Volume %d", volume_ids.get_and_step() );
           else
             PRINT_INFO( " Successfully healed Volume %d", volume_ids.get_and_step() );

           if( volume_ids.size() > 1 )
             for( kk=volume_ids.size()-1; kk--;)
             PRINT_INFO( ", %d", volume_ids.get_and_step() );
           PRINT_INFO( "!\n" ); 
         }
         else
         {
            PRINT_ERROR( "Error healing Volume %d", volume_ids.get_and_step() );
            if( volume_ids.size() > 1 )
              for( kk=volume_ids.size()-1; kk--;)
              PRINT_INFO( ", %d", volume_ids.get_and_step() );
            PRINT_INFO( "\n" ); 

            end_BODY_for_healing( copied_BODY_ptr );
				//The copied BODY is no longer needed, so delete it
            AcisModifyEngine::instance()->get_acis_query_engine()->delete_ACIS_BODY( copied_BODY_ptr, CUBIT_TRUE );
            status = CUBIT_FAILURE;
            if( logfile_ptr )
               fprintf( logfile_ptr, " ERROR: Volume(s) not healed (unable to auto-heal)\n" );
            continue;
         }
      }
      
      // Modify or make body in CUBIT
      Body* temp_body = NULL;
      if( rebuild )
         temp_body = AcisToolUtil::get_new_Body( body_ptr, BODY_ptr,
                                            new_BODY_ptr, keep_old );
      else
         temp_body = AcisToolUtil::get_new_Body( body_ptr, BODY_ptr,
                                            copied_BODY_ptr, keep_old );

      Body* new_body_ptr = CAST_TO(temp_body, Body);

      //get the ids of the volumes in this body

      if( new_body_ptr!=NULL && new_body_ptr!=body_ptr )
      {
        volumes_of_body.clean_out();
        new_body_ptr->ref_volumes( volumes_of_body );
        volume_ids.clean_out();
        for(kk=volumes_of_body.size(); kk--;) 
          volume_ids.append( volumes_of_body.get_and_step()->id() );

         PRINT_INFO( "Created new volume(s) %d", volume_ids.get_and_step() );
         if( volume_ids.size() > 1 )
           for( kk=volume_ids.size()-1; kk--;)
           PRINT_INFO( ", %d", volume_ids.get_and_step() );
         PRINT_INFO( "\n" ); 

         new_body_list.append( new_body_ptr );
      }
      else if( new_body_ptr!=NULL && new_body_ptr==body_ptr )
      {
         PRINT_INFO( "Modified volume %d", volume_ids.get_and_step() );
         if( volume_ids.size() > 1 )
           for( kk=volume_ids.size()-1; kk--;)
           PRINT_INFO( ", %d", volume_ids.get_and_step() );
         PRINT_INFO( "\n" ); 
      }

      // Figure-out which bodies we want to show for.  If the keep_old option
      // is on, we show for that body too.  The only caveat is that if
      // the old body was never analyzed, the user will get a warning
      // that it is not analyzed.  This is deemed better than forcing it
      // to be analyzed as that would take extra time.
      if( keep_old )
        show_body_list.append( body_ptr );
      if( new_body_ptr != NULL )
        show_body_list.append_unique( new_body_ptr );
      
      // Print out warnings, if settings call for it
      if( checkCurveOnHeal )
      {
        DLIList<RefEdge*> short_curve_list;
        DLIList<RefEntity*> ref_entities;
        new_body_ptr->ref_edges( short_curve_list );
        CAST_LIST_TO_PARENT(short_curve_list, ref_entities);
        measure_filter(ref_entities, curveLimit, 1e-14 );
        CAST_LIST(ref_entities, short_curve_list, RefEdge);
        if( short_curve_list.size() )
        {
          //determine the volumes containing short curves
          DLIList<RefVolume*> vols_with_short_curves, tmp_vols;
          int kk;
          for( kk=short_curve_list.size(); kk--;)
          {
            short_curve_list.get_and_step()->ref_volumes( tmp_vols );
            vols_with_short_curves += tmp_vols;
          }
          vols_with_short_curves.uniquify_unordered();

          char pre[100];
          sprintf( pre, "WARNING: Volume %d ", vols_with_short_curves.get_and_step()->id() );
          if( vols_with_short_curves.size() > 1 )
            for( kk=vols_with_short_curves.size()-1; kk--;)
              sprintf( pre, ", %d" );
          sprintf( pre, " contains %d short curves: ", short_curve_list.size() ); 

          DLIList<CubitEntity*> temp_list;
          CAST_LIST_TO_PARENT(short_curve_list, temp_list);    
          CubitUtil::list_entity_ids( pre, temp_list );
        }
      }
      if( checkSurfOnHeal )
      {
        DLIList<RefFace*> small_surf_list;
        new_body_ptr->ref_faces( small_surf_list );
        DLIList<RefEntity*> ref_entities;
        CAST_LIST_TO_PARENT(small_surf_list, ref_entities);
        measure_filter(ref_entities, surfLimit, 1e-14 );
        CAST_LIST(ref_entities, small_surf_list, RefFace);
        if( small_surf_list.size() )
        {
          //determine the volumes containing short curves
          DLIList<RefVolume*> vols_with_small_surfs, tmp_vols;
          int kk;
          for( kk=small_surf_list.size(); kk--;)
          {
            small_surf_list.get_and_step()->ref_volumes( tmp_vols );
            vols_with_small_surfs += tmp_vols;
          }
          vols_with_small_surfs.uniquify_unordered();

          char pre[100];
          sprintf( pre, "WARNING: Volume %d ", vols_with_small_surfs.get_and_step() );
          if( vols_with_small_surfs.size() > 1 )
            for( kk=vols_with_small_surfs.size()-1; kk--;)
              sprintf( pre, ", %d" );
          sprintf( pre, " contains %d small surfs: ", small_surf_list.size() ); 
          
          DLIList<CubitEntity*> temp_list;
          CAST_LIST_TO_PARENT( small_surf_list , temp_list);    
          CubitUtil::list_entity_ids( pre, temp_list);
        }
      }
      if( checkVolOnHeal )
      {
        // Warn the user if the volume is small
        if( new_body_ptr != NULL && new_body_ptr->measure() <= volLimit )
          PRINT_WARNING( "Volume %d ", volume_ids.get_and_step() );
        if( volume_ids.size() > 1 )
          for( kk=volume_ids.size()-1; kk--;)
            PRINT_INFO( ", %d", volume_ids.get_and_step() );
          PRINT_INFO(" has zero volume\n");
      }
   }

   if( showOnHeal )
     get_badgeom( CUBIT_FALSE, show_body_list, CUBIT_TRUE, logfile_ptr, bad_geometry );

   if( cleanAtt )
     clean_attributes( show_body_list );

   return status;
}

void AcisHealerTool::list_tolerances( DLIList<Body*> &body_list/*, CubitBoolean analyze*/ )
{  
   DLIList<Body*> bodies_not_analyzed;

   PRINT_INFO( "                       Stitch               Geombuild\n" );
   PRINT_INFO( "__Body__ Simplify   Min      Max   Geombuild Analytic Isospline\n" );
   int cnt = 0;
   for( int i=0; i<body_list.size(); i++ )
   {
      Body* body_ptr = body_list.get_and_step();
      BODY* BODY_ptr =  BodyACIS::get_BODY_ptr( body_ptr );

      if( BODY_ptr == NULL )
      {
        PRINT_ERROR( "Couldn't find ACIS BODY from CUBIT Body %d\n", body_ptr->id() );
        continue;
      }

      if( is_analyzed( BODY_ptr ) == CUBIT_FALSE )
      {
         bodies_not_analyzed.append( body_ptr );
         continue;
      }

      cnt++;

      PRINT_INFO( "%6d   %7.5f  %7.5f  %7.5f  %7.5f  %7.5f  %7.5f\n",
         body_ptr->id(),
         get_simplify_tol( BODY_ptr ),
         get_stitch_min_tol( BODY_ptr ),
         get_stitch_max_tol( BODY_ptr ),
         get_geombuild_tol( BODY_ptr ),
         get_analytic_tol( BODY_ptr ),
         get_isospline_tol( BODY_ptr ) );
   }

   if( cnt == 0 )
     PRINT_WARNING(" No bodies have healer tolerance settings.\n"
              " Hint: try 'healer analyze' with 'CleanAtt' option set to 'Off'\n" );
   else
   {
     PRINT_INFO( "\n" );
     DLIList<CubitEntity*> temp_list;
     CAST_LIST_TO_PARENT( bodies_not_analyzed , temp_list);    
     CubitUtil::list_entity_ids( "Bodies without tolerance settings: ", temp_list);
   }
}

CubitBoolean AcisHealerTool::is_analyzed( BODY* BODY_ptr )
{
//   if( is_geombuild_initialized( BODY_ptr ) == CUBIT_FALSE )
//      return CUBIT_FALSE;

   // This is kind-of a hack, but seems to work
//   ENTITY_LIST ent_list;
//   outcome result = api_hh_get_bad_shells( BODY_ptr, ent_list );
//   if( !result.ok() )
//      return CUBIT_FALSE;

  ATTRIB_HH_AGGR_GEOMBUILD* att = find_aggr_geombuild( BODY_ptr );
  if( att == NULL )
    return CUBIT_FALSE;

   return CUBIT_TRUE;
}

void AcisHealerTool::clean_attributes( DLIList<Body*>& body_list )
{
  for( int i=0; i<body_list.size(); i++ )
  {
    Body* body_ptr = body_list.get_and_step();
    BODY* BODY_ptr =  BodyACIS::get_BODY_ptr( body_ptr );

    if( BODY_ptr == NULL )
    {
      PRINT_ERROR( "Couldn't find ACIS BODY from CUBIT Body %d\n", 
                    body_ptr->id() );
      continue;
    }
    
    if( end_BODY_for_healing( BODY_ptr ) == CUBIT_FAILURE )
    {
      PRINT_ERROR( "Unable to cleanup healer attributes on Body %d\n", 
                   body_ptr->id() );
      continue;
    }
  }
}

void AcisHealerTool::get_curves_from_coedges( DLIList<CoEdge*> &coedge_list, 
                                             DLIList<RefEdge*> &curve_list )
{
  RefEdge* ref_edge_ptr;
  CoEdge* coedge_ptr;
  for( int i=0; i<coedge_list.size(); i++ )
  {
    coedge_ptr = coedge_list.get_and_step();
    RefEdge *temp_edge = coedge_ptr->get_ref_edge_ptr();
    ref_edge_ptr = CAST_TO(temp_edge, RefEdge);
    curve_list.append_unique( ref_edge_ptr );
  }
}
   
void AcisHealerTool::get_curves_from_loops( DLIList<Loop*> &loop_list, 
                                           DLIList<RefEdge*> &curve_list )
{
  Loop* loop_ptr;
  for( int i=0; i<loop_list.size(); i++ )
  {
    loop_ptr = loop_list.get_and_step();
    DLIList<RefEdge*> tmp_curve_list;
    loop_ptr->ordered_ref_edges( tmp_curve_list );
    curve_list.merge_unique( tmp_curve_list );
  }
}

void AcisHealerTool::get_surfaces_from_shells( DLIList<Shell*> &shell_list, 
                                               DLIList<RefFace*> &surface_list )
{
  Shell* shell_ptr;
  for( int i=0; i<shell_list.size(); i++ )
  {
    shell_ptr = shell_list.get_and_step();
    DLIList<RefFace*> tmp_surface_list;
    shell_ptr->ref_faces( tmp_surface_list );
    surface_list.merge_unique( tmp_surface_list );
  }
}

void AcisHealerTool::print_none( FILE* file_ptr, const char* str )
{
  if( file_ptr )
    fprintf( file_ptr, "%s", str );
  PRINT_INFO( "%s", str );
}

void AcisHealerTool::print_vertices( FILE* logfile_ptr, DLIList<RefVertex*> &vertex_list,  
                                    const char* str_none, const char* pre_str )
{
  if( vertex_list.size() == 0 )
  {
    print_none( logfile_ptr, str_none );
    return;
  }
  char pre[100];
  sprintf( pre, pre_str, vertex_list.size() );
  DLIList<CubitEntity*> temp_list;
  CAST_LIST_TO_PARENT( vertex_list , temp_list);    
  CubitUtil::list_entity_ids( pre, temp_list);
}

void AcisHealerTool::print_curves( FILE* logfile_ptr, DLIList<RefEdge*> &curve_list,  
                                    const char* str_none, const char* pre_str )
{
  if( curve_list.size() == 0 )
  {
    print_none( logfile_ptr, str_none );
    return;
  }
  char pre[100];
  sprintf( pre, pre_str, curve_list.size() );
  DLIList<CubitEntity*> temp_list;
  CAST_LIST_TO_PARENT( curve_list , temp_list);    
  CubitUtil::list_entity_ids( pre, temp_list);
}

void AcisHealerTool::print_coedges( FILE* logfile_ptr, DLIList<CoEdge*> &coedge_list,  
                                    const char* str_none, const char* pre_str )
{
  if( coedge_list.size() == 0 )
  {
    print_none( logfile_ptr, str_none );
    return;
  }
  char pre[100];
  sprintf( pre, pre_str, coedge_list.size() );
  DLIList<RefEdge*> curve_list;
  get_curves_from_coedges( coedge_list, curve_list );
  DLIList<CubitEntity*> temp_list;
  CAST_LIST_TO_PARENT( curve_list , temp_list);    
  CubitUtil::list_entity_ids( pre, temp_list);
}

void AcisHealerTool::print_loops( FILE* logfile_ptr, DLIList<Loop*> &loop_list,  
                                  const char* str_none, const char* pre_str )
{
  if( loop_list.size() == 0 )
  {
    print_none( logfile_ptr, str_none );
    return;
  }
  char pre[100];
  sprintf( pre, pre_str, loop_list.size() );
  DLIList<RefEdge*> curve_list;
  get_curves_from_loops( loop_list, curve_list );
  DLIList<CubitEntity*> temp_list;
  CAST_LIST_TO_PARENT( curve_list , temp_list);    
  CubitUtil::list_entity_ids( pre, temp_list);
}

void AcisHealerTool::print_surfaces( FILE* logfile_ptr, DLIList<RefFace*> &surface_list,  
                                     const char* str_none, const char* pre_str )
{
  if( surface_list.size() == 0 )
  {
    print_none( logfile_ptr, str_none );
    return;
  }
  char pre[100];
  sprintf( pre, pre_str, surface_list.size() );
  DLIList<CubitEntity*> temp_list;
  CAST_LIST_TO_PARENT( surface_list , temp_list);    
  CubitUtil::list_entity_ids( pre, temp_list);
}

void AcisHealerTool::print_shells( FILE* logfile_ptr, DLIList<Shell*> &shell_list,  
                                   const char* str_none, const char* pre_str )
{
  if( shell_list.size() == 0 )
  {
    print_none( logfile_ptr, str_none );
    return;
  }
  char pre[100];
  sprintf( pre, pre_str, shell_list.size() );
  DLIList<RefFace*> surface_list;
  get_surfaces_from_shells( shell_list, surface_list );
  DLIList<CubitEntity*> temp_list;
  CAST_LIST_TO_PARENT( surface_list , temp_list);    
  CubitUtil::list_entity_ids( pre, temp_list);
}

void AcisHealerTool::print_volumes( FILE* logfile_ptr, DLIList<RefVolume*> &volume_list,  
                                    const char* str_none, const char* pre_str )
{
  if( volume_list.size() == 0 )
  {
    print_none( logfile_ptr, str_none );
    return;
  }
  char pre[100];
  sprintf( pre, pre_str, volume_list.size() );
  DLIList<CubitEntity*> temp_list;
  CAST_LIST_TO_PARENT( volume_list , temp_list);    
  CubitUtil::list_entity_ids( pre, temp_list);
}

void AcisHealerTool::print_bodies( FILE* logfile_ptr, DLIList<Body*> &body_list,  
                                   const char* str_none, const char* pre_str )
{
  if( body_list.size() == 0 )
  {
    print_none( logfile_ptr, str_none );
    return;
  }
  char pre[100];
  sprintf( pre, pre_str, body_list.size() );
  DLIList<CubitEntity*> temp_list;
  CAST_LIST_TO_PARENT( body_list , temp_list);    
  CubitUtil::list_entity_ids( pre, temp_list);
}

CubitStatus AcisHealerTool::get_badgeom( CubitBoolean analyze_flag,
                                          DLIList<Body*> &body_list,
                                          CubitBoolean after_heal,
                                          FILE* logfile_ptr,
                                          DLIList<TopologyEntity*> &bad_geometry )
{
   if( body_list.size() == 0 )
     return CUBIT_SUCCESS;

   Body* body_ptr;
   BODY* BODY_ptr;

   DLIList<Body*> summary_body_list = body_list;

   int i, j;

   int percent_goodgeom;
   DLIList<TopologyEntity*> temp_top_list;

   // Summary output
   DLIList<RefVertex*> vertex_list, bad_vertex_list;
   DLIList<RefEdge*> curve_list, bad_curve_list, bad_curve_from_coedge_list, bad_curve_from_loop_list;
   DLIList<CoEdge*> coedge_list, bad_coedge_list;
   DLIList<Loop*> loop_list, bad_loop_list;
   DLIList<RefFace*> surface_list, bad_surface_list, bad_surface_from_shell_list;
   DLIList<Shell*> shell_list, bad_shell_list;
   DLIList<RefVolume*> volume_list, bad_volume_list;
   DLIList<Body*> bad_body_list;

   // Advanced output
   DLIList<CoEdge*> coedges_not_on_faces;
   DLIList<RefEdge*> curve_from_coedge_not_on_faces;
   DLIList<CoEdge*> coedges_no_partner;
   DLIList<RefEdge*> curve_from_coedges_no_partner;
   DLIList<RefVertex*> vertices_not_on_faces;
   DLIList<RefVertex*> vertices_not_on_edges;
   DLIList<RefVertex*> vertices_edges_dont_meet;
   DLIList<RefEdge*> discontinuous_curves;
   DLIList<RefEdge*> degenerate_curves;
   DLIList<RefEdge*> self_intersecting_curves;
   DLIList<RefEdge*> periodic_curves;
   DLIList<RefEdge*> closed_curves;
   DLIList<RefEdge*> short_edges;
   DLIList<RefEdge*> tangent_edges;
   DLIList<RefEdge*> convex_edges;
   DLIList<RefEdge*> concave_edges;
   DLIList<Loop*> loops_not_on_faces;
   DLIList<RefEdge*> curve_from_loops_not_on_faces;
   DLIList<Loop*> loops_disoriented;
   DLIList<RefEdge*> curve_from_loops_disoriented;
   DLIList<Loop*> loops_gaps;
   DLIList<RefEdge*> curve_from_loops_gaps;
   DLIList<Loop*> loops_open;
   DLIList<RefEdge*> curve_from_loops_open;
   DLIList<RefFace*> discontinuous_surfaces;
   DLIList<RefFace*> degenerate_surfaces;
   DLIList<RefFace*> self_intersecting_surfaces;
   DLIList<RefFace*> periodic_surfaces;
   DLIList<RefFace*> closed_surfaces;

   int show = 1;
/* CDE
   // ACIS bugs prevent this right now.  SRS 4-28-99
   if( showBadLoops || showBadSurfaces || showBadShells || showBadVolumes )
   {
     PRINT_WARNING( "currently you can't analyze bad loops, surfaces, shells and volumes.\n"
       "         This will be implemented in a future release.\n" );
   }
   */

#if 0
   if( showBadShells || showBadVolumes )
   {
     PRINT_WARNING( "currently you can't analyze bad shells and volumes.\n"
       "         This will be implemented in a future release.\n" );
   }
#endif

   if( !body_list.size() )
      return CUBIT_SUCCESS;

   PRINT_INFO("\n" );
   if (showMethod == 0)
     show = 0;

   if( logfile_ptr ) 
   {
      hh_set_bhl_log_file (logfile_ptr);
      CubitUtil::set_file_ptr(logfile_ptr);
   }
   
   for( i=0; i<body_list.size(); i++ )
   {
      if (cubit_intr)
      {
        PRINT_WARNING("Analysis interrupted.\n");
        break;
      }
      
      body_ptr = body_list.get_and_step();
      BODY_ptr =  BodyACIS::get_BODY_ptr(body_ptr);

      if( BODY_ptr == NULL )
      {
        PRINT_ERROR( "Couldn't find ACIS BODY from CUBIT Body %d\n", body_ptr->id() );
        continue;
      }

      //get the ids of the volumes in this body
      DLIList<RefVolume*> volumes_of_body;
      body_ptr->ref_volumes( volumes_of_body );
      DLIList<int> volume_ids;
      int kk;
      for(kk=volumes_of_body.size(); kk--;) 
        volume_ids.append( volumes_of_body.get_and_step()->id() );

      // Dump detailed information to the logfile if required
      if( logfile_ptr && !after_heal )
      {
         fprintf( logfile_ptr, "***************************************************************\n" );
         fprintf( logfile_ptr, "Body %d:\n", body_ptr->id() );
         fprintf( logfile_ptr, "Body %d containes volume %d", volume_ids.get_and_step() );
         if( volume_ids.size() > 1 )
           fprintf( logfile_ptr, ", %d", volume_ids.get_and_step() );
         fprintf( logfile_ptr, "\n");
         fprintf( logfile_ptr, "***************************************************************\n" );
      }

      if( analyze_flag )
      {
        // Force an analysis
        if( init_BODY_for_healing( BODY_ptr ) == CUBIT_FAILURE )
        {
          PRINT_ERROR( "Unable to initialize Volume %d", volume_ids.get_and_step() );
          if( volume_ids.size() > 1 )
            for( kk=volume_ids.size()-1; kk--; )
            PRINT_INFO( ", %d", volume_ids.get_and_step() );
          PRINT_INFO(" for healing\n");
          continue;
        }
        if( analyze_BODY_for_healing( BODY_ptr ) == CUBIT_FAILURE )
        {
          PRINT_ERROR( "Unable to analyze Volume %d", volume_ids.get_and_step() );
          if( volume_ids.size() > 1 )
            for( kk=volume_ids.size()-1; kk--; )
            PRINT_INFO( ", %d", volume_ids.get_and_step() );
          PRINT_INFO(" for healing\n");

          end_BODY_for_healing( BODY_ptr );
          continue;
        }
      }
      else if( is_analyzed( BODY_ptr ) == CUBIT_FALSE )
      {
        PRINT_WARNING( "Volume %d", volume_ids.get_and_step() );
        if( volume_ids.size() > 1 )
          for( kk=volume_ids.size()-1; kk--; )
            PRINT_INFO( ", %d", volume_ids.get_and_step() );
        PRINT_INFO(" not analyzed - cannot show the bad geometry\n"); 

        if( logfile_ptr )
          fprintf( logfile_ptr, " WARNING: body %d is not analyzed\n", body_ptr->id() );
        summary_body_list.remove( body_ptr );
        continue;
      }
      //else
      //{
        // We no longer write to a logfile when just showing bad geometry.  You
        // have to do the whole analysis for that.
        //
        // Explicitly write to the logfile.  Normally this is done during 
        // analysis.  Here we have to do it ourselves since we are not
        // actually analyzing the geometry.
        //if( logfile_ptr )
          //analysis_printout( BODY_ptr, logfile_ptr );
      //}
   
      // Grab percentage of good geometry and print it out
      percentage_goodgeom_after( BODY_ptr, percent_goodgeom ); // 5/3/99 - changed from before to after
      if( percent_goodgeom == -1 )
      {
        PRINT_INFO("Percentage good geometry in Volume %d", volume_ids.get_and_step() ); 
        if( volume_ids.size() > 1 )
          for( kk=volume_ids.size()-1; kk--; )
            PRINT_INFO( ", %d", volume_ids.get_and_step() );
        PRINT_INFO(" : not available\n"); 
      }
      else
      {
        PRINT_INFO("Percentage good geometry in Volume %d", volume_ids.get_and_step() ); 
        if( volume_ids.size() > 1 )
          for( kk=volume_ids.size()-1; kk--; )
            PRINT_INFO( ", %d", volume_ids.get_and_step() );
        PRINT_INFO(": %d%%\n", percent_goodgeom); 
      }

      // Keep track of which bodies have problems
      if( percent_goodgeom < 100 )
         bad_body_list.append( body_ptr );

      // Vertices
      if( get_bad_vertices( BODY_ptr, vertex_list) == CUBIT_FAILURE )
      {
        PRINT_WARNING( "Unable to retrieve bad vertices on volume %d",
          volume_ids.get_and_step() ); 
        if( volume_ids.size() > 1 )
          for( kk=volume_ids.size()-1; kk--; )
            PRINT_INFO( ", %d", volume_ids.get_and_step() );
        PRINT_INFO("\n"); 
      }
      else
      {
        if( showBadVertices && show )
          for( j=vertex_list.size(); j--; ) bad_geometry.append( vertex_list.get_and_step() );
        bad_vertex_list.merge_unique( vertex_list );
        vertex_list.clean_out();
      }

      // Curves
      if( get_bad_edges( BODY_ptr, curve_list) == CUBIT_FAILURE )
      {
        PRINT_WARNING( "Unable to retrieve bad curves on volume %d",
          volume_ids.get_and_step() ); 
        if( volume_ids.size() > 1 )
          for( kk=volume_ids.size()-1; kk--; )
           PRINT_INFO( ", %d", volume_ids.get_and_step() );
        PRINT_INFO("\n"); 
      }
      else
      {
        if( showBadCurves && show )
          for( j=curve_list.size(); j--; ) bad_geometry.append( curve_list.get_and_step() );
        bad_curve_list.merge_unique( curve_list );
        curve_list.clean_out();
      }
      
      // CoEdges
      if( get_bad_coedges( BODY_ptr, coedge_list) == CUBIT_FAILURE )
      {
        PRINT_WARNING( "Unable to retrieve bad coedges on volume %d",
          volume_ids.get_and_step() ); 
        if( volume_ids.size() > 1 )
          for( kk=volume_ids.size()-1; kk--; )
           PRINT_INFO( ", %d", volume_ids.get_and_step() );
        PRINT_INFO("\n"); 
      }
      else
      {
        if( showBadCoEdges && show )
          for( j=coedge_list.size(); j--; ) bad_geometry.append( coedge_list.get_and_step() );
        bad_coedge_list.merge_unique( coedge_list );
        coedge_list.clean_out();
      }

#if 0   
      // Loops
      if( get_bad_loops( BODY_ptr, loop_list) == CUBIT_FAILURE )
      {
        PRINT_WARNING( "Unable to retrieve bad loops on body %d\n",
          body_ptr->id() );
      }
      else
      {
        if( showBadLoops && show )
          for( j=loop_list.size(); j--; ) bad_geometry.append( loop_list.get_and_step() );
        bad_loop_list.merge_unique( loop_list );
        loop_list.clean_out();
      }

      //Surfaces
      if( get_bad_faces( BODY_ptr, surface_list) == CUBIT_FAILURE )
      {
        PRINT_WARNING( "Unable to retrieve bad surfaces on body %d\n",
          body_ptr->id() );
      }
      else
      {
        if( showBadSurfaces && show )
          for( j=surface_list.size(); j--; ) bad_geometry.append( surface_list.get_and_step() );
        bad_surface_list.merge_unique( surface_list );
        surface_list.clean_out();
      }
#endif
      
#if 0
      // Shells
      if( get_bad_shells( BODY_ptr, shell_list) == CUBIT_FAILURE )
      {
        PRINT_WARNING( "Unable to retrieve bad shells on body %d\n",
          body_ptr->id() );
      }
      else
      {
        if( showBadShells && show )
          for( j=shell_list.size(); j--; ) bad_geometry.append( shell_list.get_and_step() );
        bad_shell_list.merge_unique( shell_list );
        shell_list.clean_out();
      }
      
      // Volumes
      if( get_bad_volumes( BODY_ptr, volume_list) == CUBIT_FAILURE )
      {
        PRINT_WARNING( "Unable to retrieve bad volumes on body %d\n",
          body_ptr->id() );
      }
      else
      {
        if( showBadVolumes && show )
          for( j=volume_list.size(); j--; ) bad_geometry.append( volume_list.get_and_step() );
        bad_volume_list.merge_unique( volume_list );
        volume_list.clean_out();
      }
#endif

      // Bodies
      if( percent_goodgeom < 100 )
      {
        if( showBadBodies && show )
          bad_geometry.append( body_ptr );
      }

      if( showDetails == CUBIT_TRUE )
      {
        if( get_coedges_not_on_faces(BODY_ptr, coedges_not_on_faces) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get coedges not on surfaces\n" );
        if( get_coedges_no_partner(BODY_ptr, coedges_no_partner) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get coedges no partner\n" );
        if( get_vertices_not_on_faces(BODY_ptr, vertices_not_on_faces) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get vertices not on surfaces\n" );
        if( get_vertices_not_on_edges(BODY_ptr, vertices_not_on_edges) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get vertices not on curves\n" );
        if( get_vertices_edges_dont_meet(BODY_ptr, vertices_edges_dont_meet) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get vertices curves don't meet\n" );
        if( get_discontinuous_curves(BODY_ptr, discontinuous_curves) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get discontinuous curves\n" );
        if( get_degenerate_curves(BODY_ptr, degenerate_curves) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get degenerate curves\n" );
        if( get_self_intersecting_curves(BODY_ptr, self_intersecting_curves) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get self-interesecting curves\n" );
        if( get_periodic_curves(BODY_ptr, periodic_curves) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get periodic curves\n" );
        if( get_closed_curves(BODY_ptr, closed_curves) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get closed_curves\n" );
        if( get_short_edges(BODY_ptr, short_edges) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get short curves\n" );
        if( get_tangent_edges(BODY_ptr, tangent_edges) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get tangent curves\n" );
        if( get_convex_edges(BODY_ptr, convex_edges) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get convex curves\n" );
        if( get_concave_edges(BODY_ptr, concave_edges) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get concave curves\n" );
        if( get_loops_not_on_faces(BODY_ptr, loops_not_on_faces) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get loops not on surfaces\n" );
        if( get_loops_disoriented(BODY_ptr, loops_disoriented) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get loops disoriented\n" );
        if( get_loops_gaps(BODY_ptr, loops_gaps) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get loops with gaps\n" );
        if( get_loops_open (BODY_ptr, loops_open) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get open loops\n" );
        if( get_discontinuous_surfaces(BODY_ptr, discontinuous_surfaces) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get discontinous surfaces\n" );
        if( get_degenerate_surfaces(BODY_ptr, degenerate_surfaces) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get degenerate surfaces\n" );
        if( get_self_intersecting_surfaces(BODY_ptr, self_intersecting_surfaces) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get self-intersecting surfaces\n" );
        if( get_periodic_surfaces(BODY_ptr, periodic_surfaces) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get periodic surfaces\n" );
        if( get_closed_surfaces (BODY_ptr, closed_surfaces) == CUBIT_FAILURE )
          PRINT_WARNING( "Unable to get closed surfaces\n" );
      }
   }

   if( showSummary )
   {
     char pre[100];
     
     print_none( logfile_ptr, "\nHEALER ANALYSIS SUMMARY:\n" );
     print_none( logfile_ptr, "------------------------\n" );

     if( summary_body_list.size() == 0 )
       print_none( logfile_ptr, "No volumes were analyzed for the summary\n" );
     else
     {
       //get all the volumes ids in the summary_body_list
       DLIList<RefVolume*>tmp_volumes, summary_volumes;
       int kk;
       for(kk=summary_body_list.size(); kk--;) 
       {
         tmp_volumes.clean_out();
         summary_body_list.get_and_step()->ref_volumes( tmp_volumes );
         summary_volumes += tmp_volumes;
       }

       DLIList<int> volume_ids;
       for(kk=summary_volumes.size(); kk--;) 
         volume_ids.append( summary_volumes.get_and_step()->id() );
       
       if( volume_ids.size() == 1 )
         strcpy( pre, "Analyzed 1 Volume: " );
       else
         sprintf( pre, "Analyzed %d Volumes: ", volume_ids.size() );
       DLIList<CubitEntity*> temp_list;
       CAST_LIST_TO_PARENT( summary_volumes , temp_list);    
       CubitUtil::list_entity_ids( pre, temp_list );
       
       print_vertices( logfile_ptr, bad_vertex_list, 
         "Found 0 bad Vertices.\n", "Found %d bad Vertices: " );
       
       print_curves( logfile_ptr, bad_curve_list, 
         "Found 0 bad Curves.\n", "Found %d bad Curves: " );
       
       print_coedges( logfile_ptr, bad_coedge_list, 
         "Found 0 bad CoEdges.\n",
         "Found %d bad CoEdges. The Curves are: " );
       
#if 0 
       print_loops( logfile_ptr, bad_loop_list, 
         "Found 0 bad Loops.\n",
         "Found %d bad Loops. The Curves are: " );
       
       print_surfaces( logfile_ptr, bad_surface_list, 
         "Found 0 bad Surfaces.\n",
         "Found %d bad Surfaces: " );
#endif

#if 0
       // ACIS bugs prevent us from using this for now
       print_shells( logfile_ptr, bad_shell_list, 
         "Found 0 bad Shells.\n",
         "Found %d bad Shells. The Surfaces are: " );

       print_volumes( logfile_ptr, bad_volume_list, 
         "Found 0 bad Volumes.\n",
         "Found %d bad Volumes: " );
#endif
       
       print_bodies( logfile_ptr, bad_body_list,
         "Found 0 Volumes with problems.\n",
         "Found %d Volumes with problems: " );
     }
   }

   // Advanced analysis output
   if( showDetails )
   {
     print_none( logfile_ptr, "\nHEALER DETAILED OUTPUT:\n" );
     print_none( logfile_ptr, "-----------------------\n" );

     // NOTE: need output for curves not on faces/vertices.  I
     //       think this is direct output of api_hh_get_bad_edges
     
     print_coedges( logfile_ptr, coedges_not_on_faces,
       "Found 0 CoEdges not on Surfaces.\n",
       "Found %d CoEdges not on Surfaces. Curves are: " );
     
     print_coedges( logfile_ptr, coedges_no_partner,
       "Found 0 CoEdges with no Curve partner.\n",
       "Found %d CoEdges with no Curve partner.  Curves are: " );
     
     print_vertices( logfile_ptr, vertices_not_on_faces,
       "Found 0 Vertices not on Surfaces.\n",
       "Found %d Vertices not on Surfaces.  Vertices are: " );
     
     print_vertices( logfile_ptr, vertices_not_on_edges,
       "Found 0 Vertices not on Curves.\n",
       "Found %d Vertices not on Curves.  Vertices are: " );
     
     print_vertices( logfile_ptr, vertices_edges_dont_meet,
       "Found 0 Vertices where Curves don't meet.\n",
       "Found %d Vertices where Curves don't meet.  Vertices are: " );
     
     print_curves( logfile_ptr, discontinuous_curves,
       "Found 0 discontinuous Curves.\n",
       "Found %d discontinuous Curves: " );
     
     print_curves( logfile_ptr, degenerate_curves,
       "Found 0 degenerate Curves.\n",
       "Found %d degenerate Curves: " );
     
     print_curves( logfile_ptr, self_intersecting_curves,
       "Found 0 self-intersecting Curves.\n",
       "Found %d self-intersecting Curves: " );
     
     print_curves( logfile_ptr, periodic_curves,
       "Found 0 periodic Curves.\n",
       "Found %d periodic Curves: " );
     
     print_curves( logfile_ptr, closed_curves,
       "Found 0 closed Curves.\n",
       "Found %d closed Curves: " );
     
     print_curves( logfile_ptr, short_edges,
       "Found 0 short Curves.\n",
       "Found %d short Curves: " );
     
     print_curves( logfile_ptr, tangent_edges,
       "Found 0 tangent Curves.\n",
       "Found %d tangent Curves: " );
     
     print_curves( logfile_ptr, convex_edges,
       "Found 0 convex Curves.\n",
       "Found %d convex Curves: " );
     
     print_curves( logfile_ptr, concave_edges,
       "Found 0 concave Curves.\n",
       "Found %d concave Curves: " );
     
     print_loops( logfile_ptr, loops_not_on_faces,
       "Found 0 Loops that are not on Surfaces.\n",
       "Found %d Loops not on Surfaces. Curves are: " );
     
     print_loops( logfile_ptr, loops_disoriented,
       "Found 0 disoriented Loops.\n",
       "Found %d disoriented Loops. Curves are: " );
     
     print_loops( logfile_ptr, loops_gaps,
       "Found 0 Loops with gaps.\n",
       "Found %d Loops with gaps. Curves are: " );
     
     print_loops( logfile_ptr, loops_open,
       "Found 0 open Loops.\n",
       "Found %d open Loops. Curves are: " );
     
     print_surfaces( logfile_ptr, discontinuous_surfaces,
       "Found 0 discontinuous Surfaces.\n",
       "Found %d discontinuous Surfaces: " );
     
     print_surfaces( logfile_ptr, degenerate_surfaces,
       "Found 0 degenerate Surfaces.\n",
       "Found %d degenerate Surfaces: " );
     
     print_surfaces( logfile_ptr, self_intersecting_surfaces,
       "Found 0 self-intersecting Surfaces.\n",
       "Found %d self-intersecting Surfaces: " );
     
     print_surfaces( logfile_ptr, periodic_surfaces,
       "Found 0 periodic Surfaces.\n",
       "Found %d periodic Surfaces: " );
     
     print_surfaces( logfile_ptr, closed_surfaces,
       "Found 0 closed Surfaces.\n",
       "Found %d closed Surfaces: " );
   }
   
   // Reset the logfile
   if( logfile_ptr )
   {
      hh_reset_bhl_log_file();
      CubitUtil::reset_file_ptr();
   }

   return CUBIT_SUCCESS;
}

void AcisHealerTool::list_incremental()
{
  char on_off[2][4];
  strcpy( on_off[0], "off" );
  strcpy( on_off[1], "on" );
  char override[2][40];
  strcpy( override[0], "" );
  strcpy( override[1], " (but overridden by Geombuild 'off')" );

  PRINT_INFO( "Healer 'Incremental' Settings:\n" );
//  if( incPreprocess == CUBIT_FALSE )
//    PRINT_INFO( " PreProcess:      %s (this is not recommended)\n", on_off[incPreprocess] );
//  else
    PRINT_INFO( " PreProcess:      on\n" );
  PRINT_INFO( " Simplify:        %s\n", on_off[incSimplify] );
  PRINT_INFO( " Stitch:          %s\n", on_off[incStitch] );
  PRINT_INFO( " Geombuild:       %s\n", on_off[incGeombuild] );
  
  if( incGeombuild == CUBIT_FALSE )
  {
    PRINT_INFO( "    Analytic:        %s%s\n", on_off[incAnalytic], override[incAnalytic] );
    PRINT_INFO( "    Isospline:       %s%s\n", on_off[incIsospline], override[incIsospline] );
    //PRINT_INFO( "    Reblend:         %s%s\n", on_off[incReblend], override[incReblend] );
    PRINT_INFO( "    Sharp Edge:      %s%s\n", on_off[incSharpedge], override[incSharpedge] );
    PRINT_INFO( "    Generic Spline:  %s%s\n", on_off[incGenericspline], override[incGenericspline] );
  }
  else
  {
    PRINT_INFO( "    Analytic:        %s\n", on_off[incAnalytic] );
    PRINT_INFO( "    Isospline:       %s\n", on_off[incIsospline] );
    //PRINT_INFO( "    Reblend:         %s\n", on_off[incReblend] );
    PRINT_INFO( "    Sharp Edge:      %s\n", on_off[incSharpedge] );
    PRINT_INFO( "    Generic Spline:  %s\n", on_off[incGenericspline] );
  }

//  if( incWrapup == CUBIT_FALSE )
//    PRINT_INFO( "    Wrapup:          %s (this is not recommended)\n", on_off[incWrapup] );
//  else
    PRINT_INFO( "    Wrapup:          on\n" );

//  if( incPostprocess == CUBIT_FALSE )
//    PRINT_INFO( " PostProcess:     %s (this is not recommended)\n", on_off[incPostprocess] );
//  else
    PRINT_INFO( " PostProcess:     on\n" );
}

void AcisHealerTool::list_onshow_flgs()
{
  char on_off[2][4];
  strcpy( on_off[0], "off" );
  strcpy( on_off[1], "on" );

  PRINT_INFO( "Healer 'OnShow' Settings:\n" );
  switch (showMethod)
  {
  case 0:
    PRINT_INFO( "  Show Method:        none (won't highlight or draw bad entities)\n" );
    break;
  case 1:
    PRINT_INFO( "  Show Method:        highlight\n" );
    break;
  case 2:
    PRINT_INFO( "  Show Method:        draw\n" );
    break;
  }

  //PRINT_INFO( "  Show Details: %s\n", on_off[showDetails] );
  PRINT_INFO( "  Show Summary:       %s\n", on_off[showSummary] );
  PRINT_INFO( "  Show After Healing: %s\n", on_off[showOnHeal] );
  PRINT_INFO( "  Show Bad Vertices:  %s\n", on_off[showBadVertices] );
  PRINT_INFO( "  Show Bad Curves:    %s\n", on_off[showBadCurves] );
  PRINT_INFO( "  Show Bad CoEdges:   %s\n", on_off[showBadCoEdges] );
  //PRINT_INFO( "  Show Bad Loops:     %s\n", on_off[showBadLoops] );
  //PRINT_INFO( "  Show Bad Surfaces:  %s\n", on_off[showBadSurfaces] );
  //PRINT_INFO( "  Show Bad Shells:    %s\n", on_off[showBadShells] );
  //PRINT_INFO( "  Show Bad Volumes:   %s\n", on_off[showBadVolumes] );
  PRINT_INFO( "  Show Bad Bodies:    %s\n", on_off[showBadBodies] );
}

void AcisHealerTool::incremental_presummary( DLIList<Body*>& body_list, FILE* logfile_ptr,
                                             CubitBoolean keep_old )
{
  char tmp[200];
  double tol;
  double tol1, tol2;

  CubitUtil::set_file_ptr( logfile_ptr );
  DLIList<CubitEntity*> temp_list;
  CAST_LIST_TO_PARENT(body_list, temp_list);    
  if( keep_old == CUBIT_FALSE )
    CubitUtil::list_entity_ids( "Healing steps to be performed on Bodies ", 
                                temp_list, 80, ":\n" );
  else
    CubitUtil::list_entity_ids( "Healing steps to be performed on copies of Bodies ", 
                                temp_list, 80, ":\n" );
  CubitUtil::reset_file_ptr();

  if( doPreprocess == TRUE )
    print_none( logfile_ptr, "* Preprocess\n" );

  // Simplify
  if( doSimplify == TRUE )
  {
    tol = get_simplify_tol( NULL );
    if( tol == -1.0 ) // Healer calculates tolerance
      print_none( logfile_ptr, "* Simplify\n" );
    else
    {
      sprintf( tmp, "* Simplify (at tolerance %f)\n", tol );
      print_none( logfile_ptr, tmp );
    }
  }

  // Stitch
  
  if( doStitch == TRUE )
  {
    tol1 = get_stitch_min_tol( NULL );
    tol2 = get_stitch_max_tol( NULL );
    if( tol1 == -1.0 && tol2 == -1.0 )
      print_none( logfile_ptr, "* Stitch\n" );
    else if( tol1 != -1.0 && tol2 == -1.0 )
    {
      sprintf( tmp, "* Stitch (at Stitch Min Tol %f)\n", tol1 );
      print_none( logfile_ptr, tmp );
    }
    else if( tol1 == -1.0 && tol2 != -1.0 )
    {
      sprintf( tmp, "* Stitch (at Stitch Max Tol %f)\n", tol2 );
      print_none( logfile_ptr, tmp );
    }
    else
    {
      sprintf( tmp, "* Stitch (at Stitch Min / Max Tol %f / %f)\n", tol1, tol2 );
      print_none( logfile_ptr, tmp );
    }
  }

  if( doGeomBuild == TRUE )
  {
    tol = get_geombuild_tol( NULL );
    if( tol == -1.0 )
      print_none( logfile_ptr, "* Geombuild\n" );
    else
    {
      sprintf( tmp, "* Geombuild (at tolerance %f)\n", tol );
      print_none( logfile_ptr, tmp );
    }
  
    if( doAnalytic == TRUE )
    {
      tol = get_analytic_tol( NULL );
      if( tol == -1.0 )
        print_none( logfile_ptr, "  * Analytic\n" );
      else
      {
        sprintf( tmp, "  * Analytic (at tolerance %f)\n", tol );
        print_none( logfile_ptr, tmp );
      }
    }
    
    if( doIsospline == TRUE )
    {
      tol = get_isospline_tol( NULL );
      if( tol == -1.0 )
        print_none( logfile_ptr, "  * Isospline\n" );
      else
      {
        sprintf( tmp, "  * Isospline (at tolerance %f)\n", tol );
        print_none( logfile_ptr, tmp );
      }
    }
    
    //  if( doReblend == TRUE )
    //    print_none( logfile_ptr, "  * Reblend\n" );
    
    if( doSharpEdge == TRUE )
      print_none( logfile_ptr, "  * Sharp Edge\n" );
    
    if( doGenericSpline == TRUE )
      print_none( logfile_ptr, "  * Generic Spline\n" );
    
    if( doWrapup == TRUE )
      print_none( logfile_ptr, "  * Wrapup\n" );
  }

  if( doPostprocess == TRUE )
    print_none( logfile_ptr, "* PostProcess\n" );

  print_none( logfile_ptr, "\n" );
}

void AcisHealerTool::measure_filter(DLIList<RefEntity*> &ref_entities,
                                    double measure_value,
                                    double tolerance)
{
  int i;
  DLIList<RefEntity*> temp_entities;
  for (i = ref_entities.size(); i > 0; i--) {
    RefEntity *this_entity = ref_entities.get_and_step();
   if(fabs(this_entity->measure() - measure_value) < tolerance)
     temp_entities.append(this_entity);
  }
  ref_entities = temp_entities;
}

CubitStatus
AcisHealerTool::force_simplify_to_plane( DLIList<RefFace*> &ref_face_list, 
                                         DLIList<Body*>& new_body_list, 
                                         CubitBoolean keep_old_body )
{
  return force_simplify( 1, ref_face_list, new_body_list, keep_old_body );
}

CubitStatus
AcisHealerTool::force_simplify_to_cylinder( DLIList<RefFace*> &ref_face_list, 
                                            DLIList<Body*>& new_body_list,
                                            CubitBoolean keep_old_body )
{
  return force_simplify( 2, ref_face_list, new_body_list, keep_old_body );
}

CubitStatus
AcisHealerTool::force_simplify_to_cone( DLIList<RefFace*> &ref_face_list, 
                                        DLIList<Body*>& new_body_list,
                                        CubitBoolean keep_old_body )
{
  return force_simplify( 3, ref_face_list, new_body_list, keep_old_body );
}

CubitStatus
AcisHealerTool::force_simplify_to_sphere( DLIList<RefFace*> &ref_face_list, 
                                          DLIList<Body*>& new_body_list,
                                          CubitBoolean keep_old_body )
{
  return force_simplify( 4, ref_face_list, new_body_list, keep_old_body );
}

CubitStatus
AcisHealerTool::force_simplify_to_torus( DLIList<RefFace*> &ref_face_list, 
                                         DLIList<Body*>& new_body_list,
                                         CubitBoolean keep_old_body )
{
  return force_simplify( 5, ref_face_list, new_body_list, keep_old_body );
}

CubitStatus
AcisHealerTool::force_simplify( int simplify_type, DLIList<RefFace*> &ref_face_list, 
                               DLIList<Body*>& new_body_list,
                               CubitBoolean keep_old_body )
{
  outcome result;

  CubitBoolean delete_attribs =
     (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  // Copy the incoming ref_face_list since we will be pulling
  // surfaces out of it.
  DLIList<RefFace*> copied_ref_face_list = ref_face_list;

  int i;

  Body* body_ptr;
  RefFace* ref_face_ptr;

  BODY* BODY_ptr;
  BODY* copied_BODY_ptr;
  FACE *FACE_ptr;
 
  copied_ref_face_list.reset();
  while( copied_ref_face_list.size() )
  {
    int success = 0;
    DLIList<RefFace*> simplify_face_list;
    DLIList<FACE*> FACE_list;
    if( AcisToolUtil::get_copied_FACES_of_body( copied_ref_face_list, FACE_list, 
      simplify_face_list, copied_BODY_ptr ) == CUBIT_FAILURE )
      break;
    
    // Get original Body and BODY
    body_ptr = AcisToolUtil::get_body_of_ENTITY( copied_BODY_ptr );
    BODY_ptr =  BodyACIS::get_BODY_ptr(body_ptr);

	 // Now cleanout the owner attributes from the copied BODY, if required
    if( delete_attribs )
       AcisModifyEngine::instance()->get_acis_query_engine()->remove_cubit_owner_attrib_in_BODY( copied_BODY_ptr );
    
    // Now, simplify the surfaces on this body
    FACE_list.reset();
    simplify_face_list.reset();
    for( i=0; i<FACE_list.size(); i++ )
    {
      FACE_ptr = FACE_list.get_and_step();
      ref_face_ptr = simplify_face_list.get_and_step();

      int surface_type = (&(FACE_ptr->geometry()->equation()))->type();

      switch (simplify_type)
      {
      case 1: // PLANE
        if( surface_type == plane_type )
          PRINT_WARNING( "Surface %d in body %d is already a plane - no simplification required\n",
            ref_face_ptr->id(), body_ptr->id() );
        else
        {
          result = api_hh_force_simplify_to_plane( FACE_ptr );
          if( !result.ok() )
          {
            AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error(result);
            PRINT_ERROR( "Unable to force surface %d in body %d to be planar\n", 
              ref_face_ptr->id(), body_ptr->id() );
            continue;
          }
          else
          {
            // Check to see if the surface was indeed changed to a plane
            surface_type = (&(FACE_ptr->geometry()->equation()))->type();
            if( surface_type == plane_type )
            {
              PRINT_INFO( "Successfully changed surface %d on body %d into a plane\n", 
                          ref_face_ptr->id(), body_ptr->id() );
              success++;
            }
            else
              PRINT_ERROR( "Unable to change surface %d on body %d into a plane\n",
                           ref_face_ptr->id(), body_ptr->id() );
          }
        }
        break;
      case 2: // CYLINDER
        if( surface_type == cone_type )
          PRINT_WARNING( "Surface %d in body %d is already cylindrical/conic - no simplification required\n",
            ref_face_ptr->id(), body_ptr->id() );
        else
        {
          result = api_hh_force_simplify_to_cylinder( FACE_ptr );
          if( !result.ok() )
          {
            AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error(result);
            PRINT_ERROR( "Unable to force surface %d in body %d to be cylindrical\n", 
              ref_face_ptr->id(), body_ptr->id() );
            continue;
          }
          else
          {
            // Check to see if the surface was indeed changed to a cylinder/cone
            surface_type = (&(FACE_ptr->geometry()->equation()))->type();
            if( surface_type == cone_type )
            {
              PRINT_INFO( "Successfully changed surface %d on body %d into a cylinder\n", 
                ref_face_ptr->id(), body_ptr->id() );
              success++;
            }
            else
              PRINT_ERROR( "Unable to change surface %d on body %d into a cylinder\n",
              ref_face_ptr->id(), body_ptr->id() );
          }
        }
        break;
      case 3: // CONE
        if( surface_type ==  cone_type )
          PRINT_WARNING( "Surface %d in body %d is already conic - no simplification required\n",
            ref_face_ptr->id(), body_ptr->id() );
        else
        {
          result = api_hh_force_simplify_to_cone( FACE_ptr );
          if( !result.ok() )
          {
            AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error(result);
            PRINT_ERROR( "Unable to force surface %d in body %d to be conic\n", 
              ref_face_ptr->id(), body_ptr->id() );
            continue;
          }
          else
          {
            // Check to see if the surface was indeed changed to a cone
            surface_type = (&(FACE_ptr->geometry()->equation()))->type();
            if( surface_type == cone_type )
            {
              PRINT_INFO( "Successfully changed surface %d on body %d into a cone\n", 
                          ref_face_ptr->id(), body_ptr->id() );
              success++;
            }
            else
              PRINT_ERROR( "Unable to change surface %d on body %d into a cone\n",
                           ref_face_ptr->id(), body_ptr->id() );
          }
        }
        break;
      case 4: // SPHERE
        if( surface_type ==  sphere_type )
          PRINT_WARNING( "Surface %d in body %d is already spherical - no simplification required\n",
            ref_face_ptr->id(), body_ptr->id() );
        else
        {
          result = api_hh_force_simplify_to_sphere( FACE_ptr );
          if( !result.ok() )
          {
            AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error(result);
            PRINT_ERROR( "Unable to force surface %d in body %d to be spherical\n", 
              ref_face_ptr->id(), body_ptr->id() );
            continue;
          }
          else
          {
            // Check to see if the surface was indeed changed to a sphere
            surface_type = (&(FACE_ptr->geometry()->equation()))->type();
            if( surface_type == sphere_type )
            {
              PRINT_INFO( "Successfully changed surface %d on body %d into a sphere\n", 
                          ref_face_ptr->id(), body_ptr->id() );
              success++;
            }
            else
              PRINT_ERROR( "Unable to change surface %d on body %d into a sphere\n",
                           ref_face_ptr->id(), body_ptr->id() );
          }
        }
        break;
      case 5: // TORUS
        if( surface_type ==  torus_type )
          PRINT_WARNING( "Surface %d in body %d is already torroidal - no simplification required\n",
            ref_face_ptr->id(), body_ptr->id() );
        else
        {
          result = api_hh_force_simplify_to_torus( FACE_ptr );
          if( !result.ok() )
          {
            AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error(result);
            PRINT_ERROR( "Unable to force surface %d in body %d to be torroidal\n", 
              ref_face_ptr->id(), body_ptr->id() );
            continue;
          }
          else
          {
            // Check to see if the surface was indeed changed to a torus
            surface_type = (&(FACE_ptr->geometry()->equation()))->type();
            if( surface_type == torus_type )
            {
              PRINT_INFO( "Successfully changed surface %d on body %d into a torus\n", 
                          ref_face_ptr->id(), body_ptr->id() );
              success++;
            }
            else
              PRINT_ERROR( "Unable to change surface %d on body %d into a torus\n",
                           ref_face_ptr->id(), body_ptr->id() );
          }
        }
        break;
      default:
        PRINT_ERROR( "undefined simplification type specified\n" );
      }
    }

    if( success == 0 )
    {
      PRINT_WARNING( "No surfaces on body %d were simplified\n", body_ptr->id() );
      api_delent( copied_BODY_ptr );
      continue;
    }
    
    // If we've made it this far, the copied_BODY has been
    // modified and we can update it in CUBIT
    Body* new_body_ptr = AcisToolUtil::get_new_Body( body_ptr, BODY_ptr, 
      copied_BODY_ptr, keep_old_body );
    
    if( new_body_ptr!=NULL && new_body_ptr!=body_ptr )
    {
      PRINT_INFO( "Created new body %d (simplified %d surfaces)\n", 
        new_body_ptr->id(), success );
      new_body_list.append( new_body_ptr );
    }
    else if( new_body_ptr!=NULL && new_body_ptr==body_ptr )
    {
      PRINT_INFO( "Modified body %d (simplified %d surfaces)\n", 
        body_ptr->id(), success );
    }
    else
      PRINT_WARNING( "Body %d was not modified\n", body_ptr->id() );
  }
  
  return CUBIT_SUCCESS;
}

#endif

