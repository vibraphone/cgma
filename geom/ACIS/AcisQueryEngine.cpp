//-------------------------------------------------------------------------
// Filename      : AcisQueryEngine.cpp
//
// Purpose       : Performs all query-type operations on acis geometry
//
// Special Notes :
//
// Creator       : Tim Tautges
//
// Creation Date : 2/01
//
// Owner         : Tim Tautges
//-------------------------------------------------------------------------

// This macro can be used to eliminate compiler warnings, in the case
// where the healer is not being linked in.  Stands for use-healer-argument.
#ifdef ACIS_HEALER
  #define USEHARG(s) s
#else
  #define USEHARG(s)
#endif

#ifdef NT
  #include <direct.h>
#ifndef PATH_MAX
  #define PATH_MAX _MAX_PATH
#endif
#else
  #include <unistd.h>
#endif


// ********** BEGIN STANDARD INCLUDES         **********
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
// ********** END STANDARD INCLUDES           **********

// ********** BEGIN ACIS INCLUDES             **********
#include "acis.hxx"
#include "version.hxx"
#include "cstrapi.hxx"
#include "eulerapi.hxx"
#include "intrapi.hxx"
#include "api.hxx"
#include "kernapi.hxx"
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
#include "allsfdef.hxx"
#include "curdef.hxx"
#include "sweeping.hxx"
#include "tensor.hxx"
#include "alltop.hxx"
#include "wire.hxx"
#include "transent.hxx"
#include "af_api.hxx"
#include "q_wire.hxx"
#include "sgquery.hxx"
#include "chk_stat.hxx"
#include "chk.hxx"
#include "boolean.hxx"
#include "boolapi.hxx"
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
#include "law.hxx"
#include "wire_qry.hxx"
#include "intcucu.hxx"
#include "entwray.hxx"
#include "wire_utl.hxx"
#include "curextnd.hxx"
#include "extend.hxx"
#include "acistype.hxx"
#include "nms_attr.hxx"

#include "warp_api.hxx"
#include "at_str.hxx"

#include "ptcont.hxx"
#include "body.hxx"
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
#include "getowner.hxx"

#if CUBIT_ACIS_VERSION >= 1600
#include "clash_bodies.hxx"
#endif

#ifdef ACIS_IGES_TRANSLATOR
#include "errorcode.err"
#include "acisiges_api.hxx"
#include "SPAXUnit.h"
#include "acis_unit_api.h"
#endif

#ifdef ACIS_PROE_TRANSLATOR
   #include "proe/proehusk/api/proeapi.hxx"
#endif
#ifdef ACIS_CATIA_TRANSLATOR
   #include "cathusk/chl_api/apiinit.hxx"
   #include "cathusk/chl_api/apiread.hxx"
   #include "cathusk/chl_api/apiwrite.hxx" 
#endif
#ifdef ACIS_STEP_TRANSLATOR
#include "acisstep_api.hxx"
#endif

#if defined( ACIS_IGES_TRANSLATOR) || defined( ACIS_STEP_TRANSLATOR )
  #include "SPAXBase.h"
  #include "SPAXBoolean.h"
  #include "SPAXProgressReportCB.h"
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
#include "AcisModifyEngine.hpp"
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
#include "MergeTool.hpp"

#include "AcisQueryEngine.hpp"
#include "AcisFacetManager.hpp"
#include "AcisHealerTool.hpp"

#include "DoubleListItem.hpp"

#include "attrib_cubit_owner.hpp"
#include "attrib_snl_simple.hpp"

#include "CubitVector.hpp"

#include "AnalyticGeometryTool.hpp"
#include "GeometryUtil.hpp"


#include "DLIList.hpp"

#include "CastTo.hpp"

#include "GMem.hpp"
#include "GfxDebug.hpp"
#include "AcisFeatureEngine.hpp"


const int AcisQueryEngine::MAX_NUM_CURVE_POINTS  = 750;
const int AcisQueryEngine::AQE_SUBMINOR_VERSION = 0;

AcisQueryEngine* AcisQueryEngine::instance_ = 0;

AcisQueryEngine::~AcisQueryEngine()
{
   api_terminate_generic_attributes();
   api_terminate_spline();
   api_terminate_kernel();
   api_terminate_intersectors();
   api_terminate_constructors();
   api_terminate_euler_ops();
#ifdef ACIS_IGES_TRANSLATOR
   api_terminate_xiges();
#endif
#ifdef ACIS_STEP_TRANSLATOR
   api_terminate_xstep();
#endif
#ifdef ACIS_PROE_TRANSLATOR
   api_terminate_proe();
#endif
#ifdef  ACIS_CATIA_TRANSLATOR
   api_terminate_catia();
#endif
     //api_terminate_faceter();
   delete facetManager;
   api_stop_modeller();

    terminate_base();

   instance_ = 0;
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

const type_info& AcisQueryEngine::entity_type_info() const
{
   return typeid(AcisQueryEngine);
}


CubitStatus AcisQueryEngine::save_temp_geom_file(DLIList<TopologyBridge*> &ref_entity_list,
                                                 const char *filename,
                                                 const CubitString &cubit_version,
                                                 CubitString &created_file,
                                                 CubitString &created_file_type) 
{

  int size_before = ref_entity_list.size();
  CubitString temp_filename(filename);
  temp_filename += ".sat";
  
  if( export_solid_model( ref_entity_list, temp_filename.c_str(), "ACIS_SAT",
                          cubit_version ) == CUBIT_FAILURE )
  {
    PRINT_ERROR( "Error occured while trying to export temporary ACIS file\n");
    return CUBIT_FAILURE;
  }

  int size_after = ref_entity_list.size();

  if( size_before > size_after )
  {
    created_file +=  temp_filename;
    created_file_type += "ACIS_SAT";  
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Export the current CUBIT geometry (everything in the 
//                 Model) to a solid model format. Valid file types are:
//
//                 Valid file types are:
//                       "ACIS_SAT"    --  ACIS ASCII (SAT) file format
//                       "ACIS_SAB"    --  ACIS BINARY (SAB) file format
//                       "ACIS_DEBUG"  --  ACIS DEBUG file format
//                       "IGES"        --  IGES file
//                       "STEP"        --  STEP file
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 03/17/99
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::export_solid_model( DLIList<TopologyBridge*>& bridge_list,
                                                    const char* filename,
                                                    const char* filetype,
                                                    const CubitString &cubit_version,
                                                    const char* logfile_name) 
{
     // Check to make sure we are exporting valid types
   if ( strcmp(filetype, "ACIS_SAT") != 0 &&
        strcmp(filetype, "ACIS_SAB") != 0 &&
        strcmp(filetype, "ACIS_DEBUG") != 0 &&
        strcmp(filetype, "IGES") != 0 &&
#ifdef ACIS_CATIA_TRANSLATOR
        strcmp(filetype, "CATIA") != 0 &&
#endif
        strcmp(filetype, "STEP") != 0 )
   {
      //PRINT_ERROR("The specified file type, %s, is not supported!\n", filetype );
      return CUBIT_FAILURE;
   }

   outcome result;
   
   ENTITY* entity_ptr = NULL;
   ENTITY* copied_entity_ptr = NULL;

   FILE* file_ptr = NULL;

   int body_count = 0, lump_count = 0,
      face_count = 0,
      edge_count = 0,
      vertex_count = 0;
   
   ENTITY_LIST entity_list;
   ENTITY_LIST copied_entity_list; 
     
   int i;

   // Remove all healing attributes from all bodies.  Done here so that if user
   // is exporting individual entities it still cleans everything off.
   DLIList<Body*> all_bodies;
   GeometryQueryTool::instance()->bodies(all_bodies);
   CubitBoolean handled = CUBIT_FALSE;

   // Now actually write the entities
   for (i = bridge_list.size(); i > 0; i--) 
   {
      TopologyBridge* bridge = bridge_list.get();
      handled = CUBIT_TRUE;
      
      if( BodyACIS* body = CAST_TO( bridge, BodyACIS ) )
      {
         BODY *BODY_ptr =  body->get_BODY_ptr();
         
         // Make sure there is a valid associated ACIS BODY. If not, then
         // proceed to the next RefEntity.
         if (BODY_ptr == NULL)
           handled = CUBIT_FALSE;
         else
         {
           // If there are any healing attributes on the BODY, clean them off
           entity_list.add(BODY_ptr);
           body_count++;
           
           for( LUMP* LUMP_ptr = BODY_ptr->lump();
                LUMP_ptr != NULL; 
                LUMP_ptr = LUMP_ptr->next() )
           {
             lump_count++;
           }
         }
      }
      else if( SurfaceACIS* surf = CAST_TO( bridge, SurfaceACIS ) )
      { 
         FACE* FACE_ptr = surf->get_FACE_ptr();
         
         // Make sure there is a valid associated ACIS FACE. If not, then
         // proceed to the next RefEntity.
         if (FACE_ptr == NULL)
           handled = CUBIT_FALSE;
         else
         {
            // Since the user is requesting this FACE be exported freely,
            // copy it to a new entity (in case it's not really a free face)
            copy_single_entity((ENTITY*)FACE_ptr, copied_entity_ptr);
            entity_list.add(copied_entity_ptr);
            copied_entity_list.add(copied_entity_ptr);
            face_count++;
         }
      }
      else if( CurveACIS* curve = CAST_TO( bridge, CurveACIS ) )
      {
         EDGE* EDGE_ptr = curve->get_EDGE_ptr();
         
         // Make sure there is a valid associated ACIS EDGE. If not, then
         // proceed to the next RefEntity.
         if (EDGE_ptr == NULL)
           handled = CUBIT_FALSE;
         else
         {
            // Since the user is requesting this EDGE be exported freely,
            // copy it to a new entity (in case it's not really a free edge)
            copy_single_entity((ENTITY*)EDGE_ptr, copied_entity_ptr);
            entity_list.add(copied_entity_ptr);
            copied_entity_list.add(copied_entity_ptr);
            edge_count++;
         }
      }
      else if( PointACIS* point = CAST_TO( bridge, PointACIS ) )
      {
         VERTEX* VERTEX_ptr = point->get_VERTEX_ptr();
         
         // Make sure there is a valid associated ACIS VERTEX. If not, then
         // proceed to the next VERTEX.
         if (VERTEX_ptr == NULL)
           handled = CUBIT_FALSE;
         else
         {
            // Since the user is requesting this VERTEX be exported freely,
            // copy it to a new entity (in case it's not really a free vertex)
            copy_single_entity((ENTITY*)VERTEX_ptr, copied_entity_ptr);
            entity_list.add(copied_entity_ptr);
            copied_entity_list.add(copied_entity_ptr);
            vertex_count++;
         }
      }
      else
        handled = CUBIT_FALSE;

      if (CUBIT_TRUE == handled)
          // handled this entity - nullify in the list
        bridge_list.change_to(NULL);

        // either way, step the list
      bridge_list.step();
   }

   if (entity_list.count() == 0) 
      return CUBIT_SUCCESS;
   
     // ok, now compress the list, assuming we nullified something in it
   bridge_list.remove_all_with_value(NULL);
   
   // We now have the Acis ENTITY_LIST to be exported.  Write out
   // the correct type of file.
   
   // Export an ACIS SAT file
   if (strcmp(filetype, "ACIS_SAT") == 0 ||
      strcmp(filetype, "ACIS_SAB") == 0)
   {
      if( filename != NULL && strcmp( filename, "" ) )
      {
#ifdef NT  // NT requires "wb" type for binary - Unix doesn't
         if( strcmp(filetype, "ACIS_SAT") == 0 )
            file_ptr = fopen(filename, "w");  // ASCII output
         else
            file_ptr = fopen(filename, "wb"); // BINARY output
#else
         file_ptr = fopen(filename, "w");
#endif // NT
      }
      else
         file_ptr = NULL;
      if (file_ptr)
      {
         CubitString version = "Cubit ";
         version += cubit_version;
         FileInfo info;
         info.set_product_id(version.c_str());

         info.set_units(1.0);
	     api_set_file_info((FileId | FileUnits), info);
         if(strcmp(filetype, "ACIS_SAT") == 0)
            result = api_save_entity_list (file_ptr,
            CUBIT_TRUE,
            entity_list );
         else
            result = api_save_entity_list (file_ptr,
            CUBIT_FALSE,
            entity_list );
         
         // Make sure the save was successful
         if (!result.ok())
         {
            PRINT_ERROR("Acis could not save to the file: '%s'\n",
               filename);
            entity_list.clear();
            copied_entity_list.init();
            while ( (entity_ptr=copied_entity_list.next()) )
               api_delent( entity_ptr );
            copied_entity_list.clear();
            return CUBIT_FAILURE;
         }
         fclose (file_ptr);
      }
      
      else
      {
         PRINT_ERROR("Cannot create file: '%s'\n",
            filename);
         entity_list.clear();
         copied_entity_list.init();
         while ( (entity_ptr=copied_entity_list.next()) )
            api_delent( entity_ptr );
         return CUBIT_FAILURE;
      }
   }
   
     // Export an ACIS DEBUG file
   else if (strcmp(filetype, "ACIS_DEBUG") == 0)
   {
        // Create the file in "append" mode
      file_ptr = fopen(filename, "a");
      if (file_ptr)
      {
         for (int jj = 0; jj < entity_list.count(); jj++)
         {
            debug_entity( (BODY *) (entity_list)[jj], file_ptr);
         }
         fclose (file_ptr);
      }
      else
      {
         PRINT_ERROR("Cannot open file '%s' for appending debug data\n",
                     filename);
         entity_list.clear();
         copied_entity_list.init();
         while ( (entity_ptr=copied_entity_list.next()) )
            api_delent( entity_ptr );
         return CUBIT_FAILURE;
      }
   }

   else if( strcmp( filetype, "IGES" ) == 0 )
   {
#ifndef ACIS_IGES_TRANSLATOR
      copied_entity_list.init();
      while ( (entity_ptr=copied_entity_list.next()) )
         api_delent( entity_ptr );
      PRINT_ERROR( "The IGES translator is not licensed for this installation\n" );
      return CUBIT_FAILURE;
#else
      CubitString version = "Cubit ";
      version += cubit_version;

      // Author, organization, sending system, receiving system
      result = api_xiges_set_header( getenv("USERNAME"), getenv("USERDOMAIN"),
                                       version.c_str(), NULL, 1, 1.0 );
      if( !result.ok() )
      {
         ACIS_API_error(result);
         PRINT_WARNING("unable to define IGES header information for '%s'\n", filename);
      }

      char* logfilename = NULL;
      if( !logfile_name || !strcmp( logfile_name, "" ) )
         strcpy(logfilename, "iges_export.log");
      else
         logfilename = (char *)logfile_name;
   
      result = api_xiges_write( entity_list, (char *)filename, logfilename );

      if( !result.ok() )
      {
         PRINT_ERROR("Acis could not export to the IGES file: '%s'\n", filename);
         ACIS_API_error(result);
         entity_list.clear();
         copied_entity_list.init();
         while ( (entity_ptr=copied_entity_list.next()) )
            api_delent( entity_ptr );
         return CUBIT_FAILURE;
      }
#endif //ACIS_IGES_TRANSLATOR
   }
#ifdef ACIS_CATIA_TRANSLATOR
   else if( strcmp( filetype, "CATIA" ) == 0 )
   {
#ifndef ACIS_CATIA_TRANSLATOR
      copied_entity_list.init();
      while ( (entity_ptr=copied_entity_list.next()) )
         api_delent( entity_ptr );
      PRINT_ERROR( "The CATIA translator is not licensed for this installation\n" );
      return CUBIT_FAILURE;
#else
      CubitString version = "Cubit ";
      version += cubit_version;

      char* logfilename = NULL;
      if( !logfile_name || !strcmp( logfile_name, "" ) )
         strcpy(logfilename, "catia_export.log");
      else
         logfilename = (char *)logfile_name;
    
      result = api_catia_convert_acisentlist_to_catia(
         entity_list, (char *)filename, logfilename );
      if( !result.ok() )
      {
         PRINT_ERROR("Acis could not export to the CATIA file: '%s'\n", filename);
         ACIS_API_error(result);
         entity_list.clear();
         copied_entity_list.init();
         while ( (entity_ptr=copied_entity_list.next()) )
            api_delent( entity_ptr );
         return CUBIT_FAILURE;
      }
#endif //ACIS_CATIA_TRANSLATOR
   }
#endif
   else if( strcmp( filetype, "STEP" ) == 0 )
   {
#ifndef ACIS_STEP_TRANSLATOR
      entity_list.clear();
      copied_entity_list.init();
      while ( (entity_ptr=copied_entity_list.next()) )
         api_delent( entity_ptr );
#if defined(SGI)
      PRINT_ERROR("The STEP translator is not available for the 64-bit SGI platform\n");
#else 
      PRINT_ERROR( "The STEP translator is not licensed for this installation\n" );
#endif
      return CUBIT_FAILURE;
#else     
      if ( !stepInitialized ) {
        PRINT_ERROR("The STEP translater has not been properly initialized.\n");
        entity_list.clear();
        copied_entity_list.init();
        while ( (entity_ptr = copied_entity_list.next()) )
          api_delent( entity_ptr );
      
        return CUBIT_FAILURE;
      }
      
      char* logfilename = NULL;
      if( !logfile_name || !strcmp( logfile_name, "" ) )
         strcpy(logfilename, "step_export.log");
      else
         logfilename = (char *)logfile_name;

      result = api_xstep_write( entity_list, (char *)filename, logfilename ); 

//    BODY* block;
//    result = api_make_cuboid(10.0, 10.0, 10.0, block);
//    ENTITY_LIST elist;
//    elist.init();
//    elist.add(block);
//    result = api_step_convert_acisentlist_to_step(elist, "block.step", "block_step_write.log", NULL);
//    api_delent(block);
	
      if( !result.ok() )
      {
         PRINT_ERROR("ACIS could not export to the STEP file: '%s'\n", filename);
         ACIS_API_error(result);
         entity_list.clear();
         copied_entity_list.init();
         while ( (entity_ptr=copied_entity_list.next()) )
            api_delent( entity_ptr );
         return CUBIT_FAILURE;
      }
#endif
   }
   
   else
   {
      entity_list.clear();
      copied_entity_list.init();
      while ( (entity_ptr=copied_entity_list.next()) )
         api_delent( entity_ptr );
      PRINT_ERROR("Invalid ACIS filetype '%s'\n", filetype);
      return CUBIT_FAILURE;
   }
   entity_list.clear();
   copied_entity_list.init();
   while ( (entity_ptr=copied_entity_list.next()) )
      api_delent( entity_ptr );

   if( body_count || face_count || edge_count || vertex_count )
      PRINT_INFO( "\nExported:" );

   int flg = 0;

   if( body_count )
   {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( DEBUG_FLAG( 153 ) )
      {
        if( body_count == 1 )
           PRINT_INFO( "%4d Body\n", body_count );
        else
           PRINT_INFO( "%4d Bodies\n", body_count );
      }

      if( body_count == 1 )
         PRINT_INFO( "%4d Volume\n", lump_count );
      else
         PRINT_INFO( "%4d Volumes\n", lump_count );
      
   }
   if( face_count )
   {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( face_count == 1 )
         PRINT_INFO( "%4d Surface\n", face_count );
      else
         PRINT_INFO( "%4d Surfaces\n", face_count );
   }
   if( edge_count )
   {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( edge_count == 1 )
         PRINT_INFO( "%4d Curve\n", edge_count );
      else
         PRINT_INFO( "%4d Curves\n", edge_count );
   }
   if( vertex_count )
   {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( vertex_count == 1 )
         PRINT_INFO( "%4d Vertex\n", vertex_count );
      else
         PRINT_INFO( "%4d Vertices\n", vertex_count );
   }
   PRINT_INFO( "\n" );

   return CUBIT_SUCCESS;
}


CubitStatus 
AcisQueryEngine::import_temp_geom_file(FILE* file_ptr, 
                                       const char* /*file_name*/,
                                       const char* file_type,
                                       DLIList<TopologyBridge*> &bridge_list )
{
  CubitStatus status;
  ENTITY_LIST entity_list;

  //make sure that file_type == "ACIS_SAT"
  if( !strcmp( file_type,"ACIS_SAT") )
    status = read_acis_file( file_ptr, file_type, entity_list ); 
  else 
    return CUBIT_FAILURE;

  PRINT_INFO("Read %d ACIS Entities from the input file\n",
            entity_list.count());
  
    // create topology bridges from the ACIS entities
  CubitBoolean print_results = CUBIT_TRUE;
  status = restore_entity_list(entity_list, bridge_list, print_results); // lots of default argument values

  entity_list.clear();

#ifdef CUBIT_GUI
#ifndef NO_USAGE_TRACKING
  if(status)
      GUIInterface::instance()->log_import_acis_op();
#endif
#endif // CUBIT_GUI

  return status;
}

CubitStatus
AcisQueryEngine::remove_refinements( ENTITY_LIST& list )
{
  list.init();
  while (ENTITY* ent = list.next())
    api_set_entity_refinement( ent, NULL, TRUE );
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Combine two ENTITY_LISTs
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/15/04
//-------------------------------------------------------------------------
void AcisQueryEngine::listcat( ENTITY_LIST& result,
                               const ENTITY_LIST& other )
{
  other.init();
  while (ENTITY* ptr = other.next())
    result.add( ptr );
}

#if defined( ACIS_IGES_TRANSLATOR) || defined( ACIS_STEP_TRANSLATOR )
void progress_report_function (double iProgress, SPAXBoolean& oAbort)
{
  AppUtil::instance()->progress_tool()->percent( iProgress );
  return;
}
#endif

//-------------------------------------------------------------------------
// Purpose       : Reads in geometry in geometry and creates
//                 the necessary Reference entities associated with the 
//                 input geometry. Has ability to selectively read
//                 in either bodies, free surfaces, free curves or free 
//                 vertices from the file.
//
//                 Valid file types are:
//                       "ACIS_SAT"    --  ACIS ASCII (SAT) file format
//                       "ACIS_SAB"    --  ACIS BINARY (SAB) file format
//                       "IGES"        --  IGES file
//                       "STEP"        --  STEP file
//
// Special Notes :
//
// Creator       : Steve Storm
//
// Creation Date : 03/17/99
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::import_solid_model(
  const char* file_name,
  const char* file_type,
  DLIList<TopologyBridge*> &imported_entities,
  CubitBoolean print_results,
  const char* logfile_name,
  CubitBoolean heal_step,
  CubitBoolean import_bodies,
  CubitBoolean import_surfaces,
  CubitBoolean import_curves,
  CubitBoolean import_vertices,
  CubitBoolean free_surfaces) 
{
  int i;
  CpuTimer import_solid_model_timer;
  CubitBoolean step_import = CUBIT_FALSE;
  CubitBoolean iges_import = CUBIT_FALSE;


    // Check to make sure the requested type is supported
  if (strcmp(file_type, "ACIS_SAT") != 0 &&
      strcmp(file_type, "ACIS_SAB") != 0 &&
      strcmp(file_type, "IGES") != 0 &&
      strcmp(file_type, "STEP") != 0 &&
#ifdef ACIS_CATIA_TRANSLATOR
      strcmp(file_type, "CATIA") != 0 &&
#endif
      strcmp(file_type, "PROE") !=0)
  {
    PRINT_ERROR("The specified file type, %s, is not supported!\n", file_type);
    return CUBIT_FAILURE;
  }


  ENTITY_LIST entity_list, tmp_list;
  outcome result;
  CubitStatus status;

  // read acis sat and sab files
  if( strcmp(file_type, "ACIS_SAT") == 0 ||
      strcmp(file_type, "ACIS_SAB") == 0 )
  {
    status = read_acis_file(file_name, file_type, entity_list);
    if (CUBIT_FAILURE == status)
      return status;
  }
  // read an iges file
  else if(  strcmp(file_type, "IGES") == 0 )
  {
#ifndef ACIS_IGES_TRANSLATOR
    PRINT_ERROR( "The IGES translator is not licensed for this installation\n" );
    return CUBIT_FAILURE;
#else
	  iges_import = CUBIT_TRUE;
    status = read_iges_file(file_name, logfile_name, entity_list, heal_step);
    if (CUBIT_FAILURE == status)
      return status;

#endif //ACIS_IGES_TRANSLATOR
  }
#ifdef ACIS_CATIA_TRANSLATOR
  // read a catia file
  else if(  strcmp(file_type, "CATIA") == 0 )
  {
	iges_import = CUBIT_TRUE;
#ifndef ACIS_CATIA_TRANSLATOR
    PRINT_ERROR( "The CATIA translator is not licensed for this installation\n" );
    return CUBIT_FAILURE;
#else
    status = read_catia_file(file_name, logfile_name, entity_list);
    if (CUBIT_FAILURE == status)
      return status;
#endif //ACIS_CATIA_TRANSLATOR
  }
#endif
  // read a step file
  else if( strcmp( file_type, "STEP") == 0 )
  {
#ifndef ACIS_STEP_TRANSLATOR
#if defined(SGI)
    PRINT_ERROR("The STEP translator is not available for the 64-bit SGI platform\n");
#else 
    PRINT_ERROR( "The STEP translator is not licensed for this installation\n" );
#endif
    //char* dummy2 = (char *)logfile_name;
    return CUBIT_FAILURE;
#else
    step_import = CUBIT_TRUE;
    ENTITY_LIST my_list;

    // read the file
    status = read_step_file(file_name, logfile_name, my_list, heal_step);
    if (CUBIT_FAILURE == status)
      return status;

    // get part name attributes stored in the STEP file
    read_step_part_names(my_list);

    // Parts cut in Pro/E Mechanica and exported as STEP come in as a single body.
    // so this will auto separate the body 
    auto_separate_step_body(my_list, entity_list);
    my_list.clear();
#endif //ACIS_STEP_TRANSLATOR
  }
  // read a pro/e file
  else if( strcmp( file_type, "PROE") == 0 )
  {
#ifndef ACIS_PROE_TRANSLATOR
    PRINT_ERROR( "the PROE translator is not licensed for this installation\n" );
    //char* dummy2 = (char *)logfile_name;
    return CUBIT_FAILURE;
#else
    status = read_proe_file(file_name, entity_list);
    if (CUBIT_FAILURE == status)
      return status;
#endif //ACIS_PROE_TRANSLATOR
  }
   

  if( entity_list.count() == 0 )
  {
    PRINT_ERROR( "did not read anything from the input file '%s'\n", file_name );
     PRINT_INFO("You may be trying to read a file that is from an Acis version\n"
                "newer than you are using.\n");
    return CUBIT_FAILURE;
  }
  else
     PRINT_INFO("Read %d ACIS Entities from the input file\n",
                entity_list.count());
  
    // create topology bridges from the ACIS entities
  status = restore_entity_list(entity_list, imported_entities,
                                           print_results,
                                           step_import, 
                                           iges_import, 			
                                           import_bodies,
                                           import_surfaces, import_curves,
                                           import_vertices, free_surfaces );

    // Get the CPU time for the entire import solid model operation
  PRINT_DEBUG_3( "CPU time taken for importing the "
                 "solid model: %f secs\n",
                 import_solid_model_timer.cpu_secs()) ;  

  entity_list.clear();

#ifdef CUBIT_GUI
#ifndef NO_USAGE_TRACKING
  if(status)
      GUIInterface::instance()->log_import_acis_op();
#endif
#endif // CUBIT_GUI

  return status;
}

CubitStatus AcisQueryEngine::restore_entity_list(
                                      ENTITY_LIST &entity_list,
                                      DLIList<TopologyBridge*> &bridge_list,
                                      bool print_results,
                                      bool step_import,
                                      bool iges_import,
                                      bool import_bodies,
                                      bool import_surfaces,
                                      bool import_curves,
                                      bool import_vertices,
                                      bool free_surfaces ) const
{
  BodySM*  body_ptr;
  Surface* surface_ptr;
  Curve*   curve_ptr;
  Point*   point_ptr;
  int percent_after = 0,
     number_splines_simplified = 0;
  CubitStatus status = CUBIT_SUCCESS;
  outcome result;

    // Initialize the correct position in the list for next().
  entity_list.init();

    // show a progress bar for more than 2 entities
  if( entity_list.count() > 2 )
  {
    char message[128];
    sprintf(message, "Processing %d ACIS Entities", entity_list.count() );
    AppUtil::instance()->progress_tool()->start(0, entity_list.count(),"Progress", 
                                                message, TRUE, TRUE );
  }


  int i;
  for ( i = 0; (i < entity_list.count() && status && !CubitMessage::instance()->Interrupt()); i++ )
  {
    ENTITY *entity_ptr = entity_list[i];
    if (IS_ENTITY_TYPE( entity_ptr, BODY )) 
    {
      if( import_bodies == CUBIT_FALSE )
      {
        api_delent( entity_list[i] );
        if( entity_list.count() > 2 )
            AppUtil::instance()->progress_tool()->step();
        continue;
      }
      
        // Cast the ENTITY to a BODY
      BODY *BODY_ptr = (BODY *)entity_list[i];
      LUMP *body_lump;
      SHELL *body_shell;

      if( (body_lump = BODY_ptr->lump()) &&
          (body_shell = body_lump->shell()) &&
          body_shell->first_face())
      {
        if( body_shell->wire() )
          PRINT_WARNING("Ignoring wire(s) found in ACIS body\n");

        result = api_apply_transf ( BODY_ptr, scale_transf(1.0));
        //if (!result.ok()) 
        //{
            //ACIS_GeometryQueryTool_error (result, "Applying Transform");
        //}
        result = api_change_body_trans ( BODY_ptr, NULL, CUBIT_FALSE );
        //if (!result.ok()) 
        //{
            //ACIS_GeometryQueryTool_error (result, "Changing Body Transform");
        //}
        // Build a VGI Body (i.e, the entire VGI structure) from this 
        // ACIS BODY.
        body_ptr = populate_topology_bridges(BODY_ptr);
        bridge_list.append(body_ptr);
        
        // Make sure that a valid VGI Body was created
        if( !body_ptr )
        {
          // A valid VGI Body was not created...get out of here!
          if (print_results) PRINT_ERROR("Could not process ACIS BODY \n"
            "       Stopped processing the input file.\n");
          status = CUBIT_FAILURE;
          break;
        }
        else
        {
          PRINT_DEBUG_90( "Processed Body %p\n", body_ptr );
        }
      }
     
      //wire body case
      else if( (body_lump = BODY_ptr->lump()) &&
               (body_shell = body_lump->shell()) &&
                body_shell->wire())
      {
        
        ENTITY_LIST edge_list;
        api_get_edges( body_shell, edge_list );
        ENTITY *tmp_ent;
        edge_list.init();
        while((tmp_ent=edge_list.next())!=NULL) 
        {
          EDGE *tmp_edge = static_cast<EDGE*>(tmp_ent); 
          // Copy the EDGE to avoid having free entities that share
          // vertices and curves.  CUBIT currently doesn't like this.
          // This can occur when reading from a file.
          ENTITY* copied_entity_ptr;
          copy_single_entity( tmp_edge, copied_entity_ptr);
          EDGE* copied_EDGE_ptr = (EDGE *)copied_entity_ptr;

          curve_ptr = populate_topology_bridges(copied_EDGE_ptr);
          if( curve_ptr )
            bridge_list.append( curve_ptr );
        } 
        api_delent( BODY_ptr ); 
        if (!result.ok())
        {
          if (print_results) PRINT_ERROR("Problems deleting an ACIS WIRE BODY.\n");
          ACIS_API_error(result);
        }
      }
    }
    else if( IS_ENTITY_TYPE( entity_ptr, FACE ) )
    {
      if( import_surfaces == CUBIT_FALSE )
      {
        api_delent( entity_list[i] );
        if( entity_list.count() > 2 )
          AppUtil::instance()->progress_tool()->step();
        continue;
      }
      //         FACE* FACE_ptr = (FACE *)entity_ptr;
      
      // Copy the FACE to avoid having free entities that share
      // vertices and curves.  CUBIT currently doesn't like this.
      // This can occur when reading from a file.
      ENTITY* copied_entity_ptr;
      copy_single_entity(entity_ptr, copied_entity_ptr);
      FACE* copied_FACE_ptr = (FACE *)copied_entity_ptr;
      ENTITY_LIST entity_list2;
      entity_list2.add(copied_entity_ptr);
      // Delete the original FACE
      result = api_delent(entity_ptr);
      if (!result.ok())
      {
        if (print_results) PRINT_ERROR("Problems deleting an ACIS FACE.\n");
        ACIS_API_error(result);
      }

      if (free_surfaces == CUBIT_FALSE )
      {
        FACE *face_list[1];
        face_list[0] = copied_FACE_ptr;
        BODY *sheet_body = NULL;
        result = api_sheet_from_ff( 1, face_list, sheet_body );
        if (!result.ok() || sheet_body == NULL || sheet_body->lump() == NULL
          || sheet_body->lump()->shell() == NULL ||
          sheet_body->lump()->shell()->first_face() == NULL )
        {
          if (print_results) PRINT_ERROR("Problem with converting face to sheet body.\n");
          continue;
        }
        result =  api_body_to_2d( sheet_body );
        if (!result.ok())
        {
          if (print_results) PRINT_ERROR("Problem with converting face to sheet body.\n");
          continue;
        }
        body_ptr = populate_topology_bridges(sheet_body);
        
        if (body_ptr) 
          bridge_list.append(body_ptr);
      }
      else
      {
        surface_ptr = populate_topology_bridges(copied_FACE_ptr);
        if( surface_ptr )
          bridge_list.append( surface_ptr );
      }
    }
    else if( IS_ENTITY_TYPE( entity_ptr, EDGE) )
    {
      if( import_curves == CUBIT_FALSE )
      {
        api_delent( entity_list[i] );
        if( entity_list.count() > 2 )
          AppUtil::instance()->progress_tool()->step();
        continue;
      }
      // Copy the EDGE to avoid having free entities that share
      // vertices and curves.  CUBIT currently doesn't like this.
      // This can occur when reading from a file.
      ENTITY* copied_entity_ptr;
      copy_single_entity(entity_ptr, copied_entity_ptr);
      EDGE* copied_EDGE_ptr = (EDGE *)copied_entity_ptr;
      
      // Delete the original FACE
      result = api_delent(entity_ptr);
      if (!result.ok())
      {
        if (print_results) PRINT_ERROR("Problems deleting an ACIS EDGE.\n");
        ACIS_API_error(result);
      }
      
      curve_ptr = populate_topology_bridges(copied_EDGE_ptr);
      if( curve_ptr )
        bridge_list.append( curve_ptr );
    }
    else if( IS_ENTITY_TYPE( entity_ptr, CURVE ) )
    {
      if( import_curves == CUBIT_FALSE )
      {
        api_delent( entity_list[i] );
        if( entity_list.count() > 2 )
          AppUtil::instance()->progress_tool()->step();
        continue;
      }
      CURVE* CURVE_ptr = (CURVE*)entity_ptr;
      const curve& curve = CURVE_ptr->equation();
      EDGE* EDGE_ptr = 0;
      outcome result = api_make_edge_from_curve( &curve, EDGE_ptr );
      api_delent( entity_ptr );
      
      if (result.ok() && EDGE_ptr)
      {
        curve_ptr = populate_topology_bridges( EDGE_ptr );
        if (curve_ptr)
          bridge_list.append( curve_ptr );
      }
      else
      {
        ACIS_API_error( result, "Creating EDGE from INTCURVE" );
      }
    }
    else if( IS_ENTITY_TYPE( entity_ptr, VERTEX ) )
    {
      if( import_vertices == CUBIT_FALSE )
      {
        api_delent( entity_list[i] );
        if( entity_list.count() > 2 )
          AppUtil::instance()->progress_tool()->step();
        continue;
      }
      // Copy the VERTEX to avoid having free entities that share
      // vertices and curves.  CUBIT currently doesn't like this.
      // This can occur when reading from a file.
      ENTITY* copied_entity_ptr;
      copy_single_entity(entity_ptr, copied_entity_ptr);
      VERTEX* copied_VERTEX_ptr = (VERTEX *)copied_entity_ptr;
      
      // Delete the original FACE
      result = api_delent(entity_ptr);
      if (!result.ok())
      {
        if (print_results) PRINT_ERROR("Problems deleting an ACIS VERTEX.\n");
        ACIS_API_error(result);
      }
      
      point_ptr = populate_topology_bridges(copied_VERTEX_ptr);
      if( point_ptr )
        bridge_list.append( point_ptr );
    }
    
    if( entity_list.count() > 2 )
      AppUtil::instance()->progress_tool()->step();
  }
  
  // Added by J.Kraftcheck, July 2001.
  // If the above loop aborted for any reason
  // (particularly because cubit_intr was 
  // true), destroy any acis data that was
  // read from the file but has not yet been
  // created in cubit.
  for ( ; i < entity_list.count(); i++ )
    api_delent( entity_list[i] );

  if ( status == CUBIT_SUCCESS )
    if (print_results) PRINT_INFO("\n");

    // end progress reporting
  if( entity_list.count() > 2 )
      AppUtil::instance()->progress_tool()->end();

  return status;
}

//-------------------------------------------------------------------------
// Purpose       : To delete the input ACIS ENTITY_LIST.
//                               
// Special Notes : Re-written by J.Kraftcheck, 2003-12-3.
//
// Creator       : Wes Gill (Tim Tautges)
//
// Creation Date : 5-15-01
//-------------------------------------------------------------------------
void AcisQueryEngine::delete_ACIS_ENTITY(ENTITY_LIST& to_be_deleted) const 
{
  const int num_types = 8;
  const int type_list[num_types] = { VERTEX_TYPE,
                                     EDGE_TYPE,
                                     COEDGE_TYPE, 
                                     LOOP_TYPE, 
                                     FACE_TYPE, 
                                     SHELL_TYPE, 
                                     LUMP_TYPE, 
                                     BODY_TYPE };
                                     
  const int level_list[num_types] = { VERTEX_LEVEL,
                                      EDGE_LEVEL,
                                      COEDGE_LEVEL, 
                                      LOOP_LEVEL, 
                                      FACE_LEVEL, 
                                      SHELL_LEVEL, 
                                      LUMP_LEVEL, 
                                      BODY_LEVEL };
                                      
    // Remove duplicate entities or entities for which the parent entity is in the list.
  DLIList<ENTITY*> delete_list(to_be_deleted.count());
  int i, j, k, m;
  
  for (i = 0; i < num_types; ++i) // for each type
  {
    for (j = 0; j < to_be_deleted.count(); ++j) // for each entity of type
    {
      ENTITY* entity = to_be_deleted[j];
      if (entity->identity(level_list[i]) != type_list[i])
        continue;
      
        // check for duplicates/parents
      bool delete_this_entity = true;
      for (k = 0; delete_this_entity && k < to_be_deleted.count(); ++k)
      {
        if (k == j) 
          continue;
        ENTITY* entity2 = to_be_deleted[k];
        
          // duplicate?
        if (k < j && entity == entity2)
          delete_this_entity = false;
        
          // parent in list?
        for (m = i + 1; m < num_types; ++m)
          if (entity2->identity(level_list[m]) == type_list[m] && is_Related( entity, entity2 ))
            delete_this_entity = false;
      }
      
      if (delete_this_entity)
        delete_list.append(entity);
    }
  }
  
  
  // Delete the Acis ENTITYs in the input ENTITY_LIST
      
  CubitStatus status = CUBIT_FAILURE;
      
    // Now that the list contains no related ENTITYs, delete everything in the list
  while (delete_list.size())
  {
    ENTITY *entity = delete_list.pop();
    status = delete_ACIS_ENTITY(entity);
    if (status == CUBIT_FAILURE)
    {
      PRINT_ERROR("In AGE::delete_ACIS_ENTITY\n"
                      "       Could not delete an ENTITY in the input list.\n"
                      "       Continuing to delete other items in the list. ..\n");
    }  
  }
}

//-------------------------------------------------------------------------
// Purpose       : To delete the input ACIS ENTITY and unhook from the VGI.
//                                              
//
// Special Notes :
//
// Creator       : Wes Gill (Tim Tautges)
//
// Creation Date : 5-15-01
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::delete_ACIS_ENTITY(ENTITY* inputAcisEntity) const 
{
  
  CubitStatus status = CUBIT_SUCCESS;

  if (inputAcisEntity != NULL)
  { 
     BODY *BODYPtr = this->get_BODY_of_ENTITY(inputAcisEntity);   
     if (BODYPtr != NULL)
     {  
       this->delete_ACIS_BODY(BODYPtr, CUBIT_FALSE);   
       return status;
     }
     
     this->unhook_ENTITY_from_VGI(inputAcisEntity);
   
     outcome result = api_delent(inputAcisEntity);
     if (!result.ok())
     {
       PRINT_ERROR("Problems deleting the ACIS ENTITY.\n");  
       ACIS_API_error(result);
       status = CUBIT_FAILURE;
     }
  }
  else
  {
    PRINT_ERROR("BUG: In AGE::delete_ACIS_ENTITY\n"     
                "     Input ENTITY is NULL.\n"
                "     This is a BUG! Please report in.\n"); 
    assert (inputAcisEntity != NULL);
    status = CUBIT_FAILURE;
  }
  return status; 
}

//-------------------------------------------------------------------------
// Purpose       : "Carefully" delete the input list of ACIS BODYs.
//                 If so instructed, also clean out the VGI
//                 datastructures that refer to the ENTITYs in the
//                 input list.
//
// Special Notes :
//
// Creator       : 
//
// Creation Date : 
//-------------------------------------------------------------------------
void AcisQueryEngine::delete_ACIS_BODY(ENTITY_LIST& BODY_list) const 
{
     // Delete the Acis BODYs in the input ENTITY_LIST
   CubitStatus status = CUBIT_FAILURE;
   ENTITY* entity_ptr = NULL;
   BODY_list.init();
   while ( (entity_ptr = BODY_list.next()) != NULL)
   {
        // Delete the Acis BODY "carefully". This procedure also makes
        // sure that no RefEntities are left pointing to any of the Acis
        // ENTITYs in the deleted BODY.
      if (IS_ENTITY_TYPE( entity_ptr, BODY ))
      {
         status = delete_ACIS_BODY ((BODY*)entity_ptr);
         
         if (status == CUBIT_FAILURE)
         {
            PRINT_ERROR("In AGE::delete_ACIS_BODY\n"
                        "       Could not delete a BODY in the input list.\n"
                        "       Continuing to delete other items in the list. ..\n");
         }
         
      }
      
      else
      {
         PRINT_ERROR("In AGE::delete_ACIS_BODY\n"
                     "       One of the ENTITYs in the input list to "
                     "be deleted is not a BODY.\n"
                     "       This item was not deleted.\n"
                     "       Continuing to delete other items in the list...\n");
      }
   }
   
   return;
}


//-------------------------------------------------------------------------
// Purpose       : "Carefully" delete the input ACIS BODY.
//                 If so instructed, also clean out the VGI
//                 datastructures that refer to the ENTITYs in the
//                 input BODY.
//
// Special Notes :
//
// Creator       : 
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::delete_ACIS_BODY(BODY* inputAcisBody,
                                                 CubitBoolean skip_unhook) const
{
    // This function deletes an ACIS BODY and its associated data structures.
    // Also, it removes any references within the Cubit RefEntity datastructure 
    // to any of the entities in inputAcisBody, if any.
  
  CubitStatus status = CUBIT_SUCCESS;
  
    // Now delete the input Acis BODY, inputAcisBody
  if (inputAcisBody != NULL)
  {
      // Remove all links between VGI entities and the ENTITYs in this BODY,
      // if they exist.
    if (!skip_unhook)
      unhook_ENTITY_from_VGI(inputAcisBody) ;
    
    outcome result = api_del_entity(inputAcisBody);
    if (!result.ok())
    {
      PRINT_ERROR("Problems deleting the ACIS BODY.\n");
      ACIS_API_error(result);
      status = CUBIT_FAILURE;
    }
  }
  
  else
  {
    PRINT_ERROR("BUG: In AGE::delete_ACIS_BODY\n"
                "     Input BODY is NULL.\n"
                "  This is a BUG! Please report it.\n");
    assert ( inputAcisBody != NULL );
    status = CUBIT_FAILURE;
  }
  
  return status;
}

// Used by webcut_across_translate
CubitStatus AcisQueryEngine::get_surfs_on_plane( 
                                          BodySM* body_ptr, 
                                          const CubitVector& plane_point, 
                                          const CubitVector& plane_normal,
                                          DLIList<Surface*>& surf_list ) const
{
   DLIList<Surface*> body_surf_list;
   Surface* surf_ptr;
   double face_orig[3], face_norm[3], pln_orig[3], pln_norm[3];
   CubitVector surf_orig, surf_norm;
   plane_point.get_xyz(pln_orig);
   plane_normal.get_xyz(pln_norm);

   body_ptr->surfaces( body_surf_list );

   for( int i=0; i<body_surf_list.size(); i++ )
   {
      surf_ptr = body_surf_list.get_and_step();

      if( surf_ptr->get_point_normal( surf_orig, surf_norm ) == CUBIT_FAILURE )
         continue;

      surf_orig.get_xyz( face_orig ); surf_norm.get_xyz( face_norm );

      if( AnalyticGeometryTool::instance()->is_pln_on_pln( face_orig, face_norm,
                                                           pln_orig, pln_norm ) )
      {
         surf_list.append( surf_ptr );
      }
   }

   return CUBIT_SUCCESS;
}

CubitBoolean AcisQueryEngine::about_spatially_equal (const SPAposition &pos1, 
                                                     const SPAposition &pos2,
                                                     double tolerance_factor) const
{
    // Determines if two positions are within a tolerance value.
  CubitVector v1(pos1.x(), pos1.y(), pos1.z());
  CubitVector v2(pos2.x(), pos2.y(), pos2.z());
  return GeometryQueryTool::instance()->about_spatially_equal(v1, v2, tolerance_factor);
}

CubitBoolean AcisQueryEngine::about_spatially_equal (VERTEX* V1, 
                                                        VERTEX* V2, 
                                                        double tolerance_factor) const
{
    // Determines if two VERTEX's are within a tolerance value.
  const SPAposition &pos1 = V1->geometry()->coords();
  const SPAposition &pos2 = V2->geometry()->coords();
  return about_spatially_equal(pos1, pos2, tolerance_factor);
}

//-------------------------------------------------------------------------
// Purpose       : Unhook the ENTITY's references in VGI. Call the next 
//                 level of the series of functions to allow other 
//                 ENTITYs do their own unhooking.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/12/96
//-------------------------------------------------------------------------
void AcisQueryEngine::unhook_ENTITY_from_VGI(ENTITY *ENTITY_ptr) const
{
   if ( IS_ENTITY_TYPE( ENTITY_ptr, BODY ) )
   {
      unhook_ENTITY_from_VGI((BODY*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, LUMP ) )
   {
      unhook_ENTITY_from_VGI((LUMP*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, SHELL ) )
   {
      unhook_ENTITY_from_VGI((SHELL*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, FACE ) )
   {
      unhook_ENTITY_from_VGI((FACE*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, LOOP ) )
   {
      unhook_ENTITY_from_VGI((LOOP*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, COEDGE ) )
   {
      unhook_ENTITY_from_VGI((COEDGE*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, EDGE ) )
   {
      unhook_ENTITY_from_VGI((EDGE*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, VERTEX ) )
   {
      unhook_ENTITY_from_VGI((VERTEX*)ENTITY_ptr);
   }
   else {
     assert(0);
   }
}

void AcisQueryEngine::unhook_ENTITY_from_VGI(BODY* BODY_ptr) const
{
    // Don't proceed if the input is NULL ;
  if ( BODY_ptr == NULL)
    return;
  
    // First, get the BodyACIS using this BODY
  AcisBridge* bridge = ATTRIB_CUBIT_OWNER::cubit_owner(BODY_ptr);
  
    // If there is no BodyACIS associated with the input BODY, then
    // assume the bridge is being reused by another BODY.  Not all
    // of the sub-entities may have been replaced, so keep checking.
  if (bridge != NULL)
  {
      // Now that we have a bridge, convert it to a BodyACIS.
    BodyACIS* bodyACISPtr = CAST_TO(bridge, BodyACIS);
    
      // Make sure that the pointer to BodyACIS is valid
    if ( bodyACISPtr == NULL )
    {
      PRINT_ERROR("In AcisQueryEngine::unhook_ENTITY_from_VGI()\n"
                  "       The VGI BODY must be associated with a BodyACIS\n");
      assert(bodyACISPtr != NULL) ;
    }
    
      // Remove the links between BodyACIS and BODY
    delete bodyACISPtr;
  }
  
    // Now traverse the LUMPs
  for(LUMP* LUMP_ptr = BODY_ptr->lump();
      LUMP_ptr != NULL;
      LUMP_ptr = LUMP_ptr->next())
  {
    this->unhook_ENTITY_from_VGI(LUMP_ptr) ;
  }
}

void AcisQueryEngine::unhook_ENTITY_from_VGI(LUMP* LUMP_ptr) const
{
    //Don't proceed if the input is NULL ;
  if ( LUMP_ptr == NULL)
    return;
  
    // First get the LumpACIS associated with this LUMP
  AcisBridge* bridge = ATTRIB_CUBIT_OWNER::cubit_owner(LUMP_ptr);
  
    // If there is no TB associated with the input ENTITY, then
    // assume the bridge is being reused by another ENTITY.  Not all
    // of the sub-entities may have been replaced, so keep checking.
  if (bridge != NULL)
  {
    LumpACIS* lumpACISPtr = CAST_TO(bridge, LumpACIS) ;
      // Make sure that the pointer to the LumpACIS is valid
    if ( lumpACISPtr == NULL )
    {
      PRINT_ERROR("In AcisQueryEngine::unhook_ENTITY_from_VGI()\n"
                  "       The RefVolume must be associated with a LumpACIS\n") ;
      assert(lumpACISPtr != NULL) ;
    }
    
      // Remove the links betwee LumpACIS and LUMP
    delete lumpACISPtr;
  }
  
    // Now traverse the SHELLs
  for(SHELL* SHELL_ptr = LUMP_ptr->shell();
      SHELL_ptr != NULL;
      SHELL_ptr = SHELL_ptr->next())
  {
    this->unhook_ENTITY_from_VGI(SHELL_ptr) ;
  }
}

void AcisQueryEngine::unhook_ENTITY_from_VGI(SHELL* SHELL_ptr) const
{
    //Don't proceed if the input is NULL ;
  if ( SHELL_ptr == NULL)
    return;
  
    // First get the ShellACIS associated with this SHELL
  AcisBridge* bridge = ATTRIB_CUBIT_OWNER::cubit_owner(SHELL_ptr);
  
    // If there is no TB associated with the input ENTITY, then
    // assume the bridge is being reused by another ENTITY.  Not all
    // of the sub-entities may have been replaced, so keep checking.
  if (bridge != NULL)
  {
    ShellACIS* shellACISPtr = CAST_TO(bridge, ShellACIS);
    
      // Make sure that the pointer to ShellACIS is valid
    if ( shellACISPtr == NULL )
    {
      PRINT_ERROR("In AcisQueryEngine::unhook_ENTITY_from_VGI()\n"
                  "       The VGI Shell must be associated with a ShellACIS\n") ;
      assert(shellACISPtr != NULL) ;
    }
    
      // Remove the links between ShellACIS and SHELL
    delete shellACISPtr;
  }
  
    // Now traverse the FACEs
  for(FACE* FACE_ptr = SHELL_ptr->first_face();
      FACE_ptr != NULL;
      FACE_ptr = FACE_ptr->next())
  {
    this->unhook_ENTITY_from_VGI(FACE_ptr) ;
  }
}

void AcisQueryEngine::unhook_ENTITY_from_VGI(FACE* FACE_ptr) const
{
    //Don't proceed if the input is NULL ;
  if ( FACE_ptr == NULL)
    return;
  
    // First get the SurfaceACIS associated with this FACE
  AcisBridge* bridge = ATTRIB_CUBIT_OWNER::cubit_owner(FACE_ptr);
  
    // If there is no TB associated with the input ENTITY, then
    // assume the bridge is being reused by another ENTITY.  Not all
    // of the sub-entities may have been replaced, so keep checking.
  if (bridge != NULL)
  {
    SurfaceACIS* surfaceACISPtr = CAST_TO(bridge, SurfaceACIS) ;
    
      // Make sure that the pointer to SurfaceACIS is valid
    if ( surfaceACISPtr == NULL )
    {
      PRINT_ERROR("In AcisQueryEngine::unhook_ENTITY_from_VGI()\n"
                  "       The RefFace must be associated with a SurfaceACIS\n") ;
      assert(surfaceACISPtr != NULL);
    }
    
      // Remove the links between SurfaceACIS and FACE
    delete surfaceACISPtr;
  }
  
    // Once we are through with the FACE itself, Traverse the LOOPs
  for(LOOP* LOOP_ptr = FACE_ptr->loop();
      LOOP_ptr != NULL;
      LOOP_ptr = LOOP_ptr->next())
  {
    this->unhook_ENTITY_from_VGI(LOOP_ptr) ;
  }
}

void AcisQueryEngine::unhook_ENTITY_from_VGI(LOOP* LOOP_ptr) const
{
    //Don't proceed if the input is NULL ;
  if ( LOOP_ptr == NULL)
    return;
  
    // First get the VGI Shell associated with this LOOP
  AcisBridge* bridge = ATTRIB_CUBIT_OWNER::cubit_owner(LOOP_ptr);
  
    // If there is no TB associated with the input ENTITY, then
    // assume the bridge is being reused by another ENTITY.  Not all
    // of the sub-entities may have been replaced, so keep checking.
  if (bridge != NULL)
  {
    LoopACIS* loopACISPtr = CAST_TO(bridge, LoopACIS) ;
    
      // Make sure that the pointer to LoopACIS is valid
    if ( loopACISPtr == NULL )
    {
      PRINT_ERROR("In AcisQueryEngine::unhook_ENTITY_from_VGI()\n"
                  "       The VGI Loop must be associated with a LoopACIS\n") ;
      assert(loopACISPtr != NULL) ;
    }
    
      // Remove the links between LoopACIS and LOOP
    delete loopACISPtr;
  }
  
    // Now traverse the COEDGEs
  COEDGE* COEDGE_ptr = LOOP_ptr->start();
  if (COEDGE_ptr != NULL)
  {
    do
    {
      this->unhook_ENTITY_from_VGI(COEDGE_ptr);
      COEDGE_ptr = COEDGE_ptr->next();
    } while ( COEDGE_ptr != NULL && COEDGE_ptr != LOOP_ptr->start() );
  }
}

void AcisQueryEngine::unhook_ENTITY_from_VGI(COEDGE* COEDGE_ptr) const
{
    //Don't proceed if the input is NULL ;
  if ( COEDGE_ptr == NULL)
    return;
  
    // First get the VGI CoEdge associated with this COEDGE
  AcisBridge* bridge = ATTRIB_CUBIT_OWNER::cubit_owner(COEDGE_ptr);
  
    // If there is no TB associated with the input ENTITY, then
    // assume the bridge is being reused by another ENTITY.  Not all
    // of the sub-entities may have been replaced, so keep checking.
  if (bridge != NULL)
  {
    CoEdgeACIS* coedgeACISPtr = CAST_TO(bridge, CoEdgeACIS) ;
    
      // Make sure that the pointer to CoEdgeACIS is valid
    if ( coedgeACISPtr == NULL )
    {
      PRINT_ERROR("In AcisQueryEngine::unhook_ENTITY_from_VGI()\n"
                  "       The VGI CoEdge must be associated with a CoEdgeACIS\n") ;
      assert(coedgeACISPtr != NULL) ;
    }
    
      // Remove the links between CoEdgeACIS and COEDGE
    delete coedgeACISPtr;
  }
  
    // Once we are through with the COEDGE itself, unhook the EDGE
  this->unhook_ENTITY_from_VGI(COEDGE_ptr->edge()) ;
}

void AcisQueryEngine::unhook_ENTITY_from_VGI(EDGE* EDGE_ptr,
                                                CubitBoolean remove_lower_entities) const
{
    //Don't proceed if the input is NULL ;
  if (EDGE_ptr == NULL)
    return;
  
    // First get the CurveACIS associated with this EDGE
  AcisBridge* bridge = ATTRIB_CUBIT_OWNER::cubit_owner(EDGE_ptr);
  
    // If there is no TB associated with the input ENTITY, then
    // assume the bridge is being reused by another ENTITY.  Not all
    // of the sub-entities may have been replaced, so keep checking.
  if (bridge != NULL)
  {
    CurveACIS* curveACISPtr = CAST_TO(bridge, CurveACIS) ;
    
      // Make sure that the pointer to CurveACIS is valid
    if (curveACISPtr == NULL)
      return;
    
      // Remove the links between CurveACIS and EDGE
    delete curveACISPtr;
  }
  
    // Now traverse the VERTEXes
  if (remove_lower_entities)
  {
      // Once we are through with the EDGE itself, unhook the VERTEXes.
    this->unhook_ENTITY_from_VGI(EDGE_ptr->start()) ;
    this->unhook_ENTITY_from_VGI(EDGE_ptr->end()) ;
  }
}

void AcisQueryEngine::unhook_ENTITY_from_VGI(VERTEX* VERTEX_ptr) const
{
    //Don't proceed if the input is NULL ;
  if (VERTEX_ptr == NULL)
    return;
  
    // First get the PointACIS associated with this VERTEX
  AcisBridge* bridge = ATTRIB_CUBIT_OWNER::cubit_owner(VERTEX_ptr);
  PointACIS* pointACISPtr = CAST_TO(bridge, PointACIS) ;
  
    // Make sure that the pointer to PointACIS is valid
  if (pointACISPtr == NULL)
    return;
  
    // Remove the links between PointACIS and VERTEX
  delete pointACISPtr;
}

//-------------------------------------------------------------------------
// Purpose       : Save the input ENTITY to a SAT file
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 
//-------------------------------------------------------------------------
int AcisQueryEngine::save_ENTITY_as_sat_file ( ENTITY* entity_ptr, 
                                                  const char* filename,
                                                  const char* update_mode ) const
{
   FILE* ENTITY_file_ptr = fopen(filename, update_mode);
   if (ENTITY_file_ptr == NULL)
   {
      PRINT_ERROR ("Cannot open file. ENTITY not saved.\n");
      return CUBIT_FALSE;
   }
   
   else
   {
      CubitString version = "Cubit ";
      CubitString cubit_version("9.1b");
      version += cubit_version;
      FileInfo info;
      info.set_product_id(version.c_str());

      info.set_units(1.0);
      api_set_file_info((FileId | FileUnits), info);
     
      ENTITY_LIST entity_list;
      entity_list.add(entity_ptr);
      outcome result = api_save_entity_list(ENTITY_file_ptr,
                                    CUBIT_TRUE, entity_list);
      fclose(ENTITY_file_ptr);
      if (!result.ok())
      {
         ACIS_API_error(result);
         PRINT_ERROR ("Problem saving ENTITY to sat file.\n");
         entity_list.clear();
         return CUBIT_FALSE;
      }
      entity_list.clear();
      return CUBIT_TRUE;
   }
}

//-------------------------------------------------------------------------
// Purpose       : Return a bounding box that encompasses all the
//                 ENTITYs in the input ENTITY_LIST.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::create_super_acis_bounding_box( 
    ENTITY_LIST& entity_list, SPAbox &super_box ) const
{
     // Make sure there are some ENTITYs in the input list
   if (entity_list.count() == 0)
   {
      PRINT_ERROR ("Cannot generate bounding box.\n"
                   "       The input list of ACIS ENTITYs is empty.\n");
      SPAbox empty_box;
      super_box = empty_box;
      return CUBIT_SUCCESS;
   }
   
   ENTITY* entity_ptr = NULL;
   SPAbox entity_box;
   entity_list.init();
   while ( (entity_ptr = entity_list.next()) != NULL)
   {
        // Get the ACIS bounding box for this ENTITY. Note that this procedure
        // uses the transformation associated with the owning BODY
        // of the entity to transform the bounding box of the entity in
        // its local coordinate system into the global coordinate system
        // (of the owning BODY)
      entity_box = get_acis_entity_bounding_box(entity_ptr);
      
        // "Concatenate" this box with the super_box, creating a bounding
        // box that bounds the ENTITYs (from the list), processed so far.
      super_box = super_box | entity_box;
   }
   
   return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Return a bounding box that encompasses all the
//                 Body's in the input list of Body's.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::create_super_acis_bounding_box( 
    DLIList<BodySM*>& body_list,
    SPAbox &super_box) const
{
     // Make sure there are some Body's in the input list
   if (body_list.size() == 0)
   {
      PRINT_ERROR ("Cannot generate bounding box.\n"
                   "       The input list of Bodies is empty.\n");
      SPAbox empty_box;
      super_box = empty_box;
      return CUBIT_SUCCESS;
   }
   
   BodyACIS* body_ptr = NULL;
   ENTITY* entity_ptr = NULL;
   SPAbox entity_box;
   
   body_list.reset();
   for( int i = 0; i< body_list.size(); i++)
   {
      body_ptr = dynamic_cast<BodyACIS*>(body_list.get_and_step());
      if (!body_ptr)
      {
        PRINT_ERROR("Non-ACIS BodySM at %s:%d\n", __FILE__, __LINE__ );
        return CUBIT_FAILURE;
      }
      
        //get the BODY from the Body
      entity_ptr = body_ptr->get_BODY_ptr();
      if( entity_ptr == NULL ) 
      {
         return CUBIT_FAILURE;
      }
      
        // Get the ACIS bounding box for this ENTITY.
      entity_box = get_acis_entity_bounding_box(entity_ptr);
      
        // "Concatenate" this box with the super_box, creating a bounding
        // SPAbox that bounds the ENTITYs (from the list), processed so far.
      super_box = super_box | entity_box;
   }
   
   return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Generating a CubitBox from an ACIS SPAbox
//
// Special Notes : 
//
// Creator       : Malcolm Panthaki
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------
CubitBox AcisQueryEngine::bounding_box (SPAbox const & a_box)
{
   CubitVector minimum(a_box.low().x(),  a_box.low().y(),  a_box.low().z());
   CubitVector maximum(a_box.high().x(), a_box.high().y(), a_box.high().z());
   return CubitBox(minimum, maximum);
}
SPAbox AcisQueryEngine::bounding_box(const CubitBox& cubit_box)
{
  CubitVector min = cubit_box.minimum();
  CubitVector max = cubit_box.maximum();
  return SPAbox( SPAinterval( min.x(), max.x() ),
                 SPAinterval( min.y(), max.y() ),
                 SPAinterval( min.z(), max.z() ) );
}
CubitBox AcisQueryEngine::bounding_box (BodySM* body) const
{
  BODY* BODY_ptr = get_BODY(body);
  if (!BODY_ptr)
    return CubitBox();
  
  return bounding_box( bounding_box(BODY_ptr) );
}

//-------------------------------------------------------------------------
// Purpose       : Generating ACIS bounding boxes given ACIS ENTITYs
//
// Special Notes : If the ENTITY does not already have a SPAbox associated
//                 with it, it is computed.
//
// Creator       : 
//
// Creation Date : 
//-------------------------------------------------------------------------
SPAbox AcisQueryEngine::bounding_box(BODY* body_ptr) const
{
     // Get the pointer to the SPAtransf object associated with this BODY
   const SPAtransf* body_transf_ptr = (SPAtransf *)NULL;
   TRANSFORM* t_ptr = body_ptr->transform();
   if (t_ptr != NULL) body_transf_ptr = & (t_ptr->transform());
   
   SPAbox *tmp_box_ptr = body_ptr->bound();
   SPAbox temp_box;
   if( tmp_box_ptr == NULL )
   {
     API_BEGIN;
     temp_box = get_body_box(body_ptr, body_transf_ptr);
     API_END;
   }
   else
     temp_box |= *tmp_box_ptr;
   
   return temp_box;
}

SPAbox AcisQueryEngine::bounding_box(LUMP* lump_ptr) const
{
   const SPAtransf* lump_transf_ptr = (SPAtransf *)NULL;
   
     // Get the BODY that owns this LUMP
   BODY* BODY_ptr = this->get_BODY_of_ENTITY(lump_ptr) ;
   
     // If the LUMP is owned by a BODY, get its transform and use it
     // to compute the bounding box. Otherwise, use the default NULL
     // transform.
   if ( BODY_ptr != NULL )
   {
      TRANSFORM* t_ptr = BODY_ptr->transform();
      if (t_ptr != NULL) 
      {
         lump_transf_ptr = & (t_ptr->transform());
      }
   }

   SPAbox *tmp_box_ptr = lump_ptr->bound();
   SPAbox temp_box;
   if( tmp_box_ptr == NULL )
   {
     API_BEGIN;
     temp_box = get_lump_box(lump_ptr, lump_transf_ptr);
     API_END;
   }
   else
     temp_box |= *tmp_box_ptr;

   return temp_box;
}

SPAbox AcisQueryEngine::bounding_box(FACE* face_ptr) const
{
   if( GeometryQueryTool::instance()->get_facet_bbox() )
   {
      SPAposition box_min(0.0, 0.0, 0.0);
      SPAposition box_max(0.0, 0.0, 0.0);
      CubitStatus status = bounding_box_from_facets(face_ptr,
                                                    box_min, box_max );
      if (status == CUBIT_SUCCESS) {
         SPAbox fbbox(box_min,box_max);
         return fbbox;
      }
      // if there was a failure or no facets have been generated
      // yet then use the standard method
   }

   const SPAtransf* face_transf_ptr = (SPAtransf *)NULL;
   
     // Get the BODY that owns this FACE
   BODY* BODY_ptr = this->get_BODY_of_ENTITY(face_ptr) ;
   
     // If the FACE is owned by a BODY, get its transform and use it
     // to compute the bounding box. Otherwise, use the default NULL
     // transform.
   if ( BODY_ptr != NULL )
   {
      TRANSFORM* t_ptr = BODY_ptr->transform();
      if (t_ptr != NULL) 
      {
         face_transf_ptr = & (t_ptr->transform());
      }
   }
   
   SPAbox *tmp_box_ptr = face_ptr->bound();
   SPAbox temp_box;
   if( tmp_box_ptr == NULL )
   {
     API_BEGIN;
     temp_box = get_face_box(face_ptr, face_transf_ptr);
     API_END;
   }
   else
     temp_box |= *tmp_box_ptr;

   return temp_box;
}

SPAbox AcisQueryEngine::bounding_box(EDGE* edge_ptr) const
{
   const SPAtransf* edge_transf_ptr = (SPAtransf *)NULL;
   
     // Get the BODY that owns this EDGE
   BODY* BODY_ptr = this->get_BODY_of_ENTITY(edge_ptr) ;
     // If the EDGE is owned by a BODY, get its transform and use it
     // to compute the bounding box. Otherwise, use the default NULL
     // transform.
   if ( BODY_ptr != NULL )
   {
      TRANSFORM* t_ptr = BODY_ptr->transform();
      if (t_ptr != NULL) 
      {
         edge_transf_ptr = & (t_ptr->transform());
      }
   }

   SPAbox *tmp_box_ptr = edge_ptr->bound();
   SPAbox temp_box;
   if( tmp_box_ptr == NULL )
   {
     API_BEGIN;
     temp_box = get_edge_box(edge_ptr, edge_transf_ptr);
     API_END;
   }
   else
     temp_box |= *tmp_box_ptr;

   return temp_box;
}
//-------------------------------------------------------------------------
// Purpose       : Return the Acis bounding box based on the facetted
//                 representation of the surface
//
// Special Notes : called from bounding_box only if "set facet bbox on"
//                 has been issued previously
//
// Creator       : Steven J. Owen
//
// Creation Date : 12/15/99
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::bounding_box_from_facets(FACE *face_ptr, 
                                                         SPAposition &box_min,
                                                         SPAposition &box_max ) const
{
   int number_triangles, number_points, number_facets;
   GMem gMem;

   CubitStatus status = 
      get_graphics(face_ptr, number_triangles, number_points,
                          number_facets, &gMem, 
                          15, 0);
   if (status != CUBIT_SUCCESS || number_facets == 0) {
      status = CUBIT_FAILURE;
      return status;
   }

   double xmin, ymin, zmin, xmax, ymax, zmax;
   xmin = ymin = zmin = DBL_MAX;
   xmax = ymax = zmax = -DBL_MAX;
   GPoint* plist = gMem.point_list();
   for (int i = 0; i<number_points; i++) {
      if (plist[i].x > xmax) xmax = plist[i].x; 
      if (plist[i].y > ymax) ymax = plist[i].y; 
      if (plist[i].z > zmax) zmax = plist[i].z; 
      if (plist[i].x < xmin) xmin = plist[i].x; 
      if (plist[i].y < ymin) ymin = plist[i].y; 
      if (plist[i].z < zmin) zmin = plist[i].z; 
   }

   // add 10% to the SPAbox

   double xrange = xmax - xmin; 
   double yrange = ymax - ymin; 
   double zrange = zmax - zmin; 
   double range = (xrange > yrange) ? xrange : yrange;
   range = (range > zrange) ? range : zrange;
   xmin -= range*0.1;
   ymin -= range*0.1;
   zmin -= range*0.1;
   xmax += range*0.1;
   ymax += range*0.1;
   zmax += range*0.1;

   box_max.set_x( xmax ); 
   box_max.set_y( ymax ); 
   box_max.set_z( zmax ); 
   box_min.set_x( xmin ); 
   box_min.set_y( ymin ); 
   box_min.set_z( zmin ); 
   
   return CUBIT_SUCCESS;

}

//-------------------------------------------------------------------------
// Purpose       : Return the bounding SPAbox of the input ENTITY
//
// Special Notes : If the bounding_box does not exist, one is 
//                 created for the ENTITY.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 
//-------------------------------------------------------------------------
SPAbox AcisQueryEngine::get_acis_entity_bounding_box(ENTITY* entity_ptr) const
{
     // This procedure currently works only for ACIS BODY, LUMP, FACE and 
     // EDGE ENTITYs
     //
     // MJP NOTE:
     // Note that the get_xxx_box procedures were used because these check
     // to make sure that a bounding box exists, and if not, create one before
     // returning the SPAbox object. The Direct Interface member functions
     // (e.g., lump_ptr->bound() will return NULL if the bounding SPAbox doesn't exist.
     //
     // MJP NOTE:
     // In each of the blocks of code, one of the steps is to get the pointer to
     // the transformation object associated with the owning BODY. Presently,
     // this is done using the Direct Interface (member functions of the
     // various topology entities are called). The API associated with the
     // Geometry Husk of the ACIS 3DToolkit has a function, api_get_owner,
     // which ought to be used instead. The change should be made if the
     // ACIS 3DToolkit is licensed for Cubit.
     //
     // For example, instead of:
     //      SPAtransf* face_transf_ptr = &( ( (face_ptr->shell()->lump()->body())-> 
     //                                        transform() )->
     //                                      transform() );
     // we would have:
     //      ENTITY* owner = NULL;
     //      outcome result = api_get_owner (face_ptr, owner);
     //      if (!result.ok()) 
     //         {...};
     //      else if (owner->identity() == BODY_TYPE)
     //      {
     //         BODY* body_ptr = (BODY *) owner;
     //         SPAtransf* face_transf_ptr = &(body_ptr->transform()->transform());
     //      }
     //   else if (owner->identity() != BODY_TYPE)
     //         {...}
     //
     // Much longer....but cleaner and safer (:-) especially since the owning 
     // ENTITY of any given ACIS ENTITY may not always be a BODY.
   
     // Get the bounding SPAbox of a BODY
   if (IS_ENTITY_TYPE( entity_ptr, BODY ))
   {
      return bounding_box( (BODY *) entity_ptr );
   }
   
     // Get the bounding SPAbox of a LUMP
   else if (IS_ENTITY_TYPE( entity_ptr, LUMP ))
   {
      return bounding_box ( (LUMP *) entity_ptr );
   }
   
     // Get the bounding SPAbox of a FACE
   else if (IS_ENTITY_TYPE( entity_ptr, FACE ))
   {
      return bounding_box( (FACE *) entity_ptr );
   }
   
     // Get the bounding SPAbox of an EDGE
   else if (IS_ENTITY_TYPE( entity_ptr, EDGE ))
   {
      return bounding_box( (EDGE *) entity_ptr );
   }
   
     // The input ENTITY type is not supported yet. Crash and burn :-(
   else
   {
      PRINT_ERROR("Cannot compute the ACIS bounding box.\n"
                  "       Input ENTITY type not supported yet in "
                  "AcisQueryEngine::get_acis_entity_bounding_box\n");
      assert (0);
      return SPAbox();
   }
}


//-------------------------------------------------------------------------
// Purpose       : Get the longest side of the input ACIS bounding SPAbox.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 
//-------------------------------------------------------------------------
double AcisQueryEngine::get_max_size_of_acis_box ( 
    const SPAbox& acis_box ) const
{
     // MJP NOTE:
     // This is not the most efficient implementation, but it is more readable
     // than one that generates fewer temporaries, and I don't forsee this
     // being called umpteen times :-)
   
     // Get the position of the end points of the longest diagonal of the SPAbox
   SPAposition max_box_coords = acis_box.high();
   SPAposition min_box_coords = acis_box.low();
   
   double sizeX = fabs(max_box_coords.x() - min_box_coords.x());
   double sizeY = fabs(max_box_coords.y() - min_box_coords.y());
   double sizeZ = fabs(max_box_coords.z() - min_box_coords.z());
   
   return (  (sizeX > sizeY) ?
             ( (sizeX > sizeZ) ? sizeX : sizeZ)
             : (sizeY > sizeZ) ? sizeY : sizeZ );   
}

void AcisQueryEngine::get_all_cubit_owners(BODY *BODY_ptr,
                             DLIList<TopologyBridge*> &tb_list) const
{
  
    //- get all the te's owned by this entity and its children

  AcisBridge *tb_ptr;
  
  tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(BODY_ptr);
  if (tb_ptr != NULL) tb_list.append(CAST_TO(tb_ptr, TopologyBridge));
  
    // Get the LUMPs in this BODY
  for( LUMP* LUMP_ptr = BODY_ptr->lump();
       LUMP_ptr != NULL; 
       LUMP_ptr = LUMP_ptr->next() )
  {
    tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(LUMP_ptr);
    if (tb_ptr != NULL) {
      tb_list.append(CAST_TO(tb_ptr, TopologyBridge));
    }
        
      // For each LUMP, traverse its SHELLs
    for( SHELL* SHELL_ptr = LUMP_ptr->shell();
         SHELL_ptr != NULL;
         SHELL_ptr = SHELL_ptr->next() )
    {
      tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(SHELL_ptr);
      if (tb_ptr != NULL) tb_list.append(CAST_TO(tb_ptr, TopologyBridge));
      
        // Get the FACEs of this SHELL
      for( FACE* FACE_ptr = SHELL_ptr->first_face();
           FACE_ptr != NULL;
           FACE_ptr = FACE_ptr->next_face() )
      {
        tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(FACE_ptr);
        if (tb_ptr != NULL) {
          tb_list.append(CAST_TO(tb_ptr, TopologyBridge));
        }
        
          // Get the LOOPs of this FACE
        for( LOOP* LOOP_ptr = FACE_ptr->loop();
             LOOP_ptr != NULL;
             LOOP_ptr = LOOP_ptr->next() )
        {
          tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(LOOP_ptr);
          if (tb_ptr != NULL) tb_list.append(CAST_TO(tb_ptr, TopologyBridge));
          
            // Get the COEDGEs in this LOOP
          COEDGE* COEDGE_ptr = LOOP_ptr->start();
          if (COEDGE_ptr == NULL)  // Accounting for possible incomplete 
              // topology
            continue;
          
            // Get its underlying EDGE and search it!
            // MJP Note: 
            // This search algorithm is inefficient in that it will
            // search an EDGE as many times as there are associated
            // COEDGES. But it will suffice for now. The 3DToolkit has
            // a new set of topology routines that will make the
            // process of getting all EDGEs in a BODY easier. We can then
            // save a few microseconds :-)
            // The same inefficiency holds for VERTEXes.
          EDGE* EDGE_ptr = NULL;
          VERTEX* VERTEX_ptr = NULL;
          do
          {
            tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(COEDGE_ptr);
            if (tb_ptr != NULL) tb_list.append(CAST_TO(tb_ptr, TopologyBridge));

              // Get the underlying EDGE
            EDGE_ptr = COEDGE_ptr->edge();

            tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(EDGE_ptr);
            if (tb_ptr != NULL) {
              tb_list.append_unique(CAST_TO(tb_ptr, TopologyBridge));
            }

              // Get the start VERTEX
            VERTEX_ptr = COEDGE_ptr->start();
            tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(VERTEX_ptr);
            if (tb_ptr != NULL) tb_list.append_unique(CAST_TO(tb_ptr, TopologyBridge));

              // Get the end VERTEX
            VERTEX_ptr = COEDGE_ptr->end();
            tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(VERTEX_ptr);
            if (tb_ptr != NULL) tb_list.append_unique(CAST_TO(tb_ptr, TopologyBridge));

              // Move to the next COEDGE in this LOOP
            COEDGE_ptr = COEDGE_ptr->next();

          } while (COEDGE_ptr != NULL && COEDGE_ptr != LOOP_ptr->start());
        }
      }
    }
  }
}

//-------------------------------------------------------------------------
// Purpose       : Remove the ATTRIB_CUBIT_OWNER attributes from all the
//                 ENTITYs of the input BODY.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 
//-------------------------------------------------------------------------
void 
AcisQueryEngine::remove_cubit_owner_attrib_in_BODY(BODY* BODY_ptr) const
{
    // Some checks
  if (BODY_ptr == NULL)
    return;
  
    // Traverse the entire BODY and remove the cubit owner attributes from
    // all the ENTITYs
  
    // Remove the attribute
  ATTRIB_CUBIT_OWNER::remove_cubit_owner(BODY_ptr);
  
    // Get the LUMPs in this BODY
  LUMP* LUMP_ptr = BODY_ptr->lump();
  while (LUMP_ptr != NULL)
  {
      // Remove the attribute
    ATTRIB_CUBIT_OWNER::remove_cubit_owner(LUMP_ptr);
    
      // For each LUMP, traverse its SHELLs
    SHELL* SHELL_ptr = LUMP_ptr->shell();
    while (SHELL_ptr != NULL)
    {
        // Remove the attribute
      ATTRIB_CUBIT_OWNER::remove_cubit_owner(SHELL_ptr);
      
        // Get the FACEs of this SHELL
      FACE* FACE_ptr = SHELL_ptr->first_face();
      while (FACE_ptr != NULL)
      {
          // Remove the attribute
        ATTRIB_CUBIT_OWNER::remove_cubit_owner(FACE_ptr);
        
          // Get the LOOPs of this FACE
        LOOP* LOOP_ptr = FACE_ptr->loop();
        while (LOOP_ptr != NULL)
        {
            // Remove the attribute
          ATTRIB_CUBIT_OWNER::remove_cubit_owner(LOOP_ptr);
          
            // Get the COEDGEs in this LOOP
          COEDGE* COEDGE_ptr = LOOP_ptr->start();
          if (COEDGE_ptr != NULL)  // Accounting for possible incomplete 
              // topology
          {
              // Get its underlying EDGE and search it!
              // MJP Note: 
              // This search algorithm is inefficient in that it will
              // search an EDGE as many times as there are associated
              // COEDGES. But it will suffice for now. The 3DToolkit has
              // a new set of topology routines that will make the
              // process of getting all EDGEs in a BODY easier. We can then
              // save a few microseconds :-)
              // The same inefficiency holds for VERTEXes.
            EDGE* EDGE_ptr = NULL;
            VERTEX* VERTEX_ptr = NULL;
            do
            {
                // Remove the attribute
              ATTRIB_CUBIT_OWNER::remove_cubit_owner(COEDGE_ptr);
              
                // Get the underlying EDGE
              EDGE_ptr = COEDGE_ptr->edge();
              
                // Remove the attribute
              ATTRIB_CUBIT_OWNER::remove_cubit_owner(EDGE_ptr);
              
                // Get the start VERTEX
              VERTEX_ptr = COEDGE_ptr->start();
              
                // Remove the attribute
              ATTRIB_CUBIT_OWNER::remove_cubit_owner(VERTEX_ptr);
              
                // Get the end VERTEX
              VERTEX_ptr = COEDGE_ptr->end();
              
                // Remove the attribute
              ATTRIB_CUBIT_OWNER::remove_cubit_owner(VERTEX_ptr);
              
                // Move to the next COEDGE in this LOOP
              COEDGE_ptr = COEDGE_ptr->next();
              
            } while (COEDGE_ptr != NULL && COEDGE_ptr != LOOP_ptr->start());
          }
          
            // Move to the next LOOP in this FACE
          LOOP_ptr = LOOP_ptr->next();
        }
        
          // Move to the next FACE in this SHELL
        FACE_ptr = FACE_ptr->next_face();
      }
      
        // Move to the next SHELL in this LUMP
      SHELL_ptr = SHELL_ptr->next();
    }
    
      // Move to the next LUMP in this BODY
    LUMP_ptr = LUMP_ptr->next();
  }
}

BODY* AcisQueryEngine::get_BODY_of_ENTITY(BODY* BODY_ptr) const
{
  return BODY_ptr ;
}

BODY* AcisQueryEngine::get_BODY_of_ENTITY(LUMP* LUMP_ptr) const
{
  if (LUMP_ptr == NULL) return (BODY *)NULL;
  return LUMP_ptr->body();
}

BODY* AcisQueryEngine::get_BODY_of_ENTITY(WIRE* WIRE_ptr) const
{
  if (WIRE_ptr == NULL) return (BODY *)NULL;
  return WIRE_ptr->body();
}

BODY* AcisQueryEngine::get_BODY_of_ENTITY(SHELL* SHELL_ptr) const
{
  if (SHELL_ptr == NULL) return (BODY *)NULL;
  return this->get_BODY_of_ENTITY(SHELL_ptr->lump());
}

BODY* AcisQueryEngine::get_BODY_of_ENTITY(FACE* FACE_ptr) const
{
  if (FACE_ptr == NULL) return (BODY *)NULL;
  return this->get_BODY_of_ENTITY(FACE_ptr->shell());
}

BODY* AcisQueryEngine::get_BODY_of_ENTITY(LOOP* LOOP_ptr) const
{
  if (LOOP_ptr == NULL) return (BODY *)NULL;
  return this->get_BODY_of_ENTITY(LOOP_ptr->face());
}

BODY* AcisQueryEngine::get_BODY_of_ENTITY(COEDGE* COEDGE_ptr) const
{
  if (COEDGE_ptr == NULL) return (BODY *)NULL;
    // The owner of a COEDGE could be a LOOP, a WIRE or a SHELL. Find out
    // the right type of owner and then find the BODY that it belongs to.
  if ( COEDGE_ptr->loop() != NULL )
  {
      // The owner's a LOOP
    return this->get_BODY_of_ENTITY(COEDGE_ptr->loop()) ;
  }
  else if ( COEDGE_ptr->shell() != NULL )
  {
      // The owner's a SHELL
    return this->get_BODY_of_ENTITY(COEDGE_ptr->shell()) ;
  }
  else if ( COEDGE_ptr->wire() != NULL )
  {
      // The owner's a WIRE
    return this->get_BODY_of_ENTITY(COEDGE_ptr->wire()) ;
  }
  else
  {
      // The COEDGE is not owned by anybody. Should be an error???
      // Not sure. Return a NULL.
    return NULL ;
  }
}

BODY* AcisQueryEngine::get_BODY_of_ENTITY(EDGE* EDGE_ptr) const
{
  if (EDGE_ptr == NULL) return (BODY *)NULL;
  return this->get_BODY_of_ENTITY(EDGE_ptr->coedge());
}

BODY* AcisQueryEngine::get_BODY_of_ENTITY(VERTEX* VERTEX_ptr) const
{
  if (VERTEX_ptr == NULL) return (BODY *)NULL;
  return this->get_BODY_of_ENTITY(VERTEX_ptr->edge(0));
}

BODY* AcisQueryEngine::get_BODY(BodySM* bodysm_ptr)
{
    // Get the BodyACIS object associated with the input BodySM
  BodyACIS* BodyACIS_ptr = CAST_TO(bodysm_ptr, BodyACIS);
  
  if (BodyACIS_ptr == NULL)
  {
    PRINT_ERROR("In AcisQueryEngine::get_Body\n"
                "       Body is not a BodyACIS.");
    assert(BodyACIS_ptr != NULL);
    return NULL;
  }
  
    // Now get the LUMP associated with the BodyACIS object
  BODY* BODY_ptr = BodyACIS_ptr->get_BODY_ptr();
  
  if (BODY_ptr == NULL)
  {
    PRINT_ERROR("In AcisQueryEngine::get_Body\n"
                "       There is no BODY object associated"
                " with the BodyACIS object of RefVolume.\n");
    assert(BODY_ptr != NULL);
    return NULL;
  }
  
  return BODY_ptr;
}
  

LUMP* AcisQueryEngine::get_LUMP(Lump* lump_ptr)
{
    // Get the LumpACIS object associated with the input RefVolume
  LumpACIS* LumpACIS_ptr = CAST_TO(lump_ptr, LumpACIS);
  
  if (LumpACIS_ptr == NULL)
  {
    PRINT_ERROR("In AcisQueryEngine::get_Lump\n"
                "       Lump is not a LumpACIS.");
    assert(LumpACIS_ptr != NULL);
    return NULL;
  }
  
    // Now get the LUMP associated with the LumpACIS object
  LUMP* LUMP_ptr = LumpACIS_ptr->get_LUMP_ptr();
  
  if (LUMP_ptr == NULL)
  {
    PRINT_ERROR("In AcisQueryEngine::get_Lump\n"
                "       There is no LUMP object associated"
                " with the LumpACIS object of RefVolume.\n");
    assert(LUMP_ptr != NULL);
    return NULL;
  }
  
  return LUMP_ptr;
}
  

BODY* AcisQueryEngine::get_BODY_of_entity(Lump* lump_ptr) const
{
    // Now get the owning BODY of this LUMP and return it
  return get_BODY_of_ENTITY(get_LUMP(lump_ptr));
}

FACE* AcisQueryEngine::get_FACE(Surface* surface_ptr)
{
    // Get the SurfaceACIS object associated with the input Surface
  SurfaceACIS* SurfaceACIS_ptr = CAST_TO(surface_ptr, SurfaceACIS);
  
  if (SurfaceACIS_ptr == NULL)
  {
    PRINT_ERROR("In AcisQueryEngine::get_FACE\n"
                "       Surface is not a SurfaceACIS\n" );
    assert(SurfaceACIS_ptr != NULL);
    return NULL;
  }
  
    // Now get the FACE associated with the SurfaceACIS object
  FACE* FACE_ptr = SurfaceACIS_ptr->get_FACE_ptr();
  
  if (FACE_ptr == NULL)
  {
    PRINT_ERROR("In AcisQueryEngine::get_FACE\n"
                "       There is no FACE object associated"
                " with the SurfaceACIS object.\n");
    assert(FACE_ptr != NULL);
    return NULL;
  }

  return FACE_ptr;
}

BODY* AcisQueryEngine::get_BODY_of_entity(Surface* surface_ptr) const
{
    // Now get the owning BODY of this FACE and return it
  return get_BODY_of_ENTITY(get_FACE(surface_ptr));
}

EDGE* AcisQueryEngine::get_EDGE(Curve* curve_ptr)
{
    // Get the CurveACIS object associated with the input Curve
  CurveACIS* CurveACIS_ptr = CAST_TO(curve_ptr, CurveACIS);
  
  if (CurveACIS_ptr == NULL)
  {
    PRINT_ERROR("In AcisQueryEngine::get_EDGE\n"
                "       Curve is not a CurveACIS\n" );
    assert(CurveACIS_ptr != NULL);
    return NULL;
  }
  
    // Now get the EDGE associated with the CurveACIS object
  EDGE* EDGE_ptr = CurveACIS_ptr->get_EDGE_ptr();
  
  if (EDGE_ptr == NULL)
  {
    PRINT_ERROR("In AcisQueryEngine::get_EDGE\n"
                "       There is no EDGE object associated"
                " with the CurveACIS object.\n");
    assert(EDGE_ptr != NULL);
    return NULL;
  }
  
  return EDGE_ptr;
}

BODY* AcisQueryEngine::get_BODY_of_entity(Curve* curve_ptr) const
{
    // Now get the owning BODY of this FACE and return it
  return get_BODY_of_ENTITY(get_EDGE(curve_ptr));
}

VERTEX* AcisQueryEngine::get_VERTEX(Point* point_ptr)
{
    // Get the PointACIS object associated with the input Point
  PointACIS* PointACIS_ptr = CAST_TO(point_ptr, PointACIS);
  
  if (PointACIS_ptr == NULL)
  {
    PRINT_ERROR("In AcisQueryEngine::get_VERTEX\n"
                "       Point is not a PointACIS\n" );
    assert(PointACIS_ptr != NULL);
    return NULL;
  }
  
    // Now get the EDGE associated with the PointACIS object
  VERTEX* VERTEX_ptr = PointACIS_ptr->get_VERTEX_ptr();
  
  if (VERTEX_ptr == NULL)
  {
    PRINT_ERROR("In AcisQueryEngine::get_VERTEX\n"
                "       There is no VERTEX object associated"
                " with the PointACIS object.\n");
    assert(VERTEX_ptr != NULL);
    return NULL;
  }
  
    // Now get the owning BODY of this FACE and return it
  return VERTEX_ptr;
}

BODY* AcisQueryEngine::get_BODY_of_entity(Point* point_ptr) const
{
    // Now get the owning BODY of this FACE and return it
  return get_BODY_of_ENTITY(get_VERTEX(point_ptr));
}

ENTITY* AcisQueryEngine::get_ENTITY_of_entity( TopologyBridge *entity_ptr ) const
{
  if( BodyACIS *body_ptr = CAST_TO( entity_ptr, BodyACIS ) )
  {
    BODY *BODY_ptr = body_ptr->get_BODY_ptr();
    
    if (BODY_ptr == NULL)
    {
      PRINT_ERROR("BodyACIS without BODY at %s:%d.\n", __FILE__, __LINE__ );
      return NULL;
    }
    
    return (ENTITY *)BODY_ptr;
  }
  else if( LumpACIS *lumpl_ptr = CAST_TO( entity_ptr, LumpACIS ) )
  {
    // Now get the LUMP associated with the LumpACIS object
    LUMP* LUMP_ptr = lumpl_ptr->get_LUMP_ptr();
    
    if (LUMP_ptr == NULL)
    {
      PRINT_ERROR("LumpACIS without LUMP at %s:%d.\n", __FILE__, __LINE__ );
      return NULL;
    }
    
    return (ENTITY *)LUMP_ptr;
    
  }
  else if( SurfaceACIS *surface_ptr = CAST_TO( entity_ptr, SurfaceACIS ) )
  {
    FACE* FACE_ptr = surface_ptr->get_FACE_ptr();
    
    if (FACE_ptr == NULL)
    {
      PRINT_ERROR("SurfaceACIS without FACE at %s:%d.\n", __FILE__, __LINE__ );
      return NULL;
    }
    
    return (ENTITY *)FACE_ptr;
  }
  else if( CurveACIS *curve_ptr = CAST_TO( entity_ptr, CurveACIS ) )
  {
    EDGE* EDGE_ptr = curve_ptr->get_EDGE_ptr();
    
    if (EDGE_ptr == NULL)
    {
      PRINT_ERROR("CurveACIS without EDGE at %s:%d.\n", __FILE__, __LINE__ );
      return NULL;
    }
    
    return (ENTITY *)EDGE_ptr;
  }
  else if( PointACIS *point_ptr = CAST_TO( entity_ptr, PointACIS ) )
  {
    VERTEX* VERTEX_ptr = point_ptr->get_VERTEX_ptr();
    
    if (VERTEX_ptr == NULL)
    {
      PRINT_ERROR("PointACIS without VERTEX at %s:%d.\n", __FILE__, __LINE__ );
      return NULL;
    }
    
    return (ENTITY *)VERTEX_ptr;
  }

  PRINT_ERROR("Non-ACIS TopologyBridge at %s:%d.\n", __FILE__, __LINE__ );
  return NULL;
}

ENTITY* AcisQueryEngine::get_ENTITY_of_entity( BODY* BODY_ptr, 
                                               TopologyBridge* bridge_ptr ) const
{
    // Traverse the input Acis BODY and return the ENTITY who's
    // parent ATTRIBUTE is the same as the input TopologyEntity
  
    // Make sure we got a BODY
  if (BODY_ptr == NULL) 
  {
    PRINT_ERROR("In AcisQueryEngine::"
                "get_ENTITY_of_entity\n"
                "       Input BODY pointer is NULL.\n");
    return NULL;
  }   
  
    // Make sure we got a TopologyEntity
  AcisBridge* tb_ptr = CAST_TO(bridge_ptr, AcisBridge);
  if (tb_ptr == NULL) 
  {
    PRINT_ERROR("In AcisQueryEngine::"
                "get_ENTITY_of_entity\n"
                "       Input TopologyEntity pointer is NULL.\n");
    return NULL;
  }
  
  AcisBridge* attrib_value = NULL;
  
    // Traverse the LUMPs
  LUMP* lump = BODY_ptr->lump();
  while (lump != NULL)
  {
      // Test the LUMP and if its parent attribute matches 
      // refentity_ptr, return the LUMP
    attrib_value = ATTRIB_CUBIT_OWNER::cubit_owner(lump);
    if (attrib_value == tb_ptr)
    {
      return lump;
    }
    
      // For each LUMP, traverse its SHELLs
    SHELL* shell = lump->shell();
    while (shell != NULL)
    {
        // Test the SHELL and if its parent attribute matches 
        // refentity_ptr, return the SHELL
      attrib_value = ATTRIB_CUBIT_OWNER::cubit_owner(shell);
      if (attrib_value == tb_ptr)
      {
        return shell;
      }
      
        // Get the FACEs of this SHELL
      FACE* face = shell->first_face();
      while (face != NULL)
      {
          // Test the FACE and if its parent attribute matches 
          // refentity_ptr, return the FACE
        attrib_value = ATTRIB_CUBIT_OWNER::cubit_owner(face);
        if (attrib_value == tb_ptr)
        {
          return face;
        }
        
          // Get the LOOPs of this FACE
        LOOP* loop = face->loop();
        while (loop != NULL)
        {
            // Test the LOOP and if its parent attribute matches 
            // refentity_ptr, return the LOOP
          attrib_value = ATTRIB_CUBIT_OWNER::cubit_owner(loop);
          
          if (attrib_value == tb_ptr)
          {
            return loop;
          }
          
          COEDGE* coedge = loop->start();       
            // Accounting for possible incomplete topology
          while( coedge != NULL )
          {
              // Test the COEDGE -- if its parent attribute matches
              // refentity_ptr, return the COEDGE
            attrib_value = ATTRIB_CUBIT_OWNER::cubit_owner(coedge);
            if (attrib_value == tb_ptr)
            {
              return coedge;
            }
            
            EDGE* edge = coedge->edge();
            if( edge != NULL )
            {
                // Test the EDGE and if its parent attribute
                // matches refentity_ptr, return the EDGE
              attrib_value = ATTRIB_CUBIT_OWNER::cubit_owner(edge);
              if (attrib_value == tb_ptr)
              {
                return edge;
              }
            }
            
            VERTEX* vertex = coedge->start();
            if( vertex != NULL )
            {
                // Test the VERTEX and if its parent attribute
                // matches refentity_ptr, return the VERTEX
              attrib_value = ATTRIB_CUBIT_OWNER::cubit_owner(vertex);
              if (attrib_value == tb_ptr)
              {
                return vertex;
              }
            }
            
            vertex = coedge->end();
            if( vertex != NULL )
            {
                // Test the VERTEX and if its parent attribute
                // matches refentity_ptr, return the VERTEX
              attrib_value = ATTRIB_CUBIT_OWNER::cubit_owner(vertex);
              if (attrib_value == tb_ptr)
              {
                return vertex;
              }
            }
            
              // Move to the next COEDGE in this LOOP
            coedge = coedge->next();
            if( coedge == loop->start() )
              break;
          }
            // Move to the next LOOP in this FACE
          loop = loop->next();
        }
          // Move to the next FACE in this SHELL
        face = face->next_face();
      }
        // Move to the next SHELL in this LUMP
      shell = shell->next();
    }
      // Move to the next LUMP in this BODY
    lump = lump->next();
  }
  
    // Couldn't find a corresponding parent attribute
  return NULL;   
}


void AcisQueryEngine::acis_debug_body ( BODY *local_body ) const
{
  debug_entity ( (ENTITY *) local_body, stdout);
}

void AcisQueryEngine::acis_debug_face ( FACE *local_face ) const
{
  debug_entity ( (ENTITY *) local_face, stdout);
}

void AcisQueryEngine::acis_debug_edge ( EDGE *local_edge ) const
{
  debug_entity ( (ENTITY *) local_edge, stdout);
}

void AcisQueryEngine::acis_debug_coedge ( COEDGE *local_coedge ) const
{
  debug_entity ( (ENTITY *) local_coedge, stdout);
}

void AcisQueryEngine::acis_debug_loop ( LOOP *local_loop ) const
{
  debug_entity ( (ENTITY *) local_loop, stdout);
}

void AcisQueryEngine::ACIS_API_error ( outcome result , 
                                          char const* /*message*/ ) const
{
  if ( !result.ok() )
  {
    const char* this_error = find_err_ident(result.error_number());
    if (strncmp(this_error, "NON_POS_WIDTH", 100)  &&
        result.error_number() != API_FAILED)
    {
      err_mess_type err_no = result.error_number();
      PRINT_ERROR("ACIS API error number %d \n"
                  "ACIS API message = %s \n", 
                  result.error_number(), find_err_mess (err_no) );
    }
  }
}

//-------------------------------------------------------------------------
// Purpose       : Deletes all solid model entities associated with the 
//                 Bodies in the input list. 
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 11/07/96
//-------------------------------------------------------------------------
void AcisQueryEngine::delete_solid_model_entities(
    DLIList<BodySM*>& BodyList) const 
{
   BodySM* BodyPtr = NULL;
   for (int i = 0; i < BodyList.size(); i++ )
   {
      BodyPtr = BodyList.get_and_step();
      this->delete_solid_model_entities(BodyPtr);
   }
   
   return;
}

CubitStatus AcisQueryEngine::delete_solid_model_entities(
    GeometryEntity* ref_entity_ptr,
    bool remove_lower_entities) const
{
     // Surface
   Surface* ref_face_ptr = CAST_TO(ref_entity_ptr, Surface);
   if (ref_face_ptr != NULL)
   {
      return ( this->delete_solid_model_entities(ref_face_ptr) );      
   }
   
     // Curve
   Curve* ref_edge_ptr = CAST_TO(ref_entity_ptr, Curve);
   if (ref_edge_ptr != NULL)
   {
      return ( this->delete_solid_model_entities(ref_edge_ptr, 
                                                 remove_lower_entities) );      
   }
   
     // Point
   Point* ref_vertex_ptr = CAST_TO(ref_entity_ptr, Point);
   if (ref_vertex_ptr != NULL)
   {
      return ( this->delete_solid_model_entities(ref_vertex_ptr) );      
   }
   
     // Oops!
   PRINT_ERROR("In AcisQueryEngine::delete_solid_model_entities\n"
               "       Can only delete solid model entities underlying \n"
               "RefFaces, RefEdges and RefVertices.\n");
   return CUBIT_FAILURE;
}


CubitStatus AcisQueryEngine::delete_solid_model_entities(
    BodySM* Body_ptr) const
{
  BODY* BODY_ptr = 0;
  BodyACIS* bodyACISPtr = dynamic_cast<BodyACIS*>(Body_ptr);
  if( bodyACISPtr )
    BODY_ptr = bodyACISPtr->get_BODY_ptr();
  
    // If the BODY pointer is non-NULKL, delete it
  if (BODY_ptr != NULL)
  {
       // Do a "careful" deletion of the ACIS BODY and its contents :-)
     return this->delete_ACIS_BODY(BODY_ptr);
  }
  
  else
  {
       // This Body has no ACIS BODY asscociated with it, so return a 
       // successful completion status
     return CUBIT_SUCCESS;
  }
}

CubitStatus AcisQueryEngine::delete_solid_model_entities(
  Surface* surface_ptr ) const
{
  CubitStatus status = CUBIT_SUCCESS;
  
  SurfaceACIS* sa = CAST_TO(surface_ptr, SurfaceACIS);
  if (sa)
  {
    FACE* FACE_ptr = sa->get_FACE_ptr();
    
      // Remove the links between ACIS and Cubit
    unhook_ENTITY_from_VGI(FACE_ptr);
    
      // Delete the FACE and its contained ENTITIES.
    outcome result = api_delent(FACE_ptr);
    if (!result.ok())
    {
      PRINT_ERROR("Problems deleting an ACIS FACE.\n");
      ACIS_API_error(result);
      status = CUBIT_FAILURE;
    }
  }
  return status;
}

CubitStatus AcisQueryEngine::delete_solid_model_entities(
    Curve* curve_ptr ) const
{
  CubitStatus status = CUBIT_SUCCESS;
  
  CurveACIS* ca = CAST_TO(curve_ptr, CurveACIS);
  if (ca)
  {
    EDGE* EDGE_ptr = ca->get_EDGE_ptr();
    
      // Remove the links between ACIS and Cubit
    unhook_ENTITY_from_VGI(EDGE_ptr);
    
      // Delete the EDGE and its contained ENTITIES.
    outcome result = api_delent(EDGE_ptr);
    if (!result.ok())
    {
      PRINT_ERROR("Problems deleting an ACIS EDGE.\n");
      ACIS_API_error(result);
      status = CUBIT_FAILURE;
    }
  }

  return status;
}

CubitStatus AcisQueryEngine::delete_solid_model_entities(
    Point* point_ptr ) const
{
  CubitStatus status = CUBIT_SUCCESS;
   
  PointACIS* pa = CAST_TO(point_ptr, PointACIS);
  if (pa)
  {
    VERTEX* VERTEX_ptr = pa->get_VERTEX_ptr();
    
      // Remove the links between ACIS and Cubit
    unhook_ENTITY_from_VGI(VERTEX_ptr);
    
      // Delete the VERTEX
    outcome result = api_delent(VERTEX_ptr);
    if (!result.ok())
    {
      PRINT_ERROR("Problems deleting an ACIS VERTEX.\n");
      ACIS_API_error(result);
      status = CUBIT_FAILURE;
    }
  }
 
  return status;
}

//-------------------------------------------------------------------------
// Purpose       : This function clears(i.e. sets to NULL) the bounding 
//                 boxes of a BODY and all the other underlying ENTITYs.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 01/06/97
//-------------------------------------------------------------------------
void AcisQueryEngine::clear_bounding_box(BODY* BODYPtr) const 
{
     // First clear the bounding SPAbox of the BODY
   BODYPtr->set_bound(NULL) ;
   
     // Now clear the bounding boxes of its LUMPs
   LUMP* LUMPPtr = BODYPtr->lump() ;
   
   while(LUMPPtr != NULL)
   {
      this->clear_bounding_box(LUMPPtr) ;
      LUMPPtr = LUMPPtr->next() ;
      
        // Make sure that we don't cycle through all the LUMPs
      if ( LUMPPtr == BODYPtr->lump() )
      {
         LUMPPtr = NULL ;
      }
   }
}

//-------------------------------------------------------------------------
// Purpose       : This function clears(i.e. sets to NULL) the bounding 
//                 boxes of a LUMP and all the other underlying ENTITYs.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 01/06/97
//-------------------------------------------------------------------------
void AcisQueryEngine::clear_bounding_box(LUMP* LUMPPtr) const 
{
     // First clear the bounding SPAbox of the LUMP
   LUMPPtr->set_bound(NULL) ;
   
     // Now clear the bounding boxes of its SHELLs
   SHELL* SHELLPtr = LUMPPtr->shell() ;
   
   while(SHELLPtr != NULL)
   {
      this->clear_bounding_box(SHELLPtr) ;
      SHELLPtr = SHELLPtr->next() ;
      
        // Make sure that we don't cycle through all the SHELLs
      if ( SHELLPtr == LUMPPtr->shell() )
      {
         SHELLPtr = NULL ;
      }
   }
}

//-------------------------------------------------------------------------
// Purpose       : This function clears(i.e. sets to NULL) the bounding 
//                 boxes of a SHELL and all the other underlying ENTITYs.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 01/06/97
//-------------------------------------------------------------------------
void AcisQueryEngine::clear_bounding_box(SHELL* SHELLPtr) const 
{
     // First clear the bounding SPAbox of the SHELL
   SHELLPtr->set_bound(NULL) ;
   
     // Now clear the bounding boxes of its FACEs
   FACE* FACEPtr = SHELLPtr->face() ;
   
   while(FACEPtr != NULL)
   {
      this->clear_bounding_box(FACEPtr) ;
      FACEPtr = FACEPtr->next() ;
      
        // Make sure that we don't cycle through all the FACEs
      if ( FACEPtr == SHELLPtr->face() )
      {
         FACEPtr = NULL ;
      }
   }
}

//-------------------------------------------------------------------------
// Purpose       : This function clears(i.e. sets to NULL) the bounding 
//                 boxes of a FACE and all the other underlying ENTITYs.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 01/06/97
//-------------------------------------------------------------------------
void AcisQueryEngine::clear_bounding_box(FACE* FACEPtr) const 
{
     // First clear the bounding SPAbox of the FACE
   FACEPtr->set_bound(NULL) ;
   
     // Now clear the bounding boxes of its LOOPs
   LOOP* LOOPPtr = FACEPtr->loop() ;
   
   while(LOOPPtr != NULL)
   {
      this->clear_bounding_box(LOOPPtr) ;
      LOOPPtr = LOOPPtr->next() ;
      
        // Make sure that we don't cycle through all the LOOPs
      if ( LOOPPtr == FACEPtr->loop() )
      {
         LOOPPtr = NULL ;
      }
   }
}

//-------------------------------------------------------------------------
// Purpose       : This function clears(i.e. sets to NULL) the bounding 
//                 boxes of a LOOP and all the other underlying ENTITYs.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 01/06/97
//-------------------------------------------------------------------------
void AcisQueryEngine::clear_bounding_box(LOOP* LOOPPtr) const 
{
     // First clear the bounding SPAbox of the LOOP
   LOOPPtr->set_bound(NULL) ;
   
     // Now clear the bounding boxes of its COEDGEs
   COEDGE* COEDGEPtr = LOOPPtr->start() ;
   
   while(COEDGEPtr != NULL)
   {
      this->clear_bounding_box(COEDGEPtr->edge()) ;
      COEDGEPtr = COEDGEPtr->next() ;
      
        // Make sure that we don't cycle through all the COEDGEs
      if ( COEDGEPtr == LOOPPtr->start() )
      {
         COEDGEPtr = NULL ;
      }
   }
}

//-------------------------------------------------------------------------
// Purpose       : This function clears(i.e. sets to NULL) the bounding 
//                 boxes of an EDGE.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 01/06/97
//-------------------------------------------------------------------------
void AcisQueryEngine::clear_bounding_box(EDGE* EDGEPtr) const 
{
   if ( EDGEPtr == NULL )
       return;
   
   EDGEPtr->set_bound(NULL) ;
}

//-------------------------------------------------------------------------
// Purpose       : This function queries ACIS for the necessary HOOPS 
//                 information needed in facetting a RefFace.  This 
//                 information is stored and output in gMem.  The 
//                 number of triangles, points and facets are also 
//                 output.
//
// Special Notes : If a FACE is already faceted, just return the existing
//                 facets.
//
// Creator       : Stephen J. Verzi
//
// Creation Date : 03/04/97
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::get_graphics(
  Surface* surface_ptr,
  int& number_triangles, int& number_points, int& facet_list_size,
  GMem* g_mem, unsigned short normal_tolerance, double distance_tolerance,
  double max_edge_length ) const
{
  //  get the FACE
  FACE* face_ptr = get_FACE(surface_ptr);
  if (face_ptr == NULL)
    return CUBIT_FAILURE;
  return get_graphics ( face_ptr, number_triangles, number_points,
                               facet_list_size, g_mem, normal_tolerance,
                               distance_tolerance, max_edge_length );
}

CubitStatus AcisQueryEngine::get_graphics(
  FACE* face_ptr,
  int& number_triangles, int& number_points, int& facet_list_size,
  GMem* g_mem, unsigned short normal_tolerance, double distance_tolerance, 
  double max_edge_length ) const
{
    // Because this may be unnecessarily called twice,
    // say there is one triangle.
  if (!g_mem)
  {
    number_triangles = 1;
    number_points = 3;
    facet_list_size = 4;
    return CUBIT_SUCCESS;
  }

#if CUBIT_ACIS_VERSION >= 800
  // distance tolerance ignore value was changed
  // handle it here for everyone in cubit
  if(distance_tolerance == 0.0)
    distance_tolerance = -1.;
#endif

    // First, tell the mesh manager about g_mem
  facetManager->set_gmem(g_mem);
  
    // Update the default refinement to reflect tolerances
  REFINEMENT* def_ref = NULL;
  outcome result = api_get_default_refinement(def_ref);
  if (!result.ok())
  {
    PRINT_ERROR("Unable to set graphics tolerances.\n"
                "%s may not display properly.\n",
                face_ptr->type_name());
    this->ACIS_API_error (result);
  }
  else
  {
    API_BEGIN;
    def_ref->set_normal_tol(normal_tolerance);
    def_ref->set_surface_tol(distance_tolerance);
    def_ref->set_max_edge_length( max_edge_length );
    API_END;
  }
  
    // Now facet the face (this will add the facets to g_mem).
  result = api_facet_entity ((ENTITY*)face_ptr);
  if (!result.ok())
  {
    PRINT_ERROR("In AcisQueryEngine::get_graphics_facets\n"
                "  Acis unable to generate facets for %s.\n",
                face_ptr->type_name());
    this->ACIS_API_error (result);
    return CUBIT_FAILURE;
  }
  
    // Fill in the return variables and return.
    // These variables really aren't needed, because they
    // are returned as part of g_mem.
  number_points    = g_mem->pointListCount;
  facet_list_size  = g_mem->fListCount;
  number_triangles = facetManager->polygon_count();
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : This function queries ACIS for the edge information
//                 needed in facetting a RefEdge.  This information is
//                 stored and output in g_mem.
//
// Special Notes :
//
// Creator       : Stephen J. Verzi
//
// Creation Date : 03/04/97
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::get_graphics( Curve* curve_ptr,
                                           int& num_points,
                                           GMem* g_mem,
                                           double tolerance ) const
{
   CubitStatus rv;
   EDGE* the_EDGE = get_EDGE(curve_ptr);
     // NULL pointer is OK

   if( tolerance == 0.0 )
   {
     if( SPAresabs >= SPAresfit )
       tolerance = SPAresabs * 2;
     else
       tolerance = SPAresfit;
   }

   rv = this->facet_EDGE(the_EDGE, num_points, g_mem, (double)tolerance);
   if (rv == CUBIT_FAILURE)
      PRINT_ERROR("Unable to facet curve\n");
   return rv;
}

CubitStatus AcisQueryEngine::facet_EDGE(const EDGE* EDGE_ptr, int& num_points,
                                           GMem* g_mem, double tolerance) const
{
     // Why 'resfit'? -- that is what acis test harness uses.
     // Why resfit? -- Seemed to give good results with slightly
     // less points.
   if (tolerance <= 0)
      tolerance = (double)(2.0 * SPAresfit);
   
   if (EDGE_ptr == NULL || EDGE_ptr->geometry() == NULL ) 
   {
      num_points = 0;
      return CUBIT_SUCCESS;
   }
   
   double tLow  = (double) EDGE_ptr->start_param();
   double tHigh = (double) EDGE_ptr->end_param();
   
   if (EDGE_ptr->sense() == REVERSED) 
   {
      tLow  = -tLow;
      tHigh = -tHigh;
   }
   curve const* curvePtr = &(EDGE_ptr->geometry()->equation());

     int current_max = MAX_NUM_CURVE_POINTS;
     int old_max = CUBIT_INT_MAX;
     
     double old_tolerance = tolerance;
     int old_num_points = CUBIT_INT_MAX;     
     SPAposition *curvePts = NULL;
     double *curvePrms = NULL;
     
      // This is for faceting

   do 
   {
     if(old_max != current_max)
     {
       curvePts = new SPAposition[current_max];
       curvePrms = new double[current_max];
       old_max = current_max;
     }
     
     
        // Facet the curve with a maximum of 'MAX_NUM_CURVE_POINTS'. If
        // more points than that are generated, double the tolerance
        // and refacet. 
      outcome result = api_facet_curve(*curvePtr, tLow, tHigh, tolerance,
                                       current_max,
                                       num_points, curvePts, curvePrms);
      if (!result.ok()) 
      {
         ACIS_API_error(result);
         return CUBIT_FAILURE;
      }
      tolerance *= 2.0;
      if(tolerance > 500 && old_num_points - num_points == 0)
      {
          //Acis is refusing to lessen the number of points to facet
          //the curve, we probably have a hideous curve
        tolerance = old_tolerance;
        current_max *= 2;
      }
      
      old_num_points = num_points;
      
        // NOTE: There seems to be a bug in ACIS 1.7. It is supposed to
        // return num_points = MAX_NUM_CURVE_POINTS+1 if there are too
        // many facet points, but it seems to return
        // num_points == MAX_NUM_CURVE_POINTS instead.
   } while (num_points >= old_max);
   
     // Make sure there is enough space to store points
   g_mem->allocate_polylines(num_points-1);
   
     // Copy the facets into the array provided by the calling code
   for (int i = num_points; i--; ) 
   {
      g_mem->point_list()[i].x = (float)curvePts[i].x();
      g_mem->point_list()[i].y = (float)curvePts[i].y();
      g_mem->point_list()[i].z = (float)curvePts[i].z();
   }
   g_mem->pointListCount = num_points;

   delete [] curvePts;
   delete [] curvePrms;
 
   return CUBIT_SUCCESS;
}


CubitStatus AcisQueryEngine::get_isoparametric_points(
   Surface* surface_ptr, int &nu, int &nv,
   GMem *&g_mem) const
{
   SPAtransf ftrans;
   ENTITY_LIST edge_list;
   outcome rc;
   int i, num_points;
   
     // Get the FACE we'll be using
   FACE* facePtr = get_FACE(surface_ptr);
   if (facePtr == NULL)
   {
      nu = nv = 0;
      return CUBIT_FAILURE;
   }
     // Ask Acis to make some iso curves
   rc = api_face_nu_nv_isolines (nu, nv, facePtr, ftrans, &edge_list);
   if (!rc.ok())
   {
      this->ACIS_API_error (rc);
      nu = nv = 0;
      return CUBIT_FAILURE;
   }

     // Now facet each curve
   nu = nv = edge_list.count();
   g_mem = new GMem[nu];
   for (i = nu; i--; )
   {
        // Facet curve
      this->facet_EDGE((EDGE*)(edge_list[i]), num_points, &(g_mem[i]));
      
        // Delete curve
      rc = api_delent(edge_list[i]);
      if (!rc.ok())
      {
         PRINT_ERROR ("Unable to delete temporary curve.\n");
      }
   }
   edge_list.clear();
   
   return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::get_u_isoparametric_points(
   Surface* surface_ptr, double v, int &n,
   GMem *&g_mem) const
{
   SPAtransf ftrans;
   ENTITY_LIST edge_list;
   outcome rc;
   int i, num_points, nu;
   
     // Get the FACE we'll be using
   FACE* facePtr = get_FACE(surface_ptr);
   if (facePtr == NULL)
   {
      return CUBIT_FAILURE;
   }
     // Ask Acis to make some iso curves
   rc = api_face_u_iso(v, facePtr, ftrans, &edge_list);
   if (!rc.ok())
   {
      this->ACIS_API_error (rc);
      return CUBIT_FAILURE;
   }

     // Now facet each curve
   nu = n = edge_list.count();
   g_mem = new GMem[nu];
   for (i = nu; i--; )
   {
        // Facet curve
      this->facet_EDGE((EDGE*)(edge_list[i]), num_points, &(g_mem[i]));
      
        // Delete curve
      rc = api_delent(edge_list[i]);
      if (!rc.ok())
      {
         PRINT_ERROR ("Unable to delete temporary curve.\n");
      }
   }
   edge_list.clear();
   
   return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::get_v_isoparametric_points(
   Surface* surface_ptr, double u, int& n,
   GMem *&g_mem) const
{
   SPAtransf ftrans;
   ENTITY_LIST edge_list;
   outcome rc;
   int i, num_points, nv;
   
     // Get the FACE we'll be using
   FACE* facePtr = get_FACE(surface_ptr);
   if (facePtr == NULL)
   {
      return CUBIT_FAILURE;
   }
     // Ask Acis to make some iso curves
   rc = api_face_v_iso(u, facePtr, ftrans, &edge_list);
   if (!rc.ok())
   {
      this->ACIS_API_error (rc);
      return CUBIT_FAILURE;
   }

     // Now facet each curve
   nv = n = edge_list.count();
   g_mem = new GMem[nv];
   for (i = nv; i--; )
   {
        // Facet curve
      this->facet_EDGE((EDGE*)(edge_list[i]), num_points, &(g_mem[i]));
      
        // Delete curve
      rc = api_delent(edge_list[i]);
      if (!rc.ok())
      {
         PRINT_ERROR ("Unable to delete temporary curve.\n");
      }
   }
   edge_list.clear();
   
   return CUBIT_SUCCESS;
}

CubitStatus
AcisQueryEngine::transform_vec_position(CubitVector const& position_vector,
                                           BodySM *OSME_ptr,
                                           CubitVector &transform_vector) const
{
     // Get the BodyACIS part of the OSMEPtr
   BodyACIS* body_ACIS_ptr = CAST_TO(OSME_ptr, BodyACIS) ;
   BODY* body_acis = body_ACIS_ptr->get_BODY_ptr();
   
     //Get the transformation SPAmatrix of the body.
   SPAtransf last_tranf = body_acis->transform()->transform();
   
     //-Apply transformation to new node
   SPAposition old_position, new_position;
   
   old_position.set_x( position_vector.x() );
   old_position.set_y( position_vector.y() );
   old_position.set_z( position_vector.z() );
   
   new_position = old_position * last_tranf;
   
   transform_vector.x( new_position.x());
   transform_vector.y( new_position.y());
   transform_vector.z( new_position.z());
   
   return CUBIT_SUCCESS;
}

AcisQueryEngine::AcisQueryEngine()
  : stepInitialized(false)
{
  assert( !instance_ );

    // add this engine to geometryquerytool
  GeometryQueryTool::instance()->add_gqe(this);

  // make a default configuration to initialize base with.
   base_configuration base_config;
#ifdef NT
   // for windows, we'll make a couple modifications
   base_config.enable_freelists = FALSE;
   base_config.enable_audit_leaks = FALSE;
   base_config.enable_audit_logs = FALSE;
#endif
   initialize_base(&base_config);

  
   outcome error = api_start_modeller ( CUBIT_FALSE);
   if ( !error.ok() )
   {
      PRINT_ERROR("Problems initializing ACIS...\n");
      exit (1);
   }
//   API_BEGIN;
     //initialize the correct modules
   api_initialize_generic_attributes();
   api_initialize_spline();
   api_initialize_kernel();
   api_initialize_intersectors();
   api_initialize_constructors();
   api_initialize_euler_ops();

#ifdef ACIS_IGES_TRANSLATOR
     error = api_initialize_xiges();
#endif //end IGES TRANSLATOR


#ifdef ACIS_STEP_TRANSLATOR
     
   outcome my_result = api_initialize_xstep();
   ACIS_API_error(my_result);
   stepInitialized = true;

   // Convert assembly parts to separate bodies, as opposed as 
   // one body with multiple lumps.
   api_set_int_option( "stp_assemb_part_as_body", 1 );

#endif /* ACIS_STEP_TRANSLATOR */

#ifdef ACIS_PROE_TRANSLATOR
   api_initialize_proe();
#endif
#ifdef ACIS_CATIA_TRANSLATOR
   api_initialize_catia();
#endif
   
     // start the faceter if we want graphics
     // Start the ACIS faceter. 
   facetManager = new AcisFacetManager;

#if CUBIT_ACIS_VERSION < 1505
   api_initialise_faceter(facetManager);
#else
   api_set_mesh_manager( facetManager );
#endif
   
     // force tight bounding boxes for spheres and tori
     // this can be critical for feature consolidation - RRL
   option_header *this_option;
   
   this_option = find_option ( "tight_sphere_box" );
   this_option -> set ( TRUE );
   
   this_option = find_option ( "tight_torus_box" );
   this_option -> set ( TRUE );

#ifdef ACIS_IGES_TRANSLATOR
   api_xiges_Ig2Ac_ConvFreeCurves_Set( TRUE );

   //reads in file based on units in file....prevents scaling
   api_xacis_set_application_unit( spaxUnit_unknown );
#endif

   this_option = find_option ( "compress_bb" );
   this_option -> set ( FALSE );

	// We can have a different refinement for each body passed in.  For
	// example, we may want to adjust the refinement according to the body
	// size.  But for demonstration, we will just use the same one for all.
	REFINEMENT *refinement_ptr;
	// Create a refinement.
	outcome result = api_create_refinement(refinement_ptr);
    if (!result.ok())
    {
        ACIS_API_error(result);
    }
	// Set the refinement to obtain:
	// 1. Initial grid in the interior.  Intersection between grid lines and
	//    model edges will not cause points to inserted into the edge.
	// 2. Triangulate everywhere.
	// 3. Smooth all the triangles.
	refinement_ptr->set_triang_mode(AF_TRIANG_ALL);
	refinement_ptr->set_grid_mode(AF_GRID_INTERIOR);
	refinement_ptr->set_adjust_mode(AF_ADJUST_ALL);

    // Set the 'relative' tolerances - it is independent of model size.
    refinement_ptr->set_normal_tol(15);// Angle between two adjacent node normals
    // Set the 'absolute' tolerances - it may depend on model size.
    refinement_ptr->set_surface_tol(0.1);//Distance between facet and surface
    //R->set_max_edge_length(hmax);

    result = api_set_default_refinement(refinement_ptr);
    if (!result.ok())
    {
        ACIS_API_error(result);
    }

     // The way graphics are done now, we don't want to
     // store facets on a face to be retrieved later.
   api_mark_faceted_faces(FALSE);
// Done with stuff from faceting
   
//   API_END;
  api_logging( false );
	// This fixes a geometry file corruption problem
	// in api_save_entity_list()
    int version;
    version = get_allint_version();
    set_export_allint_version(version);

#ifdef CUBIT_LINUX
  error_harden();
#endif

  exportVersion = CUBIT_INT_MAX;

  AcisFeatureEngine::instance();
}

BodySM *AcisQueryEngine::populate_topology_bridges(BODY *body_ptr) const
{
  assert ( body_ptr != NULL);
  
  BodyACIS* bodyACISPtr = NULL;
  
    // First check to make sure we haven't already created a BodyACIS
    // from this ACIS BODY.
  AcisBridge* acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(body_ptr);
  bodyACISPtr = CAST_TO(acis_bridge, BodyACIS);
  if (bodyACISPtr == NULL)
  {
    assert(acis_bridge == NULL);
    bodyACISPtr = new BodyACIS(body_ptr) ;
  }
  else
    bodyACISPtr->set_BODY_ptr(body_ptr);
  
  assert(bodyACISPtr != NULL);
  
  ENTITY_LIST temp_list;
  get_ENTITY_from_ENTITY(LUMP_TYPE, body_ptr, temp_list);
  
  for (int i = 0; i < temp_list.count(); i++)
  {
    LUMP *this_lump = (LUMP *) temp_list[i];
    populate_topology_bridges(this_lump);
  }
  
  return bodyACISPtr;
}

Lump *AcisQueryEngine::populate_topology_bridges(LUMP *lump_ptr) const
{
  assert ( lump_ptr != NULL);
  
  LumpACIS* lumpACISPtr = NULL;
   
    // First check to make sure we haven't already created a LumpACIS
    // from this ACIS LUMP.
  AcisBridge* acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(lump_ptr);
  lumpACISPtr = CAST_TO(acis_bridge, LumpACIS);
  if (lumpACISPtr == NULL)
  {
    assert(acis_bridge == NULL);
    lumpACISPtr = new LumpACIS(lump_ptr) ;
  }
  else
    lumpACISPtr->set_LUMP_ptr(lump_ptr);
  
  assert(lumpACISPtr != NULL);

  ENTITY_LIST temp_list;
  get_ENTITY_from_ENTITY(SHELL_TYPE, lump_ptr, temp_list);

  int i;
  for (i = 0; i < temp_list.count(); i++)
  {
    SHELL *this_shell = (SHELL *) temp_list[i];

      // First check to make sure we haven't already created a ShellACIS
      // from this ACIS SHELL.
    acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(this_shell);
    ShellACIS *shellACISPtr = CAST_TO(acis_bridge, ShellACIS);
    if (shellACISPtr == NULL) {
      assert(acis_bridge == NULL);
      shellACISPtr = new ShellACIS(this_shell) ;
    }
    else
      shellACISPtr->set_SHELL_ptr(this_shell);
  }

  temp_list.clear();
  get_ENTITY_from_ENTITY(FACE_TYPE, lump_ptr, temp_list);

  for (i = 0; i < temp_list.count(); i++) {
    FACE *this_face = (FACE *) temp_list[i];
    populate_topology_bridges(this_face);
  }
  
  return lumpACISPtr;
}

Surface *AcisQueryEngine::populate_topology_bridges(FACE *face_ptr) const
{
  assert ( face_ptr != NULL);
  
  SurfaceACIS* surfaceACISPtr = NULL;
   
    // First check to make sure we haven't already created a SurfaceACIS
    // from this ACIS FACE.
  AcisBridge* acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(face_ptr);
  surfaceACISPtr = CAST_TO(acis_bridge, SurfaceACIS);
  if (surfaceACISPtr == NULL) {
    assert(acis_bridge == NULL);
    surfaceACISPtr = new SurfaceACIS(face_ptr) ;
  }
  else
    surfaceACISPtr->set_FACE_ptr(face_ptr);
  
  assert(surfaceACISPtr != NULL);

    // now loops
  ENTITY_LIST temp_list;
  get_ENTITY_from_ENTITY(LOOP_TYPE, face_ptr, temp_list);

  int i;
  for (i = 0; i < temp_list.count(); i++) {
    LOOP *this_loop = (LOOP *) temp_list[i];
  
      // First check to make sure we haven't already created a LoopACIS
      // from this ACIS LOOP.
    acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(this_loop);
    LoopACIS *loopACISPtr = CAST_TO(acis_bridge, LoopACIS);
    if (loopACISPtr == NULL) {
      assert(acis_bridge == NULL);
      loopACISPtr = new LoopACIS(this_loop) ;
    }
    else
      loopACISPtr->set_LOOP_ptr(this_loop);
  }

    // now coedges
  temp_list.clear();
  get_ENTITY_from_ENTITY(COEDGE_TYPE, face_ptr, temp_list);

  for (i = 0; i < temp_list.count(); i++) {
    COEDGE *this_coedge = (COEDGE *) temp_list[i];

      // First check to make sure we haven't already created a CoEdgeACIS
      // from this ACIS COEDGE.
    acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(this_coedge);
    CoEdgeACIS *coEdgeACISPtr = CAST_TO(acis_bridge, CoEdgeACIS);
    if (coEdgeACISPtr == NULL) {
      assert(acis_bridge == NULL);
      coEdgeACISPtr = new CoEdgeACIS(this_coedge) ;
    }
    else
      coEdgeACISPtr->set_COEDGE_ptr(this_coedge);
  }

  temp_list.clear();
  get_ENTITY_from_ENTITY(EDGE_TYPE, face_ptr, temp_list);

  for (i = 0; i < temp_list.count(); i++) {
    EDGE *this_edge = (EDGE *) temp_list[i];
    populate_topology_bridges(this_edge);
  }
  
  return surfaceACISPtr;
}

Curve *AcisQueryEngine::populate_topology_bridges(EDGE *edge_ptr) const
{
  assert ( edge_ptr != NULL);
  
  CurveACIS* curveACISPtr = NULL;
   
    // First check to make sure we haven't already created a CurveACIS
    // from this ACIS EDGE.
  AcisBridge* acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(edge_ptr);
  curveACISPtr = CAST_TO(acis_bridge, CurveACIS);
  if (curveACISPtr == NULL) {
    assert(acis_bridge == NULL);
    curveACISPtr = new CurveACIS(edge_ptr) ;
  }
  else
    curveACISPtr->set_EDGE_ptr(edge_ptr);
  
  assert(curveACISPtr != NULL);

  populate_topology_bridges(edge_ptr->start());
  populate_topology_bridges(edge_ptr->end());
  
  return curveACISPtr;
}

Point *AcisQueryEngine::populate_topology_bridges(VERTEX *vertex_ptr) const
{
  assert ( vertex_ptr != NULL);
  
  PointACIS* pointACISPtr = NULL;
   
    // First check to make sure we haven't already created a PointACIS
    // from this ACIS VERTEX.
  AcisBridge* acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner(vertex_ptr);
  pointACISPtr = CAST_TO(acis_bridge, PointACIS);
  if (pointACISPtr == NULL) {
    assert(acis_bridge == NULL);
    pointACISPtr = new PointACIS(vertex_ptr) ;
  }
  else
    pointACISPtr->set_VERTEX_ptr(vertex_ptr);
  
  assert(pointACISPtr != NULL);
  return pointACISPtr;
}

BODY* AcisQueryEngine::get_BODY_of_ENTITY(ENTITY* ENTITY_ptr) const
{
   if (ENTITY_ptr == NULL) return (BODY *)NULL;
   
     // Get the BODY that this ENTITY belongs to
   BODY* BODY_ptr = NULL ;
   
   if ( IS_ENTITY_TYPE( ENTITY_ptr, BODY ) )
   {
      BODY_ptr = (BODY*)ENTITY_ptr ;
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, LUMP ) )
   {
      BODY_ptr = ((LUMP*)ENTITY_ptr)->body() ;
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, WIRE ) )
   {
      BODY_ptr = ((WIRE*)ENTITY_ptr)->body() ;
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, SHELL ) )
   {
      BODY_ptr = this->get_BODY_of_ENTITY((SHELL*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, FACE ) )
   {
      BODY_ptr = this->get_BODY_of_ENTITY((FACE*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, LOOP ) )
   {
      BODY_ptr = this->get_BODY_of_ENTITY((LOOP*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, COEDGE ) )
   {
      BODY_ptr = this->get_BODY_of_ENTITY((COEDGE*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, EDGE ) )
   {
      BODY_ptr = this->get_BODY_of_ENTITY((EDGE*)ENTITY_ptr);
   }
   else if (IS_ENTITY_TYPE( ENTITY_ptr, VERTEX ) )
   {
      BODY_ptr = this->get_BODY_of_ENTITY((VERTEX*)ENTITY_ptr);
   }
   else
   {
      PRINT_WARNING("In AcisQueryEngine::get_BODY_of_ENTITY\n"
                    "    Unforeseen ENTITY type in argument : \"%s\"\n",
                    ENTITY_ptr->type_name() ? ENTITY_ptr->type_name() : "(null)") ;
      return NULL;
   }

   return BODY_ptr;
}
   
BodySM* AcisQueryEngine::get_body_sm_of_ENTITY(ENTITY* ENTITY_ptr) const
{

  BODY *BODY_ptr = get_BODY_of_ENTITY(ENTITY_ptr);
  
     // Find the Body associated with this BODY
   if (BODY_ptr != NULL)
   {
      AcisBridge *entity = ATTRIB_CUBIT_OWNER::cubit_owner(BODY_ptr);
      BodySM *body = CAST_TO( entity, BodySM );
      return body;
   }
   
   else
   {
      return NULL;
   }
}

CubitStatus 
AcisQueryEngine::get_EDGEs_of_Curves(DLIList<Curve*>& curve_list, 
                                     DLIList<EDGE*>& EDGE_list) const
{
   int i;
   EDGE *EDGE_ptr;
   Curve *curve_ptr;

   curve_list.reset();
   for( i=0; i<curve_list.size(); i++ )
   {
      curve_ptr = curve_list.get_and_step();
      EDGE_ptr = get_EDGE( curve_ptr );
      if( EDGE_ptr == NULL )
      {
         PRINT_ERROR( "Unable to retrieve ACIS EDGE from Curve\n" );
         return CUBIT_FAILURE;
      }

      EDGE_list.append( EDGE_ptr );
   }

   return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::fire_ray( BodySM *body,
                                       const CubitVector &ray_point,
                                       const CubitVector &unit,
                                       DLIList<double>& ray_params,
                                       DLIList<GeometryEntity*> *entity_list) const
{
   ray_params.clean_out();
     
     // fire a ray at the specified body, returning the entities hit and
     // the parameters along the ray; return non-zero if error
   
     // since ACIS only returns a given entity once, even if it is hit
     // twice (goofy, huh?), we need to fire two rays, one from each direction,
     // and possibly more if the hits don't match
   const double SCALE_FACTOR = 1.05;
   
   SPAposition pos_low(ray_point.x(), ray_point.y(), ray_point.z());
   SPAunit_vector vector_low(unit.x(), unit.y(), unit.z());
   
     // use SPAresabs to get the surface when possible
   double ray_radius = GEOMETRY_RESABS;
   int hits_wanted = 0;
   
   ENTITY_LIST entities_low;
   double *params_low = NULL;
   outcome rc;
   GeometryEntity *geom_entity;
   BODY *body_ACIS =  get_BODY( body );
   
   
   rc = api_ray_test_body (pos_low, vector_low, ray_radius, hits_wanted,
                           body_ACIS, entities_low, params_low);
   
   
   if (!rc.ok()) {
      PRINT_ERROR("First ray fire failed for body\n");
      ACIS_API_error(rc, "Ray fire");
      entities_low.clear();
      return CUBIT_FAILURE;
   }
   else if (!entities_low.count()) 
   {
     entities_low.clear();
     return CUBIT_SUCCESS;
   }
     // get new position that's totally outside bounding SPAbox
   CubitBox box = bounding_box(bounding_box(get_BODY(body)));
     // use max range of bounding box - seems like it should be easier
     // than this!
   CubitVector ray_point2 = ray_point + SCALE_FACTOR * unit * 
       CUBIT_MAX_4(box.x_range(), box.y_range(), box.z_range(), 0);
   SPAposition pos_high(ray_point2.x(), ray_point2.y(), ray_point2.z());
   
   ENTITY_LIST entities_high;
   double *params_high = NULL;
   SPAunit_vector vector_high = -vector_low;
   hits_wanted = 0;
   
   rc = api_ray_test_body (pos_high, vector_high, ray_radius, hits_wanted,
                           body_ACIS, entities_high, params_high);
   
   
   if (!rc.ok()) {
      PRINT_ERROR("Second ray fire failed for body\n");
      ACIS_API_error(rc, "Ray fire");
      entities_low.clear();
      entities_high.clear();
      return CUBIT_FAILURE;
   }
     // now compare entities on list, marching backward on one list and
     // forward on the other; if a mismatch is found, need to fire another ray
   
   if (DEBUG_FLAG(1)) {
      CubitVector dum = ray_point;
      GfxDebug::draw_vector(dum, ray_point2, CUBIT_YELLOW);
      GfxDebug::flush();
   }
   
   int ilow = 0;
   int ihigh = entities_high.count()-1;
   ENTITY_LIST entities_hit;

   double range_orig = (pos_high - pos_low).len();
   double range = range_orig;
     //  RefEntity *ref_entity;
   double res_high = 5.0 * GEOMETRY_RESABS;
   double offset = 0.0;
   
     // step up the low list and down the high list; entity on low list is
     // always added; if it is the first hit on a repeated entity, it will 
     // not match the entity on the high list; if this is the case, the low
     // ray will be fired again and the low list counter reset; continue until
     // we are at the beginning of the high list
   TopologyEntity* topo_entity_ptr;
   
   while (ihigh >= 0) {
      ray_params.append(params_low[ilow] + offset);
      
        // would have done this later, after while loop, if ACIS lists
        // were not brain dead (no duplicates allowed!!!)
      if (entity_list)
      {
         topo_entity_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity( entities_low[ilow] );
         geom_entity = CAST_TO( topo_entity_ptr, GeometryEntity );
         if (geom_entity) entity_list->append(geom_entity);
      }
      
        // do both an entity check and a SPAposition check, in case a single
        // entity is repeated (like the side surf of a cylinder)
      if (entities_low[ilow] == entities_high[ihigh] &&
          fabs(params_low[ilow] - (range - params_high[ihigh])) < res_high)
      {
         ilow++;
         ihigh--;
      }
      
      else {
         if (DEBUG_FLAG(55)) {
              // fire another low ray; add in 0.1% of the distance between high
              // and low hit parameters to make sure pos_low isn't on the surf
            double factor = .001 * (range - params_high[ihigh] - params_low[ilow]);
            offset += params_low[ilow] + factor;
            pos_low += (params_low[ilow] + factor) * vector_low; 
            entities_low.clear();
            hits_wanted = 0;
            delete [] params_low;
            rc = api_ray_test_body (pos_low, vector_low, ray_radius, hits_wanted,
                                    body_ACIS, entities_low, params_low);
            if (!rc.ok()) {
               PRINT_ERROR("Supplemental ray fire failed for body\n");
               ACIS_API_error(rc, "Ray fire");
               
                 // finish off params array with last hit on body
               ray_params.append(range - params_high[0]);
               if (entity_list) {
                  topo_entity_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity( entities_low[ilow] );
                  geom_entity = CAST_TO( topo_entity_ptr, GeometryEntity );
                  if (geom_entity) entity_list->append(geom_entity);
               }
               if (ray_params.size() % 2 == 1) {
                  ray_params.last();
                  ray_params.append(ray_params.get());
                  if (entity_list) {
                     entity_list->last();
                     entity_list->append(entity_list->get());
                  }
               }
               return CUBIT_FAILURE;
            }
            else if (!entities_low.count()) {
               PRINT_ERROR("Imbalanced rays in ray_fire!!!\n");
               entities_low.clear();
               entities_high.clear();
               entities_hit.clear();
               return CUBIT_FAILURE;
            }
            
            ilow = 0;
            range = range_orig - offset;
         }
         else {
              // try this - don't even bother firing another low ray (it'll
              // fail anyway) - just put the first high hit on the list and
              // make sure we have an even number of hits
            ray_params.append(range - params_high[0]);
            if (entity_list) {
               topo_entity_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity( entities_low[ilow] );
               geom_entity = CAST_TO( topo_entity_ptr, GeometryEntity );
               if (geom_entity) entity_list->append(geom_entity);
            }
            if (ray_params.size() % 2 == 1) {
               ray_params.last();
               ray_params.append(ray_params.get());
               if (entity_list) {
                  entity_list->last();
                  entity_list->append(entity_list->get());
               }
            }
            if (entities_low.count() == 1 && entities_high.count() == 1)
            {
               entities_low.clear();
               entities_high.clear();
               entities_hit.clear();
               return CUBIT_SUCCESS;
            }
            else 
            {
               entities_low.clear();
               entities_high.clear();
               entities_hit.clear();
               return CUBIT_FAILURE;
            }
         }
      }
   }
   
   if (ray_params.size()%2) {
        // number of hits odd - put another hit in list in case it's a
        // tangency hit; for num_hit > 1, need a point_in_body check
      SPAposition pos;
      int index = 0;
      pos_low -= offset * vector_low; 
      int i;
      
      ray_params.reset();
      for (i = 1; i < ray_params.size(); i+=2) {
         point_containment pc;
         pos = pos_low + 0.5 * (ray_params.next(i-1) + ray_params.next(i)) * vector_low;
         rc = api_point_in_body (pos, body_ACIS, pc);
         if (!rc.ok()) {
            PRINT_ERROR("Point in body during ray fire failed for body \n");
            ACIS_API_error(rc, "Ray fire point in body");
            entities_low.clear();
            entities_high.clear();
            entities_hit.clear();
            return CUBIT_FAILURE;
         }
         
         if (pc == point_outside) index = i-1;
      }
      
        // index is the point that needs to be repeated; first bubble
        // the list, then replace
      ray_params.reset();
      ray_params.step(index);
      ray_params.insert(ray_params.get());
      if (entity_list) {
         entity_list->reset();
         entity_list->step(index);
         entity_list->insert(entity_list->get());
      }
   }
   entities_low.clear();
   entities_high.clear();
   entities_hit.clear();
   return CUBIT_SUCCESS;
}

int AcisQueryEngine::get_major_version() 
{
  return ::get_major_version();
}

int AcisQueryEngine::get_minor_version()
{
  return ::get_minor_version();
}

int AcisQueryEngine::get_subminor_version() 
{
  return AQE_SUBMINOR_VERSION;
}

int AcisQueryEngine::get_allint_version()
{
     // what's the api function for this????
    //Okay, the api follows:
//   int maj_api_version;
//   int min_api_version;
//   outcome foo = api_get_save_version(maj_api_version,min_api_version);
//But the api doesn't provide the point version, so for now we'll look
//at the other way of doing things
   return 100*get_major_version() + 10*get_minor_version() + 
    ::get_point_version();
}

CubitString AcisQueryEngine::get_engine_version_string()
{
   CubitString number_string = CubitString(get_major_version());
   number_string += ".";
   number_string += CubitString(get_minor_version());
   number_string += ".";
   number_string += CubitString(::get_point_version());
   number_string += ".";
   number_string += CubitString(get_subminor_version());
   CubitString version_string = "ACIS Version ";
   version_string += number_string;
   return version_string;
}

CubitStatus AcisQueryEngine::set_export_allint_version(const int version)
{
  if(version == CUBIT_INT_MAX)
  {
    int this_version = get_allint_version();
    if(exportVersion == CUBIT_INT_MAX)
       exportVersion = this_version;
    
    PRINT_INFO("Current User-Set Geometry Version = %d\n",exportVersion);
      //Commenting this out because the api_get_save_version always
      //returns the library version, not the user-set version
//     int maj_api_version;
//     int min_api_version;
//     outcome foo = api_get_save_version(maj_api_version,min_api_version);
//     PRINT_INFO("Current set version = %d%d\n",maj_api_version,min_api_version);
    return CUBIT_SUCCESS;
  }
  
  int this_major = version / 100;
  int this_minor = version % 100;
  
  return set_export_version(this_major, this_minor);
}

CubitStatus AcisQueryEngine::set_export_version(const int major,
                                                const int minor) 
{
  bool valid = false;

  int this_minor = minor;
  
  if ((major == 1 && 
     (this_minor == 6 || this_minor == 7)))
    valid = true;
  else if (major == 2 && 
           (this_minor == 1))
    valid = true;
  
  else if (major == 3 && 
           (this_minor == 0 || this_minor == 1)) {
    this_minor = 0;
    valid = true;
  }
    
  else if (major == 4 && 
           (this_minor == 0 || this_minor == 1 || this_minor == 2 || this_minor == 3)) {
    this_minor = 0;
    valid = true;
  }

  else if (major == 5 && 
           (this_minor == 0 || this_minor == 1 || this_minor == 2 || this_minor == 3)) {
    if (this_minor != 3) this_minor = 0;
    valid = true;
  }
  
  else if (major == 6 && 
           (this_minor == 0 || this_minor == 1 || this_minor == 2 || this_minor == 3)) {
    this_minor = 0;
    valid = true;
  }
  
  else if (major == 7 && 
           (this_minor >= 0 && this_minor <= 7)) {
    this_minor = 0;
    valid = true;
  }
  
  else if (major == 8 && this_minor == 0) {
    valid = true;
  }
  
  else if (major == 10 && this_minor == 7) {
    this_minor = 0;
    valid = true;
  }
  
  else if (major == 11 && this_minor >= 0 && this_minor <= 10) {
    this_minor = 0;
    valid = true;
  }
  
  else if (major == 12 && this_minor >= 0 && this_minor <= 7) {
    this_minor = 0;
    valid = true;
  }
  
  else if (major == 13 && this_minor >= 0 && this_minor <= 6) {
    this_minor = 0;
    valid = true;
  }

  else if (major == 14 && this_minor >= 0 && this_minor <= 7) {
    this_minor = 0;
    valid = true;
  }

  else if (major == 15 && this_minor >= 0 && this_minor <= 5) {
    this_minor = 0;
    valid = true;
  }
  else if (major == 16 && this_minor >= 0 && this_minor <= 3) {
    this_minor = 0;
    valid = true;
  }
  else if (major == 16 && this_minor >= 0 && this_minor <= 5) {
    this_minor = 0;
    valid = true;
  }
  else {
    PRINT_ERROR("Unrecognized Acis version number %d.%d.\n", major, this_minor);
    return CUBIT_FAILURE;
  }
  
  int major_version = get_major_version();
  int minor_version = get_minor_version();
   
  if (major > major_version ||
      (major == major_version && 
       this_minor > minor_version)) {
    PRINT_ERROR("Can't set engine version later than current version;"
                " current version = %d.%d\n", major_version, minor_version);
    return CUBIT_FAILURE;
  }
   
  outcome rc = api_save_version(major, this_minor);
  if (!rc.ok()) {
    ACIS_API_error(rc, "Set engine version");
    return CUBIT_FAILURE;
  }

  exportVersion = 100*major + minor;

  return CUBIT_SUCCESS;
}   

CubitStatus AcisQueryEngine::list_engine_versions(CubitString &versions)
{
   int major_version = get_major_version();
   int point_version = get_point_version();
   
   versions = "106";
   if (major_version > 1 || point_version >= 7) versions += ", 107"; // 107
   if (major_version > 2 || 
       (major_version == 2 && point_version >= 0)) versions += ", 200"; // 200
   if (major_version > 2 || 
       (major_version == 2 && point_version >= 1)) versions += ", 201"; // 201
   if (major_version > 3 || 
       (major_version == 3 && point_version >= 0)) versions += ", 300"; // 300
   if (major_version > 3 || 
       (major_version == 3 && point_version >= 1)) versions += ", 301"; // 301
   if (major_version > 4 || 
       (major_version == 4 && point_version >= 0)) versions += ", 400"; // 400
   if (major_version > 4 || 
       (major_version == 4 && point_version >= 1)) versions += ", 401"; // 401
   if (major_version > 4 || 
       (major_version == 4 && point_version >= 2)) versions += ", 402"; // 402
   if (major_version > 4 || 
       (major_version == 4 && point_version >= 3)) versions += ", 403"; // 403
   if (major_version > 5 || 
       (major_version == 5 && point_version >= 0)) versions += ", 500"; // 500
   if (major_version > 5 || 
       (major_version == 5 && point_version >= 1)) versions += ", 501"; // 501
   if (major_version > 6 || 
       (major_version == 6 && point_version >= 0)) versions += ", 600"; // 600
   if (major_version > 6 || 
       (major_version == 6 && point_version >= 1)) versions += ", 601"; // 601
   if (major_version > 6 || 
       (major_version == 6 && point_version >= 2)) versions += ", 602"; // 602
   if (major_version > 6 || 
       (major_version == 6 && point_version >= 3)) versions += ", 603"; // 603
   if (major_version > 7 || 
       (major_version == 7 && point_version >= 0)) versions += ", 700"; // 603
   if (major_version > 7 || 
       (major_version == 7 && point_version >= 7)) versions += ", 707"; // 707
   if (major_version > 8 || 
       (major_version == 8 && point_version >= 0)) versions += ", 800"; // 800
   return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Finds number of the VOLUMES associated with the input 
//                 BODY.
//
// Special Notes : Returns CUBIT_FAILURE if unable to determine number of VOLUMES
//
// Creator       : Steve Storm
//
// Creation Date : 11-13-98
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::number_VOLUMES( BODY* BODY_ptr, int& num_volumes ) const
{
   // Check for a NULL input BODY
   if (!BODY_ptr)
   {
      num_volumes = 0;
      return CUBIT_SUCCESS;
   }
   
   ENTITY_LIST VOLUMES;
   outcome result;
   
   result = api_get_lumps( BODY_ptr, VOLUMES);
   if( !result.ok() )
      return CUBIT_FAILURE;

   num_volumes = VOLUMES.count();

   return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Finds number of the FACES associated with the input 
//                 BODY.
//
// Special Notes : Returns CUBIT_FAILURE if unable to determine number of FACES
//
// Creator       : Steve Storm
//
// Creation Date : 11-13-98
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::number_FACES( BODY* BODY_ptr, int& num_faces ) const
{
   // Check for a NULL input BODY
   if (!BODY_ptr)
   {
      num_faces = 0;
      return CUBIT_SUCCESS;
   }
   
   ENTITY_LIST FACES;
   outcome result;
   
   result = api_get_faces( BODY_ptr, FACES);
   if( !result.ok() )
      return CUBIT_FAILURE;

   num_faces = FACES.count();

   return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Finds number of the EDGES associated with the input 
//                 BODY.
//
// Special Notes : Returns CUBIT_FAILURE if unable to determine number of EDGES
//
// Creator       : Steve Storm
//
// Creation Date : 11-13-98
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::number_EDGES( BODY* BODY_ptr, int& num_edges ) const
{
   // Check for a NULL input BODY
   if (!BODY_ptr)
   {
      num_edges = 0;
      return CUBIT_SUCCESS;
   }
   
   ENTITY_LIST EDGES;
   outcome result;
   
   result = api_get_edges( BODY_ptr, EDGES);
   if( !result.ok() )
      return CUBIT_FAILURE;

   num_edges = EDGES.count();

   return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Finds number of the VERTICES associated with the input 
//                 BODY.
//
// Special Notes : Returns CUBIT_FAILURE if unable to determine number of VERTICEs
//
// Creator       : Steve Storm
//
// Creation Date : 11-13-98
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::number_VERTICES( BODY* BODY_ptr, int& num_vertices ) const
{
   // Check for a NULL input BODY
   if (!BODY_ptr)
   {
      num_vertices = 0;
      return CUBIT_SUCCESS;
   }
   
   ENTITY_LIST VERTICES;
   outcome result;
   
   result = api_get_edges( BODY_ptr, VERTICES);
   if( !result.ok() )
      return CUBIT_FAILURE;

   num_vertices = VERTICES.count();

   return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Finds number of the FACES, EDGES and VERTICES associated 
//                 with the input BODY.
//
// Special Notes : Returns CUBIT_FAILURE if unable to determine number of ENTITIES
//
// Creator       : Steve Storm
//
// Creation Date : 11-13-98
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::number_ENTITIES( BODY* BODY_ptr, int &num_volumes,
                                                 int& num_faces, int& num_edges,
                                                 int& num_vertices ) const
{
   // Check for a NULL input BODY
   if (!BODY_ptr)
   {
      num_volumes = 0;
      num_faces = 0;
      num_edges = 0;
      num_vertices = 0;
      return CUBIT_SUCCESS;
   }
   
   if( number_VOLUMES( BODY_ptr, num_volumes ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

   if( number_FACES( BODY_ptr, num_faces ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

   if( number_EDGES( BODY_ptr, num_edges ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

   if( number_VERTICES( BODY_ptr, num_vertices ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

   return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::bodysms(ENTITY *entity, DLIList<BodySM*> &bodies) const
{
  ENTITY_LIST entities;
  CubitStatus status = get_ENTITY_from_ENTITY(BODY_TYPE, entity, entities);
  if (status == CUBIT_FAILURE)
    return status;
  
  DLIList<TopologyBridge*> tb_list;
  ATTRIB_CUBIT_OWNER::cubit_owner(entities, tb_list);
  CAST_LIST(tb_list, bodies, BodySM);
  return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::lumps(ENTITY *entity, DLIList<Lump*> &lumps) const
{
  ENTITY_LIST entities;
  CubitStatus status = get_ENTITY_from_ENTITY(LUMP_TYPE, entity, entities);
  if (status == CUBIT_FAILURE) return status;
  
  DLIList<TopologyBridge*> tb_list;
  ATTRIB_CUBIT_OWNER::cubit_owner(entities, tb_list);
  CAST_LIST(tb_list, lumps, Lump);
  return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::shellsms(ENTITY *entity, DLIList<ShellSM*> &shellsms) const
{
  ENTITY_LIST entities;
  CubitStatus status = get_ENTITY_from_ENTITY(SHELL_TYPE, entity, entities);
  if (status == CUBIT_FAILURE) return status;
  
  DLIList<TopologyBridge*> tb_list;
  ATTRIB_CUBIT_OWNER::cubit_owner(entities, tb_list);
  CAST_LIST(tb_list, shellsms, ShellSM);
  return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::surfaces(ENTITY *entity, DLIList<Surface*> &surfaces) const
{
  ENTITY_LIST entities;
  CubitStatus status = get_ENTITY_from_ENTITY(FACE_TYPE, entity, entities);
  if (status == CUBIT_FAILURE) return status;
  
  DLIList<TopologyBridge*> tb_list;
  ATTRIB_CUBIT_OWNER::cubit_owner(entities, tb_list);
  CAST_LIST(tb_list, surfaces, Surface);
  return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::loopsms(ENTITY *entity, DLIList<LoopSM*> &loopsms) const
{
  ENTITY_LIST entities;
  CubitStatus status = get_ENTITY_from_ENTITY(LOOP_TYPE, entity, entities);
  if (status == CUBIT_FAILURE) return status;
  
  DLIList<TopologyBridge*> tb_list;
  ATTRIB_CUBIT_OWNER::cubit_owner(entities, tb_list);
  CAST_LIST(tb_list, loopsms, LoopSM);
  return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::curves(ENTITY *entity, DLIList<Curve*> &curves) const
{
  ENTITY_LIST entities;
  CubitStatus status = get_ENTITY_from_ENTITY(EDGE_TYPE, entity, entities);
  if (status == CUBIT_FAILURE) return status;
  
  DLIList<TopologyBridge*> tb_list;
  ATTRIB_CUBIT_OWNER::cubit_owner(entities, tb_list);
  CAST_LIST(tb_list, curves, Curve);
  return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::coedgesms(ENTITY *entity, DLIList<CoEdgeSM*> &coedgesms) const
{
  ENTITY_LIST entities;
  CubitStatus status = get_ENTITY_from_ENTITY(COEDGE_TYPE, entity, entities);
  if (status == CUBIT_FAILURE) return status;
  
  DLIList<TopologyBridge*> tb_list;
  ATTRIB_CUBIT_OWNER::cubit_owner(entities, tb_list);
  CAST_LIST(tb_list, coedgesms, CoEdgeSM);
  return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::points(ENTITY *entity, DLIList<Point*> &points) const
{
  ENTITY_LIST entities;
  CubitStatus status = get_ENTITY_from_ENTITY(VERTEX_TYPE, entity, entities);
  if (status == CUBIT_FAILURE)
    return status;
  
  DLIList<TopologyBridge*> tb_list;
  ATTRIB_CUBIT_OWNER::cubit_owner(entities, tb_list);
  CAST_LIST(tb_list, points, Point);
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : To determine if two ACIS ENTITYs are related 
//               
// Special Notes : 
//
// Creator       : Wes Gill (Tim Tautges)
//
// Creation Date : 5-10-01
//-------------------------------------------------------------------------
CubitBoolean AcisQueryEngine::is_Related(ENTITY *entity1, ENTITY *entity2) const
{
  ENTITY_LIST entities;
  int i;
  //if they are the same entity, they are related
  if (entity1 == entity2)
  {   
    return CUBIT_TRUE;
  }
  
  //if they are not the same but the same topological type return CUBIT_FALSE
  else if (entity1->identity() == entity2->identity()) return CUBIT_FALSE;

  else 
  {
     get_ENTITY_from_ENTITY(entity1, entity2, entities);
     for (i = 0; i < entities.count(); i++)
     {
       //if entity2 is in the list, the entities are related
       if (entity2 == entities[i]) return CUBIT_TRUE;
     }
  }
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Query from_entity for related topology of the same 
//                 base type as target_type.
//
// Special Notes : Replaces existing implementation that was broken for
//                 subclasses of base topology types (E.g. TEDGE).
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/03/03
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::get_ENTITY_from_ENTITY( ENTITY* from_entity,
                                                     const ENTITY* target_type,
                                                     ENTITY_LIST& result_set) const
{
  outcome result;
  
  if (IS_ENTITY_TYPE( target_type, BODY ))
  {
    BODY* body = get_BODY_of_ENTITY(from_entity);
    if (body)
      result_set.add(body);
  }
  else if (IS_ENTITY_TYPE( target_type, LUMP ))
    result = api_get_lumps( from_entity, result_set );
  else if (IS_ENTITY_TYPE( target_type, SHELL ))
    result = api_get_shells( from_entity, result_set );
  else if (IS_ENTITY_TYPE( target_type, FACE ))
    result = api_get_faces( from_entity, result_set );
  else if (IS_ENTITY_TYPE( target_type, LOOP ))
    result = api_get_loops( from_entity, result_set );
  else if (IS_ENTITY_TYPE( target_type, COEDGE ))
    result = api_get_coedges( from_entity, result_set );
  else if (IS_ENTITY_TYPE( target_type, EDGE ))
    result = api_get_edges( from_entity, result_set );
  else if (IS_ENTITY_TYPE( target_type, VERTEX ))
    result = api_get_vertices( from_entity, result_set );
  else
  {
    PRINT_ERROR("Internal error at %s:%d\n\tUnexpected ENTITY type \"%s\"\n",
      __FILE__, __LINE__, target_type->type_name() ? target_type->type_name() : "(null)");
    return CUBIT_FAILURE;
  }
  
  return result.ok() ? CUBIT_SUCCESS : CUBIT_FAILURE;
}
                                                     

CubitStatus AcisQueryEngine::get_ENTITY_from_ENTITY(const int ident,
                                                       ENTITY *entity,
                                                       ENTITY_LIST &entities) const
{

    //- topology traversal of TB's; implemented based on native traversal
    //- functions in the modeler; returns lists of acis entities
  outcome result;
  BODY *body;
  
  if (ident == BODY_TYPE)
  {
    body = get_BODY_of_ENTITY(entity);
    if (body)
      entities.add(body);
  }
  else if (ident == LUMP_TYPE)
  {
    result = api_get_lumps( entity, entities);
  }
  else if (ident == LOOP_TYPE)
  {
    result = api_get_loops( entity, entities);
  }
  else if (ident == SHELL_TYPE)
  {
    result = api_get_shells( entity, entities);
  }
  else if (ident == FACE_TYPE)
  {
    result = api_get_faces( entity, entities);
  }
  else if (ident == EDGE_TYPE)
  {
    result = api_get_edges( entity, entities);
  }
  else if (ident == COEDGE_TYPE)
  {
    result = api_get_coedges( entity, entities);
  }
  else if (ident == VERTEX_TYPE)
  {
    result = api_get_vertices( entity, entities);
  }
  
  if (!result.ok())
    return CUBIT_FAILURE;
  
  return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::get_child_ENTITYs(ENTITY *from_entity,
                                               ENTITY_LIST &entities,
                                               CubitBoolean get_all) const
{
  CubitBoolean continue_traversal = CUBIT_FALSE;
  
  if (IS_ENTITY_TYPE( from_entity, BODY ))
  {
    get_ENTITY_from_ENTITY(LUMP_TYPE, from_entity, entities);
    if (get_all)
      continue_traversal = CUBIT_TRUE;
  }
  if (IS_ENTITY_TYPE( from_entity, LUMP ) || continue_traversal == CUBIT_TRUE)
  {
    get_ENTITY_from_ENTITY(SHELL_TYPE, from_entity, entities);
    if (get_all)
      continue_traversal = CUBIT_TRUE;
  }
  if (IS_ENTITY_TYPE( from_entity, SHELL ) || continue_traversal == CUBIT_TRUE)
  {
    get_ENTITY_from_ENTITY(FACE_TYPE, from_entity, entities);
    if (get_all)
      continue_traversal = CUBIT_TRUE;
  }
  if (IS_ENTITY_TYPE( from_entity, FACE ) || continue_traversal == CUBIT_TRUE)
  {
    get_ENTITY_from_ENTITY(LOOP_TYPE, from_entity, entities);
    if (get_all)
      continue_traversal = CUBIT_TRUE;
  }
  if (IS_ENTITY_TYPE( from_entity, LOOP ) || continue_traversal == CUBIT_TRUE)
  {
    get_ENTITY_from_ENTITY(COEDGE_TYPE, from_entity, entities);
    if (get_all)
      continue_traversal = CUBIT_TRUE;
  }
  if (IS_ENTITY_TYPE( from_entity, COEDGE ) || continue_traversal == CUBIT_TRUE)
  {
    get_ENTITY_from_ENTITY(EDGE_TYPE, from_entity, entities);
    if (get_all)
      continue_traversal = CUBIT_TRUE;
  }
  if (IS_ENTITY_TYPE( from_entity, EDGE ) || continue_traversal == CUBIT_TRUE)
  {
    get_ENTITY_from_ENTITY(VERTEX_TYPE, from_entity, entities);
  }
  
  return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::get_FACEs( ENTITY* ENTITY_ptr, 
                                        DLIList<FACE*>& FACE_list ) const
{
   // Check for a NULL input ENTITY
   if (!ENTITY_ptr)
      return CUBIT_FAILURE;

   ENTITY_LIST FACES;
   ENTITY *FACE_ptr;
   outcome result;
      
   result = api_get_faces( ENTITY_ptr, FACES);

   if( !result.ok() )
      return CUBIT_FAILURE;
   
   FACES.init();
   
   while( (FACE_ptr = FACES.next() ) != NULL )
   {
      FACE_list.append( (FACE *)FACE_ptr );
   }
   
   FACES.clear();
   
   return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::get_LOOPs( ENTITY* ENTITY_ptr, 
                                        DLIList<LOOP*>& LOOP_list ) const
{
   // Check for a NULL input ENTITY
   if (!ENTITY_ptr)
      return CUBIT_FAILURE;

   ENTITY_LIST LOOPS;
   ENTITY *LOOP_ptr;
   outcome result;
      
   result = api_get_loops( ENTITY_ptr, LOOPS);

   if( !result.ok() )
      return CUBIT_FAILURE;
   
   LOOPS.init();
   
   while( (LOOP_ptr = LOOPS.next() ) != NULL )
   {
      LOOP_list.append( (LOOP *)LOOP_ptr );
   }
   
   LOOPS.clear();
   
   return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::get_EDGEs( ENTITY* ENTITY_ptr, 
                                           DLIList<EDGE*>& EDGE_list ) const
{
   // Check for a NULL input ENTITY
   if (!ENTITY_ptr)
      return CUBIT_FAILURE;

   ENTITY_LIST EDGES;
   ENTITY *EDGE_ptr;
   outcome result;
      
   result = api_get_edges( ENTITY_ptr, EDGES);

   if( !result.ok() )
      return CUBIT_FAILURE;
   
   EDGES.init();
   
   while( (EDGE_ptr = EDGES.next() ) != NULL )
   {
      EDGE_list.append( (EDGE *)EDGE_ptr );
   }
   
   EDGES.clear();
   
   return CUBIT_SUCCESS;
}

CubitStatus AcisQueryEngine::get_VERTICEs( ENTITY* ENTITY_ptr, 
                                              DLIList<VERTEX*>& VERTEX_list ) const
{
   // Check for a NULL input ENTITY
   if (!ENTITY_ptr)
      return CUBIT_FAILURE;

   ENTITY_LIST VERTICES;
   ENTITY *VERTEX_ptr;
   outcome result;
      
   result = api_get_vertices( ENTITY_ptr, VERTICES);

   if( !result.ok() )
      return CUBIT_FAILURE;
   
   VERTICES.init();
   
   while( (VERTEX_ptr = VERTICES.next() ) != NULL )
   {
      VERTEX_list.append( (VERTEX *)VERTEX_ptr );
   }
   
   VERTICES.clear();
   
   return CUBIT_SUCCESS;
}

double AcisQueryEngine::volume(BODY *this_body) const
{
  
    // The ACIS mass_properties functions computes a number of properties
    // including volume and a host of moments of inertia.  Setting 
    // compute_type to 2 requests the computation of only the volume --
    // the other arguments are dummy arguments.
    
  double massprop_accuracy_requested = .01;
  SPAposition origin = SPAposition(0.0, 0.0, 0.0);
  SPAunit_vector normal = z_axis;
  double desired_accuracy = massprop_accuracy_requested;
  SPAposition ignored_cofg(0.0, 0.0, 0.0);
  tensor ignored_inertia; ignored_inertia.zero();
  double ignored_p_moments[3];
  SPAunit_vector ignored_p_axes[3];
  double accuracy_achieved = 0.0;
  int compute_type = 2;
  double volume;
    
  api_body_mass_pr(this_body, origin, normal, compute_type, desired_accuracy,
                   volume, ignored_cofg, ignored_inertia,
                   ignored_p_moments, ignored_p_axes, accuracy_achieved);
  return volume;
}

// Functions following added by CAT
double AcisQueryEngine::get_sme_resabs_tolerance() const
{
   return SPAresabs;
}

double AcisQueryEngine::set_sme_resabs_tolerance( double new_resabs )
{
   double old_resabs = SPAresabs;
   SPAresabs = new_resabs;
   return old_resabs;
}


CubitStatus
AcisQueryEngine::set_int_option( const char* opt_name, int val )
{
  outcome result;
  result = api_set_int_option( opt_name, val );
  if( result.ok() )
    return CUBIT_SUCCESS;
  PRINT_ERROR( "Unable to set ACIS integer option '%s' to %d\n", opt_name, val );
  return CUBIT_FAILURE;
}

CubitStatus
AcisQueryEngine::set_dbl_option( const char* opt_name, double val )
{
  outcome result;
  result = api_set_dbl_option( opt_name, val );
  if( result.ok() )
    return CUBIT_SUCCESS;
  PRINT_ERROR( "Unable to set ACIS double option '%s' to %f\n", opt_name, val );
  return CUBIT_FAILURE;
}

CubitStatus
AcisQueryEngine::set_str_option( const char* opt_name, const char* val )
{
  outcome result;
  result = api_set_str_option( opt_name, val );
  if( result.ok() )
    return CUBIT_SUCCESS;
  PRINT_ERROR( "Unable to set ACIS string option '%s' to '%s'\n", opt_name, val );
  return CUBIT_FAILURE;
}

CubitStatus
AcisQueryEngine::get_intersections( Curve* ref_edge1, CubitVector& point1,
                                    CubitVector& point2,
                                    DLIList<CubitVector*>& intersection_list,
                                    bool bounded,
                                    bool closest )
{
  if( closest == CUBIT_TRUE )
  {
    PRINT_ERROR( "Cannot currently find 'closest' intersections for ACIS curves\n" );
    return CUBIT_FAILURE;
  }
                                                                                
  // Get the ACIS EDGES
  EDGE* EDGE_ptr1 = get_EDGE(ref_edge1);
  if( EDGE_ptr1 == NULL )
  {
    PRINT_ERROR( "Unable to retrieve ACIS EDGE from Curve\n" );
    return CUBIT_FAILURE;
  }
                                                                                
  VERTEX* vertex1 = AcisModifyEngine::instance()->make_VERTEX(point1);
  VERTEX* vertex2 = AcisModifyEngine::instance()->make_VERTEX(point2);
  EDGE* EDGE_ptr2 = AcisModifyEngine::instance()->make_straight_EDGE(vertex1,
                                                                     vertex2);
  if( EDGE_ptr2 == NULL )
  {
    PRINT_ERROR( "Unable to create ACIS EDGE from points\n" );
    return CUBIT_FAILURE;
  }
                                                                                
  logical acis_bounded;
  if( bounded == CUBIT_TRUE )
    acis_bounded = TRUE;
  else
    acis_bounded = FALSE;
                                                                                
  if( get_intersections( EDGE_ptr1, EDGE_ptr2, intersection_list, bounded ) == CUBIT_FAILURE )
  {
     PRINT_ERROR( "unable to compute intersection of Curves\n" );
    return CUBIT_FAILURE;
  }
                                                                                
  return CUBIT_SUCCESS;
}

CubitStatus
AcisQueryEngine::get_intersections( Curve* ref_edge1, Curve* ref_edge2,
                                    DLIList<CubitVector*>& intersection_list,
                                    bool bounded,
                                    bool closest )
{
  if( closest == CUBIT_TRUE )
  {
    PRINT_ERROR( "Cannot currently find 'closest' intersections for ACIS curves\n" );
    return CUBIT_FAILURE;
  }

  // Get the ACIS EDGES
  EDGE* EDGE_ptr1 = get_EDGE(ref_edge1);
  if( EDGE_ptr1 == NULL )
  {
    PRINT_ERROR( "Unable to retrieve ACIS EDGE from Curve\n" );
    return CUBIT_FAILURE;
  }

  EDGE* EDGE_ptr2 = get_EDGE(ref_edge2);
  if( EDGE_ptr2 == NULL )
  {
    PRINT_ERROR( "Unable to retrieve ACIS EDGE from Curve\n" );
    return CUBIT_FAILURE;
  }

  logical acis_bounded;
  if( bounded == CUBIT_TRUE )
    acis_bounded = TRUE;
  else
    acis_bounded = FALSE;

  if( get_intersections( EDGE_ptr1, EDGE_ptr2, intersection_list, bounded ) == CUBIT_FAILURE )
  {
     PRINT_ERROR( "unable to compute intersection of Curves\n" );
    return CUBIT_FAILURE;
  }
  
  return CUBIT_SUCCESS;
}

CubitStatus
AcisQueryEngine::get_intersections( Curve* ref_edge, Surface* ref_face,
                                    DLIList<CubitVector*>& intersection_list,
                                    bool bounded )
{
   // Get the ACIS EDGE
  EDGE* EDGE_ptr = get_EDGE(ref_edge);
  if( EDGE_ptr == NULL )
  {
    PRINT_ERROR( "Unable to retrieve ACIS EDGE from Curve\n" );
    return CUBIT_FAILURE;
  }

  // Get the ACIS FACE
  FACE* FACE_ptr = get_FACE(ref_face);
  if( EDGE_ptr == NULL )
  {
    PRINT_ERROR( "Unable to retrieve ACIS FACE from Surface\n" );
    return CUBIT_FAILURE;
  }

  if( get_intersections( EDGE_ptr, FACE_ptr, intersection_list, bounded ) == CUBIT_FAILURE )
  {
     PRINT_ERROR( "unable to compute intersection of Curve and Surface\n" );
    return CUBIT_FAILURE;
  }
  
  return CUBIT_SUCCESS;
}

CubitStatus
AcisQueryEngine::get_intersections( EDGE *EDGE_ptr1, EDGE *EDGE_ptr2, 
                                    DLIList<CubitVector*> &intersection_list,
                                    CubitBoolean bounded ) const
{
   logical acis_bounded;
   if( bounded == CUBIT_TRUE )
      acis_bounded = TRUE;
   else
      acis_bounded = FALSE;
   
   curve_curve_int *cci = NULL;
   
   outcome result = api_intersect_curves(EDGE_ptr1, EDGE_ptr2, acis_bounded, cci);
   if( !result.ok() )
   {   
     // Let's not return an error - just return success but nothing will be in
     // the intersection list.  This function is used by other internal functions
     // so we don't want to spit out API errors, which usually just mean the
     // curves aren't interesecting.  Originally, we were trying to trap the
     // api error number that meant they weren't interesecting and return success
     // with that, failure otherwise.  But the error number seems to vary from 
     // ACIS release to ACIS release (4007, 4207, 43007....!).
     return CUBIT_SUCCESS;
   }
   
   curve_curve_int *next = NULL;
   CubitVector *vec;
   intersection_list.last();
   while( cci ) 
   {	
      next = cci->next;
      vec = new CubitVector( cci->int_point.x(), cci->int_point.y(), cci->int_point.z() );
      intersection_list.append( vec );
      delete cci;
      cci = next;	
   }
   return CUBIT_SUCCESS;
}

static SPAtransf const &get_edge_trans(EDGE *edge)
{
	// Find the body

	ENTITY *entity=NULL;
	if( edge->coedge()==NULL)
	{
		return *(SPAtransf *)NULL_REF;
	}
	else if(edge->coedge()->loop()!=NULL)
	{
		entity=edge->coedge()->loop()->face()->shell()->lump()->body();
	}
	else if(edge->coedge()->wire()!=NULL)
	{
		entity=edge->coedge()->wire()->body();
		if(!entity)
		{
			entity=edge->coedge()->wire()->shell()->lump()->body();
		}
	}
	else
	{
		return *(SPAtransf *)NULL_REF;
	}

	// Get the transform
	
	if( ((BODY *)entity)->transform() != NULL )
	{
		return ((BODY *)entity)->transform()->transform();
	}
	else
	{
		return *(SPAtransf *)NULL_REF;
	}
}

CubitStatus
AcisQueryEngine::get_intersections( EDGE *EDGE_ptr, FACE *FACE_ptr, 
                                    DLIList<CubitVector*> &intersection_list,
                                    CubitBoolean bounded ) const
{
  EDGE *int_EDGE_ptr = EDGE_ptr;
  CubitBoolean free_EDGE = CUBIT_FALSE;
  if( bounded == CUBIT_FALSE )
  {
     // Extend the curve
    double curve_length = EDGE_ptr->length();

    SPAbox super_box;
    super_box = get_acis_entity_bounding_box((ENTITY*)EDGE_ptr) |
      get_acis_entity_bounding_box((ENTITY*)FACE_ptr);

    // Try to extend the curve out some reasonable length
    double desired_length = curve_length + 
      get_max_size_of_acis_box( super_box ) * 10.0;

    double ratio = desired_length/curve_length;

    double start = EDGE_ptr->start_param();
    double end = EDGE_ptr->end_param();

    //PRINT_INFO( "Curve start = %f, Curve end = %f\n", start, end );
    
    double new_start = end + ratio*( start-end );
    double new_end = start + ratio*( end-start );

    //PRINT_INFO( "New start = %f, New end = %f\n", new_start, new_end );
   
    SPAtransf const &etrans=get_edge_trans( EDGE_ptr );
    curve *cur = EDGE_ptr->geometry()->trans_curve( etrans );
    SPAinterval newint( new_start, new_end );
    //SPAinterval test1=
      extend_curve( *cur, newint );
    
    // get ends of the edge
    SPAposition pt0=cur->eval_position( new_start );
    SPAposition pt1=cur->eval_position( new_end );
    
    // create the vertices
    VERTEX *svert= ACIS_NEW VERTEX( ACIS_NEW APOINT(pt0) );
    VERTEX *evert= ACIS_NEW VERTEX( ACIS_NEW APOINT(pt1) );	
    
    // create the curve
    CURVE *the_curve=make_curve( *cur );
#ifdef BOYD16
    ACIS_DELETE cur;
#endif
    
    // make the edge
    int_EDGE_ptr = ACIS_NEW EDGE(svert,evert,the_curve,FORWARD);	
    free_EDGE = CUBIT_TRUE;
    
    // set the param range of the edge
    SPAinterval range( new_start, new_end );
    int_EDGE_ptr->set_param_range( range );
  }

  // Note - perhaps we really should extend the FACE out too.  If we do that we'll need
  // to move this code into AcisModifyEngine as the surface extension code is there.

   ENTITY_LIST intersections;
   ENTITY_LIST* tmp_list = &intersections;
   outcome result  = api_edfa_int( int_EDGE_ptr, FACE_ptr, tmp_list );
   if( free_EDGE == CUBIT_TRUE )
     api_delent( int_EDGE_ptr );
   if (!result.ok())
   {
     // Let's not return an error - just the number of intersections in the 
     // list will be zero.
      return CUBIT_SUCCESS;
   }

   int num_ents = intersections.count(); 

   if( num_ents == 0 )
   {
     // Let's not return an error - just the number of intersections in the 
     // list will be zero.
     return CUBIT_SUCCESS;
   }

   intersections.init();
   ENTITY* entity_ptr = NULL;
   CubitVector *vec = NULL;
   intersection_list.last();
   for( int i=0; i<num_ents; i++ )
   {
      entity_ptr = intersections[i];
      if (IS_ENTITY_TYPE( entity_ptr, VERTEX )) 
      {
         VERTEX* VERTEX_ptr = (VERTEX *)entity_ptr;
         APOINT* point = VERTEX_ptr->geometry() ;
         SPAposition pos = point->coords() ;
         vec = new CubitVector( pos.x(), pos.y(), pos.z() ) ;
         intersection_list.append( vec );
      }
      api_delent( entity_ptr );
   }
   return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : entity_extrema
// Member Type: PUBLIC
// Description: Find extrema location on Cubit entities
// Author     : Steve Storm
// Date       : 11/02
//===============================================================================
CubitStatus 
AcisQueryEngine::entity_extrema( DLIList<GeometryEntity*> &entity_list, 
                                 const CubitVector *dir1, 
                                 const CubitVector *dir2,
                                 const CubitVector *dir3, 
                                 CubitVector &extrema,
                                 GeometryEntity *&extrema_entity_ptr )
{
  // Build up the ENTITY_list
  ENTITY_LIST ENTITY_list;
  ENTITY *ENTITY_ptr = NULL;

  for( int i=entity_list.size(); i--; )
  {
    GeometryEntity *entity_ptr = entity_list.get_and_step();

    if( (ENTITY_ptr = get_ENTITY_of_entity( entity_ptr )) == NULL )
    {    
      PRINT_ERROR( "problem occured getting ACIS entity for extrema calculation.\n"
        "       Aborting.\n" );
      return CUBIT_FAILURE;
    }

    ENTITY_list.add( ENTITY_ptr );
  }

  ENTITY *extrema_ENTITY_ptr = NULL;

  if( entity_extrema( ENTITY_list, dir1, dir2, dir3,
                      extrema, extrema_ENTITY_ptr ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Get extrema_ENTITY_ptr back to Cubit entity
  AcisBridge *tb_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( ENTITY_ptr );
  extrema_entity_ptr = tb_ptr ? CAST_TO(tb_ptr, GeometryEntity) : NULL;

  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : entity_extrema
// Member Type: PRIVATE
// Description: Find extrema location on ACIS entities
// Author     : Steve Storm
// Date       : 11/02
//===============================================================================
CubitStatus 
AcisQueryEngine::entity_extrema( ENTITY_LIST &ENTITY_list, 
                                 const CubitVector *dir1, 
                                 const CubitVector *dir2,
                                 const CubitVector *dir3, 
                                 CubitVector &extrema,
                                 ENTITY *&extrema_ENTITY_ptr )
{
  int num_dir = 1;
  SPAvector acis_dirs[3];
  acis_dirs[0].set_x( dir1->x() );
  acis_dirs[0].set_y( dir1->y() );
  acis_dirs[0].set_z( dir1->z() );

  if( dir2 )
  {
    num_dir++;
    acis_dirs[1].set_x( dir2->x() );
    acis_dirs[1].set_y( dir2->y() );
    acis_dirs[1].set_z( dir2->z() );

    if( dir3 )
    {
      num_dir++;
      acis_dirs[2].set_x( dir3->x() );
      acis_dirs[2].set_y( dir3->y() );
      acis_dirs[2].set_z( dir3->z() );
    }
  }

  SPAposition acis_pos;
	param_info info;
	outcome result = api_entity_extrema( ENTITY_list, num_dir, acis_dirs, 
                                       acis_pos, info );
  if( !result.ok() )
	{
    ACIS_API_error( result );
    PRINT_ERROR( "Unable to get entity extrema.\n" );
    return CUBIT_FAILURE;
	}

  extrema.x( acis_pos.x() );
  extrema.y( acis_pos.y() );
  extrema.z( acis_pos.z() );

  extrema_ENTITY_ptr = info.entity();

  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : entity_entity_distance
// Member Type: PUBLIC
// Description: Find distance between entities
// Author     : Steve Storm
// Date       : 11/02
//===============================================================================
CubitStatus
AcisQueryEngine::entity_entity_distance( GeometryEntity *ref_entity_ptr1,
                                         GeometryEntity *ref_entity_ptr2,
                                         CubitVector &pos1, CubitVector &pos2,
                                         double &distance )
{
  ENTITY *ENTITY_ptr1;
  if( (ENTITY_ptr1 = get_ENTITY_of_entity( ref_entity_ptr1 )) == NULL )
  {    
    PRINT_ERROR( "problem occured getting ACIS entity.\n"
      "       Aborting.\n" );
    return CUBIT_FAILURE;
  }

  ENTITY *ENTITY_ptr2;
  if( (ENTITY_ptr2 = get_ENTITY_of_entity( ref_entity_ptr2 )) == NULL )
  {    
    PRINT_ERROR( "problem occured getting ACIS entity.\n"
      "       Aborting.\n");
    return CUBIT_FAILURE;
  }
  
  return entity_entity_distance( ENTITY_ptr1, ENTITY_ptr2, pos1, pos2, distance );
}

//===============================================================================
// Function   : entity_entity_distance
// Member Type: PRIVATE
// Description: Find distance between ACIS entities
// Author     : Steve Storm
// Date       : 11/02
//===============================================================================
CubitStatus 
AcisQueryEngine::entity_entity_distance( ENTITY *ENTITY_ptr1,
                                         ENTITY *ENTITY_ptr2,
                                         CubitVector &pos1, CubitVector &pos2,
                                         double &distance )
{
  SPAposition acis_pos1, acis_pos2;

  outcome result = api_entity_entity_distance( ENTITY_ptr1, ENTITY_ptr2,
                                               acis_pos1, acis_pos2, distance );

  if( !result.ok() )
	{
    ACIS_API_error( result );
    PRINT_ERROR( "Unable to get distance between entities.\n");
    return CUBIT_FAILURE;
	}

  // Get positions back to Cubit format
  pos1.x( acis_pos1.x() );
  pos1.y( acis_pos1.y() );
  pos1.z( acis_pos1.z() );

  pos2.x( acis_pos2.x() );
  pos2.y( acis_pos2.y() );
  pos2.z( acis_pos2.z() );

  return CUBIT_SUCCESS;
}

CubitBoolean
AcisQueryEngine::is_curve_app_straight( Curve *ref_edge_ptr )
{
   EDGE* EDGE_ptr = get_EDGE(ref_edge_ptr);
   if( EDGE_ptr == NULL )
   {
      PRINT_ERROR( "Unable to retrieve ACIS EDGE from Curve\n" );
      return CUBIT_FALSE;
   }

   const curve* acis_curve = &EDGE_ptr->geometry()->equation();
   
   if( acis_curve->type() == straight_type )
      return CUBIT_TRUE;
   else if( acis_curve->type() == ellipse_type )
      return CUBIT_FALSE;
   else
   { 
      // Check the curve - it may be a spline that is about straight
      CubitVector start, end, one_qtr, center, three_qtr;
      ref_edge_ptr->position_from_fraction( .00, start );
      ref_edge_ptr->position_from_fraction( .25, one_qtr );
      center = ref_edge_ptr->center_point();
      ref_edge_ptr->position_from_fraction( .75, three_qtr );
      ref_edge_ptr->position_from_fraction( 1.0, end );
      
      AnalyticGeometryTool *agt = AnalyticGeometryTool::instance();
      
      double start_pnt[3],end_pnt[3],one_qtr_pnt[3],center_pnt[3],three_qtr_pnt[3];
      start.get_xyz( start_pnt ); end.get_xyz( end_pnt );
      one_qtr.get_xyz( one_qtr_pnt ); center.get_xyz( center_pnt );
      three_qtr.get_xyz( three_qtr_pnt );
      
      if( agt->is_pnt_on_ln_seg( center_pnt, start_pnt, end_pnt ) &&
         agt->is_pnt_on_ln_seg( one_qtr_pnt, start_pnt, end_pnt ) &&
         agt->is_pnt_on_ln_seg( three_qtr_pnt, start_pnt, end_pnt ) )
         return CUBIT_TRUE;
   }

   return CUBIT_FALSE;   
}

CubitStatus 
AcisQueryEngine::get_EDGE_normal( EDGE *EDGE_ptr, CubitVector &norm )
{
   /* //- Some edges can't be converted to bounded curves
   SPAunit_vector normal;
   const SPAtransf ftrans;
   bounded_curve* bnd_crv_ptr = new bounded_curve ( EDGE_ptr, &ftrans );
   normal = bnd_crv_ptr->get_normal();
   delete bnd_crv_ptr;
   
   norm.set( normal.x(), normal.y(), normal.z() );
   if( norm.x() == 0.0 && norm.y == 0.0 && norm.z() == 0.0 )
      return CUBIT_FAILURE;
  */

   EDGE *copied_EDGE_ptr = NULL;
   
   // Need to create a wire from the EDGE
   
   // Copy EDGE first
   
   
   api_edge( EDGE_ptr, copied_EDGE_ptr );  

   // Remove the owner attribute from the copied edge & children
   ATTRIB_CUBIT_OWNER::remove_cubit_owner(copied_EDGE_ptr, CUBIT_TRUE);
   
   
   // The copied EDGE is put into a wire so don't delete it.
   BODY* wire_BODY = NULL;
   outcome result = api_make_ewire( 1, &copied_EDGE_ptr, wire_BODY );
   
   if( !result.ok() || wire_BODY==NULL )
   {
      ACIS_API_error(result);
      PRINT_ERROR( "unable to make ACIS wire body from curve\n" );
      if( wire_BODY )
         api_delent( (ENTITY*)wire_BODY );
      else
      {
         api_delent( (ENTITY*)copied_EDGE_ptr );
      }
      return CUBIT_FAILURE;
   }
   
   WIRE* this_wire = wire_BODY->wire() ?
      wire_BODY->wire() : wire_BODY->lump()->shell()->wire(); 
   
   SPAposition centroid;
   SPAunit_vector normal;
   if(!is_planar_wire(this_wire, centroid, normal)) 
   {
      norm.set( 0.0, 0.0, 0.0 );
      api_delent( (ENTITY*)wire_BODY );
      return CUBIT_FAILURE;
   }
   
   api_delent( (ENTITY*)wire_BODY );
   
   norm.set( normal.x(), normal.y(), normal.z() );
   
   return CUBIT_SUCCESS;
}
#ifdef ACIS_STEP_TRANSLATOR

void
AcisQueryEngine::step_import_header_data_callback(step_header& /*header*/ )
{

}
#endif


//-------------------------------------------------------------------------
// Purpose       : Transform a body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::transform( BodySM* body, const SPAtransf& xform )
{
  BODY* BODY_ptr = get_BODY( body );
  if (!BODY_ptr)
    return CUBIT_FAILURE;
  
  return transform( BODY_ptr, xform );
}
  
CubitStatus AcisQueryEngine::transform( BODY *BODY_ptr, const SPAtransf& xform )
{
  outcome result = api_apply_transf( BODY_ptr, xform );
  if (!result.ok())
  {
    ACIS_API_error( result, "transforming BODY" );
    return CUBIT_FAILURE;
  }
  
  TRANSFORM* identity = new TRANSFORM( scale_transf(1.0) );
  result = api_change_body_trans( BODY_ptr, identity, FALSE );
  identity->lose();
  if (!result.ok())
  {
     ACIS_API_error( result, "transforming BODY" );
    return CUBIT_FAILURE;
  }

  clear_bounding_box( BODY_ptr );
  bounding_box( BODY_ptr );
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Remove applied transforms
//
// Special Notes : Moved from BodyACIS
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/25/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::restore_transform( BodySM* body )
{
  BODY* BODY_ptr = get_BODY( body );
  if (!BODY_ptr)
    return CUBIT_FAILURE;
  
   // Get the current transformation SPAmatrix of the BODY
  if ( BODY_ptr->transform() == NULL )
  {
    PRINT_INFO("Sorry, restore data was lost.\n");
    return CUBIT_FAILURE;
  }
  SPAtransf transformation = BODY_ptr->transform()->transform() ; 
   
  // If there has been any shear, the BODY cannot be restored.
  if (transformation.shear()) 
  {
     PRINT_ERROR("Can't restore body that has been sheared\n");
     return CUBIT_FAILURE;
  }
 
  SPAtransf inverse_transformation = transformation.inverse();
   
  if ( transform(body, inverse_transformation) )
  {
     return CUBIT_SUCCESS;
  }
   
  else
  {
    PRINT_ERROR("Problem restoring the BODY transform.\n");
    return CUBIT_FAILURE;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Transform a Surface, Curve, or Vertex
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::transform( GeometryEntity* ptr,
                                        const SPAtransf& xform )
{
  AcisBridge* bridge_ptr = dynamic_cast<AcisBridge*>(ptr);
  if (!bridge_ptr || !bridge_ptr->ENTITY_ptr())
  {
    PRINT_ERROR("Non-ACIS entity in AcisQueryEngine::transform.\n");
    return CUBIT_FAILURE;
  }
  
  ENTITY* entity = bridge_ptr->ENTITY_ptr();
  outcome result = api_apply_transf( entity, xform );
  if (!result.ok())
  {
    ACIS_API_error( result, "transforming entity" );
    return CUBIT_FAILURE;
  }
  
  if (IS_ENTITY_TYPE( entity, LUMP ))
  {
    clear_bounding_box( (LUMP*)entity );
    bounding_box( (LUMP*)entity );
  }
  else if (IS_ENTITY_TYPE( entity, FACE ))
  {
    clear_bounding_box( (FACE*)entity );
    bounding_box( (FACE*)entity );
  }
  else if (IS_ENTITY_TYPE( entity, EDGE ))
  {
    clear_bounding_box( (EDGE*)entity );
    bounding_box( (EDGE*)entity );
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Transform a body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::translate( BodySM* body, const CubitVector& d )
{
  return transform( body, translate_transf( SPAvector(d.x(),d.y(),d.z()) ) );
}

//-------------------------------------------------------------------------
// Purpose       : Transform a body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::rotate( BodySM* body, const CubitVector& v, double a )
{
  return transform( body,  rotate_transf( DEGREES_TO_RADIANS( a ),
                                         SPAvector(v.x(),v.y(),v.z()) ) );
}

//-------------------------------------------------------------------------
// Purpose       : Transform a body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::scale( BodySM* body, double f )
{ 
  return scale( body, CubitVector( f, f, f ) );
}

//-------------------------------------------------------------------------
// Purpose       : Transform a body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::scale( BodySM* body, const CubitVector& f )
{
#if defined(WIN32) || defined(MACOSX)
  api_initialize_warp();
#endif
//    BODY *original_body = NULL;

  CubitStatus result = transform( body, scale_transf( f.x(), f.y(), f.z() ) );

#if defined(WIN32) || defined(MACOSX)
  api_terminate_warp();
#endif

  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Transform a body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::reflect( BodySM* body, const CubitVector& v )
{
  return transform( body, reflect_transf( SPAvector( v.x(), v.y(), v.z() ) ) );
}



//-------------------------------------------------------------------------
// Purpose       : Transform a Surface, Curve, or Vertex
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::translate( GeometryEntity* geom, const CubitVector& d )
{
  return transform( geom, translate_transf( SPAvector(d.x(),d.y(),d.z()) ) );
}

//-------------------------------------------------------------------------
// Purpose       : Transform a Surface, Curve, or Vertex
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::rotate( GeometryEntity* geom, const CubitVector& v, double a )
{
  return transform( geom,  rotate_transf( DEGREES_TO_RADIANS( a ),
                                          SPAvector(v.x(),v.y(),v.z()) ) );
}

//-------------------------------------------------------------------------
// Purpose       : Transform a Surface, Curve, or Vertex
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::scale( GeometryEntity* geom, double f )
{ 
  return scale ( geom, CubitVector( f, f, f ) );
}

//-------------------------------------------------------------------------
// Purpose       : Transform a Surface, Curve, or Vertex
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::scale( GeometryEntity* geom, const CubitVector& f )
{
#if defined(WIN32)
  api_initialize_warp();
#endif
  
  CubitStatus result = transform( geom, scale_transf( f.x(), f.y(), f.z() ) );

#if defined(WIN32)
  api_terminate_warp();
#endif

  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Transform a Surface, Curve, or Vertex
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/24/04
//-------------------------------------------------------------------------
CubitStatus AcisQueryEngine::reflect( GeometryEntity* geom, const CubitVector& v )
{
  return transform( geom, reflect_transf( SPAvector( v.x(), v.y(), v.z() ) ) );
}

CubitBoolean AcisQueryEngine::volumes_overlap( Lump *lump1,
                                               Lump *lump2 ) const
{
  LUMP *acis_lump1 = AcisQueryEngine::get_LUMP( lump1 );
  LUMP *acis_lump2 = AcisQueryEngine::get_LUMP( lump2 );

  //Determine that the two bodies are overlapping by calling the
  //intersection function.
  CubitBox box_1 = bounding_box( bounding_box(acis_lump1) );
  CubitBox box_2 = bounding_box( bounding_box(acis_lump2) );
  if ( !box_1.overlap( GEOMETRY_RESABS, box_2 ) )
    return CUBIT_FALSE;

  //copy the lumps
  ENTITY *new_lump;
  api_copy_entity_contents( (ENTITY*)acis_lump1, new_lump );
  LUMP *copied_lump1 = (LUMP*)new_lump;
  api_copy_entity_contents( (ENTITY*)acis_lump2, new_lump );
  LUMP *copied_lump2 = (LUMP*)new_lump;
 
  //make new bodies
  BODY *acis_body_1 = new BODY( copied_lump1 );
  BODY *acis_body_2 = new BODY( copied_lump2 );

  //remove the attributes
  remove_cubit_owner_attrib_in_BODY(acis_body_1);
  remove_cubit_owner_attrib_in_BODY(acis_body_2);

  CubitBoolean bodies_overlap = CUBIT_FALSE;

#if CUBIT_ACIS_VERSION < 1600 
  bodies_overlap = CUBIT_TRUE;
  BOOL_TYPE bool_type = INTERSECTION;
  outcome result = api_boolean(acis_body_1, acis_body_2, bool_type );
  if (acis_body_2 == NULL ||
      acis_body_2->lump() == NULL && acis_body_2->wire() == NULL ||
      !result.ok())
  {
    bodies_overlap = CUBIT_FALSE;
  }
#else
  bodies_overlap = CUBIT_FALSE;
  body_clash_result clash_result;
  outcome result = api_clash_bodies( acis_body_1, acis_body_2, 
                                     clash_result, CLASH_CLASSIFY_BODIES );
  if (!result.ok())
  {
    PRINT_ERROR("Intersection operation failed.\n");
    AcisQueryEngine::instance()->ACIS_API_error(result, "intersect Bodies");
    return CUBIT_FALSE;
  }

  if( clash_result.clash_type() == CLASH_CONTAINED ||
      clash_result.clash_type() == CLASH_CONTAINED_ABUTS ||
      clash_result.clash_type() == CLASH_COINCIDENT ||
      clash_result.clash_type() == CLASH_INTERLOCK )
  {
    bodies_overlap = CUBIT_TRUE;
  }
#endif
  
  if ( acis_body_1 != NULL )
    AcisQueryEngine::instance()->delete_ACIS_BODY(acis_body_1);
  if ( acis_body_2 != NULL )
    AcisQueryEngine::instance()->delete_ACIS_BODY(acis_body_2);

  return bodies_overlap;
}

//-------------------------------------------------------------------------
// Purpose       : Determine if the two bodies overlap each other. touching
//                 is not considered overlaping.
// Creator       : David White
// Creation Data : 8/2/2003
//-------------------------------------------------------------------------
CubitBoolean AcisQueryEngine::bodies_overlap( BodySM *body_ptr_1,
                                              BodySM *body_ptr_2 ) const
{
  //Determine that the two bodies are overlapping by calling the
  //intersection function.
  CubitBox box_1 = bounding_box(body_ptr_1);
  CubitBox box_2 = bounding_box(body_ptr_2);
  if ( !box_1.overlap( GEOMETRY_RESABS, box_2 ) )
    return CUBIT_FALSE;

    //Get the acis bodies so we can do the intersection.
  BODY *acis_body_1_orig = AcisQueryEngine::get_BODY(body_ptr_1);
  BODY *acis_body_2_orig = AcisQueryEngine::get_BODY(body_ptr_2);

  CubitBoolean bodies_overlap = CUBIT_FALSE;

#if CUBIT_ACIS_VERSION < 1600 
  BODY *acis_body_1, *acis_body_2; 
  api_copy_body(acis_body_1_orig, acis_body_1);
  api_copy_body(acis_body_2_orig, acis_body_2);

  //remove the attributes
  remove_cubit_owner_attrib_in_BODY(acis_body_1);
  remove_cubit_owner_attrib_in_BODY(acis_body_2);

  bodies_overlap = CUBIT_TRUE;
  BOOL_TYPE bool_type = INTERSECTION;
  outcome result = api_boolean(acis_body_1, acis_body_2, bool_type );
  if (acis_body_2 == NULL ||
      acis_body_2->lump() == NULL && acis_body_2->wire() == NULL ||
      !result.ok())
  {
    bodies_overlap = CUBIT_FALSE;
  }
  if ( acis_body_2 != NULL )
    AcisQueryEngine::instance()->delete_ACIS_BODY(acis_body_2);
#else
  bodies_overlap = CUBIT_FALSE;
  body_clash_result clash_result;
  outcome result = api_clash_bodies( acis_body_1_orig, acis_body_2_orig, 
                                     clash_result, CLASH_CLASSIFY_BODIES );
  if (!result.ok())
  {
    PRINT_ERROR("Intersection operation failed.\n");
    AcisQueryEngine::instance()->ACIS_API_error(result, "intersect Bodies");
    return CUBIT_FALSE;
  }

  if( clash_result.clash_type() == CLASH_CONTAINED ||
      clash_result.clash_type() == CLASH_CONTAINED_ABUTS ||
      clash_result.clash_type() == CLASH_COINCIDENT ||
      clash_result.clash_type() == CLASH_INTERLOCK )
  {
    bodies_overlap = CUBIT_TRUE;
  }
#endif
  return bodies_overlap;
}

CubitStatus AcisQueryEngine::read_acis_file(const char* file_name, const char* file_type, ENTITY_LIST &entity_list)
{
  FILE *file_ptr = 0;
  // open the file - open as binary for .sab files on Windows
#ifdef NT
  if( strcmp(file_type, "ACIS_SAB") == 0 ) 
    file_ptr = fopen(file_name, "rb");
  else
#endif
    file_ptr = fopen(file_name, "r");
  if (!file_ptr)
  {
    PRINT_ERROR("Cannot open %s file: %s (%s)\n", file_type, file_name, strerror(errno));
    return CUBIT_FAILURE;
  }

  CubitStatus status = read_acis_file(file_ptr, file_type, entity_list);
  fclose(file_ptr);
  return status;
}

CubitStatus AcisQueryEngine::read_acis_file(FILE* file_ptr, const char* file_type, ENTITY_LIST &entity_list)
{
  logical ascii = ('T' == file_type[7]);

  ENTITY_LIST tmp_list;
  outcome result;
  do {
    listcat( entity_list, tmp_list );
    tmp_list.clear();
    result = api_restore_entity_list( file_ptr, ascii, tmp_list );
  } while (result.ok());
  remove_refinements( entity_list );

  return CUBIT_SUCCESS;
}

#ifdef ACIS_IGES_TRANSLATOR
CubitStatus AcisQueryEngine::read_iges_file(const char* file_name, const char* logfile_name, 
                                            ENTITY_LIST &entity_list, bool heal )
{
  char* logfilename = NULL;
  if( !logfile_name || !strcmp( logfile_name, "" ) )
    strcpy(logfilename, "iges_import.log");
  else
    logfilename = (char *)logfile_name;
  
  PRINT_INFO( "Reading IGES file '%s'...\n", file_name );

  FILE *logfile_ptr = fopen( logfilename, "w" );
  if( logfile_ptr == NULL )
  {
    PRINT_ERROR( "unable to open iges logfile '%s'\n", logfilename );
    return CUBIT_FAILURE;
  }

  SPAXProgressReport_fct = progress_report_function;
  AppUtil::instance()->progress_tool()->start(0, 100, "Importing IGES File:" );
  outcome result = api_xiges_read(entity_list, (char *)file_name, (char *)logfilename, heal);

  //convert all bodies to 2d sheet bodies
  ENTITY *tmp_ent;
  entity_list.init();
  while( (tmp_ent = entity_list.next()) != NULL) 
  {
    if (is_BODY( tmp_ent ) )
      api_body_to_2d( (BODY*)tmp_ent );
  }

  fclose( logfile_ptr );
  
  PRINT_INFO( "Done\n" );

  if( !result.ok() )
  {
    ACIS_API_error( result );
    PRINT_ERROR( "problem reading from IGES file '%s'\n", file_name );
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}
#endif

#ifdef ACIS_CATIA_TRANSLATOR
CubitStatus AcisQueryEngine::read_catia_file(const char* file_name, const char* logfile_name, ENTITY_LIST &entity_list)
{
  char* logfilename = NULL;
  if( !logfile_name || !strcmp( logfile_name, "" ) )
    strcpy(logfilename, "catia_import.log");
  else
    logfilename = (char *)logfile_name;
  
  int success; // Currently api_iges_read does not use this variable
  PRINT_INFO( "Reading CATIA file '%s'...\n", file_name );

  FILE *logfile_ptr = fopen( logfilename, "w" );
  if( logfile_ptr == NULL )
  {
    PRINT_ERROR( "unable to open CATIA logfile '%s'\n", logfilename );
    return CUBIT_FAILURE;
  }
  char temp_logfile_name[128];
  strcpy(temp_logfile_name, logfile_name);
  
  outcome result = api_catia_convert_catia_to_acisentlist(
      (char*)file_name, entity_list, temp_logfile_name);
  fclose( logfile_ptr );
  
  PRINT_INFO( "Done\n" );

  if( !result.ok() )
  {
    ACIS_API_error( result );
    PRINT_ERROR( "problem reading from CATIA file '%s'\n", file_name );
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}
#endif

#ifdef ACIS_STEP_TRANSLATOR
CubitStatus AcisQueryEngine::read_step_file(const char* file_name, const char* logfile_name, 
                                            ENTITY_LIST &entity_list, bool heal )
{
  if ( !stepInitialized ) {
    PRINT_ERROR("The STEP translater has not been properly initialized.\n");
    return CUBIT_FAILURE;
  }

  char* logfilename = NULL;
  if( !logfile_name || !strcmp( logfile_name, "" ) )
      strcpy(logfilename, "step_import.log");
  else
      logfilename = (char *)logfile_name;

  PRINT_INFO( "Reading STEP file '%s'...\n", file_name );
  
  SPAXProgressReport_fct = progress_report_function;

  api_xstep_set_option("step_to_acis_attrib_transf", TRUE);
  AppUtil::instance()->progress_tool()->start(0, 100, "Importing STEP File:" );
  outcome result = api_xstep_read( entity_list, (char *)file_name, (char *)logfilename, heal ); 

  if( !result.ok() )
  {
    if( result.error_number() == SPAX_ACIS_E_BOX_TOO_BIG )
    {
      const char* warning = find_err_ident(result.error_number());
      PRINT_WARNING("%s\n", find_err_mess( result.error_number() ) );
      PRINT_INFO("Your model is very large. Operations on your model could be \n" 
                 "numerically unstable.  Consider scaling your model down.\n", warning);
    }
    else if( result.error_number() == SPAX_ACIS_E_BOX_TOO_SMALL )
    {
      const char* warning = find_err_ident(result.error_number());
      PRINT_WARNING("%s\n", find_err_mess( result.error_number() ) );
      PRINT_INFO("Your model is very small. Operations on your model could be \n" 
                 "numerically unstable.  Consider scaling your model up.\n", warning);
    }
    else
    {
      ACIS_API_error( result );
      PRINT_ERROR( "problem reading from STEP file '%s'\n", file_name );
      AppUtil::instance()->progress_tool()->end();
      return CUBIT_FAILURE;
    }
  }
  AppUtil::instance()->progress_tool()->end();
  return CUBIT_SUCCESS;
}

void AcisQueryEngine::read_step_part_names(ENTITY_LIST &entity_list)
{
  ENTITY *currItem = NULL;
  entity_list.init();
  
  //extract the name attributes off the step parts
  while ((currItem = entity_list.next()) != NULL)
  {
    ATTRIB_GEN_STRING * label_attrib = NULL;
    outcome result = api_find_named_attribute( currItem, "ATTRIB_XSTEP_PRODUCT_ID",(ATTRIB_GEN_NAME*&)(label_attrib) );
    
    if (label_attrib != NULL)
    {
      const char *part_name = label_attrib->value();
      CubitSimpleAttrib *tmp_attrib = new CubitSimpleAttrib( "ENTITY_NAME", part_name );
      new ATTRIB_SNL_SIMPLE( currItem, tmp_attrib );
      delete tmp_attrib;
    }          
  }
}

// Parts cut in Pro/E Mechanica and exported as STEP come in as a single body.
// so this will auto separate the body 
void AcisQueryEngine::auto_separate_step_body(ENTITY_LIST &pre_separate_list, ENTITY_LIST &post_separate_list)
{
  pre_separate_list.init();

  BODY **new_BODY_list;
  int i;
  for ( i=0; i<pre_separate_list.count(); i++ )
  {
    if (is_solid_body (pre_separate_list[i]) )
    {
      int n =0,n_body =0;
      outcome result = api_separate_body( (BODY*)(pre_separate_list[i]), n_body, new_BODY_list );
      if ((n_body > 1) && result.ok())  //can't do much error checking on result since it returns
      {                                 //an error if there was only 1 body 
        PRINT_INFO( "Auto Separated Imported STEP body\n" );
        for (n=0; n<n_body; n++)
          post_separate_list.add(new_BODY_list[n] );
        continue;// next i, so we don't add it to entity_list twice
      }
    }
    post_separate_list.add( pre_separate_list[i] );
  }
}
#endif

#ifdef ACIS_PROE_TRANSLATOR
CubitStatus AcisQueryEngine::read_proe_file(const char* file_name, ENTITY_LIST &entity_list)
{
    // open the file
  FILE *file_ptr = fopen(file_name, "r");
  if (!file_ptr)
  {
    PRINT_ERROR("AcisQueryEngine::read_proe_file - Cannot open file: %s (%s)\n", file_name, strerror(errno));
    return CUBIT_FAILURE;
  }

    //If assembly.
  PRINT_INFO( "Reading PROE file '%s'...\n", file_name );
      outcome result = api_gss_xproe_read_part (file_ptr, entity_list); 
  PRINT_INFO( "Done\n" );
  
  if( !result.ok() )
  {
      ACIS_API_error( result );
      PRINT_ERROR( "problem reading from PROE Part file '%s'\n", file_name );
      PRINT_INFO("Trying to read it as a Assembly\n");
      result = api_gss_xproe_read_assembly(
        file_ptr,(char *)file_name, entity_list);
//         result = api_gss_xproe_read_part (file_ptr, entity_list); 

  }
  if( !result.ok() )
  {
      ACIS_API_error( result );
      PRINT_ERROR( "problem reading from PROE Part file '%s'\n", file_name );
      return CUBIT_FAILURE;
  }
  fclose(file_ptr);
  return CUBIT_SUCCESS;
}
#endif

