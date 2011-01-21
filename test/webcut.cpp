/**
 * \file mergechk.cpp
 *
 * \brief mergechk, another simple C++ driver for CGM
 *
 * This program acts as a simple driver for CGM.  It reads in a geometry,
 * performs imprints between all the bodies, merges them, and writes information
 * on the results.  It also performs pairwise intersections between the
 * bodies to check for overlaps.  Results are written to stardard output.
 *
 */

#include "CpuTimer.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "InitCGMA.hpp"
#include "CubitCompat.hpp"

#define STRINGIFY(S) XSTRINGIFY(S)
#define XSTRINGIFY(S) #S

// forward declare some functions used and defined later
void read_geometry(int, char **);
void read_model6( Body*& brick, Body*& sheet );
CubitStatus evaluate_overlaps();
CubitStatus imprint_bodies();
CubitStatus print_unmerged_surfaces();
CubitStatus webcut_with_brick();
CubitStatus webcut_with_cylinder();
CubitStatus webcut_with_planar_sheet();
CubitStatus webcut_with_sweep_curves_rotated();
CubitStatus webcut_with_sweep_curves();
CubitStatus webcut_with_sweep_along_curve();
CubitStatus webcut_with_curve_loop();
CubitStatus webcut_with_extended_surf();
CubitStatus webcut_with_sweep_surfaces_rotated();
CubitStatus webcut_with_sweep_surfaces();
CubitStatus webcut_with_sweep_surfaces_along_curve();
CubitStatus webcut_with_sweep_surfaces_perp();
// macro for printing a separator line
#define PRINT_SEPARATOR   PRINT_INFO("===================%s:%d===================\n",__FILE__,__LINE__);


// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  CubitStatus result = InitCGMA::initialize_cgma();
  if (CUBIT_SUCCESS != result) return 1;
  
  CubitStatus status = CUBIT_SUCCESS;

  //Do brick webcut.
  status = webcut_with_brick();
  status = webcut_with_cylinder();
  status = webcut_with_planar_sheet();
  status = webcut_with_sweep_curves_rotated();
  status = webcut_with_sweep_curves();
  status = webcut_with_sweep_along_curve();
  status = webcut_with_curve_loop();
  status = webcut_with_extended_surf();
  status = webcut_with_sweep_surfaces_rotated();
  status = webcut_with_sweep_surfaces();
  status = webcut_with_sweep_surfaces_along_curve();
  status = webcut_with_sweep_surfaces_perp();
  if (status == CUBIT_FAILURE) 
     PRINT_INFO("Operation Failed");

  int ret_val = ( CubitMessage::instance()->error_count() );
  if ( ret_val > 0 )
  {
    PRINT_ERROR("Errors found during Mergechk session.\n");
  }
  return ret_val;
  
}

/// attribs module: list, modify attributes in a give model or models
/// 
/// Arguments: file name(s) of geometry files in which to look
///
void read_geometry(int num_files, const char **argv) 
{
  CubitStatus status = CUBIT_SUCCESS;
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  assert(gti);
  int i;
  
    // For each file, open and read the geometry
  for (i = 0; i < num_files; i++) {
    status = CubitCompat_import_solid_model(argv[i], "ACIS_SAT");
    if (status != CUBIT_SUCCESS) {
      PRINT_ERROR("Problems reading geometry file %s.\n", argv[i]);
      abort();
    }
  }

  if (gti->num_bodies() == 0) {
    PRINT_ERROR("No bodies read; aborting.\n");
    abort();
  }
}

void read_model6( Body*& volume, Body*& sheet )
{
  const char *filename = STRINGIFY(SRCDIR) "/model6.sat";
  read_geometry( 1, &filename );
  
  DLIList<Body*> bodies;
  GeometryQueryTool::instance()->bodies(bodies);
  volume = sheet = 0;
  for (int i = 0; i < bodies.size(); ++i, bodies.step()) {
    DLIList<RefVolume*> vols;
    bodies.get()->ref_volumes( vols );
    if (vols.size() != 1) {
      PRINT_ERROR("Unexpected number of volumes in body read from file: %s\n", filename);
      abort();
    }
    if (vols.get()->is_sheet()) {
      if (sheet) {
        PRINT_ERROR("Unexpected number of sheet bodies read from file: %s\n",filename);
        abort();
      }
      sheet = bodies.get();
    }
    else {
      if (volume) {
        PRINT_ERROR("Unexpected number of non-sheet bodies read from file: %s\n",filename);
        abort();
      }
      volume = bodies.get();
    }
  }
}

CubitStatus webcut_with_brick()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read in the geometry from files specified on the command line
  const char *argv = STRINGIFY(SRCDIR) "/model2.sat";
  PRINT_SEPARATOR;
  read_geometry(1, &argv);
  
  //int num_bodies = gti->num_bodies();
  DLIList<Body*> old_bodies, new_bodies,junk;
  gti->bodies(old_bodies);
  old_bodies.reset();
  //old_bodies.remove();
  CubitVector center(4.5,4.5,4.5);
  CubitVector axes[3];
  axes[0].set(1.,0.,0.);
  axes[1].set(0.,1.,0.);
  axes[2].set(0.,0.,1.);
  CubitVector extension(0.5,0.5,0.5);
  CubitStatus rsl= gmti->webcut_with_brick(old_bodies,center,axes,extension,new_bodies,junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_brick.sat";
  const char * filetype = "ACIS_SAT";
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype, 
                                 num_ents_exported, cubit_version);

  gti->delete_geometry();
  return rsl;
}

CubitStatus webcut_with_cylinder()
{

  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read in the geometry from files specified on the command line
  const char *argv = STRINGIFY(SRCDIR) "/model2.sat";
  PRINT_SEPARATOR;
  read_geometry(1, &argv);
  
  DLIList<Body*> old_bodies, new_bodies, junk;
  gti->bodies(old_bodies);
  old_bodies.reset();
  //old_bodies.remove();
  CubitVector center(4.,4.,4.);
  CubitVector axis(1.,0.,0.);
  double radius = 1.0;
  CubitStatus rsl= gmti->webcut_with_cylinder(old_bodies,radius,axis, center,new_bodies,junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_cylinder.sat";
  const char * filetype = "ACIS_SAT";
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->delete_geometry();
  return rsl;
}

CubitStatus webcut_with_planar_sheet()
{

  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read in the geometry from files specified on the command line
  const char *argv = STRINGIFY(SRCDIR) "/model2.sat";
  PRINT_SEPARATOR;
  read_geometry(1, &argv);

  DLIList<Body*> old_bodies, new_bodies, junk;
  gti->bodies(old_bodies);
  old_bodies.reset();

  CubitVector center(4.0,4.0,4.0);
  CubitVector axes[2];
  axes[0].set(1.,1.,0.);
  axes[1].set(0.,1.,1.);  
  
  CubitStatus rsl= gmti->webcut_with_planar_sheet(old_bodies,center,axes,40.,40.,new_bodies,junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_sheet.sat";
  const char * filetype = "ACIS_SAT";
  rsl =  CubitCompat_export_solid_model(ref_entity_list, filename, filetype, 
                                 num_ents_exported, cubit_version);
  
  gti->delete_geometry();
  return rsl;
 
}


CubitStatus webcut_with_sweep_curves_rotated()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read in the geometry from files specified on the command line
  const char *argv = STRINGIFY(SRCDIR) "/huge.sat";
  PRINT_SEPARATOR;
  read_geometry(1, &argv);
  
  DLIList<Body*> old_bodies, new_bodies, junk;
  gti->bodies(old_bodies);
  DLIList<RefEntity*> free_entities;
  gti->get_free_ref_entities(free_entities);
 
  DLIList<RefEdge*> curves;
  for(int i = 0; i < free_entities.size(); i++)
  {
    RefEdge* free_edge = CAST_TO(free_entities.get_and_step(), RefEdge);
    curves.append(free_edge);
  }

  old_bodies.reset();

  //CubitVector center(4.0,14.0,14.0);
  CubitVector center(-50.0,50.0,50.0);
  CubitVector axis;
  axis.set(0.,-1.,0.);

  CpuTimer webcut_BODYs_timer;
  CubitStatus rsl= gmti->webcut_with_sweep_curves_rotated(old_bodies,curves,center,axis,1.58,NULL,new_bodies,junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  PRINT_INFO( "CPU time taken to WEBCUT this body: " 
               "%f secs\n", webcut_BODYs_timer.cpu_secs() ) ;

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_sweep_rotational.sat";
  const char * filetype = "ACIS_SAT";
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->delete_geometry();
  return rsl;
}

CubitStatus webcut_with_sweep_curves()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read in the geometry from files specified on the command line
  const char *argv = STRINGIFY(SRCDIR) "/model3.sat";
  PRINT_SEPARATOR;
  read_geometry(1, &argv);

  DLIList<Body*> old_bodies, new_bodies, junk;
  gti->bodies(old_bodies);
  DLIList<RefEntity*> free_entities;
  gti->get_free_ref_entities(free_entities);

  DLIList<RefEdge*> curves;
  for(int i = 0; i < free_entities.size(); i++)
  {
    RefEdge* free_edge = CAST_TO(free_entities.get_and_step(), RefEdge);
    curves.append(free_edge);
  }

  old_bodies.reset();

  CubitVector axis;
  axis.set(1.,0.,0.);

  CubitStatus rsl= gmti->webcut_with_sweep_curves(old_bodies,curves,axis,true,NULL,NULL,new_bodies,junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_sweep_curves.sat";
  const char * filetype = "ACIS_SAT";
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->delete_geometry();
  return rsl;
}

CubitStatus webcut_with_sweep_along_curve()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read in the geometry from files specified on the command line
  const char *argv = STRINGIFY(SRCDIR) "/model3.sat";
  PRINT_SEPARATOR;
  read_geometry(1, &argv);
  
  DLIList<Body*> old_bodies, new_bodies, junk;
  gti->bodies(old_bodies);
  DLIList<RefEntity*> free_entities;
  gti->get_free_ref_entities(free_entities);

  DLIList<RefEdge*> curves;
  for(int i = 0; i < free_entities.size(); i++)
  {
    RefEdge* free_edge = CAST_TO(free_entities.get_and_step(), RefEdge);
    curves.append(free_edge);
  }

  old_bodies.reset();

  //use curve 2 as the sweep along curve
  DLIList<RefEdge*> temp_curves;
  old_bodies.get()->ref_edges(temp_curves);
  RefEdge * edge_to_sweep_along = NULL;
  // This model happen to be the forth edge wwhich is suitable to sweep along
  for (int i = 0; i < 4; i ++)
  {
     edge_to_sweep_along = temp_curves.get_and_step();
  } 
        
  CubitVector axis;
  axis.set(1.,0.,0.);

  CubitStatus rsl= gmti->webcut_with_sweep_curves(old_bodies,curves,axis,true,NULL,
                                                  edge_to_sweep_along,new_bodies, junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_sweep_along_curve.sat";
  const char * filetype = "ACIS_SAT";
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->delete_geometry();
  return rsl;
}

CubitStatus webcut_with_curve_loop()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read in the geometry from files specified on the command line
  const char *argv = STRINGIFY(SRCDIR) "/model4.sat";
  PRINT_SEPARATOR;
  read_geometry(1, &argv);
  
  DLIList<Body*> old_bodies, new_bodies, junk;
  gti->bodies(old_bodies);
  DLIList<RefEntity*> free_entities;
  gti->get_free_ref_entities(free_entities);

  DLIList<RefEdge*> curves;
  for(int i = 0; i < free_entities.size(); i++)
  {
    RefEdge* free_edge = CAST_TO(free_entities.get_and_step(), RefEdge);
    curves.append(free_edge);
  }

  old_bodies.reset();

  CubitStatus rsl= gmti->webcut_with_curve_loop(old_bodies,curves,new_bodies, junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_curve_loop.sat";
  const char * filetype = "ACIS_SAT";
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->delete_geometry();
  return rsl;
}

CubitStatus webcut_with_extended_surf()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read in the geometry from files specified on the command line
  const char *argv = STRINGIFY(SRCDIR) "/model5.sat";
  PRINT_SEPARATOR;
  read_geometry(1, &argv);

  DLIList<Body*> old_bodies, new_bodies, junk;
  gti->bodies(old_bodies);

  old_bodies.reset();
  old_bodies.step();

  DLIList<RefFace*> ref_faces;
  old_bodies.remove()->ref_faces(ref_faces);
  RefFace *refface = ref_faces.get();
  ref_faces.clean_out();
  ref_faces.append(refface);
  int num_cut = 0;
  PRINT_INFO("Webcut body %d with surface extended from %d\n", old_bodies.get()->id(), refface->id() );
  CubitStatus rsl= gmti->webcut_with_extended_sheet(old_bodies,ref_faces,new_bodies, junk, num_cut);
  if (rsl== CUBIT_FAILURE) {
     gti->delete_geometry();
     return rsl;
  }

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_extended_surf.sat";
  const char * filetype = "ACIS_SAT";
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->delete_geometry();
  return rsl;
}

CubitStatus webcut_with_sweep_surfaces_rotated()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read in the geometry from files specified on the command line
  const char *argv = STRINGIFY(SRCDIR) "/model7.sat";
  PRINT_SEPARATOR;
  read_geometry(1, &argv);
  
  DLIList<Body*> old_bodies, new_bodies, junk;
  gti->bodies(old_bodies);

  old_bodies.reset();
  old_bodies.step();

  DLIList<RefFace*> faces;
  old_bodies.remove()->ref_faces(faces);

  CubitVector center(4.0,14.0,14.0);
  CubitVector axis(0.,-1.,0.);

  // set 7th parameter to be false to indicate of to_next_surf = false
  CubitStatus rsl= gmti->webcut_with_sweep_surfaces_rotated(old_bodies,faces,center,axis,1.6,NULL,true,new_bodies, junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_sweep_surfaces_rotated.sat";
  const char * filetype = "ACIS_SAT";
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->delete_geometry();
  return rsl;
}

CubitStatus webcut_with_sweep_surfaces()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read file containing one block and one sheet body
  PRINT_SEPARATOR;
  Body *vol, *sheet;
  read_model6( vol, sheet );
  
  DLIList<Body*> old_bodies, new_bodies, junk;
  old_bodies.append( vol );
  DLIList<RefFace*> faces;
  sheet->ref_faces( faces );

  CubitVector axis(1.,0.,0.);

  CubitStatus rsl= gmti->webcut_with_sweep_surfaces(old_bodies,faces,axis,false,true,false, false,NULL,NULL,new_bodies, junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_sweep_surfaces.sat";
  const char * filetype = "ACIS_SAT";
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->delete_geometry();
  return rsl;
}

CubitStatus webcut_with_sweep_surfaces_along_curve()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read file containing one block and one sheet body
  PRINT_SEPARATOR;
  Body *vol, *sheet;
  read_model6( vol, sheet );
  
  DLIList<Body*> old_bodies, new_bodies, junk;
  old_bodies.append( vol );
  DLIList<RefFace*> faces;
  sheet->ref_faces( faces );

  //use curve 2 as the sweep along curve
  DLIList<RefEdge*> temp_curves;
  old_bodies.get()->ref_edges(temp_curves);
  RefEdge * edge_to_sweep_along = NULL;
  //This model happen to be the second curve which is suitable for sweepng along
  for (int i = 0; i < 2; i ++)
  {
     edge_to_sweep_along = temp_curves.get_and_step();
  }

  CubitVector axis(1.,0.,0.);

  CubitStatus rsl= gmti->webcut_with_sweep_surfaces(old_bodies,faces,axis,false, true,false,false,NULL,
                                                  edge_to_sweep_along,new_bodies, junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_sweep_surfaces_along_curve.sat";
  const char * filetype = "ACIS_SAT";
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->delete_geometry();
  return rsl;
}

CubitStatus webcut_with_sweep_surfaces_perp()
{
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  GeometryModifyTool *gmti = GeometryModifyTool::instance();

  // Read file containing one block and one sheet body
  PRINT_SEPARATOR;
  Body *vol, *sheet;
  read_model6( vol, sheet );
  
  DLIList<Body*> old_bodies, new_bodies, junk;
  old_bodies.append( vol );
  DLIList<RefFace*> faces;
  sheet->ref_faces( faces );

  CubitVector axis(1.,0.,0.);

  CubitStatus rsl= gmti->webcut_with_sweep_surfaces(old_bodies,faces,axis,true,true,false, false,NULL,NULL,new_bodies,junk);
  if (rsl== CUBIT_FAILURE)
     return rsl;

  DLIList<RefEntity*> ref_entity_list;
  int num_ents_exported=0;
  const CubitString cubit_version="10.2";
  const char * filename = "webcut_with_sweep_surfaces_perp.sat";
  const char * filetype = "ACIS_SAT";
  rsl = CubitCompat_export_solid_model(ref_entity_list, filename, filetype,
                                 num_ents_exported, cubit_version);

  gti->delete_geometry();
  return rsl;
}

