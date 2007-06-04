//-------------------------------------------------------------------------
// Filename      : FacetBody.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 7/18/00
//
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********


// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitPoint.hpp"
#include "CastTo.hpp"
#include "BodySM.hpp"
#include "Body.hpp"
#include "FacetBody.hpp"
#include "CubitSimpleAttrib.hpp"
#include "FacetQueryEngine.hpp"
#include "DLIList.hpp"
#include "FacetEvalTool.hpp"
#include "Surface.hpp"
#include "FacetSurface.hpp"
#include "CubitTransformMatrix.hpp"
#include "FacetPoint.hpp"
#include "FacetCurve.hpp"
#include "FacetCoEdge.hpp"
#include "FacetLoop.hpp"
#include "FacetShell.hpp"
#include "FacetLump.hpp"
#include "CubitFacetEdge.hpp"
#include "FacetModifyEngine.hpp"
#include "FacetAttrib.hpp"

//-------------------------------------------------------------------------
// Purpose       : A constructor with a list of lumps that are attached.
//
// Special Notes :
//
//-------------------------------------------------------------------------
FacetBody::FacetBody(DLIList<Lump*>& my_lumps)
{
  myLumps += my_lumps;
}
FacetBody::~FacetBody() 
{
    //Not sure what to do..
}

GeometryQueryEngine* FacetBody::get_geometry_query_engine() const
{
  return FacetQueryEngine::instance();
}

void FacetBody::append_simple_attribute_virt(CubitSimpleAttrib *csa)
  { attribSet.append_attribute(csa); }
  
void FacetBody::remove_simple_attribute_virt(CubitSimpleAttrib *csa)
  { attribSet.remove_attribute(csa); }

void FacetBody::remove_all_simple_attribute_virt()
  { attribSet.remove_all_attributes(); }
  
CubitStatus FacetBody::get_simple_attribute(DLIList<CubitSimpleAttrib*>& csa_list)
  { return attribSet.get_attributes(csa_list); }

CubitStatus FacetBody::get_simple_attribute( const CubitString& name,
                                          DLIList<CubitSimpleAttrib*>& csa_list )
  { return attribSet.get_attributes( name, csa_list ); }

CubitStatus FacetBody::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes( file_ptr); }

CubitStatus FacetBody::restore_attribs( FILE *file_ptr, unsigned int endian )
  { return attribSet.restore_attributes( file_ptr, endian ); }



//----------------------------------------------------------------
// Function: copy
// Description: create a new copy of the body.
// Author: sjowen
//----------------------------------------------------------------
BodySM* FacetBody::copy()
{
  CubitStatus rv;

  // ----------Copy the points on the body------------------------

  int ii;
  DLIList<FacetPoint*> point_list;
  get_points(point_list);
  DLIList<Point*> copy_points;
  point_list.reset();
  Point *point_ptr;
  Point *point_copy; 
  for(ii=0; ii<point_list.size(); ii++)
  {
    point_ptr = point_list.get_and_step();
    CubitVector temp_vector = point_ptr->coordinates();
    
    rv = FacetModifyEngine::instance()->make_facet_point(temp_vector,
                                                           point_copy);
    if(rv != CUBIT_SUCCESS)
    {
      PRINT_ERROR("Couldn't copy points");
      return (BodySM *)NULL;
    }
    copy_points.append( point_copy );
  }

  // ------------------Copy the curves-------------------------

  int jj;
  DLIList<FacetCurve*> curve_list;
  get_curves( curve_list );
  DLIList<Curve*> copy_curves;
  curve_list.reset();
  Curve *curve_ptr, *curve_copy;
  Curve *curvsm_copy = NULL;
  FacetCurve *fcurve;
  Point *ptsm_ptr;
  Point *start_ptr, *end_ptr, *copy_start = NULL, *copy_end = NULL;
  for (ii=0; ii<curve_list.size(); ii++)
  {
    curve_ptr = curve_list.get_and_step();
    fcurve = CAST_TO( curve_ptr, FacetCurve );
    start_ptr = fcurve->start_point();
    end_ptr = fcurve->end_point();
    int found0 = 0;
    int found1 = 0;

    // find the end points

    point_list.reset();
    copy_points.reset();
    for (jj=0; jj<point_list.size() && (!found0 || !found1); jj++)
    {
      point_ptr = point_list.get_and_step();
      ptsm_ptr = CAST_TO(point_ptr, Point);
      point_copy = copy_points.get_and_step();
      if (ptsm_ptr == start_ptr)
      {
        copy_start = point_copy;
        found0 = 1;
      }
      if (ptsm_ptr == end_ptr)
      {
        copy_end = point_copy;
        found1 = 1;
      }
    }

    // create the new curve and update the points

    rv = FacetModifyEngine::instance()->make_facet_curve(copy_start,
                                                           copy_end, 
                                                           curve_copy);
    if(rv != CUBIT_SUCCESS)
    {
      PRINT_ERROR("Couldn't copy curves");
      return (BodySM *)NULL;
    }
    curvsm_copy = CAST_TO(curve_copy, Curve);
    copy_curves.append( curvsm_copy );
  }

  // ------------------copy coedges-----------------------

  DLIList<FacetCoEdge*> coedge_list;
  get_coedges( coedge_list );
  DLIList<CoEdgeSM*> copy_coedges;
  coedge_list.reset();
  Curve *curvsm_ptr;
  CoEdgeSM *coedge_ptr, *coedge_copy;
  FacetCoEdge *fcoedge;
  for (ii=0; ii<coedge_list.size(); ii++)
  {
    coedge_ptr = coedge_list.get_and_step();
    fcoedge = CAST_TO( coedge_ptr, FacetCoEdge );
    Curve *curve_at_coedge = fcoedge->curve();
    int found = 0;

    // find the associated curve

    curve_list.reset();
    copy_curves.reset();
    for (jj=0; jj<curve_list.size() && !found; jj++)
    {
      curve_ptr = curve_list.get_and_step();
      curvsm_ptr = CAST_TO(curve_ptr, Curve);
      curvsm_copy = copy_curves.get_and_step();
      if (curve_at_coedge == curvsm_ptr)
      {
        found = 1;
      }
    }

    // create the new coedge

    CubitSense sense = fcoedge->sense();
    rv = FacetModifyEngine::instance()->make_facet_coedge(curvsm_copy,
                                                    sense, coedge_copy);
    if(rv != CUBIT_SUCCESS)
    {
      PRINT_ERROR("Couldn't copy coedge");
      return (BodySM *)NULL;
    }
    copy_coedges.append( coedge_copy );
  }

  // ----------------------copy loops--------------------------

  int kk;
  DLIList<FacetLoop*> loop_list;
  get_loops( loop_list );
  DLIList<LoopSM*> copy_loops;
  loop_list.reset();
  LoopSM *loop_ptr, *loop_copy;
  FacetLoop *floop;
  for (ii=0; ii<loop_list.size(); ii++)
  {
    floop = loop_list.get_and_step();
    DLIList<FacetCoEdge *>coedges_on_loop;
    floop->get_coedges(coedges_on_loop);
    DLIList<CoEdgeSM *>copy_coedges_on_loop;

    // find all associated coedges on the loop

    for(kk=0; kk<coedges_on_loop.size(); kk++)
    {
      int found = 0;
      coedge_list.reset();
      copy_coedges.reset();
      CoEdgeSM *coedge_on_loop = coedges_on_loop.get_and_step();
      for (jj=0; jj<coedge_list.size() && !found; jj++)
      {
        coedge_ptr = coedge_list.get_and_step();
        coedge_copy = copy_coedges.get_and_step();     
        if (coedge_on_loop == coedge_ptr)
        {
          found = 1;
          copy_coedges_on_loop.append(coedge_copy);
        }
      }
    }

    // create the new loop

    rv = FacetModifyEngine::instance()->make_facet_loop(copy_coedges_on_loop,
                                                          loop_copy);
    if(rv != CUBIT_SUCCESS)
    {
      PRINT_ERROR("Couldn't copy loops");
      return (BodySM *)NULL;
    }
    copy_loops.append( loop_copy );
  }

  // ----------------------copy surfaces--------------------------

  DLIList<FacetSurface*> surface_list;
  get_surfaces(surface_list);
  DLIList<Surface*> copy_surfaces;
  surface_list.reset();
  Surface *surface_ptr, *surface_copy;
  FacetSurface *fsurface;
  for (ii=0; ii<surface_list.size(); ii++)
  {
    fsurface = surface_list.get_and_step();
    DLIList<FacetLoop *>loops_on_surface;
    fsurface->get_loops(loops_on_surface);
    DLIList<LoopSM *>copy_loops_on_surface;

    // find all associated loops on the surface

    for(kk=0; kk<loops_on_surface.size(); kk++)
    {
      int found = 0;
      loop_list.reset();
      copy_loops.reset();
      LoopSM *loop_on_surface = loops_on_surface.get_and_step();
      for (jj=0; jj<loop_list.size() && !found; jj++)
      {
        loop_ptr = loop_list.get_and_step();
        loop_copy = copy_loops.get_and_step();     
        if (loop_on_surface == loop_ptr)
        {
          found = 1;
          copy_loops_on_surface.append(loop_copy);
        }
      }
    }

    // create the new surface

    DLIList<CubitFacet*>facet_list;
    DLIList<CubitPoint*>cpoint_list;
    rv = fsurface->copy_facets( facet_list, cpoint_list );
    if (rv != CUBIT_SUCCESS)
    {
      return (BodySM *)NULL;
    }
    int interp_order = fsurface->interp_order();
    double min_dot = fsurface->min_dot();
    const CubitEvaluatorData *eval_data = fsurface->evaluator_data();
    CubitBoolean use_point_addresses = CUBIT_FALSE;
    rv = FacetModifyEngine::instance()->make_facet_surface(eval_data,
                                                           facet_list,
                                                           cpoint_list, 
                                                           copy_loops_on_surface,
                                                           interp_order,
                                                           min_dot,
                                                           surface_copy, 
                                                           use_point_addresses);
    if(rv != CUBIT_SUCCESS)
    {
      PRINT_ERROR("Couldn't copy surfaces");
      return (BodySM *)NULL;
    }
    copy_surfaces.append( surface_copy );
  }

  // ----------------------copy shells--------------------------

  DLIList<FacetShell*> shell_list;
  get_shells(shell_list);
  DLIList<ShellSM*> copy_shells;
  shell_list.reset();
  ShellSM *shell_ptr, *shell_copy;
  FacetShell *fshell;
  for (ii=0; ii<shell_list.size(); ii++)
  {
    fshell = shell_list.get_and_step();
    DLIList<FacetSurface *>surfaces_on_shell;
    fshell->get_surfaces(surfaces_on_shell);
    DLIList<Surface *>copy_surfaces_on_shell;

    // find all associated loops on the surface

    for(kk=0; kk<surfaces_on_shell.size(); kk++)
    {
      int found = 0;
      surface_list.reset();
      copy_surfaces.reset();
      Surface *surface_on_shell = surfaces_on_shell.get_and_step();
      for (jj=0; jj<surface_list.size() && !found; jj++)
      {
        surface_ptr = surface_list.get_and_step();
        surface_copy = copy_surfaces.get_and_step();     
        if (surface_on_shell == surface_ptr)
        {
          found = 1;
          copy_surfaces_on_shell.append(surface_copy);
        }
      }
    }

    // create the new shell

    rv = FacetModifyEngine::instance()->make_facet_shell(copy_surfaces_on_shell,
                                                           shell_copy);
    if(rv != CUBIT_SUCCESS)
    {
      PRINT_ERROR("Couldn't copy shell");
      return (BodySM *)NULL;
    }

    // set the sense of the surfaces on the shell - copy from the old shell

    FacetShell *fshell_copy = CAST_TO(shell_copy, FacetShell);
    surfaces_on_shell.reset();
    copy_surfaces_on_shell.reset();
    for (kk=0; kk<surfaces_on_shell.size(); kk++)
    {
      Surface *surface_on_shell = surfaces_on_shell.get_and_step();
      Surface *copy_surface_on_shell = copy_surfaces_on_shell.get_and_step();
      fsurface = CAST_TO( surface_on_shell, FacetSurface );
      CubitSense sense = fsurface->get_shell_sense(fshell);
      FacetSurface *copy_fsurface = CAST_TO( copy_surface_on_shell, FacetSurface );
      copy_fsurface->set_shell_sense( fshell_copy, sense ); 
    }
    copy_shells.append( shell_copy );
  }

  // ----------------------copy lumps--------------------------

  DLIList<FacetLump*> lump_list;
  get_lumps(lump_list);
  DLIList<Lump*> copy_lumps;
  lump_list.reset();
  Lump *lump_copy;
  FacetLump *flump;
  for (ii=0; ii<lump_list.size(); ii++)
  {
    flump = lump_list.get_and_step();
    DLIList<FacetShell *>shells_on_lump;
    flump->get_shells(shells_on_lump);
    DLIList<ShellSM *>copy_shells_on_lump;

    // find all associated loops on the surface

    for(kk=0; kk<shells_on_lump.size(); kk++)
    {
      int found = 0;
      shell_list.reset();
      copy_shells.reset();
      ShellSM *shell_on_lump = shells_on_lump.get_and_step();
      for (jj=0; jj<shell_list.size() && !found; jj++)
      {
        shell_ptr = shell_list.get_and_step();
        shell_copy = copy_shells.get_and_step();     
        if (shell_on_lump == shell_ptr)
        {
          found = 1;
          copy_shells_on_lump.append(shell_copy);
        }
      }
    }

    // create the new lump

    rv = FacetModifyEngine::instance()->make_facet_lump(copy_shells_on_lump,
                                                          lump_copy);
    if(rv != CUBIT_SUCCESS)
    {
      PRINT_ERROR("Couldn't copy lump");
      return (BodySM *)NULL;
    }
    copy_lumps.append( lump_copy );
  }

  // ----------------------copy body--------------------------

  BodySM *body_copy;
  rv = FacetModifyEngine::instance()->make_facet_body(copy_lumps,
                                                        body_copy);
  if(rv != CUBIT_SUCCESS)
  {
    PRINT_ERROR("Couldn't copy lump");
    return (BodySM *)NULL;
  }

  return (BodySM*)body_copy;
}
//---------------------------------------------------------------- 
// Function: can_be_deleted 
// Description: determine if the body can be deleted 
// 
// Author: sjowen 
//---------------------------------------------------------------- 
CubitBoolean FacetBody::can_be_deleted( DLIList <Body*> &body_list ) 
{ 
  CubitBoolean delete_ok = CUBIT_TRUE; 
  DLIList<FacetSurface *>surf_list; 
  get_surfaces(surf_list); 
  int ii; 
  for (ii=0; ii<surf_list.size() && delete_ok; ii++) 
  { 
    FacetSurface *surf_ptr = surf_list.get_and_step(); 
    DLIList<FacetBody*>my_body_list; 
    surf_ptr->get_bodies(my_body_list); 
    int jj; 
    if (my_body_list.size() >= 2) 
    { 
      for (jj=0; jj<my_body_list.size() && delete_ok; jj++) 
      { 
        BodySM *my_body_ptr = my_body_list.get_and_step(); 
        if (my_body_ptr != this) 
        { 
          int kk; 
          int found = 0; 
          for (kk=0; kk<body_list.size() && !found; kk++) 
          { 
            Body *body_ptr = body_list.get_and_step(); 
            FacetBody* fbody_ptr = CAST_TO(body_ptr->get_body_sm_ptr(), FacetBody); 
            if (fbody_ptr) 
            { 
              if (my_body_ptr == fbody_ptr) 
                found = 1; 
            } 
          } 
          if (!found) 
          { 
            delete_ok = CUBIT_FALSE; 
            PRINT_ERROR("Body cannot be deleted because it is merged with adjacent Body\n"); 
            PRINT_INFO("    Mesh Based Geometry entities cannot be unmerged.\n" 
              "    Try using the no_merge option when importing the mesh\n"); 
          } 
        } 
      } 
    } 
  } 
  return delete_ok; 
} 
    
//----------------------------------------------------------------
// Function: move
// Description: translate the body and its child entities
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus FacetBody::move(double dx, double dy, double dz)
{
  CubitTransformMatrix tfmat;
  tfmat.translate( dx, dy, dz );

  CubitStatus stat = transform( tfmat, CUBIT_FALSE );

  if (stat == CUBIT_SUCCESS)
    myTransforms.translate( dx, dy, dz );

  return stat;
}


//----------------------------------------------------------------
// Function: rotate
// Description: rotate the body and its child entities
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus FacetBody::rotate( double x, double y, double z, 
                               double angle_in_degrees )
{

  CubitTransformMatrix rotmat;
  CubitVector axis( x, y, z );
  rotmat.rotate( angle_in_degrees, axis );

  CubitStatus stat = transform( rotmat, CUBIT_TRUE );

  if (stat == CUBIT_SUCCESS)
    myTransforms.rotate( angle_in_degrees, axis );

  return stat;
}

//----------------------------------------------------------------
// Function: scale
// Description: scale the body and its child entities
//              use a constant scale factor
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus FacetBody::scale(double scale_factor )
{
  return scale(scale_factor,scale_factor,scale_factor);
}

//----------------------------------------------------------------
// Function: scale
// Description: scale the body and its child entities
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus FacetBody::scale(double scale_factor_x,
                             double scale_factor_y,
                             double scale_factor_z )
{
  CubitTransformMatrix scalemat;
  scalemat.scale_about_origin( scale_factor_x, 
                               scale_factor_y, 
                               scale_factor_z );

  CubitStatus stat = transform( scalemat, CUBIT_FALSE );

  if (stat == CUBIT_SUCCESS)
    myTransforms.scale_about_origin( scale_factor_x, 
                                     scale_factor_y, 
                                     scale_factor_z );

  // scale the facetcurve

  DLIList<FacetCurve *> curve_list;
  get_curves(curve_list); 
  Curve *curv_ptr;
  for (int ii=0; ii<curve_list.size(); ii++)
  {
    curv_ptr = curve_list.get_and_step();
    FacetCurve *fcurve = CAST_TO( curv_ptr, FacetCurve );
    if (fcurve)
    {
      fcurve->reset_length();
    }
  }

  return stat;
}

//----------------------------------------------------------------
// Function: restore
// Description: restore the body and its child entities
//              to its original coordinates using the inverse
//              transformation matrix
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus FacetBody::restore()
{
  // invert the transformation matrix and apply to entities 
  // (assumes an orthogonal matrix (ie. no shear or non-uniform scaling)

  CubitTransformMatrix inverse_mat;
  inverse_mat = myTransforms.inverse();

  CubitStatus stat = transform( inverse_mat, CUBIT_TRUE );

  if (stat == CUBIT_SUCCESS)
    myTransforms.set_to_identity();

  return stat;
}

//----------------------------------------------------------------
// Function: reflect
// Description: reflect the body about a exis
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus FacetBody::reflect( double reflect_axis_x,
                                double reflect_axis_y,
                                double reflect_axis_z )
{
  CubitTransformMatrix reflectmat;
  CubitVector reflect_vector( reflect_axis_x, 
                              reflect_axis_y, 
                              reflect_axis_z );
  reflectmat.reflect( reflect_vector );

  CubitStatus stat = transform( reflectmat, CUBIT_TRUE );

  if (stat == CUBIT_SUCCESS)
    myTransforms.reflect( reflect_vector );

  return stat;
}

//----------------------------------------------------------------
// Function: transform
// Description: transform the body based on a transformation matrix
//              main function for applying transformations to 
//              facet-based bodies
//
// Author: sjowen
//----------------------------------------------------------------
CubitStatus FacetBody::transform( CubitTransformMatrix &tfmat, 
                                  CubitBoolean is_rotation )
{
  int ii;

  // get the list of surfaces on the body

  DLIList<FacetSurface *>surf_list;
  get_surfaces( surf_list );
  Surface *surf;
  FacetSurface *fsurf;
  FacetEvalTool *ftool;
#ifdef BOYD17 
  CubitBox box;
#endif
  //CubitVector min, max;

  // go through all the surfaces and collect the list of all points.
  // (some may be listed on multiple surfaces)

  DLIList<CubitPoint *>point_list;
  for (ii=0; ii<surf_list.size(); ii++)
  {
    surf = surf_list.get_and_step();
    fsurf = CAST_TO( surf, FacetSurface );
    fsurf->get_my_points( point_list );
  }

  // unmark all the points so we can keep track of the ones that have
  // already been transformed

  CubitPoint *cp;
  for (ii=0; ii<point_list.size(); ii++)
  {
    cp = point_list.get_and_step();
    cp->marked( 0 );
  }

  // transform the points

  //CubitVector norm, du, dv;
  for (ii=0; ii<point_list.size(); ii++)
  {
    cp = point_list.get_and_step();
    if (!cp->marked())
    {
      cp->transform( tfmat );
      if (is_rotation)
        cp->rotate_normal( tfmat );
      cp->marked( 1 );   
    }
  }

  // check the vertices - make sure they are transformed

  FacetPoint *fpt;
  Point *pt;
  DLIList<FacetPoint*>gpoint_list;
  get_points(gpoint_list);
  for (ii=0; ii<gpoint_list.size(); ii++)
  {
    pt = gpoint_list.get_and_step();
    fpt = CAST_TO( pt, FacetPoint );

    // only transform the point if it isn't already part of the facets
    // (they could be points by themselves)

    cp = fpt->get_cubit_point();
    if (cp->num_adj_facets() == 0)
    {
      cp->transform( tfmat );
      if (is_rotation)
        cp->rotate_normal( tfmat );
    }
  }

  // reset the bounding box and update the facet normal and plane

  // init flags on edges to 0
  DLIList<Surface*> tmp_surf_list( surf_list.size() );
  CAST_LIST_TO_PARENT( surf_list, tmp_surf_list );
  init_edge_flags( tmp_surf_list, 0 );
  for (ii=0; ii<surf_list.size(); ii++)
  {
    surf = surf_list.get_and_step();
    fsurf = CAST_TO( surf, FacetSurface );

    // if we are using a bspline representation, then we also need to 
    // transform the control points on the edges and facets

    ftool = fsurf->get_eval_tool();
    if (ftool->interp_order() == 4)
    {
      ftool->transform_control_points( tfmat );
    }

    DLIList<CubitFacet *>flist;
    DLIList<CubitPoint *>plist;
    fsurf->get_my_facets( flist, plist);
    int jj;
    CubitFacet *facet_ptr;
    for (jj=0; jj<flist.size(); jj++)
    {
      facet_ptr = flist.get_and_step();
      facet_ptr->update_plane();
      facet_ptr->reset_bounding_box();
    }
      //must be done after facet's are reset...
    ftool->reset_bounding_box();
    // if this facet surface has a primitive evaluator, then we need
    // to tell it about the transformation also.
    fsurf->add_transformation( tfmat );

      //re-calculate the area of the surface in case it changed
    fsurf->update_measurement();
  }
  init_edge_flags( tmp_surf_list, 0 );

    // Some transforms (those incorporating reflections)
    // invert the geometry.  Correct for it.
    // -- jason k.
  if ( tfmat.sub_matrix(3,3).determinant() < 0.0 )
  {
      // Flip CoFace senses
    DLIList<FacetShell*> shells;
    get_shells( shells );
    //modified.  mbrewer.  doing both a reverse and a 
    //reverse_surfaces.  the latter actually changes the 
    // underlying surfaces so that the normals can all
    // still be outward pointing.  It automatically changes 
    // the sense, therefore we also still need excplicity
    // call reverse so that sense is corrected.
    while (shells.size()){
      shells.get()->reverse_surfaces();
      shells.pop()->reverse();
    }
  }

  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Function: init_edge_flags
// Description: set the flags on the facet edges
// Note:  Only done for facet eval tool with order 4 interpolation
// Author: sjowen
//----------------------------------------------------------------
void FacetBody::init_edge_flags( DLIList<Surface *>&surf_list, 
                                 int )
{
  int ii, jj;
  Surface *surf;
  FacetSurface *fsurf;
  FacetEvalTool *ftool;
  CubitFacetEdge *edge_ptr;

  for (ii=0; ii<surf_list.size(); ii++)
  {  
    DLIList<CubitFacetEdge*>edge_list;
    surf = surf_list.get_and_step();
    fsurf = CAST_TO( surf, FacetSurface );
    ftool = fsurf->get_eval_tool();
    if (ftool->interp_order() == 4)
    {
      ftool->get_edges( edge_list );
      for (jj=0; jj<edge_list.size(); jj++)
      {
        edge_ptr = edge_list.get_and_step();
        edge_ptr->set_flag( 0 );
      }
    }
  }
}

CubitStatus FacetBody::get_transforms( CubitTransformMatrix &tfm ) 
{
  tfm = myTransforms;
  return CUBIT_SUCCESS;
}

CubitStatus FacetBody::set_transforms( CubitTransformMatrix tfm ) 
{
  myTransforms = tfm;
  return CUBIT_SUCCESS;
}

int FacetBody::validate(const CubitString &, DLIList <TopologyEntity*>&)
{
  PRINT_ERROR("This option is not available for mesh defined geometry.\n");
  return 0;
}

void FacetBody::get_parents_virt( DLIList<TopologyBridge*>& ) 
  {}
  
void FacetBody::get_children_virt( DLIList<TopologyBridge*>& lumps )
  { CAST_LIST_TO_PARENT( myLumps, lumps ); }

void FacetBody::get_lumps( DLIList<FacetLump*>& result_list )
{
  myLumps.reset();
  for ( int i = 0; i < myLumps.size(); i++ )
    if ( FacetLump* lump = dynamic_cast<FacetLump*>(myLumps.next(i)) )
      result_list.append(lump);
}

void FacetBody::get_shells( DLIList<FacetShell*>& result_list )
{
  DLIList<FacetLump*> lump_list;
  get_lumps( lump_list );
  lump_list.reset();
  for ( int i = 0; i < lump_list.size(); i++ )
    lump_list.next(i)->get_shells( result_list );
}

void FacetBody::get_surfaces( DLIList<FacetSurface*>& result_list )
{
  DLIList<FacetShell*> shell_list;
  DLIList<FacetSurface*> tmp_list;
  get_shells(shell_list);
  shell_list.reset();
  for ( int i = 0; i < shell_list.size(); i++ )
  {
    tmp_list.clean_out();
    shell_list.next(i)->get_surfaces( tmp_list );
    result_list.merge_unique( tmp_list );
  }
}

void FacetBody::get_loops( DLIList<FacetLoop*>& result_list )
{
  DLIList<FacetSurface*> surface_list;
  get_surfaces( surface_list );
  surface_list.reset();
  for ( int i = 0; i < surface_list.size(); i++ )
    surface_list.next(i)->get_loops( result_list );
}

void FacetBody::get_coedges( DLIList<FacetCoEdge*>& result_list )
{
  DLIList<FacetSurface*> surface_list;
  get_surfaces( surface_list );
  surface_list.reset();
  for ( int i = 0; i < surface_list.size(); i++ )
    surface_list.next(i)->get_coedges( result_list );
}

void FacetBody::get_curves( DLIList<FacetCurve*>& result_list )
{
  DLIList<FacetCoEdge*> coedge_list;
  get_coedges( coedge_list );
  coedge_list.reset();
  for ( int i = coedge_list.size(); i--; )
  {
    FacetCoEdge* coedge = coedge_list.get_and_step();
    FacetCurve* curve = dynamic_cast<FacetCurve*>(coedge->curve());
    if (curve)
      result_list.append_unique(curve);
  }
}

void FacetBody::get_points( DLIList<FacetPoint*>& result_list )
{
  DLIList<FacetCurve*> curve_list;
  get_curves( curve_list );
  curve_list.reset();
  for ( int i = curve_list.size(); i--; )
  {
    FacetCurve* curve = curve_list.get_and_step();
    FacetPoint* point = dynamic_cast<FacetPoint*>(curve->start_point());
    if (point)
      result_list.append_unique(point);
    point = dynamic_cast<FacetPoint*>(curve->end_point());
    if (point)
      result_list.append_unique(point);
  }
}

void FacetBody::add_lump( FacetLump *lump_to_add )
{
  Lump* lump = dynamic_cast<Lump*>(lump_to_add);
  if (lump)
  {
    lump_to_add->add_body(this);
    myLumps.append( lump );
  }
}

void FacetBody::remove_lump( FacetLump *lump_to_remove )
{
  FacetLump* lump = dynamic_cast<FacetLump*>(lump_to_remove);
  if (lump)
  {
    assert(lump_to_remove->get_body() == this);
    lump_to_remove->remove_body();
    myLumps.remove( lump );
  }
}


//-------------------------------------------------------------------------
// Purpose       : Tear down topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
void FacetBody::disconnect_all_lumps()
{
  myLumps.reset();
  for (int i = myLumps.size(); i--; )
  {
    Lump* sm_ptr = myLumps.get_and_step();
    FacetLump* lump = dynamic_cast<FacetLump*>(sm_ptr);
    if (lump)
    {
      assert(lump->get_body() == this);
      lump->remove_body();
    }
  }
  myLumps.clean_out();
}

//-------------------------------------------------------------------------
// Purpose       : Find centroid
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitStatus FacetBody::mass_properties( CubitVector& centroid, 
                                        double& volume )
{
  centroid.set( 0.0, 0.0, 0.0 );
  volume = 0.0;
  
  DLIList<FacetLump*> lumps (myLumps.size());
  CAST_LIST( myLumps, lumps, FacetLump );
  assert( myLumps.size() == lumps.size() );
  for (int i = lumps.size(); i--; )
  {
    CubitVector cent;
    double vol;
    if (CUBIT_SUCCESS != lumps.get_and_step()->mass_properties(cent,vol))
      return CUBIT_FAILURE;
    centroid += vol*cent;
    volume += vol;
  }
  if (volume > CUBIT_RESABS)
  {
    centroid /= volume;
  }
  else
  {
    centroid.set( 0.0, 0.0, 0.0 );
    volume = 0.0;
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Used to be FacetQueryEngine::is_point_in_body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitPointContainment FacetBody::point_containment( const CubitVector &point )
{
  CubitPointContainment pc_value; 
  FacetLump *facet_lump;

  int i;
  for(i=myLumps.size(); i--;)
  {
    facet_lump = dynamic_cast<FacetLump*>(myLumps.get_and_step()); 
    pc_value = facet_lump->point_containment( point );
    if( pc_value == CUBIT_PNT_INSIDE )
      return CUBIT_PNT_INSIDE;
    else if( pc_value == CUBIT_PNT_BOUNDARY )
      return CUBIT_PNT_BOUNDARY;
  }

  return CUBIT_PNT_OUTSIDE;
}


