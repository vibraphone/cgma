//- Class: SurfaceOverlapTool
//- Description: Utilities to debug imprinting/merging problems
//- Owner: Steve Storm
//- Created: 22 October 1999
//- Overhauled: January 2003

#include "RefEntityFactory.hpp"
#include "SurfaceOverlapTool.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "Body.hpp"
#include "CoEdge.hpp"
#include "BodySM.hpp"
#include "RefGroup.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "TopologyEntity.hpp"
#include "Surface.hpp"
#include "Curve.hpp"

#include "CubitBox.hpp"
#include "CubitUtil.hpp"

#include "DLIList.hpp"
#include "ProgressTool.hpp"
#include "AppUtil.hpp"
#include "SurfaceOverlapFacet.hpp"
#include "CurveOverlapFacet.hpp"
#include "TDSurfaceOverlap.hpp"
#include "RTree.hpp"
#include "AbstractTree.hpp"

#include "GMem.hpp"
#include "SettingHandler.hpp"

#include "GfxPreview.hpp"
#include "GfxDebug.hpp"
#include "CpuTimer.hpp"

#define NO_FACETS_FOR_ABSTRACTTREE 10

SurfaceOverlapTool* SurfaceOverlapTool::instance_ = 0;
double SurfaceOverlapTool::gapMin = 0.0;
double SurfaceOverlapTool::gapMax = 0.01;
double SurfaceOverlapTool::angleMin = 0.0;
double SurfaceOverlapTool::angleMax = 5.0;
int SurfaceOverlapTool::normalType = 1; // 1=any, 2=opposite, 3=same
double SurfaceOverlapTool::overlapTolerance = .001;
CubitBoolean SurfaceOverlapTool::groupResults = CUBIT_TRUE;
CubitBoolean SurfaceOverlapTool::listPairs = CUBIT_TRUE;
CubitBoolean SurfaceOverlapTool::displayPairs = CUBIT_TRUE;
CubitBoolean SurfaceOverlapTool::imprintResults = CUBIT_FALSE;
double SurfaceOverlapTool::facetAbsTol = 0.0; // Use default
unsigned short SurfaceOverlapTool::facetAngTol = 15; // Seems to work pretty good
CubitBoolean SurfaceOverlapTool::checkWithinBodies = CUBIT_FALSE;
CubitBoolean SurfaceOverlapTool::checkAcrossBodies = CUBIT_TRUE;
bool SurfaceOverlapTool::skipFacingSurfaces = CUBIT_FALSE; // skip the pair if the normals pass through eachother (normalType must == 2)

// Constructor
SurfaceOverlapTool* SurfaceOverlapTool::instance()
{
   if( instance_ == NULL )
       instance_ = new SurfaceOverlapTool;
   return instance_;
}

SurfaceOverlapTool::SurfaceOverlapTool()
{
  facetAbsTol = 0.0; // Use default
  facetAngTol = 10;  // Seems to work pretty good
  gapMin = 0.0;
  gapMax = 0.01;
  angleMin = 0.0;
  angleMax = 5.0;
  normalType = 1; // 1=any, 2=opposite, 3=same
  groupResults = CUBIT_TRUE;
  listPairs = CUBIT_TRUE;
  displayPairs = CUBIT_TRUE;
  overlapTolerance = .001;
  imprintResults = CUBIT_FALSE;
  checkWithinBodies = CUBIT_FALSE;
  checkAcrossBodies = CUBIT_TRUE;
  skipFacingSurfaces = CUBIT_FALSE;
}

// Destructor
SurfaceOverlapTool::~SurfaceOverlapTool()
{
   delete instance_;
}

CubitBoolean SurfaceOverlapTool::draw_overlapping_surface_pair( RefFace *ref_face_1,
                                                        RefFace *ref_face_2)
{
  CubitBoolean abort = CUBIT_FALSE;
  CubitBoolean draw_overlap = CUBIT_TRUE;
  CubitBoolean overlap = check_overlap( ref_face_1, ref_face_2, 
                                        abort, draw_overlap );

  return overlap;
}

CubitStatus 
SurfaceOverlapTool::find_overlapping_surfaces( DLIList<RefFace*> &ref_face_list, 
                                               DLIList<RefEntity*> &faces_to_draw,
                                               bool filter_slivers)
{
   DLIList<RefFace*> list1, list2;
   return find_overlapping_surfaces( ref_face_list, list1, list2, faces_to_draw, CUBIT_TRUE,
     filter_slivers);
}

CubitStatus
SurfaceOverlapTool::find_overlapping_surfaces( DLIList<Body*> &body_list,
                                               DLIList<RefEntity*> &faces_to_draw,
                                               bool filter_slivers)
{
  DLIList<RefFace*> list1, list2;
  return find_overlapping_surfaces( body_list, list1, list2, faces_to_draw, CUBIT_TRUE,
    filter_slivers);
}


CubitStatus
SurfaceOverlapTool::find_overlapping_surfaces( DLIList<BodySM*> &body_list,
                                              DLIList<Surface*> &surface_list1,
                                              DLIList<Surface*> &surface_list2,
                                              bool filter_slivers)
{
  //collect all the surfaces
  DLIList<Surface*> surface_list;
  int i;
  for( i=body_list.size(); i--; )
  {
    BodySM *body_sm = body_list.get_and_step();
    DLIList<Surface*> surfs;
    body_sm->surfaces_ignore_virtual( surfs, false );
//    body_sm->surfaces( surfs );
    surface_list.merge_unique( surfs );
  }

  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;

  // Populate the Surface AbstractTree
  AbstractTree<Surface*> *a_tree = new RTree<Surface*>( tolerance );
  surface_list.reset();
  for( i=surface_list.size(); i--; )
  {
    Surface *surface = surface_list.get_and_step();
    a_tree->add( surface );
  }

  std::map<Surface*, DLIList<SurfaceOverlapFacet*>* > surface_facet_map;
  std::map<Surface*, double > surface_to_area_map;
  std::map<Surface*, AbstractTree<SurfaceOverlapFacet*>* > a_tree_map;

  surface_list.reset();
  for( i=surface_list.size(); i--; )
  {
    Surface *surf1 = surface_list.get_and_step();

    BodySM *surf1_body = surf1->bodysm();

    // Remove this surface from AbstractTree so it is not found and never
    // found again
    a_tree->remove( surf1 );

    // Find RefFaces from AbstractTree that are within range of this surface
    CubitBox surf1_box = surf1->bounding_box();
    DLIList<Surface*> close_surfaces;
    a_tree->find( surf1_box, close_surfaces );

    int j;
    for( j=close_surfaces.size(); j--; )
    {
      Surface *surf2 = close_surfaces.get_and_step();
      BodySM *surf2_body = surf2->bodysm();

      //don't check for overlapping surfaces within bodies
      if( surf1_body == surf2_body )
        continue;

      if( check_overlap( surf1, surf2, 
                         &surface_facet_map, 
                         &surface_to_area_map, 
                         &a_tree_map  ) == CUBIT_TRUE )
      {
        surface_list1.append( surf1 );
        surface_list2.append( surf2 );
      }
    }
  }

  //clean up maps;
  std::map<Surface*, AbstractTree<SurfaceOverlapFacet*>* >::iterator tree_iter; 
  tree_iter = a_tree_map.begin();
  for(; tree_iter != a_tree_map.end(); tree_iter++ )
    delete tree_iter->second;

  std::map<Surface*, DLIList<SurfaceOverlapFacet*>* >::iterator sof_iter; 
  sof_iter = surface_facet_map.begin();
  for(; sof_iter != surface_facet_map.end(); sof_iter++)
  {
    DLIList<SurfaceOverlapFacet*> *tmp_list = sof_iter->second;

    //delete contents of list
    for( i=tmp_list->size(); i--; )
      delete tmp_list->get_and_step();

    //delete the list
    delete tmp_list;
  }

  delete a_tree;

  return CUBIT_SUCCESS;
}

CubitStatus
SurfaceOverlapTool::find_overlapping_surfaces( DLIList<Body*> &body_list,
                                              DLIList<RefFace*> &ref_face_list1,
                                              DLIList<RefFace*> &ref_face_list2,
                                              DLIList<RefEntity*> &faces_to_draw,
                                              CubitBoolean show_messages,
                                              bool filter_slivers)
{
  CubitStatus status = CUBIT_SUCCESS;

  CubitBoolean group_results = CUBIT_FALSE;
  CubitBoolean list_pairs = CUBIT_FALSE;
  CubitBoolean display_pairs = CUBIT_FALSE;
  CubitBoolean imprint_results = CUBIT_FALSE;

  if( show_messages == CUBIT_TRUE )
  {
     group_results = groupResults;
     list_pairs = listPairs;
     display_pairs = displayPairs;
     imprint_results = imprintResults;
  }

  // Handle the special case of finding overlapping surfaces within a given 
  // body - we can do this MUCH faster than the general case (this is a 
  // rare case in general but at Cat we have an application for this!).

  // The usual case
  if( checkWithinBodies == CUBIT_FALSE ||
      (checkWithinBodies == CUBIT_TRUE && checkAcrossBodies == CUBIT_TRUE) )
  {
    // Utilize a straight surface list
    DLIList<RefFace*> ref_face_list;
    body_list.reset();
    int i;
    for( i=body_list.size(); i--; )
    {
      Body* body_ptr = body_list.get_and_step();
      DLIList<RefFace*> body_face_list;
      body_ptr->ref_faces( body_face_list );
      ref_face_list.merge_unique( body_face_list );
    }

    int prog_step = 0;
    if( show_messages )
    {
      prog_step = 10;
      PRINT_INFO( "Finding surface overlap...\n" );
      if( ref_face_list.size() > prog_step )
      {
        char message[128];
        sprintf(message, "Finding Surface Overlap On %d Surfaces", ref_face_list.size() );
        AppUtil::instance()->progress_tool()->start(0, ref_face_list.size(),
          "Progress", message, TRUE, TRUE );
      }
    }

    status = find_overlapping_surfaces( ref_face_list, ref_face_list1,
      ref_face_list2, faces_to_draw, list_pairs, prog_step, filter_slivers );

    if( show_messages && ref_face_list.size() > prog_step )
      AppUtil::instance()->progress_tool()->end();
  }

  // Special case - checking within bodies
  else
  {
    int prog_step = 5;
    if( show_messages )
    {
      PRINT_INFO( "Finding surface overlap...\n" );
      if( body_list.size() > prog_step )
      {
        char message[128];
        sprintf(message, "Finding Surface Overlap On %d Bodies", body_list.size() );
        AppUtil::instance()->progress_tool()->start(0, body_list.size(),
          "Progress", message, TRUE, TRUE );
      }
    }

    int i;
    body_list.reset();
    for( i=body_list.size(); i--; )
    {
      Body* body_ptr = body_list.get_and_step();
      DLIList<RefFace*> body_face_list;
      body_ptr->ref_faces( body_face_list );

      status = find_overlapping_surfaces( body_face_list, ref_face_list1,
        ref_face_list2, faces_to_draw, list_pairs, -1, filter_slivers );

      if( show_messages && body_list.size() > prog_step )
        AppUtil::instance()->progress_tool()->step();

      if( status == CUBIT_FAILURE )
        break;
    }

    if( show_messages && body_list.size() > prog_step )
        AppUtil::instance()->progress_tool()->end();
  }

  if( faces_to_draw.size() )
  {
    if( group_results == CUBIT_TRUE )
    {
      RefGroup *new_group = RefEntityFactory::instance()->construct_RefGroup( "surf_overlap" );
      new_group->add_ref_entity( faces_to_draw );
      CubitString name = new_group->entity_name();
      PRINT_INFO( "Found %d overlapping surface pairs (added to group '%s')\n", 
        ref_face_list1.size(), name.c_str() );
    }
    else if( show_messages )
      PRINT_INFO( "Found %d overlapping surface pairs\n", ref_face_list1.size() );

    if ( imprint_results )
    {
      CubitStatus stat = imprint(ref_face_list1, ref_face_list2);
      if ( stat != CUBIT_SUCCESS )
      {
        PRINT_WARNING("Imprinting overlaps was unsuccessful\n");
      }
    }
  }
  else if (show_messages )
    PRINT_INFO( "Found 0 overlapping surface pairs\n" );

  return status;
}

CubitStatus
SurfaceOverlapTool::find_overlapping_surfaces( DLIList<RefFace*> &ref_face_list,
                                              DLIList<RefFace*> &ref_face_list1,
                                              DLIList<RefFace*> &ref_face_list2,
                                              DLIList<RefEntity*> &faces_to_draw,
                                              CubitBoolean show_messages, 
                                              bool filter_slivers)
{
  CubitStatus status = CUBIT_SUCCESS;

  CubitBoolean group_results = CUBIT_FALSE;
  CubitBoolean list_pairs = CUBIT_FALSE;
  CubitBoolean display_pairs = CUBIT_FALSE;
  CubitBoolean imprint_results = CUBIT_FALSE;

  if( show_messages == CUBIT_TRUE )
  {
     group_results = groupResults;
     list_pairs = listPairs;
     display_pairs = displayPairs;
     imprint_results = imprintResults;
  }

  int prog_step = 10;

  if( show_messages )
  {
    PRINT_INFO( "Finding surface overlap...\n" );
    if( ref_face_list.size() > prog_step )
    {
      char message[128];
      sprintf(message, "Finding Surface Overlap On %d Surfaces", ref_face_list.size() );
      AppUtil::instance()->progress_tool()->start(0, ref_face_list.size(),
        "Progress", message, TRUE, TRUE );
    }
  }
  else
    prog_step = -1;

  status = find_overlapping_surfaces( ref_face_list, ref_face_list1,
    ref_face_list2, faces_to_draw, list_pairs, prog_step, filter_slivers );

  if( show_messages && ref_face_list.size() > prog_step )
    AppUtil::instance()->progress_tool()->end();

  if( faces_to_draw.size() )
  {
    if( group_results == CUBIT_TRUE )
    {
      RefGroup *new_group = RefEntityFactory::instance()->construct_RefGroup( "surf_overlap" );
      new_group->add_ref_entity( faces_to_draw );
      CubitString name = new_group->entity_name();
      PRINT_INFO( "Found %d overlapping surface pairs (added to group '%s')\n", 
        ref_face_list1.size(), name.c_str() );
    }
    else if( show_messages )
      PRINT_INFO( "Found %d overlapping surface pairs\n", ref_face_list1.size() );

    if ( imprint_results )
    {
      CubitStatus stat = imprint(ref_face_list1, ref_face_list2);
      if ( stat != CUBIT_SUCCESS )
      {
        PRINT_WARNING("Imprinting overlaps was unsuccessful\n");
      }
    }
  }
  else if (show_messages )
    PRINT_INFO( "Found 0 overlapping surface pairs\n" );

  return status;
}

CubitStatus
SurfaceOverlapTool::find_overlapping_surfaces( DLIList<RefFace*> &ref_face_list,
                                              DLIList<RefFace*> &ref_face_list1,
                                              DLIList<RefFace*> &ref_face_list2,
                                              DLIList<RefEntity*> &pair_list,
                                              CubitBoolean list_pairs,
                                              int prog_step,
                                              bool filter_slivers)
{
  int number_pairs = 0;
  CubitBoolean abort = CUBIT_FALSE;
  RefEntity* ref_entity;

  // Check each surface with each one later in the list
  RefFace *ref_face_ptr1, *ref_face_ptr2;

  // Populate the RefFace AbstractTree
  AbstractTree<RefFace*> *a_tree = new RTree<RefFace*>( gapMax );
  int i;
  ref_face_list.reset();
  for( i=ref_face_list.size(); i--; )
  {
    RefFace *ref_face_ptr = ref_face_list.get_and_step();
    a_tree->add( ref_face_ptr );
  }

  // Main loop for finding overlapping surfaces
  ref_face_list.reset();
  for( i=ref_face_list.size(); i--; )
  {
    // Cancel button pushed or cntrl-C
    if (CubitMessage::instance()->Interrupt()) 
    {
      PRINT_INFO("Find overlap operation aborted.\n");
      goto done;
    }

    ref_face_ptr1 = ref_face_list.get_and_step();

    // Remove this surface from AbstractTree so it is not found and never
    // found again
    a_tree->remove( ref_face_ptr1 );

    // Find RefFaces from AbstractTree that are within range of this surface
    CubitBox ref_face1_box = ref_face_ptr1->bounding_box();
    DLIList<RefFace*> close_ref_faces;
    a_tree->find( ref_face1_box, close_ref_faces );

    int j;
    for( j=close_ref_faces.size(); j--; )
    {
      // Cancel button pushed or cntrl-C
      if (CubitMessage::instance()->Interrupt()) 
      {
         PRINT_INFO("Find overlap operation aborted.\n");
         goto done;
      }

      ref_face_ptr2 = close_ref_faces.get_and_step();

      bool overlap = check_overlap( ref_face_ptr1, ref_face_ptr2, abort );
      if(overlap == CUBIT_TRUE && filter_slivers)
      {
        RefFace *f1, *f2;
        DLIList<RefEdge*> f1_edges, f2_edges, *f1_edge_list, *f2_edge_list;
        ref_face_ptr1->ref_edges(f1_edges);
        ref_face_ptr2->ref_edges(f2_edges);
        if(f1_edges.size() > f2_edges.size())
        {
          f1 = ref_face_ptr2;
          f2 = ref_face_ptr1;
          f1_edge_list = &f2_edges;
          f2_edge_list = &f1_edges;
        }
        else
        {
          f1 = ref_face_ptr1;
          f2 = ref_face_ptr2;
          f1_edge_list = &f1_edges;
          f2_edge_list = &f2_edges;
        }
        int b;
        for(b=f1_edge_list->size(); b>0 && overlap; b--)
        {
          RefEdge *cur_edge = f1_edge_list->get_and_step();
          if(cur_edge->is_merged())
          {
            DLIList<RefFace*> face_list;
            cur_edge->ref_faces(face_list);
            if(face_list.is_in_list(f2))
            {
              CubitVector mid_pt;
              cur_edge->mid_point(mid_pt);
              CubitVector f1_norm = f1->normal_at(mid_pt); 
              CubitVector f2_norm = f2->normal_at(mid_pt); 
              DLIList<CoEdge*> f1_coedges, f2_coedges;
              cur_edge->get_co_edges(f1_coedges, f1);
              cur_edge->get_co_edges(f2_coedges, f2);
              if(f1_coedges.size() == 1 && f2_coedges.size() == 1)
              {
                CoEdge *ce1 = f1_coedges.get();
                CoEdge *ce2 = f2_coedges.get();
                CubitVector curve_dir;

                if( cur_edge->get_curve_ptr()->geometry_type() == STRAIGHT_CURVE_TYPE )
                  cur_edge->get_point_direction(mid_pt, curve_dir);
                else
                  cur_edge->tangent( mid_pt, curve_dir );

                CubitVector ce1_dir, ce2_dir;
                if(ce1->get_sense() == CUBIT_REVERSED)
                  ce1_dir = -curve_dir;
                else
                  ce1_dir = curve_dir;
                if(ce2->get_sense() == CUBIT_REVERSED)
                  ce2_dir = -curve_dir;
                else
                  ce2_dir = curve_dir;
                CubitVector in_dir1 = f1_norm * ce1_dir;
                CubitVector in_dir2 = f2_norm * ce2_dir;
                if(in_dir1 % in_dir2 < 0.0)
                  overlap = false;
              }
            }
          }
        }
      }
      if(overlap)
      {
        if( abort == CUBIT_TRUE )
          goto done;

        if( list_pairs == CUBIT_TRUE )
          PRINT_INFO( " Surface %d and %d overlap\n", ref_face_ptr1->id(),
          ref_face_ptr2->id() );
        
        number_pairs++;
        
        ref_entity = CAST_TO(ref_face_ptr1,RefEntity);
        pair_list.append_unique( ref_entity );
        
        ref_entity = CAST_TO(ref_face_ptr2,RefEntity);
        pair_list.append_unique( ref_entity );
        
        ref_face_list1.append( ref_face_ptr1 );
        ref_face_list2.append( ref_face_ptr2 );
      }

      if( abort == CUBIT_TRUE )
        goto done;
    }

    // Free memory, since this surface will never be accessed again.  This
    // helps to reduce memory required.
    ref_face_ptr1->delete_TD( &TDSurfaceOverlap::is_surface_overlap );

    if( prog_step>0 && ref_face_list.size() > prog_step )
      AppUtil::instance()->progress_tool()->step();
  }

done:  

  // Make sure all tool datas are deleted
  for( i=ref_face_list.size(); i--; )
    ref_face_list.get_and_step()->delete_TD( &TDSurfaceOverlap::is_surface_overlap );

  delete a_tree;

  return CUBIT_SUCCESS;
}

//Currently this function is only called when using tolerant imprinting.
//It does not use the settings controlled by the user for this tool

CubitBoolean
SurfaceOverlapTool::check_overlap( Surface *surface1, Surface *surface2,
              std::map<Surface*, DLIList<SurfaceOverlapFacet*>* > *facet_map,
              std::map<Surface*, double > *area_map,
              std::map<Surface*, AbstractTree<SurfaceOverlapFacet*>* > *a_tree_map )
{
  if( surface1 == surface2 )
    return CUBIT_FALSE;

  //if surfaces are not splines and are not of the same type, 
  //they won't overlap
  if( (surface1->geometry_type() != SPLINE_SURFACE_TYPE  &&
       surface2->geometry_type() != SPLINE_SURFACE_TYPE) &&
      (surface1->geometry_type() != surface2->geometry_type() )) 
    return CUBIT_FALSE;

  int i, j;
  AnalyticGeometryTool::instance();
  double opp_low = 180.0 - angleMax;
  double opp_high = 180.0 - angleMin;

  std::map<Surface*, DLIList<SurfaceOverlapFacet*>* >::iterator facet_iterator;

  DLIList<SurfaceOverlapFacet*> *facet_list1 = NULL;
  DLIList<SurfaceOverlapFacet*> *facet_list2 = NULL;

  //see if surface is in map...if not we have to create faceting for it.
  if( facet_map )
    facet_iterator = facet_map->find( surface1 );

  if( facet_map == NULL || facet_iterator == facet_map->end() )
  {
    //for non-planar surfaces, facet wrt the smallest curve of the surface 
    double min_edge_length = 0.0;

    if( surface1->geometry_type() != PLANE_SURFACE_TYPE )
    {
      DLIList<Curve*> tmp_curves;
      surface1->curves( tmp_curves );
      CubitBox surface_box = surface1->bounding_box();

      //ignore curves that are really small
      double min_length_threshold = (surface_box.diagonal().length())*0.01;
      if( tmp_curves.size() )
      {
        min_edge_length = CUBIT_DBL_MAX; 
        double tmp_length; 
        for( i=tmp_curves.size(); i--; )
        {
          tmp_length = tmp_curves.get_and_step()->measure();
          if( tmp_length > min_length_threshold && tmp_length < min_edge_length )
            min_edge_length = tmp_length;
        }
      }
      if( min_edge_length == CUBIT_DBL_MAX )
        min_edge_length = 0.0;
      else
        min_edge_length *= 2;
    }

    facet_list1 = new DLIList<SurfaceOverlapFacet*>;
    GMem surface_graphics;
    surface1->get_geometry_query_engine()->get_graphics( surface1, &surface_graphics, 
                                facetAngTol, facetAbsTol, min_edge_length );   

    GPoint* plist = surface_graphics.point_list();
    int* facet_list = surface_graphics.facet_list();

    GPoint p[3];
    for (i = 0; i < surface_graphics.fListCount; )
    {
      int sides = facet_list[i++];
      if (sides != 3)
      {
        PRINT_WARNING("Skipping n-sided polygone in triangle list"
                      " in TDSurfaceOverlap.\n");
        i += sides;
      }
      else
      {
        p[0] = plist[facet_list[i++]];
        p[1] = plist[facet_list[i++]];
        p[2] = plist[facet_list[i++]];
     
        SurfaceOverlapFacet *facet = new SurfaceOverlapFacet( p );
        facet_list1->append( facet );
      }
    }

    if( facet_map )
      facet_map->insert( std::map<Surface*, 
        DLIList<SurfaceOverlapFacet*>*>::value_type( surface1, facet_list1 )); 
  }
  else
    facet_list1 = facet_iterator->second;

  //see if surface is in map...if not we have to create faceting for it.
  if( facet_map )
    facet_iterator = facet_map->find( surface2 );

  if( facet_map == NULL || facet_iterator == facet_map->end() )
  {
    //for non-planar surfaces, facet wrt the smallest curve of the surface 
    double min_edge_length = 0.0;

    if( surface2->geometry_type() != PLANE_SURFACE_TYPE )
    {
      DLIList<Curve*> tmp_curves;
      surface2->curves( tmp_curves );
      CubitBox surface_box = surface2->bounding_box();

      //ignore curves that are really small
      double min_length_threshold = (surface_box.diagonal().length())*0.01;
      if( tmp_curves.size() )
      {
        min_edge_length = CUBIT_DBL_MAX; 
        double tmp_length; 
        for( i=tmp_curves.size(); i--; )
        {
          tmp_length = tmp_curves.get_and_step()->measure();
          if( tmp_length > min_length_threshold && tmp_length < min_edge_length )
            min_edge_length = tmp_length;
        }
      }
      if( min_edge_length == CUBIT_DBL_MAX )
        min_edge_length = 0.0;
      else
        min_edge_length *= 2;
    }

    facet_list2 = new DLIList<SurfaceOverlapFacet*>;
    GMem surface_graphics;
    surface2->get_geometry_query_engine()->get_graphics( surface2, &surface_graphics, 
                                facetAngTol, facetAbsTol, min_edge_length );

    GPoint* plist = surface_graphics.point_list();
    int* facet_list = surface_graphics.facet_list();

    GPoint p[3];
    for (i = 0; i < surface_graphics.fListCount; )
    {
      int sides = facet_list[i++];
      if (sides != 3)
      {
        PRINT_WARNING("Skipping n-sided polygone in triangle list"
                      " in TDSurfaceOverlap.\n");
        i += sides;
      }
      else
      {
        p[0] = plist[facet_list[i++]];
        p[1] = plist[facet_list[i++]];
        p[2] = plist[facet_list[i++]];
     
        SurfaceOverlapFacet *facet = new SurfaceOverlapFacet( p );
        facet_list2->append( facet );
      }
    }

    if( facet_map )
      facet_map->insert( std::map<Surface*, 
        DLIList<SurfaceOverlapFacet*>*>::value_type( surface2, facet_list2 )); 
  }           
  else
    facet_list2 = facet_iterator->second;

  // Compare facets
  int num_tri1 = facet_list1->size();
  if( !num_tri1 )
  {
    PRINT_WARNING( "Unable to facet surface\n" );
    return CUBIT_FALSE;
  }
  
  int num_tri2 = facet_list2->size();
  if( !num_tri2 )
  {
    PRINT_WARNING( "Unable to facet surface\n" );
    return CUBIT_FALSE;
  }
  
  Surface *tmp_surf1 = surface1;
  Surface *tmp_surf2 = surface2;

  // Compare least to most - possibly switch the lists
  if( facet_list1->size() > facet_list2->size() )
  {
    DLIList<SurfaceOverlapFacet*> *temp_list = facet_list1;
    facet_list1 = facet_list2;
    facet_list2 = temp_list;
    tmp_surf1 = surface2;
    tmp_surf2 = surface1;
  }

  // Possibly use an AbstractTree for facet_list2
  AbstractTree<SurfaceOverlapFacet*> *a_tree = NULL;
  if( facet_list2->size() > NO_FACETS_FOR_ABSTRACTTREE )
  { 
    std::map<Surface*, AbstractTree<SurfaceOverlapFacet*>* >::iterator iter1, iter2; 
    // If the same size, use the existing AbstractTree if one has
    // one and the other doesn't.  This probably won't 
    // save time, but there is a miniscule chance it might.
    if( facet_list1->size() == facet_list2->size() )
    {
      //see if both are in the a_tree_map
      iter1 = a_tree_map->find( tmp_surf1 );  
      iter2 = a_tree_map->find( tmp_surf2 );  

      if( !(iter2 != a_tree_map->end()) && (iter1 != a_tree_map->end()) ) 
      {
        // Switch them, just to use the existing AbstractTree
        DLIList<SurfaceOverlapFacet*> *temp_list = facet_list1;
        facet_list1 = facet_list2;
        facet_list2 = temp_list;
        Surface *tmp_surf = tmp_surf1;
        tmp_surf1 = tmp_surf2;
        tmp_surf2 = tmp_surf;
        std::map<Surface*, AbstractTree<SurfaceOverlapFacet*>* >::iterator tmp_iter;
        tmp_iter = iter1;
        iter1 = iter2;
        iter2 = tmp_iter; 
      }
    }
  
    //populate tree if necessary
    iter2 = a_tree_map->find( tmp_surf2 ); 
    if( iter2 == a_tree_map->end() ) 
    {
      a_tree = new RTree<SurfaceOverlapFacet*>( gapMax );

      for( i=facet_list2->size(); i--; )
        a_tree->add( facet_list2->get_and_step() ); 
      
      a_tree_map->insert( std::map<Surface*, AbstractTree<SurfaceOverlapFacet*>*>::value_type(
                                tmp_surf2, a_tree ));
    }
    else
      a_tree = iter2->second;
  }

  double area = 0;

  //set ovelap tolerance to one thousandth of the smaller area.
  double face1_area, face2_area;
  if( area_map )
  {
    std::map<Surface*, double>::iterator area_iter;
    area_iter = area_map->find( tmp_surf1);
    if( area_iter == area_map->end() )
    {
      face1_area = tmp_surf1->measure();
      area_map->insert( std::map<Surface*, double>::value_type( tmp_surf1, face1_area ) ); 
    }
    else
      face1_area = area_iter->second; 

    area_iter = area_map->find( tmp_surf2);
    if( area_iter == area_map->end() )
    {
      face2_area = tmp_surf2->measure();
      area_map->insert( std::map<Surface*, double>::value_type( tmp_surf2, face2_area )); 
    }
    else
      face2_area = area_iter->second; 
  }
  else
  {
    face1_area = tmp_surf1->measure();
    face2_area = tmp_surf2->measure();
  }

  bool sample_deviation = false; 
  Surface *larger_surface = NULL;
  //if the surfaces are non planar and one is much bigger than the other one...
  //adjust the tolerance some
  if( (tmp_surf1->geometry_type() != PLANE_SURFACE_TYPE  &&
       tmp_surf2->geometry_type() != PLANE_SURFACE_TYPE) )
  {
    double length1 = tmp_surf1->bounding_box().diagonal().length();
    double length2 = tmp_surf2->bounding_box().diagonal().length();
    
    double ratio = 0;
    if(length1 > length2)
    {
      ratio = length1/length2;
      larger_surface = tmp_surf1;
    }
    else
    {
      ratio = length2/length1;
      larger_surface = tmp_surf2;
    }
   
    if( ratio > 50 )
      sample_deviation = true;
  }

  double tmp_overlapTolerance; 
  if( face1_area < face2_area )
    tmp_overlapTolerance = face1_area * 0.001;
  else
    tmp_overlapTolerance = face2_area * 0.001;


  //If you're comparing overlap between a large surface and a very small one, 
  //you might have to adjust tolerances because of the faceting.  Here we try 
  //to determine the deviation of the faceting of the larger surface.  I 
  //take the mid point of the smallest edge on the facet, thinking the smallest 
  //edge is approximating high curvature so it would deviate the most and thus
  //provide the safest tolerance.
  double facet_compare_tol = 0;
  if( sample_deviation )
  {
    DLIList<SurfaceOverlapFacet*> *tmp_facet_list = NULL;
    if( larger_surface == tmp_surf1 )
      tmp_facet_list = facet_list1;
    else
      tmp_facet_list = facet_list2;

    int num_samples = 0;
    for( i=tmp_facet_list->size(); i--; )
    {
      SurfaceOverlapFacet *tmp_facet = tmp_facet_list->get_and_step();
        //we used to be lucky, where one seventh of the facets would give us a 
        //decent tolerance...this doesn't work now....not robust enough
//      if( (i%7) == 0)  // commenting this out....sometimes you haveto 
//      {
        CubitVector point = tmp_facet->smallest_edge_midpoint(); 
        //determine distance between surface and centroid 

        CubitVector tmp_point;
        larger_surface->closest_point_trimmed( point, tmp_point );  
        double tmp_distance = point.distance_between( tmp_point ); 
        if( tmp_distance > facet_compare_tol )
        {
          facet_compare_tol = tmp_distance;
        }
//      }
    }
  }
  else
  {
    for( i=facet_list1->size(); i--; )
      facet_compare_tol += facet_list1->get_and_step()->bounding_box().diagonal().length();

    for( i=facet_list2->size(); i--; )
      facet_compare_tol += facet_list2->get_and_step()->bounding_box().diagonal().length();

    if( facet_list1->size() || facet_list2->size() )
    {
      facet_compare_tol /= (facet_list1->size() + facet_list2->size() );
      facet_compare_tol *= 0.0025;
    }
  }
  
  //if the merge tolerance is a larger, use it
  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;
  if( tolerance > facet_compare_tol )
    facet_compare_tol = tolerance;

  facet_list1->reset();
  for( i=facet_list1->size(); i--; )
  {
    SurfaceOverlapFacet *facet1 = facet_list1->get_and_step();
    
    CubitBox facet1_bbox = facet1->bounding_box();
    DLIList<SurfaceOverlapFacet*> close_facets;
    if( a_tree )
    {
      a_tree->find( facet1_bbox, close_facets );
      facet_list2 = &close_facets;
    }
    
    facet_list2->reset();
    for( j=facet_list2->size(); j--; )
    {
      SurfaceOverlapFacet *facet2 = facet_list2->get_and_step();
      
      // Check angle between triangles, must be within criteria
      double ang = facet1->angle( *facet2 ); 
      
      // Allow overlap for angles close to 180 and 0 degrees
      
      // normalType - 1=any, 2=opposite, 3=same
      //ang>=180.0-angt || ang<angt
      if( ((normalType==1 && ang>=opp_low && ang<=opp_high) ||
          (normalType==1 && ang>=angleMin && ang<=angleMax) ||
          (normalType==2 && ang>=opp_low && ang<=opp_high) ||
          (normalType==3 && ang>=angleMin && ang<=angleMax)) && 
          (normalType != 2 || !skipFacingSurfaces || !facet1->facing( *facet2 )))// check to make sure the surfaces are not facing 
      {
          // If triangle bounding boxes don't intersect - no overlap
          if( facet1->bbox_overlap( facet_compare_tol, *facet2 ) )
          {
              // Check distance between triangles, must be within criteria
              double dist = facet1->distance( *facet2 );
              if( dist >= gapMin && dist <= facet_compare_tol )
              {
                  // Check for projected overlap
                  // We want sum of area of ALL overlapping facets
                  area += facet1->projected_overlap( *facet2 );
                  if( area > tmp_overlapTolerance )
                  {
                      return CUBIT_TRUE;
                  }                
              }
          }
      }
    }
  }

  return CUBIT_FALSE;
}

CubitBoolean 
SurfaceOverlapTool::check_overlap( DLIList<SurfaceOverlapFacet*> *facet_list1,
                                   DLIList<SurfaceOverlapFacet*> *facet_list2,
                                   AbstractTree<SurfaceOverlapFacet*> *a_tree, 
                                   CubitBoolean abort, 
                                   CubitBoolean draw_overlap,
                                   double *overlap_area ) 
{
  double area = 0;
  double opp_low = 180.0 - angleMax;
  double opp_high = 180.0 - angleMin;

  facet_list1->reset();
  int num_tri1 = facet_list1->size();
  int i,j;
  for( i=num_tri1; i--; )
  {
    // Cancel button pushed or cntrl-C
    if (CubitMessage::instance()->Interrupt()) 
    {
      PRINT_INFO("Find overlap operation aborted.\n");
      abort = CUBIT_TRUE;
      return CUBIT_FALSE;
    }

    SurfaceOverlapFacet *facet1 = facet_list1->get_and_step();

    DLIList<SurfaceOverlapFacet*> close_facets;
    if( a_tree )
    {
      a_tree->find( facet1->bounding_box(), close_facets );
      facet_list2 = &close_facets;
    }
    
    facet_list2->reset();
    for( j=facet_list2->size(); j--; )
    {
      // Cancel button pushed or cntrl-C
      if (CubitMessage::instance()->Interrupt()) 
      {
        PRINT_INFO("Find overlap operation aborted.\n");
        abort = CUBIT_TRUE;
        return CUBIT_FALSE;
      }

      SurfaceOverlapFacet *facet2 = facet_list2->get_and_step();
      
      // Check angle between triangles, must be within criteria
      double ang = facet1->angle( *facet2 ); 

      // Allow overlap for angles close to 180 and 0 degrees
      
      // normalType - 1=any, 2=opposite, 3=same
      //ang>=180.0-angt || ang<angt
      if( ((normalType==1 && ang>=opp_low && ang<=opp_high) ||
        (normalType==1 && ang>=angleMin && ang<=angleMax) ||
        (normalType==2 && ang>=opp_low && ang<=opp_high) ||
        (normalType==3 && ang>=angleMin && ang<=angleMax) ) &&
        (normalType != 2 || !skipFacingSurfaces || !facet1->facing( *facet2 )))// check to make sure the surfaces are not facing )
      {
        // If triangle bounding boxes don't intersect - no overlap
        if( facet1->bbox_overlap( gapMax, *facet2 ) )
        {
          // Check distance between triangles, must be within criteria
          double dist = facet1->distance( *facet2 );
          if( dist >= gapMin && dist <= gapMax )
          {
            // Check for projected overlap
            // We want sum of area of ALL overlapping facets
            area += facet1->projected_overlap( *facet2, draw_overlap );
            if( area > overlapTolerance &&
               (draw_overlap == CUBIT_FALSE && overlap_area == NULL ) )
            {
              return CUBIT_TRUE;
            }                
          }
        }
      }
    }
  }  
  if( draw_overlap == CUBIT_TRUE ) 
  {
    PRINT_INFO("Total overlapping area = %f\n", area );

    if( area > 0.0 )
      return CUBIT_TRUE;
  }

  if( overlap_area )
  {
    if( area > 0.0 )
    {
      *overlap_area = area;
      return CUBIT_TRUE;
    }
  }

  return CUBIT_FALSE;
}

CubitBoolean
SurfaceOverlapTool::check_overlap( RefFace *ref_face_ptr1, RefFace *ref_face_ptr2, 
                                   CubitBoolean abort, 
                                   CubitBoolean draw_overlap,
                                   double *overlap_area ) 
{
  if( ref_face_ptr1 == ref_face_ptr2 )
    return CUBIT_FALSE;

  //if surfaces are not splines and are not of the same type, 
  //they won't overlap
  if( (ref_face_ptr1->get_surface_ptr()->geometry_type() != SPLINE_SURFACE_TYPE  &&
       ref_face_ptr2->get_surface_ptr()->geometry_type() != SPLINE_SURFACE_TYPE) &&
      (ref_face_ptr1->get_surface_ptr()->geometry_type() != 
       ref_face_ptr2->get_surface_ptr()->geometry_type() )) 
    return CUBIT_FALSE;
 
  AnalyticGeometryTool::instance();
  abort = CUBIT_FALSE;

  // Check for overlap between the found surfaces using the facets
  TDSurfaceOverlap *tdso_1;
  tdso_1 = (TDSurfaceOverlap *)ref_face_ptr1->get_TD(&TDSurfaceOverlap::is_surface_overlap);
  if( !tdso_1 )
  {
    ref_face_ptr1->add_TD( new TDSurfaceOverlap( ref_face_ptr1, facetAngTol, facetAbsTol,
                                                 gapMax ) );
    tdso_1 = (TDSurfaceOverlap *)ref_face_ptr1->get_TD(&TDSurfaceOverlap::is_surface_overlap);
  }

  TDSurfaceOverlap *tdso_2;
  tdso_2 = (TDSurfaceOverlap *)ref_face_ptr2->get_TD(&TDSurfaceOverlap::is_surface_overlap);
  if( !tdso_2 )
  {
    ref_face_ptr2->add_TD( new TDSurfaceOverlap( ref_face_ptr2, facetAngTol, facetAbsTol,
                                                 gapMax ) );
    tdso_2 = (TDSurfaceOverlap *)ref_face_ptr2->get_TD(&TDSurfaceOverlap::is_surface_overlap);
  }
  
  // Check if within the same body - maybe we don't need to consider this
  if( checkWithinBodies == CUBIT_FALSE )
  {
    DLIList<Body*> *body_list_ptr1, *body_list_ptr2;
    body_list_ptr1 = tdso_1->get_body_list();
    body_list_ptr2 = tdso_2->get_body_list();
    DLIList<Body*> shared_body_list = *body_list_ptr1;
    shared_body_list.intersect( *body_list_ptr2 );
    if( shared_body_list.size() )
      return CUBIT_FALSE;
  }

  // Check if in different bodies - maybe we don't need to consider this
  if( checkAcrossBodies == CUBIT_FALSE )
  {
    DLIList<Body*> *body_list_ptr1, *body_list_ptr2;
    body_list_ptr1 = tdso_1->get_body_list();
    body_list_ptr2 = tdso_2->get_body_list();
    DLIList<Body*> shared_body_list = *body_list_ptr1;
    shared_body_list.intersect( *body_list_ptr2 );
    if( shared_body_list.size() == 0 )
      return CUBIT_FALSE;
  }
  
  // Compare facets
  DLIList<SurfaceOverlapFacet*> *facet_list1 = tdso_1->get_facet_list();
  int num_tri1 = facet_list1->size();
  if( !num_tri1 )
  {
    PRINT_WARNING( "Unable to facet surface %d\n", ref_face_ptr1->id() );
    return CUBIT_FALSE;
  }
  
  DLIList<SurfaceOverlapFacet*> *facet_list2 = tdso_2->get_facet_list();
  int num_tri2 = facet_list2->size();
  if( !num_tri2 )
  {
    PRINT_WARNING( "Unable to facet surface %d\n", ref_face_ptr2->id() );
    return CUBIT_FALSE;
  }
  
  PRINT_DEBUG_102( " Comparing %d facets in Surface %d with %d facets in Surface %d...\n",
    num_tri1, ref_face_ptr1->id(), num_tri2, ref_face_ptr2->id() );

  // Compare least to most - possibly switch the lists
  if( facet_list1->size() > facet_list2->size() )
  {
    DLIList<SurfaceOverlapFacet*> *temp_list = facet_list1;
    facet_list1 = facet_list2;
    facet_list2 = temp_list;
    TDSurfaceOverlap *temp_tdso = tdso_1;
    tdso_2 = tdso_1;
    tdso_1 = temp_tdso;
  }

  // Possibly use an AbstractTree for facet_list2
  AbstractTree<SurfaceOverlapFacet*> *a_tree = NULL;
  if( facet_list2->size() > NO_FACETS_FOR_ABSTRACTTREE )
  { 
    // If the same size, use the existing AbstractTree if one has
    // one and the other doesn't.  This probably won't 
    // save time, but there is a miniscule chance it might.
    if( facet_list1->size() == facet_list2->size() )
    {
      if( !tdso_2->has_rtree() && tdso_1->has_rtree() )
      {
        // Switch them, just to use the existing AbstractTree
        DLIList<SurfaceOverlapFacet*> *temp_list = facet_list1;
        facet_list1 = facet_list2;
        facet_list2 = temp_list;
        TDSurfaceOverlap *temp_tdso = tdso_1;
        tdso_2 = tdso_1;
        tdso_1 = temp_tdso;
      }
    }
    a_tree = tdso_2->get_facet_rtree();
  }

  bool surfs_overlap = check_overlap( facet_list1, facet_list2, a_tree,
                                      abort, draw_overlap, overlap_area );

  return surfs_overlap;
}

void SurfaceOverlapTool::list_settings()
{
  char on_off[2][4];
  strcpy( on_off[0], "Off" );
  strcpy( on_off[1], "On" );

  PRINT_INFO( "Surface Overlap Algorithm Settings:\n" );
  if( facetAbsTol == 0.0 )
    PRINT_INFO( " Facetting Absolute Tolerance: 0.0 (using default solid modeler setting)\n" );
  else
    PRINT_INFO( " Facetting Absolute Tolerance: %f\n", facetAbsTol );

  if( facetAngTol == 0 )
    PRINT_INFO( " Facetting Angle Tolerance: 0 (using default solid modeler setting)\n" );
  else
    PRINT_INFO( " Facetting Angle Tolerance: %d\n", facetAngTol );

  PRINT_INFO( " Gap Range: %f to %f\n", gapMin, gapMax );
  PRINT_INFO( " Angle Range: %.1f to %.1f degrees\n", angleMin, angleMax );
  PRINT_INFO( " Overlap Area Tolerance: %f\n", overlapTolerance );

  switch( normalType )
  {
  case 1:
    PRINT_INFO( " Pair Normals Allowed: Any\n" );
    break;
  case 2:
    PRINT_INFO( " Pair Normals Allowed: Opposite\n" );
    break;
  case 3:
    PRINT_INFO( " Pair Normals Allowed: Same\n" );
    break;
  }

  PRINT_INFO( " Display Pairs: %s\n", on_off[displayPairs] );
  PRINT_INFO( " List Pairs: %s\n", on_off[listPairs] );
  PRINT_INFO( " Group Results: %s\n", on_off[groupResults] );
  PRINT_INFO( " Across Body|Volume: %s\n", on_off[checkAcrossBodies] );
  PRINT_INFO( " Within Body|Volume: %s\n", on_off[checkWithinBodies] );
  PRINT_INFO( " Imprint: %s\n", on_off[imprintResults] );
}

CubitString SurfaceOverlapTool::get_normal_type_setting()
{
  if (normalType == 2) {
    return CubitString("opposite");
  }
  else if (normalType == 3) {
    return CubitString("same");
  }
  else {
    return CubitString("any");
  }
}

void SurfaceOverlapTool::set_normal_type_setting(CubitString type )
{
  if (CubitUtil::compare(type.c_str(), "opposite")) {
    normalType = 2;
  }
  else if (CubitUtil::compare(type.c_str(), "same")) {
    normalType = 3;
  }
  else {
    normalType = 1;
  }
}

int SurfaceOverlapTool::get_group_results_setting()
{return groupResults;}

void SurfaceOverlapTool::set_group_results_setting( int setting )
{groupResults = (setting) ? CUBIT_TRUE : CUBIT_FALSE;}

int SurfaceOverlapTool::get_list_pairs_setting()
{return listPairs;}

void SurfaceOverlapTool::set_list_pairs_setting( int setting )
{listPairs = (setting) ? CUBIT_TRUE : CUBIT_FALSE;}

int SurfaceOverlapTool::get_display_pairs_setting()
{return displayPairs;}

void SurfaceOverlapTool::set_display_pairs_setting( int setting )
{displayPairs = (setting) ? CUBIT_TRUE : CUBIT_FALSE;}

int SurfaceOverlapTool::get_facet_ang_tol_setting()
{return facetAngTol;}

void SurfaceOverlapTool::set_facet_ang_tol_setting( int val )
{facetAngTol=(unsigned short)val;}

//Imprints the bodies containing the lists of surfaces.  Also
//uses the edges imprinting rather than the normal imprint.
CubitStatus SurfaceOverlapTool::imprint(DLIList<RefFace*> &ref_face_list1,
                                        DLIList<RefFace*> &ref_face_list2)
{
  //This is a very tricky goal.  The goal is to imprint each pair for the
  //two lists.  The problem is that we need to imprint the bodies of the
  //surfaces not just the surfaces since this isn't as robust.  The problem
  //is that more than one of these surfaces may reference the same body.  So
  //in fact as we imprint one of the surfaces, the other ones could get deleted and
  //create a memory leak.
  
  //The only thing I can think of is to use the id's of the bodies as the constant.
  //On the imprint, the body id should not change.  So lets create an array of bodies,
  //and update the array as they get changed.  The array will allow constant access time.
  DLIList <Body*> body_list, temp_bodies;
  DLIList <int> body_ids1, body_ids2, modified_list;
  DLIList <DLIList<RefEdge*>*> edges_list1, edges_list2;
  DLIList <RefEdge*> *ref_edges1, *ref_edges2, tmp_edges;
  
  CubitBoolean print = CUBIT_TRUE;
  Body *tmp_body;
  int i, max_id = -1;
  
    //first get the max id of the bodies here to help us
    //get the size of the array.
  for ( i = ref_face_list1.size(); i > 0; i-- )
  {
    temp_bodies.clean_out();
    RefFace *ref_face1 = ref_face_list1.get_and_step();
    if ( ref_face1 == NULL )
    {
      PRINT_ERROR("Bad data sent to imprint overlaps.\n");
      return CUBIT_FAILURE;
    }
    ref_face1->bodies(temp_bodies);
    if ( temp_bodies.size() != 1 )
    {
      PRINT_WARNING("Entities must be in one body\n");
      PRINT_WARNING("Skipping that entity\n");
      ref_face_list2.step();
      continue;
    }
    tmp_body = temp_bodies.get();
    temp_bodies.clean_out();
    RefFace *ref_face2 = ref_face_list2.get_and_step();
    if ( ref_face2 == NULL )
    {
      PRINT_ERROR("Bad data sent to imprint overlaps.\n");
      return CUBIT_FAILURE;
    }
    ref_face2->bodies(temp_bodies);
    if ( temp_bodies.size() != 1 )
    {
      PRINT_WARNING("Entities must be in one body\n");
      PRINT_WARNING("Skipping that entity\n");
      continue;
    }
    if ( tmp_body->id() > max_id )
      max_id = tmp_body->id();
    body_list.append(tmp_body);
    body_ids1.append(tmp_body->id());
    body_list.append(temp_bodies.get());
    body_ids2.append(temp_bodies.get()->id());
    if ( temp_bodies.get()->id() > max_id )
      max_id = temp_bodies.get()->id();
      //Now set up the edge lists.
    ref_edges1 = new DLIList<RefEdge*>;
    ref_face1->ref_edges(tmp_edges);
    CubitStatus stat = copy_edges_in_list(tmp_edges, *ref_edges1);
    if ( stat != CUBIT_SUCCESS )
      return CUBIT_FAILURE;
    ref_edges2 = new DLIList<RefEdge*>;
    tmp_edges.clean_out();
    ref_face2->ref_edges(tmp_edges);
    stat = copy_edges_in_list(tmp_edges, *ref_edges2);
    if ( stat != CUBIT_SUCCESS )
      return CUBIT_FAILURE;
    edges_list1.append(ref_edges1);
    edges_list2.append(ref_edges2);
  }
    //Now create an array and store all of the bodies.
  Body **body_array = new Body* [2*max_id+1];
    //initialize the array to null.
  for ( i = 0; i < 2*max_id+1; i++ )
  {
    body_array[i] = (Body*)NULL;
  }
    //Now put the bodies into the array.
  for ( i = body_list.size(); i > 0; i-- )
  {
    tmp_body = body_list.get_and_step();
    body_array[tmp_body->id()] = tmp_body;
  }
    //Now go through and pairwise imprint the bodies.
  DLIList <Body*> bodies, tmp_modified_bodies;
  for( i = body_ids1.size(); i > 0; i-- )
  {
    tmp_modified_bodies.clean_out();
    int id1 = body_ids1.get_and_step();
    int id2 = body_ids2.get_and_step();
    ref_edges1 = edges_list1.get_and_step();
    ref_edges2 = edges_list2.get_and_step();
    
    Body *body1 = body_array[id1];
    Body *body2 = body_array[id2];
    int body1_num = num_descendants(body1);
    int body2_num = num_descendants(body2);
    bodies.clean_out();
    bodies.append(body2);
    CubitStatus stat = GeometryModifyTool::instance()->
      imprint(bodies, *ref_edges1, tmp_modified_bodies, CUBIT_FALSE, CUBIT_FALSE);
    if ( stat == CUBIT_SUCCESS && tmp_modified_bodies.size() == 1 )
    {
      body2 = tmp_modified_bodies.get();
    }
    tmp_modified_bodies.clean_out();
    bodies.clean_out();
    bodies.append(body1);
    stat = GeometryModifyTool::instance()->
      imprint(bodies, *ref_edges2, tmp_modified_bodies, CUBIT_FALSE, CUBIT_FALSE);
    if ( stat == CUBIT_SUCCESS && tmp_modified_bodies.size() == 1 )
    {
      body1 = tmp_modified_bodies.get();
    }
    body_array[id1] = body1;
    int temp_num = num_descendants(body1);
    if ( temp_num != body1_num )
    {
      modified_list.append_unique(id1);
    }
    body_array[id2] = body2;
    temp_num = num_descendants(body2);
    if ( temp_num != body2_num )
    {
      modified_list.append_unique(id2);
    }
  }
  if ( modified_list.size() && print )
  {
    if ( modified_list.size() != 1 )
      PRINT_INFO("Imprinting Modified Bodies: ");
    else
      PRINT_INFO("Imprinting Modified Body: ");
    for ( i = 0; i < modified_list.size(); i++ )
    {
      if ( i != modified_list.size()-1 )
      {
        if (i%10 == 0 && i != 0 )
          PRINT_INFO("%d,\n", modified_list.get_and_step());
        else
          PRINT_INFO("%d, ", modified_list.get_and_step());
      }
      else
        PRINT_INFO("%d.\n", modified_list.get_and_step());
    }
  }
  else if ( modified_list.size() == 0 && print )
    PRINT_INFO("Imprinting overlaps resulted in no modifications.\n");

    //clean up our memory.
  for ( i = edges_list1.size(); i > 0; i-- )
  {
    ref_edges1 =  edges_list1.remove();
    ref_edges2 = edges_list2.remove();
    delete_edges_in_list(*ref_edges1);
    delete_edges_in_list(*ref_edges2);
    delete ref_edges1;
    delete ref_edges2;
  }
  delete [] body_array;
  return CUBIT_SUCCESS;
}

int SurfaceOverlapTool::num_descendants(Body *body)
{
  int counter = 0;
  counter += body->num_ref_volumes();
  counter += body->num_ref_faces();
  counter += body->num_ref_edges();
  counter += body->num_ref_vertices();
  return counter;
}

//---------------------------------------------------------
// copy_edges_in_list:
// copy's each of the edges in the old list, and puts the
// stand alone edges in the new list.
// NOTE: THESE ARE NEWLY CREATED EDGES AND SHOULD BE DELETED
//       BY THE CALLING FUNCTION.
// Author: David R. White
// Date: 1/11/02
//---------------------------------------------------------
CubitStatus SurfaceOverlapTool::copy_edges_in_list(DLIList<RefEdge*> &old_list,
                                                   DLIList<RefEdge*> &new_list )
{
  int i;
  RefEdge *curr_edge, *new_edge;
  RefEntity *curr_entity, *new_entity;
  
  for ( i = old_list.size(); i > 0; i-- )
  {
    curr_edge = old_list.get_and_step();
    curr_entity = CAST_TO(curr_edge, RefEntity);
    assert(curr_entity != NULL );
    new_entity = GeometryModifyTool::instance()->copy_refentity(curr_entity);
    new_edge = CAST_TO(new_entity, RefEdge);
    if ( new_edge == NULL )
    {
      PRINT_ERROR("Problems copying edges for imprinting overlap.\n");
      assert(new_edge != NULL );
      return CUBIT_FAILURE;
    }
    new_list.append(new_edge);
  }
  return CUBIT_SUCCESS;
}

//Initialize all settings in this class
void SurfaceOverlapTool::initialize_settings()
{
  SettingHandler::instance()->add_setting("Overlap Facet BBox Absolute", 
    SurfaceOverlapTool::set_facet_abs_tol, 
    SurfaceOverlapTool::get_facet_abs_tol); 
  
  SettingHandler::instance()->add_setting("Overlap Facet BBox Angle", 
    SurfaceOverlapTool::set_facet_ang_tol_setting, 
    SurfaceOverlapTool::get_facet_ang_tol_setting);
  
  SettingHandler::instance()->add_setting("Overlap Minimum Gap", 
    SurfaceOverlapTool::set_gap_min, 
    SurfaceOverlapTool::get_gap_min);  
  
  SettingHandler::instance()->add_setting("Overlap Maximum Gap",
    SurfaceOverlapTool::set_gap_max,
    SurfaceOverlapTool::get_gap_max);
  
  SettingHandler::instance()->add_setting("Overlap Minimum Angle", 
    SurfaceOverlapTool::set_angle_min, 
    SurfaceOverlapTool::get_angle_min); 
  
  SettingHandler::instance()->add_setting("Overlap Maximum Angle", 
    SurfaceOverlapTool::set_angle_max, 
    SurfaceOverlapTool::get_angle_max);
  
  SettingHandler::instance()->add_setting("Overlap Normal", 
    SurfaceOverlapTool::set_normal_type_setting, 
    SurfaceOverlapTool::get_normal_type_setting);  
  
  SettingHandler::instance()->add_setting("Overlap Tolerance",
    SurfaceOverlapTool::set_overlap_tolerance,
    SurfaceOverlapTool::get_overlap_tolerance);
  
  SettingHandler::instance()->add_setting("Overlap Group", 
    SurfaceOverlapTool::set_group_results_setting, 
    SurfaceOverlapTool::get_group_results_setting); 
  
  SettingHandler::instance()->add_setting("Overlap List", 
    SurfaceOverlapTool::set_list_pairs_setting, 
    SurfaceOverlapTool::get_list_pairs_setting);
  
  SettingHandler::instance()->add_setting("Overlap Display", 
    SurfaceOverlapTool::set_display_pairs_setting, 
    SurfaceOverlapTool::get_display_pairs_setting);
  
  SettingHandler::instance()->add_setting("Overlap Within Bodies", 
    SurfaceOverlapTool::set_check_within_bodies, 
    SurfaceOverlapTool::get_check_within_bodies);

  SettingHandler::instance()->add_setting("Overlap Across Bodies", 
    SurfaceOverlapTool::set_check_across_bodies, 
    SurfaceOverlapTool::get_check_across_bodies);
}

CubitStatus SurfaceOverlapTool::delete_edges_in_list(DLIList<RefEdge*> &edge_list )
{
  int i;
  RefEdge *curr_edge;
  RefEntity *curr_entity;
  
  for ( i = edge_list.size(); i > 0; i-- )
  {
    curr_edge = edge_list.remove();
    curr_entity = CAST_TO(curr_edge, RefEntity);
    assert(curr_entity != NULL );
    CubitStatus stat = GeometryQueryTool::instance()->delete_RefEntity(curr_entity);
    if ( stat != CUBIT_SUCCESS )
    {
      PRINT_ERROR("Problems deleting edges after imprinting overlaps\n");
      assert(0);
      return CUBIT_FAILURE;
    }
  }
  return CUBIT_SUCCESS;
}

CubitStatus SurfaceOverlapTool::find_overlapping_curves( DLIList<Curve*> &curve_list,
                                DLIList< DLIList<Curve*> *> &overlapping_curve_lists,
                                std::map<Curve*, DLIList<Curve*>* > &curve_to_list_map,
                                std::multimap<BodySM*, CubitVector> &body_point_imprint_map)
{
  int i;

  //put all the curves into a tree 
  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;

  // Populate the Surface AbstractTree
  AbstractTree<Curve*> *a_tree = new RTree<Curve*>( tolerance );
  curve_list.reset();
  for( i=curve_list.size(); i--; )
  {
    Curve *curve = curve_list.get_and_step();
    a_tree->add( curve );
  }

  std::map<Curve*, DLIList<CurveOverlapFacet*>* > facet_map;
  std::map<Curve*, DLIList<Curve*>* >::iterator list_iter; 

  curve_list.reset();
  for( i=curve_list.size(); i--; )
  {
    Curve *curve1= curve_list.get_and_step();

    BodySM *curve1_body = curve1->bodysm();

    // Remove this surface from AbstractTree so it is not found and never
    // found again
    a_tree->remove( curve1 );

    // Find RefFaces from AbstractTree that are within range of this surface
    CubitBox curve1_box = curve1->bounding_box();
    DLIList<Curve*> close_curves;
    a_tree->find( curve1_box, close_curves );

    int j;
    for( j=close_curves.size(); j--; )
    {
      Curve *curve2 = close_curves.get_and_step();
      BodySM *curve2_body = curve2->bodysm();
      
      if( curve2_body == curve1_body )
        continue;

      std::multimap<BodySM*, CubitVector> tmp_body_point_imprint_map;

      if( check_overlap( curve1, curve2, &facet_map, &tmp_body_point_imprint_map  ) ) 
      {
        //check to see if the curve1 is already overlapping with another curve 
        list_iter = curve_to_list_map.find( curve1 );
        if( list_iter != curve_to_list_map.end() )
        {
          DLIList<Curve*> *tmp_curve_list = list_iter->second;
          tmp_curve_list->append( curve2 );
        }
        else
        {
          //list for curve 1 does not exist....make a new one
          DLIList<Curve*> *tmp_curve_list = new DLIList<Curve*>;
          overlapping_curve_lists.append( tmp_curve_list );
          tmp_curve_list->append( curve1 );
          tmp_curve_list->append( curve2 );
          curve_to_list_map.insert( std::map<Curve*,
               DLIList<Curve*>*>::value_type( curve1, tmp_curve_list ));
        }
      }
      else
      {
        //append what's in body_point_imprint_map to one passed in
        body_point_imprint_map.insert( tmp_body_point_imprint_map.begin(), 
                                       tmp_body_point_imprint_map.end() ); 
      }
    }
  }

  //clean up facet map
  std::map<Curve*, DLIList<CurveOverlapFacet*>* >::iterator facet_iter;
  facet_iter=facet_map.begin(); 
  for(; facet_iter != facet_map.end(); facet_iter++ )
  {
    DLIList<CurveOverlapFacet*> *co_facet_list = facet_iter->second;

    //delete all the facets in the list
    for( i=co_facet_list->size(); i--; )
      delete co_facet_list->get_and_step();
    delete co_facet_list;
  }

  delete a_tree;

  return CUBIT_SUCCESS;
}

CubitStatus SurfaceOverlapTool::find_overlapping_curves( DLIList<Surface*> &surf_list,
                                DLIList< DLIList<Curve*> *> &overlapping_curve_lists,
                                std::map<Curve*, DLIList<Curve*>* > &curve_to_list_map,
                                std::multimap<BodySM*, CubitVector> &body_point_imprint_map)
{
  //collect all the surfaces
  DLIList<Curve*> curve_list;
  int i;
  for( i=surf_list.size(); i--; )
  {
    Surface *surf = surf_list.get_and_step();
    DLIList<Curve*> curves;
    surf->curves_ignore_virtual(curves, false);
//    body_sm->curves( curves);
    curve_list.merge_unique( curves );
  }

  return find_overlapping_curves(curve_list, overlapping_curve_lists, curve_to_list_map,
    body_point_imprint_map);
}

CubitStatus SurfaceOverlapTool::find_overlapping_curves( DLIList<BodySM*> &body_list,
                                DLIList< DLIList<Curve*> *> &overlapping_curve_lists,
                                std::map<Curve*, DLIList<Curve*>* > &curve_to_list_map,
                                std::multimap<BodySM*, CubitVector> &body_point_imprint_map )
{
  //collect all the surfaces
  DLIList<Curve*> curve_list;
  int i;
  for( i=body_list.size(); i--; )
  {
    BodySM *body_sm = body_list.get_and_step();
    DLIList<Curve*> curves;
    body_sm->curves_ignore_virtual(curves, false);
//    body_sm->curves( curves);
    curve_list.merge_unique( curves );
  }

  return find_overlapping_curves(curve_list, overlapping_curve_lists, curve_to_list_map,
    body_point_imprint_map);

}

CubitStatus SurfaceOverlapTool::find_overlapping_curves( DLIList<Body*> &bodies,
                                std::multimap<RefEdge*, RefEdge*> &overlapping_edge_map,
                                double maxgap /* =-1*/)
{
  std::map<RefEdge*, DLIList<CurveOverlapFacet*>* >::iterator list_iter; 

  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;
  if (maxgap != -1 && maxgap != 0.0)
      tolerance = maxgap;

  int i,j;
  DLIList<Body*> tmp_body_list = bodies;
  for(i=tmp_body_list.size(); i--; )
  {
    Body *tmp_body1 = tmp_body_list.pop();

    tmp_body_list.reset();
    for(j=tmp_body_list.size(); j--; )
    {
      Body *tmp_body2 = tmp_body_list.get_and_step();

      //determine if bodies are close enough to have overlapping curves
      CubitBox body_box1 = tmp_body1->bounding_box();
      CubitBox body_box2 = tmp_body2->bounding_box();

      if( body_box1.overlap( tolerance, body_box2 ) )
      {
        // create list of edges on each body.
	    // put the edges in body2 in edges1 so edges1 contains
		// the list of curves contained in the body listed first in the list
        DLIList<RefEdge*> edges1, edges2;
        tmp_body1->ref_edges( edges2 );
        tmp_body2->ref_edges( edges1 );
		
		//check individual curves for bounding box interference
        if (!find_overlapping_curves(edges1, edges2, overlapping_edge_map, tolerance))
			return CUBIT_FAILURE;
      }
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus SurfaceOverlapTool::find_overlapping_curves( DLIList<RefFace*> &faces,
                        std::multimap<RefEdge*, RefEdge*> &overlapping_edge_map,
                        double maxgap /*=-1*/)
{
  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;

  if (maxgap != -1 && maxgap != 0.0)
      tolerance = maxgap;

  // Populate the Surface AbstractTree
  AbstractTree<RefFace*> *a_tree = new RTree<RefFace*>( tolerance );
  faces.reset();
  int i;
  for( i=faces.size(); i--; )
  {
    RefFace* face = faces.get_and_step();
    a_tree->add( face );
  }

  faces.reset();
  for( i=faces.size(); i--; )
  {
    RefFace* face1 = faces.get_and_step();

    Body *face1_body = face1->body();

    // Remove this face from AbstractTree so it is not found and never
    // found again
    a_tree->remove( face1 );

    // Find RefFaces from AbstractTree that are within range of this surface
    CubitBox face1_box = face1->bounding_box();
    DLIList<RefFace*> close_faces;
    a_tree->find( face1_box, close_faces );

    int j;
    for( j=close_faces.size(); j--; )
    {
      RefFace *face2 = close_faces.get_and_step();
      Body *face2_body = face2->body();

      //don't check for overlapping faces within bodies
      if( face1_body == face2_body )
              continue;

      DLIList<RefEdge*> edges1, edges2;
      face1->ref_edges( edges1 );
      face2->ref_edges( edges2 );

      if (!find_overlapping_curves(edges1, edges2, overlapping_edge_map, tolerance))
              return CUBIT_FAILURE;
            
    }
  }
  return CUBIT_SUCCESS;
}

CubitStatus SurfaceOverlapTool::find_overlapping_curves( DLIList<RefEdge*> &edgelist,
				std::multimap<RefEdge*, RefEdge*> &overlapping_edge_map, double maxgap /*=-1*/ )
{
  // for now, check every edge against every other edge
  // (the helper function will skip edges if they are in the same body)

  DLIList<RefEdge*> edgelist2 = edgelist;
  return find_overlapping_curves(edgelist, edgelist2, overlapping_edge_map, maxgap);
}

CubitStatus SurfaceOverlapTool::find_overlapping_curves( DLIList<RefEdge*> &edges1, DLIList<RefEdge*> &edges2,
                                              std::multimap<RefEdge*, RefEdge*> &overlapping_edge_map,
                                              double maxgap /*=-1*/)
{
  int i;
  std::map<RefEdge*, DLIList<CurveOverlapFacet*>* > facet_map; 
  std::map<RefEdge*, DLIList<CurveOverlapFacet*>* >::iterator list_iter; 

  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;

  if (maxgap != -1 && maxgap != 0.0)
      tolerance = maxgap;

  for( i=edges1.size(); i--; )
  {
    RefEdge *edge1 = edges1.get_and_step(); 

    int j;
    for( j=edges2.size(); j--; )
    {
      RefEdge *edge2 = edges2.get_and_step();
      
      if( edge1 == edge2 )
        continue;

      bool pair_already_in_map = false;

      std::multimap<RefEdge*, RefEdge*>::iterator it;
      it = overlapping_edge_map.find( edge2 );
      if( it != overlapping_edge_map.end() )
        continue;

      // check to see if this pair is already in the overlap map.
      //  If it is, don't need to bother checking it for overlap
     /* 
      std::multimap<RefEdge*, RefEdge*>::iterator it;
      std::pair< std::multimap<RefEdge*,RefEdge*>::iterator, std::multimap<RefEdge*,RefEdge*>::iterator > range;

      range = overlapping_edge_map.equal_range(edge1);

      for (it=range.first; it!=range.second; it++)
      {
        if ( (it->second == edge2) && (it->first == edge1) )
                pair_already_in_map = true;
      }
      if (pair_already_in_map)
        continue; 

      */

      if( check_overlap( edge1, edge2, &facet_map, &tolerance ) ) 
      {
        overlapping_edge_map.insert( std::multimap<RefEdge*, RefEdge*>
                      ::value_type( edge1, edge2));
        //PRINT_INFO("Curve %d and %d overlap.\n", edge1->id(), edge2->id());
        //PRINT_INFO("Insert Curve %d, Curve %d.\n", edge1->id(), edge2->id());
      }
    }
  }
  //clean up facet map
  list_iter=facet_map.begin(); 
  for(; list_iter != facet_map.end(); list_iter++ )
  {
    DLIList<CurveOverlapFacet*> *co_facet_list = list_iter->second;
    //delete all the facets in the list
    for( i=co_facet_list->size(); i--; )
            delete co_facet_list->get_and_step();
    delete co_facet_list;
  }

  return CUBIT_SUCCESS;
}

CubitBoolean SurfaceOverlapTool::check_overlap( Curve *curve1, Curve *curve2, 
  std::map<Curve*, DLIList<CurveOverlapFacet*>* > *facet_map,
  std::multimap<BodySM*,CubitVector> *body_point_imprint_map ) 
{
  //if surfaces are not splines and are not of the same type, 
  //they won't overlap
  GeometryType curve1_type = curve1->geometry_type();
  GeometryType curve2_type = curve2->geometry_type();

  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;

  std::map<Curve*, DLIList<CurveOverlapFacet*>* >::iterator facet_iterator;

  CubitBox curve_box1 = curve1->bounding_box(); 
  CubitBox curve_box2 = curve2->bounding_box();

  //do bounding boxes overlap
  if( curve_box1.overlap( tolerance, curve_box2 ) )
  {
    //curves must overlap at least by 100th of the smaller curve's length
    double min_overlap = 0.0;
    double curve1_length = curve1->measure();
    double curve2_length = curve2->measure();

/*
    if( curve1_length < curve2_length )
      min_overlap = curve1->measure() * 0.01;
    else 
      min_overlap = curve2->measure() * 0.01;
      */

    min_overlap = tolerance;

    DLIList<CurveOverlapFacet*> *facet_list1 = NULL;
    DLIList<CurveOverlapFacet*> *facet_list2 = NULL;

    //see if curve is in map...if not we have to create faceting for it.
    if( facet_map )
      facet_iterator = facet_map->find( curve1 );

    if( facet_map == NULL || facet_iterator == facet_map->end() )
    {
      facet_list1 = new DLIList<CurveOverlapFacet*>;

      GMem curve_graphics;
      curve1->get_geometry_query_engine()->get_graphics( curve1, 
        &curve_graphics, 0.0 );

      GPoint *points = curve_graphics.point_list();
      int num_points = curve_graphics.pointListCount;

      int kk;
      for( kk=0; kk<num_points-1; kk++ )
      {
        //create a new CurveOverlapFacets
        GPoint gpoints[2];
        gpoints[0] = points[kk];
        gpoints[1] = points[kk + 1];

        CurveOverlapFacet *tmp_facet = new CurveOverlapFacet( gpoints ); 
        facet_list1->append( tmp_facet );
      }
      
      if( facet_map )
        facet_map->insert( std::map<Curve*, 
          DLIList<CurveOverlapFacet*>*>::value_type( curve1, facet_list1 )); 
      
    }
    else
      facet_list1 = (*facet_iterator).second;

    //see if edge is in map...if not we have to create facet for it.
    if( facet_map )
      facet_iterator = facet_map->find( curve2 );

    if( facet_map == NULL || facet_iterator == facet_map->end() )
    {
      facet_list2 = new DLIList<CurveOverlapFacet*>;

      GMem curve_graphics;
      curve2->get_geometry_query_engine()->get_graphics( curve2, 
        &curve_graphics, 0.0 );

      GPoint *points = curve_graphics.point_list();
      int num_points = curve_graphics.pointListCount;

      int kk;
      for( kk=0; kk<num_points-1; kk++ )
      {
        //create a new CurveOverlapFacets
        GPoint gpoints[2];
        gpoints[0] = points[kk];
        gpoints[1] = points[kk + 1];

        CurveOverlapFacet *tmp_facet = new CurveOverlapFacet( gpoints ); 
        facet_list2->append( tmp_facet );
      }
      
      if( facet_map )
        facet_map->insert( std::map<Curve*, 
          DLIList<CurveOverlapFacet*>*>::value_type( curve2, facet_list2 )); 
    }
    else
      facet_list2 = (*facet_iterator).second;
    
    //compare smaller list to larger list...want list1 to be smaller
    if( facet_list1->size() > facet_list2->size() )
    {
      DLIList<CurveOverlapFacet*> *tmp_list = facet_list1;
      facet_list1 = facet_list2;
      facet_list2 = tmp_list;
      Curve *tmp_curve = curve1;
      curve1 = curve2;
      curve2 = tmp_curve;
    }
    
    CubitVector curv1_start_pt;
    CubitVector curv2_start_pt;
    CubitVector curv1_end_pt;
    CubitVector curv2_end_pt;
    
    facet_list1->reset();
    facet_list2->reset();
    //get start/end facet points on curves
    if( body_point_imprint_map )
    {
      curv1_start_pt = facet_list1->get()->start_point();
      curv2_start_pt = facet_list2->get()->start_point();
      facet_list1->last();
      facet_list2->last();
      curv1_end_pt = facet_list1->get()->end_point();
      curv2_end_pt = facet_list2->get()->end_point();
    }
    
    facet_list1->reset();
    facet_list2->reset();
    int kk;
    double total_overlap = 0.0;
    bool overlap = false;
    //now determine if enough curve facets overlap

    int list_size1 = facet_list1->size();
    for( kk=list_size1; kk--; )
    {
      CurveOverlapFacet *curr_facet = facet_list1->get_and_step();
      if( curr_facet->length() < GEOMETRY_RESABS )
        continue;

      //if 'curr_facet' is the first or last facet in curve1,
      //check to see if the start/end point is on curve2, and not coincident
      //with the start/end point of curve2.  If so, it's a case where you need
      //to imprint this start/end point onto the body of curve2.
      if( body_point_imprint_map )
      {
        if( kk == list_size1-1 ) //first facet in curve1
        {
          CubitVector closest_point;
          curve2->closest_point_trimmed( curr_facet->start_point(), closest_point );
          if( curv2_start_pt.distance_between( curr_facet->start_point() ) > tolerance &&  
                curv2_end_pt.distance_between( curr_facet->start_point() ) > tolerance &&  
               closest_point.distance_between( curr_facet->start_point() ) < tolerance )
          {
            BodySM *tmp_body_sm = curve2->bodysm();
            body_point_imprint_map->insert( std::multimap<BodySM*, 
                                            CubitVector>::value_type(
                                            tmp_body_sm, closest_point)); 
          }
        }
        if( kk == 0 ) //last facet in curve1 
        {
          CubitVector closest_point;
          curve2->closest_point_trimmed( curr_facet->end_point(), closest_point );
          if( curv2_start_pt.distance_between( curr_facet->end_point() ) > tolerance &&  
                curv2_end_pt.distance_between( curr_facet->end_point() ) > tolerance &&  
               closest_point.distance_between( curr_facet->end_point() ) < tolerance ) 
          {
            //insert into vertex-volume imprint map
            BodySM *tmp_body_sm = curve2->bodysm();
            body_point_imprint_map->insert( std::multimap<BodySM*, 
                                            CubitVector>::value_type(
                                            tmp_body_sm, closest_point)); 
          }
        }
      } 

      //do bounding boxes of facets overlap?
      int jj;
      int list_size2 = facet_list2->size();
      for( jj=list_size2; jj--; )
      {
        CurveOverlapFacet *other_facet = facet_list2->get_and_step();

        if( curr_facet->bbox_overlap( tolerance, other_facet ) )
        {
          double distance_between_facets = curr_facet->facet_to_facet_distance( other_facet );
          if( distance_between_facets < tolerance )
          {
            if( other_facet->length() < GEOMETRY_RESABS ) 
              continue;

            //if 'other_facet' is the first or last facet in curve1,
            //check to see if the start/end point is on curve1, and not coincident
            //with the start/end point of curve1.  If so, it's a case where you need
            //to imprint this start/end point onto the body of curve1.
            if( body_point_imprint_map  )
            {
              if( jj == list_size2-1 ) //first facet in curve1
              {
                CubitVector closest_point;
                curve1->closest_point_trimmed( other_facet->start_point(), closest_point );
                
                if( curv1_start_pt.distance_between( other_facet->start_point() ) > tolerance &&  
                      curv1_end_pt.distance_between( other_facet->start_point() ) > tolerance &&  
                     closest_point.distance_between( other_facet->start_point() ) < tolerance )
                {
                  //insert into vertex-volume imprint map
                  BodySM *tmp_body_sm = curve1->bodysm();
                  body_point_imprint_map->insert( std::multimap<BodySM*, 
                                                  CubitVector>::value_type(
                                                  tmp_body_sm, closest_point )); 
                }
              }
              if( jj == 0 ) //last facet in curve1 
              {
                CubitVector closest_point;
                curve1->closest_point_trimmed( other_facet->end_point(), closest_point );
                if( curv1_start_pt.distance_between( other_facet->end_point() ) > tolerance &&  
                      curv1_end_pt.distance_between( other_facet->end_point() ) > tolerance &&  
                     closest_point.distance_between( other_facet->end_point() ) < tolerance )
                {
                  //insert into vertex-volume imprint map
                  BodySM *tmp_body_sm = curve1->bodysm();
                  body_point_imprint_map->insert( std::multimap<BodySM*, 
                                                  CubitVector>::value_type(
                                                  tmp_body_sm, closest_point )); 
                }
              }
            }

            //get the long and short facet edge
            double curr_facet_length = curr_facet->length();
            double other_facet_length = other_facet->length();
            CurveOverlapFacet *long_facet = ( curr_facet_length > other_facet_length ?
                                             curr_facet : other_facet );
            CurveOverlapFacet *short_facet = curr_facet;
            if( long_facet == curr_facet )
              short_facet = other_facet;

            //make sure both endpoints of the smaller facet edge lie within a radius
            //or merge tolerance to an infinite line defined by longer facet edge
            CubitVector direction = long_facet->end_point() - long_facet->start_point();
            double dist1 = short_facet->start_point().distance_from_infinite_line(
                                          long_facet->start_point(), direction );
            double dist2 = short_facet->end_point().distance_from_infinite_line(
                                          long_facet->start_point(), direction );

            if( dist1 <= tolerance && dist2 <= tolerance )
            {
              double overlap_tolerance = 
                curr_facet->length() < other_facet->length() ? curr_facet->length():
                other_facet->length();
              overlap_tolerance *= 0.01;
              
              //how much of the facet overlaps?
              double tmp_overlap = curr_facet->distance_overlapping( other_facet ) ;
              if( tmp_overlap > overlap_tolerance )
                total_overlap += tmp_overlap;

              if( total_overlap > min_overlap )
              {
                overlap = true;
                break;
              }
            }
          }
        }
      }
      if( overlap == true )
        break;
    } 

    if( facet_map == NULL )
    {
      //clean up facets and list
      for( kk=facet_list1->size(); kk--; )
        delete facet_list1->get_and_step();
      for( kk=facet_list2->size(); kk--; )
        delete facet_list2->get_and_step();

      delete facet_list1; 
      delete facet_list2; 
    }

    if( overlap == false )
      return CUBIT_FALSE;

    if( curve1_type == SPLINE_CURVE_TYPE  || 
        curve2_type == SPLINE_CURVE_TYPE ) 
    {
      //measure between the 2 curves
      double dist_between_curves = 0;
      CubitVector point_on_curve1, point_on_curve2;
      GeometryQueryTool::instance()->entity_entity_distance(curve1, curve2,
        point_on_curve1, point_on_curve2, dist_between_curves);

      if( dist_between_curves > tolerance )
          return CUBIT_FALSE;

      //the curvature at the same spot on the curve could be almost the same
      //if curvature's at midpoint are not the same, they don't overlap
      CubitVector curvature1, curvature2;
      CubitVector tangent, closest_location;
      curve1->closest_point( point_on_curve1, closest_location, &tangent, &curvature1 ); 
      curve2->closest_point( point_on_curve2, closest_location, &tangent, &curvature2 ); 

      double rad1 = curvature1.length();
      double rad2 = curvature2.length();
        
      //if curvatures are more than 10% off 
      double curvature_diff = fabs( rad1 - rad2 );

      if( rad1 > GEOMETRY_RESABS || rad2 > GEOMETRY_RESABS ) 
      {
        if( curvature_diff/( (fabs(rad1) + fabs(rad2))/2 ) > 0.1 )
          return CUBIT_FALSE;
      }
    }
    else if( curve1_type == ARC_CURVE_TYPE && 
             curve2_type == ARC_CURVE_TYPE ) 
    {
      //if both curves are arcs, make sure radii are almost the same
      CubitVector mid_point1, mid_point2;
      curve1->mid_point( mid_point1 );
      curve2->mid_point( mid_point2 );
      
      CubitVector dummy_vec;
      CubitVector curvature1;
      CubitVector curvature2;
      curve1->closest_point( mid_point1, dummy_vec, NULL, &curvature1 );
      curve2->closest_point( mid_point2, dummy_vec, NULL, &curvature2 );
      
      double rad1 = 1/curvature1.length();
      double rad2 = 1/curvature2.length();

      //if( fabs( rad1 - rad2 ) > 0.001 )
      if( fabs( rad1 - rad2 ) > tolerance )
        return CUBIT_FALSE; 
    }
    else
      return CUBIT_TRUE;
  }
  else
    return CUBIT_FALSE;
  
  return CUBIT_TRUE;
}

CubitBoolean SurfaceOverlapTool::check_overlap( RefEdge *edge1, RefEdge *edge2, 
  std::map<RefEdge*, DLIList<CurveOverlapFacet*>* > *facet_map,
  double *overlap_tol ) 
{
  //if edges are not splines and are not of the same type, 
  //they won't overlap
  if( (edge1->get_curve_ptr()->geometry_type() != SPLINE_CURVE_TYPE  &&
       edge2->get_curve_ptr()->geometry_type() != SPLINE_CURVE_TYPE) &&
      (edge1->get_curve_ptr()->geometry_type() != 
       edge2->get_curve_ptr()->geometry_type() )) 
    return CUBIT_FALSE;

  // we don't want to consider edges from the same body
  if (edge1->body() == edge2->body())
	  return CUBIT_FALSE;


  double tolerance = GeometryQueryTool::get_geometry_factor()*GEOMETRY_RESABS;
  if( overlap_tol )
    tolerance = *overlap_tol;
    

  std::map<RefEdge*, DLIList<CurveOverlapFacet*>* >::iterator facet_iterator;
  CubitBox edge_box1 = edge1->bounding_box(); 
  CubitBox edge_box2 = edge2->bounding_box();

  //do bounding boxes overlap
  if( edge_box1.overlap( tolerance, edge_box2 ) )
  {
    //curves must overlap at least by 100th of the smaller curve's length
    double min_overlap = 0.0;
    double edge1_length = edge1->measure();
    double edge2_length = edge2->measure();
    if( edge1_length < edge2_length )
      min_overlap = edge1->measure() * 0.01;
    else 
      min_overlap = edge2->measure() * 0.01;

    DLIList<CurveOverlapFacet*> *facet_list1 = NULL;
    DLIList<CurveOverlapFacet*> *facet_list2 = NULL;

    //see if edge is in map...if not we have to create faceting for it.
    if( facet_map )
      facet_iterator = facet_map->find( edge1 );

    if( facet_map == NULL || facet_iterator == facet_map->end() )
    {
      facet_list1 = new DLIList<CurveOverlapFacet*>;

      GMem curve_graphics;
      edge1->get_graphics( curve_graphics );

      GPoint *points = curve_graphics.point_list();
      int num_points = curve_graphics.pointListCount;

      int kk;
      for( kk=0; kk<num_points-1; kk++ )
      {
        //create a new CurveOverlapFacets
        GPoint gpoints[2];
        gpoints[0] = points[kk];
        gpoints[1] = points[kk + 1];

        CurveOverlapFacet *tmp_facet = new CurveOverlapFacet( gpoints ); 
        facet_list1->append( tmp_facet );
      }
      
      if( facet_map )
        facet_map->insert( std::map<RefEdge*, 
          DLIList<CurveOverlapFacet*>*>::value_type( edge1, facet_list1 )); 
    }
    else
      facet_list1 = (*facet_iterator).second;

    //see if edge is in map...if not we have to create facet for it.
    if( facet_map )
      facet_iterator = facet_map->find( edge2 );

    if( facet_map == NULL || facet_iterator == facet_map->end() )
    {
      facet_list2 = new DLIList<CurveOverlapFacet*>;

      GMem curve_graphics;
      edge2->get_graphics( curve_graphics );

      GPoint *points = curve_graphics.point_list();
      int num_points = curve_graphics.pointListCount;

      int kk;
      for( kk=0; kk<num_points-1; kk++ )
      {
        //create a new CurveOverlapFacets
        GPoint gpoints[2];
        gpoints[0] = points[kk];
        gpoints[1] = points[kk + 1];

        CurveOverlapFacet *tmp_facet = new CurveOverlapFacet( gpoints ); 
        facet_list2->append( tmp_facet );
      }
      
      if( facet_map )
        facet_map->insert( std::map<RefEdge*, 
          DLIList<CurveOverlapFacet*>*>::value_type( edge2, facet_list2 )); 
    }
    else
      facet_list2 = (*facet_iterator).second;
    
    //compare smaller list to larger list...want list1 to be smaller
    if( facet_list1->size() > facet_list2->size() )
    {
      DLIList<CurveOverlapFacet*> *tmp_list = facet_list1;
      facet_list1 = facet_list2;
      facet_list2 = tmp_list;
    }

    facet_list1->reset();
    facet_list2->reset();
    int kk;
    double total_overlap = 0.0;
    bool overlap = false;
    //now determine if enough curve facets overlap
    for( kk=facet_list1->size(); kk--; )
    {
      CurveOverlapFacet *curr_facet = facet_list1->get_and_step();

      //do bounding boxes of facets overlap?
      int jj;
      for( jj=facet_list2->size(); jj--; )
      {
        CurveOverlapFacet *other_facet = facet_list2->get_and_step();

        if( curr_facet->bbox_overlap( tolerance, other_facet ) )
        {
          if( curr_facet->facet_to_facet_distance( other_facet ) < tolerance )
          {

            //are facets parallel within some tolerance angle?
            double angle = curr_facet->angle( other_facet ); 
            if( angle < 1 || fabs(180-angle ) < 1 )
            {
              double overlap_tolerance = 
                curr_facet->length() < other_facet->length() ? curr_facet->length():
                other_facet->length();
              overlap_tolerance *= 0.01;
              
              //how much of the facet overlaps?
              double tmp_overlap = curr_facet->distance_overlapping( other_facet ) ;
              if( tmp_overlap > overlap_tolerance )
                total_overlap += tmp_overlap;

              if( total_overlap > min_overlap )
              {
                if( facet_map == NULL )
                {
                  //clean up facets and list
                  int i;
                  for( i=facet_list1->size(); i--; )
                    delete facet_list1->get_and_step();
                  for( i=facet_list2->size(); i--; )
                    delete facet_list2->get_and_step();

                  delete facet_list1; 
                  delete facet_list2; 
                }
                overlap = true;
                break;
              }
            }
          }
        }
      }
      if( overlap == true )
        break;
    }
    if( overlap == false )
      return CUBIT_FALSE;

    if( edge1->get_curve_ptr()->geometry_type() == SPLINE_CURVE_TYPE  || 
        edge2->get_curve_ptr()->geometry_type() == SPLINE_CURVE_TYPE ) 
    {
      //measure between the 2 curves
      double dist_between_edges = 0;
      CubitVector point_on_edge1, point_on_edge2;
      GeometryQueryTool::instance()->entity_entity_distance( edge1, edge2, 
                                                point_on_edge1, point_on_edge2, dist_between_edges );
    
      if( dist_between_edges > tolerance )
        return CUBIT_FALSE;

      //the curvature at the same spot on the curve could be almost the same
      //if curvature's at midpoint are not the same, they don't overlap
      CubitVector curvature1, curvature2;
      CubitVector tangent, closest_location;
      edge1->closest_point( point_on_edge1, closest_location, &tangent, &curvature1 ); 
      edge2->closest_point( point_on_edge2, closest_location, &tangent, &curvature2 ); 

      double rad1 = curvature1.length();
      double rad2 = curvature2.length();
        
      //if curvatures are more than 10% off 
      double curvature_diff = fabs( rad1 - rad2 );

      if( (rad1 || rad2 ) && 
         (rad1 > GEOMETRY_RESABS ||   //radii must be of significance even look at  
          rad2 > GEOMETRY_RESABS ))
      {
        if( curvature_diff/( (fabs(rad1) + fabs(rad2))/2 ) > 0.1 )
          return CUBIT_FALSE;
      }
    }
    else
      return CUBIT_TRUE;
  }
  else
    return CUBIT_FALSE;
  
  return CUBIT_TRUE;
}
