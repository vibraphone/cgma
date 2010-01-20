//-------------------------------------------------------------------------
// Filename      : FacetorTool.cpp
//
// Purpose       : 2D Delaunay mesher, template that can use both mesh and geom entities, not yet tested with geometric
//                 entities
//
// Creator       : Christopher Hynes
//
// Creation Date : 5/31/2002
//
// Owner         : Steve Owen
//-------------------------------------------------------------------------

#ifdef INLINE_TEMPLATES
#define MY_INLINE inline
#else
#define MY_INLINE
#endif

#include "FacetorTool.hpp"
#include "FacetorUtil.hpp"
#include "TDDelaunay.hpp"
#include "BoundaryConstrainTool.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "FacetEvalTool.hpp"
#include "ParamTool.hpp"

#define DETERM(p1,q1,p2,q2,p3,q3)\
     ((q3)*((p2)-(p1)) + (q2)*((p1)-(p3)) + (q1)*((p3)-(p2)))
#define SQR(x) ((x) * (x))
#define FT_INSIDE_TOL 1.0e-6
#define QUALITY_ANGLE 0.361283155162   /* 20.7 degrees */
#define ALPHA 0.70228615

//-------------------------------------------------------------------------
// Function:    FacetorTool
// Description: constructor
// Author:      chynes
// Date:        6/5/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE 
FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::
FacetorTool( SURF *ref_face_ptr, 
			DLIList<NODE *> &bounding_nodes, 
			NODE **boundary_edge_start_nodes,
			NODE **boundary_edge_end_nodes, 
			int num_boundary_edges, 
			SIZEFUNC *sizing_function,
			ParamTool *p_tool )
{
  //update private variables
  pTool = p_tool;
  refFacePtr = ref_face_ptr;
  boundingNodes = &bounding_nodes;
  boundaryEdgeStartNodes = boundary_edge_start_nodes;
  boundaryEdgeEndNodes = boundary_edge_end_nodes;
  numBoundaryEdges = num_boundary_edges;
  sizingFunction = sizing_function;

  curVisitFlag = INT_MIN+1;
  boxNodes[0] = boxNodes[1] = boxNodes[2] = boxNodes[3] = NULL;
}


//-------------------------------------------------------------------------
// Function:    ~FacetorTool
// Description: destructor
// Author:      chynes
// Date:        5/31/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::~FacetorTool()
{
  numBoundaryEdges = 0;
}


//-------------------------------------------------------------------------
// Function:    mesh_surf
// Description: mesh the surface
// Author:      chynes
// Date:        5/31/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::mesh_surf(DLIList<TRI *> &facet_list)
{
  CubitStatus  ret_value  = CUBIT_SUCCESS;
  
  // update private variables
  facetList = &facet_list;

  // create two initial bounding triangles

  if (ret_value == CUBIT_SUCCESS)
  {
    ret_value = init_box();
  }

  // insert the boundary nodes

  if (ret_value == CUBIT_SUCCESS)
  {
      //PRINT_INFO("Inserting boundary nodes \n");
    ret_value = insert_nodes( boundingNodes );
  }

  // constrain the boundary

  if (ret_value == CUBIT_SUCCESS)
  {
      //PRINT_INFO("Constraining the mesh to the boundary\n");
    ret_value = constrain_boundary();
  }

  // delete the triangles on the outside of the boundary

  if (ret_value == CUBIT_SUCCESS)
  {
      //PRINT_INFO("Deleting the exterior entities\n");
    ret_value = delete_exterior();
  }

  // insert the hards

  // insert interior points

  if (ret_value == CUBIT_SUCCESS) {
      //PRINT_INFO("Refining the interior mesh\n");
    ret_value = refine_interior(*facetList);
  }

    // clean up
    //PRINT_INFO("Cleaning up some data\n");
  clean_up_data();
  

  return ret_value;
}

//-------------------------------------------------------------------------
// Function:    mesh_surfwoIP
// Description: mesh the surface w/o interior point refinement
// Author:      chynes
// Date:        6/20/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::mesh_surfwoIP(DLIList<TRI *> &facet_list)
{
  CubitStatus  ret_value  = CUBIT_SUCCESS; 

  //update private variables
  facetList = &facet_list;

  // create two initial bounding triangles

  if (ret_value == CUBIT_SUCCESS)
  {
    ret_value = init_box();
  }

  // insert the boundary nodes

  if (ret_value == CUBIT_SUCCESS)
  {
    ret_value = insert_nodes( boundingNodes );
  }

  // constrain the boundary

  if (ret_value == CUBIT_SUCCESS)
  {
    ret_value = constrain_boundary();
  }

  // delete the triangles on the outside of the boundary

  if (ret_value == CUBIT_SUCCESS)
  {
    ret_value = delete_exterior();
  }

  // insert the hard points

  // clean up

  clean_up_data();

  return ret_value;
}

//-------------------------------------------------------------------------
// Function:    mesh_surfwoIPBC
// Description: mesh the surface w/o interior point refinement
// Author:      chynes
// Date:        6/20/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::mesh_surfwoIPBC(DLIList<TRI *> &facet_list)
{
  CubitStatus  ret_value  = CUBIT_SUCCESS;

  // update private variables

  facetList = &facet_list;

  // create two initial bounding triangles

  if (ret_value == CUBIT_SUCCESS)
  {
    ret_value = init_box();
  }

  // insert the boundary nodes

  if (ret_value == CUBIT_SUCCESS)
  {
    ret_value = insert_nodes( boundingNodes );
  }

  // delete the triangles on the outside of the boundary

  if (ret_value == CUBIT_SUCCESS){
    ret_value = delete_exterior();
  }

  // insert the hard points

  // clean up

  clean_up_data();

  return ret_value;
}

//-------------------------------------------------------------------------
// Function:    init_box
// Description: create two initial bounding triangles
// Author:      chynes
// Date:        5/31/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::init_box(void)
{
  // find the bounding box

  int ii;
  CubitVector min_box, max_box;
  min_box.x( CUBIT_DBL_MAX );min_box.y( CUBIT_DBL_MAX );
  max_box.x( -CUBIT_DBL_MAX );max_box.y( -CUBIT_DBL_MAX );

  NODE *node_ptr;
  CubitVector node_loc;
  for (ii=0; ii<boundingNodes->size(); ii++)
  {
    node_ptr = boundingNodes->get_and_step();
    node_loc = node_ptr->coordinates();
    if (node_loc.x() < min_box.x()) min_box.x( node_loc.x() );    
    if (node_loc.y() < min_box.y()) min_box.y( node_loc.y() );
    if (node_loc.x() > max_box.x()) max_box.x( node_loc.x() );    
    if (node_loc.y() > max_box.y()) max_box.y( node_loc.y() );
  }

  // expand the box by 10%

  double dx = max_box.x() - min_box.x();
  double dy = max_box.y() - min_box.y();
  double expand;
  if (dx > dy)
    expand = 0.1 * dx;
  else 
    expand = 0.1 * dy;

  min_box.x( min_box.x() - expand );
  min_box.y( min_box.y() - expand );
  max_box.x( max_box.x() + expand );
  max_box.y( max_box.y() + expand );
 
  
   // create four new nodes

  boxNodes[0] = (NODE *) new NODECHILD(min_box.x(), min_box.y(), 0.0, refFacePtr);
  boxNodes[1] = (NODE *) new NODECHILD(max_box.x(), min_box.y(), 0.0, refFacePtr);
  boxNodes[2] = (NODE *) new NODECHILD(max_box.x(), max_box.y(), 0.0, refFacePtr);
  boxNodes[3] = (NODE *) new NODECHILD(min_box.x(), max_box.y(), 0.0, refFacePtr);

  
  // create the two triangles

  TRI *new_facet1 = (TRI *) new TRICHILD( boxNodes[0], boxNodes[1], 
                                     boxNodes[3], refFacePtr );
  TRI *new_facet2 = (TRI *) new TRICHILD( boxNodes[1], boxNodes[2], 
                                     boxNodes[3], refFacePtr );

  facetList->append(new_facet1);
  facetList->append(new_facet2);

  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Function:    insert_nodes
// Description: insert nodes into Delaunay mesh
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::
insert_nodes(DLIList<NODE *> *&bounding_nodes)
{
  CubitStatus rv = CUBIT_SUCCESS;
  int ii;
  NODE *point_ptr;
  TRI* start_tri = NULL;
  for (ii=0; ii<bounding_nodes->size() && rv == CUBIT_SUCCESS; ii++){
    point_ptr = bounding_nodes->get_and_step();
    rv = FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::insert_node(
                                                                point_ptr, *facetList,
                                                                refFacePtr, curVisitFlag,
                                                                start_tri);
  }

  
  return rv;
}


//-------------------------------------------------------------------------
// Function:    constrain_boundary
// Description: recover the boundary edges from the mesh
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::constrain_boundary(void)
{
  CubitStatus rv = CUBIT_SUCCESS;

  BoundaryConstrainTool<SURF,TRI,EDGE,NODE,TRICHILD> bctool( refFacePtr );
  EDGE *edge_ptr = NULL;
  int ii;
  for (ii=0; ii<numBoundaryEdges && rv == CUBIT_SUCCESS; ii++)
  {
    rv = bctool.recover_edge( boundaryEdgeStartNodes[ii],
			      boundaryEdgeEndNodes[ii],
			      edge_ptr,
			      facetList);
    //assert(edge_ptr != NULL);
    if (edge_ptr == NULL)
    {
      PRINT_ERROR("Could not recover a boundary in Delaunay facetor.\n"
                  "Indicates potential problems with input loop geometry\n");
      rv = CUBIT_FAILURE;
    }
  }
  return rv;
}


//-------------------------------------------------------------------------
// Function:    delete_exterior
// Description: delete the triangles on the outside of the boundary and on
//              the interior of holes
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::delete_exterior(void)
{
  CubitStatus rv = CUBIT_SUCCESS;

  // mark edges at the boundary

  EDGE *edge_ptr;
  int ii;
  for (ii=0; ii<numBoundaryEdges; ii++)
  {
    edge_ptr = boundaryEdgeStartNodes[ii]->shared_edge( boundaryEdgeEndNodes[ii] );
    if(!edge_ptr){
      PRINT_ERROR("Boundary edges were not successfully recovered.\n");
      return CUBIT_FAILURE;
    }
    edge_ptr->marked( CUBIT_TRUE );
  }

  // get a tri adjacent to a bounding box node

  DLIList<TRI*> adjtris;
  boxNodes[0]->tris( adjtris );
  TRI *tri_ptr = adjtris.get();

  // mark the tris on the outside starting from this triangle

  int num_edges_found=0;
  rv = mark_tris( tri_ptr, num_edges_found );
  if (rv != CUBIT_SUCCESS)
    return rv;

  // delete the exterior triangles

  int marked_flag = 1;
  rv = delete_tris( *facetList, marked_flag );

  // check to see if we've found all the edges
  //PRINT_INFO("num_edges_found = %i, numBoundaryEdges = %i\n",num_edges_found,numBoundaryEdges);
  
  if (num_edges_found != numBoundaryEdges)
  {

    // if not, then there are interior loops(holes) that need to have their
    // triangles removed

    // find a triangle at the boundary

    TRI *adjtri;
    int jj;
    CubitBoolean found = CUBIT_FALSE;
    if(!facetList->size()){
      PRINT_ERROR("Problem marking exterior facets.\n");
      return CUBIT_FAILURE;
    }
    
    for (ii=0; ii<facetList->size() && !found; ii++)
    {
      tri_ptr = facetList->get_and_step();
      for (jj=0; jj<3 && !found; jj++)
      {
	int kk = jj;
	adjtri = tri_ptr->adjacent( kk, refFacePtr );
        if (adjtri == NULL)
	  {
          found = CUBIT_TRUE;
        }
      }
    }

    // mark all the tris we want to keep
    if(!found){
      PRINT_WARNING("Possible problem.  No boundary edge found.\n");
    }
    
    rv = mark_tris( tri_ptr, num_edges_found );

    // delete all the rest

    marked_flag = 0;
    rv = delete_tris( *facetList, marked_flag );
  }

  return rv;
}


//-------------------------------------------------------------------------
// Function:    delete_tris
// Description: delete triangles in the list with the specified marked flag
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::delete_tris(DLIList<TRI *> &tri_list, int marked_flag)
{
  CubitStatus rv = CUBIT_SUCCESS;

  int ii, jj;
  EDGE *edge_ptr;
  TRI *tri_ptr;
  DLIList<EDGE *> edge_list;
  DLIList<TRI *> new_list;
  int ntri = tri_list.size();
  for(ii=0; ii<ntri; ii++)
  {
    tri_ptr = tri_list.pop();
    if (tri_ptr->marked() == marked_flag)
    {
      edge_list.clean_out();
      tri_ptr->edges( edge_list );
      delete tri_ptr;
      for(jj=0; jj<3; jj++)
      {
        edge_ptr = edge_list.get_and_step();
        if(edge_ptr->number_tris() == 0 && edge_ptr->number_faces() == 0)
          delete edge_ptr;
      }
    }
	else
	  new_list.append(tri_ptr);
  }
  tri_list+=new_list;
  return rv;
}


//-------------------------------------------------------------------------
// Function:    mark_tris
// Description: recursive function to mark all the triangles we want to keep
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::mark_tris(TRI *tri_ptr, 
                                       int &num_edges_found)
{
  CubitStatus rv = CUBIT_SUCCESS;
  tri_ptr->marked( CUBIT_TRUE );
  int ii;
  EDGE *edge_ptr;
  TRI *adjtri;
  for (ii=0; ii<3 && rv == CUBIT_SUCCESS; ii++)
  {
    int jj = ii;
    adjtri = tri_ptr->adjacent( jj, refFacePtr );
    edge_ptr = tri_ptr->edge( jj );
    if(edge_ptr->marked()) {
      num_edges_found++;
    }
    else {
      if (adjtri != NULL && !adjtri->marked())
        rv = mark_tris( adjtri, num_edges_found );
    }
  }
  return rv;
}


//-------------------------------------------------------------------------
// Function:    refine_interior
// Description: generate nodes on the interior of the tiangulation to
//              improve quality
// Notes:       this inserts nodes at the circumcenters of triangles roughly
//              based on the algorithm by Jonathon Shewchuk
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::refine_interior(DLIList<TRI *> &tri_list)
{
  CubitStatus rv = CUBIT_SUCCESS;

  // classify the triangles based on their minimum angle

  //rv =
  const int num_lists = 64;
  const double interval = QUALITY_ANGLE / double( num_lists - 1 );
  DLIList<TRI*>* tri_sort_array = new DLIList<TRI *> [num_lists];
  classify_triangles(tri_list, tri_sort_array, num_lists, interval );

  // process each of the triangles until done

  TRI *tri_ptr;
  TRI* start_tri = NULL;
  while ((tri_ptr = next_triangle(tri_sort_array, num_lists)) != NULL && rv == CUBIT_SUCCESS) 
  {
    rv = FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::insert_at_circumcenter(tri_ptr,
                                                                                             *facetList,
                                                                                             start_tri,
                                                                                             curVisitFlag,
                                                                                             refFacePtr,
                                                                                             tri_sort_array,
                                                                                             num_lists,
                                                                                             interval,
                                                                                             QUALITY_ANGLE,
                                                                                             sizingFunction,
                                                                                             pTool);
  }

  delete [] tri_sort_array;
  tri_sort_array = NULL;

  return rv;
}

//-------------------------------------------------------------------------
// Function:    classify_triangles
// Description: order the triangles in the current mesh by their worst angle
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::classify_triangles(DLIList<TRI *> &tri_list,
                                                                                                  DLIList<TRI*>* sorted_lists,
                                                                                                  const int num_lists,
                                                                                                  const double interval)
{
  CubitStatus rv = CUBIT_SUCCESS;

  // Create an array of lists that will hold triangles as they are
  // created (and deleted).  Each list will hold triangles whose angles 
  // fall within a specified interval.  Triangles will be processed
  // from smallest angle to best angle

  // classify each trriangle and place into sort lists

  TRI *tri_ptr;
  int ii;
  for(ii=0; ii<tri_list.size() && rv == CUBIT_SUCCESS; ii++)
  {
    tri_ptr = tri_list.get_and_step();
    //rv = 
    FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::classify_tri_by_angle(
                                    tri_ptr, sorted_lists,  num_lists, interval, QUALITY_ANGLE);
  }
   
  return rv;
}


//-------------------------------------------------------------------------
// Function:    next_triangle
// Description: get the next triangle to process
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
TRI *FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::next_triangle(DLIList<TRI*>* sorted_lists,
                                                                                      const int num_lists)
{
  TRI *tri_ptr = NULL;
  int ii;
  for( ii = 1; ii < num_lists && tri_ptr == NULL; ii++)
  {
    if (sorted_lists[ii].size() > 0)
      tri_ptr = sorted_lists[ii].remove();
  }
  if (tri_ptr == NULL)
  {
    if (sorted_lists[0].size() > 0)
      tri_ptr = sorted_lists[0].remove();
  }
  return tri_ptr;
}

//-------------------------------------------------------------------------
// Function:    clean_up_data
// Description: clean up any data we've allocated
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
void FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::clean_up_data(void)
{

  EDGE *del_edge_ptr;
  DLIList<EDGE *> del_edge_list;
  TRI *tri_ptr;
  DLIList<TRI *> tri_list;
  
  
  int ii, jj, kk;

    //loop over facets and remove marks
  for(ii=0; ii<facetList->size(); ii++)
  {	
    tri_ptr = facetList->get_and_step();
    tri_ptr->marked(CUBIT_FALSE);
    ToolData *td = tri_ptr->remove_TD( TDDelaunay< TRI, NODE >::is_delaunay );
    if (td != NULL)
    {
      TDDelaunay< TRI, NODE > *td_del = dynamic_cast<TDDelaunay< TRI, NODE >*> (td);
      delete td_del;
    }
  }

    //loop over the box ndes
    //if meshing was successful, some of the below is unnecessary
  for(ii=0; ii<4; ii++)
  {
    if (boxNodes[ii] != NULL){
      tri_list.clean_out();
        //if there are tris attached, delete them.
      boxNodes[ii]->tris(tri_list);
      for (jj=0; jj<tri_list.size(); jj++)
      {
        tri_ptr = tri_list.get_and_step();
        del_edge_list.clean_out();
        tri_ptr->edges( del_edge_list );
        delete tri_ptr;
        
          // also deltet the attached, unused edges
        for (kk=0; kk<del_edge_list.size(); kk++)
        {
          del_edge_ptr = del_edge_list.get_and_step();
          if (del_edge_ptr->number_tris() == 0)
            delete del_edge_ptr;
        }
      }
    }
      //finally delete the nodes
    delete boxNodes[ii];
  }

  // clean off the edge marks

  EDGE *edge_ptr;
  for (ii=0; ii<numBoundaryEdges; ii++)
  {
    edge_ptr = boundaryEdgeStartNodes[ii]->shared_edge( boundaryEdgeEndNodes[ii] );
    if (edge_ptr)
      edge_ptr->marked( CUBIT_FALSE );
  }
}

// EOF
 
