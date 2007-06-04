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
#define NUM_SORT_LISTS 64
#define QUALITY_ANGLE 0.361283155162   /* 20.7 degrees */
#define INTERVAL (QUALITY_ANGLE/(NUM_SORT_LISTS-1))
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
  lastTriLocated = NULL;
  boxNodes[0] = boxNodes[1] = boxNodes[2] = boxNodes[3] = NULL;
  triSortArray = NULL;
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
  if (triSortArray != NULL)
    delete [] triSortArray;
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
// Function:    circumcenter
// Description: get the circumcenter of the triangle
// Author:      chynes
// Date:        6/6/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitVector &FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::circumcenter( TRI *tri_ptr )
{
	ToolData *td = tri_ptr->get_TD( TDDelaunay< TRI, NODE >::is_delaunay );
	TDDelaunay< TRI, NODE > *td_del = dynamic_cast<TDDelaunay< TRI, NODE >*> (td);
	if(td_del == NULL) {
		td_del = new TDDelaunay<TRI, NODE>();
		tri_ptr->add_TD( td_del );
	}
	return td_del->circumcenter2d( tri_ptr );
}


//-------------------------------------------------------------------------
// Function:    radius
// Description: get the radius squared of the triangle circumcircle
// Author:      chynes
// Date:        6/6/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
double FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::radius( TRI *tri_ptr )
{
  ToolData *td = tri_ptr->get_TD( TDDelaunay< TRI, NODE >::is_delaunay );
  TDDelaunay< TRI, NODE > *td_del = dynamic_cast<TDDelaunay< TRI, NODE >*> (td);
  if (td_del == NULL)
  {
    td_del = new TDDelaunay< TRI, NODE >();
    tri_ptr->add_TD( td_del );
  }
  return td_del->radius2d( tri_ptr );
}

//-------------------------------------------------------------------------
// Function:    tri_visited
// Description: set the tri_visited flag
// Author:      chynes
// Date:        6/6/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitBoolean FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::tri_visited( TRI *tri_ptr )
{
  ToolData *td = tri_ptr->get_TD( TDDelaunay< TRI, NODE >::is_delaunay );
  TDDelaunay< TRI, NODE > *td_del = dynamic_cast<TDDelaunay< TRI, NODE >*> (td);
  if (td_del == NULL)
  {
    td_del = new TDDelaunay< TRI, NODE >();
    tri_ptr->add_TD( td_del );
  }
  return (td_del->visit_flag() == curVisitFlag);
}

//-------------------------------------------------------------------------
// Function:    tri_visited
// Description: set the tri_visited flag
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
void FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::
tri_visited( TRI *facet_ptr, CubitBoolean visited )
{
  ToolData *td = facet_ptr->get_TD( TDDelaunay< TRI, NODE >::is_delaunay );
  TDDelaunay< TRI, NODE > *td_del = dynamic_cast<TDDelaunay< TRI, NODE >*> (td);
  if (td_del == NULL)
  {
    td_del = new TDDelaunay< TRI, NODE >();
    facet_ptr->add_TD( td_del );
  }
  if (visited)
    td_del->visit_flag(curVisitFlag);
  else
    td_del->visit_flag(INT_MIN);
}


//-------------------------------------------------------------------------
// Function:    tri_sort_list
// Description: set the tri sort list index
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
void FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::tri_sort_list( TRI *facet_ptr, 
                                     int sort_list_index )
{
  ToolData *td = facet_ptr->get_TD( TDDelaunay< TRI, NODE >::is_delaunay );
  TDDelaunay< TRI, NODE > *td_del = dynamic_cast<TDDelaunay< TRI, NODE >*> (td);
  if (td_del == NULL)
  {
    td_del = new TDDelaunay<TRI, NODE>();
    facet_ptr->add_TD( td_del );
  }
  td_del->tri_sort_list(sort_list_index);
}


//-------------------------------------------------------------------------
// Function:    tri_sort_list
// Description: get the tri sort list index
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
int FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::tri_sort_list( TRI *facet_ptr )
{
  ToolData *td = facet_ptr->get_TD( TDDelaunay< TRI, NODE >::is_delaunay );
  TDDelaunay< TRI, NODE > *td_del = dynamic_cast<TDDelaunay< TRI, NODE >*> (td);
  if (td_del == NULL)
  {
    td_del = new TDDelaunay<TRI, NODE>();
    facet_ptr->add_TD( td_del );
  }
  return td_del->tri_sort_list();
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
  for (ii=0; ii<bounding_nodes->size() && rv == CUBIT_SUCCESS; ii++){
    point_ptr = bounding_nodes->get_and_step();
    rv = insert_node( point_ptr );
  }

  
  return rv;
}


//-------------------------------------------------------------------------
// Function:    insert_node
// Description: insert one node into Delaunay mesh
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::insert_node( NODE *node_ptr )
{
  CubitStatus rv = CUBIT_SUCCESS;

  // get a list of all triangles whose circumcircle contain the point
  
  DLIList<TRI *> tri_list;
  CubitVector the_point = node_ptr->coordinates();
  rv = natural_neighbor_tris( the_point, tri_list );


  // Use a Bowyer-Watson insertion 
  
  if (rv == CUBIT_SUCCESS)
  {
    rv = bowyer_watson_insert( node_ptr, tri_list );
  }

  return rv;
}


//-------------------------------------------------------------------------
// Function:    bowyer_watson_insert
// Description: Bowyer-Watson insertion into an existing Delaunay Mesh
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::bowyer_watson_insert( NODE *point_ptr,
					       DLIList<TRI *> &tri_list)
{
  CubitStatus rv = CUBIT_SUCCESS;

   // mark the tris in the list so we can distinguish them from their 
  // neighbors

  curVisitFlag++;
  int ii, jj;
  TRI *tri_ptr;
  for (ii=0; ii<tri_list.size(); ii++)
  {
    tri_ptr = tri_list.get_and_step();
    tri_visited( tri_ptr, CUBIT_TRUE );
  }
    //MBREWER:  This is not an optimal test.  But, we need need to 
    //do some tests to try to ensure that the void is valid for what
    //we need.  This is attempting to avoid crashes by not allowing nodes
    //to be inserted when the mesh starts diverging from the Delaunay
    //criteria.
  rv = valid_void( point_ptr, tri_list );
  if(!rv)
    return rv;
  
  // find all edges at the boundary of the visited triangles and create
  // new triangles with them

  // create a new triangle with this edge and the node
  TRI *adjtri_ptr;
  TRI *new_tri = NULL;
  EDGE *edge_ptr;
  DLIList<EDGE *> edge_list;
  for (ii=0; ii<tri_list.size(); ii++)
  {
    tri_ptr = tri_list.get_and_step();
    for (jj=0; jj<3; jj++){
      
      int kk = jj;
        // - if TRI == CubitFacet
        // - kk will be corrected in adjacent() to 
        // - give the correct EDGE index
      adjtri_ptr = tri_ptr->adjacent( kk, refFacePtr );
      if (!adjtri_ptr || !tri_visited( adjtri_ptr ))
      {
        edge_ptr = tri_ptr->edge(kk);
        assert(edge_list.append_unique( edge_ptr ));
        if(tri_ptr->sense(kk) == CUBIT_FORWARD)
          new_tri = (TRI *) new TRICHILD( edge_ptr->start_node(), edge_ptr->end_node(), point_ptr, refFacePtr);
        else
          new_tri = (TRI *) new TRICHILD( edge_ptr->end_node(), edge_ptr->start_node(), point_ptr, refFacePtr);
        facetList->append(new_tri);
      }
    }
  }
  lastTriLocated = new_tri;

  // delete the triangles in the original triangle list

  EDGE *del_edge_ptr;
  DLIList<EDGE *> del_edge_list;
  for (ii=0; ii<tri_list.size(); ii++)
  {
    tri_ptr = tri_list.get_and_step();
    del_edge_list.clean_out();
    facetList->move_to(tri_ptr);
    facetList->extract();
    tri_ptr->edges( del_edge_list );
    delete tri_ptr;

      // delete the unused edges
    for (jj=0; jj<del_edge_list.size(); jj++)
    {
      del_edge_ptr = del_edge_list.get_and_step();
      if (del_edge_ptr->number_tris() == 0 && del_edge_ptr->number_faces() == 0 )
        delete del_edge_ptr;
    }
  }

  return rv;
}
/************************************************************/
//author:  mbrewer
//This function performs tests on the void created by a point's
//insertion to make sure that invalid connectivities are not
//being created.  
/************************************************************/
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::valid_void( NODE * /*point_ptr*/,
					       DLIList<TRI *> &tri_list)
{
  int temp_i, temp_j;
  DLIList<EDGE*> boundary_edge_list;
  DLIList<NODE*> boundary_node_list;
  TRI *adjtri_ptr;
  TRI *tri_ptr;
  EDGE *edge_ptr;
  DLIList<EDGE *> edge_list;
    //loop over the tri's in tri_list and find all of the curves
    //on the boundary of the set (ie, on the boundary of the void). 
  for (temp_i=0; temp_i<tri_list.size(); temp_i++)
  {
    tri_ptr = tri_list.get_and_step();
      //check each edge to see whether it is a boundary edge or not
    for (temp_j=0; temp_j<3; temp_j++){
      
      int kk = temp_j;
        // - if TRI == CubitFacet
        // - kk will be corrected in adjacent() to 
        // - give the correct EDGE index
      adjtri_ptr = tri_ptr->adjacent( kk, refFacePtr );
      if (!adjtri_ptr || !tri_visited( adjtri_ptr ))
      {
        edge_ptr = tri_ptr->edge(kk);
        boundary_edge_list.append(edge_ptr);
      }
    }
  }
  int list_size = boundary_edge_list.size();
    //uniquify the boundary edge list
  boundary_edge_list.uniquify_unordered();
    //the list shouldn't have changed size during the uniquify.
    //if it did, there is a problem.
  if(list_size != boundary_edge_list.size()){
    PRINT_WARNING("Unexpected result.  Edge was duplicated on boundary.\n");
    return CUBIT_FAILURE;
  }
    //loop over the boundary edges and get the nodes in the boundary loop
  for(temp_i=0; temp_i<list_size; ++temp_i){
    edge_ptr=boundary_edge_list.get_and_step();
    boundary_node_list.append(edge_ptr->end_node());
    boundary_node_list.append(edge_ptr->start_node());
  }
  list_size = boundary_node_list.size();
    //each node should be in exactly two edges.  First make sure that there
    //are the correct number of nodes.
  if(list_size%2){
    PRINT_WARNING("Unexpected result.  Node not listed twice.\n");
    return CUBIT_FAILURE;
  }
    //now uniquify to make sure that the nodes were listed the correct number
    //of times.  
  boundary_node_list.uniquify_unordered();
  if( (list_size/2) != boundary_node_list.size()){
    PRINT_WARNING("Unexpected result.  Node was listed an incorrect number of times.\n");
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
  
}

//-------------------------------------------------------------------------
// Function:    natural_neighbor_tris
// Description: get a list of all triangles whose circumcircle contain 
//              the point
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::natural_neighbor_tris(
  CubitVector &the_point,
  DLIList <TRI *> &tri_list)
{
  CubitStatus rv = CUBIT_SUCCESS;

  // find the triangle the point is contained in

  //CubitVector areacoord;
  TRI *tri_ptr;
  rv = locate_point( the_point, tri_ptr );

  // keep track of visitation to triangles by incrementing curVisitFlag
  // and comparing with the visit flag stored with the triangle

  curVisitFlag++;

  // Recursively search, (starting with the tri the point is in)
  // search for all tris whose circumcircle contain the point and place 
  // in the tri_list 

  if (rv == CUBIT_SUCCESS)
  {
    tri_list.append( tri_ptr );
    tri_visited( tri_ptr, CUBIT_TRUE );
    int iedge;
    TRI *adjtri_ptr;
    for (iedge=0; iedge<3 && rv == CUBIT_SUCCESS; iedge++)
	{
      int ii = iedge;
      adjtri_ptr = tri_ptr->adjacent( ii, refFacePtr );
      if (adjtri_ptr != NULL && !tri_visited( adjtri_ptr )){
		rv = point_in_circumcircle( the_point, adjtri_ptr, tri_list );
      }
    }
  }
  return rv;
}


//-------------------------------------------------------------------------
// Function:    locate_point
// Description: return the triangle the point is located in
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::locate_point( CubitVector &the_point,
				       TRI *&tri_ptr )
{
  CubitStatus rv = CUBIT_SUCCESS;

  // start with the last one found

  if (lastTriLocated != NULL)
    tri_ptr = lastTriLocated;

  // otherwise use the first one on the list

  else
  {
    tri_ptr = facetList->get();
  }


  // loop until we find something

  NODE *n0, *n1, *n2;
  double aa, bb, cc;
  CubitBoolean found = CUBIT_FALSE;

  //avoiding infinite loop
  int counter = 0;
  int max_count = facetList->size();
  
  while(!found && rv == CUBIT_SUCCESS)
  {
    tri_ptr->tri_nodes( n0, n1, n2 );
    aa = DETERM(the_point.x(), the_point.y(),
                n1->coordinates().x(), n1->coordinates().y(),
                n2->coordinates().x(), n2->coordinates().y());
    bb = DETERM(n0->coordinates().x(), n0->coordinates().y(),
                the_point.x(), the_point.y(),
                n2->coordinates().x(), n2->coordinates().y());
    cc = DETERM(n0->coordinates().x(), n0->coordinates().y(),
                n1->coordinates().x(), n1->coordinates().y(),
                the_point.x(), the_point.y());
    if (aa > -FT_INSIDE_TOL &&
        bb > -FT_INSIDE_TOL &&
        cc > -FT_INSIDE_TOL)
    {
      found = CUBIT_TRUE;  // this is the one
    }
    else
    {
      // set up to check the next logical neighbor
      if (aa <= bb && aa <= cc) 
	  {
		int edge_index = 1;
		tri_ptr = tri_ptr->adjacent( edge_index, refFacePtr );
      }
      else if (bb <= aa && bb <= cc) 
	  {
		int edge_index = 2;
		tri_ptr = tri_ptr->adjacent( edge_index, refFacePtr );
      }
      else 
	  {
		int edge_index = 0;
		tri_ptr = tri_ptr->adjacent( edge_index, refFacePtr );
      }
      // check to see if we've left the triangulation
      // also make sure that we are not stuck in a cycle
      if (tri_ptr == NULL || counter > max_count)
      {
        if(counter>max_count){
          PRINT_WARNING("Encountered problem locating a triangle; going to exhaustive search.\n");
        }
        rv = exhaustive_locate_point(the_point, tri_ptr );
        found = CUBIT_TRUE;
      }
    }
    ++counter;
  }
  
  lastTriLocated = tri_ptr;

  return rv;
}

//-------------------------------------------------------------------------
// Function:    exhaustive_locate_point
// Description: return the triangle the point is located in by checking 
//              all triangles
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::exhaustive_locate_point(
  CubitVector &the_point,
  TRI *&tri_ptr )
{
  CubitStatus rv = CUBIT_SUCCESS;

  
  // loop until we find something

  int ii;
  NODE *n0, *n1, *n2;
  double aa, bb, cc;
  CubitBoolean found = CUBIT_FALSE;
  for (ii=0; ii<facetList->size() && !found; ii++)
  {
    tri_ptr = facetList->get_and_step();
    tri_ptr->tri_nodes( n0, n1, n2 );
    aa = DETERM(the_point.x(), the_point.y(),
                n1->coordinates().x(), n1->coordinates().y(),
                n2->coordinates().x(), n2->coordinates().y());
    bb = DETERM(n0->coordinates().x(), n0->coordinates().y(),
                the_point.x(), the_point.y(),
                n2->coordinates().x(), n2->coordinates().y());
    cc = DETERM(n0->coordinates().x(), n0->coordinates().y(),
                n1->coordinates().x(), n1->coordinates().y(),
                the_point.x(), the_point.y());
    if (aa > -FT_INSIDE_TOL &&
        bb > -FT_INSIDE_TOL &&
        cc > -FT_INSIDE_TOL)
    {
      found = CUBIT_TRUE;  // this is the one
    }
  }
  if (!found)
  {
    rv = CUBIT_FAILURE;
    tri_ptr = NULL;
  }

  return rv;
}

//-------------------------------------------------------------------------
// Function:    are_nodes_colinear
// Description: determine if the TRI is valid
// Author:      mbrewer
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitBoolean FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::are_nodes_colinear( TRI *tri_ptr )
{
  NODE *nodes[3];
  tri_ptr->tri_nodes( nodes[0], nodes[1], nodes[2] );
  double det = DETERM( nodes[0]->coordinates().x(),
                       nodes[0]->coordinates().y(),
                       nodes[1]->coordinates().x(),
                       nodes[1]->coordinates().y(),
                       nodes[2]->coordinates().x(),
                       nodes[2]->coordinates().y());
    //PRINT_INFO("Det = %f\n",det);
  
  if(fabs(det) > CUBIT_DBL_MIN){
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
                       
}

//-------------------------------------------------------------------------
// Function:    point_in_circumcircle
// Description: determine if the point is inside the circumcircle of the
//              triangle and recurse to the adjacent triangles
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::point_in_circumcircle(
  CubitVector &the_point,
  TRI *tri_ptr,
  DLIList <TRI *> &tri_list)
{
  CubitStatus rv = CUBIT_SUCCESS;

  // check this triangle.  If the nodes are colinear do not try to calculate
    //the circumcenter.
  if(!are_nodes_colinear(tri_ptr))
  {
    PRINT_ERROR("Can't evaluate center of circumcircle\n");
    return CUBIT_FAILURE;
  }
     
  CubitVector cc = circumcenter( tri_ptr );
  tri_visited( tri_ptr, CUBIT_TRUE );
  double dist2 = SQR(the_point.x() - cc.x()) + SQR(the_point.y() - cc.y());
  double r2 = radius( tri_ptr );
  double tol_factor = CUBIT_MAX(CUBIT_MAX(tri_ptr->edge(0)->length(),
                                          tri_ptr->edge(1)->length()),
                                tri_ptr->edge(2)->length());
    //PRINT_INFO("Tolerance factor = %f\n", tol_factor);
  
  
  
  if (r2-dist2 > -(tol_factor*FT_INSIDE_TOL*FT_INSIDE_TOL))// inside or on circle
  {
    tri_list.append( tri_ptr );

    // go to its neighbors

    int iedge;
    TRI *adjtri_ptr;
    for (iedge=0; iedge<3 && rv == CUBIT_SUCCESS; iedge++){
     
      int ii = iedge;
      adjtri_ptr = tri_ptr->adjacent( ii, refFacePtr );
      if (adjtri_ptr != NULL && !tri_visited( adjtri_ptr ))
      {
        rv = point_in_circumcircle( the_point, adjtri_ptr, tri_list );
      }
    }
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
  classify_triangles(tri_list);

  // process each of the triangles until done

  TRI *tri_ptr;
  while ((tri_ptr = next_triangle()) != NULL && rv == CUBIT_SUCCESS) 
  {
    rv = insert_at_circumcenter( tri_ptr );
  }

  delete [] triSortArray;
  triSortArray = NULL;

  return rv;
}

//-------------------------------------------------------------------------
// Function:    classify_triangles
// Description: order the triangles in the current mesh by their worst angle
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::classify_triangles(DLIList<TRI *> &tri_list)
{
  CubitStatus rv = CUBIT_SUCCESS;

  // Create an array of lists that will hold triangles as they are
  // created (and deleted).  Each list will hold triangles whose angles 
  // fall within a specified interval.  Triangles will be processed
  // from smallest angle to best angle

  triSortArray = new DLIList<TRI *> [NUM_SORT_LISTS];

  // classify each trriangle and place into sort lists

  TRI *tri_ptr;
  int ii;
  for(ii=0; ii<tri_list.size() && rv == CUBIT_SUCCESS; ii++)
  {
    tri_ptr = tri_list.get_and_step();
    //rv = 
    classify_tri_by_angle( tri_ptr );
  }
   
  return rv;
}

//-------------------------------------------------------------------------
// Function:    classify_tri_by_angle
// Description: compute the angles at the triangle vertices and classify
//              by its worst triangle
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::classify_tri_by_angle( TRI *tri_ptr )
{
  //CubitStatus rv = CUBIT_SUCCESS;

  // Determine the minimum angle

  NODE *nodes[3];
  tri_ptr->tri_nodes( nodes[0], nodes[1], nodes[2] );
  double x0 = nodes[0]->coordinates().x();
  double y0 = nodes[0]->coordinates().y();
  double x1 = nodes[1]->coordinates().x();
  double y1 = nodes[1]->coordinates().y();
  double x2 = nodes[2]->coordinates().x();
  double y2 = nodes[2]->coordinates().y();

  double ax = x1 - x0;
  double ay = y1 - y0;
  double bx = x2 - x0;
  double by = y2 - y0;
  double dot = ax*bx + ay*by;
  double a_mag = sqrt( ax*ax + ay*ay );
  double b_mag = sqrt( bx*bx + by*by );
  double angle0 = dot / ( a_mag * b_mag );
  angle0 = acos( angle0 );

  ax = -ax;
  ay = -ay;
  bx = x2 - x1;
  by = y2 - y1;
  dot = ax*bx + ay*by;
  b_mag = sqrt( bx*bx + by*by );
  double angle1 = dot / ( a_mag * b_mag );
  angle1 = acos( angle1 );

  double angle2 = CUBIT_PI - angle0 - angle1;

  double min_angle = CUBIT_MIN( CUBIT_MIN(angle0,angle1), 
                                CUBIT_MIN(angle1,angle2) );
  if (min_angle < 0.0) {
    assert(0);
    return CUBIT_FAILURE;
  }

  // If it is greater than the QUALITY_ANGLE then place it in
  // the triSortArray[0], otherwise place it in one of the other deques
  // depending upon its minimum angle

  // Determine which list

  int index;
  if (min_angle >= QUALITY_ANGLE) {
    index = 0;
  }
  else {
    index = 1 + (int)(min_angle/INTERVAL);
    if (index < 1) index = 1;
    if (index > NUM_SORT_LISTS-1) index = NUM_SORT_LISTS-1;
  }

  // Place it on a list

  triSortArray[index].append( tri_ptr );
  tri_sort_list( tri_ptr, index );

  return CUBIT_SUCCESS;

}


//-------------------------------------------------------------------------
// Function:    next_triangle
// Description: get the next triangle to process
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
TRI *FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::next_triangle(void)
{
  TRI *tri_ptr = NULL;
  int ii;
  for(ii=1; ii<NUM_SORT_LISTS && tri_ptr == NULL; ii++)
  {
    if (triSortArray[ii].size() > 0)
      tri_ptr = triSortArray[ii].remove();
  }
  if (tri_ptr == NULL)
  {
    if (triSortArray[0].size() > 0)
      tri_ptr = triSortArray[0].remove();
  }
  return tri_ptr;
}
//-------------------------------------------------------------------------
// Function:    get_size
// Description: get the distortion factor for point inside tri, if one exists
//				otherwise return 1;
// Author:      chynes
// Date:        7/24/02
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
double FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::get_size(CubitVector &cc, TRI *tri_ptr)
{
	//extract data
	NODE *n0,*n1,*n2;
	CubitVector area_coord;
	tri_ptr->tri_nodes(n0,n1,n2);
	
	if    (n0->coordinates().z() - 1.0 < fabs(FT_INSIDE_TOL) 
		&& n1->coordinates().z() - 1.0 < fabs(FT_INSIDE_TOL)
		&& n2->coordinates().z() - 1.0 < fabs(FT_INSIDE_TOL)  )
		return 1.0;
	else 
	{
		//get vectors
		CubitVector v0 = n0->coordinates();
		CubitVector v1 = n1->coordinates();
		CubitVector v2 = n2->coordinates();

		//set z direction
		v0.z(cc.z());
		v1.z(cc.z());
		v2.z(cc.z());

		//create points
		CubitPoint *p0 = (CubitPoint*) new CubitPointData(v0);
		CubitPoint *p1 = (CubitPoint*) new CubitPointData(v1);
		CubitPoint *p2 = (CubitPoint*) new CubitPointData(v2);

		//create facet
		CubitFacet *temp_facet = (CubitFacet*) new CubitFacetData(p0,p1,p2);

		FacetEvalTool::facet_area_coordinate(temp_facet, cc, area_coord);
		
		delete p0;
		delete p1;
		delete p2;
		delete temp_facet;

		return	  (area_coord.x()*n0->coordinates().z()) 
				+ (area_coord.y()*n1->coordinates().z()) 
				+ (area_coord.z()*n2->coordinates().z());
	}

}

//-------------------------------------------------------------------------
// Function:    insert_at_circumcenter
// Description: insert a node at the circumcenter of a tri
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus FacetorTool<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::insert_at_circumcenter(TRI *tri_ptr)
{
  CubitStatus rv = CUBIT_SUCCESS;

  // find the cicumcenter of the triangle and the target size there
    //if nodes are colinear do not try to find circumenter
  if(!are_nodes_colinear(tri_ptr))
  {
    PRINT_ERROR("Can't evaluate center of circumcircle\n");
    return CUBIT_FAILURE;
  }
  CubitVector cc = circumcenter( tri_ptr );
  
  //translate cc into 3D space
  CubitVector cc_xyz;
  pTool->transform_to_xyz(cc_xyz, cc);
  // get target length in 3D space
  double target_length = sizingFunction->size_at_point( cc_xyz );
  // get new size
  double size = get_size(cc, tri_ptr);
  // update size
  cc.z(size);
  // update target_length
  target_length = target_length*size;

  // Determine if we should now insert the point.  Allow insertions
  // in the general case down to circumcircle size of ALPHA times the
  // interpolated target edge length size.  For tris with small
  // angles, allow additional inserts to improve the quality down
  // to 1/2 the target size
  if(!are_nodes_colinear(tri_ptr))
  {
    PRINT_ERROR("Can't evaluate radius of circumcircle\n");
    return CUBIT_FAILURE;
  }
     
    
  double r2 = radius( tri_ptr );
  CubitBoolean insert = CUBIT_FALSE;
  int tsindex = tri_sort_list( tri_ptr );
  assert(tsindex > -1);
  if (tsindex == 0) 
  {
    if (r2 > SQR(ALPHA*target_length)) 
    {
      insert = CUBIT_TRUE;
    }
  }
  else 
  {
    if (r2 > SQR(0.5*ALPHA*target_length)) 
    {
      insert = CUBIT_TRUE;
    }
  }
  if (insert) 
  {

    // Determine the tris that will be affected by the insertion

    lastTriLocated = tri_ptr;
    DLIList <TRI *> tri_list;
    rv = natural_neighbor_tris( cc, tri_list );
    // If it was outside, then we are done with it

    if (tri_list.size() == 0) 
    {
      return CUBIT_SUCCESS;
    }
	if (rv != CUBIT_SUCCESS) {
      return rv;
    }

    // See if we are too close to a boundary

    double x0, y0, x1, y1, cx, cy, edge_radius, dist;
    TRI *nntri_ptr;
    EDGE *edge_ptr;
    int ii, iedge;
    for (ii=0; ii<tri_list.size(); ii++) 
    {
      nntri_ptr = tri_list.get_and_step();
      for (iedge=0; iedge<3; iedge++) 
      {
		edge_ptr = tri_ptr->edge( iedge );

        // An edge encroaches if the distance from the prospective
        // new point to the midpoint of the edge is less than 
        // half the length of the edge

        if (edge_ptr->marked()) // on the boundary?
        {
          x0 = (edge_ptr->start_node())->coordinates().x();
          y0 = (edge_ptr->start_node())->coordinates().y();
          x1 = (edge_ptr->end_node())->coordinates().x();
          y1 = (edge_ptr->end_node())->coordinates().y();
          cx = (x0 + x1) * 0.5;
          cy = (y0 + y1) * 0.5;
          edge_radius = sqrt(SQR(x1-x0) + SQR(y1-y0)) * 0.5;       
          dist = sqrt( SQR(cx-cc.x()) + SQR(cy-cc.y()) );

          // Edge encroaches: don't insert, return now

          if (dist - edge_radius < FT_INSIDE_TOL) 
          {
            return CUBIT_SUCCESS;
          }
        }
      }
    }

    // Before inserting, remove all the tris on the neighbor
    // tri_list from the lists 
   
    int index;
    for (ii=0; ii<tri_list.size(); ii++) 
    {
      nntri_ptr = tri_list.get_and_step();
      index = tri_sort_list( nntri_ptr );
      assert(index >= 0);
      triSortArray[index].remove( nntri_ptr );
    }

    // Create the new node

    NODE *new_node_ptr = (NODE *)new NODECHILD( cc, refFacePtr );

    // Insert the new node into the mesh

    rv = bowyer_watson_insert( new_node_ptr, tri_list );
    if (rv != CUBIT_SUCCESS)
      return rv;

    // get the new tris at the node and classify them

    tri_list.clean_out();
    new_node_ptr->tris( tri_list );
    for (ii=0; ii<tri_list.size() && rv == CUBIT_SUCCESS; ii++) 
    {
      tri_ptr = tri_list.get_and_step();
      rv = classify_tri_by_angle( tri_ptr );
    }
  }

  return rv;

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
 
