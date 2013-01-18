
#ifdef INLINE_TEMPLATES
#define MY_INLINE inline
#else
#define MY_INLINE
#endif

#include "FacetorUtil.hpp"
#include "TDDelaunay.hpp"
#include "FacetEvalTool.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "ParamTool.hpp"

#define DETERM(p1,q1,p2,q2,p3,q3)\
     ((q3)*((p2)-(p1)) + (q2)*((p1)-(p3)) + (q1)*((p3)-(p2)))
#define FT_INSIDE_TOL 1.0e-6
#define SQR(x) ((x) * (x))
#define ALPHA 0.70228615

template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::locate_point(
                                      CubitVector& the_point,
                                      DLIList<TRI*>& facet_list,
                                      TRI* starting_tri,
                                      SURF* owning_surface,
                                      TRI*& tri_ptr)
{
  CubitStatus rv = CUBIT_SUCCESS;

    // start with the last one found

  if (starting_tri != NULL)
      tri_ptr = starting_tri;

    // otherwise use the first one on the list

  else
  {
    tri_ptr = facet_list.get();
  }


    // loop until we find something

  NODE *n0, *n1, *n2;
  double aa, bb, cc;
  CubitBoolean found = CUBIT_FALSE;

    //avoiding infinite loop
  int counter = 0;
  int max_count = facet_list.size();
  
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
        tri_ptr = tri_ptr->adjacent( edge_index, owning_surface );
      }
      else if (bb <= aa && bb <= cc) 
      {
        int edge_index = 2;
        tri_ptr = tri_ptr->adjacent( edge_index, owning_surface );
      }
      else 
      {
        int edge_index = 0;
        tri_ptr = tri_ptr->adjacent( edge_index, owning_surface );
      }
        // check to see if we've left the triangulation
        // also make sure that we are not stuck in a cycle
      if (tri_ptr == NULL || counter > max_count)
      {
        if(counter>max_count){
          PRINT_WARNING("Encountered problem locating a triangle; going to exhaustive search.\n");
        }

        rv = exhaustive_locate_point( the_point, facet_list, tri_ptr );
        found = CUBIT_TRUE;
      }
    }
    ++counter;
  }
  
  return rv;
}

template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::exhaustive_locate_point(
                                                 CubitVector& the_point,
                                                 DLIList<TRI*>& facet_list,
                                                 TRI*& tri_ptr)
{
  CubitStatus rv = CUBIT_SUCCESS;
  
  // loop until we find something

  int ii;
  NODE *n0, *n1, *n2;
  double aa, bb, cc;
  CubitBoolean found = CUBIT_FALSE;
  for (ii=0; ii<facet_list.size() && !found; ii++)
  {
    tri_ptr = facet_list.get_and_step();
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
CubitBoolean
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::are_nodes_colinear(TRI* tri_ptr)
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
// Function:    circumcenter
// Description: get the circumcenter of the triangle
// Author:      chynes
// Date:        6/6/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitVector& 
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::circumcenter(TRI* tri_ptr)
{
	ToolData *td = tri_ptr->get_TD( TDDelaunay<TRI,NODE>::is_delaunay );
	TDDelaunay< TRI, NODE > *td_del = dynamic_cast<TDDelaunay< TRI, NODE >*> (td);
	if(td_del == NULL) {
		td_del = new TDDelaunay<TRI, NODE>();
		tri_ptr->add_TD( td_del );
	}
	return td_del->circumcenter2d( tri_ptr );
}

//-------------------------------------------------------------------------
// Function:    tri_visited
// Description: set the tri_visited flag
// Author:      chynes
// Date:        6/6/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitBoolean
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::tri_visited(
                                      TRI *tri_ptr,
                                      int curr_visit_flag)
{
  ToolData *td = tri_ptr->get_TD( TDDelaunay< TRI, NODE >::is_delaunay );
  TDDelaunay< TRI, NODE > *td_del = dynamic_cast<TDDelaunay< TRI, NODE >*> (td);
  if (td_del == NULL)
  {
    td_del = new TDDelaunay< TRI, NODE >();
    tri_ptr->add_TD( td_del );
  }
  return (td_del->visit_flag() == curr_visit_flag);
}

//-------------------------------------------------------------------------
// Function:    tri_visited
// Description: set the tri_visited flag
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
void
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::tri_visited(
                               TRI *facet_ptr,
                               CubitBoolean visited,
                               int curr_visit_flag)
{
  ToolData *td = facet_ptr->get_TD( TDDelaunay< TRI, NODE >::is_delaunay );
  TDDelaunay< TRI, NODE > *td_del = dynamic_cast<TDDelaunay< TRI, NODE >*> (td);
  if (td_del == NULL)
  {
    td_del = new TDDelaunay< TRI, NODE >();
    facet_ptr->add_TD( td_del );
  }
  if (visited)
    td_del->visit_flag(curr_visit_flag);
  else
    td_del->visit_flag(INT_MIN);
}

//-------------------------------------------------------------------------
// Function:    radius
// Description: get the radius squared of the triangle circumcircle
// Author:      chynes
// Date:        6/6/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
double
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::radius(TRI* tri_ptr)
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
// Function:    point_in_circumcircle
// Description: determine if the point is inside the circumcircle of the
//              triangle and recurse to the adjacent triangles
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::point_in_circumcircle(
                                               CubitVector& the_point,
                                               TRI* tri_ptr,
                                               DLIList<TRI*>& tri_list,
                                               SURF* surface_ptr,
                                               int curr_visit_flag)
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
  tri_visited( tri_ptr, CUBIT_TRUE, curr_visit_flag );
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
      adjtri_ptr = tri_ptr->adjacent( ii, surface_ptr );
      if (adjtri_ptr != NULL && !tri_visited( adjtri_ptr, curr_visit_flag ))
      {
        rv = point_in_circumcircle( the_point, adjtri_ptr, tri_list,
                                                   surface_ptr, curr_visit_flag );
      }
    }
  }
  return rv;
}

//-------------------------------------------------------------------------
// Function:    natural_neighbor_tris
// Description: get a list of all triangles whose circumcircle contain 
//              the point
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::natural_neighbor_tris(
                                               CubitVector& the_point,
                                               DLIList<TRI*>& facet_list,
                                               TRI*& start_tri,
                                               SURF* surface_ptr,
                                               int& curr_visit_flag,
                                               DLIList<TRI*>& tri_list)
{
  CubitStatus rv = CUBIT_SUCCESS;

  // find the triangle the point is contained in

  //CubitVector areacoord;
  TRI *tri_ptr;
  rv = locate_point( the_point, facet_list, start_tri, surface_ptr, tri_ptr );
  start_tri = tri_ptr;
  
  // keep track of visitation to triangles by incrementing curr_visit_flag
  // and comparing with the visit flag stored with the triangle

  curr_visit_flag++;

  // Recursively search, (starting with the tri the point is in)
  // search for all tris whose circumcircle contain the point and place 
  // in the tri_list 

  if (rv == CUBIT_SUCCESS)
  {
    tri_list.append( tri_ptr );
    tri_visited( tri_ptr, CUBIT_TRUE, curr_visit_flag );
    int iedge;
    TRI *adjtri_ptr;
    for (iedge=0; iedge<3 && rv == CUBIT_SUCCESS; iedge++)
	{
      int ii = iedge;
      adjtri_ptr = tri_ptr->adjacent( ii, surface_ptr );
      if (adjtri_ptr != NULL && !tri_visited( adjtri_ptr, curr_visit_flag )){
        rv = point_in_circumcircle( the_point, adjtri_ptr,
                                                   tri_list, surface_ptr,
                                                   curr_visit_flag);
      }
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
CubitStatus
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::valid_void(
                                     NODE * /*point_ptr*/,
                                     DLIList<TRI *> &tri_list,
                                     SURF* surface_ptr,
                                     int curr_visit_flag)
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
      adjtri_ptr = tri_ptr->adjacent( kk, surface_ptr );
      if (!adjtri_ptr || !tri_visited( adjtri_ptr, curr_visit_flag ))
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
// Function:    bowyer_watson_insert
// Description: Bowyer-Watson insertion into an existing Delaunay Mesh
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::bowyer_watson_insert(
                                              NODE* point_ptr,
                                              DLIList<TRI*>& tri_list,
                                              DLIList<TRI*>& facet_list,
                                              int& curr_visit_flag,
                                              SURF* surface_ptr,
                                              TRI*& last_tri)
{
  CubitStatus rv = CUBIT_SUCCESS;

    // mark the tris in the list so we can distinguish them from their 
    // neighbors

  curr_visit_flag++;
  int ii, jj;
  TRI *tri_ptr;
  for (ii=0; ii<tri_list.size(); ii++)
  {
    tri_ptr = tri_list.get_and_step();
    tri_visited( tri_ptr, CUBIT_TRUE, curr_visit_flag );
  }
    //MBREWER:  This is not an optimal test.  But, we need need to 
    //do some tests to try to ensure that the void is valid for what
    //we need.  This is attempting to avoid crashes by not allowing nodes
    //to be inserted when the mesh starts diverging from the Delaunay
    //criteria.
  rv = valid_void( point_ptr, tri_list, surface_ptr, curr_visit_flag );
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
      adjtri_ptr = tri_ptr->adjacent( kk, surface_ptr );
      if (!adjtri_ptr || !tri_visited( adjtri_ptr, curr_visit_flag ))
      {
        edge_ptr = tri_ptr->edge(kk);
        assert(edge_list.append_unique( edge_ptr ));
        if(tri_ptr->sense(kk) == CUBIT_FORWARD)
            new_tri = (TRI *) new TRICHILD( edge_ptr->start_node(), edge_ptr->end_node(), point_ptr, surface_ptr);
        else
            new_tri = (TRI *) new TRICHILD( edge_ptr->end_node(), edge_ptr->start_node(), point_ptr, surface_ptr);
        facet_list.append(new_tri);
      }
    }
  }
  last_tri = new_tri;

    // delete the triangles in the original triangle list

  EDGE *del_edge_ptr;
  DLIList<EDGE *> del_edge_list;
  for (ii=0; ii<tri_list.size(); ii++)
  {
    tri_ptr = tri_list.get_and_step();
    del_edge_list.clean_out();
    facet_list.move_to(tri_ptr);
    facet_list.extract();
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

//-------------------------------------------------------------------------
// Function:    insert_node
// Description: insert one node into Delaunay mesh
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::insert_node(
                                     NODE *node_ptr,
                                     DLIList<TRI*>& facet_list,
                                     SURF* surface_ptr,
                                     int& curr_visit_flag,
                                     TRI*& last_tri)
{
  CubitStatus rv = CUBIT_SUCCESS;

  // get a list of all triangles whose circumcircle contain the point
  
  DLIList<TRI *> tri_list;
  CubitVector the_point = node_ptr->coordinates();
  rv = natural_neighbor_tris( the_point, facet_list,
                              last_tri, surface_ptr,
                              curr_visit_flag, tri_list );


  // Use a Bowyer-Watson insertion 
  
  if (rv == CUBIT_SUCCESS)
  {
    rv = bowyer_watson_insert( node_ptr, tri_list,
                               facet_list, curr_visit_flag,
                               surface_ptr, last_tri);
  }

  return rv;
}

//-------------------------------------------------------------------------
// Function:    get_size
// Description: get the distortion factor for point inside tri, if one exists
//              otherwise return 1;
// Author:      chynes
// Date:        7/24/02
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
double
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::get_size(CubitVector &cc, TRI *tri_ptr)
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
// Function:    tri_sort_list
// Description: set the tri sort list index
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
void
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::tri_sort_list(
                                TRI *facet_ptr,
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
int
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::tri_sort_list( TRI *facet_ptr )
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
// Function:    classify_tri_by_angle
// Description: compute the angles at the triangle vertices and classify
//              by its worst triangle
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::classify_tri_by_angle(
                                               TRI* tri_ptr,
                                               DLIList<TRI*>* sorted_lists,
                                               const int num_lists,
                                               const double interval,
                                               const double quality_angle)
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
  if (min_angle >= quality_angle) {
    index = 0;
  }
  else {
    index = 1 + (int)(min_angle/interval);
    if (index < 1) index = 1;
    if (index > num_lists-1) index = num_lists-1;
  }

  // Place it on a list

  sorted_lists[index].append( tri_ptr );
  tri_sort_list( tri_ptr, index );

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Function:    insert_at_circumcenter
// Description: insert a node at the circumcenter of a tri
// Author:      chynes
// Date:        6/3/2002
//-------------------------------------------------------------------------
template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC> MY_INLINE
CubitStatus
FacetorUtil<SURF, TRI, EDGE, NODE, TRICHILD, NODECHILD, SIZEFUNC>::insert_at_circumcenter(
                                                TRI* tri_ptr,
                                                DLIList<TRI*>& facet_list,
                                                TRI*& start_tri,
                                                int& curr_visit_flag,
                                                SURF* surface_ptr,
                                                DLIList<TRI*>* sorted_lists,
                                                const int num_lists,
                                                const double interval,
                                                const double quality_angle,
                                                SIZEFUNC* sizing_function,
                                                ParamTool* p_tool)
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
  p_tool->transform_to_xyz(cc_xyz, cc);
  // get target length in 3D space
  double target_length = sizing_function->size_at_point( cc_xyz );
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

    start_tri = tri_ptr;
    DLIList <TRI *> tri_list;
    rv = natural_neighbor_tris( cc, facet_list,
                                               start_tri, surface_ptr,
                                               curr_visit_flag, tri_list );
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
      sorted_lists[index].remove( nntri_ptr );
    }

      // Create the new node

    NODE *new_node_ptr = (NODE *)new NODECHILD( cc, surface_ptr );

      // Insert the new node into the mesh

    rv = bowyer_watson_insert( new_node_ptr, tri_list,
                               facet_list, curr_visit_flag,
                               surface_ptr, start_tri);
    if (rv != CUBIT_SUCCESS)
        return rv;

      // get the new tris at the node and classify them

    tri_list.clean_out();
    new_node_ptr->tris( tri_list );
    for (ii=0; ii<tri_list.size() && rv == CUBIT_SUCCESS; ii++) 
    {
      tri_ptr = tri_list.get_and_step();
      rv = classify_tri_by_angle( tri_ptr, sorted_lists, num_lists, interval, quality_angle );
    }
  }

  return rv;
}
