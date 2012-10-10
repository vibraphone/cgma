//- Class: FacetEvalTool
//- Description:  The FacetEvalTool is a tool to perform generic geometric 
//-               operations on a set of facets.
//- Assumptions:  All of the facets are connected topologically correct.
//-
//- Owner: Steve J. Owen
//- Checked by: 

#include <float.h>
#include "CastTo.hpp"
#include "CubitVector.hpp"
#include "CubitPoint.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CpuTimer.hpp"
#include "CubitMessage.hpp"
#include "CubitBox.hpp"
#include "DLIList.hpp"
#include "FacetEvalTool.hpp"
#include "GeometryDefines.h"
#include "CubitTransformMatrix.hpp"
#include "debug.hpp"
#include "TDFacetBoundaryEdge.hpp"
#include "TDFacetBoundaryPoint.hpp"
#include "GfxDebug.hpp"
//#include "FacetParamTool.hpp"
#include "KDDTree.hpp"
#include "AbstractTree.hpp"
#include "CubitFileIOWrapper.hpp"
#include "FacetDataUtil.hpp"

double FacetEvalTool::timeGridSearch = 0.0;
double FacetEvalTool::timeFacetProject = 0.0;
int FacetEvalTool::numEvals = 0;
#define GRID_SEARCH_THRESHOLD 20

//===========================================================================
//Function Name: FacetEvalTool
//
//Member Type:  PUBLIC
//Description:  constructor
//===========================================================================
FacetEvalTool::FacetEvalTool( 
  DLIList<CubitFacet*> &facet_list,
  DLIList<CubitPoint*> &point_list,  // if this isn't filled in, I'll do it here
  int interp_order, // 0 = linear or 4 = Bezier only for now
  double min_dot )  // from feature angle
{
  static int counter = 1;
  toolID = counter++;
  myFacetList = facet_list;
  myPointList = point_list;
  isParameterized = CUBIT_FALSE;
  minDot = min_dot;
  output_id = -1;

  // mark all of the facets on this surface with the toolID

  int i;
  for (i = 0; i< myFacetList.size(); i++)
    myFacetList.get_and_step()->set_tool_id( toolID );

  // generate the list of points, if not already defined

  if (myPointList.size() == 0)
  {
    get_points_from_facets( facet_list, myPointList );
  }
  interpOrder = 0;  // default to project to linear facet

  // set the bounding box and compareTol and setup grid search
  myBBox = NULL;
  aTree = NULL;
  
  reset_bounding_box();
  
             
  lastFacet = NULL;
  isFlat = -999;
  isFlat = is_flat();

  if (interp_order == -1) {
    interpOrder = 0;  // use linear as default interp order
  }
  else {
    interpOrder = interp_order;
  }

  // generate edges 

  get_edges_from_facets( facet_list, myEdgeList );

  // generate loops
  
  get_loops_from_facets( myEdgeList, myLoopList );

  if (interpOrder > 0) {
    init_gradient();
    if (interpOrder == 3) { // least_squares
      init_quadrics();
    }
    else if(interpOrder == 4) { // spline
      init_bezier_surface();
    }
  }

  // compute the area of all facets

  myArea = -1.0;
  myArea = area();

  
  int mydebug = 0;
  if (mydebug)
  {
    debug_draw_eval_bezier_facet( myFacetList.get_and_step() );
    debug_draw_facets();
    debug_draw_facet_normals();
    if (interpOrder > 0) debug_draw_point_normals();
    GfxDebug::flush();
    GfxDebug::mouse_xforms();
  }

  //parameterize();
}

//===========================================================================
//Function Name: FacetEvalTool
//Member Type:  PUBLIC
//Descriptoin:  constructor.  This one is an empty constructor
//              sed with restore method
//===========================================================================
FacetEvalTool::FacetEvalTool()
{
  static int counter = 1;
  toolID = counter++;
  myBBox = NULL;
  aTree = NULL;
  lastFacet = NULL;
  isFlat = -999;
  myArea = -1.0;
  interpOrder = 0;
  minDot = 0;
}

//===========================================================================
//Function Name: ~FacetEvalTool
//
//Member Type:  PUBLIC
//Descriptoin:  destructor
//===========================================================================
FacetEvalTool::~FacetEvalTool()
{
  if ( DEBUG_FLAG(110) )
  {
    PRINT_INFO("num facets       = %d\n",myFacetList.size());
    PRINT_INFO("numEvals         = %d\n",numEvals);
    PRINT_INFO("timeFacetProject = %f\n",timeFacetProject);
    PRINT_INFO("timeGridSearch   = %f\n",timeGridSearch);
  }

  if (aTree != NULL)
    delete aTree;

  if (myBBox) delete myBBox;
  destroy_facets();

  myLoopList.reset();
  int i;
  for (i=myLoopList.size(); i>0; i--)
  {
    delete myLoopList.get_and_step();
  }
  myLoopList.clean_out();
}

void FacetEvalTool::remove_facets( DLIList<CubitFacet*>& facets )
{
  facets = myFacetList;
  for ( int i = facets.size(); i--; )
    facets.step_and_get()->set_tool_id(0);
  myFacetList.clean_out();
  myEdgeList.clean_out();
  myPointList.clean_out();
}

CubitStatus FacetEvalTool::reverse_facets( )
{
  int i;
  CubitFacet* temp_facet;
  
  for(i=0; i<myFacetList.size(); i++){
    temp_facet=myFacetList.get_and_step();
    if(!temp_facet){
      PRINT_ERROR("Unexpected NULL pointer for facet.\n");
      return CUBIT_FAILURE;
    }
    temp_facet->flip();
  }

  return CUBIT_SUCCESS;
}


void FacetEvalTool::replace_facets( DLIList<CubitFacet *> &facet_list )
{
  myFacetList.clean_out();
  myEdgeList.clean_out();
  myPointList.clean_out();
  FacetDataUtil::get_points( facet_list, myPointList );
  FacetDataUtil::get_edges( facet_list, myEdgeList );
  myFacetList = facet_list;
  for ( int i = myFacetList.size(); i--; )
    myFacetList.step_and_get()->set_tool_id(toolID);
}
  

//===========================================================================
//Function Name: set_up_grid_search
//
//Member Type:  PUBLIC
//Descriptoin:  set up grid search if we have a lot of facets
//===========================================================================
void FacetEvalTool::set_up_grid_search(
  double geom_factor )
{
  if(aTree)
    delete aTree;
  
  if (myFacetList.size() < GRID_SEARCH_THRESHOLD)
  {
    aTree = NULL;
  }
  else
  {
    aTree = new KDDTree<CubitFacet*> (GEOMETRY_RESABS*geom_factor, false/*self-balancing off*/);
    for ( int ii = myFacetList.size(); ii > 0; ii--  )
    {
      aTree->add(myFacetList.get_and_step());
    }
    aTree->balance();
  }
}

//===========================================================================
//Function Name: is_flat
//
//Member Type:  PUBLIC
//Descriptoin:  Determine if the facets are all in the same plane
//===========================================================================
int FacetEvalTool::is_flat()
{
  int ii;
  if (isFlat != -999) {
    return isFlat;
  }
  else {
    isFlat = CUBIT_TRUE;
    CubitVector firstnormal = myFacetList.get_and_step()->normal();
    firstnormal.normalize();
    CubitVector normal;
    for ( ii = myFacetList.size(); ii > 1 &&isFlat; ii-- ){
      normal = myFacetList.get_and_step()->normal();
      normal.normalize();
      if (fabs(normal%firstnormal) < 0.9987654321) {
        isFlat = CUBIT_FALSE;
      }
    }
  }
  return isFlat;
}

//===========================================================================
//Function Name: init_gradient
//
//Member Type:  PRIVATE
//Descriptoin:  initialize the gradients for order 1 interpolation
//===========================================================================
CubitStatus FacetEvalTool::init_gradient()
{
  int i,j;

  // retrieve all faces attached to the points in point_list

  for (i = 0; i < myPointList.size(); i++)
  {
    CubitPoint* point = myPointList.get_and_step();

    DLIList<CubitFacet*> adj_facet_list;
    point->facets(adj_facet_list);
    if (adj_facet_list.size() > 0) {
      CubitVector avg_normal(0.0e0, 0.0e0, 0.0e0);
      double totangle = 0.0e0;

      // weight the normal by the spanning angle at the point

      for (j = 0; j < adj_facet_list.size(); j++)
      {
        CubitFacet* facet = adj_facet_list.get_and_step();
        double angle = facet->angle( point );
        facet->weight( angle );
        totangle += angle;
      }
      for (j = 0; j < adj_facet_list.size(); j++)
      {
        CubitFacet* facet = adj_facet_list.get_and_step();
        CubitVector normal = facet->normal();
        normal.normalize();
        avg_normal += (facet->weight() / totangle) * normal;
      }
      avg_normal.normalize();
      point->normal(avg_normal);
      double coefd = -(point->coordinates()%avg_normal);
      point->d_coef( coefd );
    }
  }
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: init_quadrics
// NOTE: I (Roshan) fixed couple of bugs on Aug 02, 2010.  See BUGFIX comments below.  I didn't test this function; however, I have tested CMLSmoothTool::init_quadric(). 
//Member Type:  PRIVATE
//Descriptoin:  initialize the quadrics at the facet vertices for order 2 
//              interpolation
//===========================================================================
CubitStatus FacetEvalTool::init_quadrics()
{
  // use the normal at the vertices as a local space with which to do 
  // interpolation.  A quadric will be approximated from the surrounding
  // vertices.  Interpolation will provide a "z" deviation from the tangent
  // plane at the vertex

  // define a basis set of vectors at each point (assumes the gradients
  // have already been approximated

  int i,j;
  CubitPoint* point;
  for (i = 0; i < myPointList.size(); i++)
  {
    point = myPointList.get_and_step();
    point->define_tangent_vectors(); 
  }

  // set up for least squares

  CubitStatus status;
#define MAX_CLOSE_POINTS 100
  CubitPoint *close_points[MAX_CLOSE_POINTS];
  CubitVector coords[MAX_CLOSE_POINTS], cp;
  double weight[MAX_CLOSE_POINTS];
  int num_close;
  CubitMatrix lhs(5,5);
  CubitMatrix rhs(5,1);
  CubitMatrix coef(5,1);
  for (i = 0; i < myPointList.size(); i++)
  {
    point = myPointList.get_and_step();
    status = get_close_points( point, close_points, num_close, 
                               MAX_CLOSE_POINTS, 5 );
    if (status != CUBIT_SUCCESS) {
      return status;
    }

    // transform to local system in x-y
    // determine weights based on inverse distance

    weight[0] = 0.0e0;
    double maxdist = -1e100;
    double totweight = 0.0e0;
    for(j=0; j<num_close; j++) {
      cp = close_points[j]->coordinates();
      point->transform_to_local( cp, coords[j] );
      weight[j] = sqrt( sqr(coords[j].x()) + sqr(coords[j].y()) );
      if (weight[j] > maxdist) maxdist = weight[j];
    }
    maxdist *= 1.1e0;
    for (j=0; j<num_close; j++) {
      weight[j] = sqr((maxdist-weight[j])/(maxdist*weight[j]));
      totweight += weight[j];
    }

    // fill up the matrices

    //lhs.set_to_identity();
    //rhs.set_to_identity();
    //coef.set_to_identity();

    double dx, dy, wjdx, wjdy, dx2, dy2, dxdy, dz;
    for (j=0; j<num_close; j++) {
      weight[j] /= totweight;
      weight[j] = 1; // BUGFIX: ignore weights for now and reset weights to 1
      dx = /*-*/ coords[j].x(); //BUGFIX: Why we need -ve coords?
      dy = /*-*/ coords[j].y(); //BUGFIX: Why we need -ve coords?
      wjdx = weight[j] * dx;
      wjdy = weight[j] * dy;
      dx2 = sqr( dx );
      dy2 = sqr( dy );
      dxdy = dx * dy;
      dz = coords[j].z(); 
      
      lhs.add( 0, 0, wjdx * dx );
      lhs.add( 0, 1, wjdx * dy );
      lhs.add( 0, 2, wjdx * dx2 );
      lhs.add( 0, 3, wjdx * dxdy );
      lhs.add( 0, 4, wjdx * dy2 );
      rhs.add( 0, 0, wjdx * dz ); // BUGFIX: dz was missing
      
      lhs.add( 1, 1, wjdy * dy );
      lhs.add( 1, 2, wjdy * dx2 );
      lhs.add( 1, 3, wjdy * dxdy );
      lhs.add( 1, 4, wjdy * dx * dy2 );
      rhs.add( 1, 0, wjdy * dz ); // BUGFIX: dz was missing
      
      lhs.add( 2, 2, wjdx * dx2 * dx );
      lhs.add( 2, 3, wjdx * dx2 * dy );
      lhs.add( 2, 4, wjdx * dx * dy2 ); 
      rhs.add( 2, 0, wjdx * dx * dz );// BUGFIX: dz was missing

      lhs.add( 3, 3, wjdx * dx * dy2 );
      lhs.add( 3, 4, wjdx * dy * dy2 );
      rhs.add( 3, 0, wjdx * dy * dz ); // BUGFIX: dz was missing
      
      lhs.add( 4, 4, wjdy * dy * dy2 );
      rhs.add( 4, 0, wjdy * dy * dz ); // BUGFIX: dz was missing
    }
    lhs.set( 1, 0, lhs.get(0,1) );
    lhs.set( 2, 0, lhs.get(0,2) );
    lhs.set( 2, 1, lhs.get(1,2) );
    lhs.set( 3, 0, lhs.get(0,3) );
    lhs.set( 3, 1, lhs.get(1,3) );
    lhs.set( 3, 2, lhs.get(2,3) );
    lhs.set( 4, 0, lhs.get(0,4) );
    lhs.set( 4, 1, lhs.get(1,4) );
    lhs.set( 4, 2, lhs.get(2,4) );
    lhs.set( 4, 3, lhs.get(3,4) );
    
    // solve the system

    status = lhs.solveNxN( rhs, coef );
    if (status != CUBIT_SUCCESS) {
      return status;
    }

    // store the quadric coefficents with the point

    point->coef_vector( coef );
  }
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: get_close_points
//
//Member Type:  PRIVATE
//Descriptoin:  get a list of points close to the current point on the facets
//              return a minimum of min_close and a maximum of max_close points
//===========================================================================
CubitStatus FacetEvalTool::get_close_points(
   CubitPoint *point,
   CubitPoint **close_point,
   int &num_close,
   int max_close,
   int min_close )
{

  // get the points immediately adjacent

  DLIList<CubitFacet*> adj_facet_list;
  point->facets(adj_facet_list);
  if (adj_facet_list.size() > max_close) {
    return CUBIT_FAILURE;
  }
  point->adjacent_points(close_point, num_close);
  
  // if we don't have enough yet, then go to the next level

  if (num_close < min_close) {
    CubitPoint *cpoint[100];
    CubitPoint *cur_point;
    CubitBoolean found;
    int nclose, cnclose;
    nclose = num_close;
    for(int i=0; i<num_close; i++) {
      cur_point = close_point[i];
      DLIList<CubitFacet*> cadj_facet_list;
      cur_point->facets(cadj_facet_list);
      if (cadj_facet_list.size() + nclose > max_close) {
        return CUBIT_FAILURE;
      }
      cur_point->adjacent_points(cpoint,cnclose);
      for(int k=0; k<cnclose; k++) {

        // check that it is not already on the list

        found = CUBIT_FALSE;
        for(int l=0; l<nclose && !found; l++) {
          if (close_point[l] == cpoint[k] ||
              cpoint[k] == point) {
            found = CUBIT_TRUE;
          }
        }
        if (!found) {

          // make sure the normal at this point is not more than 90 degrees 
          // from the normal at the point
          
          double dot = cpoint[k]->normal() % point->normal();
          if (dot > GEOMETRY_RESABS) {

            // add the point to the list

            close_point[nclose++] = cpoint[k];
          }
        }
      }
    }
    num_close = nclose;
    if (num_close < min_close) {
      return CUBIT_FAILURE;
    }
  }
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: mark edge pairs
//
//Member Type:  PRIVATE
//Descriptoin:  For the given point, we loop over the attached boundary
//              edges and figure out which pairs should be C1 continuous.
//              A given edge is paired with at most two other edges.  The
//              other edges (or the other edge) is stored on a tool data
//              to be used later.
//===========================================================================
CubitStatus FacetEvalTool::mark_edge_pairs(CubitPoint* point)
{
  int i,j;
  DLIList<CubitFacetEdge*> edge_list;
  DLIList<CubitFacetEdge*> edge_list_init;
  CubitPoint* prev_point = NULL;
  CubitPoint* next_point = NULL;
  CubitFacetEdge* prev_edge = NULL;
  CubitFacetEdge* next_edge = NULL;
  CubitVector e0, e1;
  double current_dot;
  point->edges(edge_list_init);
  TDFacetBoundaryEdge *td_fbe = NULL;
  
    // make a list with only the boundary edges
  for(i=0; i< edge_list_init.size(); i++){
    prev_edge = edge_list_init.get_and_step();
    td_fbe = TDFacetBoundaryEdge::get_facet_boundary_edge(prev_edge);
    if(td_fbe)
      edge_list.append(prev_edge);
  }
  prev_edge=NULL;
    //if there is one or none, we don't need to bother looking
    // for pairs.
  if(edge_list.size() < 2)
    return CUBIT_SUCCESS;
  int num_edges=edge_list.size();
  int* other_index = new int[num_edges];
  double* dot_array = new double[num_edges];  
  bool pair_exists = false;
    // loop over the edges.  For each edge, find the edge that would
    // make the best other edge to pair this one with
  for(i = 0; i<num_edges; i++){
    other_index[i]=-1;
    dot_array[i]=-1.0;
    
    prev_edge=edge_list[i];
    prev_point = prev_edge->other_point(point);
      // now figure out which other edge would be the best pair with
      // this one.  This could be sped up by only looking at the edges
      // i+1 to num_edges, but ever this more exhaustive search is not
      // guaranteed to make an optimal pairing... we take the longer
      // approach hoping it will give slightly better results for
      // hard problems...
    for(j = 0; j<num_edges; j++){
      if(j!=i){
        next_edge = edge_list[j];
        next_point = next_edge->other_point(point);
        e0 = point->coordinates() - next_point->coordinates();
        e1 = prev_point->coordinates() - point->coordinates();
        e0.normalize();
        e1.normalize();
        current_dot = e0%e1;
          //if the current dot satisfies the feature angle criterion
          // and is better than any other we've seen so far for this
          // given edge, save it.
        if(current_dot >= minDot && current_dot > dot_array[i]){
          dot_array[i]=current_dot;
          other_index[i]=j;
          pair_exists=true;//keep track of whether a pair has been saved
        }
      }
    }
  }
    //if there aren't any pairs, don't bother moving forward.
  if(!pair_exists)
    return CUBIT_SUCCESS;
    //now find the best pair.  That is the pair with biggest
    // dot product.  Then find the next, and so on until we
    // are done.
  double best_this_time = CUBIT_DBL_MAX;
  int best_index=-1;
    // given num_edges > 2 and each edge is paired with at most
    // one other edge at this node, there can't be more than
    // num_edges - 1 pairs.  Actually, num_edges / 2, but
    // it is just a safety check anyway, so num_edges-1 will work.
  for(i=0;i<num_edges-1 && best_this_time >= minDot; i++){
    best_this_time = -1.0;
    best_index = -1;
      //loop over and find the biggest dot
    for(j=0;j<num_edges; j++){
      if(dot_array[j] > best_this_time){
        best_this_time = dot_array[j];
        best_index = j;
      }
    }
      //if we found a pair that we can make C1
    if(best_index >= 0){
        //Don't let the above loop find either of these again
        // (unless they are the 'other' in the pair).
      dot_array[best_index] = -1.0;
      dot_array[other_index[best_index]] = -1.0;

        //First, make sure the other in the pair hasn't already
        // been used.
      CubitFacetEdge* edge_1 = edge_list[best_index];
      CubitFacetEdge* edge_2 = edge_list[other_index[best_index]];
      td_fbe = TDFacetBoundaryEdge::get_facet_boundary_edge( edge_2 );
      if(!td_fbe){
        PRINT_ERROR("Expected a tool data.\n");
        return CUBIT_FAILURE;
      }
        //figure out whether it should be stored as prev or next
      if(point == edge_2->point(0)){
          //if something has already been stored here, then
          // need to skip this pair
        if(td_fbe->get_prev_edge())
          edge_1 = NULL;
        else
          td_fbe->set_prev_edge(edge_1);
      }
      else{
        if(td_fbe->get_next_edge())
          edge_1 = NULL;
        else
          td_fbe->set_next_edge(edge_1);
      }
        //edge_1 will be NULL if we decided above to skip this pair
      if(edge_1){//otherwise save edge_2 in the appropriate spot
          // for edge_1's tool data.
        td_fbe = TDFacetBoundaryEdge::get_facet_boundary_edge( edge_1 );
        if(!td_fbe){
          PRINT_ERROR("Expected another tool data.\n");
          return CUBIT_FAILURE;
        }
        
        if(point == edge_1->point(0))
          td_fbe->set_prev_edge(edge_2);
        else
          td_fbe->set_next_edge(edge_2);
      }
    }
    
  }
  return CUBIT_SUCCESS;
}

      
  
//===========================================================================
//Function Name: pair_edges
//
//Member Type:  PRIVATE
//Descriptoin:  Basically, loops over the points and calls
//              mark_edge_pairs
//===========================================================================
CubitStatus FacetEvalTool::pair_edges()
{
  CubitStatus status = CUBIT_SUCCESS;
  int i;
  CubitPoint* point = NULL;
  for( i=0;i<myPointList.size(); i++){
    point = myPointList.get_and_step();
    status = mark_edge_pairs(point);
    if(status!=CUBIT_SUCCESS)
      return status;
  }
  return status;
}
 
//===========================================================================
//Function Name: init_bezier_surface
//
//Member Type:  PRIVATE
//Descriptoin:  compute the surface control points
//===========================================================================
CubitStatus FacetEvalTool::init_bezier_surface()
{
  // initialize the edges

  int i;
  CubitStatus status = CUBIT_SUCCESS;
  CubitFacetEdge *edge;
    // figure out which edges should be paired for C1 continuity.
  status = pair_edges();
  if(status!=CUBIT_SUCCESS)
    return status;
  
  for (i=0; i<myEdgeList.size() && status == CUBIT_SUCCESS; i++) {
    edge = myEdgeList.get_and_step();
    status = init_bezier_edge( edge, minDot );
    if (status != CUBIT_SUCCESS)
      return status;
  }
  int mydebug = 0;
  if (mydebug)
  {
    debug_draw_bezier_edges();
    dview();
  }
  
  // initialize the facets

  if (status == CUBIT_SUCCESS) {
    CubitFacet *facet;
    for (i=0; i<myFacetList.size() && status == CUBIT_SUCCESS; i++) {
      facet = myFacetList.get_and_step();
      status = init_bezier_facet( facet );
    }
  }
  if(status != CUBIT_SUCCESS){
      PRINT_ERROR("Problem initializing bezier facet.\n");
      return status;
  }

  // reset the bounding box to account for new control points

  reset_bounding_box();

  //draw_eval_bezier_facet( myFacetList.get_and_step() );
  // draw_bezier_facet( myFacetList.get_and_step() );

  return status;
}

//===========================================================================
//Function Name: init_bezier_edge
//
//Member Type:  PRIVATE
//Descriptoin:  compute the control points for an edge
//===========================================================================
CubitStatus FacetEvalTool::init_bezier_edge( CubitFacetEdge *edge,
                                             double min_dot )
{
  CubitStatus stat = CUBIT_SUCCESS;

  TDFacetBoundaryEdge *td_bfe =
    TDFacetBoundaryEdge::get_facet_boundary_edge( edge );
  if (td_bfe)
  {
    stat = td_bfe->init_control_points( min_dot );
    if (stat != CUBIT_SUCCESS)
      return stat;
    stat = td_bfe->merge_control_points();
  }
  else
  {
    CubitVector ctrl_pts[3];
    CubitVector P0 = edge->point(0)->coordinates();
    CubitVector P3 = edge->point(1)->coordinates();
    CubitVector N0 = edge->point(0)->normal( edge );
    CubitVector N3 = edge->point(1)->normal( edge );
    CubitVector T0, T3;
    if (edge->num_adj_facets() <= 1)
    {
      stat = compute_curve_tangent( edge, min_dot, T0, T3 );
      if (stat != CUBIT_SUCCESS)
        return stat;
    }
    else
    {
      T0 = P3 - P0;
      T0.normalize();
      T3 = T0;
    }
    stat = init_edge_control_points( P0, P3, N0, N3, T0, T3, ctrl_pts );
    if (stat != CUBIT_SUCCESS)
      return stat;
    edge->control_points( ctrl_pts, 4 );
  }
  return stat;
}


//===========================================================================
//Function Name: init_edge_control_points
//
//Member Type:  PRIVATE
//Descriptoin:  compute the control points for an edge
//===========================================================================
CubitStatus FacetEvalTool::init_edge_control_points( CubitVector &P0,
                                                     CubitVector &P3,
                                                     CubitVector &N0,
                                                     CubitVector &N3,
                                                     CubitVector &T0,
                                                     CubitVector &T3,
                                                     CubitVector Pi[3] )
{
  CubitVector Vi[4];
  Vi[0] = P0;
  Vi[3] = P3;
  double di = P0.distance_between( P3 );
  double ai = N0 % N3;
  double ai0 = N0 % T0;
  double ai3 = N3 % T3;
  double denom = 4 - sqr(ai);
  if (fabs(denom) < 1e-20) {
    return CUBIT_FAILURE;
  }
  double row = 6.0e0 * (2.0e0 * ai0 + ai * ai3) / denom;
  double omega = 6.0e0 * (2.0e0 * ai3 + ai * ai0) / denom;
  Vi[1] = Vi[0] + (di * (((6.0e0 * T0) - ((2.0e0 * row) * N0) + (omega * N3)) / 18.0e0));
  Vi[2] = Vi[3] - (di * (((6.0e0 * T3) + (row * N0) - ((2.0e0 * omega) * N3)) / 18.0e0));
  CubitVector Wi[3];
  Wi[0] = Vi[1] - Vi[0];
  Wi[1] = Vi[2] - Vi[1];
  Wi[2] = Vi[3] - Vi[2];

  Pi[0] = 0.25 * Vi[0] + 0.75 * Vi[1];
  Pi[1] = 0.50 * Vi[1] + 0.50 * Vi[2];
  Pi[2] = 0.75 * Vi[2] + 0.25 * Vi[3];

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: init_edge_control_points_single
//
//Member Type:  PRIVATE
//Descriptoin:  compute the control points for an edge without normals
//===========================================================================
CubitStatus FacetEvalTool::init_edge_control_points_single( 
                                     CubitVector &P0,
                                     CubitVector &P3,
                                     CubitVector &T0,
                                     CubitVector &T3,
                                     CubitVector Pi[3] )
{
  CubitVector Vi[4];
  Vi[0] = P0;
  Vi[3] = P3;
  double di = P0.distance_between( P3 );
  double ai = T0 % T3;
  double denom = 4 - sqr(ai);
  if (fabs(denom) < 1e-20) {
    return CUBIT_FAILURE;
  }
  Vi[1] = Vi[0] + (di * (((6.0e0 * T0)) / 18.0e0));
  Vi[2] = Vi[3] - (di * (((6.0e0 * T3)) / 18.0e0));

  Pi[0] = 0.25 * Vi[0] + 0.75 * Vi[1];
  Pi[1] = 0.50 * Vi[1] + 0.50 * Vi[2];
  Pi[2] = 0.75 * Vi[2] + 0.25 * Vi[3];

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: init_bezier_facet
//
//Member Type:  PRIVATE
//Descriptoin:  compute the control points for a facet
//===========================================================================
CubitStatus FacetEvalTool::init_bezier_facet( CubitFacet *facet )
{

  CubitStatus stat = CUBIT_SUCCESS;
  CubitVector P[3][5];
  CubitVector N[6], G[6];
  stat = facet->get_edge_control_points( P );
  if (stat != CUBIT_SUCCESS)
    return stat;

  // retreive the normals from the edge points.  Note we duplicate the
  // pointer data here only because of the edge_use

  int mydebug = 0;
  if (mydebug)
  {
    dcolor(CUBIT_RED);
    dfdraw(facet);
    dview();
  }
  for (int i=0; i<3; i++) {
    CubitFacetEdge *edge = facet->edge(i);
    if (facet->edge_use(i) == 1) {
      N[i*2] = edge->point(0)->normal( facet );
      N[i*2+1] = edge->point(1)->normal( facet );
    }
    else {
      N[i*2] = edge->point(1)->normal( facet );
      N[i*2+1] = edge->point(0)->normal( facet );
    }
  }

  // init the facet control points.

  stat = init_facet_control_points( N, P, G );
  if (stat != CUBIT_SUCCESS)
    return stat;
  facet->set_control_points ( G );
  facet->update_bezier_bounding_box();
  return stat;
}

//===========================================================================
//Function Name: init_facet_control_points
//
//Member Type:  PRIVATE
//Descriptoin:  compute the control points for a facet
//===========================================================================
CubitStatus FacetEvalTool::init_facet_control_points( 
  CubitVector N[6],     // vertex normals (per edge)
  CubitVector P[3][5],  // edge control points
  CubitVector G[6] )    // return internal control points
{
  CubitVector Di[4], Ai[3], N0, N3, Vi[4], Wi[3];
  double denom;
  double lambda[2], mu[2];

  CubitStatus stat = CUBIT_SUCCESS;

  for (int i=0; i<3; i++) {
    N0 = N[i*2];
    N3 = N[i*2+1];
    Vi[0] = P[i][0];
    Vi[1] = (P[i][1] - 0.25 * P[i][0]) / 0.75;
    Vi[2] = (P[i][3] - 0.25 * P[i][4]) / 0.75;
    Vi[3] = P[i][4];
    Wi[0] = Vi[1] - Vi[0];
    Wi[1] = Vi[2] - Vi[1];
    Wi[2] = Vi[3] - Vi[2];
    Di[0] = P[(i+2)%3][3] - 0.5*(P[i][1] + P[i][0]);
    Di[3] = P[(i+1)%3][1] - 0.5*(P[i][4] + P[i][3]);
    Ai[0] = (N0 * Wi[0]) / Wi[0].length();
    Ai[2] = (N3 * Wi[2]) / Wi[2].length();
    Ai[1] = Ai[0] + Ai[2];
    denom = Ai[1].length();
    Ai[1] /= denom;
    lambda[0] = (Di[0] % Wi[0]) / (Wi[0] % Wi[0]);
    lambda[1] = (Di[3] % Wi[2]) / (Wi[2] % Wi[2]); 
    mu[0] = (Di[0] % Ai[0]);
    mu[1] = (Di[3] % Ai[2]);
    G[i*2] = 0.5 * (P[i][1] + P[i][2]) + 
              0.66666666666666 * lambda[0] * Wi[1] +
              0.33333333333333 * lambda[1] * Wi[0] +
              0.66666666666666 * mu[0] * Ai[1] +
              0.33333333333333 * mu[1] * Ai[0];
    G[i*2+1] = 0.5 * (P[i][2] + P[i][3]) +
              0.33333333333333 * lambda[0] * Wi[2] +
              0.66666666666666 * lambda[1] * Wi[1] +
              0.33333333333333 * mu[0] * Ai[2] +
              0.66666666666666 * mu[1] * Ai[1];
  }
  return stat;
}

//===========================================================================
//Function Name: eval_bezier_patch
//
//Member Type:  PRIVATE
//Descriptoin:  evaluate the bezier patch defined at a facet
//===========================================================================
CubitStatus FacetEvalTool::eval_bezier_patch( CubitFacet *facet, 
                                              CubitVector &areacoord,
                                              CubitVector &pt)
{
  // interpolate internal control points

  CubitVector gctrl_pts[6];
  facet->get_control_points( gctrl_pts );
  CubitVector P_facet[3];
  if (fabs(areacoord.y() + areacoord.z()) < 1.0e-6) {
    pt = facet->point(0)->coordinates();
    return CUBIT_SUCCESS;
  }
  if (fabs(areacoord.x() + areacoord.z()) < 1.0e-6) {
    pt = facet->point(1)->coordinates();
    return CUBIT_SUCCESS;
  }
  if (fabs(areacoord.x() + areacoord.y()) < 1.0e-6) {
    pt = facet->point(2)->coordinates();
    return CUBIT_SUCCESS;
  }

  //2,1,1
  P_facet[0] = (1.0e0 / (areacoord.y() + areacoord.z())) *
               (areacoord.y() * gctrl_pts[3] +
                areacoord.z() * gctrl_pts[4]);
  //1,2,1
  P_facet[1] = (1.0e0 / (areacoord.x() + areacoord.z())) *
               (areacoord.x() * gctrl_pts[0] +
                areacoord.z() * gctrl_pts[5]);
  //1,1,2
  P_facet[2] = (1.0e0 / (areacoord.x() + areacoord.y())) *
               (areacoord.x() * gctrl_pts[1] +
                areacoord.y() * gctrl_pts[2]);

  // sum the contribution from each of the control points

  pt.set(0.0e0, 0.0e0, 0.0e0);
  CubitFacetEdge *edge;
  edge = facet->edge( 2 );
  CubitVector ctrl_pts[5];
  edge->control_points(facet, ctrl_pts);

  //i=4; j=0; k=0;
  double B = quart(areacoord.x());
  pt += B * ctrl_pts[0];

  //i=3; j=1; k=0;
  B = 4.0 * cube(areacoord.x()) * areacoord.y();
  pt += B * ctrl_pts[1];

  //i=2; j=2; k=0;
  B = 6.0 * sqr(areacoord.x()) * sqr(areacoord.y());
  pt += B * ctrl_pts[2];

  //i=1; j=3; k=0;
  B = 4.0 * areacoord.x() * cube(areacoord.y());
  pt += B * ctrl_pts[3];

  edge = facet->edge( 0 );
  edge->control_points(facet, ctrl_pts);

  //i=0; j=4; k=0;
  B = quart(areacoord.y());
  pt += B * ctrl_pts[0];

  //i=0; j=3; k=1;
  B = 4.0 * cube(areacoord.y()) * areacoord.z();
  pt += B * ctrl_pts[1];

  //i=0; j=2; k=2;
  B = 6.0 * sqr(areacoord.y()) * sqr(areacoord.z());
  pt += B * ctrl_pts[2];

  //i=0; j=1; k=3;
  B = 4.0 * areacoord.y() * cube(areacoord.z());
  pt += B * ctrl_pts[3];

  edge = facet->edge( 1 );
  edge->control_points(facet, ctrl_pts);

  //i=0; j=0; k=4;
  B = quart(areacoord.z());
  pt += B * ctrl_pts[0];

  //i=1; j=0; k=3;
  B = 4.0 * areacoord.x() * cube(areacoord.z());
  pt += B * ctrl_pts[1];

  //i=2; j=0; k=2;
  B = 6.0 * sqr(areacoord.x()) * sqr(areacoord.z());
  pt += B * ctrl_pts[2];

  //i=3; j=0; k=1;
  B = 4.0 * cube(areacoord.x()) * areacoord.z();
  pt += B * ctrl_pts[3];

  //i=2; j=1; k=1;
  B = 12.0 * sqr(areacoord.x()) * areacoord.y() * areacoord.z();
  pt += B * P_facet[0];

  //i=1; j=2; k=1;
  B = 12.0 * areacoord.x() * sqr(areacoord.y()) * areacoord.z();
  pt += B * P_facet[1];

  //i=1; j=1; k=2;
  B = 12.0 * areacoord.x() * areacoord.y() * sqr(areacoord.z());
  pt += B * P_facet[2];

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: eval_bezier_patch_normal
//
//Member Type:  PRIVATE
//Descriptoin:  evaluate the bezier patch defined at a facet
//===========================================================================
CubitStatus FacetEvalTool::eval_bezier_patch_normal( CubitFacet *facet, 
                                                     CubitVector &areacoord,
                                                     CubitVector &normal )
{
  // interpolate internal control points

  CubitVector gctrl_pts[6];
  facet->get_control_points( gctrl_pts );
  if (fabs(areacoord.y() + areacoord.z()) < 1.0e-6) {
    normal = facet->point(0)->normal(facet);
    return CUBIT_SUCCESS;
  }
  if (fabs(areacoord.x() + areacoord.z()) < 1.0e-6) {
    normal = facet->point(1)->normal(facet);
    return CUBIT_SUCCESS;
  }
  if (fabs(areacoord.x() + areacoord.y()) < 1.0e-6) {
    normal = facet->point(2)->normal(facet);
    return CUBIT_SUCCESS;
  }

  // compute the hodograph of the quartic Gregory patch

  CubitVector Nijk[10];
  hodograph(facet,areacoord,Nijk);

  // sum the contribution from each of the control points

  normal.set(0.0e0, 0.0e0, 0.0e0);

  //i=3; j=0; k=0;
  double Bsum = 0.0;
  double B = cube(areacoord.x());
  Bsum += B;
  normal += B * Nijk[0];

  //i=2; j=1; k=0;
  B = 3.0 * sqr(areacoord.x()) * areacoord.y();
  Bsum += B;
  normal += B * Nijk[1];

  //i=1; j=2; k=0;
  B = 3.0 * areacoord.x() * sqr(areacoord.y());
  Bsum += B;
  normal += B * Nijk[2];

  //i=0; j=3; k=0;
  B = cube(areacoord.y());
  Bsum += B;
  normal += B * Nijk[3];

  //i=2; j=0; k=1;
  B = 3.0 * sqr(areacoord.x()) * areacoord.z();
  Bsum += B;
  normal += B * Nijk[4];

    //i=1; j=1; k=1;
  B = 6.0 * areacoord.x() * areacoord.y() * areacoord.z();
  Bsum += B;
  normal += B * Nijk[5];

  //i=0; j=2; k=1;
  B = 3.0 * sqr(areacoord.y()) * areacoord.z();
  Bsum += B;
  normal += B * Nijk[6];

  //i=1; j=0; k=2;
  B = 3.0 * areacoord.x() * sqr(areacoord.z());
  Bsum += B;
  normal += B * Nijk[7];

  //i=0; j=1; k=2;
  B = 3.0 * areacoord.y() * sqr(areacoord.z());
  Bsum += B;
  normal += B * Nijk[8];

  //i=0; j=0; k=3;
  B = cube(areacoord.z());
  Bsum += B;
  normal += B * Nijk[9];

  assert(fabs(Bsum - 1.0) < 1e-9);

  normal.normalize();

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: hodograph
//
//Member Type:  PUBLIC
//Description:  get the hodograph control points for the facet
//Note:  This is a triangle cubic patch that is defined by the
//       normals of quartic facet control point lattice.  Returned coordinates
//       in Nijk are defined by the following diagram
//
//
//                         *9               index  polar
//                        / \                 0     300    point(0)
//                       /   \                1     210
//                     7*-----*8              2     120    
//                     / \   / \              3     030    point(1)
//                    /   \ /   \             4     201
//                  4*----5*-----*6           5     111
//                  / \   / \   / \           6     021
//                 /   \ /   \ /   \          7     102
//                *-----*-----*-----*         8     012
//                0     1     2     3         9     003    point(2)
//
//===========================================================================
CubitStatus FacetEvalTool::hodograph( CubitFacet *facet, 
                                      CubitVector &areacoord,
                                      CubitVector Nijk[10] )
{

  // compute the control points on the interior of the patch based on areacoord

  CubitVector gctrl_pts[6];
  facet->get_control_points( gctrl_pts );
  CubitVector P_facet[3];

  //2,1,1
  P_facet[0] = (1.0e0 / (areacoord.y() + areacoord.z())) *
               (areacoord.y() * gctrl_pts[3] +
                areacoord.z() * gctrl_pts[4]);
  //1,2,1
  P_facet[1] = (1.0e0 / (areacoord.x() + areacoord.z())) *
               (areacoord.x() * gctrl_pts[0] +
                areacoord.z() * gctrl_pts[5]);
  //1,1,2
  P_facet[2] = (1.0e0 / (areacoord.x() + areacoord.y())) *
               (areacoord.x() * gctrl_pts[1] +
                areacoord.y() * gctrl_pts[2]);

  // corner control points are just the normals at the points
  
  //3, 0, 0
  Nijk[0] = facet->point(0)->normal(facet);
  //0, 3, 0
  Nijk[3] = facet->point(1)->normal(facet);
  //0, 0, 3
  Nijk[9] = facet->point(2)->normal(facet);

  // fill in the boundary control points.  Define as the normal to the local
  // triangle formed by the quartic control point lattice
    
  CubitFacetEdge *edge;
  edge = facet->edge( 2 );
  CubitVector ctrl_pts[5];
  edge->control_points(facet, ctrl_pts);

  //2, 1, 0
  Nijk[1] = (ctrl_pts[2] - ctrl_pts[1]) * (P_facet[0] - ctrl_pts[1]);
  Nijk[1].normalize();

  //1, 2, 0
  Nijk[2] = (ctrl_pts[3] - ctrl_pts[2]) * (P_facet[1] - ctrl_pts[2]);
  Nijk[2].normalize();

  edge = facet->edge( 0 );
  edge->control_points(facet, ctrl_pts);

  //0, 2, 1
  Nijk[6] = (ctrl_pts[1] - P_facet[1]) * (ctrl_pts[2] - P_facet[1]);
  Nijk[6].normalize();

  //0, 1, 2
  Nijk[8] = (ctrl_pts[2] - P_facet[2]) * (ctrl_pts[3] - P_facet[2]);
  Nijk[8].normalize();

  edge = facet->edge( 1 );
  edge->control_points(facet, ctrl_pts);

  //1, 0, 2
  Nijk[7] = (P_facet[2] - ctrl_pts[2]) * (ctrl_pts[1] - ctrl_pts[2]);
  Nijk[7].normalize();

  //2, 0, 1
  Nijk[4] = (P_facet[0] - ctrl_pts[3]) * (ctrl_pts[2] - ctrl_pts[3]);
  Nijk[4].normalize();

  //1, 1, 1
  Nijk[5] = (P_facet[1] - P_facet[0]) * (P_facet[2] - P_facet[0]);
  Nijk[5].normalize();

  int mydebug = 0;
  if (mydebug)
  {
#define ONE_THIRD 0.33333333333333333333333333333333333
#define TWO_THIRDS .666666666666666666666666666666666666667
    int ii;
    CubitVector ac, pt;
    for(ii=0; ii<10; ii++)
    {
      switch(ii)
      {
      case 0: ac.set(1.0,        0.0,        0.0       ); break;
      case 1: ac.set(TWO_THIRDS, ONE_THIRD,  0.0       ); break;
      case 2: ac.set(ONE_THIRD,  TWO_THIRDS, 0.0       ); break;
      case 3: ac.set(0.0,        1.0,        0.0       ); break;
      case 4: ac.set(TWO_THIRDS, 0.0,        ONE_THIRD ); break;
      case 5: ac.set(ONE_THIRD,  ONE_THIRD,  ONE_THIRD ); break;
      case 6: ac.set(0.0,        TWO_THIRDS, ONE_THIRD ); break;
      case 7: ac.set(ONE_THIRD,  0.0,        TWO_THIRDS); break;
      case 8: ac.set(0.0,        ONE_THIRD,  TWO_THIRDS); break;
      case 9: ac.set(0.0,        0.0,        1.0       ); break;
      }
      eval_bezier_patch(facet, ac, pt);
      dray(pt,Nijk[ii],1.0);
    }
    dview();
  }

  return CUBIT_SUCCESS;
}


//===========================================================================
//Function Name: project_to_patch
//
//Member Type:  PUBLIC
//Descriptoin:  Project a point to a bezier patch. Pass in the area
//              of the point projected to the linear facet.  Function
//              assumes that the point is contained within the patch -
//              if not, it will project to one of its edges.
//===========================================================================
CubitStatus FacetEvalTool::project_to_patch(
  CubitFacet *facet,     // (IN) the facet where the patch is defined                                                
  CubitVector &ac,       // (IN) area coordinate initial guess (from linear facet)
  const CubitVector &pt,       // (IN) point we are projecting to patch
  CubitVector &eval_pt,  // (OUT) The projected point
  CubitVector *eval_norm, // (OUT) normal at evaluated point
  CubitBoolean &outside, // (OUT) the closest point on patch to pt is on an edge
  double compare_tol,    // (IN) comparison tolerance
  int edge_id )          // (IN) only used if this is to be projected to one
                         //      of the edges.  Otherwise, should be -1
{
  CubitStatus status = CUBIT_SUCCESS;

  // see if we are at a vertex

#define INCR 0.01  
  const double tol = compare_tol;

  if (is_at_vertex( facet, pt, ac, compare_tol, eval_pt, eval_norm ))
  {
    outside = CUBIT_FALSE;
    return CUBIT_SUCCESS;
  }

  // check if the start ac is inside the patch -if not, then move it there
  
  int nout = 0;
  const double atol = 0.001;
  if(move_ac_inside( ac, atol ))
    nout++;

  int diverge = 0;
  int iter = 0;
  CubitVector newpt;
  eval_bezier_patch(facet, ac, newpt);
  CubitVector move = pt - newpt;
  double lastdist = move.length();
  double bestdist = lastdist;
  CubitVector bestac = ac;
  CubitVector bestpt = newpt;
  CubitVector bestnorm;

  // If we are already close enough, then return now

  if (lastdist <= tol && !eval_norm && nout == 0) {
    eval_pt = pt;
    outside = CUBIT_FALSE;
    return status;
  }
  
  double ratio, mag, umove, vmove, det, distnew, movedist;
  CubitVector lastpt = newpt;
  CubitVector lastac = ac;
  CubitVector norm;
  CubitVector xpt, ypt, zpt, xac, yac, zac, xvec, yvec, zvec;
  CubitVector du, dv, newac;
  CubitBoolean done = CUBIT_FALSE;
  while (!done) {

    // We will be locating the projected point within the u,v,w coordinate
    // system of the triangular bezier patch.  Since u+v+w=1, the components
    // are not linearly independant.  We will choose only two of the 
    // coordinates to use and compute the third.

    int system;
    if (lastac.x() >= lastac.y() && lastac.x() >= lastac.z()) {
      system = 0;
    }
    else if (lastac.y() >= lastac.z()) {
      system = 1;
    }
    else {
      system = 2;
    }

    // compute the surface derivatives with respect to each 
    // of the barycentric coordinates

    
    if (system == 1 || system == 2) {
      xac.x( lastac.x() + INCR );
      if (lastac.y() + lastac.z() == 0.0)
        return CUBIT_FAILURE;
      ratio = lastac.z() / (lastac.y() + lastac.z());
      xac.y( (1.0 - xac.x()) * (1.0 - ratio) );
      xac.z( 1.0 - xac.x() - xac.y() );
      eval_bezier_patch(facet, xac, xpt);
      xvec = xpt - lastpt;
      xvec /= INCR;
    }
    if (system == 0 || system == 2) {
      yac.y( lastac.y() + INCR );
      if (lastac.x() + lastac.z() == 0.0)
        return CUBIT_FAILURE;
      ratio = lastac.z() / (lastac.x() + lastac.z());
      yac.x( (1.0 - yac.y()) * (1.0 - ratio) );
      yac.z( 1.0 - yac.x() - yac.y() );
      eval_bezier_patch(facet, yac, ypt);
      yvec = ypt - lastpt;
      yvec /= INCR;
    }
    if (system == 0 || system == 1) {
      zac.z( lastac.z() + INCR );
      if (lastac.x() + lastac.y() == 0.0)
        return CUBIT_FAILURE;
      ratio = lastac.y() / (lastac.x() + lastac.y());
      zac.x( (1.0 - zac.z()) * (1.0 - ratio) );
      zac.y( 1.0 - zac.x() - zac.z() );
      eval_bezier_patch(facet, zac, zpt);
      zvec = zpt - lastpt;
      zvec /= INCR;
    }

    // compute the surface normal

    switch (system) {
    case 0:
      du = yvec;
      dv = zvec;
      break;
    case 1:
      du = zvec;
      dv = xvec;
      break;
    case 2:
      du = xvec;
      dv = yvec;
      break;
    }
    norm = du * dv;
    mag = norm.length();
    if (mag < DBL_EPSILON) {
      return CUBIT_FAILURE;  
      // do something else here (it is likely a flat triangle - 
      // so try evaluating just an edge of the bezier patch)
    }
    norm /= mag;
    if (iter == 0)
      bestnorm = norm;

    // project the move vector to the tangent plane

    move = (norm * move) * norm;

    // compute an equivalent u-v-w vector

    CubitVector absnorm( fabs(norm.x()), fabs(norm.y()), fabs(norm.z()) );
    if (absnorm.z() >= absnorm.y() && absnorm.z() >= absnorm.x()) {
      det = du.x() * dv.y() - dv.x() * du.y();
      if (fabs(det) <= DBL_EPSILON) {
        return CUBIT_FAILURE;  // do something else here
      }
      umove = (move.x() * dv.y() - dv.x() * move.y()) / det;
      vmove = (du.x() * move.y() - move.x() * du.y()) / det;
    }
    else if (absnorm.y() >= absnorm.z() && absnorm.y() >= absnorm.x()) {
      det = du.x() * dv.z() - dv.x() * du.z();
      if (fabs(det) <= DBL_EPSILON) {
        return CUBIT_FAILURE;
      }
      umove = (move.x() * dv.z() - dv.x() * move.z()) / det;
      vmove = (du.x() * move.z() - move.x() * du.z()) / det;
    }
    else {
      det = du.y() * dv.z() - dv.y() * du.z();
      if (fabs(det) <= DBL_EPSILON) {
        return CUBIT_FAILURE;
      }
      umove = (move.y() * dv.z() - dv.y() * move.z()) / det;
      vmove = (du.y() * move.z() - move.y() * du.z()) / det;
    }

    /* === compute the new u-v coords and evaluate surface at new location */

    switch (system) {
    case 0:
      newac.y( lastac.y() + umove );
      newac.z( lastac.z() + vmove );
      newac.x( 1.0 - newac.y() - newac.z() );
      break;
    case 1:
      newac.z( lastac.z() + umove );
      newac.x( lastac.x() + vmove );
      newac.y( 1.0 - newac.z() - newac.x() );
      break;
    case 2:
      newac.x( lastac.x() + umove );
      newac.y( lastac.y() + vmove );
      newac.z( 1.0 - newac.x() - newac.y() );
      break;
    }

    // Keep it inside the patch

    if ( newac.x() >= -atol && 
         newac.y() >= -atol && 
         newac.z() >= -atol) {
      nout = 0;
    }
    else {
      if (move_ac_inside( newac, atol ) == CUBIT_TRUE)
        nout++;
    }

    // Evaluate at the new location

    if (edge_id != -1)
      ac_at_edge( newac, newac, edge_id );  // move to edge first
    eval_bezier_patch(facet, newac, newpt);

    // Check for convergence

    distnew = pt.distance_between(newpt);
    move = newpt - lastpt;
    movedist = move.length();
    if (movedist < tol || 
        distnew < tol ) {
      done = CUBIT_TRUE;
      if (distnew < bestdist)
      {
        bestdist = distnew;
        bestac = newac;
        bestpt = newpt;
        bestnorm = norm;
      }
    }
    else {

      // don't allow more than 30 iterations

      iter++;
      if (iter > 30) {
        //if (movedist > tol * 100.0) nout=1;
        done = CUBIT_TRUE;
      }

      // Check for divergence - don't allow more than 5 divergent
      // iterations

      if (distnew > lastdist) {
        diverge++;
        if (diverge > 10) {
          done = CUBIT_TRUE;
          //if (movedist > tol * 100.0) nout=1;
        }
      }

      // Check if we are continuing to project outside the facet.  
      // If so, then stop now

      if (nout > 3) {
        done = CUBIT_TRUE;
      }

      // set up for next iteration

      if (!done) {
        if (distnew < bestdist)
        {
          bestdist = distnew;
          bestac = newac;
          bestpt = newpt;
          bestnorm = norm;
        }
        lastdist = distnew;
        lastpt = newpt;
        lastac = newac;
        move = pt - lastpt;
      }
    }
  }

  eval_pt = bestpt;
  if (eval_norm) {
    *eval_norm = bestnorm;
  }
  outside = (nout > 0) ? CUBIT_TRUE : CUBIT_FALSE;
  ac = bestac;

  return status;
}

//===========================================================================
//Function Name: is_at_vertex
//
//Member Type:  PRIVATE
//Description:  determine if the point is at one of the facet's vertices
//===========================================================================
CubitBoolean FacetEvalTool::is_at_vertex( 
  CubitFacet *facet,  // (IN) facet we are evaluating 
  const CubitVector &pt,    // (IN) the point
  CubitVector &ac,    // (IN) the ac of the point on the facet plane
  double compare_tol, // (IN) return TRUE of closer than this
  CubitVector &eval_pt, // (OUT) location at vertex if TRUE
  CubitVector *eval_norm_ptr ) // (OUT) normal at vertex if TRUE
{
  double dist;
  CubitVector vert_loc;
  const double actol = 0.1;
  if (fabs(ac.x()) < actol && fabs(ac.y()) < actol) {
    vert_loc = facet->point(2)->coordinates();
    dist = pt.distance_between( vert_loc );
    if (dist <= compare_tol)
    {
      eval_pt = vert_loc;
      if (eval_norm_ptr) 
      {
        *eval_norm_ptr = facet->point(2)->normal( facet );
      }
      return CUBIT_TRUE;
    }
  }

  if (fabs(ac.x()) < actol && fabs(ac.z()) < actol) {
    vert_loc = facet->point(1)->coordinates();
    dist = pt.distance_between( vert_loc );
    if (dist <= compare_tol)
    {
      eval_pt = vert_loc;
      if (eval_norm_ptr) 
      {
        *eval_norm_ptr = facet->point(1)->normal( facet );
      }
      return CUBIT_TRUE;
    }
  }

  if (fabs(ac.y()) < actol && fabs(ac.z()) < actol) {
    vert_loc = facet->point(0)->coordinates();
    dist = pt.distance_between( vert_loc );
    if (dist <= compare_tol)
    {
      eval_pt = vert_loc;
      if (eval_norm_ptr) 
      {
        *eval_norm_ptr = facet->point(0)->normal( facet );
      }
      return CUBIT_TRUE;
    }
  }

  return CUBIT_FALSE;
}

//===========================================================================
//Function Name: ac_at_edge
//
//Member Type:  PRIVATE
//Description:  determine the area coordinate of the facet at the edge
//===========================================================================
void FacetEvalTool::ac_at_edge( CubitVector &fac,  // facet area coordinate
                                CubitVector &eac,   // edge area coordinate
                                int edge_id )      // id of edge
{
  double u, v, w;
  switch (edge_id)
  {
  case 0:
    u = 0.0;
    v = fac.y() / (fac.y() + fac.z());
    w = 1.0 - v;
    break;
  case 1:
    u = fac.x() / (fac.x() + fac.z());
    v = 0.0;
    w = 1.0 - u;
    break;
  case 2:
    u = fac.x() / (fac.x() + fac.y());
    v = 1.0 - u;
    w = 0.0;
    break;
  default:
    assert(0);
    break;
  }
  eac.set(u, v, w);
}

//===========================================================================
//Function Name: move_ac_inside
//
//Member Type:  PRIVATE
//Description:  find the closest area coordinate to the boundary of the
//              patch if any of its components are < 0
//              Return if the ac was modified.
//===========================================================================
CubitBoolean FacetEvalTool::move_ac_inside( CubitVector &ac, double tol )
{
  int nout = 0;
  if (ac.x() < -tol) {
    ac.x( 0.0 );
    ac.y( ac.y() / (ac.y() + ac.z()) );
    ac.z( 1.0 - ac.y() );
    nout++;
  }
  if (ac.y() < -tol) {
    ac.y( 0.0 );
    ac.x( ac.x() / (ac.x() + ac.z()) );
    ac.z( 1.0 - ac.x() );
    nout++;
  }
  if (ac.z() < -tol) {
    ac.z( 0.0 );
    ac.x( ac.x() / (ac.x() + ac.y()) );
    ac.y( 1.0 - ac.x() );
    nout++;
  }
  return (nout > 0) ? CUBIT_TRUE : CUBIT_FALSE;
}

bool FacetEvalTool::have_data_to_calculate_bbox(void)
{
  if(myPointList.size() > 0 ||
     (interpOrder == 4 && 
         (myEdgeList.size() > 0 || 
          myFacetList.size() > 0)))
  {
    return true;
  }
  return false;
}

//===========================================================================
//Function Name: reset_bounding_box
//
//Member Type:  PUBLIC
//Descriptoin:  Calculates the bounding box of the surface (also sets
// the compareTol and grid search).
//===========================================================================
void FacetEvalTool::reset_bounding_box()
{
  if(have_data_to_calculate_bbox())
  {
    if (myBBox != NULL)
    {
      delete myBBox;
      myBBox = NULL;
    }
    
    bounding_box();

    double diag = sqrt(sqr(myBBox->x_range()) + 
                       sqr(myBBox->y_range()) + 
                       sqr(myBBox->z_range()));
    compareTol = 1.0e-3 * diag;

    set_up_grid_search( diag );
  }
}


//===========================================================================
//Function Name: bounding_box
//
//Member Type:  PUBLIC
//Descriptoin:  Calculates the bounding box of the surface.
//===========================================================================
CubitBox FacetEvalTool::bounding_box()
{
  if ( !myBBox )
  {
    int ii;
    CubitPoint *point_ptr;
    double x, y, z;
    double x_min = CUBIT_DBL_MAX, x_max = -CUBIT_DBL_MAX;
    double y_min = CUBIT_DBL_MAX, y_max = -CUBIT_DBL_MAX;
    double z_min = CUBIT_DBL_MAX, z_max = -CUBIT_DBL_MAX;
    for ( ii = myPointList.size(); ii > 0; ii-- )
    {
      point_ptr = myPointList.get_and_step();
      x = point_ptr->x();
      y = point_ptr->y();
      z = point_ptr->z();
      if ( x < x_min )
        x_min = x;
      if ( y < y_min )
        y_min = y;
      if ( z < z_min )
        z_min = z;
      if ( x > x_max )
        x_max = x;
      if ( y > y_max )
        y_max = y;
      if ( z > z_max )
        z_max = z;
    }
    if (interpOrder == 4)
    {
      CubitFacetEdge *edge_ptr;
      CubitVector ctrl_pts[6];
      int jj;
      for ( ii = myEdgeList.size(); ii > 0; ii-- )
      {
        edge_ptr = myEdgeList.get_and_step();
        edge_ptr->control_points(ctrl_pts);
        for (jj=1; jj<4; jj++)
        {
          x = ctrl_pts[jj].x();
          y = ctrl_pts[jj].y();
          z = ctrl_pts[jj].z();
          if ( x < x_min )
            x_min = x;
          if ( y < y_min )
            y_min = y;
          if ( z < z_min )
            z_min = z;
          if ( x > x_max )
            x_max = x;
          if ( y > y_max )
            y_max = y;
          if ( z > z_max )
            z_max = z;
        }
      }
      CubitFacet *facet_ptr;
      for ( ii = myFacetList.size(); ii > 0; ii-- )
      {
        facet_ptr = myFacetList.get_and_step();
        facet_ptr->get_control_points(ctrl_pts);
        for (jj=0; jj<6; jj++)
        {
          x = ctrl_pts[jj].x();
          y = ctrl_pts[jj].y();
          z = ctrl_pts[jj].z();
          if ( x < x_min )
            x_min = x;
          if ( y < y_min )
            y_min = y;
          if ( z < z_min )
            z_min = z;
          if ( x > x_max )
            x_max = x;
          if ( y > y_max )
            y_max = y;
          if ( z > z_max )
            z_max = z;
        }
      }
    }
    CubitVector min_v(x_min, y_min, z_min );
    CubitVector max_v(x_max, y_max, z_max );
    myBBox = new CubitBox( min_v, max_v );
  }
  return *(myBBox);
}

//===========================================================================
//Function Name: closest_point
//
//Member Type:  PUBLIC
//Description:  Finds the closest point from the vector (this_point) to the
//              set of facets that lies on the set of facets.  If the point
//              lies outside this set, the closest point will be on the plane
//              of the closest facet.  The closest_point is set to be that point.
//===========================================================================
CubitStatus FacetEvalTool::closest_point(CubitVector &this_point,
                                         CubitVector *closest_point_ptr,
                                         CubitVector *normal_ptr)
{
  CubitBoolean trim = CUBIT_FALSE;
  CubitBoolean outside;
  CubitStatus rv = CUBIT_SUCCESS;
  static int count = 0;  count++;

  int mydebug = 0;
  if (mydebug)
  {
    debug_draw_vec( this_point, CUBIT_RED );
  }

  rv = project_to_facets( this_point, trim, &outside, 
                          closest_point_ptr, normal_ptr );

  if (DEBUG_FLAG(49))
  {
    if (closest_point_ptr)
    {
      double dist = closest_point_ptr->distance_between( this_point );
      if (dist > 1.0)
      {
        PRINT_ERROR("Appears to be bad projection in FacetEvalTool::project_to_facets\n");
      }
    }
  }

  if (mydebug && closest_point_ptr)
  {
    debug_draw_vec( *closest_point_ptr, CUBIT_GREEN );
  }
  return rv;
}


CubitFacet* FacetEvalTool::closest_facet( const CubitVector& point )
{
  CubitVector non_const(point);
  CubitBoolean junk;
  CubitStatus result = project_to_facets( non_const, true, &junk, 0, 0 ); 
  return result ? lastFacet : 0;
}

//===========================================================================
//Function Name: facets_from_search_grid
//
//Member Type:  PRIVATE
//Description:  find the closest facets to the point in the search grid
// by starting with a default tolerance and expanding until facets are found
//===========================================================================
void FacetEvalTool::facets_from_search_grid( 
  CubitVector &this_point,
  DLIList<CubitFacet *> &facet_list,
  double &tol_used)
{
  double tol = compareTol * 10;
  while (facet_list.size() == 0)
  {
    CubitVector ptmin( this_point.x() - tol, 
                       this_point.y() - tol, 
                       this_point.z() - tol );
    CubitVector ptmax( this_point.x() + tol, 
                       this_point.y() + tol, 
                       this_point.z() + tol );
    CubitBox ptbox(ptmin,ptmax);
    aTree->find(ptbox,facet_list);

    if (0 == facet_list.size())
      tol *= 2.0;
  }

  tol_used = tol;
}

//===========================================================================
//Function Name: facets_from_search_grid
//
//Member Type:  PRIVATE
//Description:  find the closest facets to the point in the search grid
// for a given tolerance
//===========================================================================
void FacetEvalTool::facets_from_search_grid( 
  CubitVector &this_point,
  double compare_tol,
  DLIList<CubitFacet *> &facet_list )
{
  double tol = compare_tol;

  CubitVector ptmin( this_point.x() - tol, 
                     this_point.y() - tol, 
                     this_point.z() - tol );
  CubitVector ptmax( this_point.x() + tol, 
                     this_point.y() + tol, 
                     this_point.z() + tol );
  CubitBox ptbox(ptmin,ptmax);
  aTree->find(ptbox,facet_list);
}

//===========================================================================
//Function Name: project_to_facets
//
//Member Type:  PRIVATE
//Description:  Project a point to the facets.  Use the interpOrder.
//              if trim is set, then trim the point to the boundary.
//              This is a non-static version.  it uses the facets
//              in the evaltool
//===========================================================================
CubitStatus 
FacetEvalTool::project_to_facets(CubitVector &this_point,
                                 CubitBoolean trim,
                                 CubitBoolean *outside,
                                 CubitVector *closest_point_ptr,
                                 CubitVector *normal_ptr)
{
  CpuTimer function_time;
  if (DEBUG_FLAG(110))
  {
    function_time.cpu_secs();
    numEvals++;
  }

  CubitStatus rv = CUBIT_SUCCESS;

  // if there are a lot of facets on this surface - use the grid search first 
  // to narrow the selection

  if (aTree != NULL)
  {
    DLIList<CubitFacet *> facet_list;
    double search_tol = DBL_MAX;
    facets_from_search_grid( this_point, facet_list, search_tol );
    if ( DEBUG_FLAG(110) )
      timeGridSearch += function_time.cpu_secs();

    if (facet_list.size())
    {
      CubitVector grid_close_pt;
      rv = project_to_facets(facet_list,lastFacet,interpOrder,compareTol,
                             this_point,trim,outside, &grid_close_pt,
                             normal_ptr);

      if (CUBIT_SUCCESS == rv)
      {
        if (closest_point_ptr)
        {
          *closest_point_ptr = grid_close_pt;
        }

        // when we do the projection, if we end up with the closest point being farther
        // away than the grid search tolerance, we may have missed a closer facet
        // so redo the grid search using the distance as a tolerance
        double distance = grid_close_pt.distance_between( this_point );
        if (distance > search_tol)
        {
          DLIList<CubitFacet*> facets_within_distance;
          CubitVector grid_close_pt2;
          facets_from_search_grid(this_point, distance, facets_within_distance);
          if (facets_within_distance.size() &&
              (facets_within_distance != facet_list) )
          {
            rv = project_to_facets(facets_within_distance, lastFacet, interpOrder, compareTol,
                                   this_point, trim, outside, &grid_close_pt2,
                                   normal_ptr);
            if (CUBIT_SUCCESS == rv)
            {
              double distance2 = grid_close_pt2.distance_between( this_point );
              if (distance2 < distance)
              {
                if (closest_point_ptr)
                {
                  *closest_point_ptr = grid_close_pt2;
                }
              }
            }
          }
        }
      }

      if ( DEBUG_FLAG(110) )
      {
        timeFacetProject += function_time.cpu_secs();
        if (closest_point_ptr)
        {
          double dist = closest_point_ptr->distance_between( this_point );
          if (dist > compareTol * 100)
          {
            PRINT_ERROR("Appears to be bad projection in FacetEvalTool::project_to_facets\n");
            dcolor(CUBIT_GREEN);
            dpoint(this_point);
            dcolor(CUBIT_RED);
            dpoint(*closest_point_ptr);
            dcolor(CUBIT_YELLOW);
            dfldraw(facet_list);
            dview();
            rv = CUBIT_FAILURE;
          }
        }
      }
    }
  }

  // otherwise just use the complete list of facets

  else
  {
    rv = project_to_facets(myFacetList,lastFacet,interpOrder,compareTol,
                           this_point,trim,outside,closest_point_ptr,
                           normal_ptr);
    if ( DEBUG_FLAG(110) )
      timeFacetProject += function_time.cpu_secs();
  }

  return rv;
}

//===========================================================================
//Function Name: project_to_facets
//
//Member Type:  PRIVATE
//Description:  Project a point to the facets.  Use the interpOrder.
//              if trim is set, then trim the point to the boundary.
//              This is a static version of the above.  Any list of facets
//              can be passed to this function
//===========================================================================
CubitStatus 
FacetEvalTool::project_to_facets(
  DLIList <CubitFacet *> &facet_list,  // (IN) facets that we can project to
  CubitFacet *&last_facet,             // (IN/OUT) last facet projected to - 
                                       //          it will try this one first
  int interp_order,                    // (IN) 0 = linear facets, 
                                       //      4 = b-spline patches
  double compare_tol,                  // (IN) tolerance for projection - 
                                       //      should be about 1e-3*edge
  const CubitVector &this_point,       // (IN) point we are projecting
  CubitBoolean trim,                   // (IN) trim to facet (always trimmed 
                                       //      if b-spline patch)
  CubitBoolean *outside,               // (OUT) TRUE if projected outside 
                                       //       the facet
  CubitVector *closest_point_ptr,      // (OUT) resulting projection point 
                                       //       (NULL if only want normal)
  CubitVector *normal_ptr)             // (OUT) resulting normal at projected 
                                       //       point (NULL if not required)
{
  int ncheck, ii, nincr=0;
  static int calls=0;
  static int nncheck=0;
  static int ntol=0;
  static int mydebug=0;
  CubitBoolean outside_facet, best_outside_facet;
  CubitVector close_point, best_point, best_areacoord;
  CubitFacet *best_facet = NULL; 
  CubitFacet *facet;
  assert (facet_list.size() > 0);
  double big_dist = compare_tol * 1.0e3;

  // set the first facet to be checked as the last one located

  if (last_facet) {
    if (!facet_list.move_to(last_facet)) {
      facet_list.reset();
    }
  }
  else {
    facet_list.reset();
  }

  // so we don't evaluate a facet more than once - mark the facets
  // as we evaluate them.  Put the evaluated facets on a used_facet_list
  // so we clear the marks off when we are done.  Note: this assumes
  // theat marks are initially cleared.
  
  DLIList<CubitFacet *>used_facet_list;
  for(ii=0; ii<facet_list.size(); ii++)
    facet_list.get_and_step()->marked(0);

  int nfacets = facet_list.size();
  int nevald = 0;
  double tol = compare_tol * 10;
  const double atol = 0.001;
  double mindist = CUBIT_DBL_MAX;
  CubitBoolean eval_all = CUBIT_FALSE;
  CubitBoolean done = CUBIT_FALSE;
  best_outside_facet = CUBIT_TRUE;
  
  while(!done) {

    // define a bounding box around the point

    double ptmin_x = this_point.x() - tol;
    double ptmin_y = this_point.y() - tol;
    double ptmin_z = this_point.z() - tol;
    double ptmax_x = this_point.x() + tol;
    double ptmax_y = this_point.y() + tol;
    double ptmax_z = this_point.z() + tol;

    ncheck = 0;
    for ( ii = facet_list.size(); ii > 0 && !done; ii-- ) 
    {
      facet = facet_list.get_and_step();
      if (facet->marked())
        continue;

      // Try to trivially reject this facet with a bounding box test
      // (Does the bounding box of the facet intersect with the 
      // bounding box of the point)

      if (!eval_all)
      {
        const CubitBox &bbox = facet->bounding_box();
        if (ptmax_x < bbox.min_x() ||
            ptmin_x > bbox.max_x()) {
          continue;
        }
        if (ptmax_y < bbox.min_y() ||
            ptmin_y > bbox.max_y()) {
          continue;
        }
        if (ptmax_z < bbox.min_z() ||
            ptmin_z > bbox.max_z()) {
          continue;
        }
      }

      // skip zero area facets
      if(facet->area() <= 0.0)
        continue;

      // Only facets that pass the bounding box test will get past here!

      // Project point to plane of the facet and determine its area coordinates

      ncheck++; nevald++;
      CubitVector pt_on_plane;
      double dist_to_plane;
      project_to_facet_plane( facet, this_point, pt_on_plane, dist_to_plane );

      CubitVector areacoord;
      facet_area_coordinate( facet, pt_on_plane, areacoord );
      if (interp_order != 0)
      {

        // modify the areacoord - project to the bezier patch- snaps to the 
        // edge of the patch if necessary

        if (project_to_facet( facet, this_point, areacoord, 
                        close_point, outside_facet, compare_tol ) != CUBIT_SUCCESS) 
        {
          return CUBIT_FAILURE;
        }
      }
      else
      {

        // If sign of areacoords are all positive then its inside the triangle
        // and we are done - go interpolate the point. (use an absolute
        // tolerance since the coordinates arenormalized)
        if (areacoord.x() > -atol && 
            areacoord.y() > -atol && 
            areacoord.z() > -atol) {
          if (dist_to_plane < compare_tol) {
            close_point = this_point;
          }
          else
          {
            close_point = pt_on_plane;
          }
          outside_facet = CUBIT_FALSE;
        }

        // otherwise find the closest vertex or edge to the projected point

        else if (areacoord.x() < atol) {
          outside_facet = CUBIT_TRUE;
          if (areacoord.y() < atol) {
            if (eval_point( facet, 2, close_point ) 
              != CUBIT_SUCCESS) {
              return CUBIT_FAILURE;
            }
          }
          else if(areacoord.z() < atol) {
            if (eval_point( facet, 1, close_point ) 
              != CUBIT_SUCCESS) {
              return CUBIT_FAILURE;
            }
          }
          else {
            if (project_to_facetedge( facet, 1, 2, this_point, pt_on_plane, 
                           close_point, outside_facet, trim ) !=CUBIT_SUCCESS) {
              return CUBIT_FAILURE;
            }
          }
        }
        else if (areacoord.y() < atol) {
          outside_facet = CUBIT_TRUE;
          if (areacoord.z() < atol) {
            if (eval_point( facet, 0, close_point ) 
              != CUBIT_SUCCESS) {
              return CUBIT_FAILURE;
            }
          }
          else {
            if (project_to_facetedge( facet, 2, 0, this_point, pt_on_plane, 
                           close_point, outside_facet, trim ) !=CUBIT_SUCCESS) {
              return CUBIT_FAILURE;
            }
          }
        }
        else {
          outside_facet = CUBIT_TRUE;
          if (project_to_facetedge( facet, 0, 1, this_point, pt_on_plane, 
                         close_point, outside_facet, trim ) !=CUBIT_SUCCESS) {
            return CUBIT_FAILURE;
          }
        }
      }
      
      // keep track of the minimum distance

      double dist = close_point.distance_between( this_point );
      if ((best_outside_facet == outside_facet && dist < mindist) ||
          (best_outside_facet && !outside_facet && (dist < big_dist || !best_facet)) )
      {
        mindist = dist;
        best_point = close_point;
        best_facet = facet;
        best_areacoord = areacoord;
        best_outside_facet = outside_facet;

        if (dist < compare_tol) {
          done = CUBIT_TRUE;
        }
        big_dist = 10.0 * mindist;
      }
      facet->marked(1);
      used_facet_list.append(facet);
    }

    // We are done if we found at least one triangle.  Otherwise
    // increase the tolerance and try again

    if (!done)
    {
      if (nevald == nfacets)
      {
        done = CUBIT_TRUE;
      }
      else
      {
        nincr++;
        if (ncheck > 0) {
          if (best_outside_facet) {
            if (nincr < 10)
            {
              tol *= 2.0;
              ntol++;
            }
            else
            // getting here means that the compare_tol probably is too small
            // just try all the remaining facets
            {
              eval_all = CUBIT_TRUE;
            }
          }
          else
          {
            done = CUBIT_TRUE;
          }
        }
        else {
          if (nincr >= 10)
          {
            eval_all = CUBIT_TRUE;
          }
          else
          {
            tol *= 2.0e0;
            ntol++;
          }
        }
      }
    }
  }  // while(!done)
  if(best_facet == NULL){
      PRINT_ERROR("Unable to determine facet correctly.\n");
      return CUBIT_FAILURE;
  }
  // make sure we actually got something
  assert(best_facet != NULL);

  // if the closest point is outside of a facet, then evaluate the point
  // on the facet using its area coordinates (otherwise it would be 
  // trimmed to an edge or point)


  if ( !trim && best_outside_facet && interp_order != 4) {
    if (project_to_facet( best_facet, this_point, best_areacoord, 
      best_point, best_outside_facet, compare_tol ) 
      != CUBIT_SUCCESS) {
      return CUBIT_FAILURE;
    }
    
    // see if its really outside (it could just be on an edge where the
    // curvature is convex)

    best_outside_facet = is_outside( best_facet, best_areacoord );
  }

  // evaluate the normal if required

  if (normal_ptr) {
    CubitVector normal;
    if (eval_facet_normal( best_facet, best_areacoord, normal ) 
      != CUBIT_SUCCESS) {
      return CUBIT_FAILURE;
    }
    *normal_ptr = normal;
  }

  if (closest_point_ptr) {
    *closest_point_ptr = best_point;
  }
  
  *outside = best_outside_facet;
  last_facet = best_facet;


  // clear the marks from the used facets

  for (ii=0; ii<used_facet_list.size(); ii++)
  {
    facet = used_facet_list.get_and_step();
    facet->marked( 0 );
  }

  if (mydebug) {
    nncheck+= ncheck;
    calls++;
    if (calls%100==0){
      PRINT_INFO("calls = %d, ckecks = %d, ntol = %d\n",calls,nncheck,ntol);
    }
  }
  
  return CUBIT_SUCCESS;
}


//===========================================================================
//Function Name: eval_facet
//
//Member Type:  PUBLIC
//Descriptoin:  Evaluate the location and normal of a set of area coordinates  
//              on a facet.  Use the interpOrder to evaluate.
//              Static function 
//===========================================================================
CubitStatus FacetEvalTool::eval_facet( CubitFacet *facet, 
                                       CubitVector &areacoord,
                                       CubitVector *eval_point,
                                       CubitVector *eval_normal )
{
  CubitStatus rv = CUBIT_SUCCESS;
  int mydebug = 0;
  if (mydebug)
  {
    dcolor( CUBIT_RED );
    dfdraw( facet );
    dview();
    dcolor( CUBIT_YELLOW );
  }

  CubitPoint *pt0 = facet->point(0);
  CubitPoint *pt1 = facet->point(1);
  CubitPoint *pt2 = facet->point(2);
  CubitVector close_point;

  CubitVector this_point;
  this_point.x( areacoord.x() * pt0->x() +
                areacoord.y() * pt1->x() +
                areacoord.z() * pt2->x() );
  this_point.y( areacoord.x() * pt0->y() +
                areacoord.y() * pt1->y() +
                areacoord.z() * pt2->y() );
  this_point.z( areacoord.x() * pt0->z() +
                areacoord.y() * pt1->z() +
                areacoord.z() * pt2->z() );

  int eval_order = facet->eval_order();
  if (eval_order != 0)
  {
    if (facet->is_flat())
      eval_order = 0;
  }

  switch(eval_order) {
  case 0:
    if (eval_point)
      *eval_point = this_point;
    if (eval_normal)
      eval_facet_normal(facet, areacoord, *eval_normal);
    break;

  case 4:
    eval_bezier_patch( facet, 
                       areacoord,
                       close_point );
    //project_to_patch(facet, areacoord, this_point, close_point, 
    //                 eval_normal, outside );

    //for now over-ride the normal from project_to_patch -- it is bogus
    if (eval_normal)
      eval_facet_normal(facet, areacoord, *eval_normal);
    if (eval_point)
      *eval_point = close_point;

    break;
  default:
    // the interpolation order for now is limited to 0 and 4
    // something other than that is being attempted (or eval_order is not
    // returning the correct value)
    assert(0);
    break;
  }

  return rv;
}

//===========================================================================
//Function Name: eval_facet
//
//Member Type:  PUBLIC
//Descriptoin:  Evaluate the location of a set of area coordinates  
//              on a facet.  Use the interpOrder to evaluate. 
//===========================================================================
CubitStatus FacetEvalTool::eval_facet( CubitFacet *facet, 
                                       CubitVector &pt,
                                       CubitVector &areacoord, 
                                       CubitVector &close_point,
                                       CubitBoolean &outside_facet )
{

  CubitPoint *pt0 = facet->point(0);
  CubitPoint *pt1 = facet->point(1);
  CubitPoint *pt2 = facet->point(2);

  int interp_order = facet->eval_order();
  if (interp_order != 0 && facet->is_flat())
  {
    interp_order = 0;
  }

  switch(interp_order) {
  case 0:
    close_point.x( areacoord.x() * pt0->x() +
                   areacoord.y() * pt1->x() +
                   areacoord.z() * pt2->x() );
    close_point.y( areacoord.x() * pt0->y() +
                   areacoord.y() * pt1->y() +
                   areacoord.z() * pt2->y() );
    close_point.z( areacoord.x() * pt0->z() +
                   areacoord.y() * pt1->z() +
                   areacoord.z() * pt2->z() );
    outside_facet = CUBIT_FALSE;
    break;
  case 1:
    {
      CubitVector tp0 = pt0->project_to_tangent_plane( pt );
      CubitVector tp1 = pt1->project_to_tangent_plane( pt );
      CubitVector tp2 = pt2->project_to_tangent_plane( pt );
      close_point.x( areacoord.x() * tp0.x() +
                     areacoord.y() * tp1.x() +
                     areacoord.z() * tp2.x() );
      close_point.y( areacoord.x() * tp0.y() +
                     areacoord.y() * tp1.y() +
                     areacoord.z() * tp2.y() );
      close_point.z( areacoord.x() * tp0.z() +
                     areacoord.y() * tp1.z() +
                     areacoord.z() * tp2.z() );
      outside_facet = CUBIT_FALSE;
    }
    break;
  case 2:
    {
      CubitVector qp0, qp1, qp2, qn0, qn1, qn2;
      eval_quadratic( facet, 0, pt, qp0, qn0 );
      eval_quadratic( facet, 1, pt, qp1, qn1 );
      eval_quadratic( facet, 2, pt, qp2, qn2 );

      close_point.x( areacoord.x() * qp0.x() +
                     areacoord.y() * qp1.x() +
                     areacoord.z() * qp2.x() );
      close_point.y( areacoord.x() * qp0.y() +
                     areacoord.y() * qp1.y() +
                     areacoord.z() * qp2.y() );
      close_point.z( areacoord.x() * qp0.z() +
                     areacoord.y() * qp1.z() +
                     areacoord.z() * qp2.z() );
      outside_facet = CUBIT_FALSE;
    }
    break;
  case 3:
    {
      CubitVector qp0, qp1, qp2;
      eval_quadric( facet, 0, pt, qp0 );
      eval_quadric( facet, 1, pt, qp1 );
      eval_quadric( facet, 2, pt, qp2 );

      close_point.x( areacoord.x() * qp0.x() +
                     areacoord.y() * qp1.x() +
                     areacoord.z() * qp2.x() );
      close_point.y( areacoord.x() * qp0.y() +
                     areacoord.y() * qp1.y() +
                     areacoord.z() * qp2.y() );
      close_point.z( areacoord.x() * qp0.z() +
                     areacoord.y() * qp1.z() +
                     areacoord.z() * qp2.z() );
      outside_facet = CUBIT_FALSE;
    }
    break;
  case 4:
    {
      //CubitStatus stat = eval_bezier_patch( facet, areacoord, pt );
      //close_point = pt;
      //outside_facet = CUBIT_FALSE;
      double compare_tol = /*sqrt(facet->area()) **/ 1.0e-3;
      int edge_id = -1;
      project_to_patch(facet, areacoord, pt, close_point, NULL, 
                       outside_facet, compare_tol, edge_id );
    }
    break;
  }

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: project_to_facet
//
//Member Type:  PUBLIC
//Description:  project to a single facet.  Uses the input areacoord as
//              a starting guess.
//===========================================================================
CubitStatus FacetEvalTool::project_to_facet( CubitFacet *facet, 
                                       const CubitVector &pt,
                                       CubitVector &areacoord, 
                                       CubitVector &close_point,
                                       CubitBoolean &outside_facet,
                                       double compare_tol)
{

  CubitStatus stat = CUBIT_SUCCESS;
  CubitPoint *pt0 = facet->point(0);
  CubitPoint *pt1 = facet->point(1);
  CubitPoint *pt2 = facet->point(2);

  int interp_order = facet->eval_order();
  if (facet->is_flat())
  {
    interp_order = 0;
  }

  switch(interp_order) {
  case 0:
    close_point.x( areacoord.x() * pt0->x() +
                   areacoord.y() * pt1->x() +
                   areacoord.z() * pt2->x() );
    close_point.y( areacoord.x() * pt0->y() +
                   areacoord.y() * pt1->y() +
                   areacoord.z() * pt2->y() );
    close_point.z( areacoord.x() * pt0->z() +
                   areacoord.y() * pt1->z() +
                   areacoord.z() * pt2->z() );
    outside_facet = CUBIT_FALSE;
    break;
  case 1:
  case 2:
  case 3:
    assert(0);  // not available from this function 
    break;
  case 4:
    {
      int edge_id = -1;
      stat = project_to_patch(facet, areacoord, pt, close_point, NULL, 
                              outside_facet, compare_tol, edge_id );
    }
    break;
  }

  return stat;
}

//===========================================================================
//Function Name: project_to_facetedge
//
//Member Type:  PUBLIC
//Description:  Project a point to the facet edge.  
//              Use the interpOrder to evaluate. 
//===========================================================================
CubitStatus FacetEvalTool::project_to_facetedge( CubitFacet *facet, 
                                      int vert0, int vert1,
                                      const CubitVector &the_point,
                                      CubitVector &pt_on_plane, 
                                      CubitVector &close_point,
                                      CubitBoolean &outside_facet,
                                      CubitBoolean must_be_on_edge)
{
  CubitPoint *pt0 = facet->point(vert0);
  CubitPoint *pt1 = facet->point(vert1);

  // the edge vector

  CubitVector e0 ( pt1->x() - pt0->x(),
                   pt1->y() - pt0->y(),
                   pt1->z() - pt0->z() );
  e0.normalize();
  
  // vector from vert0 to point

  CubitVector v0 ( pt_on_plane.x() - pt0->x(),
                   pt_on_plane.y() - pt0->y(),
                   pt_on_plane.z() - pt0->z() );

  CubitVector v1 ( pt_on_plane.x() - pt1->x(),
                   pt_on_plane.y() - pt1->y(),
                   pt_on_plane.z() - pt1->z() );
  
  // project to edge

  double projection1 = v0%e0;
  double projection2 = v1%(-e0);

  if( !must_be_on_edge || (projection1 > 0 && projection2 > 0 ))
  {
    close_point.x ( pt0->x() + e0.x() * projection1 );
    close_point.y ( pt0->y() + e0.y() * projection1 );
    close_point.z ( pt0->z() + e0.z() * projection1 );
  }
  else //we are closer to one a facet vertex than to the edge
  {
    if(   the_point.distance_between( pt0->coordinates() ) 
        < the_point.distance_between( pt1->coordinates())) 
      close_point = pt0->coordinates();
    else 
      close_point = pt1->coordinates();
  }

  // project the point on the facet (if the order is higher than 0)

  outside_facet = CUBIT_TRUE;
  if (facet->eval_order() > 0 && !facet->is_flat()) {
    CubitVector areacoord;
    facet_area_coordinate( facet, close_point, areacoord ); 
    int edge_id;
    if ((vert0 == 0 && vert1 == 1) ||
        (vert0 == 1 && vert1 == 0)) 
    {
      edge_id = 2;
    }
    else if ((vert0 == 1 && vert1 == 2) ||
             (vert0 == 2 && vert1 == 1))
    {
      edge_id = 0;
    }
    else if ((vert0 == 2 && vert1 == 0) ||
             (vert0 == 0 && vert1 == 2))
    {
      edge_id = 1;
    }
    else 
    {
      assert(0);  //edge_id wasn't set
    }

    double compare_tol = projection1 * 1.0e-3;
    if (project_to_patch( facet, areacoord, the_point, close_point, 
                          NULL, outside_facet, compare_tol, edge_id )!= CUBIT_SUCCESS) {
      return CUBIT_FAILURE;     
    }
  }

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: eval_edge
//
//Member Type:  PUBLIC
//Description:  Evaluate the location of a point projected to a 
//              linear edge. 
//===========================================================================
CubitStatus FacetEvalTool::eval_edge( CubitFacet *facet, 
                                      int vert0, int vert1,
                                      CubitVector &pt_on_plane, 
                                      CubitVector &close_point )
{
  CubitPoint *pt0 = facet->point(vert0);
  CubitPoint *pt1 = facet->point(vert1);

  // the edge vector

  CubitVector e0 ( pt1->x() - pt0->x(),
                   pt1->y() - pt0->y(),
                   pt1->z() - pt0->z() );
  e0.normalize();
  
  // vector from vert0 to point

  CubitVector v0 ( pt_on_plane.x() - pt0->x(),
                   pt_on_plane.y() - pt0->y(),
                   pt_on_plane.z() - pt0->z() );
  
  // project to edge

  double len = v0%e0;
  close_point.x ( pt0->x() + e0.x() * len );
  close_point.y ( pt0->y() + e0.y() * len );
  close_point.z ( pt0->z() + e0.z() * len );

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: eval_point
//
//Member Type:  PUBLIC
//Descriptoin:  Evaluate the location and normal of a set of area coordinates  
//              on a facet. 
//===========================================================================
CubitStatus FacetEvalTool::eval_point( CubitFacet *facet, 
                                       int vertex_id, 
                                       CubitVector &close_point )
{
  close_point.x (facet->point(vertex_id)->x());
  close_point.y (facet->point(vertex_id)->y());
  close_point.z (facet->point(vertex_id)->z());

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: eval_facet_normal
//
//Member Type:  PUBLIC
//Descriptoin:  Evaluate the normal of the facet (use the interpOrder)
//              return normalized normal 
//===========================================================================
CubitStatus FacetEvalTool::eval_facet_normal( CubitFacet *facet,
                                              CubitVector &areacoord,
                                              CubitVector &normal )
{
  switch(facet->eval_order()) {
  case 0:
    normal = facet->normal();
    break;
  case 1: case 2: case 3:
    {
      CubitVector norm0 = facet->point(0)->normal( facet );
      CubitVector norm1 = facet->point(1)->normal( facet );
      CubitVector norm2 = facet->point(2)->normal( facet );
      normal.x( areacoord.x() * norm0.x() +
                areacoord.y() * norm1.x() +
                areacoord.z() * norm2.x() );
      normal.y( areacoord.x() * norm0.y() +
                areacoord.y() * norm1.y() +
                areacoord.z() * norm2.y() );
      normal.z( areacoord.x() * norm0.z() +
                areacoord.y() * norm1.z() +
                areacoord.z() * norm2.z() );
      normal.normalize();
    }
    break;
  case 4:
    if(facet->is_flat())
    {
      normal = facet->normal();
    }
    else
    {
      eval_bezier_patch_normal(facet, areacoord, normal);

      // check for reasonableness of the normal

#if 0
      if (DEBUG_FLAG(110))
      {
        CubitVector norm0 = facet->point(0)->normal( facet );
        CubitVector norm1 = facet->point(1)->normal( facet );
        CubitVector norm2 = facet->point(2)->normal( facet );
        CubitVector lin_normal;
        lin_normal.x( areacoord.x() * norm0.x() +
                  areacoord.y() * norm1.x() +
                  areacoord.z() * norm2.x() );
        lin_normal.y( areacoord.x() * norm0.y() +
                  areacoord.y() * norm1.y() +
                  areacoord.z() * norm2.y() );
        lin_normal.z( areacoord.x() * norm0.z() +
                  areacoord.y() * norm1.z() +
                  areacoord.z() * norm2.z() );
        lin_normal.normalize();

        PRINT_INFO("(facet %4d) ac=%7.5lf %7.5lf %7.5lf\n", 
          facet->id(), areacoord.x(), areacoord.y(), areacoord.z());
        PRINT_INFO("            bn=%7.5lf %7.5lf %7.5lf\n",
          normal.x(), normal.y(), normal.z());
        PRINT_INFO("            ln=%7.5lf %7.5lf %7.5lf\n",
          lin_normal.x(), lin_normal.y(), lin_normal.z());

        const double tol = 1e-2;
        if (fabs(lin_normal.x() - normal.x()) > tol ||
            fabs(lin_normal.y() - normal.y()) > tol ||
            fabs(lin_normal.z() - normal.z()) > tol)
        {
          int mydebug = 0;
          if (mydebug)
          {
            CubitVector pt;
            eval_bezier_patch(facet, areacoord, pt);
            dcolor(CUBIT_GREEN);
            dray(pt,normal,1.0);
            dcolor(CUBIT_RED);
            dray(pt,lin_normal,1.0);
            dview();
          }
          
          PRINT_INFO("^=============^==============^=============^=============^\n");
        }
      }
#endif
    }

    break;
  }

  return CUBIT_SUCCESS;
}


//===========================================================================
//Function Name: eval_quadratic
//
//Member Type:  PRIVATE
//Descriptoin:  Evaluate the point based on a quadratic approximation 
//===========================================================================
CubitStatus FacetEvalTool::eval_quadratic( CubitFacet *facet, 
                                           int pt_idx, 
                                           CubitVector &eval_pt,
                                           CubitVector &qpoint,
                                           CubitVector &qnorm )
{
  // interpolate a point on a circle that is defined by two points and
  // two normals.  The first pont is a vertex on the facet, the second
  // is a point on the opposite edge.  The point to be interpolated lies
  // somewhere between the two

  CubitVector normA = facet->point(pt_idx)->normal( facet );
  CubitVector ptA = facet->point(pt_idx)->coordinates();
  int idx0 = -1, idx1 = -1;
  switch(pt_idx) {
  case 0:
    idx0 = 1;
    idx1 = 2;
    break;
  case 1:
    idx0 = 2;
    idx1 = 0;
    break;
  case 2:
    idx0 = 0;
    idx1 = 1;
    break;
  }
  CubitVector ptB, normB, pt0, pt1, norm0, norm1;
  pt0 = facet->point(idx0)->coordinates();
  pt1 = facet->point(idx1)->coordinates();
  norm0 = facet->point(idx0)->normal( facet );
  norm1 = facet->point(idx1)->normal( facet ); 
  on_circle( pt0, norm0, pt1, norm1,
             eval_pt, ptB, normB );
  on_circle( ptA, normA, ptB, normB,
             eval_pt, qpoint, qnorm );

  return CUBIT_SUCCESS;

}

//===========================================================================
//Function Name: on_circle
//
//Member Type:  PRIVATE
//Description:  compute a projection of a point onto a circle.  The circle
//              is defined by two points and two normals.  Return the 
//              projected point and its normal 
//===========================================================================
void FacetEvalTool::on_circle( CubitVector &ptA,
                               CubitVector &normA,
                               CubitVector &ptB,
                               CubitVector &normB,
                               CubitVector &eval_pt,
                               CubitVector &pt_on_circle,
                               CubitVector &normal )
{
  // angle between the normals
  
  double cosang = normA%normB;

  // check for flat surfaces - project to the segment

  if (cosang >= 0.99999 || cosang <= -0.99999) {
    CubitVector vAB = ptB - ptA;
    vAB.normalize();
    CubitVector vAeval_pt = eval_pt - ptA;
    double len = vAB%vAeval_pt;
    pt_on_circle = ptA + len * vAB;
    if (cosang <= -0.99999) {  // this is bad! (facet spans 180 degrees)
      normal = normA;
    }
    else {
      normal = 0.5e0 * (normA + normB);
    }
  }
  else {

    // curved surface
    // define a common plane at eval_pt

    CubitVector pnorm = normA*normB;
    pnorm.normalize();
    double pcoefd = -pnorm%eval_pt;

    // project everything to common plane

    double pdist = pnorm%ptA + pcoefd;
    CubitVector pptA = ptA - pnorm * pdist;
    pdist = pnorm%ptB + pcoefd;
    CubitVector pptB = ptB - pnorm * pdist; 
    
    double angle = acos(cosang);
    double dist = pptA.distance_between(pptB);

    // kradius is the radius of curvature
    // center is the center of curvature
    // centerA and centerB should be the same within tol

    double kradius = dist / (2.0e0 * sin( angle * 0.5e0 ));
    CubitVector centerA = pptA - kradius * normA;
    CubitVector centerB = pptB - kradius * normB;
    CubitVector center = (centerA + centerB) * 0.5e0;

    normal = eval_pt - center;
    normal.normalize();
    pt_on_circle = center + kradius * normal;
  }
}

//===========================================================================
//Function Name: eval_quadric
//
//Member Type:  PRIVATE
//Description:  evaluate a point on an interpolated quadric surface  
//===========================================================================
void FacetEvalTool::eval_quadric( CubitFacet *facet,
                                  int pt_index,
                                  CubitVector &eval_pt,
                                  CubitVector &qpt )
{
  // transform the point to the local system

  CubitPoint *point = facet->point(pt_index);
  CubitVector loc_eval_pt;
  point->transform_to_local( eval_pt, loc_eval_pt );
//  point->transform_to_local( point->coordinates(), loc_pt );
  
  // interpolate a "z" value in the local coordinate system

  double *coef = point->coef_vector();
  loc_eval_pt.z( coef[0] * loc_eval_pt.x() +
                 coef[1] * loc_eval_pt.y() +
                 coef[2] * sqr(loc_eval_pt.x()) +
                 coef[3] * loc_eval_pt.x() * loc_eval_pt.y() +
                 coef[4] * sqr(loc_eval_pt.y()) );
  
  
//  loc_eval_pt.z( point->z() -
//                 coef[0] * loc_eval_pt.x() -
//                 coef[1] * loc_eval_pt.y() +
//                 coef[2] * sqr(loc_eval_pt.x()) +
//                 coef[3] * loc_eval_pt.x() * loc_eval_pt.y() +
//                 coef[4] * sqr(loc_eval_pt.y()) );

  // transform back to global system

  point->transform_to_global( loc_eval_pt, qpt );

}

//===========================================================================
//Function Name: project_to_facet_plane
//
//Member Type:  PUBLIC
//Descriptoin:  Project a point to the plane of a facet 
//===========================================================================
void FacetEvalTool::project_to_facet_plane(
  CubitFacet *facet,
  const CubitVector &pt,
  CubitVector &point_on_plane,
  double &dist_to_plane)
{
  CubitVector normal = facet->normal();
  normal.normalize();
  CubitPoint *facPoint = facet->point(0);
  double coefd = -(normal.x() * facPoint->x() +
                   normal.y() * facPoint->y() +
                   normal.z() * facPoint->z());
  double dist = normal.x() * pt.x() +
                normal.y() * pt.y() +
                normal.z() * pt.z() + coefd;
  dist_to_plane = fabs(dist);

  point_on_plane.set( pt.x() - normal.x() * dist,
                      pt.y() - normal.y() * dist,
                      pt.z() - normal.z() * dist );

}

//===========================================================================
//Function Name: facet_area_coordinate
//
//Member Type:  PUBLIC
//Descriptoin:  Determine the area coordinates of a point on the plane 
//              of a facet 
//===========================================================================
void FacetEvalTool::facet_area_coordinate( 
  CubitFacet *facet,
  const CubitVector &pt_on_plane,
  CubitVector &areacoord )
{
  double area2;
  CubitVector normal = facet->normal();
  CubitPoint *pt0 = facet->point(0);
  CubitPoint *pt1 = facet->point(1);
  CubitPoint *pt2 = facet->point(2);
  double tol = GEOMETRY_RESABS * 1.e-5;
  CubitVector v1( pt1->x() - pt0->x(),
                  pt1->y() - pt0->y(),
                  pt1->z() - pt0->z());//(*p1-*p0);
  CubitVector v2( pt2->x() - pt0->x(),
                  pt2->y() - pt0->y(),
                  pt2->z() - pt0->z());// = (*p2-*p0);
  area2 = (v1*v2).length_squared();
  if(area2 < 100 * tol){
      tol = .01 * area2;
  }
  CubitVector absnorm( fabs(normal.x()), fabs(normal.y()), fabs(normal.z()) );
  
  // project to the closest coordinate plane so we only have to do this in 2D

  if (absnorm.x() >= absnorm.y() && absnorm.x() >= absnorm.z()) {
    area2 = determ3(pt0->y(), pt0->z(),
                    pt1->y(), pt1->z(),
                    pt2->y(), pt2->z());
    if (fabs(area2) < tol) {
      areacoord.set( -CUBIT_DBL_MAX, -CUBIT_DBL_MAX, -CUBIT_DBL_MAX );
    }
    else if ( pt_on_plane.within_tolerance( pt0->coordinates(), GEOMETRY_RESABS ) )
    {
        areacoord.set( 1.0, 0.0, 0.0 );
    }
    else if ( pt_on_plane.within_tolerance( pt1->coordinates(), GEOMETRY_RESABS ) )
    {
        areacoord.set( 0.0, 1.0, 0.0 );
    }
    else if ( pt_on_plane.within_tolerance( pt2->coordinates(), GEOMETRY_RESABS ) )
    {
        areacoord.set( 0.0, 0.0, 1.0 );
    }
    else {
      areacoord.x( 
        determ3( pt_on_plane.y(), pt_on_plane.z(),
                 pt1->y(), pt1->z(), pt2->y(), pt2->z() ) / area2 );
      areacoord.y( 
        determ3( pt0->y(), pt0->z(),
                 pt_on_plane.y(), pt_on_plane.z(),
                 pt2->y(), pt2->z() ) / area2 );
      areacoord.z( 
        determ3( pt0->y(), pt0->z(), pt1->y(), pt1->z(),
                 pt_on_plane.y(), pt_on_plane.z() ) / area2 );
    }
  }
  else if(absnorm.y() >= absnorm.x() && absnorm.y() >= absnorm.z()) {
    area2 = determ3(pt0->x(), pt0->z(),
                    pt1->x(), pt1->z(),
                    pt2->x(), pt2->z());
    if (fabs(area2) < tol) {
      areacoord.set( -CUBIT_DBL_MAX, -CUBIT_DBL_MAX, -CUBIT_DBL_MAX );
    }
    else if ( pt_on_plane.within_tolerance( pt0->coordinates(), GEOMETRY_RESABS ) )
    {
        areacoord.set( 1.0, 0.0, 0.0 );
    }
    else if ( pt_on_plane.within_tolerance( pt1->coordinates(), GEOMETRY_RESABS ) )
    {
        areacoord.set( 0.0, 1.0, 0.0 );
    }
    else if ( pt_on_plane.within_tolerance( pt2->coordinates(), GEOMETRY_RESABS ) )
    {
        areacoord.set( 0.0, 0.0, 1.0 );
    }
    else {
      areacoord.x( 
        determ3( pt_on_plane.x(), pt_on_plane.z(),
                 pt1->x(), pt1->z(), pt2->x(), pt2->z() ) / area2 );
      areacoord.y( 
        determ3( pt0->x(), pt0->z(),
                 pt_on_plane.x(), pt_on_plane.z(),
                 pt2->x(), pt2->z() ) / area2 );
      areacoord.z( 
        determ3( pt0->x(), pt0->z(), pt1->x(), pt1->z(),
                 pt_on_plane.x(), pt_on_plane.z() ) / area2 );
    }
  }
  else { 
    area2 = determ3(pt0->x(), pt0->y(),
                    pt1->x(), pt1->y(),
                    pt2->x(), pt2->y());
    if (fabs(area2) < tol) {
      areacoord.set( -CUBIT_DBL_MAX, -CUBIT_DBL_MAX, -CUBIT_DBL_MAX );
    }
    else if ( pt_on_plane.within_tolerance( pt0->coordinates(), GEOMETRY_RESABS ) )
    {
        areacoord.set( 1.0, 0.0, 0.0 );
    }
    else if ( pt_on_plane.within_tolerance( pt1->coordinates(), GEOMETRY_RESABS ) )
    {
        areacoord.set( 0.0, 1.0, 0.0 );
    }
    else if ( pt_on_plane.within_tolerance( pt2->coordinates(), GEOMETRY_RESABS ) )
    {
        areacoord.set( 0.0, 0.0, 1.0 );
    }
    else {
      areacoord.x( 
        determ3( pt_on_plane.x(), pt_on_plane.y(),
                 pt1->x(), pt1->y(), pt2->x(), pt2->y() ) / area2 );
      areacoord.y( 
        determ3( pt0->x(), pt0->y(),
                 pt_on_plane.x(), pt_on_plane.y(),
                 pt2->x(), pt2->y() ) / area2 );
      areacoord.z( 
        determ3( pt0->x(), pt0->y(), pt1->x(), pt1->y(),
                 pt_on_plane.x(), pt_on_plane.y() ) / area2 );
    }
  }
}

//===========================================================================
//Function Name: is_outside
//
//Member Type:  PRIVATE
//Descriptoin:  Determines if the areacoord is actually outside the range
//              of the surface facets.
//===========================================================================
CubitBoolean FacetEvalTool::is_outside(
  CubitFacet *facet, 
  CubitVector &areacoord)
{
  if (areacoord.x() < -GEOMETRY_RESABS) {
    if (NULL == facet->shared_facet_on_surf( facet->point(1), 
                                             facet->point(2), 
                                             facet->tool_id() )) {
      return CUBIT_TRUE;
    }
  }
  if (areacoord.y() < -GEOMETRY_RESABS) {
    if (NULL == facet->shared_facet_on_surf( facet->point(2), 
                                             facet->point(0), 
                                             facet->tool_id() )) {
      return CUBIT_TRUE;
    }
  }
  if (areacoord.z() < -GEOMETRY_RESABS) {
    if (NULL == facet->shared_facet_on_surf( facet->point(0), 
                                             facet->point(1), 
                                             facet->tool_id() )) {
      return CUBIT_TRUE;
    }
  }
  return CUBIT_FALSE;
}

//===========================================================================
//Function Name: closest_point_trimmed
//
//Member Type:  PUBLIC
//Descriptoin:  Finds the closest point from the vector (this_point) to the
//              set of facets that lies on the set of facets.  If the point
//              lies outside this set, the closest point will be on the edge
//              of the closest facet.  
//              This function also determines if the point is outside or
//              inside the current set of facets.  You should call
//              a bounding box test first before this...
//===========================================================================
CubitStatus FacetEvalTool::closest_point_trimmed(
  CubitVector &this_point,
  CubitVector *closest_point_ptr,
  CubitBoolean &lies_outside,
  CubitVector *normal_ptr)
{

  CubitBoolean trim = CUBIT_TRUE;
  return project_to_facets ( this_point, trim, &lies_outside,
                             closest_point_ptr, normal_ptr );
}

//===========================================================================
//Function Name: destroy_facets
//
//Member Type:  PRIVATE
//Descriptoin:  Deletes the points and facets.
//===========================================================================
void FacetEvalTool::destroy_facets()
{
  int i, j;
  CubitPoint *point;
  CubitFacet *facet;
  CubitFacetEdge *edge;
  myFacetList.last();
  CubitFacetEdgeData *cfe_data = NULL;
  for( i = myFacetList.size(); i > 0; i-- )
  {
    facet = myFacetList.pop();
    for (j = 0; j<3; j++)
    {
      point = facet->point( j );
      point->remove_facet( facet );
      edge = facet->edge( j );
      cfe_data = CAST_TO( edge, CubitFacetEdgeData );
      if (cfe_data)
        cfe_data->remove_facet( facet );
    }
    delete facet;
  }
  for( i = myEdgeList.size(); i > 0; i-- )
  {
    edge = myEdgeList.pop();
    if (edge && edge->num_adj_facets() == 0)
    {
      delete edge;
    }
  }
  for( i = myPointList.size(); i > 0; i-- )
  {
    point = myPointList.pop();
    if (point && point->num_adj_facets() == 0)
    {
      delete point;
    }
  }
}
  
//===========================================================================
//Function Name: get_points_from_facets
//
//Member Type:  PRIVATE
//Descriptoin:  Gets the set of points contained in the list of facets.  
//              Populates the point_list with those points.
//===========================================================================
CubitStatus FacetEvalTool::get_points_from_facets(
  DLIList<CubitFacet*> &facet_list,
  DLIList<CubitPoint*> &point_list )
{
  int i, j;

  DLIList<CubitPoint*> all_points;
  
  for ( i = 0; i < facet_list.size(); i++)
  {
    CubitFacet* facet = facet_list.get_and_step();    
    facet->points(all_points);
  }
  for ( i = 0; i < all_points.size(); i++)
  {
    all_points.get_and_step()->marked( CUBIT_FALSE );
  }

  for ( j = 0; j < all_points.size(); j++)
  {
    CubitPoint* point = all_points.get_and_step();
    if (!point->marked())
    {
      point->marked(CUBIT_TRUE);
      point_list.append(point);
    }
  }

  // unmark the found points
  
  for (i = 0; i < point_list.size(); i++)
  {
    point_list.get_and_step()->marked(CUBIT_FALSE);
  }  
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: get_loops_from_facets
//
//Member Type:  PRIVATE
//Descriptoin:  generate an ordered list of facetedge lists representing the
//              boundary loops for this set of facets
//===========================================================================
CubitStatus FacetEvalTool::get_loops_from_facets( 
  DLIList<CubitFacetEdge*> &all_edge_list,    // all the edges on the facets
  DLIList<DLIList<CubitFacetEdge*>*> &loop_list )  // return a list of edge lists
{
  int i;
  int mydebug = 0;
  CubitFacetEdge* edge, *startedge;
  CubitFacet *facet, *nextfacet;
  CubitBoolean found;

  if (mydebug)
  {
    GfxDebug::clear();
    for (i = 0; i< myFacetList.size(); i++)
      myFacetList.get_and_step()->debug_draw();
    GfxDebug::flush();
    GfxDebug::mouse_xforms();
  } 

  for ( i = 0; i < all_edge_list.size(); i++)
  {
    startedge = all_edge_list.get_and_step();
    if (startedge->get_flag() == 0) {
      startedge->set_flag( 1 );

      // Find an edge without a neighboring facet or a facet that
      // is not on the current surface to start a loop

      if (startedge->num_adj_facets_on_surf(toolID) == 1)
      {     
        
        // Start a new loop

        DLIList<CubitFacetEdge*> *edge_list = new DLIList<CubitFacetEdge*>;
        loop_list.append( edge_list );
        edge_list->append( startedge );
        if (mydebug) debug_draw_edge(startedge);

        // find the next ccw edge on the loop.  Edges are placed on the
        // list based on ccw order wrt the orientation of the facets
        // (ie. same orientation as facets)

        edge = startedge;
        facet = edge->adj_facet_on_surf(toolID);       
        do {
          found = CUBIT_FALSE;
          do {
            edge = facet->next_edge( edge );
            assert(edge != NULL);
            edge->set_flag( 1 );
            nextfacet = edge->other_facet_on_surf(facet);
            if (nextfacet == NULL) {
              if (edge != startedge) {
                edge_list->append( edge );
                if (mydebug) 
                {
                  debug_draw_edge(edge);
                  GfxDebug::mouse_xforms();
                }

              }
              found = CUBIT_TRUE;
            }
            else {
              facet = nextfacet;
            }
          } while (!found);
        } while (edge != startedge);
      }
    }
  }

  // reset the flags

  for ( i = 0; i < all_edge_list.size(); i++) {
    all_edge_list.get_and_step()->set_flag( 0 );
  }
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: get_edges_from_facets
//
//Member Type:  PRIVATE
//Descriptoin:  Populates the edge list from the list of facets  
//===========================================================================
CubitStatus FacetEvalTool::get_edges_from_facets(
  DLIList<CubitFacet*> &facet_list,
  DLIList<CubitFacetEdge*> &edge_list )
{
  int i, j, k;
  CubitPoint *p0, *p1;
  CubitFacet *facet, *adj_facet;
  CubitFacetEdge *edge;
  for ( i = 0; i < facet_list.size(); i++) {
    facet = facet_list.get_and_step();
    for (j=0; j<3; j++) {
      if (!(edge = facet->edge(j))) 
      {
        facet->get_edge_pts( j, p0, p1 );
        edge = NULL;
        k = -1;
        adj_facet = facet->shared_facet_on_surf( p0, p1, toolID );
        if (adj_facet) {
          edge = adj_facet->edge_from_pts(p0, p1, k);
        }
        if (!edge) {
          edge = (CubitFacetEdge *) 
            new CubitFacetEdgeData( p0, p1, facet, adj_facet, j, k );
          
        }
      }
      edge->set_flag( 0 );
    }
  }
    // create a unique list of edges

  for ( i = 0; i < facet_list.size(); i++)
  {
    facet = facet_list.get_and_step();
    for (j=0; j<3; j++)
    {
      edge = facet->edge(j);
      if (0 == edge->get_flag())
      {
        edge->set_flag( 1 );
        edge_list.append( edge );
      }
    }
  }

  // reset the flags on the edges

  for ( i = 0; i < edge_list.size(); i++) 
  {
    edge = edge_list.get_and_step();
    edge->set_flag( 0 );
  }
  return CUBIT_SUCCESS;
}


//===========================================================================
//Function Name: check_faceting
//
//Member Type:  PRIVATE
//Descriptoin:  check the edge/face orientations  
//===========================================================================
void FacetEvalTool::check_faceting()
{
  int ii;
  for (ii = 0; ii< myFacetList.size(); ii++)
      myFacetList.get_and_step()->debug_draw( CUBIT_YELLOW );
  GfxDebug::flush();
  GfxDebug::mouse_xforms();

  CubitFacet *facet = myFacetList.get_and_step();
  CubitVector snorm = facet->normal();
  snorm.normalize();
  for (ii=1; ii<myFacetList.size(); ii++)
  {
    facet = myFacetList.get_and_step();
    CubitVector norm = facet->normal();
    norm.normalize();
    if (norm % snorm < 0.8)
    {
      facet->debug_draw( CUBIT_RED );
    }
    else
    {
      facet->debug_draw( CUBIT_BLUE );
    }
  }
}


//===========================================================================
//Function Name: draw_facets
//
//Member Type:  PRIVATE
//Descriptoin:  draw the facets for debug 
//===========================================================================
void FacetEvalTool::debug_draw_facets( int color )
{
  int ii;
  if ( color == -1 )
    color = CUBIT_YELLOW;
  CubitBoolean flush_it = CUBIT_FALSE;
  for ( ii = myFacetList.size(); ii > 0; ii-- )
  {
    myFacetList.get_and_step()->debug_draw(color, flush_it);
  }
  GfxDebug::flush();
}

//===========================================================================
//Function Name: draw_vec
//
//Member Type:  PRIVATE
//Descriptoin:  draw a single point 
//===========================================================================
void FacetEvalTool::debug_draw_vec( CubitVector &vec, int color )
{
  if ( color == -1 )
    color = CUBIT_YELLOW;
  GfxDebug::draw_point(vec, color);
  GfxDebug::flush();
}

//===========================================================================
//Function Name: draw_facet_normals
//
//Member Type:  PRIVATE
//Descriptoin:  the the normal at the centroid of the facets 
//===========================================================================
void FacetEvalTool::debug_draw_facet_normals(int color)
{
  int ii,jj;
  if ( color == -1 )
    color = CUBIT_RED;
  for ( ii = myFacetList.size(); ii > 0; ii-- )
  {
    CubitFacet *facet = myFacetList.get_and_step();
    CubitVector center(0.0e0, 0.0e0, 0.0e0);
    for (jj = 0; jj < 3; jj++) 
    {
      center += facet->point(jj)->coordinates();
    }
    center /= 3.0;
    CubitVector normal = facet->normal();
    double mag = sqrt(sqr(normal.x()) + sqr(normal.y()) + sqr(normal.z()));
    mag = sqrt(mag);
    normal.normalize();
    CubitVector end = center + normal*mag;
    GfxDebug::draw_line(center.x(), center.y(), center.z(),
                        end.x(), end.y(), end.z(), color);
  }
  GfxDebug::flush();
}

//===========================================================================
//Function Name: draw_point_normals
//
//Member Type:  PRIVATE
//Descriptoin:  the the normal at the points 
//===========================================================================
void FacetEvalTool::debug_draw_point_normals(int color)
{
  int ii;
  if ( color == -1 )
    color = CUBIT_RED;
  double len = 0.1 * sqrt(sqr(myBBox->x_range()) + 
	                        sqr(myBBox->y_range()) + 
	                        sqr(myBBox->z_range()));
  for ( ii = myPointList.size(); ii > 0; ii-- )
  {
    CubitPoint *point = myPointList.get_and_step();
    CubitVector normal = point->normal();
    CubitVector end = point->coordinates() + normal*len;
    GfxDebug::draw_line(point->x(), point->y(), point->z(),
                        end.x(), end.y(), end.z(), color);
  }
  GfxDebug::flush();
}

//===========================================================================
//Function Name: draw_edges
//
//Member Type:  PRIVATE
//Descriptoin:  draw the facet edges 
//===========================================================================
void FacetEvalTool::debug_draw_edges(int color)
{
  int ii;
  if ( color == -1 )
    color = CUBIT_BLUE;
  for ( ii = myEdgeList.size(); ii > 0; ii-- )
  {
    CubitFacetEdge *edge = myEdgeList.get_and_step();
    CubitPoint *begin_point = edge->point(0);
    CubitPoint *end_point = edge->point(1);
    GfxDebug::draw_line(begin_point->x(), 
                        begin_point->y(), 
                        begin_point->z(),
                        end_point->x(), 
                        end_point->y(), 
                        end_point->z(), 
                        color);
  }
  GfxDebug::flush();
}

//===========================================================================
//Function Name: draw_edge
//
//Member Type:  PRIVATE
//Descriptoin:  draw the facet edge 
//===========================================================================
void FacetEvalTool::debug_draw_edge(CubitFacetEdge *edge, int color)
{
  if ( color == -1 ) {
    color = CUBIT_YELLOW;
  }
  CubitPoint *begin_point = edge->point(0);
  CubitPoint *end_point = edge->point(1);
  GfxDebug::draw_line(begin_point->x(), 
                      begin_point->y(), 
                      begin_point->z(),
                      end_point->x(), 
                      end_point->y(), 
                      end_point->z(), 
                      color);
  GfxDebug::flush();
}

//===========================================================================
//Function Name: draw_bezier_edges
//
//Member Type:  PRIVATE
//Descriptoin:  draw the control polygons from the bezier control points on
//               the edges 
//===========================================================================
void FacetEvalTool::debug_draw_bezier_edges(int color)
{
  int ii, i;
  CubitVector ctrl_pts[5], begin, end;
  if ( color == -1 )
    color = CUBIT_RED;
  for ( ii = myEdgeList.size(); ii > 0; ii-- )
  {
    CubitFacetEdge *edge = myEdgeList.get_and_step();
    edge->control_points( ctrl_pts );
    for (i=1; i<5; i++) {
      begin = ctrl_pts[i-1];
      end = ctrl_pts[i];
      GfxDebug::draw_line(begin.x(), begin.y(), begin.z(),
                          end.x(), end.y(), end.z(), color);
    }
  }
  GfxDebug::flush();
}

//===========================================================================
//Function Name: draw_bezier_facets
//
//Member Type:  PRIVATE
//Descriptoin:  draw the control polygons from the bezier control points on
//               the facets 
//===========================================================================
void FacetEvalTool::debug_draw_bezier_facets(int color)
{
  int ii;
  if ( color == -1 )
    color = CUBIT_WHITE;
  for ( ii = myFacetList.size(); ii > 0; ii-- )
  {
    debug_draw_bezier_facet( myFacetList.get_and_step(), color );
  }
  
}

//===========================================================================
//Function Name: draw_bezier_facet
//
//Member Type:  PRIVATE
//Descriptoin:  draw the control polygons from the bezier control points on
//               a facet 
//===========================================================================
void FacetEvalTool::debug_draw_bezier_facet(CubitFacet *facet, int color)
{
  CubitVector areacoord( 0.3333333333333333333, 
                         0.3333333333333333333, 
                         0.3333333333333333333 );

  // interpolate internal control points

  CubitVector gctrl_pts[6];
  facet->get_control_points( gctrl_pts );
  CubitVector P_facet[3];

  //2,1,1
  P_facet[0] = (1.0e0 / (areacoord.y() + areacoord.z())) *
               (areacoord.y() * gctrl_pts[3] +
                areacoord.z() * gctrl_pts[4]);
  //1,2,1
  P_facet[1] = (1.0e0 / (areacoord.x() + areacoord.z())) *
               (areacoord.x() * gctrl_pts[0] +
                areacoord.z() * gctrl_pts[5]);
  //1,1,2
  P_facet[2] = (1.0e0 / (areacoord.x() + areacoord.y())) *
               (areacoord.x() * gctrl_pts[1] +
                areacoord.y() * gctrl_pts[2]);

  CubitVector cp0[5], cp1[5], cp2[5];

  facet->edge(0)->control_points(facet,cp0);
  facet->edge(1)->control_points(facet,cp1);
  facet->edge(2)->control_points(facet,cp2);

  color = CUBIT_WHITE;
  debug_draw_line(cp0[1], cp2[3], color);
  debug_draw_line(cp0[2], P_facet[1], color);
  debug_draw_line(P_facet[1], cp2[2], color);
  debug_draw_line(cp0[3], P_facet[2], color);
  debug_draw_line(P_facet[2], P_facet[0], color);
  debug_draw_line(P_facet[0], cp2[1], color);

  color = CUBIT_GREEN;
  debug_draw_line(cp1[1], P_facet[2], color);
  debug_draw_line(P_facet[2], P_facet[1], color);
  debug_draw_line(P_facet[1], cp2[3], color);
  debug_draw_line(cp1[2], P_facet[0], color);
  debug_draw_line(P_facet[0], cp2[2], color);
  debug_draw_line(cp1[3], cp2[1], color);

  color = CUBIT_BLUE;
  debug_draw_line(cp0[1], P_facet[1], color);
  debug_draw_line(P_facet[1], P_facet[0], color);
  debug_draw_line(P_facet[0], cp1[3], color);
  debug_draw_line(cp0[2], P_facet[2], color);
  debug_draw_line(P_facet[2], cp1[2], color);
  debug_draw_line(cp0[3], cp1[1], color);
}

//===========================================================================
//Function Name: draw_eval_bezier_facet
//
//Member Type:  PRIVATE
//Descriptoin:  draw points on the evaluated bezier patch 
//===========================================================================
void FacetEvalTool::debug_draw_eval_bezier_facet( CubitFacet *facet )
{
  CubitVector areacoord, pt, loc;
  double u, v, w;
  CubitBoolean outside;
#if 0 
  for (int i=0; i<=10; i++) {
    v = w = 0.5 * (double)i/10.0;
    u = 1.0 - v - w;
    areacoord.set(u, v, w);
    eval_facet(facet, pt, areacoord, loc, outside);
    draw_location(loc);
    v = 0.25 * (double)i/10.0;
    w = 0.75 * (double)i/10.0;
    u = 1.0 - v - w;
    areacoord.set(u, v, w);
    eval_facet(facet, pt, areacoord, loc, outside);
    draw_location(loc);
    w = 0.25 * (double)i/10.0;
    v = 0.75 * (double)i/10.0;
    u = 1.0 - v - w;
    areacoord.set(u, v, w);
    eval_facet(facet, pt, areacoord, loc, outside);
    debug_draw_location(loc);
  }
#endif
  for (int j=0; j<=10; j++) {
    for (int i=0; i<=20; i++) {
      u = ((double)j/10.0) * (double)i/20.0;
      w = (1.0 - ((double)j/10.0)) * (double)i/20.0;
      v = 1.0 - u - w;
      areacoord.set(u, v, w);
  eval_facet(facet, pt, areacoord, loc, outside);
      debug_draw_location(loc);
    }
  }
 
#if 0
  for (i=0; i<=10; i++) {
    u = v = 0.5 * (double)i/10.0;
    w = 1.0 - u - v;
    areacoord.set(u, v, w);
    eval_facet(facet, pt, areacoord, loc, outside);
    draw_location(loc);
    u = 0.25 * (double)i/10.0;
    v = 0.75 * (double)i/10.0;
    w = 1.0 - u - v;
    areacoord.set(u, v, w);
    eval_facet(facet, pt, areacoord, loc, outside);
    draw_location(loc);
    v = 0.25 * (double)i/10.0;
    u = 0.75 * (double)i/10.0;
    w = 1.0 - u - v;
    areacoord.set(u, v, w);
    eval_facet(facet, pt, areacoord, loc, outside);
    draw_location(loc);
  }
#endif

}

//===========================================================================
//Function Name: draw_line
//
//Member Type:  PRIVATE 
//===========================================================================
void FacetEvalTool::debug_draw_line(CubitVector &begin, CubitVector &end, int color)
{
  GfxDebug::draw_line(begin.x(), begin.y(), begin.z(),
                      end.x(), end.y(), end.z(), color);
  GfxDebug::flush();
}

//===========================================================================
//Function Name: draw_location
//
//Member Type:  PRIVATE 
//===========================================================================
void FacetEvalTool::debug_draw_location(CubitVector &loc, int color )
{
  if ( color == -1 )
    color = CUBIT_YELLOW;
  GfxDebug::draw_point(loc, color);
  GfxDebug::flush();
}


//===========================================================================
//Function Name: write_loops
//
//Member Type:  PRIVATE 
//===========================================================================
void FacetEvalTool::write_loops()
{
  int ii, jj;
  for (ii=0; ii<myLoopList.size(); ii++)
  {
    DLIList<CubitFacetEdge*> *loop = myLoopList.get_and_step();
    PRINT_INFO("======= Loop %d =========\n", ii);
    for (jj=0; jj<loop->size(); jj++)
    {
      CubitFacetEdge *edge = loop->get_and_step();
      CubitPoint *point0 = edge->point( 0 );
      CubitPoint *point1 = edge->point( 1 );
      PRINT_INFO("  (%d) %f, %f, %f   (%d) %f, %f, %f\n",
        point0->id(), point0->x(), point0->y(), point0->z(),
        point1->id(), point1->x(), point1->y(), point1->z());
    }
  }
}

//===========================================================================
//Function Name: transform_control_points
//
//Member Type:  PUBLIC
//===========================================================================
void FacetEvalTool::transform_control_points( CubitTransformMatrix &tfmat )
{
  if (interpOrder != 4)
    return;

  CubitVector control_points[6];
  CubitFacet *facet;
  int ii, jj;
  for (ii=0; ii<myFacetList.size(); ii++)
  {
    facet = myFacetList.get_and_step();
    facet->get_control_points( control_points );
    for (jj=0; jj<6; jj++)
    {
      control_points[jj] = tfmat * control_points[jj];
    }
    facet->set_control_points( control_points );
  }

  CubitFacetEdge *edge;
  for (ii=0; ii<myEdgeList.size(); ii++)
  { 
    edge = myEdgeList.get_and_step();
    if (!edge->get_flag())
    {
      edge->set_flag(1);   
      edge->control_points( control_points );
    
      // assumes the end control points (the vertices) have already 
      // been transformed

      control_points[0] = tfmat * control_points[1];
      control_points[1] = tfmat * control_points[2];
      control_points[2] = tfmat * control_points[3];
      edge->control_points( control_points, 4 );
    }
  }
}

//===========================================================================
//Function Name: area
//
//Member Type:  PUBLIC
//===========================================================================
double FacetEvalTool::area()
{
  if (myArea < 0.0)
    calculate_area();
  
  return myArea;
}
//===========================================================================
//Function Name: calcualte_area
//
//Member Type:  Public
//===========================================================================
void FacetEvalTool::calculate_area()
{
    myArea = 0.0;
    int ii;
    CubitFacet *facet;
    for (ii=0; ii<myFacetList.size(); ii++)
    {
      facet = myFacetList.get_and_step();
      myArea += facet->area();
    }
}

//===========================================================================
//Function Name: parameterize
//Description: compute the parameterization of the facetted representation
//Member Type:  PUBLIC
//===========================================================================
CubitStatus FacetEvalTool::parameterize()
{ 

  if (myLoopList.size() != 1)
  {
    PRINT_WARNING("Cannot parameterize surface.  Multiple loops detected\n");
    isParameterized = CUBIT_FALSE;
    return CUBIT_FAILURE;
  }

  // make arrays out of the points and facets

  double *points = new double [3 * myPointList.size()];
  int *facets = new int [3 * myFacetList.size()];
  if (!points || !facets)
  {
    PRINT_ERROR("Could not define parameterization for surface (out of memory)\n");
    return CUBIT_FAILURE;
  }
  int ii;
  CubitPoint *pt;
  myPointList.reset();
  for (ii=0; ii<myPointList.size(); ii++)
  {
    pt = myPointList.get_and_step();
    points[ii*3] = pt->x();
    points[ii*3+1] = pt->y();
    points[ii*3+2] = pt->z();
    pt->marked(ii);
  }
  CubitFacet *facet;
  CubitPoint *pts[3];
  for (ii=0; ii<myFacetList.size(); ii++)
  {
    facet = myFacetList.get_and_step();    
    facet->points( pts[0], pts[1], pts[2] );
    facets[ii*3]   = pts[0]->marked();
    facets[ii*3+1] = pts[1]->marked();
    facets[ii*3+2] = pts[2]->marked();
  }

  // do the parameterization

// comment out for now
// Note to sjowen:  this depends on FacetParamTool and facetParamTool is a ParamTool
// (which is in geom directory). We ned a solution that breaks that dependency.
//  FacetParamTool facetparamtool( myPointList.size(), myFacetList.size(),
//                                 points, facets );
  if(1)//!facetparamtool.flatten())
  {
    PRINT_ERROR("Surface Parameterizer Failed\n");
    isParameterized = CUBIT_FALSE;
  }
  else
  {
    double *sizes = NULL;
    double *uv = NULL;//facetparamtool.get_uvs_sizing( ratio, sizes ); 

    // update the points with uv values

    TDFacetBoundaryPoint *td_fbp;
    CubitBoolean on_internal_boundary;
    myPointList.reset();
    for (ii=0; ii<myPointList.size(); ii++)
    {
      pt = myPointList.get_and_step();
      DLIList <CubitFacet *> facet_list;
      pt->facets_on_surf( toolID, facet_list, on_internal_boundary );
      if (on_internal_boundary)
      {
        td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( pt );
        if (!td_fbp)
        {
          TDFacetBoundaryPoint::add_facet_boundary_point( pt );
          td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( pt );
          td_fbp->add_surf_facets( facet_list );
          td_fbp->set_uv( facet_list.get(), uv[ii*2], uv[ii*2+1] );
        }
        else
        {
          if (td_fbp->set_uv( facet_list.get(),
              uv[ii*2], uv[ii*2+1] ) != CUBIT_SUCCESS)
          {
            td_fbp->add_surf_facets( facet_list );
            td_fbp->set_uv( facet_list.get(), uv[ii*2], uv[ii*2+1] );
          }
        }
      }
      else
      {
        pt->set_uv( uv[ii*2], uv[ii*2+1] );
      }
    }
    isParameterized = CUBIT_TRUE;
    PRINT_INFO("Surface Parameterization succeeded\n");
    delete [] sizes;
    delete [] uv;
  }

  // clean up

  delete [] points;
  delete [] facets;
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: compute_curve_tangent
//
//Member Type:  PRIVATE
//Descriptoin:  compute the tangents to the endpoints of a boundary edge.
//===========================================================================
CubitStatus FacetEvalTool::compute_curve_tangent( 
  CubitFacetEdge *edge,
  double min_dot,
  CubitVector &T0,
  CubitVector &T3 )
{

  CubitPoint *p0 = edge->point( 0 );
  CubitPoint *p1 = edge->point( 1 );
  CubitFacetEdge *prev_edge = next_boundary_edge( edge, p0 );
  if (prev_edge == NULL)  // could be end of a hard line
  {
    T0 = p1->coordinates() - p0->coordinates();  
    T0.normalize();  
  }
  else
  {
    CubitPoint *p2 = prev_edge->other_point( p0 );
    CubitVector e0 = p0->coordinates() - p2->coordinates();
    CubitVector e1 = p1->coordinates() - p0->coordinates();
    e0.normalize();
    e1.normalize();
    if (e0 % e1 >= min_dot)
    {
      T0 = (p0->coordinates() - p2->coordinates()) + 
           (p1->coordinates() - p0->coordinates());
      T0.normalize();
    }
    else
    {
      T0 = e1;
    }
  }
  

  CubitFacetEdge *next_edge = next_boundary_edge( edge, p1 );
  if (next_edge == NULL)  // could be end of a hard line
  {
    T3 = p1->coordinates() - p0->coordinates();
    T3.normalize();
  }
  else
  {
    CubitPoint *p2 = next_edge->other_point( p1 );
    CubitVector e0 = p2->coordinates() - p1->coordinates();
    CubitVector e1 = p1->coordinates() - p0->coordinates();
    e0.normalize();
    e1.normalize();
    if (e0 % e1 >= min_dot)
    {
      T3 = (p2->coordinates() - p1->coordinates()) + 
           (p1->coordinates() - p0->coordinates());
      T3.normalize();
    }
    else
    {
      T3 = e1;
    }
  }
  
  return CUBIT_SUCCESS;
}


//===========================================================================
//Function Name: next_boundary_edge
//
//Member Type:  PRIVATE
//Descriptoin:  given a facet boundary edge and one of its nodes, find the
//              next edge on the same surface  
//===========================================================================
CubitFacetEdge *FacetEvalTool::next_boundary_edge( 
  CubitFacetEdge *this_edge,
  CubitPoint *p0 )
{
  CubitFacetEdge *next_edge = NULL;

  DLIList<CubitFacetEdge*> edge_list;
  p0->edges( edge_list );
  int ii;

  CubitFacetEdge *edge_ptr = NULL;
  for (ii=0; ii<edge_list.size() && next_edge == NULL; ii++)
  {
    edge_ptr = edge_list.get_and_step();
    if (edge_ptr != this_edge)
    {
      if (edge_ptr->num_adj_facets() <= 1)
      {
        next_edge = edge_ptr;
      }
    }
  }

  return next_edge;
}

//===========================================================================
//Function Name: save
//
//Member Type:  PRIVATE
//Description:  save the facet eval tool to a cubit file  
// Assumption:  contained facets have been previuosly saved.  This function
//              saves only the facet ids.
//===========================================================================
CubitStatus FacetEvalTool::save( 
  FILE *fp )
{
  NCubitFile::CIOWrapper cio(fp);
  typedef NCubitFile::UnsignedInt32 int32;

  cio.Write(reinterpret_cast<int32*>(&interpOrder), 1);
  cio.Write(reinterpret_cast<int32*>(&isFlat), 1);
  cio.Write(reinterpret_cast<int32*>(&isParameterized), 1);

  cio.Write(&myArea, 1);
  cio.Write(&minDot, 1);

   // write ids of facets that belong to this tool 
  int ii;
  CubitFacet *facet_ptr;
  int nfacets = myFacetList.size();
  int32* facet_id = new int32 [nfacets];
  myFacetList.reset();
  for (ii=0; ii<nfacets; ii++)
  {
    facet_ptr = myFacetList.get_and_step();
    facet_id[ii] = facet_ptr->id();
  }
  cio.Write(reinterpret_cast<int32*>(&nfacets), 1);
  if (nfacets > 0)
  {
    cio.Write(facet_id, nfacets);
  }
  delete [] facet_id;

   // write ids of edges that belong to this tool 
  CubitFacetEdge *edge_ptr;
  int32 nedges = myEdgeList.size();
  int32* edge_id = new int32 [nedges];
  myEdgeList.reset();
  for (ii=0; ii<nedges; ii++)
  {
    edge_ptr = myEdgeList.get_and_step();
    edge_id[ii] = edge_ptr->id();
    volatile int test_int = edge_id[ii] + 1;   // this is a test line to look for uninitialized data rjm
  }
  cio.Write( &nedges, 1);
  if (nedges > 0)
  {
    cio.Write(edge_id, nedges);
  }
  delete [] edge_id;

   // write ids of points that belong to this tool
  CubitPoint *point_ptr;
  int npoints = myPointList.size();
  int32* point_id = new int32 [npoints];
  myPointList.reset();
  for (ii=0; ii<npoints; ii++)
  {
    point_ptr = myPointList.get_and_step();
    point_id[ii] = point_ptr->id();
  }
  cio.Write(reinterpret_cast<int32*>(&npoints), 1);
  if (npoints > 0)
  {
    cio.Write(point_id, npoints);
  }
  delete [] point_id;

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: restore
//
//Member Type:  PRIVATE
//Description:  restore a facetevaltool from a CUB file  
//===========================================================================
CubitStatus FacetEvalTool::restore( 
  FILE *fp,
  unsigned int endian,
  int num_facets, 
  int num_edges, 
  int num_points,
  CubitFacet **facets, 
  CubitFacetEdge **edges, 
  CubitPoint **points )
{
  NCubitFile::CIOWrapper cio(endian, fp);
  typedef NCubitFile::UnsignedInt32 int32;

  // read interpOrder, isFlat, isParameterized 
  int int_data[3];
  cio.Read(reinterpret_cast<int32*>(int_data), 3);
  interpOrder = int_data[0];
  isFlat = int_data[1];
  isParameterized = (int_data[2] != 0);

  // read myArea, minDot
  double double_data[2];
  cio.Read( double_data, 2);
  myArea = double_data[0];
  minDot  = double_data[1];

  lastFacet = NULL;

  // read the facet ids and assign facets to this eval tool
  int nfacets; 
  cio.Read(reinterpret_cast<int32*>(&nfacets), 1);
  int32 *facet_id = NULL;
  if (nfacets > 0)
  {
    facet_id = new int32 [nfacets];
    cio.Read(facet_id, nfacets);
    int ii, id;
    for(ii=0; ii<nfacets; ii++)
    {
      id = facet_id[ii];
      if (ii < 0 || ii >= num_facets)
      {
        delete [] facet_id;
        return CUBIT_FAILURE;
      }
      myFacetList.append( facets[id] );
      facets[id]->set_tool_id( toolID );
    }
    delete [] facet_id;
  }

  // read the edges
  int nedges; 
  cio.Read(reinterpret_cast<int32*>(&nedges), 1);
  int32 *edge_id = NULL;
  if (nedges > 0)
  {
    edge_id = new int32 [nedges];
    cio.Read(edge_id, nedges);
    int ii, id;
    for(ii=0; ii<nedges; ii++)
    {
      id = edge_id[ii];
      if (ii < 0 || ii >= num_edges)
      {
        delete [] edge_id;
        return CUBIT_FAILURE;
      }
      myEdgeList.append( edges[id] );
    }
    delete [] edge_id;
  }

  // read the points
  int npoints; 
  cio.Read(reinterpret_cast<int32*>(&npoints), 1);
  int32 *point_id = NULL;
  if (npoints > 0)
  {
    point_id = new int32 [npoints];
    cio.Read(point_id, npoints);
    int ii, id;
    for(ii=0; ii<npoints; ii++)
    {
      id = point_id[ii];
      if (ii < 0 || ii >= num_points)
      {
        delete [] point_id;
        return CUBIT_FAILURE;
      }
      myPointList.append( points[id] );
    }
    delete [] point_id;
  }

  bounding_box();
  double diag = sqrt(sqr(myBBox->x_range()) + 
	                   sqr(myBBox->y_range()) + 
	                   sqr(myBBox->z_range()));
  compareTol = 1.0e-3 * diag;

  return CUBIT_SUCCESS;
}

CubitStatus FacetEvalTool::get_intersections(CubitVector point1,
                                             CubitVector point2,
                                             DLIList<CubitVector*>& intersection_list,
                                             bool bounded)
{
    //Find the points were the line intersects the bounding box.
  CubitVector intersect_1;
  CubitVector intersect_2;
  CubitBox bbox = *myBBox;

    //Increase the size of the box in each direction.
    //Don't use scale because the box may be too thin (planar surface).
  double offset = 2.0 * (point1 - point2).length();
  
  CubitVector min;
  min.x( bbox.min_x() - offset );
  min.y( bbox.min_y() - offset );
  min.z( bbox.min_z() - offset );
  CubitVector max;
  max.x( bbox.max_x() + offset );
  max.y( bbox.max_y() + offset );
  max.z( bbox.max_z() + offset );

  bbox.reset( min, max );
  int box_intersections = FacetDataUtil::get_bbox_intersections( point1, point2, bbox,
                                                                 intersect_1, intersect_2 );

    //The bounding box is larger than the surface we are checking.
    //This means that if there are less than two intersections
    //the line will not intersect the surface.
  if( 2 > box_intersections )
      return CUBIT_SUCCESS;

  bbox.reset( intersect_1 );
  bbox |= intersect_2;
  
    //Find the facets that are intersected by the bbox that was just created.
  DLIList<CubitFacet*> search_facets;
  if( aTree )
  {
      //Get the facets from the tree.
    aTree->find( bbox, search_facets );
  }
  else
      search_facets = myFacetList;

  int ii;
  for( ii = search_facets.size(); ii > 0; ii-- )
  {
    CubitFacet* test_facet = search_facets.get_and_step();

    CubitVector intersection;
    CubitVector area_coord;
    CubitPointContainment contain = FacetDataUtil::intersect_facet( intersect_1, intersect_2, test_facet,
                                                                    intersection, area_coord );

    if( bounded )
    {
      CubitVector dir1 = point2 - point1;
      CubitVector dir2 = intersection - point1;

      if( dir2.length_squared() > (GEOMETRY_RESABS * GEOMETRY_RESABS) )
      {
        if( dir1 % dir2 < 0 ||
            ( ( dir2.length_squared() - dir1.length_squared() ) >
              ( GEOMETRY_RESABS * GEOMETRY_RESABS ) ) )
        {
            //The inserction point is not between the two end points.
          contain = CUBIT_PNT_OUTSIDE;
        }
      }
    }
    
    if( CUBIT_PNT_BOUNDARY == contain ||
        CUBIT_PNT_INSIDE == contain )
    {
        //The point intersects the facets.
      CubitVector* new_intersection = new CubitVector;
      *new_intersection = intersection;
      intersection_list.append( new_intersection );
    }
  }

    //Remove duplicate intersections.
  intersection_list.reset();
  for( ii = 0; ii < intersection_list.size(); ii++ )
  {
    CubitVector* base_vec = intersection_list.next(ii);
    if( NULL == base_vec )
        continue;
    
    int jj;
    for( jj = ii+1; jj < intersection_list.size(); jj++ )
    {
      CubitVector* compare_vec = intersection_list.next(jj);
      if( NULL != compare_vec )
      {
        if( base_vec->distance_between_squared( *compare_vec ) < GEOMETRY_RESABS * GEOMETRY_RESABS )
        {
          intersection_list.step(jj);
          delete intersection_list.get();
          intersection_list.change_to( NULL );
          intersection_list.reset();
        }
      }
    }
  }
  intersection_list.remove_all_with_value( NULL );
  
  return CUBIT_SUCCESS;
}

int FacetEvalTool::intersect_ray( CubitVector &origin, CubitVector &direction, CubitFacet* facet, CubitVector* point, double &distance )
{
	// This algorithm can be found at http://geometryalgorithms.com/

	CubitVector n;           // triangle vectors
    CubitVector w0, w;       // ray vectors
    double a, b;             // params to calc ray-plane intersect

    // get triangle edge vectors and plane normal
	CubitVector normal = facet->normal();
	CubitPoint *pt0 = facet->point(0);
	CubitPoint *pt1 = facet->point(1);
	CubitPoint *pt2 = facet->point(2);
	double tol = GEOMETRY_RESABS;
	
	CubitVector u( pt1->x() - pt0->x(),
					pt1->y() - pt0->y(),
					pt1->z() - pt0->z()); //(*p1-*p0);
	CubitVector v( pt2->x() - pt0->x(),
					pt2->y() - pt0->y(),
					pt2->z() - pt0->z()); // = (*p2-*p0);

	//u = T.V1 - T.V0;
    //v = T.V2 - T.V0;
    n = u * v;             // cross product
    if (n.length_squared() == 0)   // triangle is degenerate
        return -1;                 // do not deal with this case

    //dir = R.P1 - R.P0;             // ray direction vector
    //w0 = R.P0 - T.V0;
	w0 = CubitVector(origin.x() - pt0->x(),
		origin.y() - pt0->y(),
		origin.z() - pt0->z());

    a = -(n%w0);
    b = (n%direction);
    if (fabs(b) < tol) {     // ray is parallel to triangle plane
        if (a == 0)                // ray lies in triangle plane
            return 2;
        else return 0;             // ray disjoint from plane
    }

    // get intersect point of ray with triangle plane
    distance = a / b;
    if (distance < 0.0)                   // ray goes away from triangle
        return 0;                  // => no intersect
    // for a segment, also test if (r > 1.0) => no intersect

    point->set(origin + distance * direction);           // intersect point of ray and plane

    // the distance we want to return is real distance, not distance/magnitude
    distance *= direction.length();

    // is point inside facet?
    double uu, uv, vv, wu, wv, D;
    uu = u%u;
    uv = u%v;
    vv = v%v;
    //w = *I - T.V0;
	w = CubitVector(point->x() - pt0->x(),
					point->y() - pt0->y(),
					point->z() - pt0->z());
    wu = w%u;
    wv = w%v;
    D = uv * uv - uu * vv;

    // get and test parametric coords
    double s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)        // point is outside facet
        return 0;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // point is outside facet
        return 0;

    return 1;                      // point is in facet

}

int FacetEvalTool::intersect_ray( CubitVector &origin, CubitVector &direction, CubitFacetEdge* facet_edge, CubitVector* point, double &hit_distance )
{
	// This algorithm can be found at http://geometryalgorithms.com/
	double tol = GEOMETRY_RESABS;

	CubitPoint* p0 = facet_edge->point(0);
	CubitPoint* p1 = facet_edge->point(1);
	CubitVector u0 = CubitVector(p0->x(), p0->y(), p0->z());
    CubitVector u1 = CubitVector(p1->x(), p1->y(), p1->z());
	
	CubitVector u = CubitVector(u0, u1);
	CubitVector v = direction;
	v.normalize();

	CubitVector w = CubitVector(origin, u0);

	double sc, tc;         // sc is fraction along facet edge, tc is distance along ray
	
	double a = u%u;        // always >= 0
    double b = u%v;
    double c = v%v;        // always >= 0
    double d = u%w;
    double e = v%w;
    double D = a*c - b*b;  // always >= 0

    // compute the line parameters of the two closest points
    if (D < tol)
	{
		// the lines are almost parallel
        sc = 0.0;
        tc = (b>c ? d/b : e/c);   // use the largest denominator
    }
    else
	{
        sc = (b*e - c*d) / D;
        tc = (a*e - b*d) / D;
    }

    // get the difference of the two closest points
    CubitVector dP = CubitVector(w + (sc * u) - (tc * v));  // = <0 0 0> if intersection

    double distance = sqrt(dP % dP); // return the closest distance (0 if intersection)

	point->set(u0 + (sc * u));
	hit_distance = tc; //distance from origin to intersection point

	if (distance < tol)
	{
		//check if parallel (infinite intersection)
		if (D < tol)
			return 2;
		//check if on edge
		if (sc <= 1.0 && sc >= 0.0)
			return 1;
		else
			return 0;
	}

	return 0;
}
