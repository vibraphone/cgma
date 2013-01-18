//-Class:       BoundaryConstrainTool
//-Description: recover edges from a triangle mesh.  Currently 2D meshes only
//-Owner:       Steve Owen
//-Checked by:


#include "GeometryDefines.h"
#include "CubitVector.hpp"
#include "BoundaryConstrainTool.hpp"

#define ON_LINE_TOL 1.0e-7
#ifndef SQR
#define SQR(x) ((x) * (x))
#endif

//-----------------------------------------------------------------------------
// Function:    BoundaryConstrainTool
// Type:        Public
// Description: constructor
// Author:      sjowen
// Altered by:  jitken(07/08/02)
// Date:        2/17/02
//-----------------------------------------------------------------------------
template <class SURF, class TRI, class EDGE, class NODE, class TRICHILD> 
BoundaryConstrainTool<SURF,TRI,EDGE,NODE,TRICHILD>::BoundaryConstrainTool(SURF *ref_face_ptr)
{
  refFacePtr = ref_face_ptr;
}


//-----------------------------------------------------------------------------
// Function:    BoundaryConstrainTool
// Type:        Public
// Description: constructor
// Author:      sjowen
// Altered by:  jitken(07/08/02)
// Date:        2/17/02
//-----------------------------------------------------------------------------
template <class SURF, class TRI, class EDGE, class NODE, class TRICHILD> 
BoundaryConstrainTool<SURF,TRI,EDGE,NODE,TRICHILD>::BoundaryConstrainTool()
{
  refFacePtr = NULL;
}

//-----------------------------------------------------------------------------
// Function:    ~BoundaryConstrainTool
// Type:        Public
// Description: desctructor
// Author:      sjowen
// Date:        2/17/02
//-----------------------------------------------------------------------------
template <class SURF, class TRI, class EDGE, class NODE, class TRICHILD> 
BoundaryConstrainTool<SURF,TRI,EDGE,NODE,TRICHILD>::~BoundaryConstrainTool()
{
}

//-----------------------------------------------------------------------------
// Function:    ~BoundaryConstrainTool
// Type:        Public
// Description: desctructor
// Author:      sjowen
// Date:        2/17/02
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function:    recover_edge
// Type:        Public
// Description: flip edges in the triangulation to recover the edge between
//              n0_ptr and n1_ptr. Assumes 2D (x-y-0)
// Author:      sjowen
// Altered by:  jitken(07/08/02)
// Date:        2/17/02
//-----------------------------------------------------------------------------
template <class SURF, class TRI, class EDGE, class NODE, class TRICHILD> 
CubitStatus BoundaryConstrainTool<SURF,TRI,EDGE,NODE,TRICHILD>::recover_edge(
  NODE *n0_ptr,
  NODE *n1_ptr,
  EDGE *&recovered_edge_ptr,
  DLIList <TRI *> *facet_list)
{
  assert(n0_ptr && n1_ptr);
  CubitStatus rv = CUBIT_SUCCESS;
  
  //Update the private Variable
  facetList = facet_list;

  // TRIvial case: edge is already in the triangulation.  Return
  // now with the edge
  
  recovered_edge_ptr = n0_ptr->shared_edge(n1_ptr);

  if ( recovered_edge_ptr != NULL ) 
    return CUBIT_SUCCESS;
  
 
  // Define the vector from node0 to node1.  This is the direction
  // of the edge we want to recover
  
  
  CubitVector edge_vec = n1_ptr->coordinates() - n0_ptr->coordinates();
  double length = edge_vec.normalize( );
  if (length == 0.0)
    return CUBIT_FAILURE;
  
  zeroTol = length * CUBIT_RESABS;
  
  // Initialize the edge cross queue
  
  edgeCrossQueue.clean_out();
  
  // Fill in the Queue with edges that cross the edge_vec
  
  rv = get_crossing_edges( edge_vec, n0_ptr, n1_ptr );
  if (rv != CUBIT_SUCCESS)
    return rv;
  
  // Process each edge in the queue and swap edges on faces
  // adjacent to the edges crossing the recover vector
    
  EDGE *edge_ptr;
  SwapStatus swap_stat;
  while (edgeCrossQueue.size() > 0)
  {
    edge_ptr = edgeCrossQueue.remove();
    swap_stat = swap_edge( edge_ptr );
    switch (swap_stat) 
    {
      
    case SWAP_FAILURE:
      return CUBIT_FAILURE;
      
    case SWAP_INVALID:
      
      // If flip was not valid, then put it back at the end of
      // the queue to be processed later
      
      edgeCrossQueue.append( edge_ptr );
      break;
      
    case SWAP_SUCCESS:
      
      // Check if the edge just created by the flip also
      // crosses the recover vector.  If it does, then
      // add these at the end of the queue to be processed later
      
      if (edgeCrossQueue.size() > 0)
      {
        rv = edge_intersected( edge_vec, edge_ptr, n0_ptr, n1_ptr );
        if (rv != CUBIT_SUCCESS)
          return rv;
      }
    }
  }
  
  // Retrieve the recovered edge to pass back
  
  recovered_edge_ptr = n0_ptr->shared_edge(n1_ptr);
  if ( recovered_edge_ptr == NULL ) 
    rv = CUBIT_FAILURE;
  
  return rv;
}



//-----------------------------------------------------------------------------
// Function:    get_crossing_edges
// Type:        Private
// Description: generate the initial list of all edges that cross the edge_vec
//              between n0 and n1  
// Author:      sjowen
// Date:        2/17/02
//-----------------------------------------------------------------------------
template <class SURF, class TRI, class EDGE, class NODE, class TRICHILD> 
CubitStatus BoundaryConstrainTool<SURF,TRI,EDGE,NODE,TRICHILD>::get_crossing_edges( 
  CubitVector &edge_vec, 
  NODE *n0_ptr, 
  NODE *n1_ptr )
{
  CubitStatus status = CUBIT_FAILURE;
  IntersectionType isect_type = AT_BEGIN;
  IntersectionType last_isect_type = NO_ISECT;
  NODE *node_ptr;
  EDGE *edge_ptr;
  TRI *tri_ptr;
  CubitBoolean done = CUBIT_FALSE;
  
  do 
  {
    switch(isect_type) 
    {
    case AT_BEGIN:
      last_isect_type = AT_BEGIN;
      isect_type = intersect_from_node( edge_vec, n0_ptr, 
                                        tri_ptr, edge_ptr, node_ptr );
      break;
    case AT_EDGE:
      last_isect_type = AT_EDGE;
      edgeCrossQueue.append( edge_ptr );
      isect_type = intersect_from_edge( edge_vec, n0_ptr, n1_ptr,
                                        tri_ptr, edge_ptr, node_ptr );
      break;
    case AT_MID:
      status = node_at_mid( edge_vec, tri_ptr, node_ptr );
      if (status != CUBIT_SUCCESS) {
        return status;
      }
      isect_type = last_isect_type;
      break;
    case AT_END:
      done = CUBIT_TRUE;
      status = CUBIT_SUCCESS;
      break;
    case NO_ISECT:
      done = 1;
      return CUBIT_FAILURE;
    }
  } while (!done);

  return status;
}

//-----------------------------------------------------------------------------
// Function:    intersect_from_node
// Type:        Private
// Description: determine next intersection assuming start from n0
// Author:      sjowen
// Date:        2/17/02
//-----------------------------------------------------------------------------
template <class SURF, class TRI, class EDGE, class NODE, class TRICHILD> 
IntersectionType BoundaryConstrainTool<SURF,TRI,EDGE,NODE,TRICHILD>::intersect_from_node(
  CubitVector &edge_vec, 
  NODE *n0_ptr, 
  TRI *&tri_ptr, 
  EDGE *&edge_ptr, 
  NODE *&node_ptr ) 
{
  IntersectionType isect_type = NO_ISECT;

  // get all tris adjacent to n0_ptr
  
  DLIList<TRI *> tri_list;
  n0_ptr->tris( tri_list );
  CubitBoolean found = CUBIT_FALSE;
  
  int itri;
  for (itri=0; itri<tri_list.size() && !found; itri++) 
  {
    tri_ptr = tri_list.get_and_step();
      
    // If the dot product of edges radiating from the node
    // are both greater than zero, then it is a candidate 
      
    NODE *n1_ptr = tri_ptr->next_node( n0_ptr );
    assert( n1_ptr != NULL );
    edge_ptr = n0_ptr->shared_edge( n1_ptr );
    assert( edge_ptr != NULL );

    // define the normal vector to the edge

    CubitVector vec = n1_ptr->coordinates() - n0_ptr->coordinates();
    vec.z( vec.x() );
    vec.x( -vec.y() );
    vec.y( vec.z() );
    vec.z( 0.0 );
    double length = sqrt(SQR(vec.x()) + SQR(vec.y()));
    if (length == 0.0)
    {
      node_ptr = NULL;
      edge_ptr = NULL;
      tri_ptr = NULL;
      return NO_ISECT;
    }
    vec.x( vec.x() / length );
    vec.y( vec.y() / length );

    // dot with the recover vector

    double dot = edge_vec.x() * vec.x() + edge_vec.y() * vec.y();
         
    if (dot > -ON_LINE_TOL) {

      // Dot product within tolerance of zero means the vector
      // follows along edge of face and passes through a node

      if (fabs(dot) < ON_LINE_TOL) {
	vec = n1_ptr->coordinates() - n0_ptr->coordinates();
	double direction = edge_vec.x() * vec.x() + edge_vec.y() * vec.y();
	//if vector n1-n0 is in the same direction as edge_vec
	if(direction > 0){
	  vec = n1_ptr->coordinates() - n0_ptr->coordinates();
	  node_ptr = n1_ptr;
	  isect_type = AT_MID;
	  return isect_type;
	}
      }

      // do the same check on the other edge
      
      NODE *n2_ptr = tri_ptr->next_node( n1_ptr );
      assert( n2_ptr != NULL );
      edge_ptr = n0_ptr->shared_edge( n2_ptr );
      assert( edge_ptr != NULL );
      vec = n0_ptr->coordinates() - n2_ptr->coordinates();
      vec.z( vec.x() );
      vec.x( -vec.y() );
      vec.y( vec.z() );
      vec.z( 0.0 );
      length = sqrt(SQR(vec.x()) + SQR(vec.y()));
      if (length == 0.0)
      {
        node_ptr = NULL;
        edge_ptr = NULL;
        tri_ptr = NULL;
        return NO_ISECT;
      }
      vec.x( vec.x() / length );
      vec.y( vec.y() / length );
      dot = edge_vec.x() * vec.x() + edge_vec.y() * vec.y();
      if (fabs(dot) < ON_LINE_TOL) {
	node_ptr = n2_ptr;
	isect_type = AT_MID;
	return isect_type;
      }
      
      if (dot > -ON_LINE_TOL){
	node_ptr = NULL;
	edge_ptr = n1_ptr->shared_edge( n2_ptr );
	isect_type = AT_EDGE;
	return isect_type;
      }
    }
  }
  
  if (!found)
  {
    isect_type = NO_ISECT;
    node_ptr = NULL;
    edge_ptr = NULL;
    tri_ptr = NULL;
  }

  return isect_type;
}

//-----------------------------------------------------------------------------
// Function:    intersect_from_edge
// Type:        Private
// Description: determine next intersection assuming start from edge_ptr
// Author:      sjowen
// Date:        2/17/02
//-----------------------------------------------------------------------------
template <class SURF, class TRI, class EDGE, class NODE, class TRICHILD> 
IntersectionType BoundaryConstrainTool<SURF,TRI,EDGE,NODE,TRICHILD>::intersect_from_edge( 
  CubitVector &edge_vec, 
  NODE *n0_ptr, 
  NODE *n1_ptr,
  TRI * &tri_ptr, 
  EDGE *&edge_ptr, 
  NODE *&node_ptr )
{
  IntersectionType isect_type = NO_ISECT;
  assert( edge_ptr != NULL && tri_ptr != NULL);

  // Get the adjacent face to the edge - there must be exactly two faces
  // next to this edge otherwise we've left the triangulation
  
  DLIList<TRI *> adjtris;
  edge_ptr->tris( adjtris );
  if (adjtris.size() != 2) 
    return NO_ISECT;
  TRI * nexttri_ptr = adjtris.get_and_step();
  if (nexttri_ptr == tri_ptr)
  {
    nexttri_ptr = adjtris.get();
    assert(nexttri_ptr != NULL && nexttri_ptr != tri_ptr);
  }
  tri_ptr = nexttri_ptr;
  
  // Check if we've arrived at the end (Does the this triangle contain n1)
  
  int ii;
  NODE *n_ptr[3];
  nexttri_ptr->tri_nodes( n_ptr[0], n_ptr[1], n_ptr[2] );
  CubitBoolean found = CUBIT_FALSE;
  for (ii=0; ii<3 && !found; ii++) {
    if (n_ptr[ii] == n1_ptr) {
      found = CUBIT_TRUE;
      edge_ptr = NULL;
      node_ptr = n1_ptr;
      tri_ptr = nexttri_ptr;
      isect_type = AT_END;
    }
  }
  
  // Determine which edge (or node) the vector intersects
  
  if (!found) {
    NODE *tn0_ptr, *tn1_ptr, *tn2_ptr;
    tn0_ptr = edge_ptr->start_node();
    tn1_ptr = edge_ptr->end_node();
    tn2_ptr = tri_ptr->next_node(tn1_ptr);
    //Checking for CCW order. Correct it if it is wrong
    if(tn2_ptr == tn0_ptr){
      tn0_ptr = tn1_ptr;
      tn1_ptr = tn2_ptr;
      tn2_ptr = tri_ptr->next_node(tn1_ptr);
    }
    
    // Determine vector from n0 to tn2
    
    CubitVector vec( tn2_ptr->coordinates().x() - n0_ptr->coordinates().x(),
                     tn2_ptr->coordinates().y() - n0_ptr->coordinates().y(),
                     0.0 );
    double len = sqrt(SQR(vec.x()) + SQR(vec.y()));
    if (len == 0.0)
      return NO_ISECT;
    vec.x( vec.x()/len );
    vec.y( vec.y()/len );
    
    // Compute the dot product of the normal to vec (above)
    // with the edgevec. 
    // One of the following will result:
    // Dot > 0:  edge defined by tn2 - tn0 is intersected
    // Dot < 0:  edge defined by tn1 - tn2 is intersected
    // Dot = 0:  node tn2 is intersected  */

    vec.z( vec.x() );
    vec.x( -vec.y() );
    vec.y( vec.z() );
    vec.z( 0.0 );

    double dot = vec.x() * edge_vec.x() + vec.y() * edge_vec.y();

    if (fabs(dot) < ON_LINE_TOL) {
      isect_type = AT_MID;
      edge_ptr = NULL;
      node_ptr = tn2_ptr;
    }
    else if (dot > 0) {
      isect_type = AT_EDGE;
      edge_ptr = tn0_ptr->shared_edge( tn2_ptr );
      node_ptr = NULL;
    }
    else {
      isect_type = AT_EDGE;
      edge_ptr = tn1_ptr->shared_edge( tn2_ptr );
      node_ptr = NULL;
    }
  }
  
  return isect_type;
}

//-----------------------------------------------------------------------------
// Function:    node_at_mid
// Type:        Private
// Description: the vector hit a node on its way to the n1
// Author:      sjowen
// Date:        2/17/02
//-----------------------------------------------------------------------------
template <class SURF, class TRI, class EDGE, class NODE, class TRICHILD> 
CubitStatus BoundaryConstrainTool<SURF,TRI,EDGE,NODE,TRICHILD>::node_at_mid( 
  CubitVector &edge_vec, 
  TRI * tri_ptr, 
  NODE *node_ptr )
{
  // don't know what to do here yet
  return CUBIT_FAILURE;
}

//-----------------------------------------------------------------------------
// Function:    swap_edge
// Type:        Private
// Description: attempt to swap a single edge.  Check for valid swap 
//              before doing so
// Author:      sjowen
// Date:        2/17/02
//-----------------------------------------------------------------------------
template <class SURF, class TRI, class EDGE, class NODE, class TRICHILD> 
SwapStatus BoundaryConstrainTool<SURF,TRI,EDGE,NODE,TRICHILD>::swap_edge(  
  EDGE *&edge_ptr )
{

  DLIList<TRI *>adjtris;
  edge_ptr->tris(refFacePtr, adjtris);
  if(adjtris.size() != 2)
    return SWAP_FAILURE;

  // check the potential new triangles

  TRI * tri0 = adjtris.get_and_step();
  TRI * tri1 = adjtris.get();
  NODE *n0, *n1, *n2, *n3;
  n0 = edge_ptr->start_node();
  n1 = edge_ptr->end_node();
  n2 = tri0 ->next_node(n1);
  //if direction is not in CCW order, correct it.
  if(n2 == n0){
    n0 = n1;
    n1 = n2;
    n2 = tri0->next_node(n1);
  }
  n3 = tri1->next_node(n0);

  CubitVector e1( n1->coordinates().x() - n3->coordinates().x(),
                  n1->coordinates().y() - n3->coordinates().y(),
                  0.0 );
  double len = sqrt(SQR(e1.x()) + SQR(e1.y()));
  if (len == 0.0)
    return SWAP_FAILURE;
  e1.x( e1.x()/len ); e1.y( e1.y()/len );
  CubitVector diag( n2->coordinates().x() - n3->coordinates().x(),
                    n2->coordinates().y() - n3->coordinates().y(),
                    0.0 );
  len = sqrt(SQR(diag.x()) + SQR(diag.y()));
  if (len == 0.0)
    return SWAP_FAILURE;
  diag.x( diag.x()/len ); diag.y( diag.y()/len );

  CubitVector e2( n0->coordinates().x() - n3->coordinates().x(),
                  n0->coordinates().y() - n3->coordinates().y(),
                  0.0 );
  len = sqrt(SQR(e2.x()) + SQR(e2.y()));
  if (len == 0.0)
    return SWAP_FAILURE;
  e2.x( e2.x()/len ); e2.y( e2.y()/len );

  // return now if we would create invalid triangles

  double cross = e1.x() * diag.y() - e1.y() * diag.x();
  if (cross < ON_LINE_TOL)
    return SWAP_INVALID;

  cross = diag.x() * e2.y() - diag.y() * e2.x();
  if (cross < ON_LINE_TOL)
    return SWAP_INVALID;

  // Do the swap
  // Delete the old triangles

  facetList->move_to(tri0);
  facetList->extract();
  facetList->move_to(tri1);
  facetList->extract();
  delete tri0;
  delete tri1;
  delete edge_ptr;

  // Create the new triangles

  tri0 = (TRI *) new TRICHILD( n1, n2, n3, refFacePtr);
  tri1 = (TRI *) new TRICHILD( n0, n3, n2, refFacePtr);
  facetList->append(tri0);
  facetList->append(tri1);
  edge_ptr = n2->shared_edge( n3 );

  return SWAP_SUCCESS;
}

//-----------------------------------------------------------------------------
// Function:    edge_intersected
// Type:        Private
// Description: check to see if the edge we just generated crosses the
//              edge_vec.  If so, the add it to the queue 
// Author:      sjowen
// Date:        2/17/02
//-----------------------------------------------------------------------------
template <class SURF, class TRI, class EDGE, class NODE, class TRICHILD> 
CubitStatus BoundaryConstrainTool<SURF,TRI,EDGE,NODE,TRICHILD>::edge_intersected( 
   CubitVector &edge_vec, 
   EDGE *edge_ptr, 
   NODE *n0_ptr, 
   NODE *n1_ptr )
{

  // if the end nodes match the end nodes of the edge we are checking, then
  // we are OK - no intersection

  NODE *en0_ptr = edge_ptr->start_node();
  NODE *en1_ptr = edge_ptr->end_node();
  if (en0_ptr == n0_ptr || en1_ptr == n0_ptr ||
      en0_ptr == n1_ptr || en1_ptr == n1_ptr)
    return CUBIT_SUCCESS;

  CubitVector vec1( en0_ptr->coordinates().x() - n0_ptr->coordinates().x(),
                    en0_ptr->coordinates().y() - n0_ptr->coordinates().y(),
                    0.0 );
  double len = sqrt( SQR( vec1.x() ) + SQR( vec1.y() ) );
  if (len == 0.0)
    return CUBIT_FAILURE;
  vec1.x( vec1.x() / len );  vec1.y( vec1.y() / len );

  CubitVector vec2( en1_ptr->coordinates().x() - n0_ptr->coordinates().x(),
                    en1_ptr->coordinates().y() - n0_ptr->coordinates().y(),
                    0.0 );
  len = sqrt( SQR( vec2.x() ) + SQR( vec2.y() ) );
  if (len == 0.0)
    return CUBIT_FAILURE;
  vec2.x( vec2.x() / len );  vec2.y( vec2.y() / len );

  double cross1 = edge_vec.x() * vec1.y() - edge_vec.y() * vec1.x();
  double cross2 = edge_vec.x() * vec2.y() - edge_vec.y() * vec2.x();

  // if the en0 and en1 are both on opposite sides of edge_vec then
  // we know that it intersects somewhere

  if (cross1 * cross2 <= 0.0)
  {
    edgeCrossQueue.append( edge_ptr );
  }

  return CUBIT_SUCCESS;
}
// EOF
