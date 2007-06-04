//-------------------------------------------------------------------
//- Class:          AcisFacetManager
//- Description:    Receives the facet information as Acis
//-                 generates it.  This is faster than copying
//-                 the facets from an Acis-supplied mesh manager.
//-
//- Owner:          Darryl Melander
//-------------------------------------------------------------------
#include "AcisFacetManager.hpp"
#include "CubitMessage.hpp"
#include "GMem.hpp"
#include <assert.h>

//#define ACIS_FACET_MAN_DEBUG

#if CUBIT_ACIS_VERSION < 1100
#include "baseutil/vector/position.hxx"
#else
#include "position.hxx"
#endif

// #define ACIS_FACET_MAN_DEBUG

// This is called to say the 'entity' is about to
// be faceted w/ refinement 'app_ref', using output
// format 'format'
void AcisFacetManager::begin_mesh_output(
   ENTITY *, //entity,
   ENTITY *, //app_ref,
   ENTITY *  //format)
   )
{
   if (!gMem)
   {
      PRINT_ERROR("Unable to display surface.  Nowhere to "
                  "store facets!\n");
      assert(0);
   }
   gMem->clean_out();
   gMem->points_consolidated(CUBIT_TRUE);

#ifdef ACIS_FACET_MAN_DEBUG
   PRINT_INFO("Called AFM::begin_mesh_output()\n");
#endif
}

// This is the first thing called, and says how many polygons,
// nodes, and node references there are.
// nref is how many nodes there would be if each
// polygon had its own copy of each node.
void AcisFacetManager::announce_counts(int npoly, int nnode, int nref)
{
     // Make sure the arrays are big enough
   if (gMem->point_list_size() < nnode)
      gMem->replace_point_list(new GPoint[nnode], 0, nnode);
   if (gMem->facet_list_size() < npoly + nref)
      gMem->replace_facet_list(new int[npoly+nref], 0, npoly+nref);
     // Make sure arrays are empty
   gMem->clean_out();
     // store how many polygons there are.
   numPoly = npoly;
   
#ifdef ACIS_FACET_MAN_DEBUG
   PRINT_INFO("Called AFM::announce_counts(%d poly, %d points, %d refs)\n",
              npoly, nnode, nref);
#endif
}

// Tells the facet manager about a node,
// including its index, parametric location,
// world coordinate location, and normal SPAvector.
void *AcisFacetManager::announce_indexed_node(
      int inode,               // Index in array
      const SPApar_pos&,          // parametric coords
      const SPAposition& X,       // world location coords
      const SPAunit_vector&)      // normal SPAvector
{
     // Store the coordinates in the gMem->points array
   gMem->point_list()[inode].x = (float)X.x();
   gMem->point_list()[inode].y = (float)X.y();
   gMem->point_list()[inode].z = (float)X.z();

   gMem->pointListCount++;
   
#ifdef ACIS_FACET_MAN_DEBUG
   PRINT_INFO("Called AFM::announce_indexed_node(index %d)\n",
              inode);
#endif
   
     // Return a pointer to the data so it can be
     // retrieved when forming the polygons
   return (void*)(gMem->point_list() + inode);
}

// Says we are about to start a new polygon with
// index 'ipoly' containing 'npolynode' nodes.
void AcisFacetManager::start_indexed_polygon(int ipoly, int npolynode, int)
{ 
  start_indexed_polygon(ipoly, npolynode); 
}

void AcisFacetManager::start_indexed_polygon(int, int npolynode)
{
     // Store how many nodes are in the polygon
   gMem->facet_list()[gMem->fListCount++] = npolynode;
#ifdef ACIS_FACET_MAN_DEBUG
   PRINT_INFO("Called AFM::start_indexed_polygon(%d points in poly)\n",
              npolynode);
#endif
}

// Tells the facet manager that the 'i'th node of
// polygon 'ipoly' is pointed at by 'pnode'.
void AcisFacetManager::announce_indexed_polynode(int, // ipoly,
                                                 int i,
                                                 void* pnode)
{
     // Get the index
   i = ((GPoint*)pnode - gMem->point_list());
     // Store the index
   gMem->facet_list()[gMem->fListCount++] = i;

   
#ifdef ACIS_FACET_MAN_DEBUG
   PRINT_INFO("Called AFM::annonce_indexed_polynode(index %d)\n",
              i);
#endif
}

