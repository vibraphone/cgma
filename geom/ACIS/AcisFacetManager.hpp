//---------------------------------------------------------------------------
//- Class:          AcisFacetManager
//- Description:    Receives the facet information as Acis
//-                 generates it.  This is faster than copying
//-                 the facets from an Acis-supplied mesh manager.
//-
//- Owner:          Darryl Melander
//---------------------------------------------------------------------------
#ifndef ACISFACETMANAGER_HPP
#define ACISFACETMANAGER_HPP

#if CUBIT_ACIS_VERSION < 1100
#include "faceter/meshmgr/meshmg.hxx"
#else
#include "meshmg.hxx"
#endif

#include "AcisTypes.h"

class GMem;

class AcisFacetManager : public MESH_MANAGER
{   
public:
     // Constructor just initializes the pointer
   AcisFacetManager() : gMem(NULL), numPoly(0)
      {}
   
     // Destructor shouldn't do anything because
     // it never allocates its own members
   ~AcisFacetManager() 
      {}
   
     // This is how the calling code sets where to store the data
   void set_gmem(GMem* new_gmem) 
      { gMem = new_gmem; }
   
     // This tells the calling code how many polygons there are.
     // You can't tell from gMem unless all polygons are triangles
     // (which they are if you set the refinement the way we usually do).
   int polygon_count() const
      { return numPoly; }
   
     // All the rest of the functions are invoked by the faceter
     // to construct the mesh.

     // This is called to say the 'entity' is about to
     // be faceted w/ refinement 'app_ref', using output
     // format 'format'
   virtual void begin_mesh_output(
      ENTITY *entity,
      ENTITY *app_ref,
      ENTITY *format);
   
     // This says we are done faceting 'entity'
   virtual void end_mesh_output(
      ENTITY *, //entity,
      ENTITY *, //app_ref,
      ENTITY *) //format) 
      {} // We don't do anything.
   
     // This says we want to attach the facets to 'entity'
   virtual void save_mesh_output(
      ENTITY *, //entity,
      ENTITY *, //app_ref,
      ENTITY *) //format) 
      {} // We don't do anything.
   
     // This is the first thing called, and says how many polygons,
     // nodes, and node references there are.
   virtual void announce_counts(int npoly,int nnode,int nref);

     // This says to recompute the points along curves
     // when the attached face is faceted a second time.
     // Otherwise, wierd things happen.
   virtual AF_EDGE_DIRECTIVE check_edge_refinement(EDGE *, double &, double &, double &, int, int)
      { return AF_EDGE_RECOMPUTE; }
   virtual logical need_edge_grading(double &)
      { return FALSE; }
     // This returns true because we want to use the indexed mode
   virtual logical need_indexed_polygons()
      { return TRUE; }
     // We don't need duplicates
   virtual logical need_duplicate_indexed_nodes_on_surface_seams()
      { return FALSE; }
   virtual logical need_indexed_polyedges()
      { return FALSE; }
   
     // Tells the facet manager about a node,
     // including its index, parametric location,
     // world coordinate location, and normal vector.
   virtual void *announce_indexed_node(
      int inode,               // Index in array
      const SPApar_pos &uv,       // parametric coords
      const SPAposition &X,       // world location coords
      const SPAunit_vector &N);   // normal vector
   
     // Says we are about to start a new polygon with
     // index 'ipoly' containing 'npolynode' nodes.
   void start_indexed_polygon(int ipoly, int npolynode, int ishare);
   void start_indexed_polygon(int ipoly, int npolynode);
     // Comment - Between Acis 3.x and 4.x, they added a parameter, but they
     // didn't leave the old version there!!!  The result is that we have to
     // have both, one calling the other.
   
#ifdef BOYD14
     // This is just to avoid hiding the base classes' function
   virtual void announce_indexed_polynode(ENTITY *E,int i1,int i2,void* v)
    { MESH_MANAGER::announce_indexed_polynode(E, i1, i2, v); }
#endif

  virtual void announce_indexed_polynode(ENTITY* E, int i1, int i2, void* v,
                                         const double& d, const SPApar_pos& pp,
                                         const SPAposition& p, const SPAunit_vector& uv)
    { MESH_MANAGER::announce_indexed_polynode(E, i1, i2, v,
                                              d, pp, p, uv); }
   
     // Tells the facet manager that polygon 'ipoly' has node 'i'
     // found in location 'pnode' in the global array.
   virtual void announce_indexed_polynode(int ipoly,int i,void *pnode);

     // Tell it that a Null pointer means invalid node
   virtual void* null_node_id()
      { return (void*)0; }
   
private:
   GMem* gMem;
   int numPoly;
};

#endif
