#include "AppUtil.hpp"
#include "CGMApp.hpp"
#include "GeometryQueryTool.hpp"
#include "FacetModifyEngine.hpp"
#include "FacetQueryEngine.hpp"
#include "CubitFacetData.hpp"
#include "CubitPointData.hpp"
#include "CubitQuadFacetData.hpp"
#include "DLIList.hpp"
#include "Surface.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "Body.hpp"
#include "BodySM.hpp"
#include "RefFace.hpp"
#include "RefVertex.hpp"

#include <iostream>

bool generateFacetModel(
  DLIList<TopologyBridge*>& model)
{
  FacetModifyEngine* fme = FacetModifyEngine::instance();

  CubitStatus stat;
  int i;

  // Build a tetrahedron, with 1 triangular facet per faceted-surface:
  double pts[4][3] = {
      { 0., 0., 0. },
      { 1., 0., 0. },
      { 1., 1., 0. },
      { 0., 1., 1. },
  };

  int tconn[4][3] = {
      {0, 2, 1},
      {0, 3, 2},
      {0, 1, 3},
      {1, 2, 3}
  };
  DLIList<Surface*> all_surfs;
  for (int k = 0; k < 2; ++k)
    {
    std::cout << "Build faceted surface " << k << "\n";
    DLIList<CubitFacet*> tflist;
    DLIList<CubitQuadFacet*> qflist;
    DLIList<Surface*> slist;
    DLIList<CubitPoint*> facet_plist;
    std::map<int,int> ptmap;
    int nxt = 0;
    // We will create 2 faceted surfaces with 2 triangles in each:
    for (i = 0; i < 2; ++i)
      {
      // Create shared vertices for the triangles in this faceted surface.
      // They are not shared with the other faceted surface but are shared
      // by facets (triangles) in this surface:
      for (int j = 0; j < 3; ++j)
        {
        int connEnt = tconn[2*k + i][j];
        if (ptmap.find(connEnt) == ptmap.end())
          {
          ptmap[connEnt] = nxt++;
          facet_plist.append(
            new CubitPointData(
              pts[connEnt][0], pts[connEnt][1], pts[connEnt][2]));
          }
        }
      }
    for (i = 0; i < 2; ++i)
      {
      // Create a single facet for our faceted surface:
      tflist.append(
        new CubitFacetData(
          facet_plist[ptmap[tconn[2*k + i][0]]],
          facet_plist[ptmap[tconn[2*k + i][1]]],
          facet_plist[ptmap[tconn[2*k + i][2]]]));
      }
    // Attempt to build a faceted surface:
    stat = fme->build_facet_surface(
      qflist, tflist, facet_plist, -1., 0, false, false, slist);
    if (slist.size() <= 0 || stat != CUBIT_SUCCESS)
      return false;
    all_surfs += slist;
    }

  ShellSM* shell;
  stat = fme->make_facet_shell(all_surfs, shell);
  if (!shell || stat != CUBIT_SUCCESS) return false;

  DLIList<ShellSM*> shlist;
  shlist.append(shell);
  Lump* lump;
  stat = fme->make_facet_lump(shlist, lump);
  if (!lump || stat != CUBIT_SUCCESS) return false;

  DLIList<Lump*> llist;
  BodySM* bodySM;
  Body* body;
  llist.append(lump);
  stat = fme->make_facet_body(llist, bodySM);
  body = GeometryQueryTool::instance()->make_Body(bodySM);
  if (!body || stat != CUBIT_SUCCESS) return false;

  model.append(bodySM);
  return true;
}

int main(int argc, char* argv[])
{
  // Initialize CGM
  FacetQueryEngine* fqe = FacetQueryEngine::instance();

  DLIList<TopologyBridge*> model;
  const char *arg = "multifaceted-2trisPerSurface-works.stl";
  if (generateFacetModel(model))
    {
    ModelExportOptions opts;
    fqe->export_solid_model(
      model, arg, FACET_TYPE,
      CubitString(), opts);
    std::cout << "Wrote \"" << arg << "\"\n";
    }

  return 0;
}
