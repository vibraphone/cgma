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

// forward declare some functions used and defined later
bool generateFacetModel(DLIList<TopologyBridge*>& model);

int main(int argc, char* argv[])
{
  // Initialize CGM
  FacetQueryEngine* fqe = FacetQueryEngine::instance();

  DLIList<TopologyBridge*> model;
  const char *arg = "multifaceted.stl";
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

bool generateFacetModel(
  DLIList<TopologyBridge*>& model)
{
  FacetModifyEngine* fme = FacetModifyEngine::instance();

  CubitStatus stat;
  int i;

  DLIList <CubitPoint*> plist;
  double pts[4][3] = {
      { 0., 0., 0. },
      { 1., 0., 0. },
      { 1., 1., 0. },
      { 0., 1., 1. },
  };
  for (i = 0; i < (int)(sizeof(pts)/sizeof(pts[0])); ++i)
    {
    plist.append(new CubitPointData(pts[i][0], pts[i][1], pts[i][2]));
    }

  // Build a tetrahedron, with 1 triangular facet per faceted-surface:
  int tconn[4][3] = {
      {0, 2, 1},
      {0, 3, 2},
      {0, 1, 3},
      {1, 2, 3}
  };
  DLIList<Surface*> all_surfs;
  DLIList<CubitFacet*> tflist;
  DLIList<CubitQuadFacet*> qflist;
  DLIList<CubitPoint*> facet_plist;
  for (i = 0; i < (int)(sizeof(tconn)/sizeof(tconn[0])); ++i)
    {
    std::cout << "Build faceted surface " << i << "\n";
    //DLIList<CubitFacet*> tflist;
    //DLIList<CubitQuadFacet*> qflist;
    //DLIList<CubitPoint*> facet_plist;
    // Fails with or without the next 3 lines:
    //facet_plist.append(plist[tconn[i][0]]);
    //facet_plist.append(plist[tconn[i][1]]);
    //facet_plist.append(plist[tconn[i][2]]);
    // Create a single facet for our faceted surface:
    tflist.append(
      new CubitFacetData(
        plist[tconn[i][0]], plist[tconn[i][1]], plist[tconn[i][2]]));
    }
  DLIList<Surface*> slist;
  // Attempt to build a faceted surface:
  stat = fme->build_facet_surface(
    qflist, tflist, facet_plist, -1., 0, false, false, slist);
  if (slist.size() <= 0 || stat != CUBIT_SUCCESS)
  {
    PRINT_ERROR("Make surface failed. \n");
    return false;
  }
  all_surfs += slist;

  ShellSM* shell;
  stat = fme->make_facet_shell(all_surfs, shell);
  if (!shell || stat != CUBIT_SUCCESS) 
  {
    PRINT_ERROR("Make shell failed. \n");
    return false;
  }

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

