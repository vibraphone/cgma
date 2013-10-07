
#undef NDEBUG
#include <cassert>

#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "Body.hpp"
#include "AppUtil.hpp"
#include "InitCGMA.hpp"

#ifndef TEST_ENGINE
# define TEST_ENGINE 0
#endif

#ifdef HAVE_OCC
#define export_b  0
#else
#define export_b  1
#endif 

#define BS 1.0
#define RAD 0.4
#define GQI GeometryQueryTool::instance()
#define GMI GeometryModifyTool::instance()

int main( int argc, char* argv[] )
{
    // Start up CGM
  CubitStatus result = InitCGMA::initialize_cgma(TEST_ENGINE);
  if (CUBIT_SUCCESS != result) return 1;

    // Create a brick
  Body* brick = GMI->brick(BS, BS, BS);
  assert(brick != 0);
  DLIList<Body *> bodies, single_body, all_bodies, neighbor_list, new_bodies;
  all_bodies.append(brick);
  for (int i = -1; i <= 1; i+= 2) {
    for (int j = -1; j <= 1; j+= 2) {
      for (int k = -1; k <= 1; k+= 2) {
        Body* sph = GMI->sphere(RAD);
        assert(brick != 0);
        GQI->translate(sph, CubitVector(i*.5*BS, j*.5*BS, k*.5*BS));
        bodies.append(sph);
      }
    }
  }
  CubitStatus stat = GMI->unite(bodies, single_body);
  assert(CUBIT_SUCCESS == stat && single_body.size() == 1);
  stat = GMI->webcut_with_body(all_bodies, single_body.get(), new_bodies, neighbor_list);
  assert(CUBIT_SUCCESS == stat);
  std::cout << "Number of resulting bodies = " << new_bodies.size() << std::endl;
  DLIList<RefEntity*> re_list;
  for (int i = 0; i < new_bodies.size(); i++) re_list.append(new_bodies.get_and_step());

if(export_b == 1)
{
  int num_exp;
  CubitString vers;
  ModelExportOptions opts;
  stat = GQI->export_solid_model (re_list, "spheres.sat", ACIS_SAT_TYPE, num_exp, vers, opts);
}

  bodies.clean_out();
  GQI->bodies(bodies);
  assert (bodies.size() == 10);
  //delete all entities
  GQI->delete_Body(bodies);

  DLIList<RefEntity*> free_entities;
  GQI->get_free_ref_entities(free_entities);
  assert(free_entities.size() ==0);
  return 0;
}

