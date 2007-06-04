/**
 * \file testgeom.cpp
 *
 * \brief testgeom, a unit test for the TSTT geometry interface
 *
 */
#include "TSTTG_CGM.hh"
#include "TSTTG.hh"
#include "TSTTB.hh"
#include <iostream>
#include <set>
#include "TSTTB_SNL_SIDL_defs.h"

// define IMPLEMENTATION_CLASS to be the namespace::class of your implementation
#define IMPLEMENTATION_CLASS TSTTG_CGM::CgmGeom

typedef void* TagHandle;
typedef void* GentityHandle;
typedef void* GentitysetHandle;

bool gLoad_test(const std::string filename, TSTTG::Geometry &geom);

bool tags_test(TSTTG::Geometry &geom);
bool tag_get_set_test(TSTTG::Geometry &geom);
bool tag_info_test(TSTTG::Geometry &geom);
bool gentityset_test(TSTTG::Geometry &geom, bool /*multiset*/, bool /*ordered*/);
bool topology_adjacencies_test(TSTTG::Geometry &geom);
bool construct_test(TSTTG::Modify &geom);
bool primitives_test(TSTTG::Modify &geom);
bool transforms_test(TSTTG::Modify &geom);
bool booleans_test(TSTTG::Modify &geom);

void handle_error_code(const bool result,
                       int &number_failed,
                       int &/*number_not_implemented*/,
                       int &number_successful)
{
  if (result) {
    std::cout << "Success";
    number_successful++;
  }
  else {
    std::cout << "Failure";    
    number_failed++;
  }
}

#define CAST_TSTT_INTERFACE(var_in, var_out, iface, rvalue) \
          TSTTB::iface var_out; \
          try {var_out = var_in;}\
          catch(TSTTB::Error) {\
            std::cerr << "Error: current interface doesn't support iface." << std::endl; \
            return rvalue;}

#define CAST_GEOM_INTERFACE(var_in, var_out, iface, rvalue) \
          TSTTG::iface var_out; \
          try {var_out = var_in;}\
          catch(TSTTB::Error) {\
            std::cerr << "Error: current interface doesn't support iface." << std::endl; \
            return rvalue;}

int main( int argc, char *argv[] )
{
    // Check command line arg
  std::string filename;

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <geom_filename>" << std::endl;
    return 1;
  }
  filename = argv[1];

  bool result;
  int number_tests = 0;
  int number_tests_successful = 0;
  int number_tests_not_implemented = 0;
  int number_tests_failed = 0;

    // initialize the Mesh
  TSTTG::Geometry geom = IMPLEMENTATION_CLASS::_create();

    // Print out Header information
  std::cout << "\n\nTSTT GEOMETRY INTERFACE TEST PROGRAM:\n\n";

    // gLoad test
  std::cout << "   gLoad: ";
  result = gLoad_test(filename, geom);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  std::cout << "\n";

    // tags test
  std::cout << "   tags: ";
  result = tags_test(geom);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  std::cout << "\n";

    // gentitysets test
  std::cout << "   gentity sets: ";
  result = gentityset_test(geom, false, false);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  std::cout << "\n";

    // gentitysets test
  std::cout << "   topology adjacencies: ";
  result = topology_adjacencies_test(geom);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  std::cout << "\n";

  CAST_GEOM_INTERFACE(geom, geom_modify, Modify, false);
    // construct test
  std::cout << "   construct: ";
  result = construct_test(geom_modify);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  std::cout << "\n";

    // primitives test
  std::cout << "   primitives: ";
  result = primitives_test(geom_modify);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  std::cout << "\n";

    // transforms test
  std::cout << "   transforms: ";
  result = transforms_test(geom_modify);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  std::cout << "\n";

    // booleans test
  std::cout << "   booleans: ";
  result = booleans_test(geom_modify);
  handle_error_code(result, number_tests_failed,
                    number_tests_not_implemented,
                    number_tests_successful);
  number_tests++;
  std::cout << "\n";

    // summary

  std::cout << "\nTSTT TEST SUMMARY: \n"
            << "   Number Tests:           " << number_tests << "\n"
            << "   Number Successful:      " << number_tests_successful << "\n"
            << "   Number Not Implemented: " << number_tests_not_implemented << "\n"
            << "   Number Failed:          " << number_tests_failed 
            << "\n\n" << std::endl;
  
  return 0;
}

/*!
  @test 
  Load Mesh
  @li Load a mesh file
*/
bool gLoad_test(const std::string filename, TSTTG::Geometry &geom)
{
  CAST_GEOM_INTERFACE(geom, geom_cq, CoreQuery, false);
  CAST_GEOM_INTERFACE(geom, geom_topo, Topology, false);
  
    // load a mesh
  sidl::array<std::string> options;
  options = options.create1d(1);
  
  try {
    geom_cq.gLoad(filename, options, 1);
  }
  catch (TSTTB::Error) {
    std::cerr << "ERROR : can not load a geometry" << std::endl;
    return false;
  }

    // print out the number of entities
  std::cout << "Model contents: " << std::endl;
  const char *gtype[] = {"Gvertices: ", "Gedges: ", "Gfaces: ", "Gregions: "};
  try {
    for (int i = 0; i <= 3; i++)
      std::cout << gtype[i] 
                << geom_topo.gentitysetGetNumberGentitiesOfType(0, (TSTTG::GentityType)i) 
                << std::endl;
  } catch (TSTTB::Error) {
    std::cerr << "Error: problem getting entities after gLoad." 
              << std::endl;
    return false;
  }

  return true;
}

/*!
  @test 
  Test tag creating, reading, writing, deleting
  @li Load a mesh file
*/
bool tags_test(TSTTG::Geometry &geom)
{
  CAST_TSTT_INTERFACE(geom, geom_tag, Tag, false);

  bool success = true;

  success = tag_info_test(geom);
  if (!success) return success;
  
  success = tag_get_set_test(geom);
  if (!success) return success;

  return true;
}

bool tag_info_test(TSTTG::Geometry &geom) 
{
  CAST_TSTT_INTERFACE(geom, geom_tag, Tag, false);
  CAST_TSTT_INTERFACE(geom, geom_etag, EntTag, false);
  //CAST_TSTT_INTERFACE(geom, geom_atag, ArrTag, false);

    // create an arbitrary tag, size 4
  TagHandle this_tag, tmp_handle;
  std::string tag_name("tag_info tag"), tmp_name;
  try {
    geom_tag.createTag(tag_name, 4, TSTTB::TagValueType_BYTES, this_tag);
  }
  catch (TSTTB::Error) {
    std::cerr << "ERROR : can not create a tag." << std::endl;
    return false;
  }

  bool success = true;
  
    // get information on the tag
  try {
    tmp_name = geom_tag.getTagName(this_tag);
  }
  catch (TSTTB::Error) {
    std::cerr << "ERROR : Couldn't get tag name." << std::endl;
    success = false;
  }
  if (tmp_name != tag_name) {
    std::cerr << "ERROR: getTagName didn't return consistent result." << std::endl;
    success = false;
  }

  try {
    tmp_handle = geom_tag.getTagHandle(tag_name);
  }
  catch (TSTTB::Error) {
    std::cerr << "ERROR : Couldn't get tag handle." << std::endl;
    success = false;
  }
  if (tmp_handle != this_tag) {
    std::cerr << "ERROR: getTagHandle didn't return consistent result." << std::endl;
    success = false;
  }

  int tag_size;
  try {
    tag_size = geom_tag.getTagSizeBytes(this_tag);
  }
  catch (TSTTB::Error) {
    std::cerr << "ERROR : Couldn't get tag size." << std::endl;
    success = false;
  }
  if (tmp_handle != this_tag) {
    std::cerr << "ERROR: getTagSize didn't return consistent result." << std::endl;
    success = false;
  }

  try {
    geom_tag.destroyTag(this_tag, true);
  }
  catch (TSTTB::Error) {
    std::cerr << "ERROR : Couldn't delete a tag." << std::endl;
    success = false;
  }

    // print information about all the tags in the model
  std::set<std::string> tag_names;
  CAST_GEOM_INTERFACE(geom_tag, geom_topo, Topology, false);
  try {
    for (int dim = 0; dim <= 3; dim++) {
      sidl::array<GentityHandle> gentity_handles;
      int num_gentity_handles;
        // get all entities of this dimension
      geom_topo.gentitysetGetGentitiesOfType(0, (TSTTG::GentityType)dim, gentity_handles,
                                             num_gentity_handles);

      sidl::array<void*> tag_handles;
      int num_tag_handles;
      int num_ents = ARRAY_SIZE(gentity_handles);
        // for each entity of this dimension
      for (int i = 0; i < num_ents; i++) {
	// get the tags on this entity
        //geom_tag.gentityGetAllTagHandles(gentity_handles.get(i), tag_handles,
	//                               num_tag_handles);
	geom_etag.getAllTags(gentity_handles.get(i), tag_handles,
			     num_tag_handles);
        int num_tags = ARRAY_SIZE(tag_handles);
          // for each tag, get the name and add to the set
        for (int j = 0; j < num_tags; j++) {
          std::string this_name = geom_tag.getTagName(tag_handles.get(j));
          tag_names.insert(this_name);
        }
      }
    }
  }
  catch (TSTTB::Error) {
    std::cerr << "ERROR : Trouble getting all tags in model." << std::endl;
    success = false;
  }

  if (tag_names.empty())
    std::cout << "No tags defined on model." << std::endl;
  else {
    std::cout << "Tags defined on model: ";
    bool first = true;
    for (std::set<std::string>::iterator sit = tag_names.begin(); sit != tag_names.end(); sit++) {
      if (!first) std::cout << ", ";
      std::cout << *sit;
      first = false;
    }
    std::cout << std::endl;
  }

  return success;
}

bool tag_get_set_test(TSTTG::Geometry &geom) 
{
  CAST_TSTT_INTERFACE(geom, geom_tag, Tag, false);
  CAST_GEOM_INTERFACE(geom, geom_topo, Topology, false);
  CAST_TSTT_INTERFACE(geom, geom_atag, ArrTag, false);

    // create an arbitrary tag, size 4
  TagHandle this_tag;
  std::string tag_name("tag_get_set tag");
  try {
    geom_tag.createTag(tag_name, 4, TSTTB::TagValueType_BYTES, this_tag);
  }
  catch (TSTTB::Error) {
    std::cerr << "ERROR : can not create a tag for get_set test." << std::endl;
    return false;
  }

  bool success = true;
  
    // set this tag to an integer on each entity; keep track of total sum
  int sum = 0, num = 0, dim;
  for (dim = 0; dim <= 3; dim++) {
    sidl::array<GentityHandle> gentity_handles;
    int num_gentity_handles;
    try {
      geom_topo.gentitysetGetGentitiesOfType(0, (TSTTG::GentityType)dim, gentity_handles,
                                             num_gentity_handles);
    }
    catch (TSTTB::Error) {
      std::cerr << "ERROR : can't get entities of dimension " << dim << "." << std::endl;
      success = false;
      continue;
    }

    int num_ents = ARRAY_SIZE(gentity_handles);
    int tag_vals_size = num_ents*sizeof(int);
    sidl::array<char> tag_vals = tag_vals.create1d(tag_vals_size);
    int *tag_ptr = ARRAY_PTR(tag_vals, int);
    for (int i = 0; i < num_ents; i++) {
      tag_ptr[i] = num;
      sum += num++;
    }
      
    try {
      geom_atag.setArrData(gentity_handles, num_gentity_handles,
                           this_tag, tag_vals, tag_vals_size);
    }
    catch (TSTTB::Error) {
      std::cerr << "ERROR : can't set tag on entities of dimension " << dim 
                << "." << std::endl;
      success = false;
    }
  }
  
    // check tag values for entities now
  int get_sum = 0;
  for (dim = 0; dim <= 3; dim++) {
    sidl::array<GentityHandle> gentity_handles;
    int num_gentity_handles;
    try {
      geom_topo.gentitysetGetGentitiesOfType(0, (TSTTG::GentityType)dim, gentity_handles,
                                             num_gentity_handles);
    }
    catch (TSTTB::Error) {
      std::cerr << "ERROR : can't get entities of dimension " << dim << "." << std::endl;
      success = false;
      continue;
    }

    int num_ents = ARRAY_SIZE(gentity_handles);
    sidl::array<char> tag_vals;
    int tag_vals_size;
    try {
      geom_atag.getArrData(gentity_handles, num_gentity_handles,
                           this_tag, tag_vals, tag_vals_size);
    }
    catch (TSTTB::Error) {
      std::cerr << "ERROR : can't get tag on entities of dimension " << dim 
                << "." << std::endl;
      success = false;
    }

    int *tag_ptr = ARRAY_PTR(tag_vals, int);
    for (int i = 0; i < num_ents; i++) {
      get_sum += tag_ptr[i];
    }
  }
  if (get_sum != sum) {
    std::cerr << "ERROR: getData didn't return consistent results." << std::endl;
    success = false;
  }
  
  try {
    geom_tag.destroyTag(this_tag, true);
  }
  catch (TSTTB::Error) {
    std::cerr << "ERROR : couldn't delete tag." << std::endl;
    success = false;
  }

  return success;
}

/*!
  @test
  TSTT gentity sets test (just implemented parts for now)
  @li Check gentity sets
*/
bool gentityset_test(TSTTG::Geometry &geom, bool /*multiset*/, bool /*ordered*/)
{
  int num_type = 4;
  int num_all_gentities_super = 0;
  GentitysetHandle ges_array[num_type];
  int number_array[num_type];
  int ent_type = TSTTG::GentityType_GVERTEX;
  //CAST_INTERFACE(geom, geom_cgeso, CoreGentitysetOperations, false);
  CAST_TSTT_INTERFACE(geom, geom_eset, EntSet, false);
  CAST_GEOM_INTERFACE(geom, geom_topo, Topology, false);
  //CAST_INTERFACE(geom, geom_bso, BooleanSetOperations, false);
  CAST_TSTT_INTERFACE(geom, geom_sbo, SetBoolOps, false);
  CAST_TSTT_INTERFACE(geom, geom_srel, SetRelation, false);
  
  CAST_GEOM_INTERFACE(geom, geom_cq, CoreQuery, false);
  //CAST_INTERFACE(geom, geom_gesr, GentitysetRelations, false);

  bool success = true;
  
    // get the number of sets in the whole model
  int all_sets = 0;
  try {
    all_sets = geom_eset.getNumEntSets(0, 0);
  } catch (TSTTB::Error) {
    std::cerr << "Problem getting the number of all gentity sets in whole model."
              << std::endl;
    success = false;
  }

    // add gentities to entitysets by type
  for (; ent_type < num_type; ent_type++) {
      // initialize the entityset
    try {
      geom_eset.createEntSet(true, ges_array[ent_type]);
    } catch (TSTTB::Error) {
      std::cerr << "Problem creating entityset." << std::endl;
      success = false;
      continue;
    }

      // get entities by type in total "mesh"
    sidl::array<GentityHandle> gentities;
    int num_gentities;
    try {
      geom_topo.gentitysetGetGentitiesOfType
        (NULL, static_cast<TSTTG::GentityType>(ent_type),
         gentities, num_gentities);
    } catch (TSTTB::Error err) {
      std::cerr << "Failed to get gentities by type in gentityset_test."
                << std::endl;
      success = false;
      continue;
    }

      // add gentities into gentity set
    try {
      geom_eset.addEntArrToSet(gentities, num_gentities, ges_array[ent_type]);
    } catch (TSTTB::Error) {
      std::cerr << "Failed to add gentities in entityset_test."
                << std::endl;
      success = false;
      continue;
    }

      // Check to make sure entity set really has correct number of entities in it
    try {
      number_array[ent_type] = 
        geom_topo.gentitysetGetNumberGentitiesOfType(ges_array[ent_type],
                                                     static_cast<TSTTG::GentityType>(ent_type));
      
    } catch (TSTTB::Error) {
      std::cerr << "Failed to get number of gentities by type in entityset_test."
                << std::endl;
      return false;
    }  

      // compare the number of entities by type
    int num_type_gentity = ARRAY_SIZE(gentities);

    if (number_array[ent_type] != num_type_gentity)
    {
      std::cerr << "Number of gentities by type is not correct"
                << std::endl;
      success = false;
      continue;
    }

      // add to number of all entities in super set
    num_all_gentities_super += num_type_gentity;
  }

    // make a super set having all entitysets
  GentitysetHandle super_set = NULL;
  try {
    geom_eset.createEntSet(true, super_set);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to create a super set in gentityset_test."
              << std::endl;
    success = false;
  }

  for (int i = 0; i < num_type; i++) {
    try {
      geom_eset.addEntSet(ges_array[i], super_set);
    } catch (TSTTB::Error) {
      std::cerr << "Failed to create a super set in gentityset_test."
		<< std::endl;
      success = false;
    }
  }


    //----------TEST BOOLEAN OPERATIONS----------------//

  GentitysetHandle temp_ges1;
  try {
    geom_eset.createEntSet(true, temp_ges1);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to create a super set in gentityset_test."
              << std::endl;
    success = false;
  }

    // Subtract
    // add all EDGEs and FACEs to temp_es1
    // get all EDGE entities
  sidl::array<GentityHandle> gedges;
  int num_gedges;
  sidl::array<GentityHandle> gfaces;
  int num_gfaces;
  sidl::array<GentityHandle> temp_gentities1;
  int num_temp_gentities1;

  try {
    geom_topo.gentitysetGetGentitiesOfType
      (ges_array[TSTTG::GentityType_GEDGE], TSTTG::GentityType_GEDGE, gedges, num_gedges);

  } catch (TSTTB::Error err) {
    std::cerr << "Failed to get gedge gentities in gentityset_test."
              << std::endl;
    success = false;
  }

    // add GEDGEs to ges1
  try {
    geom_eset.addEntArrToSet(gedges, num_gedges, temp_ges1);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to add gedge gentities in gentityset_test."
              << std::endl;
    success = false;
  }

    // get all GFACE gentities
  try {
    geom_topo.gentitysetGetGentitiesOfType
      (ges_array[TSTTG::GentityType_GFACE], TSTTG::GentityType_GFACE, gfaces, num_gfaces);
  } catch (TSTTB::Error err) {
    std::cerr << "Failed to get gface gentities in gentityset_test."
              << std::endl;
    success = false;
  }

    // add FACEs to es1
  try {
    geom_eset.addEntArrToSet(gfaces, num_gfaces, temp_ges1);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to add gface gentities in gentityset_test."
              << std::endl;
    success = false;
  }

    // subtract GEDGEs
  try {
    geom_sbo.subtract(temp_ges1, ges_array[TSTTG::GentityType_GEDGE], temp_ges1);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to subtract gentitysets in gentityset_test."
              << std::endl;
    success = false;
  }

  try {
    geom_topo.gentitysetGetGentitiesOfType
      (temp_ges1, TSTTG::GentityType_GFACE, temp_gentities1, num_temp_gentities1);
  } catch (TSTTB::Error err) {
    std::cerr << "Failed to get gedge gentities in gentityset_test."
              <<std::endl;
    success = false;
  }

  if (ARRAY_SIZE(gfaces) != ARRAY_SIZE(temp_gentities1)) {
    std::cerr << "Number of entitysets after subtraction not correct \
             in gentityset_test." << std::endl;
    success = false;
  }

    // check there's nothing but gfaces in temp_ges1
  int num_gents;
  try {
    num_gents = geom_topo.gentitysetGetNumberGentitiesOfType(temp_ges1, TSTTG::GentityType_GFACE);
  } catch (TSTTB::Error err) {
    std::cerr << "Failed to get dimensions of gentities in gentityset_test." << std::endl;
    success = false;
  }
  if (num_gents != ARRAY_SIZE(temp_gentities1)) {
    std::cerr << "Wrong number of faces in subtract test." << std::endl;
    success = false;
  }

    //------------Intersect------------
    //

    // clean out the temp_ges1
  try {
    geom_eset.rmvEntArrFromSet(gfaces, num_gfaces, temp_ges1);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to remove gface gentities in gentityset_test."
              << std::endl;
    success = false;
  }

    // check if it is really cleaned out
  try {
    num_gents = geom_topo.gentitysetGetNumberGentitiesOfType(temp_ges1, TSTTG::GentityType_GFACE);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to get number of gentities by type in gentityset_test."
              << std::endl;
    success = false;
  }

  if (num_gents != 0) {
    std::cerr << "failed to remove correctly." << std::endl;
    success = false;
  }
  
    // add GEDGEs to temp ges1
  try {
    geom_eset.addEntArrToSet(gedges, num_gedges, temp_ges1);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to add gedge gentities in gentityset_test."
              << std::endl;
    success = false;
  }

    // add GFACEs to temp ges1
  try {
    geom_eset.addEntArrToSet(gfaces, num_gfaces, temp_ges1);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to add gface gentities in gentityset_test."
              << std::endl;
    success = false;
  }

    // intersect temp_ges1 with gedges set 
    // temp_ges1 entityset is altered
  try {
    geom_sbo.intersect(temp_ges1, ges_array[TSTTG::GentityType_GEDGE], temp_ges1);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to intersect in gentityset_test."
              << std::endl;
    success = false;
  }

    // try to get GFACEs, but there should be nothing but GEDGE
  try {
    num_gents = geom_topo.gentitysetGetNumberGentitiesOfType
      (temp_ges1, TSTTG::GentityType_GFACE);
  } catch (TSTTB::Error err) {
    std::cerr << "Failed to get gface gentities in gentityset_test."
              <<std::endl;
    success = false;
  }

  if (num_gents != 0) {
    std::cerr << "wrong number of gfaces." << std::endl;
    success = false;
  }


    //-------------Unite--------------

    // get all regions
  GentitysetHandle temp_ges2;
  sidl::array<GentityHandle> gregions;
  int num_gregions;

  try {
    geom_eset.createEntSet(true, temp_ges2);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to create a temp gentityset in gentityset_test."
              << std::endl;
    success = false;
  }

  try {
    geom_topo.gentitysetGetGentitiesOfType
      (ges_array[TSTTG::GentityType_GREGION], TSTTG::GentityType_GREGION, gregions, num_gregions);
  } catch (TSTTB::Error err) {
    std::cerr << "Failed to get gregion gentities in gentityset_test." << std::endl;
    success = false;
  }

    // add GREGIONs to temp es2
  try {
    geom_eset.addEntArrToSet(gregions, num_gregions, temp_ges2);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to add gregion gentities in gentityset_test." << std::endl;
    success = false;
  }

    // unite temp_ges1 and temp_ges2
    // temp_ges1 gentityset is altered
  try {
    geom_sbo.unite(temp_ges1, temp_ges2, temp_ges1);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to unite in gentityset_test." << std::endl;
    success = false;
  }

    // perform the check
  try {
    num_gents =
      geom_topo.gentitysetGetNumberGentitiesOfType(temp_ges1, TSTTG::GentityType_GREGION);
  } catch (TSTTB::Error) {
    std::cerr << "Failed to get number of gregion gentities by type in gentityset_test."
              << std::endl;
    success = false;
  }
  
  if (num_gents != number_array[TSTTG::GentityType_GREGION]) {
    std::cerr << "different number of gregions in gentityset_test." << std::endl;
    success = false;
  }


    //--------Test parent/child stuff in entiysets-----------

    // Add 2 sets as children to another
  GentitysetHandle parent_child;
  try {
    geom_eset.createEntSet(true, parent_child);
  } catch (TSTTB::Error) {
    std::cerr << "Problem creating gentityset in gentityset_test."
              << std::endl;
    success = false;
  }

  try {
    geom_srel.addPrntChld(ges_array[TSTTG::GentityType_GVERTEX], parent_child);
  } catch (TSTTB::Error) {
    std::cerr << "Problem add parent in gentityset_test."
              << std::endl;
    success = false;
  }

    // check if parent is really added
  sidl::array<void*> parents;
  int num_parents;
  try {
    geom_srel.getPrnts(parent_child, 1, parents, num_parents);
  } catch (TSTTB::Error) {
    std::cerr << "Problem getting parents in gentityset_test."
              << std::endl;
    success = false;
  }

  if (ARRAY_SIZE(parents) != 1) {
    std::cerr << "number of parents is not correct in gentityset_test."
              << std::endl;
    success = false;
  }

    // add parent and child
  //sidl::array<void*> parent_child_array = sidl::array<void*>::create1d(1);
  //int num_parent_child_array;
  //sidl::array<void*> temp_gedge_array = sidl::array<void*>::create1d(1);
  //int num_temp_gedge_array;
  //parent_child_array.set(0, parent_child);
  //temp_gedge_array.set(0, ges_array[TSTTG::GentityType_GEDGE]);
  try {
    geom_srel.addPrntChld(ges_array[TSTTG::GentityType_GEDGE], parent_child);
  } catch (TSTTB::Error) {
    std::cerr << "Problem adding parent and child in gentityset_test."
              << std::endl;
    success = false;
  }

  //sidl::array<void*> temp_gface_array = sidl::array<void*>::create1d(1);
  //int num_temp_gface_array;
  //temp_gface_array.set(0, ges_array[TSTTG::GentityType_GFACE]);
  try {
    geom_srel.addPrntChld(parent_child, ges_array[TSTTG::GentityType_GFACE]);
  } catch (TSTTB::Error) {
    std::cerr << "Problem adding parent and child in gentityset_test."
              << std::endl;
    success = false;
  }

    // add child
  try {
    geom_srel.addPrntChld(parent_child, ges_array[TSTTG::GentityType_GREGION]);
  } catch (TSTTB::Error) {
    std::cerr << "Problem adding child in gentityset_test."
              << std::endl;
    success = false;
  }

    // get the number of parent gentitysets
  num_gents = -1;
  try {
    num_gents = geom_srel.getNumPrnt(parent_child, 1);
  } catch (TSTTB::Error) {
    std::cerr << "Problem getting number of parents in gentityset_test."
              << std::endl;
    success = false;
  }

  if (num_gents != 2) {
    std::cerr << "number of parents is not correct in gentityset_test."
              << std::endl;
    success = false;
  }

    // get the number of child gentitysets
  try {
    num_gents = geom_srel.getNumChld(parent_child, 1);
  } catch (TSTTB::Error) {
    std::cerr << "Problem getting number of children in gentityset_test."
              << std::endl;
    success = false;
  }

  if (num_gents != 2) {
    std::cerr << "number of children is not correct in gentityset_test."
              << std::endl;
    success = false;
  }

  sidl::array<void*> children;
  int num_children;
  try {
    geom_srel.getChldn(parent_child, 1, children, num_children);
  } catch (TSTTB::Error) {
    std::cerr << "Problem getting children in gentityset_test."
              << std::endl;
    success = false;
  }

  if (ARRAY_SIZE(children) != 2) {
    std::cerr << "number of children is not correct in gentityset_test."
              << std::endl;
    success = false;
  }

    // remove children
  try {
    geom_srel.rmvPrntChld(parent_child, ges_array[TSTTG::GentityType_GFACE]);
  } catch (TSTTB::Error) {
    std::cerr << "Problem removing parent child in gentityset_test."
              << std::endl;
    success = false;
  }

    // get the number of child gentitysets
  try {
    num_gents = geom_srel.getNumChld(parent_child, 1);
  } catch (TSTTB::Error) {
    std::cerr << "Problem getting number of children in gentityset_test."
              << std::endl;
    success = false;
  }

  if (num_gents != 1) {
    std::cerr << "number of children is not correct in gentityset_test."
              << std::endl;
    success = false;
  }

    // parent_child and ges_array[TSTTG::GentityType_GEDGE] should be related
  try {
    if (!geom_srel.isChildOf
	(ges_array[TSTTG::GentityType_GEDGE], parent_child)) {
      std::cerr << "parent_child and ges_array[TSTTG::GentityType_GEDGE] should be related" << std::endl;
      success = false;
    }
  } catch (TSTTB::Error) {
    std::cerr << "Problem checking relation in gentityset_test."
              << std::endl;
    success = false;
  }

    // ges_array[TSTTG::GentityType_GFACE] and ges_array[TSTTG::REGION] are not related
  try {
    if (geom_srel.isChildOf
        (ges_array[TSTTG::GentityType_GFACE], ges_array[TSTTG::GentityType_GREGION])) {
      std::cerr << "ges_array[TSTTG::REGION] and ges_array[TSTTG::GentityType_GFACE] should not be related" << std::endl;
      success = false;
    }
  } catch (TSTTB::Error) {
    std::cerr << "Problem checking relation in gentityset_test."
              << std::endl;
    success = false;
  }
  

    //--------test modify and query functions-----------------------------
  
    // check the number of gentity sets in whole mesh
  sidl::array<void*> gentity_sets;
  int num_gentity_sets;
  try {
    geom_eset.getEntSets(NULL, 1, gentity_sets, num_gentity_sets);
  } catch (TSTTB::Error) {
    std::cerr << "Problem to get all gentity sets in mesh."
              << std::endl;
    success = false;
  }

  if (ARRAY_SIZE(gentity_sets) != all_sets + 8) {
    std::cerr << "the number of gentity sets in whole mesh should be 8 times of num_iter."
              << std::endl;
    success = false;
  }

    // get all gentity sets in super set
  sidl::array<void*> ges_array1;
  int num_ges_array1;
  try {
    geom_eset.getEntSets(super_set, 1, ges_array1, num_ges_array1);
  } catch (TSTTB::Error) {
    std::cerr << "Problem to get gentity sets in super set."
              << std::endl;
    success = false;
  }

  int num_super;

    // get the number of gentity sets in super set
  try {
    num_super = geom_eset.getNumEntSets(super_set, 1);
  } catch (TSTTB::Error) {
    std::cerr << "Problem to get the number of all gentity sets in super set."
              << std::endl;
    success = false;
  }

    // the number of gentity sets in super set should be same
  if (num_super != ARRAY_SIZE(ges_array1)) {
    std::cerr << "the number of gentity sets in super set should be same." << std::endl;
    success = false;
  }

    // get all entities in super set
  sidl::array<void*> all_gentities;
  int num_all_gentities = 0;
  try {
    geom_eset.getEntSets(super_set, 1, all_gentities, num_all_gentities);
//    geom_eset.getEntSets(super_set, 2, all_gentities, num_all_gentities);
  } catch (TSTTB::Error) {
    std::cerr << "Problem to get all gentities in super set."
              << std::endl;
    success = false;
  }
  
    // compare the number of all gentities in super set
  // HJK : num_hops is not implemented
  //if (num_all_gentities_super != ARRAY_SIZE(all_gentities)) {
  //std::cerr << "number of all gentities in super set should be same." << std::endl;
  //success = false;
  //}

    // test add, remove and get all entitiy sets using super set
    // check GetAllGentitysets works recursively and dosen't return
    // multi sets
  for (int k = 0; k < num_super; k++) {
      // add gentity sets of super set to each gentity set of super set
      // make multiple child super sets
    GentitysetHandle ges_k = ges_array1.get(k);

    for (int a = 0; a < num_ges_array1; a++) {
      try {
	geom_eset.addEntSet(ges_array1[a], ges_k);
      } catch (TSTTB::Error) {
	std::cerr << "Problem to add entity set."
		  << std::endl;
	success = false;
      }
    }

      // add super set to each entity set
    //    sidl::array<GentitysetHandle> superset_array
    //= sidl::array<GentitysetHandle>::create1d(1);
    //superset_array.set(0, super_set);
    //int num_superset_array;
    
    try {
      geom_eset.addEntSet(super_set, ges_k);
    } catch (TSTTB::Error) {
      std::cerr << "Problem to add super set to gentitysets." << std::endl;
      success = false;
    }

      // add one gentity sets multiple times
    // HJK: ??? how to deal this case?
    //sidl::array<GentitysetHandle> temp_array1
    //= sidl::array<GentitysetHandle>::create1d(1);
    //int num_temp_array1;
    //temp_array1.set(0, temp_ges1);

    //for (int l = 0; l < 3; l++) {
      try {
	geom_eset.addEntSet(temp_ges1, ges_k);
      } catch (TSTTB::Error) {
	std::cerr << "Problem to add temp set to gentitysets." << std::endl;
	success = false;
      }
      //}
  }

  return success;
}
  
/*!
@test
TSTTG topology adjacencies Test
@li Check topology information
@li Check adjacency
*/
// make each topological entity vectors, check their topology
// types, get interior and exterior faces of model
bool topology_adjacencies_test(TSTTG::Geometry &geom)
{
  CAST_GEOM_INTERFACE(geom, geom_topo, Topology, false);
  CAST_GEOM_INTERFACE(geom, geom_cq, CoreQuery, false);

  int top = TSTTG::GentityType_GVERTEX;
  int num_test_top = TSTTG::GentityType_GALL_TYPES;
  std::vector<GentityHandle> **gentity_vectors = 
    new std::vector<GentityHandle> *[num_test_top];

  // make the array of vectors having each topology entities
  int i;
  for (i = top; i < num_test_top; i++)
    gentity_vectors[i] = new std::vector<GentityHandle>;

  // fill the vectors of each topology entities
  // like lines vector, polygon vector, triangle vector,
  // quadrilateral, polyhedrron, tet, hex, prism, pyramid,
  // septahedron vectors
  bool success = true;
  for (i = top; i < num_test_top; i++) {
    sidl::array<GentityHandle> gentities;
    int num_gentities;
    try {
      geom_topo.gentitysetGetGentitiesOfType(NULL, 
                                             static_cast<TSTTG::GentityType>(i),
                                             gentities, num_gentities);
    } catch (TSTTB::Error err) {
      std::cerr << "Failed to get gentities in adjacencies_test." << std::endl;
      success = false;
      continue;
    }
    
    int siz_gentities = ARRAY_SIZE(gentities);
    GentityHandle *handles = ARRAY_PTR(gentities, GentityHandle);
    std::copy(handles, handles+siz_gentities, std::back_inserter(*gentity_vectors[i]));
  }

  // check number of entities for each topology
  for (i = top; i < num_test_top; i++) {
    int num_tops = 0;
    try {
      num_tops = geom_topo.gentitysetGetNumberGentitiesOfType
	(NULL, static_cast<TSTTG::GentityType>(i));
    } catch (TSTTB::Error err) {
      std::cerr << "Failed to get number of gentities in adjacencies_test." << std::endl;
      success = false;
      continue;
    }
    
    if (static_cast<int>(gentity_vectors[i]->size()) != num_tops) {
      std::cerr << "Number of gentities doesn't agree with number returned for dimension " 
                << i << std::endl;
      success = false;
      continue;
    }
  }

  // check adjacencies in both directions
  try {
    std::vector<GentityHandle>::iterator vit;
    for (i = TSTTG::GentityType_GREGION; i >= TSTTG::GentityType_GVERTEX; i--) {
      for (vit = gentity_vectors[i]->begin(); vit != gentity_vectors[i]->end(); vit++) {
        GentityHandle this_gent = *vit;

          // check downward adjacencies
        for (int j = TSTTG::GentityType_GVERTEX; j < i; j++) {
        
          sidl::array<GentityHandle> lower_ents;
          int num_lower_ents;
          geom_topo.gentityGetAdjacencies(this_gent, 
                                          static_cast<TSTTG::GentityType>(j), 
                                          lower_ents, num_lower_ents);
        
            // for each of them, make sure they are adjacent to the upward ones
          int num_lower = ARRAY_SIZE(lower_ents);
          for (int k = 0; k < num_lower; k++) {
            sidl::array<GentityHandle> upper_ents;
            int num_upper_ents;
            geom_topo.gentityGetAdjacencies(lower_ents.get(k), 
                                            static_cast<TSTTG::GentityType>(i), 
                                            upper_ents, num_upper_ents);
            GentityHandle *upper_ptr = ARRAY_PTR(upper_ents, GentityHandle);
            int num_upper = ARRAY_SIZE(upper_ents);
            if (std::find(upper_ptr, upper_ptr+num_upper, this_gent) ==
                upper_ptr+num_upper) {
              std::cerr << "Didn't find lower-upper adjacency which was supposed to be there, dims = "
                   << i << ", " << j << std::endl;
              success = false;
            }
          }
        }
      }
    }
  } catch (TSTTB::Error err) {
    std::cerr << "Bi-directional adjacencies test failed." << std::endl;
  }

  for (int i = 0; i < num_test_top; i++)
    delete gentity_vectors[i];
  delete [] gentity_vectors;

  return success;
}

/*!
@test
TSTTG construct Test
@li Check construction of geometry
*/
bool construct_test(TSTTG::Modify &geom)
{
    // query ones
  CAST_GEOM_INTERFACE(geom, geom_topo, Topology, false);
  CAST_GEOM_INTERFACE(geom, geom_shape, Shape, false);

    // modify ones
  CAST_GEOM_INTERFACE(geom, geom_construct, Construct, false);
  CAST_GEOM_INTERFACE(geom, geom_prim, Primitives, false);
  CAST_GEOM_INTERFACE(geom, geom_transforms, Transforms, false);

    // construct a cylinder, sweep it about an axis, and delete the result
  GentityHandle cyl = 0;
  try {
    geom_prim.Cylinder(1.0, 1.0, 0.0, cyl);
    
  } catch (TSTTB::Error err) {
    std::cerr << "Creating cylinder failed." << std::endl;
    return false;
  }
  
    // move it onto the y axis
  GentityHandle new_body = 0;
  try {
    geom_transforms.Move(cyl, 0.0, 1.0, -0.5);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems moving surface." << std::endl;
    return false;
  }

      // get the surface with max z
  GentityHandle max_surf = 0;
  try {
    sidl::array<GentityHandle> surfs;
    int num_surfs;
    geom_topo.gentityGetAdjacencies(cyl, TSTTG::GentityType_GFACE, surfs, num_surfs);

    sidl::array<double> max_corn, min_corn;
    int num_max, num_min;
    geom_shape.gentityBoundingBox(surfs, num_surfs, 
                                  min_corn, num_min, max_corn, num_max);
        
      // for each of them, make sure they are adjacent to the upward ones
    num_surfs = ARRAY_SIZE(surfs);
    for (int i = 0; i < num_surfs; i++) {
      if (max_corn.get(3*i+2) >= 0.0 && min_corn.get(3*i+2) >= 0.0) {
        max_surf = surfs[i];
        break;
      }
    }
  } catch (TSTTB::Error err) {
    std::cerr << "Problems getting max surf for rotation." << std::endl;
    return false;
  }
  
  if (0 == max_surf) {
    std::cerr << "Couldn't find max surf for rotation." << std::endl;
    return false;
  }
  
    // sweep it around the x axis
  try {
    geom_construct.SweepAboutAxis(max_surf, 360.0, 1.0, 0.0, 0.0,
                                  new_body);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems sweeping surface about axis." << std::endl;
    return false;
  }
  
    // now delete
  try {
    geom_construct.Delete(new_body);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems deleting cylinder or swept surface body." << std::endl;
    return false;
  }
  
    // if we got here, we were successful
  return true;
}

bool primitives_test(TSTTG::Modify &geom) 
{
    // query ones
  CAST_GEOM_INTERFACE(geom, geom_topo, Topology, false);
  CAST_GEOM_INTERFACE(geom, geom_shape, Shape, false);

    // modify ones
  CAST_GEOM_INTERFACE(geom, geom_construct, Construct, false);
  CAST_GEOM_INTERFACE(geom, geom_prim, Primitives, false);
  CAST_GEOM_INTERFACE(geom, geom_transforms, Transforms, false);

  sidl::array<GentityHandle> prims
    = sidl::array<GentityHandle>::create1d(3);

    // construct a primitive
  int failed_num = 0;
  try {
    for (int i = 0; i < 3; i++) {
      GentityHandle prim = 0;
      switch (i) {
        case 0:
          geom_prim.Brick(1.0, 2.0, 3.0, prim);
          if (0 == prim) failed_num += 1;
          break;
        case 1:
          geom_prim.Cylinder(1.0, 4.0, 2.0, prim);
          if (0 == prim) failed_num += 2;
          break;
        case 2:
          geom_prim.Torus(2.0, 1.0, prim);
          if (0 == prim) failed_num += 4;
          break;
      }
      if (0 != prim) prims.set(i, prim);
    }
  } catch (TSTTB::Error err) {
    std::cerr << "Creating primitive failed, failed_num = " << failed_num << std::endl;
    return false;
  }

    // verify the bounding boxes
  sidl::array<double> max_corn, min_corn;
  int num_max, num_min;
  try {
    geom_shape.gentityBoundingBox(prims, 3, 
                                  min_corn, num_min, max_corn, num_max);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems getting bounding boxes of primitives." << std::endl;
    return false;
  }

  failed_num = 0;
  double preset_min_corn[] = 
      // min brick corner xyz
    {-0.5, -1.0, -1.5, 
      // min cyl corner xyz
     -4.0, -2.0, -0.5,
      // min torus corner xyz
     -3.0, -3.0, -1.0
    };
  
  double preset_max_corn[] = 
      // max brick corner xyz
    {0.5, 1.0, 1.5, 
      // max cyl corner xyz
     4.0, 2.0, 0.5,
      // max torus corner xyz
     3.0, 3.0, 1.0
    };
  
      
  for (int i = 0; i < 3; i++) {
    if (min_corn.get(3*i) != preset_min_corn[3*i] ||
        min_corn.get(3*i+1) != preset_min_corn[3*i+1] ||
        min_corn.get(3*i+2) != preset_min_corn[3*i+2] ||
        max_corn.get(3*i) != preset_max_corn[3*i] ||
        max_corn.get(3*i+1) != preset_max_corn[3*i+1] ||
        max_corn.get(3*i+2) != preset_max_corn[3*i+2])
      failed_num += (i == 0 ? 1 : (i == 2 ? 2 : 4));
  }
  
  if (failed_num != 0) {
    std::cerr << "Bounding box check failed for primitive; failed_num = " 
             << failed_num << std::endl;
    return false;
  }
  
    // must have worked; delete the entities then return
  try {
    for (int i = 0; i < 3; i++)
      geom_construct.Delete(prims.get(i));
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems deleting primitive after boolean check." << std::endl;
    return false;
  }
  
  return true;
}
  
bool transforms_test(TSTTG::Modify &geom) 
{
    // query ones
  CAST_GEOM_INTERFACE(geom, geom_topo, Topology, false);
  CAST_GEOM_INTERFACE(geom, geom_shape, Shape, false);

    // modify ones
  CAST_GEOM_INTERFACE(geom, geom_construct, Construct, false);
  CAST_GEOM_INTERFACE(geom, geom_prim, Primitives, false);
  CAST_GEOM_INTERFACE(geom, geom_transforms, Transforms, false);

    // construct a brick
  GentityHandle brick = 0;
  try {
    geom_prim.Brick(1.0, 2.0, 3.0, brick);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems creating brick for transforms test." << std::endl;
    return false;
  }

    // move it, then test bounding box
  try {
    geom_transforms.Move(brick, 0.5, 1.0, 1.5);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems moving brick for transforms test." << std::endl;
    return false;
  }

  sidl::array<double> bb_min, bb_max;
  int num_min, num_max;
  sidl::array<GentityHandle> dum = sidl::array<GentityHandle>::create1d(1);
  dum.set(0, brick);
  try {
    geom_shape.gentityBoundingBox(dum, 1, bb_min, num_min, bb_max, num_max);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems getting bounding box after move." << std::endl;
    return false;
  }

  num_min = ARRAY_SIZE(bb_min);
  num_max = ARRAY_SIZE(bb_max);
  if (3 != num_min || 3 != num_max ||
      bb_min.get(0) != 0.0 || bb_min.get(1) != 0.0 || bb_min.get(2) != 0.0 ||
      bb_max.get(0) != 1.0 || bb_max.get(1) != 2.0 || bb_max.get(2) != 3.0) {
    std::cerr << "Wrong bounding box after move." << std::endl;
    return false;
  }
  
    // now rotate it about +x, then test bounding box
  try {
    geom_transforms.Rotate(brick, 90, 1.0, 0.0, 0.0);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems rotating brick for transforms test." << std::endl;
    return false;
  }

  try {
    geom_shape.gentityBoundingBox(dum, 1, bb_min, num_min, bb_max, num_max);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems getting bounding box after rotate." << std::endl;
    return false;
  }

  num_min = ARRAY_SIZE(bb_min);
  num_max = ARRAY_SIZE(bb_max);
  if (3 != num_min || 3 != num_max ||
      bb_min.get(0) != 0.0 || bb_min.get(1) != -3.0 || bb_min.get(2) != 0.0 ||
      bb_max.get(0) != 1.0 || bb_max.get(1) < -1.0e-6 || 
      bb_max.get(1) > 1.0e-6 || bb_max.get(2) != 2.0) {
    std::cerr << "Wrong bounding box after rotate." << std::endl;
    return false;
  }
  
    // now reflect through y plane; should recover original bb
  try {
    geom_transforms.Reflect(brick, 0.0, 1.0, 0.0);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems reflecting brick for transforms test." << std::endl;
    return false;
  }

  try {
    geom_shape.gentityBoundingBox(dum, 1, bb_min, num_min, bb_max, num_max);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems getting bounding box after rotate." << std::endl;
    return false;
  }

  num_min = ARRAY_SIZE(bb_min);
  num_max = ARRAY_SIZE(bb_max);
  if (3 != num_min || 3 != num_max ||
      bb_min.get(0) != 0.0 || bb_min.get(1) < -1.0e-6 
      || bb_min.get(1) > 1.0e-6 ||
      bb_min.get(2) < -1.0e-6 || bb_min.get(2) > 1.0e-6 ||
      bb_max.get(0) != 1.0 || bb_max.get(1) != 3.0 || bb_max.get(2) != 2.0) {
    std::cerr << "Wrong bounding box after reflect." << std::endl;
    return false;
  }

    // must have worked; delete the entities then return
  try {
    geom_construct.Delete(brick);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems deleting brick after transforms check." << std::endl;
    return false;
  }
  
  return true;
}


bool booleans_test(TSTTG::Modify &geom) 
{
    // query ones
  CAST_GEOM_INTERFACE(geom, geom_topo, Topology, false);
  CAST_GEOM_INTERFACE(geom, geom_shape, Shape, false);

    // modify ones
  CAST_GEOM_INTERFACE(geom, geom_construct, Construct, false);
  CAST_GEOM_INTERFACE(geom, geom_prim, Primitives, false);
  CAST_GEOM_INTERFACE(geom, geom_transforms, Transforms, false);
  CAST_GEOM_INTERFACE(geom, geom_booleans, Booleans, false);

    // construct a brick size 1, and a cylinder rad 0.25 height 2
  GentityHandle brick = 0, cyl = 0;
  try {
    geom_prim.Brick(1.0, 0.0, 0.0, brick);
    geom_prim.Cylinder(1.0, 0.25, 0.0, cyl);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems creating brick or cylinder for booleans test." << std::endl;
    return false;
  }

    // subtract the cylinder from the brick
  GentityHandle subtract_result = 0;
  try {
    geom_booleans.Subtract(brick, cyl, subtract_result);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems subtracting for booleans subtract test." << std::endl;
    return false;
  }

    // section the brick
  GentityHandle section_result = 0;
  try {
    geom_booleans.Section(subtract_result, 1.0, 0.0, 0.0, 0.25, true, section_result);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems sectioning for booleans section test." << std::endl;
    return false;
  }

    // unite the section result with a new cylinder
  try {
    geom_prim.Cylinder(1.0, 0.25, 0.0, cyl);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems creating cylinder for unite test." << std::endl;
    return false;
  }
  GentityHandle unite_results;
  sidl::array<GentityHandle> unite_input = sidl::array<GentityHandle>::create1d(2);
  unite_input.set(0, section_result);
  unite_input.set(1, cyl);
  try {
    geom_booleans.Unite(unite_input, 2, unite_results);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems uniting for booleans unite test." << std::endl;
    return false;
  }

    // must have worked; delete the entities then return
  try {
    geom_construct.Delete(unite_results);
  }
  catch (TSTTB::Error err) {
    std::cerr << "Problems deleting for booleans unite test." << std::endl;
    return false;
  }
  
  return true;
}
