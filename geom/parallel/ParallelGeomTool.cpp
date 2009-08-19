#include <cstdio>

#include "CubitString.hpp"
#include "CubitMessage.hpp"
#include "DLList.hpp"
#include "RefEntity.hpp"
#include "CubitEntity.hpp"
#include "CastTo.hpp"
#include "CubitUtil.hpp"

#include "TopologyBridge.hpp"
#include "GeometryQueryTool.hpp"
#include "ParallelGeomTool.hpp"
#include "CGMParallelConventions.h"
#include "CATag.hpp"
#include "CGMParallelComm.hpp"

const bool debug = true;

enum CGMParallelActions {PA_READ=0, PA_BROADCAST, PA_DELETE_NONLOCAL,
			 PA_SCATTER,
			 PA_CHECK_GIDS_SERIAL, PA_GET_ENTS, 
			 PA_RESOLVE_SHARED_ENTS,
			 PA_EXCHANGE_GHOSTS};

enum CGMPartitionActions {PT_GEOM_DIM=0, PT_PAR_PART};

const char *CGMParallelActionsNames[] = {
  "PARALLEL READ",
  "PARALLEL BROADCAST", 
  "PARALLEL DELETE NONLOCAL",
  "PARALLEL SCATTER",
  "PARALLEL CHECK_GIDS_SERIAL",
  "PARALLEL GET_FILESET_ENTS",
  "PARALLEL RESOLVE_SHARED_ENTS",
  "PARALLEL EXCHANGE_GHOSTS"
};

const char* ParallelGeomTool::CGMparallelOptsNames[] = { "NONE", "READ", "READ_DELETE", "BCAST", 
							 "BCAST_DELETE", "SCATTER", "SCATTER_DELETE",
							 "READ_PARALLEL", "FORMAT", "", 0 };

const char* ParallelGeomTool::CGMpartitionOptsNames[] = { "NONE", "GEOM_DIMENSION",
							  "PARARELL_PARTITION", "", 0 };

//ParallelGeomTool *ParallelGeomTool::instance_ = NULL;

ParallelGeomTool::ParallelGeomTool(CGMTagManager* impl, CGMParallelComm *pc)
  : cgmImpl(impl), myPcomm(pc) 
{
  if (!myPcomm) {
    myPcomm = CGMParallelComm::get_pcomm(impl, 0);
    if (NULL == myPcomm) myPcomm = new CGMParallelComm(cgmImpl);
  }
  //gqt = GeometryQueryTool::instfance();
}

CubitStatus ParallelGeomTool::load_file(const char *file_name,
					const char *options,
					const char* set_tag_name,
					const int* set_tag_values,
					int num_set_tag_values) 
{
  FileOptions opts(options);

  // Get parallel settings
  int parallel_mode;
  FOErrorCode result = opts.match_option("PARALLEL", CGMparallelOptsNames, 
					 parallel_mode);
  if (FO_FAILURE == result) {
    PRINT_ERROR( "Unexpected value for 'PARALLEL' option\n" );
    return CUBIT_FAILURE;
  }
  else if (FO_ENTITY_NOT_FOUND == result) {
    parallel_mode = 0;
  }

  bool distrib = false;
  bool surf_partition = false;
  bool body_partition = false;
  std::string partition_tag_name;
  std::vector<int> partition_tag_vals;

  // Get partition tag value(s), if any, and whether they're to be
  result = opts.get_ints_option("PARTITION_VAL", partition_tag_vals);

  // Get partition setting
  result = opts.get_option("PARTITION", partition_tag_name);

  if (FO_ENTITY_NOT_FOUND == result || partition_tag_name.empty()) {
    partition_tag_name = PARALLEL_PARTITION_TAG_NAME;
    distrib = true;
  }
  else {
    if (partition_tag_name == "GEOM_DIMENSION") {
      int geom_dim = 0;
      for (std::vector<int>::iterator pit = partition_tag_vals.begin(); 
	   pit != partition_tag_vals.end(); pit++) {
	geom_dim = *pit;
	if (geom_dim == 2) surf_partition = true; // body & surface distribution
	else if (geom_dim == 3) body_partition = true; // body only distribution
      }
      if (!surf_partition && !body_partition) {
	PRINT_ERROR("Geometry dimension %d is not supported.\n", geom_dim);
	return CUBIT_FAILURE;
      }
    }
    else if (partition_tag_name == "PARALLEL_PARTITION") {
    }

    // distributed or assigned
    result = opts.get_null_option("PARTITION_DISTRIBUTE");
    if (FO_SUCCESS == result) {
      distrib = true;
      body_partition = true;
    }
  }

  // get MPI IO processor rank
  int reader_rank;
  result = opts.get_int_option("MPI_IO_RANK", reader_rank);
  if (FO_ENTITY_NOT_FOUND == result)
    reader_rank = 0;
  else if (FO_SUCCESS != result) {
    PRINT_ERROR( "Unexpected value for 'MPI_IO_RANK' option\n" );
    return CUBIT_FAILURE;
  }

  // now that we've parsed all the parallel options, make an instruction
  // queue
  std::vector<int> pa_vec;
  bool is_reader = (reader_rank == (int) myPcomm->proc_config().proc_rank());
  
  switch (parallel_mode) {

  case POPT_READ:
    pa_vec.push_back(PA_READ);
    //pa_vec.push_back(PA_CHECK_GIDS_SERIAL);
    pa_vec.push_back(PA_GET_ENTS);
    break;
  
  case POPT_DEFAULT:
  case POPT_READ_DELETE:
    pa_vec.push_back(PA_READ);
    //pa_vec.push_back(PA_CHECK_GIDS_SERIAL);
    pa_vec.push_back(PA_GET_ENTS);
    pa_vec.push_back(PA_DELETE_NONLOCAL);
    break;

  case POPT_BCAST:
    if (is_reader) {
      pa_vec.push_back(PA_READ);
      //pa_vec.push_back(PA_CHECK_GIDS_SERIAL);
      pa_vec.push_back(PA_GET_ENTS);
    }
    pa_vec.push_back(PA_BROADCAST);
    if (!is_reader) pa_vec.push_back(PA_GET_ENTS);
    break;
    
  case POPT_BCAST_DELETE:
    if (is_reader) {
      pa_vec.push_back(PA_READ);
      //pa_vec.push_back(PA_CHECK_GIDS_SERIAL);
      pa_vec.push_back(PA_GET_ENTS);
    }
    pa_vec.push_back(PA_BROADCAST);
    //if (!is_reader) pa_vec.push_back(PA_GET_FILESET_ENTS);
    pa_vec.push_back(PA_DELETE_NONLOCAL);
    break;
    
  case PORT_SCATTER:
    if (is_reader) {
      pa_vec.push_back(PA_READ);
      //pa_vec.push_back(PA_CHECK_GIDS_SERIAL);
      pa_vec.push_back(PA_GET_ENTS);
    }
    pa_vec.push_back(PA_SCATTER);
    if (!is_reader) pa_vec.push_back(PA_GET_ENTS);
    break;
      
  case PORT_SCATTER_DELETE:
    if (is_reader) {
      pa_vec.push_back(PA_READ);
      //pa_vec.push_back(PA_CHECK_GIDS_SERIAL);
      pa_vec.push_back(PA_GET_ENTS);
    }
    pa_vec.push_back(PA_SCATTER);
    if (is_reader) pa_vec.push_back(PA_DELETE_NONLOCAL);
    else pa_vec.push_back(PA_GET_ENTS);
    break;
  
  case POPT_FORMAT:
    PRINT_ERROR( "Access to format-specific parallel read not implemented.\n");
    return CUBIT_FAILURE;

  case POPT_READ_PARALLEL:
    PRINT_ERROR( "Partitioning for PARALLEL=READ_PARALLEL not supported yet.\n");
    return CUBIT_FAILURE;

  default:
    return CUBIT_FAILURE;
  }

  return load_file(file_name, parallel_mode, 
                   partition_tag_name,
                   partition_tag_vals, distrib, pa_vec, opts,
                   set_tag_name, set_tag_values,
		   num_set_tag_values,
                   reader_rank, surf_partition, body_partition
		   );
}

CubitStatus ParallelGeomTool::load_file(const char *file_name,
                        		int parallel_mode, 
                        		std::string &partition_tag_name, 
					std::vector<int> &partition_tag_vals, 
					bool distrib,
					std::vector<int> &pa_vec,
					const FileOptions &opts,
					const char* set_tag_name,
					const int* set_tag_values,
					const int num_set_tag_values,
					const int reader_rank,
					const bool surf_partition,
					const bool body_partition
					)
{
  // check file type
  CubitString file_type;
  //if (strstr(file_name, ".sab")) file_type = "ACIS_SAB";
  //else if (strstr(file_name, ".sat")) file_type = "ACIS_SAT";
  if (strstr(file_name, ".stp")) file_type = "STEP";
  else if (strstr(file_name, ".igs")) file_type = "IGES";
  else if (strstr(file_name, ".occ") || strstr(file_name, ".brep")) file_type = "OCC";
  else {
    PRINT_ERROR("File type not known for file %s; skipping.\n", file_name);
    return CUBIT_FAILURE;
  }
  
  // do the work by options
  bool i_read = false;
  std::vector<int>::iterator vit;
  int i;
  DLIList<RefEntity*> surf_entity_list, body_entity_list;
  for (i = 1, vit = pa_vec.begin(); vit != pa_vec.end(); vit++, i++) {
    CubitStatus tmp_result = CUBIT_SUCCESS;
    switch (*vit) {
//==================
    case PA_READ:
      i_read = true;
      double tStart, tEnd;

      if (debug) {
	std::cout << "Reading file " << file_name << std::endl;
	tStart = MPI_Wtime();
      }
      tmp_result = GeometryQueryTool::instance()->import_solid_model(file_name, file_type.c_str());

      if (CUBIT_SUCCESS != tmp_result) {
	PRINT_ERROR("Reading file %s failed.\n", file_name);
	return CUBIT_FAILURE;
      }
      else if (debug) {
	tEnd = MPI_Wtime();
	PRINT_INFO("Read time in proc %d is %f.\n", myPcomm->proc_config().proc_size(),
		   tEnd - tStart);
	PRINT_INFO("Read done.\n");
      }
      
      break;

//==================
    case PA_GET_ENTS:
      if (debug) std::cout << "Getting entities." << std::endl;
      
      // get entities
      if (body_partition) {
	tmp_result = GeometryQueryTool::instance()->ref_entity_list("body", body_entity_list, CUBIT_FALSE);
	
	if (CUBIT_SUCCESS != tmp_result) {
	  PRINT_ERROR("Getting body entities failed.\n");
	  return CUBIT_FAILURE;
	}
      }
      if (surf_partition) {
	tmp_result = GeometryQueryTool::instance()->ref_entity_list("surface", surf_entity_list, CUBIT_FALSE);
	
	if (CUBIT_SUCCESS != tmp_result) {
	  PRINT_ERROR("Getting surf entities failed.\n");
	  return CUBIT_FAILURE;
	}
      }
      if (debug) PRINT_INFO("Getting entities done.\n");

      break;
      
//==================     
    case PA_DELETE_NONLOCAL:
      if (debug) {
	PRINT_INFO("Deleting nonlocal entities.\n");
	tStart = MPI_Wtime();
      }
     
      tmp_result = delete_nonlocal_entities(partition_tag_name, 
					    partition_tag_vals,
					    surf_entity_list,
					    body_entity_list,
					    distrib);
     
      if (CUBIT_SUCCESS != tmp_result) {
	PRINT_ERROR("Delete failed.\n");
	return CUBIT_FAILURE;
      }
      else if (debug) {
	tEnd = MPI_Wtime();
	PRINT_INFO("Delete done.\n");
	PRINT_INFO("Delete time in proc %d is %f.\n", myPcomm->proc_config().proc_size(),
		   tEnd - tStart);
      }
      break;

//==================      
    case PA_BROADCAST:
      // do the actual broadcast; if single-processor, ignore error
      if (myPcomm->proc_config().proc_size() > 1) {
	if (body_partition) {
	  if (debug) {
	    PRINT_INFO("Broadcasting body entities.\n");
	    tStart = MPI_Wtime();
	  }

	  tmp_result = myPcomm->broadcast_entities(reader_rank, body_entity_list);
	  
	  if (CUBIT_SUCCESS != tmp_result) {
	    PRINT_ERROR("Broadcasting body entities failed.\n");
	    return CUBIT_FAILURE;
	  }
	  else if (debug) {
	    tEnd = MPI_Wtime();
	    PRINT_INFO("Bcast bodies done.\n");
	    PRINT_INFO("Broadcast bodies time in proc %d is %f.\n", myPcomm->proc_config().proc_size(),
		       tEnd - tStart);
	  }
	}
	if (surf_partition) {
	  if (debug) {
	    PRINT_INFO("Broadcasting surface entities.\n");
	    tStart = MPI_Wtime();
	  }
	  tmp_result = myPcomm->broadcast_entities(reader_rank, surf_entity_list);
	  
	  if (CUBIT_SUCCESS != tmp_result) {
	    PRINT_ERROR("Broadcasting surface entities failed.\n");
	    return CUBIT_FAILURE;
	  }
	  else if (debug) {
	    tEnd = MPI_Wtime();
	    PRINT_INFO("Bcast surface done.\n");
	    PRINT_INFO("Broadcast surfaces time in proc %d is %f.\n", myPcomm->proc_config().proc_size(),
		       tEnd - tStart);
	  }
	}
      }
      
      break;

//==================      
    case PA_SCATTER:
      // do the actual scatter
      if (myPcomm->proc_config().proc_size() > 1) {
	if (body_partition) {
	  if (debug) {
	    PRINT_INFO("Scattering body entities.\n");
	    tStart = MPI_Wtime();
	  }
	  tmp_result = myPcomm->scatter_entities(reader_rank, body_entity_list);
	  
	  if (CUBIT_SUCCESS != tmp_result) {
	    PRINT_ERROR("Scattering body entities failed.\n");
	    return CUBIT_FAILURE;
	  }
	  else if (debug) {
	    tEnd = MPI_Wtime();
	    PRINT_INFO("Scatter bodies done.\n");
	    PRINT_INFO("Scatter bodies time in proc %d is %f.\n", myPcomm->proc_config().proc_size(),
		       tEnd - tStart);
	  }
	}
	if (surf_partition) {
	  if (debug) {
	    PRINT_INFO("Scattering surface entities.\n");
	    tStart = MPI_Wtime();
	  }
	  tmp_result = myPcomm->scatter_entities(reader_rank, surf_entity_list);
	  
	  if (CUBIT_SUCCESS != tmp_result) {
	    PRINT_ERROR("Scattering surf entities failed.\n");
	    return CUBIT_FAILURE;
	  }
	  else if (debug) {
	    tEnd = MPI_Wtime();
	    PRINT_INFO("Scatter surfaces done;\n");
	    PRINT_INFO("Scatter surfaces time in proc %d is %f.\n", myPcomm->proc_config().proc_size(),
		       tEnd - tStart);
	  }
	}
      }
      if (debug) PRINT_INFO("Scatter done.\n");
      
      break;

//==================    
    default:
      return CUBIT_FAILURE;
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus ParallelGeomTool::delete_nonlocal_entities(std::string &ptag_name,
						       std::vector<int> &ptag_vals,
						       DLIList<RefEntity*> surf_entity_list,
						       DLIList<RefEntity*> body_entity_list,
						       const bool distribute)
{
  unsigned int i;
  CubitStatus result;
  int proc_sz = myPcomm->proc_config().proc_size();
  int proc_rk = myPcomm->proc_config().proc_rank();
  
  if (distribute) { // Round-Robin
    // for now, require that number of partition sets be greater
    // than number of procs
    if (body_entity_list.size() < proc_sz &&
	surf_entity_list.size() < proc_sz) {
      result = CUBIT_FAILURE;
      PRINT_ERROR("Number of procs greater than number of partitions.");
      return CUBIT_FAILURE;
    }
    
    // delete unneeded surf entities
    unsigned int tot_entity = surf_entity_list.size();
    unsigned int n_entity = tot_entity / proc_sz;
    unsigned int n_entity_leftover = tot_entity % proc_sz;
    unsigned int begin_entity = 0;

    if (proc_rk < (int) n_entity_leftover) {
      n_entity++;
      begin_entity = n_entity * proc_rk;
    }
    else
      begin_entity = proc_rk * n_entity + n_entity_leftover;

    DLIList<RefEntity*> tmp_entity_list, delete_surf_list, delete_body_list;
    unsigned int end_entity = begin_entity + n_entity;
    for (i = 0; i < tot_entity; i++) {
      if (i >= begin_entity && i < end_entity) tmp_entity_list.append(surf_entity_list[i]);
      else {
	delete_surf_list.append(surf_entity_list[i]);
	//GeometryQueryTool::instance()->delete_RefEntity(surf_entity_list[i]);
      }
    }
    
    // change partition surf ref entity list
    myPcomm->partition_surf_list().clean_out();
    myPcomm->partition_surf_list() += tmp_entity_list;

    // delete unneeded body entities
    tot_entity = body_entity_list.size();
    n_entity = tot_entity / proc_sz;
    n_entity_leftover = tot_entity % proc_sz;

    if (proc_rk < (int) n_entity_leftover) {
      n_entity++;
      begin_entity = n_entity * proc_rk;
    }
    else
      begin_entity = proc_rk * n_entity + n_entity_leftover;

    end_entity = begin_entity + n_entity;
    for (i = 0; i < tot_entity; i++) {
      if (i >= begin_entity && i < end_entity) tmp_entity_list.append(body_entity_list[i]);
      else {
	delete_body_list.append(body_entity_list[i]);
	//GeometryQueryTool::instance()->delete_RefEntity(body_entity_list[i]);
      }
    }

    char pre_surf[100], pre_body[100];
    DLIList<CubitEntity*> tmp_surf_list, tmp_body_list;
    if (debug) {
      CAST_LIST_TO_PARENT(delete_surf_list, tmp_surf_list);
      sprintf( pre_surf, "Deleted %d Surfaces: ", tmp_surf_list.size() );
      CAST_LIST_TO_PARENT(delete_body_list, tmp_body_list);
      sprintf( pre_body, "Deleted %d Bodies: ", tmp_body_list.size() );
      CubitUtil::list_entity_ids( pre_surf, tmp_surf_list );
    }

    unsigned int nDelete = delete_surf_list.size();
    for (i = 0; i < nDelete; i++) {
      GeometryQueryTool::instance()->delete_RefEntity(surf_entity_list[i]);
    }
    
    if (debug) CubitUtil::list_entity_ids( pre_body, tmp_body_list );
    
    nDelete = delete_body_list.size();
    for (i = 0; i < nDelete; i++) {
      GeometryQueryTool::instance()->delete_RefEntity(body_entity_list[i]);
    }

    // change partition surf ref entity list
    myPcomm->partition_body_list().clean_out();
    myPcomm->partition_body_list() += tmp_entity_list;
  }
  else {
    myPcomm->partition_surf_list().clean_out();
    myPcomm->partition_surf_list() += surf_entity_list;
    myPcomm->partition_body_list().clean_out();
    myPcomm->partition_body_list() += body_entity_list;
  }
       
  if (debug) {
    /*
    surf_entity_list.clean_out();
    GeometryQueryTool::instance()->ref_entity_list("surface", surf_entity_list, CUBIT_FALSE);
    body_entity_list.clean_out();
    GeometryQueryTool::instance()->ref_entity_list("body", body_entity_list, CUBIT_FALSE);
    std::cerr << "My partition surf/body ref entity list size1: "
    << surf_entity_list.size() << " "
	      << body_entity_list.size() << std::endl;*/

    std::cerr << "My partition surf/body ref entity list size: "
	      << myPcomm->partition_surf_list().size() << " "
	      << myPcomm->partition_body_list().size() << std::endl;
  }    
  
  return CUBIT_SUCCESS;
}
