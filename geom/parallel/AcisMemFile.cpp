//---------------------------------------------------------------------
// Implement the AcisMemFile class for doing ACIS save and restore to a
// memory buffer.
//---------------------------------------------------------------------

#include <memory.h>
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kerndata/lists/lists.hxx"
#include "kernel/kernutil/fileio/binfile.hxx"
#include "kernel/kerndata/savres/savres.hxx"

#include "kernel/kerndata/savres/fileinfo.hxx"
#include "kernel/kernapi/api/kernapi.hxx"

#include "kernel/kernapi/api/api.hxx"
#include "kernel/kerndata/top/alltop.hxx"
#include "kernel/kerndata/data/entity.hxx"
#include "cstr/constrct/kernapi/api/cstrapi.hxx"
#else
#include "lists.hxx"
#include "binfile.hxx"
#include "savres.hxx"

#include "fileinfo.hxx"
#include "kernapi.hxx"

#include "api.hxx"
#include "alltop.hxx"
#include "entity.hxx"
#include "cstrapi.hxx"
#endif
// ********** END ACIS INCLUDES               **********

#include "AcisMemFile.hpp"
#include "GeometryQueryTool.hpp"
#include "AcisQueryEngine.hpp"
#include "TopologyEntity.hpp"
#include "AcisBridge.hpp"
#include "TopologyBridge.hpp"
#include "RefEntity.hpp"
#include "DLIList.hpp"
#include "CastTo.hpp"
#include "CubitMessage.hpp"
#include "ProcData.hpp"
#include "CubitDefines.h"

extern DECL_KERN logical save_entity_list_on_file(
            FileInterface*,     // open file pointer
            ENTITY_LIST const&  // list of ENTITIES to save
        );
extern DECL_KERN logical restore_entity_list_from_file(
            FileInterface*,     // open file pointer
            ENTITY_LIST&        // list of restored ENTITIES
        );

//---------------------------------------------------------------------
// Purpose--
//	Create the AcisMemFile object.
//---------------------------------------------------------------------

//-------------------------------------
// Write num bytes of data with no translation
void AcisMemFile::write(const void* data, size_t num, logical)
{
    // See how much space is left in the buffer
  size_t spaceLeft = m_sizeBuffer - m_currentPosition;

  if( spaceLeft < num ) {
      // There is not enough space. Allocate a bigger buffer
    unsigned char *temp_buff_ptr = (unsigned char *) malloc(num+m_sizeBuffer);
    memcpy(temp_buff_ptr, m_pBuffer, m_currentPosition);
    delete m_pBuffer;
    m_pBuffer = temp_buff_ptr;
  }

    // Copy the data into the buffer
  memcpy(m_pBuffer+m_currentPosition, data, num);
  m_currentPosition += num;

}

//-------------------------------------
// Read from the buffer
size_t AcisMemFile::read(void* buf, size_t bytesNeeded, logical)
{
    // See how much data is left in the buffer
  size_t bytesLeft = m_sizeBuffer - m_currentPosition;

	// See if we can get everything we need
  size_t bytesRead = bytesLeft >= bytesNeeded ? bytesNeeded : bytesLeft;

    // Get the data
  if( bytesRead > 0 ) {
    memcpy(buf, m_pBuffer+m_currentPosition, bytesRead);
    m_currentPosition += bytesRead;
  }

  return bytesRead;
}

//-------------------------------------
FilePosition AcisMemFile::set_mark()
{
  return (FilePosition)m_currentPosition;
}

//-------------------------------------
FilePosition AcisMemFile::goto_mark(FilePosition pos)
{
  if( pos > 0 ) {
    if( (unsigned long)pos < m_sizeBuffer ) {
      m_currentPosition = (unsigned long)pos;
    } else {
      m_currentPosition = m_sizeBuffer;
    }
  } else {
    m_currentPosition = 0;
  }

return (FilePosition)m_currentPosition;
}

size_t AcisMemFile::write_memory_buffer(ENTITY_LIST &entity_list) 
{

    // first, create and write a FileSizeComputer
  FileSizeComputer file_size_computer;

  //  BinaryFile bin_computer;
  //  cout << "AcisMemFile::write_memory_buffer1" << endl;

  FileInfo file_info;
  file_info.reset();
  file_info.set_units(1.0);
  file_info.set_product_id("test");

  api_set_file_info(FileUnits, file_info);
  api_set_file_info(FileId, file_info);

  api_save_entity_list_file(&file_size_computer, entity_list);

  //cout << "AcisMemFile::write_memory_buffer2" << endl;
  // now check the size, allocate if necessary, and write to the real buffer
  if (file_size_computer.GetSize() > m_sizeBuffer) {
    //cout << "AcisMemFile::write_memory_buffer3" << endl;
    delete m_pBuffer;
    m_pBuffer = (unsigned char *) malloc(file_size_computer.GetSize());
    if (m_pBuffer == NULL) {
      PRINT_ERROR("Couldn't allocate big enough buffer to pass model.\n");
      return 0;
    }

    PRINT_DEBUG_100("Allocated buffer of size %d bytes.\n", file_size_computer.GetSize());
    
    m_sizeBuffer = file_size_computer.GetSize();
  }
  //cout << "AcisMemFile::write_memory_buffer4" << endl;

  PRINT_DEBUG_100("Saving entity list into buffer.\n");
  //cout << "AcisMemFile::write_memory_buffer5" << endl;
  api_save_entity_list_file(this, entity_list);
  //cout << "AcisMemFile::write_memory_buffer6" << endl;

  return 1;
}

size_t AcisMemFile::append_memory_buffer(ENTITY_LIST &entity_list,
                                         int &buffer_size) 
{
    // if buffer_size is zero, need to check the size
  if (buffer_size == 0) {
    buffer_size = get_memory_buffer_size(entity_list) + m_currentPosition;
    check_size(buffer_size, CUBIT_TRUE);
  }

  PRINT_DEBUG_100("Saving entity list into buffer.\n");
  api_save_entity_list_file(this, entity_list);
  return 1;
}

size_t AcisMemFile::scatter_memory_buffer(ENTITY_LIST &entity_list,
                                          int *buffer_sizes) 
{
  // first, create and write a FileSizeComputer
  ENTITY_LIST dummy_list;
  unsigned int total_size = 0;
  for (int i = 0; i < entity_list.count(); i++) {
    FileSizeComputer file_size_computer;
    FileInfo file_info;
    file_info.reset();
    file_info.set_units(1.0);
    file_info.set_product_id("test");

    api_set_file_info(FileUnits, file_info);
    api_set_file_info(FileId, file_info);

    dummy_list.clear();
    dummy_list.add(entity_list[i]);
    dummy_list.add(entity_list[i]);
    api_save_entity_list_file(&file_size_computer, dummy_list);
    buffer_sizes[i] = file_size_computer.GetSize();
    total_size += buffer_sizes[i];
  }

  // now check the size, allocate if necessary, and write to the real buffer
  if (total_size > m_sizeBuffer) {
    delete m_pBuffer;
    m_pBuffer = (unsigned char *) malloc(total_size);
    if (m_pBuffer == NULL) {
      PRINT_ERROR("Couldn't allocate big enough buffer to pass model.\n");
      return 0;
    }

    PRINT_DEBUG_100("Allocated buffer of size %d bytes.\n", total_size);
    
    m_sizeBuffer = total_size;
  }

  //cout << "total_size=" << total_size << endl;

  PRINT_DEBUG_100("Saving entity list into buffer.\n");
  for (int i = 0; i < entity_list.count(); i++) {
    dummy_list.clear();
    //cout << "m_sizeBuffer before save" << i << '=' << m_sizeBuffer << endl; 
    //cout << "m_currentPosition before save" << i << '=' << m_currentPosition << endl;
    dummy_list.add(entity_list[i]);
    dummy_list.add(entity_list[i]);
    api_save_entity_list_file(this, dummy_list);
    //cout << "m_sizeBuffer after save" << i << '=' << m_sizeBuffer << endl; 
    //cout << "m_currentPosition after save" << i << '=' << m_currentPosition << endl;
  }

  return 1;
}

/*
unsigned char* AcisMemFile::get_memory_buffer(ENTITY_LIST &entity_list) 
{
  // first, create and write a FileSizeComputer
  FileSizeComputer file_size_computer;
  //  AcisMemFile acmFile;
  //save_entity_list_on_file(&acmFile, entity_list);
  return acmFile.m_pBuffer;
}
*/
int AcisMemFile::get_memory_buffer_size(ENTITY_LIST entity_list) 
{
  // first, create and write a FileSizeComputer
  FileSizeComputer file_size_computer;
  FileInfo file_info;
  file_info.reset();
  file_info.set_units(1.0);
  file_info.set_product_id("test");

  api_set_file_info(FileUnits, file_info);
  api_set_file_info(FileId, file_info);

  api_save_entity_list_file(&file_size_computer, entity_list);
  return file_size_computer.GetSize();
}

size_t AcisMemFile::read_memory_buffer(ENTITY_LIST &entity_list)
{
  if (m_pBuffer == NULL) return 1;
  
  PRINT_DEBUG_100("Reading entity list from buffer on proc %d.\n", ProcData::instance()->myRank);
  restore_entity_list_from_file(this, entity_list);

  return 1;
}

//---------------------------------------------------------------------

// Implementation of the FileSizeComputer class

//---------------------------------------------------------------------
FileSizeComputer::FileSizeComputer() 
    : BinaryFile()
{
  m_bytesNeeded = 0;
  num_writes = 0;
}

//-------------------------------------
FileSizeComputer::~FileSizeComputer()
{
}

//-------------------------------------
// Write num bytes of data with no translation
void FileSizeComputer::write(const void* , size_t num, logical)
{
  m_bytesNeeded += num;
}

//-------------------------------------
// Read from the buffer
size_t FileSizeComputer::read(void*, size_t, logical)
{
  return 0;
}

//-------------------------------------
FilePosition FileSizeComputer::set_mark()
{
  return (FilePosition)0;
}

//-------------------------------------
FilePosition FileSizeComputer::goto_mark(FilePosition pos)
{
  return pos;
}

//-------------------------------------
size_t
FileSizeComputer::GetSize()
{
  return m_bytesNeeded;
}

//---------------------------------------------------------------------

//-------------------------------------
int FileSizeComputer::GetWrites()
{
  return num_writes;
}

//---------------------------------------------------------------------

CubitStatus AcisMemFile::read_refentity_list(DLIList<RefEntity*> &ref_entity_list) 
{
  ENTITY_LIST entity_list;
  goto_mark(0);
  read_memory_buffer(entity_list);

  DLIList<TopologyBridge*> bridge_list;
  acisQueryEngine->restore_entity_list(entity_list, 
                                       bridge_list,
				       CUBIT_FALSE);
  /*
                                       CUBIT_FALSE, CUBIT_FALSE, CUBIT_TRUE,
				       CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE,
				       CUBIT_TRUE, CUBIT_TRUE, CUBIT_FALSE);
  */
  PRINT_INFO("bridge_list_size=%d\n", bridge_list.size());
  GeometryQueryTool::instance()->construct_refentities(bridge_list, &ref_entity_list);
  return CUBIT_SUCCESS;
}

CubitStatus AcisMemFile::write_refentity_list(DLIList<RefEntity*> &ref_entity_list)
{
    // get all the acis entities, then call the internal function
    // 

  int i;
  ENTITY_LIST entity_list;
  AcisBridge *acis_bridge;
  RefEntity *ref_entity;
  /*
  BODY* cuboid;
  api_make_cuboid(100,150,200,cuboid);
  ENTITY_LIST savelist;
  savelist.add(cuboid);
  cout << "AcisMemFile::write_refentity_list1" << endl;

  write_memory_buffer(savelist);
  */
  //cout << "AcisMemFile::write_refentity_list2" << endl;

  for (i = ref_entity_list.size(); i > 0; i--) {
    //cout << "AcisMemFile::write_refentity_list3" << endl;
    ref_entity = ref_entity_list.get_and_step();
    acis_bridge = CAST_TO(CAST_TO(ref_entity, TopologyEntity)->
                          bridge_manager()->topology_bridge(), AcisBridge);
    entity_list.add(acis_bridge->ENTITY_ptr());
  }

  //cout << "AcisMemFile::write_refentity_list4" << endl;
  return (CubitStatus) write_memory_buffer(entity_list);
}

CubitStatus AcisMemFile::append_refentity_list(DLIList<RefEntity*> &ref_entity_list,
                                               int &buffer_size)
{
    // get all the acis entities, then call the internal function
    // 

  int i;
  ENTITY_LIST entity_list;
  AcisBridge *acis_bridge;
  RefEntity *ref_entity;

  for (i = ref_entity_list.size(); i > 0; i--) {
    ref_entity = ref_entity_list.get_and_step();
    acis_bridge = CAST_TO(CAST_TO(ref_entity, TopologyEntity)->
                          bridge_manager()->topology_bridge(), AcisBridge);
    entity_list.add(acis_bridge->ENTITY_ptr());
  }

  return (CubitStatus) append_memory_buffer(entity_list, buffer_size);
}

CubitStatus AcisMemFile::get_refentity_list(DLIList<RefEntity*> &ref_entity_list)
{
    // get all the acis entities, then call the internal function
    // 

  int i;
  ENTITY_LIST entity_list;
  AcisBridge *acis_bridge;
  RefEntity *ref_entity;
  int nEntity = ref_entity_list.size();

  for (i = nEntity; i > 0; i--) {
    ref_entity = ref_entity_list.get_and_step();
    acis_bridge = CAST_TO(CAST_TO(ref_entity, TopologyEntity)->
                          bridge_manager()->topology_bridge(), AcisBridge);
    entity_list.add(acis_bridge->ENTITY_ptr());
  }

  int *bufferSizes = new int[nEntity];
  
  return (CubitStatus) scatter_memory_buffer(entity_list, bufferSizes);
}

int AcisMemFile::get_refentity_list_size(DLIList<RefEntity*> ref_entity_list)
{
  //cout << "AcisMemFile::get_refentity_list_size1"
  //   << ",ref_entity_list.size()=" << ref_entity_list.size() << endl;
  // get the size of entities in memory buffer
  int i;
  ENTITY_LIST entity_list;
  AcisBridge *acis_bridge;
  RefEntity *ref_entity;

  for (i = ref_entity_list.size(); i > 0; i--) {
    //cout << "AcisMemFile::get_refentity_list_size2" << endl;
    ref_entity = ref_entity_list.get_and_step();
    acis_bridge = CAST_TO(CAST_TO(ref_entity, TopologyEntity)->
                          bridge_manager()->topology_bridge(), AcisBridge);
    entity_list.add(acis_bridge->ENTITY_ptr());
  }

  //cout << "AcisMemFile::get_refentity_list_size3" << endl;
  return get_memory_buffer_size(entity_list);
}

/*
CubitStatus AcisMemFile::delete_refentity(DLIList<RefEntity*> &ref_entity_list)
{
  */
