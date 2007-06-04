#ifndef ACISMEMFILE
#define ACISMEMFILE

//---------------------------------------------------------------------
// Define the AcisMemFile class for doing ACIS save and restore to a memory
// buffer.  This file is intended as an example of how to create a new
// ACIS FileInterface object.
//
// This class is used to store ACIS-specific geometry in a memory file;
// generic memory file functions are implemented in CGMMemFile; ACIS-specific
// functions for gathering ACIS geometry and putting it in a memory buffer
// is implemented in AcisMemFile and in the ACIS BinaryFile class.
//
// There is also a FileSizeComputer class.  This can be used to determine
// how many bytes a list of entities will use when saved.
//---------------------------------------------------------------------

#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernutil/fileio/binfile.hxx"
#else
#include "binfile.hxx"
#endif
#include "CGMMemFile.hpp"
#include "DLIList.hpp"

class AcisQueryEngine;
class RefEntity;
class ENTITY_LIST;
//---------------------------------------------------------------------

class AcisMemFile : public BinaryFile,
                    public CGMMemFile
{

public:

  AcisMemFile(AcisQueryEngine *aqe) {acisQueryEngine = aqe;};
    //- (empty) constructor

  virtual ~AcisMemFile() {};
    //- (empty) destructor

  virtual FilePosition set_mark();
    //- The methods for positioning the file pointer must also be
    //- implemented for each derived class.  These are not normally
    //- used.

  virtual FilePosition goto_mark(FilePosition);
    //- The methods for positioning the file pointer must also be
    //- implemented for each derived class.  These are not normally
    //- used.

  size_t BytesWritten() { return (size_t)m_currentPosition; }
    //- Get the number of bytes written

protected:

  CubitStatus read_refentity_list(DLIList<RefEntity*> &ref_entity_list);
    //- read a RefEntity list from the buffer into ref_entity_list; relies on engine-
    //- specific implementation

  CubitStatus write_refentity_list(DLIList<RefEntity*> &ref_entity_list);
    //- write a RefEntity list to the buffer from ref_entity_list; relies on engine-
    //- specific implementation
  CubitStatus append_refentity_list(DLIList<RefEntity*> &ref_entity_list,
                                    int &buffer_size);

  CubitStatus get_refentity_list(DLIList<RefEntity*> &ref_entity_list);

  virtual int get_refentity_list_size(DLIList<RefEntity*> ref_entity_list);
  //- get a the size of entities

  virtual size_t read(void* buf, size_t length, logical swap);
    //- ACIS-dependent function to read data from a buffer

  virtual void write(const void* data, size_t len, logical swap);
    //- ACIS-dependent function to write data to a buffer

private:
  size_t write_memory_buffer(ENTITY_LIST &entity_list);
    //- write ACIS entity list to memory buffer

  size_t append_memory_buffer(ENTITY_LIST &entity_list,
                              int &buffer_size);

  //unsigned char* get_memory_buffer(ENTITY_LIST &entity_list);
  
  int get_memory_buffer_size(ENTITY_LIST entity_list);
    //- get a the size of entities in memory buffer
  
  size_t read_memory_buffer(ENTITY_LIST &entity_list);
    // read ACIS entity list from memory buffer
  size_t scatter_memory_buffer(ENTITY_LIST &entity_list,
                                          int *buffer_sizes);
  AcisQueryEngine *acisQueryEngine;
    //- for convenience
  
};

//---------------------------------------------------------------------
// This class can be used to compute how many bytes will be required
// to save a list of entities in binary format.  The following example
// shows how it can be used:
//
// ENTITY_LIST elist;
// ... add entities you want saved to the list
// FileSizeComputer fsc;
// outcome result = api_save_entity_list_file(&fsc, elist);
// if(result.ok()) {
//    unsigned long bytesNeeded = sc.GetSize();
// }

class FileSizeComputer : public BinaryFile {

private:
  size_t m_bytesNeeded;
  int num_writes; // number of calls to write function
  
protected:

    // These are the virtual methods which actually do the reading
    // and writing of the data.

    // Read doesn't really do anything.
  virtual size_t read(void* buf, size_t length, logical swap);
  virtual void write(const void* data, size_t len, logical swap);

public:

    // The methods for positioning the file pointer must also be
    // implemented for each derived class.  These do not do
    // anything in this class
  virtual FilePosition set_mark();
  virtual FilePosition goto_mark(FilePosition);

public:
  FileSizeComputer();
  virtual ~FileSizeComputer();
  size_t GetSize();
  int GetWrites();
};

//---------------------------------------------------------------------
#endif
