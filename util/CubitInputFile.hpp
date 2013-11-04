#ifndef CUBITINPUTFILE_HPP
#define CUBITINPUTFILE_HPP

#include "CubitString.hpp"
#include <cstdio>
#include "CubitMessage.hpp"
#include "CubitFileUtil.hpp"

struct CubitInputFile
{
  enum FileType
  {
    FILE_NORMAL=1,
    FILE_FASTQ=2,
    FILE_TEMPORARY=3
  };

  CubitString              filename;
  CubitFile                filePointer;
  int                      lineNumber;
  CubitInputFile::FileType fileType;
  int                      loopCount;
  int                      breakPoint;
  
  CubitInputFile(const CubitString& fileName,
                 FileType type=FILE_NORMAL,
                 int loop=1,
                 char *default_path=NULL);
  
  
  ~CubitInputFile();
  
};

inline CubitInputFile::CubitInputFile(const CubitString& fileName,
                                      CubitInputFile::FileType type,
                                      int loop,
                                      char *includePath)
  : breakPoint(0)
{
  CubitString file_and_path;
  filePointer.open(fileName, "r");

  if (!filePointer && includePath) {
    file_and_path = includePath;
    file_and_path += "/";
    file_and_path += fileName;
    filePointer.open(file_and_path.c_str(), "r");
  }

  if (filePointer) {
    if (includePath)
      filename  = file_and_path;
    else
      filename  = fileName;
  }
  else {
    filename = "<Invalid File>";
    PRINT_WARNING("Could not open file: %s\n", fileName.c_str() );
  }
  lineNumber  = 1;
  fileType = type;
  loopCount = --loop;
}

inline CubitInputFile::~CubitInputFile() {
  filePointer.close();
  if (fileType == FILE_TEMPORARY)
    CubitFileUtil::remove_file(filename); /* Delete file if temporary */
}

#endif 

