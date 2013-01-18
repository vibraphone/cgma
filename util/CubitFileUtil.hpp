//- Class:          CubitFileUtil
//- Description:    Class with functions to get files from a directory, etc.
//- Owner:          Steve Storm

#ifndef CUBITFILEUTIL_HPP
#define CUBITFILEUTIL_HPP

#include "DLIList.hpp"
#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitUtilConfigure.h"

// an iterator to iterate over files in a directory
class CUBIT_UTIL_EXPORT CubitDirIterator
{
public:
  CubitDirIterator(const CubitString& path,
                   const CubitString& pattern_match = "");
  virtual ~CubitDirIterator();
  void open(const CubitString& path,
            const CubitString& pattern_match = "");
  bool has_next();
  CubitString next();
protected:
  struct Helper;
  Helper* mHelper;
  CubitString mPattern;
  bool mDirs;
  bool atEnd;
  void cleanup();
};

// wrapper for C FILE*
// Includes support for non-ascii filenames.
// Automatically opens/closes the file.
class CUBIT_UTIL_EXPORT CubitFile
{
public:
  // default constructor
  CubitFile();

  // constructor that also opens the file
  CubitFile(const CubitString& file, const char* mode);

  // closes the file, if open
  virtual ~CubitFile();

  // opens a file
  bool open(const CubitString& file, const char* mode);

  // closes the file
  void close();

  // Return the FILE* handle.
  // Returns NULL if status is not 'OK'.
  FILE* file() const;

  // returns whether FILE* is valid
  operator bool () const;

  // Return status of opening the file.
  // For example ENOENT for non-existant file.
  int error();

protected:
  FILE* mFile;
  int mError;
};


class CUBIT_UTIL_EXPORT CubitFileUtil
{
public:

  CubitFileUtil();
  //- Constructor
  
   ~CubitFileUtil();
  //- Destructor

   // returns platform specific path separator
   static const char* separator();

   // returns whether a given path is a directory
   static bool is_directory(const CubitString& path);

   // return whether a given path (file or directory) exists
   static bool path_exists(const CubitString& path);
   
   // returns whether a given path is an executable
   static bool is_executable(const CubitString& path);
  
   // returns whether a path is absolute 
   static bool is_absolute(const CubitString& path);
   
   // returns information about a file (returns 0 on success and errno on failure)
   static int file_info(const CubitString& path, off_t& size, time_t& modification_time, int& mode);


   //! Find the users home directory independent of OS
   static CubitString find_home_path(const CubitString& which_user="");

   //! Get/Set current working directory
   static CubitStatus get_current_working_directory( CubitString& wd );
   static CubitStatus set_current_working_directory( const CubitString& wd );

   //! Given a path append the name onto the path correctly 
   //! \param path, also contains the returned value
   //! \param name
   static CubitString add_name_to_path( const CubitString& path, const CubitString& name);
  
   static CubitStatus create_directory( const CubitString& wd );
   // Create a directory
  
   //! remove a file or directory
   static CubitStatus remove_file( const CubitString& file );

   //! rename a file
   static CubitStatus rename_file( const CubitString& old_file, const CubitString& new_file );

   static CubitStatus get_full_path_str( const CubitString& part,
                                         CubitString &full_path );
   // Given a relative path, return the full path to the filename or 
   // directory.

   static CubitString make_path_platform_compatible( const CubitString& path );
   // Replaces slashes with backslashes or vice versa in the given path.

   static CubitString get_nice_filename( const CubitString& path );
   // Gets a nice looking filename.  All it does is return just the filename
   // without the path if the input filename is in the current working
   // directory (ie., so the filename doesn't contain the full path).

   static void split_path( const CubitString& path, CubitString& dirpart, CubitString& filepart );
   // Split a path string into directory and filename parts.  The function 
   // returns the directory part with a slash (or backslash for NT) on the end.
   // Note the logic is based entirely on looking for the last directory separator 
   // (slash or backslash) in the path - only the string is used (not lstat, for 
   // instance).  Because of this if you know you have a directory coming in put 
   // a slash on the end of your path, otherwise the last directory in your path 
   // will be parsed out as a file.

   static CubitString get_default_cubit_file_name();
   // Gets the default Cubit filename, when the user doesn't specify one (ie., 
   // at Cubit startup.  Format cubitXX.cub, where XX = integer.  Returned 
   // filename does not include full path (assumed current working directory).

   static int get_next_filenumber( const CubitString& filepattern,
                                   int &num_matches,
                                   int &first_match );
   // Gets the next next available number for a filename pattern (use absolute path).
   // The '*' character can be used to denote where a number is searched for.
   // (i.e., use a pattern of file.cub.* to return file.cub.2, if file.cub.1 already exists).
   // Also returns the number of files that matched the file.cub.* pattern in the directory,
   // and the first number found.  Handles cases like file.cub.3.4, etc..
   // (user started from a backup file).
   // Returns -1 if failure.
  
  static CubitString get_file_extension(const CubitString& file, bool remove_version = false /* remove .1, .2, ...*/ );

private:
  static CubitBoolean all_chars_are( char ch, const char *str );
  // Check if all chars from src are ch

  static CubitBoolean is_int_number( const char *str );
  // Check if the string represents an integer number
  
  static CubitBoolean contains_char( char ch, const char *str );
  // Check if the string contains the given character
};

#endif

