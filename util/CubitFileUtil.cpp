//- Class:          CubitFileUtil
//- Description:    Class with functions to get files from a directory, etc.
//- Owner:          Steve Storm

#define NOMINMAX

#include "CubitFileUtil.hpp"
#include "CubitUtil.hpp"
#include "CubitString.hpp"

#ifdef WIN32
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <windows.h>
  #ifndef PATH_MAX
    #define PATH_MAX _MAX_PATH
  #endif
  #include "shlwapi.h"
#else
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <dirent.h>
  #include <cstdlib>
  #include <sys/param.h>
  #include <unistd.h>
  #include <pwd.h>
  #include <fnmatch.h>
#endif
#include <errno.h>


#ifdef WIN32
static const char* DIR_SEP_STR = "\\";
static const char DIR_SEP_CHAR = '\\';
#else
static const char* DIR_SEP_STR = "/";
static const char DIR_SEP_CHAR = '/';
#endif
   
const char* CubitFileUtil::separator()
{
  return DIR_SEP_STR;
}

CubitStatus
CubitFileUtil::get_current_working_directory( CubitString& wd )
{
#ifdef WIN32
  wchar_t* buffer = _wgetcwd( NULL, 0 );
#else
  char* buffer = getcwd( NULL, 0 );
#endif
  if (!buffer)
  {
    PRINT_WARNING( "Unable to get new working directory\n" );
    return CUBIT_FAILURE;
  }
  else
  {
    // convert to string
#ifdef WIN32
    wd = CubitString::toUtf8(buffer);
#else
    wd = buffer;
#endif

    // Add a slash at the end, if not already there
    int wd_len = wd.length();
    if( wd.c_str()[wd_len-1] != DIR_SEP_CHAR )
    {
      wd += DIR_SEP_STR;

      free(buffer);
    }

// TODO: need to figure out the right way to do this. This assumes
// variables of length PATH_MAX which is bad!
//#ifdef WIN32
//    // Make sure format is compatible with full path format
//    CubitString full_path_str;
//    if( get_full_path_str( wd, full_path_str ) == CUBIT_SUCCESS )
//      strcpy( wd, full_path_str.c_str() );
//#endif

    return CUBIT_SUCCESS;
  }
}

CubitStatus CubitFileUtil::set_current_working_directory( const CubitString& wd )
{
#ifdef WIN32
  int ret = _wchdir(CubitString::toUtf16(wd).c_str());
#else
  int ret = chdir(wd.c_str());
#endif
  return ret == 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

CubitString CubitFileUtil::add_name_to_path( const CubitString& path, const CubitString& name )
{
  CubitString result = path;
  // Add a slash at the end of the path, if not already there
  int path_len = result.length();
  if( result.c_str()[path_len-1] != DIR_SEP_CHAR )
  {
    result += DIR_SEP_STR;
  }
  // append the name to the end of the path
  result += name;
      
  return result;
}

CubitString CubitFileUtil::find_home_path(const CubitString& which_user)
{
  CubitString home_dir;

#ifdef WIN32
  home_dir = CubitUtil::getenv("USERPROFILE");
#else
  if(which_user.length() == 0)
  {
    home_dir = CubitUtil::getenv("HOME");
    if( home_dir.length() == 0 )
    {
      struct passwd* userdata = getpwuid( getuid() );
      if( userdata )
        home_dir = userdata->pw_dir;
    }
  }
  else
  {
    struct passwd* userdata = getpwnam( which_user.c_str() );
    if(userdata)
      home_dir = userdata->pw_dir;
  }
#endif

  return home_dir;
}

CubitStatus
CubitFileUtil::create_directory( const CubitString& wd )
{
  // Create the directory
#ifdef WIN32
  if (_wmkdir(CubitString::toUtf16(wd).c_str()) == -1)
  {  
    PRINT_WARNING( "Unable to create new directory\n" );
    return CUBIT_FAILURE;
  }
#else    
  if (mkdir(wd.c_str(), 0777) == -1)
  {  
    PRINT_WARNING( "Unable to create new directory\n" );
    return CUBIT_FAILURE;
  }
#endif    
  return CUBIT_SUCCESS;
}


CubitStatus CubitFileUtil::remove_file( const CubitString& file )
{
#ifdef WIN32
  int status = _wremove(CubitString::toUtf16(file).c_str());
#else
  int status = remove(file.c_str());
#endif
  return status == 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

CubitStatus CubitFileUtil::rename_file( const CubitString& old_file, const CubitString& new_file )
{
#ifdef WIN32
  int status = _wrename(CubitString::toUtf16(old_file).c_str(), CubitString::toUtf16(new_file).c_str());
#else
  int status = rename(old_file.c_str(), new_file.c_str());
#endif
  return status == 0 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

CubitStatus
CubitFileUtil::get_full_path_str( const CubitString& part,
                                  CubitString &full_path_str )
{
  CubitString my_part = CubitFileUtil::make_path_platform_compatible(part);

#ifdef WIN32

  wchar_t* full = _wfullpath(NULL, CubitString::toUtf16(my_part).c_str(), 0);
  if(!full)
  {
    PRINT_ERROR( "problem getting full path to %s\n", part.c_str() );
    return CUBIT_FAILURE;
  }
  full_path_str = CubitString::toUtf8(full);
  free(full);
  
#else

  // we loop removing parts until realpath can resolve an existing path,
  // then add the non-existing parts back on.

  std::vector<CubitString> split_parts;
  CubitString trypart = part;

  if(!CubitFileUtil::is_absolute(trypart))
  {
    CubitString cwd;
    CubitFileUtil::get_current_working_directory(cwd);
    trypart = CubitFileUtil::add_name_to_path(cwd, trypart);
  }

  char full[PATH_MAX];
  while(trypart.length() && !realpath(trypart.c_str(), full))
  {
    CubitString split_part1, split_part2;
    CubitFileUtil::split_path(trypart, split_part1, split_part2);
    split_parts.push_back(split_part2);
    if(split_part1.length() == 0)
    {
      PRINT_ERROR( "problem getting full path to %s\n", part.c_str() );
      return CUBIT_FAILURE;
    }
    trypart = split_part1;
  }

  full_path_str = full;
  for(size_t i=0; i<split_parts.size(); i++)
  {
    full_path_str += CubitString("/") + split_parts[split_parts.size() - i - 1];
  }

#endif

  return CUBIT_SUCCESS;
}

CubitString
CubitFileUtil::make_path_platform_compatible( const CubitString& path)
{
  CubitString ret = path;
  for(size_t i=0; i<ret.length(); i++)
  {
#ifdef WIN32
  // Replace '/' with '\\'
    if(ret.get_at(i) == '/' )
      ret.put_at(i, '\\');
#else
     // Replace '\\' with '/'
    if(ret.get_at(i) == '\\' )
      ret.put_at(i, '/');
#endif
  }
  return ret;
}

CubitString
CubitFileUtil::get_nice_filename( const CubitString& path )
{
  CubitString ret_str;

  CubitString dpart_in, fpart_in;
  split_path( path, dpart_in, fpart_in );

  CubitString wd;
  get_current_working_directory( wd );

  if( dpart_in == wd )
    ret_str = fpart_in;
  else
    ret_str = path;

  return ret_str;
}

void
CubitFileUtil::split_path( const CubitString& path, CubitString& dirpart, CubitString& filepart )
{
  CubitString mypath = path;
  while(mypath.length() && mypath.get_at(mypath.length()-1) == DIR_SEP_CHAR)
  {
    mypath = mypath.substr(0, mypath.length()-1);
  }
  size_t pos = mypath.find_last(DIR_SEP_CHAR);

  // No separator - could be filename or directory.  We assume
  // it's a directory.
  if(pos == MAX_POS)
  {
    filepart = ".";
    dirpart = mypath;
  }
  else
  {
    filepart = mypath.substr(pos+1);
    dirpart = mypath.substr(0, pos);
  }

  // Add slash on end of dirpart if not already there
  if(dirpart.length() && dirpart.get_at(dirpart.length()-1) != DIR_SEP_CHAR)
    dirpart += DIR_SEP_STR;

  return;
}


CubitString
CubitFileUtil::get_file_extension(
  const CubitString& file,
  bool remove_version /* remove .1, .2, ...*/ )
{
  size_t dot_pos = file.find_last('.');
  size_t dot_pos2 = 0;

  if ( dot_pos == MAX_POS )
    return "";

  if(remove_version)
  {
	  dot_pos2 = file.find_last('.',dot_pos);
	  if ( dot_pos2 == MAX_POS )
		  remove_version = false;
	  else if(!is_int_number( file.substr(dot_pos+1).c_str() ))
		  remove_version = false;
  }

  CubitString extension;
  if(!remove_version)
	  extension = file.substr(dot_pos);
  else
	  extension = file.substr(dot_pos2,dot_pos-dot_pos2);
  
  for ( size_t i = 0; i < extension.length(); i++ )
    extension.put_at(i, tolower( extension.get_at(i) ) );
    
  return extension;
  
}  //  get_file_extension()
                     
CubitBoolean 
CubitFileUtil::all_chars_are( char ch, const char *str )
{
  while( *str )
    if( *str++ != ch )
      return CUBIT_FALSE;
    return CUBIT_TRUE;
}

CubitBoolean 
CubitFileUtil::is_int_number( const char *str )
{
  while( *str )
    if( !isdigit(*str++) )
      return CUBIT_FALSE;
    return CUBIT_TRUE;
}

CubitBoolean 
CubitFileUtil::contains_char( char ch, const char *str )
{
  while( *str )
    if( *str++ == ch )
      return CUBIT_TRUE;
    return CUBIT_FALSE;
}

CubitString
CubitFileUtil::get_default_cubit_file_name()
{
  // Get a list of all files matching 'cubit*.cub', in current directory

  CubitDirIterator dir_iter(".", "cubit*.cub");

  int max_number = 0;

  while(dir_iter.has_next())
  {
    CubitString f = dir_iter.next();

    // cut off the "cubit" at the start and the ".cub" at the end
    f = f.substr(5, f.length()-9);
    if( is_int_number(f.c_str()) )
    {
      int num = atoi(f.c_str());
      max_number = std::max(num, max_number);
    }
  }

  max_number++;

  CubitString num_str = CubitString(max_number);
  if(max_number < 10)
    num_str = CubitString("0") + num_str;

  return CubitString("cubit") + num_str + CubitString(".cub");
}

int
CubitFileUtil::get_next_filenumber( const CubitString& file_pattern,
                                    int &num_matches,
                                    int &first_id )
{
  // TODO: check this bug?
  // Continue if nothing after dot.  Also, I noticed the WIN32 will match 
  // file.cub against wildcard file.cub.*

  int max_number = 0;
  int min_number = CUBIT_INT_MAX;
  num_matches = 0;

  // get directory and file parts
  CubitString dir, file, full_file_pattern;
  CubitFileUtil::get_full_path_str(file_pattern, full_file_pattern);
  CubitFileUtil::split_path(full_file_pattern, dir, file);

  size_t wildcard_location = file.find("*");
  if(wildcard_location == MAX_POS)
    return 1;

  CubitDirIterator dir_iter(dir, file);

  while(dir_iter.has_next())
  {
    CubitString f = dir_iter.next();

    // cut off the matched part
    f = f.substr(wildcard_location, 1 + f.length() - file.length());
    if( is_int_number(f.c_str()) )
    {
      int num = atoi(f.c_str());
      max_number = std::max(num, max_number);
      min_number = std::min(num, min_number);
    }
    num_matches++;
  }

  if( min_number != CUBIT_INT_MAX )
    first_id = min_number;

  return max_number+1;
}

struct CubitDirIterator::Helper
{
  Helper()
  {
#ifdef WIN32
    mFileHandle = INVALID_HANDLE_VALUE;
#else
    mDirHandle = NULL;
    mFileHandle = NULL;
#endif
  }

#ifdef WIN32
  WIN32_FIND_DATAW mDirHandle;
  HANDLE mFileHandle;
#else
  DIR* mDirHandle;
  dirent* mFileHandle;
#endif
};
  

CubitDirIterator::CubitDirIterator(const CubitString& path, const CubitString& pattern)
{
  this->mHelper = new Helper;
  open(path, pattern);
}

CubitDirIterator::~CubitDirIterator()
{
  cleanup();
  delete this->mHelper;
}

void CubitDirIterator::open(const CubitString& path, const CubitString& pattern)
{
  cleanup();

  if(pattern.length())
    mPattern = pattern;
  else
    mPattern = "";

  this->atEnd = false;
#ifdef WIN32
  CubitString p = path;
  if (p.get_at(p.length()-1) != DIR_SEP_CHAR)
    p += "\\";
  p += mPattern.length() == 0 ? "*" : mPattern;
  this->mHelper->mFileHandle = FindFirstFileW(CubitString::toUtf16(p).c_str(), &this->mHelper->mDirHandle);
  if(this->mHelper->mFileHandle == INVALID_HANDLE_VALUE)
    this->atEnd = true;
#else
  this->mHelper->mDirHandle = opendir(path.c_str());
  if(this->mHelper->mDirHandle)
  {
    do
    {
      this->mHelper->mFileHandle = readdir(this->mHelper->mDirHandle);
    } while (this->mHelper->mFileHandle && 
             (mPattern.length() != 0 && 
             fnmatch(mPattern.c_str(), this->mHelper->mFileHandle->d_name, 0))
             );
  }
  if(!this->mHelper->mDirHandle || !this->mHelper->mFileHandle)
    this->atEnd = true;
#endif
}

void CubitDirIterator::cleanup()
{
#ifdef WIN32
  if(this->mHelper->mFileHandle != 0)
    FindClose(this->mHelper->mFileHandle);
  this->mHelper->mFileHandle = 0;
#else
  if(this->mHelper->mDirHandle)
    closedir(this->mHelper->mDirHandle);
  this->mHelper->mFileHandle = 0;
  this->mHelper->mDirHandle = 0;
#endif
}

bool CubitDirIterator::has_next()
{
  return !this->atEnd;
}

CubitString CubitDirIterator::next()
{
  CubitString file;
  if(this->atEnd)
    return file;

#ifdef WIN32
  file = CubitString::toUtf8(this->mHelper->mDirHandle.cFileName);
  BOOL result = FindNextFileW(this->mHelper->mFileHandle, &this->mHelper->mDirHandle);
  if(result == 0)
    this->atEnd = true;
#else
  file = this->mHelper->mFileHandle->d_name;
  do
  {
    this->mHelper->mFileHandle = readdir(this->mHelper->mDirHandle);
  } while (this->mHelper->mFileHandle && 
           (mPattern.length() != 0 && 
           fnmatch(mPattern.c_str(), this->mHelper->mFileHandle->d_name, 0))
           );
  if(this->mHelper->mFileHandle == 0)
    this->atEnd = true;
#endif
  return file;
}

bool CubitFileUtil::is_directory(const CubitString& path)
{
  off_t size;
  time_t time;
  int mode;
  if(0 == file_info(path, size, time, mode))
  {
#ifdef WIN32
  if( (_S_IFDIR & mode) )
#else
  if( S_ISDIR( mode ) )
#endif
    {
      return true;
    }
  }
  return false;
}

bool CubitFileUtil::path_exists(const CubitString& path)
{
  off_t size;
  time_t time;
  int mode;
  return 0 == file_info(path, size, time, mode);
}
   
bool CubitFileUtil::is_absolute(const CubitString& path)
{
#ifdef WIN32
  return !PathIsRelativeW(CubitString::toUtf16(path).c_str());
#else
  return path.c_str()[0] == '/' ? true : false;
#endif
}
   
int  CubitFileUtil::file_info(const CubitString& path, off_t& size, time_t& time, int& mode)
{
#ifdef WIN32
  // remove trailing separators
  CubitString mypath = path;
  if(mypath.get_at(mypath.length()-1) == DIR_SEP_CHAR)
    mypath = mypath.substr(0, mypath.length()-1);
  struct _stat file_info;
  int stat_result = _wstat( CubitString::toUtf16(mypath).c_str(), &file_info );
#else
  struct stat file_info;
  int stat_result = lstat( path.c_str(), &file_info );
#endif

  if(stat_result == 0)
  {
    size = file_info.st_size;
    time = file_info.st_mtime;
    mode = file_info.st_mode;
    return 0;
  }
  return errno ? errno : ENOENT;
}

CubitFile::CubitFile()
  : mFile(NULL), mError(0)
{
}

CubitFile::CubitFile(const CubitString& file, const char* mode)
  : mFile(NULL)
{
  open(file, mode);
}

CubitFile::~CubitFile()
{
  close();
}

bool CubitFile::open(const CubitString& file, const char* mode)
{
  close();

  this->mError = 0;

#ifdef WIN32
  this->mFile = _wfopen(CubitString::toUtf16(file).c_str(), CubitString::toUtf16(mode).c_str());
#else
  this->mFile = fopen(file.c_str(), mode);
#endif
  if(!this->mFile)
  {
    this->mError = errno;
  }

  return mError == 0;
}

void CubitFile::close()
{
  if(mFile)
  {
    fclose(mFile);
  }
  mFile = NULL;
  mError = 0;
}

FILE* CubitFile::file() const
{
  return this->mFile;
}

CubitFile::operator bool () const
{
  return this->mFile != NULL;
}

int CubitFile::error()
{
  return this->mError;
}
