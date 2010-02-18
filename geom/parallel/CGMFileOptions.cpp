/**\file CGMFileOptions.cpp
 *\ copied from MOAB
 *\date 2009-06-11
 */

#include "CGMFileOptions.hpp"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

const char DEFAULT_SEPARATOR = ';';

static inline bool strempty( const char* s ) { return !*s; }

CGMFileOptions::CGMFileOptions( const char* str )
  : mData(0)
{
    // if option string is null, just return
  if (!str)
    return;
  
    // check if alternate separator is specified
  char separator[2] = { DEFAULT_SEPARATOR, '\0' };
  if (*str == DEFAULT_SEPARATOR) {
    ++str;
    if (strempty(str))
      return;
    separator[0] = *str;
    ++str;
  }
  
    // don't bother allocating copy of input string if
    // input string is empty.
  if (!strempty(str))
  {
       // tokenize at separator character
    mData = strdup( str );
    for (char* i = strtok( mData, separator ); i; i = strtok( 0, separator )) 
      if (!strempty(i)) // skip empty strings
        mOptions.push_back( i );
  }
}

CGMFileOptions::CGMFileOptions( const CGMFileOptions& copy ) :
  mData(0), mOptions( copy.mOptions.size() )
{
  if (!copy.mOptions.empty()) {
    const char* last = copy.mOptions.back();
    const char* endptr = last + strlen(last) + 1;
    size_t len = endptr - copy.mData;
    mData = (char*)malloc( len );
    memcpy( mData, copy.mData, len );
    for (size_t i = 0; i < mOptions.size(); ++i)
      mOptions[i] = mData + (copy.mOptions[i] - copy.mData);
  }
}

CGMFileOptions& CGMFileOptions::operator=( const CGMFileOptions& copy )
{
  free( mData );
  mData = 0;
  mOptions.resize( copy.mOptions.size() );

  if (!copy.mOptions.empty()) {
    const char* last = copy.mOptions.back();
    const char* endptr = last + strlen(last) + 1;
    size_t len = endptr - copy.mData;
    mData = (char*)malloc( len );
    memcpy( mData, copy.mData, len );
    for (size_t i = 0; i < mOptions.size(); ++i)
      mOptions[i] = mData + (copy.mOptions[i] - copy.mData);
  }
    
  return *this;
}

CGMFileOptions::~CGMFileOptions()
{
  free( mData );
}

CGMFOErrorCode CGMFileOptions::get_null_option( const char* name ) const
{
  const char* s;
  CGMFOErrorCode rval = get_option( name, s );
  if (FO_SUCCESS != rval)
    return rval;
  return strempty(s) ? FO_SUCCESS : FO_TYPE_OUT_OF_RANGE;
}

CGMFOErrorCode CGMFileOptions::get_int_option( const char* name, int& value ) const
{
  const char* s;
  CGMFOErrorCode rval = get_option( name, s );
  if (FO_SUCCESS != rval)
    return rval;
  
    // empty string
  if (strempty(s))
    return FO_TYPE_OUT_OF_RANGE;
  
    // parse value
  char* endptr;
  long int pval = strtol( s, &endptr, 0 );
  if (!strempty(endptr)) // syntax error
    return FO_TYPE_OUT_OF_RANGE;
  
    // check for overflow (parsing long int, returning int)
  value = pval;
  if (pval != (long int)value)
    return FO_TYPE_OUT_OF_RANGE;
  
  return FO_SUCCESS;
}

CGMFOErrorCode CGMFileOptions::get_ints_option( const char* name, 
                                          std::vector<int>& values) const
{
  const char* s;
  CGMFOErrorCode rval = get_option( name, s );
  if (FO_SUCCESS != rval)
    return rval;
  
    // empty string
  if (strempty(s))
    return FO_TYPE_OUT_OF_RANGE;
  
    // parse values
  while (!strempty(s)) {
    char* endptr;
    long int sval = strtol( s, &endptr, 0 );

#define EATSPACE(a) while ((!strncmp(a, " ", 1) ||	\
			    !strncmp(a, ",", 1)) && !strempty(a)) a++;
    //EATSPACE(endptr);

    while ((!strncmp(endptr, " ", 1) ||
	    !strncmp(endptr, ",", 1)) && !strempty(endptr)) {
      endptr++;
    }

    long int eval = sval;
    if (!strcmp(endptr, "-")) {
      endptr++;
      s = endptr;
      eval = strtol(s, &endptr, 0);
      EATSPACE(endptr);
    }
  
      // check for overflow (parsing long int, returning int)
    int value = sval;
    if (sval != (long int)value)
      return FO_TYPE_OUT_OF_RANGE;
    value = eval;
    if (eval != (long int)value)
      return FO_TYPE_OUT_OF_RANGE;
  
    for (int i = sval; i <= eval; i++)
      values.push_back(i);

    s = endptr;
  }
  
  return FO_SUCCESS;
}

CGMFOErrorCode CGMFileOptions::get_real_option ( const char* name, double& value ) const
{
  const char* s;
  CGMFOErrorCode rval = get_option( name, s );
  if (FO_SUCCESS != rval)
    return rval;
  
    // empty string
  if (strempty(s))
    return FO_TYPE_OUT_OF_RANGE;
  
    // parse value
  char* endptr;
  value = strtod( s, &endptr );
  if (!strempty(endptr)) // syntax error
    return FO_TYPE_OUT_OF_RANGE;
  
  return FO_SUCCESS;
}

CGMFOErrorCode CGMFileOptions::get_str_option( const char* name, std::string& value ) const
{
  const char* s;
  CGMFOErrorCode rval = get_option( name, s );
  if (FO_SUCCESS != rval)
    return rval;
  if (strempty(s))
    return FO_TYPE_OUT_OF_RANGE;
  value = s;
  return FO_SUCCESS;
}

CGMFOErrorCode CGMFileOptions::get_option( const char* name, std::string& value ) const
{
  const char* s;
  CGMFOErrorCode rval = get_option( name, s );
  if (FO_SUCCESS != rval)
    return rval;
  
  value = s;
  return FO_SUCCESS;
}  

CGMFOErrorCode CGMFileOptions::get_option( const char* name, const char*& value ) const
{
  std::vector<const char*>::const_iterator i;
  for (i = mOptions.begin(); i != mOptions.end(); ++i) {
    const char* opt = *i;
    if (compare( name, opt )) {
      value = opt + strlen(name);
        // if compare returned true, next char after option
        // name must be either the null char or an equals symbol.
      if (*value == '=') 
        ++value;
        
      return FO_SUCCESS;
    }
  }
  
  return FO_ENTITY_NOT_FOUND;
}

CGMFOErrorCode CGMFileOptions::match_option( const char* name, 
                                       const char* value ) const
{
  int idx;
  const char* array[] = { value, NULL };
  return match_option( name, array, idx );
}

CGMFOErrorCode CGMFileOptions::match_option( const char* name, 
                                       const char* const* values, 
                                       int& index ) const
{
  const char* optval;
  CGMFOErrorCode rval = get_option( name, optval );
  if (FO_SUCCESS != rval)
    return rval;
  
  for (index = 0; values[index]; ++index)
    if (compare( optval, values[index] ))
      return FO_SUCCESS;
  
  index = -1;
  return FO_FAILURE;
}


bool CGMFileOptions::compare( const char* name, const char* option )
{
  while (!strempty(name) && toupper(*name) == toupper(*option)) {
    ++name;
    ++option;
  }
   // match if name matched option for length of name,
   // and option either matched entirely or matches up to
   // and equals sign.
  return strempty(name) && (strempty(option) || *option == '=');
}

void CGMFileOptions::get_options( std::vector<std::string>& list ) const
{
  list.clear();
  list.resize( mOptions.size() );
  std::copy( mOptions.begin(), mOptions.end(), list.begin() );
}

#ifdef TEST

#include <iostream>

#define CHECK(A) \
  if (FO_SUCCESS != (A)) { \
    std::cerr << "Failure at line " << __LINE__ << ": error code " << (A) << std::endl; \
    return 1; \
  }

#define EQUAL(A,B) \
  if (A != B) { \
    std::cerr << "Failure at line " << __LINE__ << ": expected " << (B) << " but got " << (A) << std::endl; \
    return 2; \
  }

int main()
{
  CGMFileOptions tool( "INT1=1;NUL1;STR1=ABC;DBL1=1.0;dbl2=2.0;DBL3=3.0;INT2=2;nul2;NUL3;INT3=3;str2=once upon a time;str3==fubar=;;" );

  std::string s;
  int i;
  double d;
  CGMFOErrorCodeyy rval;
  
    // test basic get_option method without deleting entry
  rval = tool.get_option( "STR1", s );
  CHECK(rval);
  EQUAL( s, "ABC" );
  
    // test basic get_option method again, this time deleting the entry
  rval = tool.get_option( "STR1", s );
  CHECK(rval);
  EQUAL( s, "ABC" );
  
    // test basig get_option method with a null option
  rval = tool.get_option( "NUL2", s );
  CHECK( rval );
  EQUAL( s.empty(), true );

  
    // test null option
  rval = tool.get_null_option( "nul1" );
  CHECK( rval );
  
    // try null option method on non-null value
  rval = tool.get_null_option( "INT1" ) ;
  EQUAL( rval, FO_TYPE_OUT_OF_RANGE) ;
  

    // test integer option
  rval = tool.get_int_option( "int1", i );
  CHECK( rval );
  EQUAL( i, 1 );
  
  rval = tool.get_int_option( "int2", i );
  CHECK( rval );
  EQUAL( i, 2 );
  
    // test integer option on non-integer value
  rval = tool.get_int_option( "dbl2", i );
  EQUAL( rval, FO_TYPE_OUT_OF_RANGE );
  
    // test integer option on null value
  rval = tool.get_int_option( "NUL3", i);
  EQUAL( rval, FO_TYPE_OUT_OF_RANGE );
  
    // test double option
  rval = tool.get_real_option( "dbl1", d );
  CHECK( rval );
  EQUAL( d, 1.0 );
  
  rval = tool.get_real_option( "dbl2", d );
  CHECK( rval );
  EQUAL( d, 2.0 );
  
  rval = tool.get_real_option( "int3", d );
  CHECK( rval );
  EQUAL( d, 3.0 );
  
    // test real option on non-real value
  rval = tool.get_real_option( "str2", d );
  EQUAL( rval, FO_TYPE_OUT_OF_RANGE );
  
  
    // test real option on null value
  rval = tool.get_real_option( "NUL3", d );
  EQUAL( rval, FO_TYPE_OUT_OF_RANGE );
  
    // test get a simple string option
  rval = tool.get_str_option( "DBL3", s );
  CHECK( rval );
  EQUAL( s, "3.0" );
  
    // test get a string with spaces
  rval = tool.get_str_option("STR2", s );
  CHECK( rval );
  EQUAL( s, "once upon a time" );
  
    // try to get a string value for a null option
  rval = tool.get_str_option( "nul3", s );
  EQUAL( rval, FO_TYPE_OUT_OF_RANGE );
  
    // test options using generic get_option method
    
  rval = tool.get_option( "NUL3", s );
  CHECK( rval );
  EQUAL( s.empty(), true );
  
  rval = tool.get_option( "STR3", s );
  CHECK( rval );
  EQUAL( s, "=fubar=" );
  
    // test size of options string
  unsigned l = tool.size();
  EQUAL( l, 12u );
  
  
    // test alternate separator
  
  CGMFileOptions tool2( ";+OPT1=ABC+OPT2=" );
  l = tool2.size();
  EQUAL( l, 2 );
  
  rval = tool2.get_option( "opt1", s );
  CHECK( rval );
  EQUAL( s, "ABC" );
  
  rval = tool2.get_option( "opt2", s );
  CHECK( rval );
  bool e = s.empty();
  EQUAL( e, true );
  
  l = tool2.size();
  EQUAL( l, 2 );
  
    
    // test empty options string
    
  CGMFileOptions tool3( ";;;;" );
  e = tool3.empty();
  EQUAL( e, true );
  l = tool3.size();
  EQUAL( l, 0 );
  
  CGMFileOptions tool4(NULL);
  e = tool4.empty();
  EQUAL( e, true );
  l = tool4.size();
  EQUAL( l, 0 );
  
  CGMFileOptions tool5(";+");
  e = tool5.empty();
  EQUAL( e, true );
  l = tool5.size();
  EQUAL( l, 0 );
  
    // test copy constructor
  
  CGMFileOptions tool6( tool2 );
  
  rval = tool6.get_option( "opt1", s );
  CHECK( rval );
  EQUAL( s, "ABC" );
  
  rval = tool6.get_option( "opt2", s );
  CHECK( rval );
  e = s.empty();
  EQUAL( e, true );
  
  l = tool6.size();
  EQUAL( l, 2 );
  
  CGMFileOptions tool7( tool5 );
  e = tool7.empty();
  EQUAL( e, true );
  l = tool7.size();
  EQUAL( l, 0 );
  
    // test assignment operator
  
  CGMFileOptions tool8( tool2 );
  tool8 = tool;
  EQUAL( tool8.size(), tool.size() );
    
  return 0;
}

#endif
