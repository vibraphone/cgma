/**\file CGMFileOptions.hpp
 *\copied from MOAB
 *\date 2009-06-11
 */

#ifndef CGM_FILE_OPTIONS_HPP
#define CGM_FILE_OPTIONS_HPP

#include <string>
#include <vector>
#include "CubitDefines.h"

/* file option type */
enum CGMFOErrorCode
{
  FO_SUCCESS = 0,
  FO_TYPE_OUT_OF_RANGE,
  FO_ENTITY_NOT_FOUND,
  FO_FAILURE
};

/**\brief Parse options string passed to file IO routines
 *
 * This is a utility class used by file-IO-related code to 
 * parse the options string passed to ParallelMeshTool::load_file
 */
class CGMFileOptions {
public:

  /*\param options_string The concatenation of a list of 
   *          options, separated either by the default separator
   *          (semicolon) with a custom separator specified at
   *          the beginning of the string (semicolon followed by
   *          destired separator character.)
   */
  CGMFileOptions( const char* option_string );
  
  CGMFileOptions( const CGMFileOptions& copy );
  CGMFileOptions& operator=( const CGMFileOptions& copy );
  
  ~CGMFileOptions();
  
  /**\brief Check for option with no value 
   *
   * Check for an option w/out a value.
   *\param name The option name
   *\return - CUBIT_SUCCESS if option is found
   *        - CUBIT_TYPE_OUT_OF_RANGE if options is found, but has value
   *        - CUBIT_ENTITY_NOT_FOUND if option is not found.
   */
  CGMFOErrorCode get_null_option( const char* name ) const;
  
  /**\brief Check for option with an integer value.
   *
   * Check for an option with an integer value
   *\param name The option name
   *\param value Output. The value.
   *\return - CUBIT_SUCCESS if option is found
   *        - CUBIT_TYPE_OUT_OF_RANGE if options is found, but does not have an integer value
   *        - CUBIT_ENTITY_NOT_FOUND if option is not found.
   */
  CGMFOErrorCode get_int_option( const char* name, int& value ) const;
  
  /**\brief Check for option with a double value.
   *
   * Check for an option with a double value
   *\param name The option name
   *\param value Output. The value.
   *\return - CUBIT_SUCCESS if option is found
   *        - CUBIT_TYPE_OUT_OF_RANGE if options is found, but does not have a double value
   *        - CUBIT_ENTITY_NOT_FOUND if option is not found.
   */
  CGMFOErrorCode get_real_option( const char* name, double& value ) const;
  
  /**\brief Check for option with any value.
   *
   * Check for an option with any value.
   *\param name The option name
   *\param value Output. The value.
   *\return - CUBIT_SUCCESS if option is found
   *        - CUBIT_TYPE_OUT_OF_RANGE if options is found, but does not have a value
   *        - CUBIT_ENTITY_NOT_FOUND if option is not found.
   */
  CGMFOErrorCode get_str_option( const char* name, std::string& value ) const;
  
  /**\brief Check for option 
   *
   * Check for an option
   *\param name The option name
   *\param value The option value, or an empty string if no value.
   *\return CUBIT_SUCCESS or CUBIT_ENTITY_NOT_FOUND
   */
  CGMFOErrorCode get_option( const char* name, std::string& value ) const;
  
  /**\brief Check the string value of an option
   *
   * Check which of a list of possible values a string option contains.
   *\param name The option name
   *\param values A NULL-terminated array of C-style strings enumerating
   *              the possible option values.
   *\param index  Output: The index into <code>values</code> for the
   *              option value.
   *\return CUBIT_SUCCESS if matched name and value.
   *        CUBIT_ENTITY_NOT_FOUND if the option was not specified
   *        CUBIT_FAILURE if the option value is not in the input <code>values</code> array.
   */
  CGMFOErrorCode match_option( const char* name, const char* const* values, int& index ) const;
  
  /**\brief Check if an option matches a string value
   *
   * Check if the value for an option is the passed string.
   *\param name The option name
   *\param value The expected value.
   *\return CUBIT_SUCCESS if matched name and value.
   *        CUBIT_ENTITY_NOT_FOUND if the option was not specified
   *        CUBIT_FAILURE if the option value doesn't match the passed string/
   */
  CGMFOErrorCode match_option( const char* name, const char* value ) const;
  
  /**\brief Check for option for which the value is a list of ints
   *
   * Check for an option which is an int list.  The value is expected to
   * be a comma-separated list of int ranges, where an int range can be 
   * either a single integer value or a range of integer values separated
   * by a dash ('-').
   *
   *\param name The option name
   *\param values Output. The list of integer values.
   *\return - CUBIT_SUCCESS if option is found
   *        - CUBIT_TYPE_OUT_OF_RANGE if options is found, but does not contain an ID list
   *        - CUBIT_ENTITY_NOT_FOUND if option is not found.
   */
  CGMFOErrorCode get_ints_option( const char* name, std::vector<int>& values) const;
  
  /** number of options */
  inline unsigned size() const 
    { return mOptions.size(); }
  
  /** true if no options */
  inline bool empty() const 
    { return mOptions.empty(); }
  
  /** Get list of options */
  void get_options( std::vector<std::string>& list ) const;
  
private:
  
  /**\brief Check for option 
   *
   * Check for an option
   *\param name The option name
   *\param value The option value, or an empty string if no value.
   *\return CUBIT_SUCCESS or CUBIT_ENTITY_NOT_FOUND
   */
  CGMFOErrorCode get_option( const char* name, const char*& value) const;

  char* mData;
  std::vector<const char*> mOptions;

    /** Case-insensitive compare of name with option value. */
  static bool compare( const char* name, const char* option );
};

#endif

