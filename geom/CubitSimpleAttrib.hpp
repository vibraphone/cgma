//- Class: CubitSimpleAttrib
//- Owner: Greg Nielson
//- Description: This class is used to carry information from
//-              the cubit attribute classes to the solid model
//-               attributes.
//- Checked By:
//- Version: $Id:

#ifndef CUBIT_SIMPLE_ATTRIB_HPP
#define CUBIT_SIMPLE_ATTRIB_HPP

#if defined(WIN32)
#pragma warning(push)
#pragma warning(disable : 4251)  // hide warnings about template dll exports
#endif

#include "CubitString.hpp"
#include "CubitDefines.h"
#include <vector>
#include "CubitGeomConfigure.h"
#include <new>

class CUBIT_GEOM_EXPORT CubitSimpleAttrib {
private:

  std::vector<CubitString> stringDataList;  // List of CubitString attribute data

  std::vector<double> doubleDataList;  // List of double attribute data

  std::vector<int> intDataList;  // List of int attribute data

  static CubitBoolean pushAttribs;
    //- if true, when saving attributes, push them to non-primary bridges
    //- on merged ref entities
  
public:

  // returns whether this CSA is valid (it is invalid if the data is empty)
  bool isEmpty() const;

  CubitSimpleAttrib();
  //- default constructor

  explicit CubitSimpleAttrib(const CubitString new_character_type,
                             const CubitString new_string_data = CubitString(),
                             const CubitString new_more_string_data = CubitString(),
                             const int* new_integer_data = NULL,
                             const double* new_double_data = NULL);
  //- constructor

  explicit CubitSimpleAttrib(const std::vector<CubitString> *string_list,
                             const std::vector<double> *double_list = NULL,
                             const std::vector<int> *int_list = NULL);
  //- constructor  

  CubitSimpleAttrib(const CubitSimpleAttrib& csa_ptr);
    //- copy constructor, sort of
  
  ~CubitSimpleAttrib();
   //- destructor

  CubitSimpleAttrib& operator=(const CubitSimpleAttrib &other);

  CubitString character_type() const;
  //- returns the type of attribute from the first item in the string list
  
  static void initialize_settings();

  const std::vector<CubitString>& string_data_list() const {return stringDataList;}
  std::vector<CubitString>& string_data_list() {return stringDataList;}

  void string_data_list( const std::vector<CubitString>& new_list);

  const std::vector<double>& double_data_list() const {return doubleDataList;}
  std::vector<double>& double_data_list() {return doubleDataList;}

  void double_data_list(const std::vector<double>& new_list);

  const std::vector<int>& int_data_list() const {return intDataList;}
  std::vector<int>& int_data_list() {return intDataList;}

  void int_data_list(const std::vector<int>& new_list);

  
  bool operator==(const CubitSimpleAttrib& other) const;
    //- return true if the two csa's are equivalent

  void print() const;
    //- print the contents of this simpleattrib

  static CubitBoolean get_push_attribs() {return pushAttribs;};
  static void set_push_attribs(CubitBoolean flag) {pushAttribs = flag;};
    //- get/set pushAttribs
  
};


#if defined(WIN32)
#pragma warning(pop)
#endif

#endif








