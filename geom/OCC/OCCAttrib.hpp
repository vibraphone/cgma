#ifndef FACET_ATTRIB_HPP
#define FACET_ATTRIB_HPP

#include "CubitDefines.h"
#include "CubitString.hpp"
class CubitSimpleAttrib;
class CubitString;

class OCCAttrib
{
public:

  OCCAttrib( CubitSimpleAttrib* );
  
  ~OCCAttrib();
  
  CubitSimpleAttrib* get_CSA() const;
  
  bool equals( CubitSimpleAttrib* ) const;
  
  CubitStatus save( FILE* file ) const;
  
  static OCCAttrib* restore( FILE* file, unsigned int endian );

  CubitString name() const
    { return stringArray[0]; }

    // Use arrays rather than std::vector so IO routines
    // can use the data without making copies.
  CubitString* stringArray;
  double* doubleArray;
  int* integerArray;
  
  int numStrings;
  int numDoubles;
  int numIntegers;
  
  OCCAttrib* listNext;




private:

  OCCAttrib( int string_count, CubitString* strings,
               int double_count, double* doubles,
               int integer_count, int* integers );

};

#endif
