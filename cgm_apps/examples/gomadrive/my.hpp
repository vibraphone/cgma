
//--------------------------
//
// File:  my.hpp
//
// description:
//     generic header file
//
// ------------------------------
#ifndef MY_DEFINES
#define MY_DEFINES

#ifndef CANT_USE_STD_IO
#include <iostream>
#include <fstream>
#else
#include <iostream.h>
#include <fstream.h>
#endif

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>



#define TRUE 1
#define FALSE 0
enum body_type
{
    BODY_NONE = 0,
    VERTEX = 1,
    EDGE = 2,
    FACE = 3,
    SPHERE = 4,
    POLYGON = 5,
    BTEXT = 6,
    FTEXT = 7
};

inline enum body_type String_to_enum(char* type)
{
    //return the type of geometry
  if(strcasecmp(type,"vertex")==0)
      return VERTEX;
  else if(strcasecmp(type,"edge")==0)
      return EDGE;
  else if(strcasecmp(type,"face")==0)
      return FACE;
  else if(strcasecmp(type,"sphere")==0)
      return SPHERE;
  else if(strcasecmp(type,"polygon")==0)
      return POLYGON;
  else
      return BODY_NONE;
}
inline enum body_type String_to_enum(const char* type)
{
    //return the type of geometry
  if(strcasecmp(type,"vertex")==0)
      return VERTEX;
  else if(strcasecmp(type,"edge")==0)
      return EDGE;
  else if(strcasecmp(type,"face")==0)
      return FACE;
  else if(strcasecmp(type,"sphere")==0)
      return SPHERE;
  else if(strcasecmp(type,"polygon")==0)
      return POLYGON;
  else
      return BODY_NONE;
}

#endif
