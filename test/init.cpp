/**
 * \file init.cpp
 *
 * \brief Tests of CGM initialization / shutdown features
 */
#include "InitCGMA.hpp"

#include <typeinfo>


// main program - initialize, then send to proper function
int main (int argc, char **argv)
{
  // first initialization, which should succeed
  CubitStatus result = InitCGMA::initialize_cgma();
  if (CUBIT_SUCCESS != result) return 1;

  // second initialization, which should also succeed
  result = InitCGMA::initialize_cgma();
  if(CUBIT_SUCCESS != result) return 1;

  // now try again, with a different argument; should fail.
  result = InitCGMA::initialize_cgma( "facet" );
  if(CUBIT_FAILURE != result) return 1;

  return 0;
  
}


    
