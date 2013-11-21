#include "FacetAttrib.hpp"
#include "CubitString.hpp"
#include "CubitSimpleAttrib.hpp"
#include "CubitFileIOWrapper.hpp"
#include <algorithm>

  // Constructor - copy from CubitSimpleAttrib
FacetAttrib::FacetAttrib( const CubitSimpleAttrib& csa ) : listNext(0)
{
  int i;

    // save counts
  numStrings = csa.string_data_list().size();
  numDoubles = csa.double_data_list().size();
  numIntegers = csa.int_data_list().size();

    // allocate arrays, but don't try to allocate zero-length arrays
  stringArray = numStrings ? new CubitString[numStrings] : NULL;
  doubleArray = numDoubles ? new double[numDoubles] : NULL;
  integerArray = numIntegers ? new int[numIntegers] : NULL;

    // copy data into arrays
  for( i = 0; i < numStrings; i++ )
    stringArray[i] = csa.string_data_list()[i];
  for( i = 0; i < numIntegers; i++ )
    integerArray[i] = csa.int_data_list()[i];
  for( i = 0; i < numDoubles; i++ )
      doubleArray[i] = csa.double_data_list()[i];
}

  // Private constructor for use by restore(FILE*)
FacetAttrib::FacetAttrib( int string_count, CubitString strings[],
                          int double_count, double doubles[],
                          int int_count, int integers[] )
: stringArray(strings), doubleArray(doubles), integerArray(integers),
  numStrings(string_count),  numDoubles(double_count), numIntegers(int_count),
  listNext(0)
{}


  // Destructor -- free arrays
FacetAttrib::~FacetAttrib()
{
    // "delete"ing NULL pointers is okay.
  delete [] integerArray;
  delete [] doubleArray;
  delete [] stringArray;
}

  // Copy this into a new CubitSimpleAttrib
CubitSimpleAttrib FacetAttrib::get_CSA() const
{
    // Set initial list size
  std::vector<CubitString> string_list(numStrings);
  std::vector<int> int_list(numIntegers);
  std::vector<double> double_list(numDoubles);

    // Don't need to 'new' objects in DLIList because
    // CSA will make copies.  Just put addresses in list.
  int i;
  for( i = 0; i < numStrings; i++ )
    string_list[i] = stringArray[i];
  for( i = 0; i < numIntegers; i++ )
    int_list[i] = integerArray[i];
  for( i = 0; i < numDoubles; i++ )
    double_list[i] = doubleArray[i];

  return CubitSimpleAttrib( &string_list, &double_list, &int_list );
}

  // compare to a CubitSimpleAttrib
bool FacetAttrib::equals( const CubitSimpleAttrib& csa ) const
{
  // compare counts
  if( csa.int_data_list().size() != (int)numIntegers ||
      csa.double_data_list().size() != (int)numDoubles ||
      csa.string_data_list().size() != (int)numStrings )
    return false;

    // compare strings first because most likely the
    // first string (the name) will differ.
  if(!std::equal(stringArray, stringArray+numStrings, csa.string_data_list().begin()))
    return false;
  if(!std::equal(doubleArray, doubleArray+numDoubles, csa.double_data_list().begin()))
    return false;
  if(!std::equal(integerArray, integerArray+numIntegers, csa.int_data_list().begin()))
    return false;

  return true;
}

  // write to a file at the current file offset
CubitStatus FacetAttrib::save(FILE *save_file) const
{
  if( save_file == NULL)
  {
    PRINT_ERROR("Problem saving MBG attributes: null FILE ptr\n");
    return CUBIT_FAILURE;
  }

  NCubitFile::CIOWrapper wrapper(save_file);

  // write a version number for the attribute data
  unsigned int Attrib_Version = 1;
  wrapper.Write(&Attrib_Version, 1);

  // write the number of strings, number of doubles, and number of integers
  int counts[3] = { numStrings, numDoubles, numIntegers };
  wrapper.Write(reinterpret_cast<unsigned int*>(counts), 3);

  // write the string data
  int i;
  for( i = 0; i < numStrings; i++ )
    wrapper.Write(stringArray[i].c_str());

  // write the doubles
  wrapper.Write(doubleArray, numDoubles);

  // write the integers
  wrapper.Write(reinterpret_cast<unsigned int*>(integerArray), numIntegers);
  
  return CUBIT_SUCCESS;
}


  // read from file starting at current file offset
FacetAttrib* FacetAttrib::restore(FILE *restore_file, unsigned int endian)
{
  if( restore_file == NULL )
    return NULL;

  NCubitFile::CIOWrapper wrapper(endian, restore_file );

  // write a version number for the attribute data
  unsigned int version;
  wrapper.Read(&version, 1);

  // haven't handled any version changes yet
  if( version != 1 )
  {
    PRINT_ERROR("Wrong FacetAttrib version : %u\n", version );
    return NULL;
  }

  // read the number of strings, number of doubles, and number of integers
  int counts[3];
  wrapper.Read(reinterpret_cast<unsigned int*>(counts), 3);
  int n_strings = counts[0];
  int n_doubles = counts[1];
  int n_ints = counts[2];
  
    // allocate arrays, but don't try to allocate zero-length array
  CubitString* strings = n_strings ? new CubitString[n_strings] : NULL;
  double *doubles = n_doubles ? new double[n_doubles] : NULL;
  int *ints = n_ints ? new int[n_ints] : NULL;
  
  // read the string data
  int i;
  for( i = 0; i < n_strings; i++ )
    strings[i] = CubitString(wrapper.Read());

  // write the doubles
  wrapper.Read(doubles, n_doubles);

  // write the integers
  wrapper.Read(reinterpret_cast<unsigned int*>(ints), n_ints);

  return new FacetAttrib(n_strings, strings, n_doubles, doubles, n_ints, ints);
}


