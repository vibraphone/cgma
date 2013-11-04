//- Filename:    CubitSimpleAttrib.cc
//- Description: Implementation of class CubitSimpleAttrib
//- Owner:       Greg Nielson
//- Checked by:
//- Version:

#include <string.h>
#include "SettingHandler.hpp"

#include "CubitSimpleAttrib.hpp"
#include "CubitString.hpp"
#include "RefEntity.hpp"
#include "RefEntityFactory.hpp"
#include "CubitFileIOWrapper.hpp"

CubitBoolean CubitSimpleAttrib::pushAttribs = CUBIT_FALSE;

//SettingHolder SettingHandler::settingObject43 = SettingHolder("Push_Attribs", CubitSimpleAttrib::set_push_attribs, CubitSimpleAttrib::get_push_attribs, CUBIT_FALSE);

CubitSimpleAttrib::CubitSimpleAttrib()
{
}

CubitSimpleAttrib::~CubitSimpleAttrib()
{
}

CubitSimpleAttrib::CubitSimpleAttrib(const CubitString new_character_type,
                                     const CubitString new_string_data,
                                     const CubitString new_more_string_data,
                                     const int* new_integer_data,
                                     const double* new_double_data)
{
  assert(new_character_type.length() > 0);

  stringDataList.push_back(new_character_type);

  if (new_string_data.length()) {
    stringDataList.push_back(new_string_data);
  }
  
  if (new_more_string_data.length()) {
    stringDataList.push_back(new_more_string_data);
  }
  
  if(new_double_data)
    doubleDataList.push_back( *new_double_data );
  if(new_integer_data)
    intDataList.push_back( *new_integer_data );

  return;
}

bool CubitSimpleAttrib::isEmpty() const
{
  return stringDataList.empty() &&
      doubleDataList.empty() &&
      intDataList.empty();
}

CubitSimpleAttrib::CubitSimpleAttrib(const CubitSimpleAttrib& csa_ptr)
{
  this->stringDataList = csa_ptr.stringDataList;
  this->doubleDataList = csa_ptr.doubleDataList;
  this->intDataList = csa_ptr.intDataList;
}

CubitSimpleAttrib::CubitSimpleAttrib(const std::vector<CubitString> *string_list,
                                     const std::vector<double> *double_list,
                                     const std::vector<int> *int_list)
{
  if(string_list)
    this->stringDataList = *string_list;
  if(double_list)
    this->doubleDataList = *double_list;
  if(int_list)
    this->intDataList = *int_list;
}


CubitSimpleAttrib &CubitSimpleAttrib::operator=(const CubitSimpleAttrib &other)
{
    this->stringDataList = other.stringDataList;
    this->intDataList = other.intDataList;
    this->doubleDataList = other.doubleDataList;
    return *this;
}

CubitString CubitSimpleAttrib::character_type() const //- returns the type of attribute
{
  if(stringDataList.size() >= 1 )
  {
    return stringDataList[0];
  }
  return CubitString();
}

void CubitSimpleAttrib::string_data_list( const std::vector<CubitString>& new_list)
{
  this->stringDataList = new_list;
} 

void CubitSimpleAttrib::double_data_list(const std::vector<double>& new_list)
{
  doubleDataList = new_list;
}

void CubitSimpleAttrib::int_data_list(const std::vector<int>& new_list)
{
  intDataList = new_list;
}

bool CubitSimpleAttrib::operator==(const CubitSimpleAttrib& other) const
{
  return this->stringDataList == other.stringDataList &&
         this->doubleDataList == other.doubleDataList &&
         this->intDataList == other.intDataList;
}


//Initialize all settings in this class
void CubitSimpleAttrib::initialize_settings()
{

  SettingHandler::instance()->add_setting("Push Attribs", 
                                         CubitSimpleAttrib::set_push_attribs, 
					 CubitSimpleAttrib::get_push_attribs); 
}

void CubitSimpleAttrib::print() const
{
  PRINT_INFO("CSA: type = %s\n", stringDataList[0].c_str());

  PRINT_INFO("String data: ");
  for ( size_t i=1; i<stringDataList.size(); i++)
    PRINT_INFO("%s;", stringDataList[i].c_str());
  if (stringDataList.size() == 0)
    PRINT_INFO("(none)");
  PRINT_INFO("\n");
  
  PRINT_INFO("Int data: ");
  for (size_t i=0; i<intDataList.size(); i++)
    PRINT_INFO("%d;", intDataList[i]);
  if (intDataList.size() == 0)
    PRINT_INFO("(none)");
  PRINT_INFO("\n");
  
  PRINT_INFO("Double data: ");
  for (size_t i=0; i<doubleDataList.size(); i++)
    PRINT_INFO("%f;", doubleDataList[i]);
  if (doubleDataList.size() == 0)
    PRINT_INFO("(none)");
  PRINT_INFO("\n");
}

