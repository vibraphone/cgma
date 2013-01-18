//-------------------------------------------------------------------------
// Filename      : OCCAttribSet.cpp
//
// Purpose       : Common attribute functionality for MBG
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/01/03
//-------------------------------------------------------------------------

#include "OCCAttribSet.hpp"
//#include "OCCAttrib.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitSimpleAttrib.hpp"
#include "CubitFileIOWrapper.hpp"

#ifndef OCC_VERSION_MINOR
#include "Standard_Version.hxx"
#endif

#if OCC_VERSION_MINOR < 5
  #include "TDataStd_Shape.hxx"
  typedef TDataStd_Shape TDataXtd_Shape;
  typedef Handle_TDataStd_Shape Handle_TDataXtd_Shape;
#else
  #include "TDataXtd_Shape.hxx"
#endif

#include "TDF_Label.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "TCollection_ExtendedString.hxx"
#include "Handle_TDataStd_Name.hxx"
#include "Handle_TDataStd_ExtStringArray.hxx"
#include "Handle_TDataStd_RealArray.hxx"
#include "Handle_TDataStd_IntegerArray.hxx"
#include "TDataStd_Name.hxx"
#include "TDataStd_ExtStringArray.hxx"
#include "TDataStd_RealArray.hxx"
#include "TDataStd_IntegerArray.hxx"
#include "TopoDS_Shape.hxx"
#include "TDF_ChildIterator.hxx"
#include <vector>
typedef std::map<int, TDF_Label>::value_type labType;

void OCCAttribSet::FindShape(TopoDS_Shape& shape,
                             TDF_Label& aLabel,
                             CubitBoolean& found) 
{
  found = CUBIT_FALSE;
  if(!OCCQueryEngine::instance()->OCCMap->IsBound(shape))
    return;

  int k = OCCQueryEngine::instance()->OCCMap->Find(shape);
  std::map<int, TDF_Label>::iterator it = 
            OCCQueryEngine::instance()->Shape_Label_Map->find(k);
  if(it != OCCQueryEngine::instance()->Shape_Label_Map->end())
  {
    aLabel = (*it).second; 
    found = CUBIT_TRUE;
  }
}

void OCCAttribSet::append_attribute( const CubitSimpleAttrib& csa, TopoDS_Shape& shape )
{
  if(shape.IsNull())
    return;

  //Add attributes on child of myLabel
  //1. add shape attribute, first check to make sure there's no shape attribute
  CubitBoolean found = CUBIT_FALSE;
  int result = 1;
  TDF_Label aLabel, lab;
  Handle_TDataStd_Name attr_name;
  CubitString type = csa.string_data_list()[0];
  TCollection_ExtendedString cstring( (Standard_CString)type.c_str() );

  FindShape(shape, aLabel, found);

  if(!found)
  { 
    aLabel = OCCQueryEngine::instance()->mainLabel.NewChild();
    Handle_TDataXtd_Shape attr_shape = TDataXtd_Shape::Set(aLabel, shape);
    //aLabel.AddAttribute(attr_shape);
    //test if the attribute has been attached
    assert(aLabel.IsAttribute(TDataXtd_Shape::GetID()));

    if(!OCCQueryEngine::instance()->OCCMap->IsBound(shape))
    {
      OCCQueryEngine::instance()->iTotalTBCreated++;
      int k = OCCQueryEngine::instance()->iTotalTBCreated;
      OCCQueryEngine::instance()->OCCMap->Bind(shape, k);
    }
 
    int k = OCCQueryEngine::instance()->OCCMap->Find(shape);
    OCCQueryEngine::instance()->Shape_Label_Map->insert(labType(k, aLabel));
  }

  //2. add type attribute , below attributes are added on child lable of myLabel
  //First check if the attributes are already exist, if so, create new data.
  else
  {
    TCollection_ExtendedString name_string;
    for (TDF_ChildIterator it(aLabel,CUBIT_FALSE); it.More(); it.Next()) 
    {
      lab = it.Value();
      result = find_attribute(lab, csa);
      if(result != 1)
        break;
    }
  }
  //if not exist,  create new child label
  Handle_TDataStd_ExtStringArray attr_string;
  const std::vector<CubitString>& strings = csa.string_data_list();
  int size = strings.size()-1;
 
  if (result == 1) 
  {
    lab = aLabel.NewChild();
    attr_name = TDataStd_Name::Set(lab, cstring);
    //lab.AddAttribute(attr_name);
    assert(lab.IsAttribute(TDataStd_Name::GetID()));
  }
     
  else if(result == 3)//nothing to do
    return;

  else if(result == 2) 
  //adding the name strings
  {
    lab.FindAttribute(TDataStd_ExtStringArray::GetID(), attr_string);
    TCollection_ExtendedString name1, name2;
    std::vector<CubitString> new_names;
    int length = 0;
    for(int j = 1; j < (int)strings.size();j++)
    { 
      CubitString name_string = strings[j];
      TCollection_ExtendedString cstring(name_string.c_str(), Standard_True );
      bool dup = false;
      length = attr_string->Length() ;

      for(int i = 0; i < length; i++)
      {
        name1 = attr_string->Value(i);
        //it could be string[0] = "volume 1", string[1] = "volume_1";
    
        if(name1.IsEqual(cstring))
        {
          dup = true;
          break;
        }
      }
      if(!dup) 
        new_names.push_back(name_string);
    }
    if(new_names.size() > 0)
    {
      lab.ForgetAttribute(attr_string);
      Handle_TDataStd_ExtStringArray attr_string_new; 
      attr_string_new = TDataStd_ExtStringArray::Set(lab, 0, 
                                       length+new_names.size()-1); 
      for(int i = 0; i < attr_string_new->Length(); i++)
      {
        if (i < length)
          attr_string_new->SetValue(i, attr_string->Value(i));
        else
        {
          TCollection_ExtendedString cstring(new_names[i-length].c_str(), Standard_True );
          attr_string_new->SetValue(i, cstring);
        }
      }
    } 
    return;
  } 
  //3. add string attribute
  if(size > 0)
  {
    //set the length of String Array .
    attr_string = TDataStd_ExtStringArray::Set(lab, 0, size-1);
  }
  for(int i = 1; i < (int)strings.size(); i++)
  {
    TCollection_ExtendedString cstring(strings[i].c_str(), Standard_True );
    attr_string->SetValue(i-1, cstring) ;
  }
  
  if (size > 0)
  {
    //lab.AddAttribute(attr_string);
    assert(lab.IsAttribute(TDataStd_ExtStringArray::GetID()));
  }
  //4. add double attribute
  const std::vector<double>& doubles = csa.double_data_list();
  Handle_TDataStd_RealArray attr_double;
  size = doubles.size();
  if(size > 0)
  {
    //set the length of double array .
    attr_double = TDataStd_RealArray::Set(lab,0, size-1);
  } 
  for(int i = 0; i < size; i++)
    attr_double->SetValue(i, doubles[i]);

  if(size > 0)
  {
    //lab.AddAttribute(attr_double);
    assert(lab.IsAttribute(TDataStd_RealArray::GetID()));
  }  
  //5. add int attribute
  const std::vector<int>& ints = csa.int_data_list();
  Handle_TDataStd_IntegerArray attr_int;
  size = ints.size();
  if(size > 0)
  {
    //set the length of int array be 11, can be extended.
    attr_int = TDataStd_IntegerArray::Set(lab, 0, size-1);
  }
  for(int i = 0; i < size; i++)
    attr_int->SetValue(i, ints[i]);

  if(size > 0)
  {
    //lab.AddAttribute(attr_int);
    assert(lab.IsAttribute(TDataStd_IntegerArray::GetID()));
  }
}

// return 1 = not found
// return 2 = found ENTITY_NAME attribute
// return 3 = found the same attribute
int  OCCAttribSet::find_attribute(TDF_Label child,
                                  const CubitSimpleAttrib& csa)
{
  const std::vector<int>& ints = csa.int_data_list();
  CubitString type = csa.character_type();
  Standard_Boolean isMultiByte = Standard_True;
  TCollection_ExtendedString cstring( type.c_str(), isMultiByte );
  const std::vector<double>& doubles = csa.double_data_list();
  const std::vector<CubitString>& strings = csa.string_data_list();

  //find the same type attribute first
  Handle_TDataStd_Name attr_name;
  TCollection_ExtendedString old_string;
  if(child.FindAttribute(TDataStd_Name::GetID(), attr_name))
  {
    old_string = attr_name->Get();
    Standard_Boolean isMultiByte = Standard_True;
    char const* string = "ENTITY_NAME";
    TCollection_ExtendedString string2( string, isMultiByte );
    if(old_string != cstring)
      return 1;

    else if(cstring.IsEqual(string2))
      return 2;  
  }

  else
    return 1;

  //continue to compare the rest attributes.
  CubitBoolean is_same = CUBIT_TRUE;
  if(ints.size()> 0 )
  {
    Handle_TDataStd_IntegerArray attr_ints;
    if(child.FindAttribute(TDataStd_IntegerArray::GetID(), attr_ints) &&
       attr_ints->Length() == (int)ints.size())
    {
      for(int i = 0; i < (int)ints.size(); i++)
        if(attr_ints->Value(i) != ints[i])
        {
          is_same = CUBIT_FALSE;
          break;
        }
    }
  }
  if(!is_same)
    return 1;

  if(doubles.size() > 0 )
  {
    Handle_TDataStd_RealArray attr_doubles;
    if(child.FindAttribute(TDataStd_RealArray::GetID(), attr_doubles) &&
       attr_doubles->Length() == (int)doubles.size())
    {
      for(int i = 0; i < (int)doubles.size(); i++)
        if(attr_doubles->Value(i) != doubles[i])
        {
          is_same = CUBIT_FALSE;
          break;
        }
    }
  }

  if(!is_same)
    return 1;

  if(strings.size() > 0 )
  {
    Handle_TDataStd_ExtStringArray attr_strings;
    if(child.FindAttribute(TDataStd_ExtStringArray::GetID(), attr_strings) &&
       attr_strings->Length() == (int)strings.size()-1)
    {
      for(int i = 1; i < (int)strings.size(); i++)
      {
        CubitString astring = strings[i];
        TCollection_ExtendedString string( astring.c_str() , Standard_True);
        if(attr_strings->Value(i-1) != string)
        {
          is_same = CUBIT_FALSE;
          break;
        }
      }
    }
    if(!is_same)
      return 1;
  }
  return 3  ;
}

void OCCAttribSet::remove_attribute( const CubitSimpleAttrib& csa)
{
  CubitString type = csa.character_type();
  if(type.length() == 0)
      return;

  TDF_Label myLabel;
  for (TDF_ChildIterator it(OCCQueryEngine::instance()->mainLabel, CUBIT_FALSE); it.More(); it.Next())
  {
    //find the same type attribute first
    myLabel = it.Value();

    //forget csa attribute from the document
    for (TDF_ChildIterator it1(myLabel, CUBIT_FALSE); it1.More(); it1.Next())
    {
      //find the same type attribute first
      TDF_Label child = it1.Value();

      int found = find_attribute(child, csa);

      if(found == 3)
        child.ForgetAllAttributes(false );
    }
  }
}

void OCCAttribSet::remove_attribute( TopoDS_Shape& shape)
{
  CubitBoolean found = CUBIT_FALSE;
  TDF_Label myLabel;

  FindShape(shape, myLabel, found);

  if(!found)
    return;
 
  myLabel.ForgetAllAttributes(false);
/*
  if(OCCQueryEngine::instance()->OCCMap->IsBound(shape)) 
  {
    int k = OCCQueryEngine::instance()->OCCMap->Find(shape);
    OCCQueryEngine::instance()->Shape_Label_Map->erase(k);
  }
*/
}

void OCCAttribSet::remove_attribute( const CubitSimpleAttrib& csa,
                                     TopoDS_Shape& shape)
{
  CubitBoolean found = CUBIT_FALSE;
  TDF_Label myLabel;

  FindShape(shape, myLabel, found);

  if(!found)
    return;

  int num_children = 0;
  int num_found = 0;
  for (TDF_ChildIterator it1(myLabel, CUBIT_FALSE); it1.More(); it1.Next())
  {
    //find the same type attribute first
    TDF_Label child = it1.Value();
    if (child.NbAttributes() == 0)
      continue;

    num_children ++;

    if(csa.isEmpty())
    {
      //only name attributes can't be forget here.
      Handle_TDataStd_Name attr_name;
      if(child.FindAttribute(TDataStd_Name::GetID(), attr_name))
      {
        TCollection_ExtendedString old_string = attr_name->Get();
        char const* string1 = "ENTITY_NAME";
        Standard_Boolean isMultiByte = Standard_True;
        TCollection_ExtendedString cstring( string1, isMultiByte );
        if( old_string == cstring)
          continue;
      }
      child.ForgetAllAttributes(false );
      num_found ++;
    }

    else
    {
      int found = find_attribute(child, csa);

      if(found == 3)
      {
        child.ForgetAllAttributes(false );
        num_found ++;
      } 
    }
  }
/*
  if(num_children == num_found)
  {
    if(OCCQueryEngine::instance()->OCCMap->IsBound(shape))
    {
      int k = OCCQueryEngine::instance()->OCCMap->Find(shape);
      OCCQueryEngine::instance()->Shape_Label_Map->erase(k);
    }
  }
*/
}

void OCCAttribSet::get_attributes(TDF_Label &lab,
                                  DLIList<CubitSimpleAttrib>& list)
{
  std::vector<CubitString> strings;
  CubitString string;
  std::vector<double> doubles;
  std::vector<int> ints;

  Handle_TDataStd_Name attr_name;
  TCollection_ExtendedString name_string;
  Handle_TDataStd_ExtStringArray attr_string;
  Handle_TDataStd_RealArray attr_doubles;
  Handle_TDataStd_IntegerArray attr_ints;

  if(lab.FindAttribute(TDataStd_Name::GetID(), attr_name))
  {
    name_string = attr_name->Get();
    int length = name_string.Length();
    std::vector<char> temp_type(length+1);
    for (int i = 1 ;  i <= length; i++)
    {
      Standard_ExtCharacter c = name_string.Value(i);
      temp_type[i-1] = ToCharacter(c);
    }
    temp_type[length] = '\0';
    string = CubitString(&temp_type[0]);
    strings.push_back(string);
  }
  if(lab.FindAttribute(TDataStd_ExtStringArray::GetID(), attr_string))
  {
    int length = attr_string->Length();
    for (int i = 0 ;  i < length; i++)
    {
      name_string = attr_string->Value(i);
      int length2 = name_string.Length();
      if (length2 == 0)
        continue;
      std::vector<char> temp_string(length2+1);
      for(int j = 1; j <= length2; j ++)
      {
        Standard_ExtCharacter c = name_string.Value(j);
        temp_string[j-1] = ToCharacter(c);
      }
      temp_string[length2] = '\0';
      const char *s = &temp_string[0];
      string = CubitString(s);
      strings.push_back(string);
    }
  }
  if(lab.FindAttribute(TDataStd_RealArray::GetID(), attr_doubles))
  {
    for(int i = 0; i < attr_doubles->Length(); i++)
    {
      double d = attr_doubles->Value(i);
      doubles.push_back(d);
    }
  }
  if(lab.FindAttribute(TDataStd_IntegerArray::GetID(), attr_ints))
  {
    for(int i = 0; i < attr_ints->Length(); i++)
    {
      int j = attr_ints->Value(i);
      ints.push_back(j);
    }
  }
  if(strings.size() > 0 || doubles.size() || ints.size())
  {
    list.append(CubitSimpleAttrib(&strings, &doubles, &ints));
  }
}

CubitStatus OCCAttribSet::get_attributes( TopoDS_Shape& shape,
                                          DLIList<CubitSimpleAttrib>& list)
{
  CubitBoolean found = CUBIT_FALSE;
  TDF_Label myLabel;

  FindShape(shape, myLabel, found);

  if(!found)
    return CUBIT_FAILURE;

  TDF_Label lab;
  for (TDF_ChildIterator it(myLabel,CUBIT_FALSE); it.More(); it.Next())
  {
    lab = it.Value();
    get_attributes(lab, list);
  }     
  return CUBIT_SUCCESS;
}

CubitStatus OCCAttribSet::get_attributes( const CubitString& name,
                                          TopoDS_Shape& shape,
                                    DLIList<CubitSimpleAttrib>& list )
{
  CubitBoolean found = CUBIT_FALSE;
  TDF_Label aLabel;

  FindShape(shape, aLabel, found);

  if(!found)
    return CUBIT_FAILURE;

  Handle_TDataStd_Name attr_name;
  TCollection_ExtendedString cstring( (Standard_CString)name.c_str() );
  TDF_Label lab;
  for (TDF_ChildIterator it(aLabel,CUBIT_FALSE); it.More(); it.Next())
  {
    lab = it.Value();

    //check if the type attribute is the same, disable it
    TCollection_ExtendedString old_string;
    if(lab.FindAttribute(TDataStd_Name::GetID(), attr_name))
    {
      old_string = attr_name->Get();

      if(old_string == cstring)
      {
        get_attributes(lab, list);
        break;
      }
    }
  }
 
  return CUBIT_SUCCESS;
}


int OCCAttribSet::attribute_count() 
{
  int count = 0;
  
  return count;
}

