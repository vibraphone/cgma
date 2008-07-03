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
#include "OCCAttrib.hpp"
#include "CubitSimpleAttrib.hpp"
#include "CubitFileIOWrapper.hpp"
#include "Handle_TDataStd_Shape.hxx"
#include "TCollection_ExtendedString.hxx"
#include "Handle_TDataStd_Name.hxx"
#include "Handle_TDataStd_ExtStringArray.hxx"
#include "Handle_TDataStd_RealArray.hxx"
#include "Handle_TDataStd_IntegerArray.hxx"
#include "TDataStd_Name.hxx"
#include "TDataStd_ExtStringArray.hxx"
#include "TDataStd_RealArray.hxx"
#include "TDataStd_IntegerArray.hxx"
#include "TDataStd_Shape.hxx"
#include "TopoDS_Shape.hxx"
#include "TDF_ChildIterator.hxx"

void OCCAttribSet::append_attribute( CubitSimpleAttrib* csa, TopoDS_Shape& shape )
{
  OCCAttrib* new_attrib = new OCCAttrib(csa);
  new_attrib->listNext = listHead;
  listHead = new_attrib;

  TDF_Label child;
  //Add attributes on child of myLabel
  //1. add shape attribute, first check to make sure there's no shape attribute
  CubitBoolean found = CUBIT_FALSE;
  for (TDF_ChildIterator it1(myLabel, CUBIT_FALSE); it1.More(); it1.Next())
  {
    //find the same type attribute first
    child = it1.Value();

    Handle_TDataStd_Shape attr_shape;
    TopoDS_Shape old_shape;
    if(TDataStd_Shape::Find(child, attr_shape))
      old_shape = attr_shape->Get(child);
    if(old_shape.IsEqual(shape)) 
    {
      found = CUBIT_TRUE;
      break;
    }
  }

  if(!found)
  { 
    child = myLabel.NewChild();
    Handle_TDataStd_Shape attr_shape = TDataStd_Shape::Set(child, shape);
    child.AddAttribute(attr_shape);
  }

  //2. add type attribute , belowing attributes are added on child lable of child
  //First create new child label
  TDF_Label lab = child.NewChild();
  
  //Then add child labels
  CubitString type = csa->character_type(); 
  if(type.length() > 0)
  {
    TCollection_ExtendedString cstring( (Standard_CString)type.c_str() );
    Handle_TDataStd_Name attr_name = TDataStd_Name::Set(lab, cstring);
    lab.AddAttribute(attr_name);
  }

  //3. add string attribute
  DLIList<CubitString*>* strings = csa->string_data_list();
  TDF_Label lab_child = lab.NewChild(); 
  Handle_TDataStd_ExtStringArray attr_string;
  if(strings && strings->size() > 0)
  {
    //set the length of String Array be 11, can be extended.
    attr_string = TDataStd_ExtStringArray::Set(lab_child, 0, 10);
  }

  for(int i = 0; strings && i < strings->size(); i++)
  {
    TCollection_ExtendedString 
       cstring((Standard_CString)strings->get_and_step()->c_str() );
    attr_string->SetValue(i, cstring) ;
  }
  if(strings && strings->size() > 0)
    lab_child.AddAttribute(attr_string);

  //4. add double attribute
  DLIList<double*>* doubles = csa->double_data_list();
  Handle_TDataStd_RealArray attr_double;
  if(doubles && doubles->size() > 0)
  {
    //set the length of double array be 11, can be extended.
    attr_double = TDataStd_RealArray::Set(lab_child,0, 10);
  } 

  for(int i = 0; doubles && i < doubles->size(); i++)
    attr_double->SetValue(i, *(doubles->get_and_step()));

  if(doubles && doubles->size() > 0)
    lab_child.AddAttribute(attr_double);
    
  //5. add int attribute
  DLIList<int*>* ints = csa->int_data_list();
  Handle_TDataStd_IntegerArray attr_int;
  if(ints && ints->size() > 0)
  {
    //set the length of int array be 11, can be extended.
    attr_int = TDataStd_IntegerArray::Set(lab_child, 0, 10);
  }

  for(int i = 0; ints && i < ints->size(); i++)
    attr_int->SetValue(i, *(ints->get_and_step()));

  if(ints && ints->size() > 0)
    lab_child.AddAttribute(attr_int);
}

void OCCAttribSet::remove_attribute( CubitSimpleAttrib* csa)
{
  if( !listHead )
    return;
    
  OCCAttrib* attrib = 0;
  if ( listHead->equals(csa) )
  {
    attrib = listHead;
    listHead = listHead->listNext;
    delete attrib;
    return;
  }
  
  for ( OCCAttrib* prev = listHead; prev->listNext; prev = prev->listNext )
  {
    if( prev->listNext->equals(csa) )
    {
      attrib = prev->listNext;
      prev->listNext = attrib->listNext;
      delete attrib;
      return;
    }
  }

  //forget csa attribute from the document
  DLIList<int*>* ints = csa->int_data_list();
  CubitString type = csa->character_type();
  TCollection_ExtendedString cstring( (Standard_CString)type.c_str() );

  if(type.length() == 0)
      return;

  DLIList<double*>* doubles = csa->double_data_list();
  DLIList<CubitString*>* strings = csa->string_data_list();  
  for (TDF_ChildIterator it1(myLabel, CUBIT_TRUE); it1.More(); it1.Next())
  {
    //find the same type attribute first
    TDF_Label child = it1.Value();

    Handle_TDataStd_Name attr_name;
    TCollection_ExtendedString old_string;
    if(child.FindAttribute(TDataStd_Name::GetID(), attr_name))
      old_string = attr_name->Get(); 

    if(old_string != cstring) 
      continue; 
  
    //continue to compare the rest attributes.
    CubitBoolean is_same = CUBIT_TRUE;
    for(TDF_ChildIterator it2(child,CUBIT_FALSE); it2.More(); it2.Next())
    {
      is_same = CUBIT_TRUE;
      TDF_Label g_child = it2.Value();
      if(ints->size() > 0 )
      {
        Handle_TDataStd_IntegerArray attr_ints;
        if(g_child.FindAttribute(TDataStd_IntegerArray::GetID(), attr_ints) && 
           attr_ints->Length() == ints->size())
        {
          for(int i = 0; i < ints->size(); i++)
            if(attr_ints->Value(i) != *ints->get_and_step())   
            {
              is_same = CUBIT_FALSE;
              break;
            } 
        }
      }
      if(!is_same)
        continue;
 
      if(doubles->size() > 0 )
      {
        Handle_TDataStd_RealArray attr_doubles;
        if(g_child.FindAttribute(TDataStd_RealArray::GetID(), attr_doubles) &&
           attr_doubles->Length() == doubles->size())
        {
          for(int i = 0; i < doubles->size(); i++)
            if(attr_doubles->Value(i) != *doubles->get_and_step())
            {
              is_same = CUBIT_FALSE;
              break;
            }
        }
      }

      if(!is_same)
        continue;

      if(strings->size() > 0 )
      {
        Handle_TDataStd_ExtStringArray attr_strings;
        if(g_child.FindAttribute(TDataStd_ExtStringArray::GetID(), attr_strings) &&
           attr_strings->Length() == strings->size())
        {
          for(int i = 0; i < strings->size(); i++)
          {
            CubitString astring = *strings->get_and_step();
            TCollection_ExtendedString string( (Standard_CString)astring.c_str() );
            if(attr_strings->Value(i) != string)
            {
              is_same = CUBIT_FALSE;
              break;
            }
          }
        }
      }
      if(!is_same)
        continue;

      child.ForgetAllAttributes( );
    }
  }
}

void OCCAttribSet::remove_all_attributes()
{
  while( listHead )
  {
    OCCAttrib* dead = listHead;
    listHead = dead->listNext;
    delete dead;
  }
  myLabel.ForgetAllAttributes( );
}

CubitStatus OCCAttribSet::get_attributes( DLIList<CubitSimpleAttrib*>& list ) const
{
  for( OCCAttrib* attrib = listHead; attrib; attrib = attrib->listNext )
    list.append( attrib->get_CSA() );
  return CUBIT_SUCCESS;
}

CubitStatus OCCAttribSet::get_attributes( const CubitString& name,
                                    DLIList<CubitSimpleAttrib*>& list ) const
{
  for( OCCAttrib* attrib = listHead; attrib; attrib = attrib->listNext )
    if( attrib->name() == name )
      list.append( attrib->get_CSA() );
  return CUBIT_SUCCESS;
}

CubitStatus OCCAttribSet::save_attributes( FILE* file_ptr ) const
{
  OCCAttrib *curr_attrib;
  CubitStatus status = CUBIT_SUCCESS;
  
  //save # attribs
  unsigned int size = attribute_count();
  NCubitFile::CIOWrapper wrapper( file_ptr );
  wrapper.Write( &size, 1 ); 

  //save each attrib
  for( curr_attrib = listHead; curr_attrib; curr_attrib = curr_attrib->listNext )
    if( !curr_attrib->save(file_ptr) )
      status = CUBIT_FAILURE;

  return status;
}
  
CubitStatus OCCAttribSet::restore_attributes( FILE* file_ptr, unsigned endian )
{
  OCCAttrib *curr_attrib;
  
  //Read # attribs
  unsigned int size;
  NCubitFile::CIOWrapper wrapper( endian, file_ptr );
  wrapper.Read( &size, 1 ); 

  for (unsigned i = 0; i < size; i++)
  {
    curr_attrib = OCCAttrib::restore( file_ptr, endian);  
    if (!curr_attrib)
    {
        // file corrupt?  don't try to read any more
      return CUBIT_FAILURE;
    }
    
    curr_attrib->listNext = listHead;
    listHead = curr_attrib;
  }

  return CUBIT_SUCCESS;
}


int OCCAttribSet::attribute_count() const
{
  int count = 0;
  for( OCCAttrib* attrib = listHead; attrib; attrib = attrib->listNext )
    count++;
  return count;
}

