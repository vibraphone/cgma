#include "iGeom.h"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <assert.h>

// helper macro for igeom
#define ASSERT(A) if (!(A)) failed(#A,__FILE__,__LINE__)
#define CHECK( STR ) if (err != iBase_SUCCESS) return Print_Error( STR, err, geom, __FILE__, __LINE__ )
int Tag_Body( iGeom_Instance geom, const iBase_TagHandle tag, const std::string id, const iBase_EntityHandle body);
int Tag_Get( iGeom_Instance geom, const iBase_TagHandle tag, const iBase_EntityHandle body);
bool Print_Error( const char* desc, 
		  int err,
		  iGeom_Instance geom,
		  const char* file,
		  int line );
int main(){
  int err;
  iGeom_Instance geom;
  double dSide = 10.0, dHeight = 10.0, dRad = 3.0;
  iBase_EntityHandle assm = NULL, cyl = NULL, tmp_cyl = NULL, tmp_new = NULL;
  std::string geomFile = "test.occ";


  iGeom_newGeom( 0, &geom, &err, 0 );

  iBase_TagHandle pin_tag = NULL, name_tag = NULL;
  char* tag_name_pin = (char*)"PIN";
  char* tag_name = (char*)"NAME";
  iGeom_getTagHandle(geom, tag_name_pin, &pin_tag, &err, 3);
  if(err == iBase_TAG_NOT_FOUND){
      
    iGeom_createTag(geom, tag_name_pin, 1, iBase_INTEGER, &pin_tag, &err,  3);

  }

  iGeom_getTagHandle(geom, tag_name, &name_tag, &err, 4);
  if(err == iBase_TAG_NOT_FOUND){
      
    iGeom_createTag(geom, tag_name, 1, iBase_INTEGER, &name_tag, &err,  4);

  }

  // creating prism
  std::cout << "\n\n\nCreating Prism\n" << std::endl;
  iGeom_createPrism(geom, dHeight, 6, 
		    dSide, dSide,
		    &assm, &err); 
  std::string one = "1";
  Tag_Body(geom, pin_tag, one, assm);
  Tag_Get (geom, pin_tag, assm);
  std::cout << "name tag: ";
  std::string a = "a";
  Tag_Body(geom, name_tag, a, assm);
  Tag_Get (geom, name_tag, assm);  

  std::cout << "\n\n\nCreating Cylinder\n" << std::endl;
  // create cylinder
  iGeom_createCylinder(geom, dHeight, dRad, dRad, &cyl, &err);
  std::string two = "2";
  Tag_Body(geom, pin_tag, two, cyl);
  Tag_Get(geom, pin_tag, cyl);

  std::cout << "name tag: ";
  std::string b = "b";
  Tag_Body(geom, name_tag, b, cyl);
  Tag_Get (geom, name_tag, cyl); 

  // Copy
  iGeom_copyEnt(geom, cyl, &tmp_cyl, &err);
  std::cout << "\n\n After copy operation\n" << std::endl;
  Tag_Get(geom, pin_tag, tmp_cyl);
  std::cout << "\ngetting name tag " << std::endl;
  Tag_Get(geom, name_tag, tmp_cyl);

  // Substract
  iGeom_subtractEnts(geom, assm, tmp_cyl, &tmp_new, &err);
  std::cout << "\n\n After subtract operation\n" << std::endl;
  Tag_Get(geom, pin_tag, tmp_new);
  std::cout << "\n getting name tag " << std::endl;
  Tag_Get(geom, name_tag, tmp_new);

  // save
  iGeom_save(geom, geomFile.c_str(), NULL, &err, geomFile.length() , 0);
  
  //check that the two single volume bodys' attributes are exported as SINGLELUMP%
  std::string search = "SINGLELUMP%";
  std::ifstream Myfile;
  Myfile.open (geomFile.c_str());
  int found = 0;
  std::string line;
  size_t offset;
  if(Myfile.is_open())
  {
    while(!Myfile.eof())
    {
      getline(Myfile,line);
      if ((offset = line.find(search, 0)) != std::string::npos)
        found ++;
    }
    Myfile.close();
  }

  assert (found == 4);

  return 0;
}
  

int Tag_Body( iGeom_Instance geom, const iBase_TagHandle tag, const std::string id, const iBase_EntityHandle body)
//---------------------------------------------------------------------------
//Function: Tag body with pin number
//Input:    none
//Output:   none
//---------------------------------------------------------------------------
{
  int err = 0;

  iGeom_setData(geom, body, tag, id.c_str(), id.size(), &err);
  std::cout<<"\nset pin tag - " << id<< " on " << body << std::endl;
  return 0;
}



int Tag_Get(iGeom_Instance geom, const iBase_TagHandle tag, const iBase_EntityHandle body)
//---------------------------------------------------------------------------
//Function: Tag get
//Input:    none
//Output:   none
//---------------------------------------------------------------------------
{
  int err = 0, bytes;
  iGeom_getTagSizeBytes( geom, tag, &bytes, &err );
  CHECK( "Can't get size" );

  std::vector<char> data(bytes);

  //just check if pin tag exist, IT'S NOT USED IN THIS FUNCTION
  iBase_TagHandle pin1 = NULL;
  char* p = (char*)"PIN";

  iGeom_getTagHandle(geom, p, &pin1, &err, 3);
  CHECK( "PIN tag doesn't exist" );

  void* ptr = &data[0];
  int junk1 = bytes, junk2;
  iGeom_getData( geom, body, tag, (void**)&ptr, &junk1, &junk2, &err );
  CHECK( "failed to getData for tag" );
  std::cout << "Able to get this tag: "<<*(void**)ptr << " on " << body <<std::endl;
  return 0;
}



// print error function definition (iGeom)
bool Print_Error( const char* desc, 
		  int err,
		  iGeom_Instance geom,
		  const char* file,
		  int line )
{
  char buffer[1024];
  iGeom_getDescription( geom, buffer, sizeof(buffer) );
  buffer[sizeof(buffer)-1] = '\0';
  
  std::cerr << "ERROR: " << desc << std::endl
	    << "  Error code: " << err << std::endl
	    << "  Error desc: " << buffer << std::endl
	    << "  At        : " << file << ':' << line << std::endl
    ;
  
  return false; // must always return false or CHECK macro will break
}
