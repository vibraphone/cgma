/*




*/

//-------------------------------------------------------------
//
// File:  main.cpp
// 
// Description:
//  This is the big enchilada
//      
// Author: Bob Kerr
//
//-------------------------------------------------------------


#include "my.hpp"
#include "gm_cgm_c_interface.h"
#include "CubitDefines.h"
#include "DLIList.hpp"
#include "RefEntity.hpp"
#include "GeometryQueryTool.cpp"

int fill_vertex(char* name, STD_IO(ifstream) &fin);

int main(int argc, char **argv)
{
  int i;
  int flag = 0;
  char filename[21];
  for(i = 0; i<argc; i++)
  {
    printf("Argument %d was %s\n",i,argv[i]);
  }
  
  
    //user wants to open a file
  if(argc < 2)
  {
    printf("Hey dummy, you forgot an input file\n");
    return 1;
  }
  strncpy(filename, argv[1], 20);
 
    //set file pointer again
  STD_IO(ifstream) fin(argv[1]);
  if(!fin.good())
  {
      //check to see if file exists
    STD_IO(cout)<<"The file named "<<filename<<" does not exist!\n";
    flag = 0;
  }
  else
     flag = 1;
  if(cpp_cgm_initialize() != 0)
  {
    STD_IO(cout)<<"Error, cgm didn't initialize"<<STD_IO(endl);
    exit(1);
  }
  
  
  char line[100];
  char* tokens1 = "= \t,\n";
  char* tag;
  char* type;
  char* currentline;
  int abortflag = 1;
 
  fin.getline(line,80);
    //read a line from the file until we're at the end
  do
  {
    currentline = line;

    type = strtok(currentline, tokens1);
    tag = strtok(NULL, tokens1);
      //Switch to whichever type of geometry we are given
    if(type)
    {
      int newtype = String_to_enum(type);
      int   numItems=0;
      switch(newtype)
      {
        case VERTEX:
           STD_IO(cout)<<"Read Vertex ";
           double items[3];
           while ( (type = strtok(NULL, tokens1)) != NULL )
           {
             if (!sscanf(type, "%lf", &items[numItems]))
                break;
             
             numItems++;
             
             currentline= NULL;
           }
             
           if (numItems == 3)
           {
             VertexHandle* vert_handle;
             
             STD_IO(cout)<<tag<<" with three coordinates "<<items[0]<<" "<<items[1]<<" "<<items[2]<<STD_IO(endl);
             if(cgm_vertex_construct(UNDEFINED_POINT_TYPE, tag, 3, items, &vert_handle) !=0)
                STD_IO(cout)<<"Error creating Vertex "<<tag<<STD_IO(endl);
           }
           else
              STD_IO(cout)<<" but didn't have three coordinates\n"<<STD_IO(endl);
           break;
        case EDGE:
           STD_IO(cout)<<"Read EDGE";
           type = strtok(NULL, tokens1);
           if(!strcmp(type, "STRAIGHT"))
           {
             char* vert1, *vert2;
             vert1 = strtok(NULL, tokens1);
             vert2 = strtok(NULL, tokens1);
             EdgeHandle* edge_handle;
             STD_IO(cout)<<" with ends "<<vert1<<" and "<<vert2<<STD_IO(endl);
             if(cgm_straight_edge_construct(
                   STRAIGHT_CURVE_TYPE, tag, vert1, vert2, &edge_handle)!=0)
                STD_IO(cout)<<"Error creating Curve "<<tag<<STD_IO(endl);
             
           }
           else if(!strcmp(type, "ELLIPSE"))
           {
             char* vert1, *vert2, *bulge;
             vert1 = strtok(NULL, tokens1);
             vert2 = strtok(NULL, tokens1);
             double coordinate[3];
             for(int j = 0; j < 3; j++)
                coordinate[j] = atof(strtok(NULL, tokens1));
             bulge = strtok(NULL, tokens1);
             EdgeHandle* edge_handle;
             STD_IO(cout)<<" elliptical, with bulge "<<bulge;
             STD_IO(cout)<<" with ends "<<vert1<<" and "<<vert2<<STD_IO(endl);
             CubitSense edge_sense = strcmp(bulge, "FORWARD") ?
                CUBIT_REVERSED : CUBIT_FORWARD;
             if(cgm_quadratic_edge_construct(ELLIPSE_CURVE_TYPE, tag, vert1, vert2,
                                             3, coordinate, edge_sense, &edge_handle)!=0)
                STD_IO(cout)<<"Error creating Curve "<<tag<<STD_IO(endl);
           }
           else if(!strcmp(type, "PARABOLA"))
           {
           }
           else if(!strcmp(type, "COMPOSITE"))
           {
           }
           

           break;
        case FACE:
           STD_IO(cout)<<"Read FACE ";
           type = strtok(NULL, tokens1);
           if(!strcmp(type, "DISK"))
           {
             double items[3];
             while ( (type = strtok(NULL, tokens1)) != NULL )
             {
               if (!sscanf(type, "%lf", &items[numItems]))
                  break;
               
               numItems++;
               currentline= NULL;
             }
             
             if (numItems == 3)
             {
               FaceHandle* face_handle;
               STD_IO(cout)<<tag<< " type DISK with x = "<<items[0]<<" y = "<<items[1]
                   <<" radius = "<<items[2]<<STD_IO(endl);
               if(cgm_face_disk_construct(PLANE_SURFACE_TYPE ,tag,items,
                                          items[2], &face_handle)!= 0)
                  STD_IO(cout)<<"Error creating face "<<tag<<STD_IO(endl);
             }
             else
                STD_IO(cout)<<" but didn't have the right inputs\n";
           }
           else
              STD_IO(cout)<<"Only DISK is available"<<STD_IO(endl);
           break;
        default:
           STD_IO(cout)<<"\nError, unable to determine element type\n";
           break;
      }
    }
    if(!abortflag)
    {
        //There were errors so we'll quit reading.
      STD_IO(cout)<<"\nERROR--can't parse input file, aborting. . .";
      fin.close();
      exit(1);
    }

      //Check for end of file
    if(!fin.getline(line,80))
       flag = FALSE;
  
  }
 
  while(flag);

  DLIList<RefEntity*> entity_list;
  const char* file_type = "ACIS_SAT";
  CubitString version = "7.0.0";
  CubitString name = "temp.sat";
  int numexport = 0;
  int status = GeometryQueryTool::instance()->
     export_solid_model(entity_list,
                        name.c_str(),
                        file_type,
                        numexport,
                        version);
  
  
// }
 
}

