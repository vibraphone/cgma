#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>

#include <vector>

#include <iGeom.h>

//#include "SimpleArray.h"
template<typename T>
class SimpleArray 
{
public:
  T *array;
  int array_alloc;
  int array_size;
  
  SimpleArray() : array(NULL), array_alloc(0), array_size(0){}
  ~SimpleArray() {reset();}

  T operator[](unsigned int lhs) {return array[lhs];}
        
  void reset() {
    if (0 != array) {
      free(array); array = NULL; array_alloc = array_size = 0;
    }
  }
  int size() {return array_size;}
  
};
#define SA_ARRAY_INOUT(a) &a.array, &a.array_alloc, &a.array_size
#define SA_ARRAY_IN(a) a.array, a.array_size

using namespace std;

int main( int argc, char **argv)
{
  assert( argc == 2 );
  string engine_opt = ";engine=OCC";
  string filename   = argv[1]; 

  int err;
  iGeom_Instance geom;
  iGeom_newGeom( engine_opt.c_str(), &geom, &err, engine_opt.length() );

  iGeom_load( geom, &filename[0], 0, &err, filename.length(), 0 );

  iBase_EntitySetHandle root_set;
  iGeom_getRootSet( geom, &root_set, &err);

  cout << "Model Contents " << endl;
  const char *gtype[] = {"vertices: ", "edges: ", "faces: ", "regions: "};
  for (int i = 0; i <= 3; ++i) {
    int count;
    iGeom_getNumOfType( geom, root_set, i, &count, &err );
    std::cout << gtype[i] << count << std::endl;
  }

  iBase_TagHandle idtag;
  iGeom_getTagHandle( geom, "GLOBAL_ID", &idtag, &err, strlen("GLOBAL_ID") );

  SimpleArray<iBase_EntityHandle> vertices;
  iGeom_getEntities( geom, root_set, iBase_VERTEX, SA_ARRAY_INOUT(vertices), &err);

  SimpleArray<double> vCoords;
  iGeom_getVtxArrCoords(geom, SA_ARRAY_IN(vertices), iBase_INTERLEAVED, SA_ARRAY_INOUT(vCoords), &err);

  SimpleArray<iBase_EntityHandle> edges;
  iGeom_getEntities( geom, root_set, iBase_EDGE, SA_ARRAY_INOUT(edges), &err);

  SimpleArray<double> edgelength;
  iGeom_measure( geom, SA_ARRAY_IN(edges), SA_ARRAY_INOUT(edgelength), &err);

  double minlength = edgelength[0];
  for( int i = 0; i < edges.size(); i++) 
    minlength = std::min(minlength, edgelength[i]);

  double ustart, uend, delta;
  double x, y, z;
  int N = 10;

  int id;
  int numNodes = 0, numEdges = 0;
  vector<double> xCoord, yCoord, zCoord;
  for( int i = 0; i < edges.size(); i++) {
    iGeom_getEntURange( geom, edges[i], &ustart, &uend, &err);
    iGeom_getIntData( geom, edges[i], idtag, &id, &err);
    delta = (uend-ustart)/(double)N;
    for( int j = 0; j < N+1; j++) {
      double u = ustart + j*delta;
      iGeom_getEntUtoXYZ( geom, edges[i], u, &x, &y, &z, &err);
      xCoord.push_back(x);
      yCoord.push_back(y);
      zCoord.push_back(z);
      numNodes++;
    }
    numEdges += N;
  }

  SimpleArray<iBase_EntityHandle> faces;
  iGeom_getEntities( geom, root_set, iBase_FACE, SA_ARRAY_INOUT(faces), &err);

  int eid, vid1, vid2;
  SimpleArray<iBase_EntityHandle>  edgenodes;
  SimpleArray<iBase_EntityHandle>  faceedges;
  SimpleArray<iBase_EntityHandle>  facenodes;

  for( int i = 0; i < faces.size(); i++) 
  {
    iGeom_getIntData( geom, faces[i], idtag, &id, &err);
    iGeom_getEntAdj( geom, faces[i], iBase_EDGE,   SA_ARRAY_INOUT( faceedges ), &err);
    cout << "Face " << id << ": " << faceedges.size() << " edges." << endl;
    for(int j = 0; j < faceedges.size(); ++j) 
    {
      iGeom_getIntData( geom, faceedges[j], idtag, &eid, &err);
      iGeom_getEntAdj( geom, faceedges[j], iBase_VERTEX, SA_ARRAY_INOUT(edgenodes), &err);
      cout << "Edge " << eid << ", " << edgenodes.size() << " vertices: ";
      assert( edgenodes.size() == 2 );
      iGeom_getIntData( geom, edgenodes[0], idtag, &vid1, &err);
      iGeom_getIntData( geom, edgenodes[1], idtag, &vid2, &err);
      cout << vid1 << ", " << vid2 << endl;
    }
  }

}
