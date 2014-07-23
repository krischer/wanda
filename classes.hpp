#include <iostream>
#include <vector>
#include "kdtree.h"

class Kernel;

class Kernel
{

public:


  // These are used for the external kernel mesh.
  float *xExt;
  float *yExt;
  float *zExt;

  // These values hold the min/max points for kernal mesh.
  float minX;
  float minY;
  float minZ;
  float maxX;
  float maxY;
  float maxZ;

  // These values hold the min/max for the entire regMesh.
  float minXBox;
  float maxXBox;
  float minYBox;
  float maxYBox;
  float minZBox;
  float maxZBox;

  // This is the internal regular mesh which we interpolate onto.
  float *regMeshArr;
  float *regX;
  float *regY;
  float *regZ;

  // This array holds the data (indices) for the kdtree.
  int    *KDdat;
  kdtree *tree;

  // This array holds the gll kernel values from one processor.
  float *kernStore;
  int    numGLL = 0;

  // These are the functions.
  void readNetcdf          ( std::string mode, std::string fname );
  void createKDtree        ( std::vector<Kernel> &allKern );
  void mergeKernels        ( std::vector<Kernel> &allKern );
  void createRegMesh       ( std::vector<Kernel> &allKern );
  void getMinMaxCartesian  ( );

};

void createRegMesh ( Kernel &kern, std::vector<Kernel> &allKern );
