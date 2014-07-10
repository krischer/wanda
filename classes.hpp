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
  float *param;

  // This is the internal regular mesh which we interpolate onto.
  void *regMesh;

  // This array holds the data (indices) for the kdtree.
  int *KDdat;

  // These are the functions.
  void interpolate   ( float &x, float &y, float &z );
  void createKDtree  ( float &x, float &y, float &z );
  void createRegMesh ( );

};
