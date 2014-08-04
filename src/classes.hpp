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
  float minC;
  float maxC;
  float minL;
  float maxL;
  float minR;
  float maxR;
  float centerX;
  float centerY;
  float centerZ;

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
  float *smoothArr;
  int    numGLL=0;
  int    NX;
  int    NY;
  int    NZ;

  // These are the functions.
  void readNetcdf          ( std::string mode, std::string fname );
  void createKDtree        ( );
  void mergeKernels        ( std::vector<Kernel> &allKern );
  void createRegMesh       ( std::vector<Kernel> &allKern );
  void quickSortCenter     ( int i1st, int i2nd );
  void quickSortPoint      ( int i1st, int i2nd, float px, float py, float pz );
  void getMinMaxCartesian  ( );
  void writeExodus         ( );
  void rotateXaxis         ( double &deg );
  void rotateYaxis         ( double &deg );
  void rotateZaxis         ( double &deg );
  double angleFromAxis      ( char axis );


  void smoothGaussian      ( float &var );
  float distFromCenter     ( float &x, float &y, float &z );
  float distFromPoint      ( float &x,  float &y,  float &z, 
                             float &px, float &py, float &pz );

  int  pivot               ( int   &i1st, int &i2nd,   
                             float &d1st, float &d2nd );
};

void createRegMesh ( Kernel &kern, std::vector<Kernel> &allKern );
void exodusErrorCheck ( int ier, std::string function );
