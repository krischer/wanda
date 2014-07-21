#include "classes.hpp"
#include "kdtree.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>

int main ()
{

  Kernel kern;

  std::vector <Kernel> allKern;

  std::cout << "Hello world." << std::endl;
  std::cout << "How many processors." << std::endl;

  int nProc;
  std::cin >> nProc;

  allKern.reserve (nProc);

  for ( int i=0; i<nProc; i++ )
  {

    Kernel kern;

    /*
     * This part could be replaced by food from python 
     */
    std::string fNameBase = "./testKernels/proc";
    std::string fNameApp  = "_reg1_betah_kernel.nc";

    std::string xyzNameBase = "./cemRequest/xyz_reg01_proc";

    std::stringstream ssPrc;
    ssPrc << std::setw (6) << std::setfill ('0');
    ssPrc << std::to_string (i);

    std::stringstream ssPrcShort;
    ssPrcShort << std::setw (4) << std::setfill ('0');
    ssPrcShort << std::to_string (i);

    fNameBase.append ( ssPrc.str() );
    fNameBase.append ( fNameApp );

    xyzNameBase.append ( ssPrcShort.str() );

    kern.readNetcdf ( "param", fNameBase );
    kern.readNetcdf ( "coordinates", xyzNameBase );
    /*
     * This part could be replaced by food from python 
     */

    kern.getMinMaxCartesian ();
    allKern.push_back       (kern);

  }

  kern.createKDtree ( allKern );
  createRegMesh     ( kern, allKern );

  return 0;

}

void createRegMesh ( Kernel &kern, std::vector<Kernel> &allKern )
{

  float DX = 10.;
  float DY = 10.;
  float DZ = 1.;

  // Initialize box values.
  kern.minXBox = allKern[0].minX;
  kern.maxXBox = allKern[0].maxX;
  kern.minYBox = allKern[0].minY;
  kern.maxYBox = allKern[0].maxY;
  kern.minZBox = allKern[0].minZ;
  kern.maxZBox = allKern[0].maxZ;

  // Find the max/min of the master box.
  for ( int i=0; i<allKern.size(); i++ )
  {

    if ( kern.minXBox > (allKern[i].minX) )
      kern.minXBox = allKern[i].minX;
    if ( kern.maxXBox < (allKern[i].maxX) )
      kern.maxXBox = allKern[i].maxX;

    if ( kern.minYBox > (allKern[i].minY) )
      kern.minYBox = allKern[i].minY;
    if ( kern.maxYBox < (allKern[i].maxY) )
      kern.maxYBox = allKern[i].maxY;

    if ( kern.minZBox > (allKern[i].minZ) )
      kern.minZBox = allKern[i].minZ;
    if ( kern.maxZBox < (allKern[i].maxZ) )
      kern.maxZBox = allKern[i].maxZ;

  }

  // Determine # of grid points in each direction, given DX.
  long NX = int ( ((kern.maxXBox - kern.minXBox) / DX) + 1);
  long NY = int ( ((kern.maxYBox - kern.minYBox) / DY) + 1);
  long NZ = int ( ((kern.maxZBox - kern.minZBox) / DZ) + 1);


  kern.regMeshArr = new float [NY*NY*NZ];
  float (*regMesh)[NX][NY][NZ] = ( (float (*)[NX][NY][NZ]) kern.regMeshArr );

  float x = 0;
  float y = 0;
  float z = 0;

//##pragma omp parallel for 
  for ( int i=0; i<NX; i++ ) {
    for ( int j=0; j<NY; j++ ) {
      for ( int k=0; k<NZ; k++ ) {

        x = kern.minXBox + i * DX;
        y = kern.minYBox + i * DY;
        z = kern.minZBox + i * DZ;

//        regMesh[i][j][k] = kd


      }
    }
  }

}




  




