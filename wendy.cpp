#include <cstdlib>
#include <netcdf>
#include "classes.hpp"

void returnRegularArray ( int NUM_X, int NUM_Y, int NUM_Z, void *testArr )
{

  Kernel kern;
  std::cout << "HELLO WORLD." << std::flush << std::endl;

  std::cout << NUM_X << ' ' << NUM_Y << ' ' << NUM_Z << ' ' << testArr << std::endl;

//  kern.regMesh = std::malloc ( 10 * 10 * 10 * sizeof(float) );
  float (*dmat)[NUM_X][NUM_Y][NUM_Z] = ( (float (*)[NUM_X][NUM_Y][NUM_Z]) testArr );
  for ( int i=0; i<NUM_X; i++ ) {
    for ( int j=0; j<NUM_Y; j++ ) {
      for ( int k=0; k<NUM_Z; k++ ) {

        *dmat[i][j][k] = 10.;

      }
    }
  }

}

void Kernel::createKDtree ( float &x, float &y, float &z )
{

  std::cout << "Creating KDTree." << std::flush << std::endl;

  // Initialize the KDtree.
  kdtree *tree;
  tree = kd_create (3);

  // Initizalize the data array.
  KDdat = new int [numGLL];

  // Populate the tree with an index (KDdat) and x, y, z.
  for ( size_t i=0; i<numGLL; i++ )
  {
    KDdat[i] = i;
    kd_insert3 ( tree, xExt[i], yExt[i], zExt[i], &KDdat[i] );

  }

}

void Kernel::readNetcdf ( std::string mode, std::string fname )
{

  using namespace netCDF;
  using namespace netCDF::exceptions;

  try
  {

    NcFile dataFile ( fname, NcFile::read );

    if ( mode == "kernel" )
    {

      NcVar gllValues = dataFile.getVar  ("param");
      NcDim dim       = gllValues.getDim (0);
      numGLL          = dim.getSize      ();

      kernStore = new float [numGLL];
    
      gllValues.getVar (kernStore);

    }
    else if ( mode == "coordinates" )
    {

      NcVar dataX = dataFile.getVar ("dataX");
      NcVar dataY = dataFile.getVar ("dataY");
      NcVar dataZ = dataFile.getVar ("dataZ");
      NcDim dim   = dataX.getDim    (0);
      numGLL      = dim.getSize     ();

      xExt = new float [numGLL];
      yExt = new float [numGLL];
      zExt = new float [numGLL];

      dataX.getVar (xExt);
      dataY.getVar (yExt);
      dataZ.getVar (zExt);

    }
    
  } 
  catch (NcException &e)
  {

    e.what();
    std::cout << "Failure reading in kernel." << std::endl;
    std::exit ( EXIT_FAILURE );

  }
}
