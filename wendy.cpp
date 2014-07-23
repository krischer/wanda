#include <cstdlib>
#include <netcdf>
#include "classes.hpp"
#include "mpi.h"

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

void Kernel::createKDtree ( std::vector<Kernel> &allKern )
{

  // Initialize the KDtree.
  tree = kd_create (3);

  // Initizalize the data array.
  KDdat = new int [numGLL];

  // Populate the tree with an index (KDdat) and x, y, z
  // We need to loop over both mesh chunks and gll points.
  if ( MPI::COMM_WORLD.Get_rank() == 0 )
    std::cout << "Creating KDtree." << std::flush << std::endl;
  for ( size_t i=0; i<numGLL; i++ )
  {
    KDdat[i] = i;
    kd_insert3 ( tree, xExt[i], yExt[i], zExt[i], &KDdat[i] );
  }

}

void Kernel::mergeKernels ( std::vector<Kernel> &allKern )
{

  // Determine total number of gll points.
  for ( size_t i=0; i<allKern.size(); i++ )
    numGLL += allKern[i].numGLL;

  // Initialize parameter array.
  kernStore = new float [numGLL];

  // Initialize coordinate arrays.
  xExt = new float [numGLL]; 
  yExt = new float [numGLL]; 
  zExt = new float [numGLL]; 

  // Copy to master parameter array.
  int totIter = 0;
  for ( size_t chunk=0; chunk<allKern.size(); chunk++ ) 
  {
    for ( size_t j=0; j<allKern[chunk].numGLL; j++ ) 
    {
      kernStore[totIter] = allKern[chunk].kernStore[j];
      xExt[totIter]      = allKern[chunk].xExt[j];
      yExt[totIter]      = allKern[chunk].yExt[j];
      zExt[totIter]      = allKern[chunk].zExt[j];
      totIter++;
    }

    // Free memory associated with the original kernal storage 
    delete [] allKern[chunk].kernStore;
    delete [] allKern[chunk].xExt;
    delete [] allKern[chunk].yExt;
    delete [] allKern[chunk].zExt;
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

      NcVar gllValues = dataFile.getVar  ("kernel");
      NcDim dim       = gllValues.getDim (0);
      numGLL          = dim.getSize      ();

      kernStore = new float [numGLL];
   
#if defined (withMPI)
      if ( MPI::COMM_WORLD.Get_rank() == 0 )
      {
        gllValues.getVar (kernStore);
      }
      MPI::COMM_WORLD.Bcast ( kernStore, numGLL, MPI_FLOAT, 0 );
#else
      gllValues.getVar (kernStore);
#endif

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
#if defined (withMPI)
      if ( MPI::COMM_WORLD.Get_rank() == 0 )
      {
        dataX.getVar (xExt);
        dataY.getVar (yExt);
        dataZ.getVar (zExt);
      }
      MPI::COMM_WORLD.Bcast ( xExt, numGLL, MPI_FLOAT, 0 );
      MPI::COMM_WORLD.Bcast ( yExt, numGLL, MPI_FLOAT, 0 );
      MPI::COMM_WORLD.Bcast ( zExt, numGLL, MPI_FLOAT, 0 );
#else
        dataX.getVar (xExt);
        dataY.getVar (yExt);
        dataZ.getVar (zExt);
#endif

    }
    
  } 
  catch (NcException &e)
  {

    std::cout << e.what() << std::endl;
    std::cout << "Failure reading in kernel." << std::endl;
    std::exit ( EXIT_FAILURE );

  }

}

void Kernel::getMinMaxCartesian ()
{

  minX = xExt[0];
  maxX = xExt[0];
  minY = yExt[0];
  maxY = yExt[0];
  minZ = zExt[0];
  maxZ = zExt[0];
  for ( int i=0; i<numGLL; i++ )
  {
    
    if ( xExt[i] < minX )
      minX = xExt[i];
    if ( xExt[i] > maxX )
      maxX = xExt[i];

    if ( yExt[i] < minY )
      minY = yExt[i];
    if ( yExt[i] > maxY )
      maxY = yExt[i];

    if ( zExt[i] < minZ )
      minZ = zExt[i];
    if ( zExt[i] > maxZ )
      maxZ = zExt[i];

  }


}
