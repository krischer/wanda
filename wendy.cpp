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

void Kernel::createKDtree ( std::vector<Kernel> &allKern )
{

  // Initialize the KDtree.
  tree = kd_create (3);

  // Determine total number of gll points.
  int totalGLL = 0;
  for ( size_t i=0; i<allKern.size(); i++ )
    totalGLL += allKern[i].numGLL;

  // Initizalize the data array.
  KDdat = new int [totalGLL];

  // Populate the tree with an index (KDdat) and x, y, z
  // We need to loop over both mesh chunks and gll points.
  int totIter = 0;
  for ( size_t chunk=0; chunk<allKern.size(); chunk++ ) 
  {

    std::cout << "Creating KDtree for chunk " << chunk << "." << std::flush << std::endl;
    for ( size_t i=0; i<allKern[chunk].numGLL; i++ )
    {
      KDdat[totIter] = totIter;
      kd_insert3 ( tree, 
                   allKern[chunk].xExt[i], 
                   allKern[chunk].yExt[i], 
                   allKern[chunk].zExt[i], 
                   &KDdat[i] );
      totIter++;
    }

    /* After we've created the KDtree, we no longer need the
     * inital coordinate arrays 
     */
    delete [] allKern[chunk].xExt;
    delete [] allKern[chunk].yExt;
    delete [] allKern[chunk].zExt;

  } 

}

void Kernel::mergeKernels ( std::vector<Kernel> &allKern )
{

  // Determine total number of gll points.
  int totalGLL = 0;
  for ( size_t i=0; i<allKern.size(); i++ )
    totalGLL += allKern[i].numGLL;

  // Initialize parameter array.
  kernStore = new float [totalGLL];

  // Copy to master parameter array.
  int totIter = 0;
  for ( size_t i=0; i<allKern.size(); i++ ) 
  {
    for ( size_t j=0; j<allKern[i].numGLL; j++ ) 
    {
      kernStore[totIter] = allKern[i].kernStore[j];
      totIter++;
    }

    /* Free memory associated with the original kernal storage
     */
    delete [] allKern[i].kernStore;
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
