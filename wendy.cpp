#include <cstdlib>
#include <netcdf>
#include <math.h>
#include "classes.hpp"
#include "mpi.h"
#include <exodusII.h>

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

void Kernel::createKDtree ( )
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

  centerX = (minX + maxX) / 2.;
  centerY = (minY + maxY) / 2.;
  centerZ = (minZ + maxZ) / 2.;

}

void Kernel::quickSortPoint ( int i1st, int i2nd,
                              float px, float py, float pz )
{

  int pivotElement;

  // Distance of two points from box center.
  float d1st = distFromPoint ( xExt[i1st], yExt[i1st], zExt[i1st], px, py, pz );
  float d2nd = distFromPoint ( xExt[i2nd], yExt[i2nd], zExt[i2nd], px, py, pz );

  if ( i1st < i2nd )
  {
    pivotElement = pivot ( i1st, i2nd, d1st, d2nd );
    quickSortPoint ( i1st, pivotElement-1, px, py, pz );
    quickSortPoint ( pivotElement+1, i2nd, px, py, pz );
  }

}

void Kernel::quickSortCenter ( int i1st, int i2nd )
{

  int pivotElement;

  // Distance of two points from box center.
  float d1st = distFromCenter ( xExt[i1st], yExt[i1st], zExt[i1st] );
  float d2nd = distFromCenter ( xExt[i2nd], yExt[i2nd], zExt[i2nd] );

  if ( i1st < i2nd )
  {
    pivotElement = pivot ( i1st, i2nd, d1st, d2nd );
    quickSortCenter ( i1st, pivotElement-1 );
    quickSortCenter ( pivotElement+1, i2nd );
  }

}

int Kernel::pivot ( int &i1st, int &i2nd, float &d1st, float &d2nd )
{

  int p              = i1st;
  float pivotElement = d1st;


  for ( int i=i1st+1; i<=i2nd; i++ )
  {
    float dTest = distFromCenter ( xExt[i], yExt[i], zExt[i] );
    if ( dTest <= pivotElement )
    {
      p++;
      std::swap ( xExt[i], xExt[p] );
      std::swap ( yExt[i], yExt[p] );
      std::swap ( zExt[i], zExt[p] );
      std::swap ( kernStore[i], kernStore[p] );
    }
  }

  std::swap ( xExt[p], xExt[i1st] );
  std::swap ( yExt[p], yExt[i1st] );
  std::swap ( zExt[p], zExt[i1st] );
  std::swap ( kernStore[p], kernStore[i1st] );

  return p;

}

float Kernel::distFromPoint ( float &x,  float &y,  float &z,
                              float &px, float &py, float &pz )
{

  float diffX = ( x - px );
  float diffY = ( y - py );
  float diffZ = ( z - pz );

  float dist  = sqrt ( diffX * diffX + diffY * diffY + diffZ * diffZ );

  return dist;

}

float Kernel::distFromCenter ( float &x, float &y, float &z )
{

  float diffX = ( x - centerX );
  float diffY = ( y - centerY );
  float diffZ = ( z - centerZ );

  float dist  = sqrt ( diffX * diffX + diffY * diffY + diffZ * diffZ );

  return dist;

}


void Kernel::writeExodus ( )
{

  int comp_ws = sizeof(float);
  int io_ws   = 0;

  int numNodes    = NX*NY*NZ;
  int *nodeNumArr = new int [numNodes];


  float *nodeCorZ = new float [numNodes];
  float *nodeCorY = new float [numNodes];
  float *nodeCorX = new float [numNodes];

  int it = 0;
  for ( int i=0; i<NX; i++ ) {
    for ( int j=0; j<NY; j++ ) {
      for ( int k=0; k<NZ; k++ ) {
    
        nodeNumArr[it] = it+1;
        nodeCorZ[it]   = regZ[k];
        nodeCorY[it]   = regY[j];
        nodeCorX[it]   = regX[i];
        it++;

      }
    }
  }

  int numElem    = (NX-1)*(NY-1)*(NZ-1);
  int *connect = new int [numElem*8];

  int count=0;
  for ( int i=0; i<NX-1; i++ ) {
    for ( int j=0; j<NY-1; j++ ) {
      for ( int k=0; k<NZ-1; k++ ) {

        connect[count]   = nodeNumArr[k+(NZ)*(j+i*NY)];
        connect[count+1] = nodeNumArr[k+(NZ)*(j+i*NY)+NZ];
        connect[count+2] = nodeNumArr[k+(NZ)*(j+i*NY)+NY*NZ+NZ];
        connect[count+3] = nodeNumArr[k+(NZ)*(j+i*NY)+NY*NZ];
        connect[count+4] = nodeNumArr[(k+1)+(NZ)*(j+i*NY)];
        connect[count+5] = nodeNumArr[(k+1)+(NZ)*(j+i*NY)+NZ];
        connect[count+6] = nodeNumArr[(k+1)+(NZ)*(j+i*NY)+NY*NZ+NZ];
        connect[count+7] = nodeNumArr[(k+1)+(NZ)*(j+i*NY)+NY*NZ];

        count=count+8;

      }
    }
  }

  int idexo = ex_create        ( "./test.ex2", EX_CLOBBER, &comp_ws, &io_ws );
  int ier   = ex_put_init      ( idexo, "Thing", 3, NX*NY*NZ, (NX-1)*(NY-1)*(NZ-1), 1, 0, 0 );
  ex_put_coord ( idexo, nodeCorX, nodeCorY, nodeCorZ );
  ex_put_elem_block ( idexo, 1, "HEX", numElem, 8, 0 );
  ex_put_node_num_map ( idexo, nodeNumArr );
  ex_put_elem_conn  ( idexo, 1, connect );
      ier   = ex_put_var_param ( idexo, "n", 1 );
      ier   = ex_put_nodal_var ( idexo, 1, 1, numNodes, regMeshArr ); 
  
  ex_close ( idexo );


}
  




