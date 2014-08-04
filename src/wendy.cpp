#include <cstdlib>
#include <netcdf>
#include <cmath>
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

void Kernel::rotateXaxis ( double &deg )
{

  double rad = deg * M_PI / 180.0;

  double rx11 = 1;
  double rx12 = 0;
  double rx13 = 0;
  double rx21 = 0;
  double rx22 = cos (rad);
  double rx23 = sin (rad) * (-1);
  double rx31 = 0;
  double rx32 = sin (rad);
  double rx33 = cos (rad);

  for ( size_t i=0; i<numGLL; i++ )
  {
    double xNew = rx11 * xExt[i] + rx12 * yExt[i] + rx13 * zExt[i];
    double yNew = rx21 * xExt[i] + rx22 * yExt[i] + rx23 * zExt[i];
    double zNew = rx31 * xExt[i] + rx32 * yExt[i] + rx33 * zExt[i];

    xExt[i] = xNew;
    yExt[i] = yNew;
    zExt[i] = zNew;
  }

}

void Kernel::rotateYaxis ( double &deg )
{

  double rad = deg * M_PI / 180.0;

  double ry11 = cos (rad);
  double ry12 = 0;
  double ry13 = sin (rad);
  double ry21 = 0;
  double ry22 = 1;
  double ry23 = 0;
  double ry31 = sin (rad) * (-1);
  double ry32 = 0;
  double ry33 = cos (rad);

  for ( size_t i=0; i<numGLL; i++ )
  {
    double xNew = ry11 * xExt[i] + ry12 * yExt[i] + ry13 * zExt[i];
    double yNew = ry21 * xExt[i] + ry22 * yExt[i] + ry23 * zExt[i];
    double zNew = ry31 * xExt[i] + ry32 * yExt[i] + ry33 * zExt[i];


    xExt[i] = xNew;
    yExt[i] = yNew;
    zExt[i] = zNew;
  }

}

void Kernel::rotateZaxis ( double &deg )
{

  double rad = deg * M_PI / 180.0;

  double  rz11 = cos (rad);
  double  rz12 = sin (rad) * (-1);
  double  rz13 = 0;
  double  rz21 = sin (rad);
  double  rz22 = cos (rad);
  double  rz23 = 0;
  double  rz31 = 0;
  double  rz32 = 0;
  double  rz33 = 1;

  for ( size_t i=0; i<numGLL; i++ )
  {

    double xNew = rz11 * xExt[i] + rz12 * yExt[i] + rz13 * zExt[i];
    double yNew = rz21 * xExt[i] + rz22 * yExt[i] + rz23 * zExt[i];
    double zNew = rz31 * xExt[i] + rz32 * yExt[i] + rz33 * zExt[i];

    xExt[i] = xNew;
    yExt[i] = yNew;
    zExt[i] = zNew;
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

double Kernel::angleFromAxis ( char axis )
{

  double v1 = 0;
  double v2 = 0;
  double v3 = 0;

  if ( axis == 'x' )
  {
    v1 = 1;
    v2 = 0;
    v3 = 0;
  }

  if ( axis == 'y' )
  {
    v1 = 0;
    v2 = 1;
    v3 = 0;
  }

  if ( axis == 'z' )
  {
    v1 = 0;
    v2 = 0;
    v3 = 1;
  }

  double uDotv = centerX * v1 + centerY * v2 + centerZ * v3;
  double magU  = sqrt ( centerX * centerX + centerY * centerY + centerZ * centerZ );
  double div   = uDotv / magU;
  double angle = acos (div) * 180.0 / M_PI;

  return angle;

}




void Kernel::getMinMaxCartesian ()
{

  minX = xExt[0];
  maxX = xExt[0];
  minY = yExt[0];
  maxY = yExt[0];
  minZ = zExt[0];
  maxZ = zExt[0];
  minC = 180;
  maxC = 0;
  minL = 180;
  maxL = -180;
  minR = 6371;
  maxR = 0;
  for ( int i=0; i<numGLL; i++ )
  {

    float rad = sqrt  ( xExt[i] * xExt[i] + yExt[i] * yExt[i] + zExt[i] * zExt[i] );
    float lon = atan2 ( yExt[i] , xExt[i] );
    float col = acos  ( zExt[i] / rad );

    if ( rad < minR )
      minR = rad;
    if ( rad > maxR )
      maxR = rad;

    if ( lon < minL )
      minL = lon;
    if ( lon > maxL )
      maxL = lon;

    if ( col < minC )
      minC = col;
    if ( col > maxC )
      maxC = col;
    
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

  // Exodus parameters.
  int  comp_ws      = sizeof(float);
  int  io_ws        = 0;
  int  nVars        = 1;
  int  nDim         = 3;
  int  nBlock       = 1;
  int  nNodeSet     = 0;
  int  nSideSet     = 0;
  int  nNodePerElem = 8;
  int  numNodes     = NX*NY*NZ;
  int  numElem      = (NX-1)*(NY-1)*(NZ-1);

  char *varNames[1];
  varNames[0] = "Sensitivity";
  
  int *nodeNumArr = new int [numNodes];
  float *nodeCorZ = new float [numNodes];
  float *nodeCorY = new float [numNodes];
  float *nodeCorX = new float [numNodes];

  // Unpack coordinate arrays.
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

  int *connect = new int [numElem*nNodePerElem];

  // Construct connectivity array.
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


  std::cout << "Writing exodus file." << std::flush << std::endl;
  int idexo = ex_create ( "./test.ex2", EX_CLOBBER, &comp_ws, &io_ws );

  exodusErrorCheck ( ex_put_init( idexo, "Kernel", nDim, numNodes, numElem, nBlock, nNodeSet, nSideSet ), "ex_put_init" );
  exodusErrorCheck ( ex_put_coord ( idexo, nodeCorX, nodeCorY, nodeCorZ ), "ex_put_coord" );
  exodusErrorCheck ( ex_put_elem_block ( idexo, 1, "HEX", numElem, nNodePerElem, 0 ), "ex_put_elem_block" );
  exodusErrorCheck ( ex_put_node_num_map ( idexo, nodeNumArr ), "ex_put_node_num_map" );
  exodusErrorCheck ( ex_put_elem_conn  ( idexo, 1, connect ), "ex_put_elem_conn" );
  exodusErrorCheck ( ex_put_var_param ( idexo, "n", nVars ), "ex_put_var_param" );
  exodusErrorCheck ( ex_put_var_names ( idexo, "n", nVars, varNames ), "ex_put_var_names" );
  exodusErrorCheck ( ex_put_nodal_var ( idexo, 1, 1, numNodes, regMeshArr ), "ex_put_nodal_var" ); 
  exodusErrorCheck ( ex_close ( idexo ), "ex_close" );

}

void exodusErrorCheck ( int ier, std::string function )
{

  if ( ier != 0 )
  {
    std::cout << "ERROR in " << function << std::flush << std::endl;
    exit (EXIT_FAILURE);
  }

}
