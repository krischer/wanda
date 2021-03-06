#include "classes.hpp"
#include "mpi.h"

#include <cstdlib>
#include <netcdf>
#include <cmath>
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

void Kernel::smoothGaussian ( float &var )
{

  // Gaussian constants.
  float denom     = pow ( sqrt (2*M_PI), 3 );// * pow ( var, 3 );
  float norm      = sqrt (var) / denom;
  float varSquare = var*var;

  // Get rank for MPI looping purposes.
  int rank = MPI::COMM_WORLD.Get_rank();

  // Array which will hold smoothed image. Initialize to zero.
  smoothArr = new float [NX*NY*NZ]();

  // Scratch gaussian array. Initialize to zero.
  float *gaussArr = new float [NX*NY*NZ]();

  // Assign initial ranks.
  size_t iStartMPI = rank;
  size_t iEndMPI   = rank+1;

  bool loop = true;
  if ( iEndMPI > NX )
    loop = false;

  clock_t begin = std::clock();

  if ( rank == 0 )
    std::cout << "Launching gaussian smoother." << std::flush << std::endl;

  while ( loop == true )
  {
    for ( size_t iSmooth=iStartMPI; iSmooth<iEndMPI; iSmooth++ ) {
      for ( size_t jSmooth=0; jSmooth<NY; jSmooth++ ) {
        for ( size_t kSmooth=0; kSmooth<NZ; kSmooth++ ) {

          int indSmooth  = kSmooth + NZ * (jSmooth + iSmooth * NY);
          float iSmoothX = regX[iSmooth];
          float iSmoothY = regY[jSmooth];
          float iSmoothZ = regZ[kSmooth];
          float fullSum  = 0.;
          for ( size_t i=0; i<NX; i++ ) {
            for ( size_t j=0; j<NY; j++ ) {
              for ( size_t k=0; k<NZ; k++ ) {

                // Calculate the indices for the smoothed, and source, mesh.
                int ind = k + j * (NZ) + i * ((NZ) * (NY));
                
                // Get the distance from arbitrary grid point to 
                // source smoothing point.
                float xDist = regX[i] - iSmoothX; 
                float yDist = regY[j] - iSmoothY; 
                float zDist = regZ[k] - iSmoothZ;

                // If we're within 5 sigma, smooth.
                gaussArr[ind] = 0.;
                float magnitude = xDist*xDist + yDist*yDist + zDist*zDist;
                if ( sqrt (magnitude) <= (3 * var) ) {

                  float shape   = (-1) * (magnitude) / (2 * varSquare);

                  gaussArr[ind] = exp (shape);
                  fullSum      += gaussArr[ind];

                }

              }
            }
          }

          for ( size_t i=0; i<NX; i++ ) {
            for ( size_t j=0; j<NY; j++ ) {
              for ( size_t k=0; k<NZ; k++ ) {

                int ind               = k + j * (NZ) + i * ((NZ) * (NY));
                smoothArr[indSmooth] += gaussArr[ind] * regMeshArr[ind] / fullSum;

              }
            }
          }
        }
      }
    }

    iStartMPI += MPI::COMM_WORLD.Get_size();
    iEndMPI    = iStartMPI + 1;
    if ( iEndMPI > NX )
      loop = false;
  }

  // Free scratch gaussian array.
  delete [] gaussArr;

  // Bring all together.
  MPI::COMM_WORLD.Barrier ();
  clock_t end = std::clock();
  double elapsed_secs = double (end - begin) / CLOCKS_PER_SEC;
  if ( MPI::COMM_WORLD.Get_rank() == 0 )
    std::cout << "The smoothing took: " << elapsed_secs
      << " seconds." << std::flush << std::endl;

  MPI::COMM_WORLD.Allreduce ( MPI_IN_PLACE, smoothArr, NX*NY*NZ, MPI_FLOAT, MPI_SUM );

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

void Kernel::crossProductZ ( )
{

  rotVec.reserve (3);

  float x1 = centerX;
  float y1 = centerY;
  float z1 = centerZ;
  float x2 = 0;
  float y2 = 0;
  float z2 = 1;

  rotVec[0] = y1 * z2 - z1 * y2;
  rotVec[1] = z1 * x2 - x1 * z2;
  rotVec[2] = x1 * y2 - y1 * x2;

  float mag = sqrt ( rotVec[0]*rotVec[0] + rotVec[1]*rotVec[1] + rotVec[2]*rotVec[2] );
  rotVec[0] = rotVec[0] / mag;
  rotVec[1] = rotVec[1] / mag;
  rotVec[2] = rotVec[2] / mag;
  
}

void Kernel::rotateArbitraryVector ( float &angle )
{

  float a = angle * M_PI / 180.;

  float x = rotVec[0];
  float y = rotVec[1];
  float z = rotVec[2];

  float rot11 = cos(a) + (x * x) * (1 - cos(a));
  float rot21 = z * sin(a) + x * y * (1 - cos(a));
  float rot31 = y * sin(a) + x * z * (1 - cos(a));
  float rot12 = x * y * (1 - cos(a)) - z * sin(a);
  float rot22 = cos(a) + (y * y) * (1 - cos(a));
  float rot32 = x * sin(a) + y * z * (1 - cos(a));
  float rot13 = y * sin(a) + x * z * (1 - cos(a));
  float rot23 = x * sin(a) + y * z * (1 - cos(a));
  float rot33 = cos(a) + (z * x) * (1 - cos(a));

  rot23 = (-1) * rot23;
  rot31 = (-1) * rot31;

  regXRot = new float [numGLL];
  regYRot = new float [numGLL];
  regZRot = new float [numGLL];

  for ( size_t i=0; i<numGLL; i++ ) {

    regXRot[i] = rot11 * xExt[i] + rot21 * yExt[i] + rot31 * zExt[i];
    regYRot[i] = rot12 * xExt[i] + rot22 * yExt[i] + rot32 * zExt[i];
    regZRot[i] = rot13 * xExt[i] + rot23 * yExt[i] + rot33 * zExt[i];
  }

  for ( size_t i=0; i<numGLL; i++ ) {
    xExt[i] = regXRot[i];
    yExt[i] = regYRot[i];
    zExt[i] = regZRot[i];
  }

  delete [] regXRot;
  delete [] regYRot;
  delete [] regZRot;

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

  centerC = (minC + maxC) / 2.;
  centerL = (minL + maxL) / 2.;
  centerR = (minR + maxR) / 2.;

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

  // WARNING WARNING TODO.
  centerX = 0.;
  centerY = 0.;
  centerZ = centerR; 

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
        count           += 8;

      }
    }
  }

  std::cout << "Writing exodus file." << std::flush << std::endl;

  // Interpolated array.
  int idexo = ex_create ( "./testInterp.ex2", EX_CLOBBER, &comp_ws, &io_ws );
  exodusErrorCheck ( ex_put_init( idexo, "Kernel", nDim, numNodes, numElem, nBlock, nNodeSet, nSideSet ), "ex_put_init" );
  exodusErrorCheck ( ex_put_coord ( idexo, nodeCorX, nodeCorY, nodeCorZ ), "ex_put_coord" );
  exodusErrorCheck ( ex_put_elem_block ( idexo, 1, "HEX", numElem, nNodePerElem, 0 ), "ex_put_elem_block" );
  exodusErrorCheck ( ex_put_node_num_map ( idexo, nodeNumArr ), "ex_put_node_num_map" );
  exodusErrorCheck ( ex_put_elem_conn  ( idexo, 1, connect ), "ex_put_elem_conn" );
  exodusErrorCheck ( ex_put_var_param ( idexo, "n", nVars ), "ex_put_var_param" );
  exodusErrorCheck ( ex_put_var_names ( idexo, "n", nVars, varNames ), "ex_put_var_names" );
  exodusErrorCheck ( ex_put_nodal_var ( idexo, 1, 1, numNodes, regMeshArr ), "ex_put_nodal_var" ); 
  exodusErrorCheck ( ex_close ( idexo ), "ex_close" );

  // Smoothed array.
  idexo = ex_create ( "./testSmooth.ex2", EX_CLOBBER, &comp_ws, &io_ws );
  exodusErrorCheck ( ex_put_init( idexo, "Kernel", nDim, numNodes, numElem, nBlock, nNodeSet, nSideSet ), "ex_put_init" );
  exodusErrorCheck ( ex_put_coord ( idexo, nodeCorX, nodeCorY, nodeCorZ ), "ex_put_coord" );
  exodusErrorCheck ( ex_put_elem_block ( idexo, 1, "HEX", numElem, nNodePerElem, 0 ), "ex_put_elem_block" );
  exodusErrorCheck ( ex_put_node_num_map ( idexo, nodeNumArr ), "ex_put_node_num_map" );
  exodusErrorCheck ( ex_put_elem_conn  ( idexo, 1, connect ), "ex_put_elem_conn" );
  exodusErrorCheck ( ex_put_var_param ( idexo, "n", nVars ), "ex_put_var_param" );
  exodusErrorCheck ( ex_put_var_names ( idexo, "n", nVars, varNames ), "ex_put_var_names" );
  exodusErrorCheck ( ex_put_nodal_var ( idexo, 1, 1, numNodes, smoothArr ), "ex_put_nodal_var" ); 
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
