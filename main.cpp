#include "classes.hpp"
#include "kdtree.h"

#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>

#include <omp.h>
#include "mpi.h"

int main ()
{

  // Initialize MPI.
#if defined (withMPI)
  MPI::Init ();
  int rank = MPI::COMM_WORLD.Get_rank ();
  int size = MPI::COMM_WORLD.Get_size ();
#else
  int rank = 0;
  int size = 1;
#endif

  // Initialize master Kernel object, and vector of individual baby chunk
  // kernels (will be read in as netcdf files).
  Kernel kern;
  std::vector <Kernel> allKern;

  // Use input for number of processors kernel is coming from. Need to
  // nProcBuf here for MPI broadcast purposes.
  int  nProc;
  int *nProcBuf = &nProc;

  // Read in number of processors from user. Make this a parameter.
  if ( rank == 0 ) 
  {
    std::cout << "How many processors." << std::endl;
    std::cin >> nProc;
  }

#if defined (withMPI)
  // Broadcast number of processor files.
  MPI::COMM_WORLD.Bcast ( nProcBuf, 1, MPI_INT, 0 );
  allKern.reserve (*nProcBuf);
#endif

  // Loop over input netcdf files.
  for ( int i=0; i<*nProcBuf; i++ )
  {

    // Baby kernel (could rename).
    Kernel kern;

    // Directory basenames (make parameter).
    std::string fNameBase   = "./testKernels/proc";
    std::string fNameApp    = "_reg1_betah_kernel.nc";
    std::string xyzNameBase = "./cemRequest/xyz_reg01_proc";

    // Get proc num name for kernel.
    std::stringstream ssPrc;
    ssPrc << std::setw (6) << std::setfill ('0');
    ssPrc << std::to_string (static_cast<long long>(i));

    // Get proc num name for coordinates.
    std::stringstream ssPrcShort;
    ssPrcShort << std::setw (4) << std::setfill ('0');
    ssPrcShort << std::to_string (static_cast<long long>(i));

    // Create full file string.
    fNameBase.append   ( ssPrc.str() );
    fNameBase.append   ( fNameApp );
    xyzNameBase.append ( ssPrcShort.str() );

    // Read in the kernels and coordiantes.
    kern.readNetcdf ( "kernel", fNameBase );
    kern.readNetcdf ( "coordinates", xyzNameBase );

    // Determine the min/max xyz coordinates in each kernel. Then, populate the
    // master kernel array.
    kern.getMinMaxCartesian ();
    allKern.push_back       (kern);

  }


  // Merge the kernel data as well.
  kern.mergeKernels ( allKern );

  // Find min/max box of full kernel.
  kern.getMinMaxCartesian ();

  clock_t begin = std::clock();
  if ( rank == 0 )
    std::cout << "minDistBeforeSort: " << kern.distFromCenter ( kern.xExt[0], kern.yExt[0], kern.zExt[0] ) << " km." << std::endl;

  // Quicksort by distance from center.
  kern.quickSortCenter ( 0, kern.numGLL-1 );

  if ( rank == 0 )
    std::cout << "minDistAfterSort:  " << kern.distFromCenter ( kern.xExt[0], kern.yExt[0], kern.zExt[0] ) << " km." << std::endl;

  clock_t end = std::clock();
  double elapsed_secs = double (end - begin) / CLOCKS_PER_SEC;
  if ( rank == 0 )
    std::cout << "hey mike, sorting [quicksort] took: " << elapsed_secs << " seconds. [scaled up]" 
      << std::flush << std::endl;

//  for ( int i=0; i<kern.numGLL-1; i++ )
//  {
//    float dist1 = kern.distFromCenter ( kern.xExt[i], kern.yExt[i], kern.zExt[i] );
//    float dist2 = kern.distFromCenter ( kern.xExt[i+1], kern.yExt[i+1], kern.zExt[i+1] );
//    if ( dist2 < dist1)
//      std::cout << "BADBADBAD" << std::endl;
//  }

  // Create a single kdtree from the baby kernels.
  kern.createKDtree ( );

  // Go ahead and create the regular mesh.
  createRegMesh     ( kern, allKern );

  return 0;

}

void createRegMesh ( Kernel &kern, std::vector<Kernel> &allKern )
{

  // Grid discretization. Make this a parameter.
  float DX = 10.;
  float DY = 10.;
  float DZ = 10.;

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

  // Report on the grid.
  int gridSize = NX*NY*NZ;
  if ( MPI::COMM_WORLD.Get_rank() == 0 )
    std::cout << "The grid is [ " << NX << ", " << NY << ", " << NZ << " ]" << std::flush
      << std::endl;

  // Allocate space for regular grid (initialized to 0).
  kern.regMeshArr = new float [gridSize]();
  kern.regX       = new float [gridSize]();
  kern.regY       = new float [gridSize]();
  kern.regZ       = new float [gridSize]();

  // Determine the parts of the mesh sent to each processor.
  float *sliceDistForward = new float [NX];
  float *sliceDistSort    = new float [NX];
  for ( int i=0; i<NX; i++ )
  {
    float x = kern.minXBox + i * DX;
    float y = kern.minYBox + NY/2 * DY;
    float z = kern.minZBox + NZ/2 * DZ;
    sliceDistForward[i] = kern.distFromCenter ( x, y, z );
    sliceDistSort[i]    = kern.distFromCenter ( x, y, z );
  }

  // Sort the array by (average) distance from center.
  std::sort    ( sliceDistSort, sliceDistSort+NX );

  int iStartMPI = 0;
  int iEndMPI   = 0;
  for ( int i=0; i<NX; i++ )
  {
    if ( sliceDistSort[MPI::COMM_WORLD.Get_rank()] == sliceDistForward[i] )
    {
      iStartMPI = i;
      iEndMPI   = i+1;
    }
  }

  clock_t begin = std::clock();

  // Loop over the 3D regular grid.
  while ( iEndMPI <= NX )
  {

#pragma omp parallel for schedule (dynamic) 
    for ( size_t i=iStartMPI; i<iEndMPI; i++ ) {

      // For the curious: a report.
      int size      = MPI::COMM_WORLD.Get_size();
      int rank      = MPI::COMM_WORLD.Get_rank();
      int threadNum = omp_get_thread_num();
      int threadSiz = omp_get_num_threads();
      if ( (rank == 0) && (threadNum == 0) ) 
        std::cout << "Hi there. I've just launched interpolation on "
          << size << " nodes, using " << threadSiz << " threads per "
          << "node. Enjoy your day!" << std::flush << std::endl;

      for ( size_t j=0; j<NY; j++ ) {
        for ( size_t k=0; k<NZ; k++ ) {

          // Calculate the position in the grid.
          float x = kern.minXBox + i * DX;
          float y = kern.minYBox + j * DY;
          float z = kern.minZBox + k * DZ;

          // Extract the closest point from the complete kdtree.
          kdres *set = kd_nearest3      ( kern.tree, x, y, z );
          void  *ind = kd_res_item_data ( set );
          int    pnt =                * ( int * ) ind;
          kd_res_free                   ( set );

          // Store both the kernel and the coordinates.
          int index = k + j * (NZ) + i * ((NZ) * (NY));
          kern.regMeshArr[index] = kern.kernStore[pnt];
          kern.regX      [index] = x;
          kern.regY      [index] = y;
          kern.regZ      [index] = z;

        }
      }
    }

    iStartMPI += MPI::COMM_WORLD.Get_size();
    iEndMPI    = iStartMPI + 1;

  }

  std::cout << "I'm done: " << MPI::COMM_WORLD.Get_rank() << ' ' << iStartMPI << ' ' << iEndMPI << std::flush << std::endl;

  MPI::COMM_WORLD.Barrier ();
  clock_t end = std::clock();
  double elapsed_secs = double (end - begin) / CLOCKS_PER_SEC;
  if ( MPI::COMM_WORLD.Get_rank() == 0 )
    std::cout << "The interpolation took: " << elapsed_secs << " seconds." 
      << std::flush << std::endl;

  MPI::COMM_WORLD.Allreduce ( MPI_IN_PLACE, kern.regMeshArr, NX*NY*NZ, MPI_FLOAT, MPI_SUM );
  MPI::COMM_WORLD.Allreduce ( MPI_IN_PLACE, kern.regX, NX*NY*NZ, MPI_FLOAT, MPI_SUM );
  MPI::COMM_WORLD.Allreduce ( MPI_IN_PLACE, kern.regY, NX*NY*NZ, MPI_FLOAT, MPI_SUM );
  MPI::COMM_WORLD.Allreduce ( MPI_IN_PLACE, kern.regZ, NX*NY*NZ, MPI_FLOAT, MPI_SUM );

}

//void writeExodus ()
//{
//
//  int exoid = ex_create ("./test.ex2", EX_CLOBBER, 4, 4 );
//  int ier   = ex_put_init ( idexo, "Title", 3, kern.totGLL, 0, 0, 0 );
//  int ier   = ex_put_var_param ( idexo, "n", 1 );
//  int ier   = ex_put_nodal_var ( idexo, 1, 1, NX*NY*NZ, 
//
//  float (*regMesh)[NX][NY][NZ] = ( (float (*)[NX][NY][NZ]) kern.regMeshArr );
  




