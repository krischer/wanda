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
    std::cout << "How many processors: ";
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
    std::string fNameApp    = "_reg1_betav_kernel.nc";
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

  // Quicksort by distance from center.
  if ( rank == 0 )
    std::cout << "Sorting." << std::flush << std::endl;
  
  kern.quickSortCenter ( 0, kern.numGLL-1 );

  // Create a single kdtree from the baby kernels.
  kern.createKDtree ( );

  // Go ahead and create the regular mesh.
  createRegMesh     ( kern, allKern );

  if ( rank == 0 )
    kern.writeExodus ( );

  MPI::Finalize();
  return 0;

}

void createRegMesh ( Kernel &kern, std::vector<Kernel> &allKern )
{

  // Grid discretization. Make this a parameter.
  float DX = 50.;
  float DY = 50.;
  float DZ = 50.;


  // Determine # of grid points in each direction, given DX.
  long NX = int ( (((kern.maxX+DX) - kern.minX) / DX) + 1);
  long NY = int ( (((kern.maxY+DY) - kern.minY) / DY) + 1);
  long NZ = int ( (((kern.maxZ+DZ) - kern.minZ) / DZ) + 1);

  // Report on the grid.
  int gridSize = NX*NY*NZ;
  if ( MPI::COMM_WORLD.Get_rank() == 0 )
    std::cout << "The grid is: [ " << NX << ", " << NY << ", " << NZ << " ]." << std::flush
      << std::endl;

  // Allocate space for regular grid (initialized to 0).
  kern.regMeshArr = new float [gridSize]();
  kern.regX       = new float [NX]();
  kern.regY       = new float [NY]();
  kern.regZ       = new float [NZ]();

  // Determine the parts of the mesh sent to each processor.
  float *sliceDistForward = new float [NX];
  float *sliceDistSort    = new float [NX];
  for ( int i=0; i<NX; i++ )
  {
    float x = kern.minX + i * DX;
    float y = kern.minY + NY/2 * DY;
    float z = kern.minZ + NZ/2 * DZ;
    sliceDistForward[i] = kern.distFromCenter ( x, y, z );
    sliceDistSort[i]    = kern.distFromCenter ( x, y, z );
  }

  // Sort the array by (average) distance from center.
  std::sort    ( sliceDistSort, sliceDistSort+NX );

  // Sort the order of MPI processes by distance (makes it much faster).
  int iStartMPI = 0;
  int iEndMPI   = 0;
  int loopRank  = MPI::COMM_WORLD.Get_rank();
  for ( size_t i=0; i<NX; i++ )
  {
    if ( sliceDistSort[loopRank] == sliceDistForward[i] )
    {
      iStartMPI = i;
      iEndMPI   = i+1;
    }
  }

  bool loop = true;
  if ( iStartMPI == 0 && iEndMPI == 0 )
    loop = false;

  // Time the interpolation.
  clock_t begin = std::clock();

  // Loop over the 3D regular grid.
  if ( loopRank == 0 )
    std::cout << "Launching interpolation." << std::flush << std::endl;

  while ( loop == true )
  {

    for ( size_t i=iStartMPI; i<iEndMPI; i++ ) {
      for ( size_t j=0; j<NY; j++ ) {
        for ( size_t k=0; k<NZ; k++ ) {

          // Calculate the position in the grid.
          float x = kern.minX + i * DX;
          float y = kern.minY + j * DY;
          float z = kern.minZ + k * DZ;

          // Use these to determine if we're in the pad box or not.
          float rad = sqrt  ( x*x + y*y + z*z );
          float lon = atan2 ( y, x );
          float col = acos  ( z / rad );

          // Assume we're outside the box, and test
          bool inRad = false;
          bool inLon = false;
          bool inCol = false;
          if ( (rad < kern.maxR) && ( rad > kern.minR) )
            inRad = true;
          if ( (lon < kern.maxL) && ( lon > kern.minL) )
            inLon = true;
          if ( (col < kern.maxC) && ( col > kern.minC) )
            inCol = true;

          // This value will be either taken from kernel (if we're inside)
          // or kept at 0 (if we're outside).
          float value = 0.;
          if ( inRad == true && inLon == true && inCol == true )
          {
            // Extract the closest point from the complete kdtree.
            kdres *set = kd_nearest3      ( kern.tree, x, y, z );
            void  *ind = kd_res_item_data ( set );
            int    pnt =                * ( int * ) ind;
            value      = kern.kernStore[pnt];
            kd_res_free ( set );
          }

          // Store both the kernel and the coordinates.
          int index = k + j * (NZ) + i * ((NZ) * (NY));
          kern.regMeshArr[index] = value;
          kern.regX[i]           = x;
          kern.regY[j]           = y;
          kern.regZ[k]           = z;

        }
      }

    }

    // Increment the loop counter.
    loopRank += MPI::COMM_WORLD.Get_size();
    if ( loopRank < NX )
    {
      for ( size_t i=0; i<NX; i++ )
      {
        if ( sliceDistSort[loopRank] == sliceDistForward[i] )
        {
          iStartMPI = i;
          iEndMPI   = i+1;
        }
      }
    }
    else
    {
      loop = false;
    }

  }

  // Barrier here for timer.
  MPI::COMM_WORLD.Barrier ();
  clock_t end = std::clock();
  double elapsed_secs = double (end - begin) / CLOCKS_PER_SEC;
  if ( MPI::COMM_WORLD.Get_rank() == 0 ) 
  {
     std::cout << "The interpolation took: " << elapsed_secs 
       << " seconds."  << std::flush << std::endl;
  }

  // Bring array back together.
  MPI::COMM_WORLD.Allreduce ( MPI_IN_PLACE, kern.regMeshArr, NX*NY*NZ, MPI_FLOAT, MPI_SUM );
  MPI::COMM_WORLD.Allreduce ( MPI_IN_PLACE, kern.regX, NX, MPI_FLOAT, MPI_SUM );

  kern.NX = NX;
  kern.NY = NY;
  kern.NZ = NZ;

}

