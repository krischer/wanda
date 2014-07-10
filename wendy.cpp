#include <cstdlib>
#include "classes.hpp"

void returnRegularArray ( int NUM_X, int NUM_Y, int NUM_Z, void *testArr )
{

  Kernel kern;
  std::cout << "HELLO WORLD." << std::flush << std::endl;

  std::cout << NUM_X << ' ' << NUM_Y << ' ' << NUM_Z << ' ' << testArr << std::endl;

  kern.regMesh = std::malloc ( 10 * 10 * 10 * sizeof(float) );
  float (*dmat)[NUM_X][NUM_Y][NUM_Z] = ( (float (*)[NUM_X][NUM_Y][NUM_Z]) kern.regMesh );
  for ( int i=0; i<NUM_X; i++ ) {
    for ( int j=0; j<NUM_Y; j++ ) {
      for ( int k=0; k<NUM_Z; k++ ) {

        *dmat[i][j][k] = 10.;
        std::cout << *dmat[i][j][k] << std::endl;

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

  // Initizalize the data array
//  KDdat = new int [x.size()];

  // Populate the tree with an index (KDdat) and x, y, z.
//  for ( size_t i=0; i<x.size(); i++ )
//  {
//    KDdat[i] = i;
//    kd_insert3 ( tree, x[i], y[i], z[i], &KDdat[i] );
//  }

}
