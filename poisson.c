#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int grid_size, int n);
void transpose (Real **transposed_grid, Real **grid, int grid_size);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);

Real f( Real x, Real y )
{
	return 1.0; //change later
}

int main(int argc, char **argv )
{
  Real *diagonal, **grid, **transposed_grid, *z;
  Real pi, point_distance, umax;
  int i, j, n, grid_size, nn, k, processors, threads, per_proc, world_size;

  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */

 if( argc < 2 ) {
    printf("need a problem size\n");
    return 1;
  }

  n = pow( 2, atoi( argv[1]) );
  processors = atoi( argv[2] );
  threads    = atoi( argv[3] );
  grid_size  = n-1;
  nn = 4*n;

  per_proc = grid_size / processors;
  omp_set_num_threads( threads );

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &processors );
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );

  diagonal = createRealArray (grid_size);
  grid     = createReal2DArray (grid_size,grid_size);
  transposed_grid   = createReal2DArray (grid_size,grid_size);
  z        = createRealArray (nn);

  point_distance    = 1./(Real)n;
  pi   = 4.*atan(1.);

  #pragma omp parallel for
  for (i=0; i < grid_size; i++) {
    diagonal[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
  }

  #pragma omp parallel for
  for (j=0; j < grid_size; j++) {
    for (i=0; i < grid_size; i++) {
      grid[j][i] = point_distance * point_distance * f(j,i);
    }
  }

  #pragma omp parallel for
  for (j=0; j < grid_size; j++) {
    fst_(grid[j], &n, z, &nn);
  }

  transpose( transposed_grid, grid, grid_size );

  #pragma omp parallel for
  for (i=0; i < grid_size; i++) {
    fstinv_( transposed_grid[i], &n, z, &nn );
  }
  
  #pragma omp parallel for
  for (j=0; j < grid_size; j++) {
    for (i=0; i < grid_size; i++) {
      transposed_grid[j][i] = transposed_grid[j][i] / (diagonal[i]+diagonal[j]);
    }
  }
  
  #pragma omp parallel for
  for (i=0; i < grid_size; i++) {
    fst_(transposed_grid[i], &n, z, &nn);
  }

  transpose( grid, transposed_grid, grid_size );

  #pragma omp parallel for
  for (j=0; j < grid_size; j++) {
    fstinv_( grid[j], &n, z, &nn );
  }

  umax = 0.0;
  #pragma omp parallel for
  for (j=0; j < grid_size; j++) {
    for (i=0; i < grid_size; i++) {
      if (grid[j][i] > umax) 
		  #pragma omp critical
		  umax = grid[j][i];
    }
  }
  printf (" umax = %e \n",umax);

  return 0;
}

void transpose (Real **transposed_grid, Real **grid, int grid_size)
{
  int i, j;
  for (j=0; j < grid_size; j++) {
    for (i=0; i < grid_size; i++) {
      transposed_grid[j][i] = grid[i][j];
    }
  }
}

Real *createRealArray (int n)
{
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

Real **createReal2DArray (int n1, int n2)
{
  int i, n;
  Real **a;
  a    = (Real **)malloc(n1   *sizeof(Real *));
  a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
  for (i=1; i < n1; i++) {
    a[i] = a[i-1] + n2;
  }
  n = n1*n2;
  memset(a[0],0,n*sizeof(Real));
  return (a);
}
