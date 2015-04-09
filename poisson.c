#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

/* function prototypes */
double *createDoubleArray (int n);
double **createDouble2DArray (int grid_size, int n);
void transpose (double **transposed_grid, double **grid, int grid_size);
void fst_(double *v, int *n, double *w, int *nn);
void fstinv_(double *v, int *n, double *w, int *nn);

double f( double x, double y )
{
	return 1.0; //change later
}

int main(int argc, char **argv )
{
  double *diagonal, **grid, **transposed_grid, *z;
  double pi, point_distance, umax;
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

  diagonal = createDoubleArray (grid_size);
  grid     = createDouble2DArray (grid_size,grid_size);
  transposed_grid   = createDouble2DArray (grid_size,grid_size);
  z        = createDoubleArray (nn);

  point_distance    = 1./(double)n;
  pi   = 4.*atan(1.);

  #pragma omp parallel for
  for (i=0; i < grid_size; i++) {
    diagonal[i] = 2.*(1.-cos((i+1)*pi/(double)n));
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

void transpose (double **transposed_grid, double **grid, int grid_size)
{
  int i, j;
  for (j=0; j < grid_size; j++) {
    for (i=0; i < grid_size; i++) {
      transposed_grid[j][i] = grid[i][j];
    }
  }
}

double *createDoubleArray (int n)
{
  double *a;
  int i;
  a = (double *)malloc(n*sizeof(double));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

double **createDouble2DArray (int n1, int n2)
{
  int i, n;
  double **a;
  a    = (double **)malloc(n1   *sizeof(double *));
  a[0] = (double  *)malloc(n1*n2*sizeof(double));
  for (i=1; i < n1; i++) {
    a[i] = a[i-1] + n2;
  }
  n = n1*n2;
  memset(a[0],0,n*sizeof(double));
  return (a);
}
