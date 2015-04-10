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
void parallel_transpose (double **transposed_grid, double **grid, int grid_size, int* rows_in_proc, int proc_id, int proc_num);
void prepare_send (double **grid, double *send_block, int grid_size, int *rows_in_proc, int proc_id, int proc_num);
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
  int i, j, n, grid_size, nn, k, proc_id, threads, proc_num;
  int *rows_in_proc;

  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */

 if( argc < 2 ) {
    printf("need a problem size\n");
    return 1;
  }

  n = pow( 2, atoi( argv[1]) );
  threads    = atoi( argv[2] );
  grid_size  = n-1;
  nn = 4*n;

  diagonal = createDoubleArray (grid_size);
  grid     = createDouble2DArray (grid_size,grid_size);
  transposed_grid   = createDouble2DArray (grid_size,grid_size);
  z        = createDoubleArray (nn);

  point_distance    = 1./(double)n;
  pi   = 4.*atan(1.);

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &proc_num );
  omp_set_num_threads( threads );

  // todo: calculate once, scatter to all 
  rows_in_proc = malloc( sizeof(int)*proc_num );
  int num_rows = grid_size / proc_num;
  int more_rows = grid_size % proc_num;
  
  #pragma omp parallel for
  for( i=0; i<more_rows; i++ ){
	  rows_in_proc[i] = num_rows + 1;
  }
  #pragma omp parallel for
  for( i=more_rows; i<proc_num; i++ ){
	  rows_in_proc[i] = num_rows;
  }

  #pragma omp parallel for
  for (i=0; i < grid_size; i++) {
    diagonal[i] = 2.*(1.-cos((i+1)*pi/(double)n));
  }

  double** test_array = malloc( 7*7*sizeof(int) );


  #pragma omp parallel for schedule(static) private(i)
  for (j=0; j < grid_size; j++) {
    for (i=0; i < grid_size; i++) {
      grid[j][i] = point_distance * point_distance * f(j,i);
    }
  }

  #pragma omp parallel for
  for (j=0; j < grid_size; j++) {
    fst_(grid[j], &n, z, &nn);
  }

  parallel_transpose( transposed_grid, grid, grid_size, rows_in_proc, proc_id, proc_num );
  MPI_Barrier( MPI_COMM_WORLD );

  #pragma omp parallel for
  for (i=0; i < grid_size; i++) {
    fstinv_( transposed_grid[i], &n, z, &nn );
  }
  
  #pragma omp parallel for schedule(static) private(i)
  for (j=0; j < grid_size; j++) {
    for (i=0; i < grid_size; i++) {
      transposed_grid[j][i] = transposed_grid[j][i] / (diagonal[i]+diagonal[j]);
    }
  }
  
  #pragma omp parallel for
  for (i=0; i < grid_size; i++) {
    fst_(transposed_grid[i], &n, z, &nn);
  }

  parallel_transpose( transposed_grid, grid, grid_size, rows_in_proc, proc_id, proc_num );
  MPI_Barrier( MPI_COMM_WORLD );

  #pragma omp parallel for
  for (j=0; j < grid_size; j++) {
    fstinv_( grid[j], &n, z, &nn );
  }

  umax = 0.0;
  // TODO: MPI scatter, gather
  for (j=0; j < grid_size; j++) {
    for (i=0; i < grid_size; i++) {
      if (grid[j][i] > umax) 
		  // TODO use f to calculate error
		  umax = grid[j][i];
    }
  }
  double global_umax = 0.0;
  MPI_Reduce( &umax, &global_umax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
  if (proc_id == 0)
  	printf ("umax = %e\n", global_umax);

  free( rows_in_proc );
  MPI_Finalize();

  return 0;
}

void prepare_send (double **grid, double *send_block, int grid_size, int *rows_in_proc, int proc_id, int proc_num)
{
	int i, j, k, l, rows, rows_before;

	rows_before = 0;
	for( i=0; i<proc_id; i++ ){
		rows_before += rows_in_proc[i];
	}
	rows = rows_in_proc[ proc_id ];
	for( j=0; j<grid_size; j++ ){
		for( i=0; i<rows; i++ ){
			send_block[ j*rows+i ] = grid[j][i+rows_before];
		}
	}
}

void transpose_block (double **grid, double *send_block, int grid_size, int *rows_in_proc, int proc_id, int proc_num)
{
	int i, j, k, l, rows, rows_before;

	rows_before = 0;
	for( i=0; i<proc_id; i++ ){
		rows_before += rows_in_proc[i];
	}
	rows = rows_in_proc[ proc_id ];
	for( j=0; j<grid_size; j++ ){
		for( i=0; i<rows; i++ ){
			send_block[ j*rows+i ] = grid[j][i+rows_before];
		}
	}
}

void parallel_transpose (double **transposed_grid, double **grid, int grid_size, int *rows_in_proc, int proc_id, int proc_num)
{
	double *send_block, *receive_block;
	int *send_sizes, *send_offset;
	int i, offset, j;
	// divide grid into equal sized blocks
	// transpose each block locally
	printf("1\n");
	send_block = malloc( sizeof(double)*rows_in_proc[ proc_id ]*grid_size );
	receive_block = malloc( sizeof(double)*rows_in_proc[ proc_id ]*grid_size );
	send_sizes = malloc( sizeof(int)*proc_num );
	send_offset = malloc( sizeof(int)*proc_num );

	printf("2\n");
	prepare_send( grid, send_block, grid_size, rows_in_proc, proc_id, proc_num );
	for( i=0; i<proc_num; i++ ){
		send_sizes[i] = grid_size*rows_in_proc[i];
	}
	printf("3\n");
	offset = 0;
	for( i=0; i<proc_num; i++ ){
		send_offset[i] = offset;
		offset += send_sizes[i];
	}

	printf("4\n");
	// all to all communication to send transposed blocks
	MPI_Alltoallv( send_block, send_sizes, send_offset, MPI_DOUBLE, \
				receive_block, send_sizes, send_offset, MPI_DOUBLE, MPI_COMM_WORLD);
	// merge to n/p x n matrix locally
	
	printf("5, %d\n", proc_id);
	// MPI_Alltoallv( *sendbuf, array of num of items sent to each process, 
	// index of first element in items sent to proc, sendtype,
	// *recvbuf, array of num of items max received from each process, 
	// result i position in recvbuffer, recvtype,
	// MPI_WORLD_COMM );
/*	
  for (j=0; j < grid_size; j++) {
    for (i=0; i < grid_size; i++) {
      transposed_grid[j][i] = grid[i][j];
    }
  }
*/
  free( send_block    );
  free( receive_block );
  free( send_sizes    );
  free( send_offset   );

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
