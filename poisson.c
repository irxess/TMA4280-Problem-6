/*
  C-program to solve the two-dimensional Poisson equation on 
  a unit square using one-dimensional eigenvalue decompositions
  and fast sine transforms

  einar m. ronquist
  ntnu, october 2000
  revised, october 2001
*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>

#include <mpi.h>

/* function prototypes */
double *createdoubleArray (int n);
double **createdouble2DArray (int m, int n);
void freedouble2DArray (double **array);

void block_transpose (double **bt, double **b, int m, int m_loc, int size, int rank);
void block_copy (double **from, double **to, int m, int offset);
void transpose (double **bt, double **b, int m, int offset);
void fst_(double *v, int *n, double *w, int *nn);
void fstinv_(double *v, int *n, double *w, int *nn);

/* global variables <3 */
static int m_loc, m_padded;
static MPI_Datatype block_type;
static double *sendbuf, *recvbuf;
static double pi;

double f(double x, double y) 
{
  return 5 * pi * pi * sin (pi*x) * sin (2*pi*y);
}

int main(int argc, char **argv )
{
  double *diag, **b, **bt, **z;
  double h, umax;
  double starttime;
  int i, j, k;
  int n, m, nn;

  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */

  if( argc < 2 ) {
    printf("need a problem size\n");
    return 1;
  }

  n  = atoi(argv[1]);
  m  = n-1;
  nn = 4*n;

  int size, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  starttime = MPI_Wtime();

  int m_loc = m / size;
  if (m % size) {
    m_loc++;
    m_padded = m_loc * size;
  } else {
    m_padded = m;
  }

  sendbuf = (double *) malloc (m_padded*m_loc*sizeof(double));
  recvbuf = (double *) malloc (m_padded*m_loc*sizeof(double));

  MPI_Type_vector(m_loc, m_loc,
		  m_padded, MPI_DOUBLE, &block_type);
  MPI_Type_commit(&block_type);

  diag = createdoubleArray (m);
  b    = createdouble2DArray (m_loc,m_padded);
  bt   = createdouble2DArray (m_loc,m_padded);
  z    = createdouble2DArray (omp_get_max_threads(), nn);

  h    = 1./(double)n;
  pi   = 4.*atan(1.);

  for (i=0; i < m; i++) {
    diag[i] = 2.*(1.-cos((i+1)*pi/(double)n));
  }
  for (j=0; (rank*m_loc + j < m) && (j < m_loc); j++) {
    for (i=0; i < m; i++) {
      b[j][i] = h*h*f(rank*m_loc*h + j*h, i*h);
    }
  }
#pragma omp parallel for
  for (j=0; j < m_loc; j++) {
    fst_(b[j], &n, z[omp_get_thread_num()], &nn);
  }

  block_transpose (bt, b, m_padded, m_loc, size, rank);

#pragma omp parallel for
  for (i=0; i < m_loc; i++) {
    fstinv_(bt[i], &n, z[omp_get_thread_num()], &nn);
  }
  
  for (j=0; j < m_loc; j++) {
    for (i=0; i < m; i++) {
      bt[j][i] = bt[j][i]/(diag[i]+diag[j + rank*m_loc]);
    }
  }
#pragma omp parallel for
  for (i=0; i < m_loc; i++) {
    fst_(bt[i], &n, z[omp_get_thread_num()], &nn);
  }

  block_transpose (b, bt, m_padded, m_loc, size, rank);
#pragma omp parallel for
  for (j=0; j < m_loc; j++) {
    fstinv_(b[j], &n, z[omp_get_thread_num()], &nn);
  }

  umax = 0.0;
  for (j=0; j < m_loc; j++) {
    for (i=0; i < m; i++) {
      double uxy = sin(pi*(rank*m_loc*h + j*h)) * sin(2*pi*i*h);
      double err = fabs(b[j][i] - uxy);
      if (err > umax) umax = err;
    }

  }

  double global_err;
  MPI_Reduce(&umax, &global_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    printf (" Highest error was %e.\n", global_err);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    printf ("Elapsed time: %f secs\n", MPI_Wtime() - starttime);
  }

  MPI_Finalize();

  return 0;
}

void block_transpose (double **bt, double **b, int m, int m_loc, int size, int rank)
{
  int i,j, pos = 0;
  for (i=0; i < size; i++) {
    MPI_Pack(b[0]+i*m_loc, 1, block_type, sendbuf, m*m_loc*sizeof(double), &pos, MPI_COMM_WORLD);
  }
  MPI_Alltoall(sendbuf, m_loc*m_loc, MPI_DOUBLE, recvbuf, m_loc*m_loc, MPI_DOUBLE, MPI_COMM_WORLD);
  pos = 0;
  for (i=0; i < size; i++) {
    MPI_Unpack(recvbuf, m*m_loc*sizeof(double), &pos, bt[0]+i*m_loc, 1, block_type, MPI_COMM_WORLD);
  }
  double **temp_block = createdouble2DArray (m_loc, m_loc);
  for (i = 0; i < size; i++) {
    transpose(temp_block, bt, m_loc, m_loc*i);
    block_copy(temp_block, bt, m_loc, m_loc*i);
  }
  freedouble2DArray (temp_block);
}

void block_copy (double **from, double **to, int m, int offset)
{
  int i,j;
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      to[j][i+offset] = from[j][i];
    }
  }
}

void transpose (double **bt, double **b, int m, int offset)
{
  int i, j;
#pragma omp parallel for private(i)
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      bt[j][i] = b[i][j+offset];
    }
  }
}

double *createdoubleArray (int n)
{
  double *a;
  int i;
  a = (double *)malloc(n*sizeof(double));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

double **createdouble2DArray (int n1, int n2)
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

void freedouble2DArray (double **array)
{
  free(array[0]);
  free(array);
}
