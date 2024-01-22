#include <fstream>
#include <math.h>
#include <mpi.h>

#define SIZE 400

using namespace std;

int main(){
  int i, j, n, rank, nproc, ht, niter, source, dest;
  double *g1, *g2, *p, *sendbuf, *recvbuf;
  ofstream img;
  MPI_Status st;

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ht = SIZE / nproc + 2;
  niter = 40000;

  g1 = new double[SIZE*ht];
  g2 = new double[SIZE*ht];
  
  for(i = 0; i < SIZE * ht ; i++) {
    g1[i] = 0.0; /* Zero our part of of grid */
  }
  
  if (rank == 0)
    {
      for(i = SIZE / 8; i < 7 * SIZE / 8; i++) {
	g1[i] = 1.0; /* Top bc */
      }
    }
  
  for(n = 0; n < niter; n++){
    /* Send our bottom row, which is neighbour’s top row */
    dest = rank + 1;
    source = rank - 1;

    if (rank == 0)
      {
	for(i = 0; i < SIZE; i++) {
	  g2[i] = g1[i]; //Top
	}
	source = MPI_PROC_NULL;
      }
    if (rank == nproc - 1)
      {
	for(i = 0; i < SIZE; i++) {
	  g2[i + SIZE * (ht - 1)] = g1[i + SIZE * (ht - 1)]; //Bottom
	}
	dest = MPI_PROC_NULL;
      }

    MPI_Sendrecv(g1 + (ht - 2) * SIZE, SIZE, MPI_DOUBLE, dest, 1, g1, SIZE, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
    
    /* Send our top row, which is neighbour's bottom row */
    dest = rank - 1;
    source = rank + 1;
    if (rank == 0)
      {
	dest = MPI_PROC_NULL;
      }
    else if (rank == nproc - 1)
      {
	source = MPI_PROC_NULL;
      }

    MPI_Sendrecv(g1 + SIZE, SIZE, MPI_DOUBLE, dest, 1, g1 + SIZE * (ht - 1), SIZE, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
    
    /* Bottom */

    for(j = 1;j < ht - 1; j++){

      /* Left */
      g2[j * SIZE] = (g1[j * SIZE + 1] +
		  g1[(j - 1) * SIZE] + g1[(j + 1) * SIZE] ) / 3.0;

      /* Body */
      for(i = 1; i < SIZE - 1; i++){
	g2[j * SIZE + i] = 0.25 * (g1[j * SIZE + i - 1] + g1[j * SIZE + i + 1] + g1[(j - 1) * SIZE + i] + g1[(j + 1) * SIZE + i]);
      }

      /* Right */
      g2[j * SIZE + SIZE - 1] = (g1[j * SIZE + SIZE - 2] +
			 g1[(j - 1) * SIZE + SIZE - 1] + g1[(j + 1) * SIZE + SIZE - 1]) / 3.0;
    }

    p = g1; g1 = g2; g2 = p;

  }

  if (rank == 0){

    img.open("laplace.pam");
    img << "P2\n" << SIZE << " " << SIZE << " 255\n";

    for(i = SIZE; i < (ht - 1) * SIZE; i++) {
      img << (int)(255*g1[i]) << "\n";
    }

    for(j = 1; j < nproc; j++){
      MPI_Recv(g1, ht * SIZE, MPI_DOUBLE,
	       j, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

      for(i = SIZE; i < (ht - 1) * SIZE; i++)
	{
	  img << (int)(255*g1[i]) << "\n";
	}
    }
    
    img.close();
  }
  else
    MPI_Send(g1, ht * SIZE, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

  MPI_Finalize();

  return 0;
}
