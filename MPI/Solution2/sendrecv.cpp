#include <iostream>
#include <mpi.h>

int main(void)
{
  int rank,nproc;
  MPI_Status st;
  
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int v = 0;

  if(rank == 0)
    {
      MPI_Send(&v, 1, MPI_INT, (rank+1) % nproc, 1, MPI_COMM_WORLD);
      MPI_Recv(&v, 1, MPI_INT, nproc-1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
    }
  else
    {
      MPI_Recv(&v, 1, MPI_INT, rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
      v += rank;
      std::cout << "Hit rank " << rank << " v = " << v << std::endl;
      MPI_Send(&v, 1, MPI_INT, (rank+1) % nproc, 1, MPI_COMM_WORLD);
    }

  if(rank == 0)
    {
      std::cout << "Final v = " << v << ", should be = " << (nproc*(nproc-1))/2 << std::endl;
    }
  
  MPI_Finalize();
}
