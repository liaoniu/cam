#include <mpi.h>
#include <iostream>

int main(void)
{
  int N, rank;
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&N);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int* matrix;
  if(rank == 0)
    {
      matrix = new int[N*N];
      std::cout << "Set matrix on process " << rank << std::endl;
      for(int j=0 ; j < N ; j++)
	{
	  for(int i=0 ; i < N ; i++)
	    {
	      matrix[i + j*N] = i*N + j;
	      std::cout << matrix[i + j*N] << " ";
	    }
	  std::cout << std::endl;
	}
    }

  int* row = new int[N];
  MPI_Scatter(matrix, N, MPI_INT, row, N, MPI_INT, 0, MPI_COMM_WORLD);

  std::cout << "On process " << rank << ":";
  int sum = 0;
  for(int i=0 ; i < N ; i++)
    {
      std::cout << row[i] << " ";
      sum += row[i];
    }

  std::cout << std::endl << "Process " << rank << " has sum = " << sum << std::endl;

  int* sums;
  if(rank == 0)
    {
      sums = new int[N];
    }

  MPI_Gather(&sum, 1, MPI_INT, sums, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(rank == 0)
    {
      std::cout << "On rank 0 sums = ";
      for(int i=0 ; i < N ; i++)
	{
	  std::cout << sums[i] << " ";
	}
      std::cout << std::endl;
    }

  MPI_Finalize();
}
