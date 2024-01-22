#include <mpi.h>
#include <iostream>

int main(void)
{
  int N, rank;

  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&N);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int* matrix = NULL;
  if(rank == 0)
    {
      matrix = new int[N*N];
      std::cout << "Set matrix on process " << rank << std::endl;
      for(int i=0 ; i < N ; i++)
	{
	  for(int j=0 ; j < N ; j++)
	    {
	      matrix[i + j*N] = i*N + j;
	      std::cout << matrix[i + j*N] << " ";
	    }
	  std::cout << std::endl;
	}
      std::cout << std::endl;
    }

  int* row = new int[N];
  int* column = new int[N];

  MPI_Scatter(matrix, N, MPI_INT, row, N, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Alltoall(MPI_IN_PLACE, 1, MPI_INT, row, 1, MPI_INT, MPI_COMM_WORLD);
  
  int* matrixT = NULL;
  if(rank == 0)
    {
      matrixT = new int[N*N];
    }
  
  MPI_Gather(row, N, MPI_INT, matrixT, N, MPI_INT, 0, MPI_COMM_WORLD);

  if(rank == 0)
    {
      // Now check for correctness
      bool correct = true;
      for(int i=0 ; i < N ; i++)
	{
	  for(int j=0 ; j < N ; j++)
	    {
	      std::cout << matrixT[i + j*N] << " ";
	      if( matrixT[i + j*N] != j*N + i )
		{
		  correct = false;
		}
	    }
	  std::cout << std::endl;
	}

      if(correct)
	{
	  std::cout << "Matrix is correctly transposed!" << std::endl;
	}
    }
  
  delete[] matrixT;
  delete[] row;
  
  MPI_Finalize();
}
