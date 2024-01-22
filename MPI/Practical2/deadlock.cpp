#include <iostream>
#include <mpi.h>

#define SIZE 10240

int main(){
  int rank;
  void *ptr;
	
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
  ptr=malloc(SIZE);
	
  if (rank == 0){
    MPI_Send(ptr, SIZE, MPI_CHAR, 1, 1, MPI_COMM_WORLD);
    MPI_Recv(ptr, SIZE, MPI_CHAR, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (rank == 1){
    MPI_Send(ptr, SIZE, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
    MPI_Recv(ptr, SIZE, MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
	
  std::cout << "Rank " << rank << "finished" << std::endl;
  MPI_Finalize();
  return 0;
}
