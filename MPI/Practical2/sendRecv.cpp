#include <iostream>
#include <mpi.h>

int main()
{
    int v = 0;
    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank == 0)
    {
        MPI_Send(&v, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD);
        MPI_Recv(&v, 1, MPI_INT, size-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    else
    {
        MPI_Recv(&v, 1, MPI_INT, (rank-1), MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        v += rank;
        MPI_Send(&v, 1, MPI_INT, (rank + 1)%size, 1, MPI_COMM_WORLD);
    }

    if (rank == 0)
    {
        std::cout << v << std::endl;
    }
    MPI_Finalize();
    return 0;

}