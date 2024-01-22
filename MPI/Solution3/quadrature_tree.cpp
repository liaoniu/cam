#include <iostream>
#include <iomanip>
#include <cmath>
#include <mpi.h>

#define STEPS 27720

int main(int argc, char *argv[])
{
  int rank, nproc, start, end, i;
  double total = 0.0;
  
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  start = rank * (STEPS / (double) nproc);
  end = (rank + 1) * (STEPS / (double) nproc);

  for(i = start; i < end; i++)
    {
      total += exp( -pow(10.0 * ((i + 0.5) / (STEPS + 1)), 2));
    }

  std::cout << "Total from process " << rank << " is " << total << std::endl;

  int iters = ceil(log2(nproc));
  for (int i = iters - 1; i >= 0; i--)
    {
      int offset = (1 << i);
      if (rank < offset && rank + offset < nproc)
	{
	  double tmp;
	  MPI_Recv(&tmp, 1, MPI_DOUBLE, rank + offset, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  total += tmp;
	}
      else if (rank >= offset && rank < 2 * offset)
	{
	  MPI_Send(&total, 1, MPI_DOUBLE, rank - offset, 1, MPI_COMM_WORLD);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }
  
  if (rank == 0)
    {
      std::cout << std::fixed << std::setprecision(16);
      std::cout << "Grand total = " << total << ", integral = " << 10.0 * total / STEPS << std::endl;
    }

  MPI_Finalize();

  return 0;
}
