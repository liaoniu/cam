#include <iostream>
#include <iomanip>
#include <math.h>
#include <mpi.h>

#define STEPS 27720

int main(int argc, char *argv[])
{
  int rank, nproc, start, end, i;
  double total = 0, result;
  
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

  if (rank == 0)
    {
      result = total;
      for (int i = 1; i < nproc; i++)
	{
	  MPI_Recv(&total, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  result += total;
	}

      std::cout << std::fixed << std::setprecision(16);
      std::cout << "Grand total = " << result << ", integral = " << 10.0 * result / STEPS << std::endl;
    }
  else
    {
      MPI_Send(&total, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

  MPI_Finalize();

  return 0;
}
