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

  MPI_Reduce(&total, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0)
    {
      std::cout << std::fixed << std::setprecision(16);
      std::cout << "Grand total = " << result << ", integral = " << 10.0 * result / STEPS << std::endl;
    }

  MPI_Finalize();

  return 0;
}