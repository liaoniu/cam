#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<cmath>
#include<mpi.h>

int main(int argc, char **argv){
  int nproc, rank;
  long int i, count, hit = 0;
  double x, y, t, result;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (rank == 0)
    {
      count = atoi(argv[1]);
    }
  MPI_Bcast(&count, 1, MPI_LONG, 0, MPI_COMM_WORLD);

  srand(rank + 1);
  
  for(i = 0; i < count; i++)
    {
      x = rand() / (double) RAND_MAX;
      y = rand() / (double) RAND_MAX;
      if (x*x + y*y < 1.0)
	{
	  hit++;
	}
    }

  t = hit / (double) count;

  MPI_Reduce(&t, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank==0)
    {
      std::cout << std::fixed << std::setprecision(16);
      std::cout << "Ratio: " << result / nproc << std::endl;
      std::cout << "Pi/4 = " << M_PI / 4.0 << "\tError = " << M_PI / 4.0 - result / nproc << std::endl;
    }
  
  MPI_Finalize();
  
  return 0;
}
