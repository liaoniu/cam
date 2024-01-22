#include<complex>
#include<iostream>
#include<fstream>
#include <vector>
#include <mpi.h>

std::vector<std::vector<int>> colors = {{68, 1, 84}, {32, 164, 134}, {68, 1, 84}, {57, 86, 140}, {31, 150, 139}, {115, 208, 85}, {253, 231, 37}, {42, 120, 142}, {255, 255, 255}};

int mandel(std::complex<double> z0, int iters)
{
  /* Repeats 320 iterations of the recurrence relation to 
     find the smallest integer, n, such that |z_n| < 2 
     
     If none is found, returns i = iters. This signifies 
     that z0 is not in Mandelbrot set */
  int i;
  std::complex<double> z;

  z = z0;
  for(i = 1; i < iters; i++){
    z = z*z + z0;
    if (abs(z) > 2.0)
      {
	break;
      }
  }
  return i;
}

int main(){
  double bndr[4];
  int i, j, res, iters, rank, nproc, r, s, msg;
  std::complex<double> z;
  MPI_Status st;
  
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == 0)
    {
      std::cout << "Enter resolution" << std::endl;
      std::cin >> res;
      do {
	std::cout << "Enter xmin" << std::endl;
	std::cin >> bndr[0];
	std::cout << "Enter xmax" << std::endl;
	std::cin >> bndr[1];
	std::cout << "Enter ymin" << std::endl;
	std::cin >> bndr[2];
	std::cout << "Enter ymax" << std::endl;
	std::cin >> bndr[3];
	if (bndr[0] >= bndr[1] || bndr[2] >= bndr[3])
	  {
	    std::cout << "Invalid domain boundaries." << std::endl;
	  }
      } while(bndr[0] >= bndr[1] || bndr[2] >= bndr[3]);

      std::cout << "Enter iterations" << std::endl;
      std::cin >> iters;
    }

  MPI_Bcast(&res, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bndr, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&iters, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int row[res];
  
  if (rank == 0)
    {
      int total_inside = 0;
      // Send first n jobs
      for (i = 1; i < nproc; i++)
	{
	  MPI_Send(&i, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}

      // Receive incoming results
      for (i = 0; i < res; i++)
	{
	  int num_inside;
	  MPI_Recv(&num_inside, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

	  r = st.MPI_TAG;
	  s = st.MPI_SOURCE;

	  total_inside += num_inside;
	  
	  if (i + nproc - 1 < res)
	    {
	      msg = i + nproc - 1;
	    }
	  else
	    {
	      msg = -1;
	    }

	  MPI_Send(&msg, 1, MPI_INT, s, 0, MPI_COMM_WORLD);
	}

      int pixels = res * res;
      double region_area = (bndr[3]-bndr[2]) * (bndr[1]- bndr[0]);
      std::cout << "Mandelbrot set area = " << total_inside / (double) pixels * region_area << std::endl;
    }
  else
    {
      for (; ;)
	{
	  MPI_Recv(&i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);

	  if (i == -1) {
	    break;
	  }

	  int num_inside = 0;
	  
	  for (j = 0; j < res; j++)
	    {
	      z = std::complex<double>(bndr[0] + j * (bndr[1]-bndr[0]) / res, bndr[3] - i * (bndr[3]-bndr[2])/res);
	      int a = mandel(z, iters);
	      if (a == iters) num_inside++;
	    }

	  MPI_Send(&num_inside, 1, MPI_INT, 0, i, MPI_COMM_WORLD);
	}
    }
  
  MPI_Finalize();
  return 0;
}
