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

int getColor(int iter, int totaliters, int index)
{
  // Returns r, g, or b color value given iter (result from calling function mandel)
  // The color map is a linear interpolation of the colors vector defined above
  
  //Gradient region
  int no_gradients = colors.size() - 1;
  double var = (double)no_gradients * (1.0 - (double)iter / (double)totaliters);
  int gr = (int) var;
  if (gr > no_gradients - 1) gr = no_gradients - 1;
  if (gr < 0) gr = 0;

  return colors[gr][index] + (int) round((colors[gr+1][index] - colors[gr][index]) * (var - (double)gr));
}

int main(){
  double bndr[4];
  int i, j, res, iters, rank, nproc;
  std::complex<double> z;
  std::ofstream img;

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

  int istart;
  int iend;
  int starts[nproc], ends[nproc], lengths[nproc], displs[nproc];
  
  if (rank == 0)
    {
      displs[0] = 0;
      for (i = 0; i < nproc; i++)
	{
	  starts[i] = (i * res) / nproc;
	  ends[i] = ((i+1) * res) / nproc;
	  lengths[i] = 3 * res * (ends[i] - starts[i]);
	  if (i > 0) displs[i] = displs[i-1] + lengths[i-1];
	}
    }
  
  MPI_Scatter(starts, 1, MPI_INT, &istart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(ends, 1, MPI_INT, &iend, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int row[res];
  unsigned char mylines[3 * (iend-istart) * res];
  unsigned char* alllines;

  if (rank == 0)
    {
      alllines = new unsigned char[3 * res * res];
    }

  for(i = istart; i < iend; i++){
    // Calculates mandel() for each row in the domain
    for(j = 0; j < res; j++){
      z.real(bndr[0] + j * ((bndr[1]-bndr[0]) / res));
      z.imag((bndr[3] - i * ((bndr[3]-bndr[2]) / res)));
      row[j] = mandel(z, iters);
    }

    // Creates colours
    for(j = 0; j < res; j++)
      {
	for (int index = 0; index < 3; index++)
	  {
	    mylines[3 * res * (i-istart) + 3 * j + index] = getColor(row[j], iters, index);
	  }
      }
  }

  MPI_Gatherv(mylines, (iend-istart) * 3 * res, MPI_CHAR, alllines, lengths, displs, MPI_CHAR, 0, MPI_COMM_WORLD);

  if (rank == 0)
    {
      // Open output file and write width and height of image
      img.open("mandel.pam", std::ios::binary);
      img << "P6\n" << res << " " << res << " 255\n";

      for (i = 0; i < res; i++)
	{
	  /* Writes row data to file. (3 * res is the total number 
	     of characters needed for colours in one row) */
	  img.write((char*)(alllines + 3 * res * i), 3 * res);
	}
        img.close();
	delete[] alllines;
    }
  
  MPI_Finalize();
  return 0;
}
