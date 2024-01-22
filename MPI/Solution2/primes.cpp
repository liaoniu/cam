#include <mpi.h>
#include <cmath>

bool isPrime(int n)
{
  if(n == 2) return true;
  if(n == 1) return false;
  if(n % 2 == 0)
    {
      return false;
    }
  else
    {
      for(int t=2 ; t <= sqrt(n) ; t++)
	{
	  if( n % t == 0 )
	    {
	      return false;
	    }
	}
      return true;
    }
}


int main(void)
{
  int rank, nproc;
  MPI_Status st;
  
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int maxP = 10000000;
  int local_start = maxP * rank / nproc;
  int local_end = maxP * (rank + 1) / nproc;
  int local_nprimes = 0;
  int* local_primes = new int[(local_end - local_start)/2];
  
  for (int i = local_start; i < local_end; i++)
    {
      if (isPrime(i))
	{
	  local_primes[local_nprimes] = i;
	  local_nprimes++;
	}
    }

  if (rank != 0)
    {
      MPI_Send(local_primes, local_nprimes, MPI_INT, 0, rank, MPI_COMM_WORLD);
    }
  else
    {
      int totalprimes = local_nprimes;
      int* primes = new int[maxP/2];
      for (int p = 0; p < local_nprimes; p++)
	{
	  primes[p] = local_primes[p];
	}

      std::cout << local_nprimes << " primes computed in rank 0" << std::endl;
      
      for (int i = 1; i < nproc; i++)
	{
	  MPI_Recv(primes + totalprimes, maxP/nproc + 1, MPI_INT, i, i, MPI_COMM_WORLD, &st);
	  int primesRecvd;
	  MPI_Get_count(&st, MPI_INT, &primesRecvd);
	  totalprimes += primesRecvd;
	  std::cout << primesRecvd << " primes received from rank " << i << std::endl;
	}
      std::cout << "There are " << totalprimes << " primes less than " << maxP << std::endl;
    }
 
  MPI_Finalize();
}
