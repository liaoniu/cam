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
  int segment = 100;

  if(nproc < 2)
    {
      std::cout << "Need at least 2 processors to run this" << std::endl;
      std::abort();
    }

  int* primes = new int[segment / 2];

  if(rank == 0)
    {
      int start = 0;
      for(int i=1 ; i < nproc ; i++)
	{
	  MPI_Send(&start, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	  start += segment;
	}

      int totalPrimes = 0;
      int finish = maxP + (nproc - 1) * segment;
      for( ; start < finish ; start += segment )
	{
	  MPI_Recv(primes, segment, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &st);
	  int numPrimes;
	  MPI_Get_count(&st, MPI_INT, &numPrimes);
	  totalPrimes += numPrimes;
	  int worker = st.MPI_SOURCE;
	  if(start < finish)
	    MPI_Send(&start, 1, MPI_INT, worker, 0, MPI_COMM_WORLD);
	}

      std::cout << "There are " << totalPrimes << " primes less than " << maxP << std::endl;
    }
  else
    {
      int numSegs = 0;
      int start;
      for( ; ; )
	{
	  MPI_Recv(&start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
	  if(start >= maxP) break;
      
	  int nprimes = 0;
	  for(int i=start ; i < start + segment ; i++)
	    {
	      if(isPrime(i))
		{
		  primes[nprimes++] = i;
		}
	    }
	  MPI_Send(primes, nprimes, MPI_INT, 0, 1, MPI_COMM_WORLD);
	  numSegs++;
	}
      std::cout << "Process " << rank << " did " << numSegs << " chunks of work" << std::endl;
    }
  
  MPI_Finalize();
}
