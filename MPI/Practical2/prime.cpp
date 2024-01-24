#include <iostream>
#include <vector>
#include <mpi.h>

using namespace std;

bool isPrime(int n)
{
    if (n == 0)
        return false;
    if (n == 1)
        return false;
    if (n == 2)
        return true;
    
    else
    {
        for (int i = 2; i <= sqrt(n); i++)
        {
            if (n % i == 0)
                return false;
        }
    }
    return true;
    
}




int main()
{
    int N = 10000000;
    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Status st;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int n = N/size;
    int i_start = rank*n;
    int i_end = ((rank+1)*n);
    int* local_primes = new int[(i_end - i_start)/2];
    int local_nprimes = 0;
    for (int i = i_start; i != i_end; i++)
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
        int *primes = new int[N/2];
        int curr_nprimes = local_nprimes;
        for (int i = 0; i != local_nprimes; i++)
        {
            primes[i] = local_primes[i];
        }

        for (int i = 1; i != size; i++)
        {
            MPI_Recv(primes + curr_nprimes, N/size + 1, MPI_INT, i, i, MPI_COMM_WORLD, &st);
            int nRecv;
            MPI_Get_count(&st, MPI_INT, &nRecv);
            curr_nprimes += nRecv;
        }


        cout << "Total " << curr_nprimes << " primes < 10000000" << endl;
    }
    
    MPI_Finalize();
    return 0;


}