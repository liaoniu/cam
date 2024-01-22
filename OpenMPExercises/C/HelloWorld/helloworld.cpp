#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


int main(){
#pragma omp parallel
    {
#pragma omp single
    {
    printf("num of threads = %d\n",omp_get_num_threads());
    }
#pragma omp critical 
     {
      printf("hello from thread %d\n",omp_get_thread_num());
     }
    }
}
