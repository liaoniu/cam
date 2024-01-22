program findprimes

  use mpi_f08

  implicit none
  
  integer :: rank, nproc, maxp, local_start, local_end, local_nprimes, totalprimes, p, primesrecvd, i
  integer, allocatable :: local_primes(:)
  integer, allocatable :: primes(:)
  logical :: isprime
  type(mpi_status) :: st
  
  call mpi_init()
  call mpi_comm_size(mpi_comm_world, nproc)
  call mpi_comm_rank(mpi_comm_world, rank)

  maxp = 10000000
  local_start = maxp * rank / nproc + 1
  local_end = maxp * (rank + 1) / nproc + 1
  local_nprimes = 0
  allocate(local_primes((local_end - local_start)/2))

  do i = local_start, local_end - 1
     call checkprime(i, isprime)
     if (isprime) then
        local_primes(local_nprimes + 1) = i
        local_nprimes = local_nprimes + 1
     end if
  end do
  
  if (rank .eq. 0) then
     totalprimes = local_nprimes
     allocate(primes(maxp/2))
     primes(1:local_nprimes) = local_primes(1:local_nprimes)

     print *, local_nprimes, ' primes computed in rank 0'
     
     do i = 1, nproc-1
        call mpi_recv(primes(totalprimes + 1), maxp/nproc + 1, mpi_integer, i, i, mpi_comm_world, st)
        call mpi_get_count(st, mpi_integer, primesrecvd)
        print *, primesrecvd, ' primes received from rank ', i
        totalprimes = totalprimes + primesrecvd
     end do

     print *, 'There are ', totalprimes, ' primes less than ', maxp 
  else
     call mpi_send(local_primes(1), local_nprimes, mpi_integer, 0, rank, mpi_comm_world)
  end if
  
  call mpi_finalize()
end program findprimes

subroutine checkprime(n, isprime)
  integer n, t
  logical isprime
  if (n .eq. 2) then
     isprime = .true.
     return 
  else if (n .eq. 1) then
     isprime = .false.
     return 
  else if (modulo(n, 2) .eq. 0) then
     isprime = .false.
     return 
  else
     t = 2
     do while (t * t .le. n)
        if (modulo(n, t) .eq. 0) then
           isprime = .false.
           return
        end if
        t = t + 1
     end do
     isprime = .true.
  end if
  return
end subroutine checkprime
