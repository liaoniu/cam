program findprimes_ms

  use mpi_f08

  implicit none
  
  integer :: rank, nproc, maxp, segment, start, finish, i, totalprimes, numprimes, worker, numsegs, nprimes
  integer, allocatable :: primes(:)
  logical :: isprime
  type(mpi_status) :: st
  
  call mpi_init()
  call mpi_comm_size(mpi_comm_world, nproc)
  call mpi_comm_rank(mpi_comm_world, rank)

  maxp = 10000000
  segment = 100

  if (nproc .lt. 2) then
     print *, 'Need at least 2 processors to run this'
     call abort
  end if
    
  allocate(primes(segment / 2))
  
  if (rank .eq. 0) then
     start = 0
     do i = 1, nproc - 1
        call mpi_send(start, 1, mpi_integer, i, 0, mpi_comm_world)
        start = start + segment
     end do

     totalprimes = 0
     finish = maxp + (nproc - 1) * segment
     do while (start < finish)
        call mpi_recv(primes, segment, mpi_integer, mpi_any_source, 1, mpi_comm_world, st)
        call mpi_get_count(st, mpi_integer, numprimes)
        totalprimes = totalprimes + numprimes
        worker = st%mpi_source
        if (start .lt. finish) then
           call mpi_send(start, 1, mpi_integer, worker, 0, mpi_comm_world)
        end if
           start = start + segment
     end do
     
     print *, 'There are ', totalprimes, ' primes less than ', maxp 
  else
     numsegs = 0
     do
        call mpi_recv(start, 1, mpi_integer, 0, 0, mpi_comm_world, st)
        if (start .ge. maxp) then
           exit
        end if

        nprimes = 0
        do i = start, start + segment - 1
           call checkprime(i, isprime)
           if (isprime) then
              primes(nprimes + 1) = i
              nprimes = nprimes + 1
           end if
        end do
        call mpi_send(primes, nprimes, mpi_integer, 0, 1, mpi_comm_world)
        numsegs = numsegs + 1
     end do
     print *, 'Process ', rank, ' did ', numsegs, ' chunks of work'
  end if
  
  call mpi_finalize()
end program findprimes_ms

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
