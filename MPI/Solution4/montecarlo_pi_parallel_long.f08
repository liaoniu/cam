program montecarlo_pi_parallel
  use mpi_f08
  implicit none

  real(kind(1d0)) :: count = 10000000000.0
  integer :: nproc, rank, j
  real(kind(1d0)) :: i, hit = 0.0
  integer, allocatable :: seed(:)
  real(kind(1d0)) :: x(2), result, t, pi

  call mpi_init()
  call mpi_comm_size(mpi_comm_world, nproc)
  call mpi_comm_rank(mpi_comm_world, rank)

  call random_seed(size = j)
  allocate(seed(j))
  j = rank + 1
  call random_seed(put = seed)

  i = 1.0
  do while (i .le. count)
     call random_number(x)
     if ((x(1)**2 + x(2)**2) < 1d0) hit = hit + 1.0
     i = i + 1.0
  enddo

  t = hit / count

  call mpi_reduce(t, result, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world)

  pi = 3.1415926535897932384626433
  
  if (rank.eq.0) then
     write(*,101) result / nproc
     write(*, 102) pi/4.0, pi/4.0 - result / nproc
  end if
  
101 format('Ratio: ', f18.16)
102 format('Pi/4 = ', f18.16, t35, 'Error = ', f19.16)

  call mpi_finalize()
end program montecarlo_pi_parallel
