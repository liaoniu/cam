program quadrature_reduce
  
  use mpi_f08

  implicit none

  integer, parameter :: steps = 27720
  integer :: nproc, rank, ierr, i, istart, iend
  real(kind(1d0)) :: total, result

  call mpi_init()

  call mpi_comm_size(mpi_comm_world, nproc)
  call mpi_comm_rank(mpi_comm_world, rank)

  total = 0.0
  istart = rank * (steps / real(nproc,kind(1d0)))
  iend = (rank + 1) * (steps / real(nproc,kind(1d0)))
  do i = istart, iend - 1
     total = total + exp(- (10d0 * (i + 0.5d0)/(steps + 1))**2)
  end do
write(*,101) rank, total
101 format('Total from process ',i3,' is ',g16.8)

call mpi_reduce(total,result,1,mpi_double_precision,mpi_sum,0, &
     mpi_comm_world)

if (rank .eq. 0) write(*,102) result, 10 * result / steps

102 format('Grand total is ',g21.16,' and integral is ',f18.16)

call mpi_finalize()

end program quadrature_reduce
