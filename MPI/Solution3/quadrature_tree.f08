program quadrature_tree
  
  use mpi_f08

  implicit none

  integer, parameter :: steps = 27720
  integer :: nproc, rank, ierr, i, istart, iend, iters, offset
  real(kind(1d0)) :: total, tmp

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

iters = ceiling(log(nproc * 1.0)/log(2.0))

do i = iters - 1, 0, -1
   offset = ishft(1, i)
   if ((rank .lt. offset) .and. (rank + offset .lt. nproc)) then
      call mpi_recv(tmp, 1, mpi_double_precision, rank + offset, 1, mpi_comm_world, mpi_status_ignore)
      total = total + tmp
   else if ((rank .ge. offset) .and. (rank .lt. 2 * offset)) then
      call mpi_send(total, 1, mpi_double_precision, rank - offset, 1, mpi_comm_world)
   end if
   call mpi_barrier(mpi_comm_world)
end do


if (rank .eq. 0) then
   write(*,102) total, 10 * total / steps
102 format('Grand total is ',g21.16,' and integral is ',f18.16)
end if

call mpi_finalize()

end program quadrature_tree
