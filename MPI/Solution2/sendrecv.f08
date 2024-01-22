program sendrecv

  use mpi_f08

  implicit none

  integer :: v, rank, size

  call mpi_init()
  call mpi_comm_size(mpi_comm_world, size)
  call mpi_comm_rank(mpi_comm_world, rank)

  v = 0

  if (rank .eq. 0) then
     call mpi_send(v, 1, mpi_integer, modulo(rank+1, size), 1, mpi_comm_world)
     call mpi_recv(v, 1, mpi_integer, size - 1, mpi_any_tag, mpi_comm_world, mpi_status_ignore)
  else
     call mpi_recv(v, 1, mpi_integer, rank - 1, mpi_any_tag, mpi_comm_world, mpi_status_ignore)
     v = v + rank
     write(*, 101) rank, v
     101 format('Hit rank ', I2, ' with v = ', I3)
     call mpi_send(v, 1, mpi_integer, modulo(rank+1,size), 1, mpi_comm_world)
 end if

 if (rank .eq. 0) then
    write (*, 102) v, size*(size-1)/2
    102 format ('Final v = ', I3, ', should be = ', I3)
    end if

    call mpi_finalize()
end program sendrecv
