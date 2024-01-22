program deadlock

  use mpi_f08

  integer :: rank, size
  character(len=:), allocatable :: msg
  
  call mpi_init()

  size=10240
  
  allocate(character(len=10240)::msg)
  
  call mpi_comm_rank(mpi_comm_world, rank)

  if (rank .eq. 0) then
     call mpi_send(msg, size, mpi_character, 1, 1, mpi_comm_world);
     call mpi_recv(msg, size, mpi_character, 1, 1, mpi_comm_world, mpi_status_ignore);
  end if
  if (rank .eq. 1) then
     call mpi_send(msg, size, mpi_character, 0, 1, mpi_comm_world);
     call mpi_recv(msg, size, mpi_character, 0, 1, mpi_comm_world, mpi_status_ignore);
  end if
  
  
  write(*, 101) rank
101 format('Rank ', I2, ' finished')

  call mpi_finalize()

end program deadlock
