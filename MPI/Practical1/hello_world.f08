program hello

  use mpi_f08

  integer :: size, rank

  call mpi_init()
  call mpi_comm_size(mpi_comm_world, size)
  call mpi_comm_rank(mpi_comm_world, rank)

  write(*, 101) rank, size
101 format('Hello parallel world, I am process ', I2, ' out of ', I2)

  call mpi_finalize()

end program hello
