program version

  use mpi_f08

  integer :: minor, major

  call mpi_get_version(major, minor)
    
  write(*, 101) major, minor
101 format('Run-time MPI version is ', I1, '.', I1)

end program version
