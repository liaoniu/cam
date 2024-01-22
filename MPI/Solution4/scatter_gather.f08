program scattergather

  use mpi_f08

  implicit none

  integer :: n, rank, i, j, sum
  integer, dimension(:), allocatable :: matrix, row, sums
  
  call mpi_init()
  call mpi_comm_size(mpi_comm_world, n)
  call mpi_comm_rank(mpi_comm_world, rank)

  if (rank .eq. 0) then
     allocate(matrix(n*n))
     print *, 'Set matrix on process 0'
     do j = 1, n
        do i = 1, n
           matrix(i + (j-1) * n) = (i-1) * n + (j-1)
           write(*,fmt="(i2, a)", advance="no") matrix(i + (j-1) * n), ' '
        end do
        print *, ''
     end do
  end if
  
  allocate(row(n))
  call mpi_scatter(matrix, n, mpi_integer, row, n, mpi_integer, 0, mpi_comm_world)

  write (*, fmt="(a, i1, a)") 'On process ', rank, ':'
  sum = 0
  do i = 1, n
     write (*,fmt="(i2, a)",advance="no") row(i), ' '
     sum = sum + row(i)
  end do
  print *
  
  write (*,fmt="(a,i1,a,i2)") 'Process ', rank, ' has sum = ', sum

  if (rank .eq. 0) then
     allocate(sums(n))
  end if

  call mpi_gather(sum, 1, mpi_integer, sums, 1, mpi_integer, 0, mpi_comm_world)

  if (rank .eq. 0) then
     write(*,fmt="(a)",advance="no") 'On rank 0 sums = '
     do i = 1, n
        write(*,fmt="(i2,a)",advance="no") sums(i), ' '
     end do
     print *
  end if

  call mpi_finalize()
  
end program scattergather
