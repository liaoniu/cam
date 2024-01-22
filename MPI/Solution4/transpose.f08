program transpose

  use mpi_f08

  implicit none

  integer :: n, rank, i, j
  integer, dimension(:), allocatable :: matrix, matrixT, row
  logical :: correct
  
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
     print *
  end if
  
  allocate(row(n))

  do i = 1, n
     call mpi_scatter(matrix(1 + (i-1) * n), 1, mpi_integer, row(1 + (i-1)), 1, mpi_integer, 0, mpi_comm_world)
  end do

  if (rank .eq. 0) then
     allocate(matrixT(n*n))
  end if

  call mpi_gather(row, n, mpi_integer, matrixT, n, mpi_integer, 0, mpi_comm_world)
  
  if (rank .eq. 0) then
     correct = .true.
     
     do j = 1, n
        do i = 1, n
           write(*, fmt="(i2, a)",advance="no") matrixT(i + (j-1)*n), ' '
           if (matrixT(i + (j-1)*n) .ne. (j-1) * n + (i-1)) then
              correct = .false.
           end if
           
        end do
        print *
     end do

     if (correct) then
        print *, 'Matrix is correctly transposed!'
     end if
  end if
  
  call mpi_finalize()
  
end program transpose
