program laplace_send

  use mpi_f08

  real (kind(1d0)), pointer, dimension(:,:) :: g1, g2, p
  integer :: i, j, k, img, ht, n, niter, rank, nproc
  integer, parameter :: size = 400
  type(mpi_status) :: st
  
  call mpi_init()
  call mpi_comm_size(mpi_comm_world, nproc)
  call mpi_comm_rank(mpi_comm_world, rank)
  
  ht = size/nproc + 2
  niter = 40000

  allocate (g1(size, ht), g2(size, ht))

  g1 = 0.0
  if (rank .eq. 0) then
     g1(size / 8 + 1 : 7 * size / 8, 1) = 1.0
  end if
  
  ! Top
  g2(:, 1) = g1(:, 1)
  ! Bottom
  g2(:, ht) = g1(:, ht)

  do n = 1, niter
     if (rank .ne. 0) then
        call mpi_recv(g1(:, 1), size, mpi_double_precision, rank - 1, mpi_any_tag, mpi_comm_world, st)
     end if
     
     if (rank .ne. nproc - 1) then
        call mpi_send(g1(:, ht - 1), size, mpi_double_precision, rank + 1, 1, mpi_comm_world)
        call mpi_recv(g1(:, ht), size, mpi_double_precision, rank + 1, mpi_any_tag, mpi_comm_world, st)
     end if

     if (rank .ne. 0) then
        call mpi_send(g1(:, 2), size, mpi_double_precision, rank - 1, 1, mpi_comm_world)
     end if

     do j = 2, ht - 1
        ! Left
        g2(1, j) = (g1(1, j - 1) + g1(1, j + 1) + g1(2, j)) / 3
        ! Body
        do i = 2, size - 1
           g2(i, j) = 0.25 * (g1(i, j - 1) + g1(i, j + 1) + g1(i - 1, j) + g1(i + 1, j))
        enddo
        ! Right
        g2(size, j) = (g1(size, j - 1) + g1(size, j + 1) + g1(size - 1, j)) / 3
     enddo

     p => g1
     g1 => g2
     g2 => p

  enddo

  if (rank .eq. 0) then 
     img = 10
     open(unit = img, file = 'laplace.pgm')
     write(img,'(a2)')"P2"
     write(img,'(I0,x,I0,x,I0)')size, size, 255
     do j = 2, ht - 1
        do i = 1, size
           write(img,'(I0)')int(255 * g1(i, j))
        enddo
     enddo
     do k = 1, nproc - 1
        call mpi_recv(g1, ht * size, mpi_double_precision, k, mpi_any_tag, mpi_comm_world, st)
        do j = 2, ht - 1
           do i = 1, size
              write(img,'(I0)')int(255 * g1(i, j))
           end do
        end do
     end do
  else
     call mpi_send(g1, ht * size, mpi_double_precision, 0, 1, mpi_comm_world)
  end if

call mpi_finalize()
end program laplace_send
