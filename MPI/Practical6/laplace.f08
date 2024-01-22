program laplace
  real (kind(1d0)), pointer, dimension(:,:) :: g1, g2, p
  integer :: i, j, k, img, ht, n, niter
  integer, parameter :: SIZE=400

  ht = SIZE + 2
  niter = 40000

  allocate (g1(SIZE, ht), g2(SIZE, ht))

  g1 = 0.0
  g1(SIZE / 8 + 1 : 7 * SIZE / 8, 1) = 1.0

  ! Top
  g2(:, 1) = g1(:, 1)
  ! Bottom
  g2(:, ht) = g1(:, ht)

  do n = 1, niter

     do j = 2, ht - 1
        ! Left
        g2(1, j) = (g1(1, j - 1) + g1(1, j + 1) + g1(2, j)) / 3
        ! Body
        do i = 2, SIZE - 1
           g2(i, j) = 0.25 * (g1(i, j - 1) + g1(i, j + 1) + g1(i - 1, j) + g1(i + 1, j))
        enddo
        ! Right
        g2(SIZE, j) = (g1(SIZE, j - 1) + g1(SIZE, j + 1) + g1(SIZE - 1, j)) / 3
     enddo

     p => g1
     g1 => g2
     g2 => p

  enddo

  img = 10
  open(unit = img, file = 'laplace.pgm')
  write(img,'(a2)')"P2"
  write(img,'(I0,x,I0,x,I0)')SIZE, SIZE, 255
  do j = 2, ht - 1
     do i = 1, SIZE
        write(img,'(I0)')int(255 * g1(i, j))
     enddo
  enddo
end program laplace
