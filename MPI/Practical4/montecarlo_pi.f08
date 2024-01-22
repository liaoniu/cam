program montecarlo_pi

  integer, parameter :: count=100000
  integer :: i, hit=0
  real(kind(1d0)) :: x(2), pi

  do i=1, count
     call random_number(x)
     if ((x(1)**2+x(2)**2)<1d0) hit = hit + 1
  enddo

  write(*,101) hit / real(count,kind(1d0))

101 format ('Ratio: ', f18.16)

  pi = 3.1415926535897932384626433

  write(*, 102) pi/4.0, pi/4.0 - hit / real(count,kind(1d0))
  102 format('Pi/4 = ', f18.16, t35, 'Error = ', f19.16)
  
end program montecarlo_pi
