program quadrature
  
  implicit none

  integer, parameter :: steps = 27720
  integer :: i
  real(kind(1d0)) :: total

  total = 0.0
  
  do i = 0, steps
     total = total + exp(- (10d0 * (i + 0.5d0)/(steps + 1))**2)
  end do

write(*,102) total, 10 * total / steps
102 format('Grand total is ',g21.16,' and integral is ',f18.16)

end program quadrature
