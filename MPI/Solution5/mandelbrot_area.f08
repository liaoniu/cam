function mand(z0,n)
  complex (kind=kind(1d0)) :: z,z0
  integer i,n,mand
  logical :: in = .false.

  z=z0
  do i=1,n-1
     z=z*z+z0
     if (abs(z).gt.2d0) then
        in = .true.
        exit
     endif
  enddo
  if (.not. in) then
     mand = n
  else
     mand = i
  end if

end function mand

program mandelbrot_area

  use mpi_f08
  
  implicit none

  complex (kind=kind(1d0)) :: z
  real (kind=kind(1d0)) :: x,y
  real :: var
  real (kind=kind(1d0)) :: bndr(4)
  integer :: ix,iy,res,i, j, a, iters, rank, nproc, r, s, msg, num_inside, total_inside
  integer mand
  type(mpi_status) :: st
  
  call mpi_init()
  call mpi_comm_size(mpi_comm_world, nproc)
  call mpi_comm_rank(mpi_comm_world, rank)
  
  if (rank .eq. 0) then 
     write(*,*)'Enter resolution'
     read(*,*) res

     do
        write(*,*)'Enter xmin'
        read(*,*) bndr(1)
        write(*,*)'Enter xmax'
        read(*,*) bndr(2)
        write(*,*)'Enter ymin'
        read(*,*) bndr(3)
        write(*,*)'Enter ymax'
        read(*,*) bndr(4)
        if ((bndr(1) .ge. bndr(2)) .or. (bndr(3) .ge. bndr(4))) then
           write(*,*)'Invalid domain boundaries'
        else
           exit
        end if
     end do

     write(*, *)'Enter iterations'
     read(*, *)iters
  end if

  call mpi_bcast(res, 1, mpi_integer, 0, mpi_comm_world)
  call mpi_bcast(bndr(1), 4, mpi_double_precision, 0, mpi_comm_world)
  call mpi_bcast(iters, 1, mpi_integer, 0, mpi_comm_world)

  if (rank .eq. 0) then

     total_inside = 0

     do i=1,nproc-1
        call mpi_send(i, 1, mpi_integer, i, 0, mpi_comm_world)
     end do

     do i = 0, res-1
        call mpi_recv(num_inside, 1, mpi_integer, mpi_any_source, mpi_any_tag, mpi_comm_world, st)

        r = st%mpi_tag
        s = st%mpi_source

        total_inside = total_inside + num_inside
        
        if (i + nproc - 1 .lt. res) then
           msg = i + nproc - 1
        else
           msg = -1
        end if

        call mpi_send(msg, 1, mpi_integer, s, 0, mpi_comm_world)
     end do

     print *, 'Mandelbrot set area = ', total_inside / ( 1.0 * res * res) * (bndr(4)-bndr(3)) * (bndr(2)-bndr(1))
     
  else
     do
        call mpi_recv(i, 1, mpi_integer, 0, 0, mpi_comm_world, st)

        if (i .eq. -1) then
           exit
        end if

        num_inside = 0
        
        y = bndr(3) + i * (bndr(4) - bndr(3)) / (res-1)
        
        do j=0, res-1
           x = bndr(1) + j * (bndr(2) - bndr(1)) / (res-1)
           a = mand(cmplx(x,y,kind(1d0)),iters)
           if (a .eq. iters) num_inside = num_inside + 1
        end do

        call mpi_send(num_inside, 1, mpi_integer, 0, i, mpi_comm_world)
     end do
  end if
  
  call mpi_finalize()
end program mandelbrot_area
