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

module pnm
contains
  subroutine ppm_write(set,unit,name,res)
    ! Writes data to file
    integer :: i,j,k,c,unit,res
    integer, allocatable :: set (:)
    character(*) :: name
    character, allocatable :: line(:)

    allocate(line(3*res))

    open(unit,file=name)
    ! NB First two bytes must be "P6"
    ! having an initial space will not do, so
    ! write(unit,*)'P3',res,res  is wrong 
    write(unit,'(''P6'')')
    write(unit,*) res,res
    write(unit,*) 255
    close(unit)

    open(unit,file=name,access='stream',position='append')
    do i=res-1,0,-1
       do j=0,res-1
          do k=0,2
             line(3*j+k+1)=char(set(3*res*i+3*j+k+1))
          enddo
       enddo
       write(unit) line
    enddo
    close(unit)

  end subroutine ppm_write
end module pnm

program mandel
  use pnm

  use mpi_f08
  
  implicit none

  complex (kind=kind(1d0)) :: z
  real (kind=kind(1d0)) :: x,y
  real :: var
  real (kind=kind(1d0)) :: bndr(4)
  integer :: ix,iy,res,i, j, k, iters, rank, nproc
  integer,  allocatable :: pallette(:,:), myset(:), allsets(:)
  integer :: colors(3, 9)
  integer mand
  
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

  allocate(pallette(3, 0:iters))
  
  ! Set color pallette using a linear interpolation of the following colors
  colors(:, 9) = (/ 68, 1, 84 /)
  colors(:, 8) = (/ 32, 164, 134 /)
  colors(:, 7) = (/ 68, 1, 84 /)
  colors(:, 6) = (/ 57, 86, 140 /)
  colors(:, 5) = (/ 31, 150, 139 /)
  colors(:, 4) = (/ 115, 208, 85 /)
  colors(:, 3) = (/ 253, 231, 37 /)
  colors(:, 2) = (/ 42, 120, 142 /)
  colors(:, 1) = (/ 255, 255, 255 /)
  
  pallette(:,0)=colors(:, 1)
  do i=1,8
     do j = 1 + ((i-1) * iters)/8, (i * iters)/8
        var = (1.0*j - (1.0 + (i-1.0) * iters/8.0)) / (iters/8.0)
        do k = 1, 3
           pallette(k, j) = colors(k, i) + nint((colors(k, i+1) - colors(k, i)) * var)
        end do
     end do
  end do

  allocate(myset(3 * (res/nproc + 1) * res))

  ! Call routine mand for each point in the domain 
  do iy=rank, res-1, nproc
     y = bndr(3) + iy * (bndr(4) - bndr(3)) / (res-1)
     do ix=0,res-1
        x = bndr(1) + ix * (bndr(2) - bndr(1)) / (res-1)
        do k=0,2
           myset(3 * res * (iy-rank)/nproc + 3 * ix + k + 1)=pallette(k+1,mand(cmplx(x,y,kind(1d0)),iters))
        end do
     enddo
  enddo

  do iy=0, res-1, nproc
     call mpi_gather(myset(3 * res * iy/nproc + 1), 3 * res, mpi_integer, allsets(3 * res * iy + 1), 3 * res, &
          mpi_integer, 0, mpi_comm_world)
  end do
  
  if (rank .eq. 0) then
     call ppm_write(allsets,10,"m.ppm", res)
  end if

  call mpi_finalize()
end program mandel
