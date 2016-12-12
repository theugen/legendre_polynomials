include 'legendre_routines.f90'
program legendre
  use, intrinsic :: iso_fortran_env !This is for using real128 
  implicit none
  integer, parameter :: qp = REAL128 !quadruple precision
  integer :: l=4, m=0, i, j
  real(kind=qp) :: x=sin(1.57d0)
  real(kind=qp), allocatable, dimension(:,:) :: ans

  print *, 'sin(pi/2)', x
  allocate(ans(l+1, l+1))
  call alegendre_increment_l(x, l, m, ans)
  call alegendre_increment_m(x, l, ans)

  do i=1, l
      print *, (ans(i,j), j=1,l+1)
  end do

end program
