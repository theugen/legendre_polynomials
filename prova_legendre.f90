include 'legendre_routines.f90'

!Function to compute the gravitational potential, given
!the gravitational constant G, Earth's mean radius a_e,
!a set of associated Legendre polynomials lpoly and relatives
!Stokes coefficient stokes_coef.
!lmax is maximum degree for which associated Legendre polynomials
!contained in lpoly have been computed, while length is the length
!of GRACE GSM file considered, stokes_coef.
function compute_v(G, a_e, lpoly, lmax, stokes_coef, length) result(v_grav)
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, parameter :: qp = REAL128
  real(kind=qp), intent(in):: G, a_e
  integer :: lmax, length
  real(kind=qp), dimension(ceiling((real(lmax**2 + 4*lmax + 1))/(2)),&
                 &3), intent(in)&
                 & :: lpoly 
  real(kind=qp), dimension(length, 3), intent(in) :: stokes_coef
  real(kind=qp) :: v_grav



end function


!Program to process and prepare for plotting GRACE level 2
!data.
program legendre
  use, intrinsic :: iso_fortran_env !This is for using real128 
  implicit none
  integer, parameter :: qp = REAL128 !quadruple precision
  integer :: l=4, m=0, i, j
  real(kind=qp) :: x=sin(1.57d0)
  real(kind=qp), allocatable, dimension(:,:) :: ans
  real(kind=qp), allocatable, dimension(:,:) :: normalized
  integer :: nrow
  !Declaring gravitational constant and Earth's mean radius
  real(kind=qp), parameter :: G=6.674d-11, a_e = 6.371d6  
  !Declaring Earth's Potential
  real(kind=qp) :: v_grav
  !Declaring return type of compute_v
  real(kind=qp) :: compute_v
  integer :: in_unit = 1 !Input device
  !Computing # of row of result
  nrow = ceiling((real(l**2 + 4*l + 1))/(2))

  
  !Reading GRACE file!
  do
    read(in_unit, *, iostat=stat)

  end do 

  print *, 'sin(pi/2)', x
  allocate(ans(l+1, l+1))
  allocate(normalized(nrow, 3))
  call alegendre_increment_l(x, l, m, ans)
  call alegendre_increment_m(x, l, ans)
  call normalize_alegendre(ans, l, normalized)

  do i=1, l
      print *, (ans(i,j), j=1,l+1)
  end do

  do i=1, nrow
    print *, normalized(i, 1), ' ', normalized(i,2)&
            &, ' ', normalized(i,3)  
  end do


end program
