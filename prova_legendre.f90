include 'legendre_routines.f90'

!Function to compute the gravitational potential, given
!the gravitational constant G, Earth's mean radius a_e,
!a set of associated Legendre polynomials lpoly and relatives
!Stokes coefficient stokes_coef.
!lmax is maximum degree for which associated Legendre polynomials
!contained in lpoly have been computed, while length is the length
!of GRACE GSM file considered, stokes_coef.
function compute_v(G, a_e, r, phi, lambda, lpoly, lmax, stokes_coef, length) result(v_grav)
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, parameter :: qp = REAL128
  real(kind=qp), intent(in):: G, a_e
  integer :: lmax, length, l, m, i, j, nrow
  real(kind=qp), dimension(ceiling((real(lmax**2 + 4*lmax + 1))/(2)),&
                 &3), intent(in)&
                 & :: lpoly 
  real(kind=qp), dimension(length, 4), intent(in) :: stokes_coef
  real(kind=qp) :: v_grav, plm, clm, slm
  real(kind=qp), intent(in) :: r, phi, lambda
  real(kind=qp) :: scala 
  nrow = ceiling((real(lmax**2 + 4*lmax + 1))/2)

  v_grav = 0d0
  scala = r

  do l=2, lmax
    do m=0, l
      do i=1, nrow 
        if(lpoly(i,1)==l .and. lpoly(i,2)==m) then
          plm = lpoly(i, 3)
        end if
      end do
      do i=1, length
        if(stokes_coef(i,1) == l .and. stokes_coef(i,2)==m) then
          clm = stokes_coef(i, 3)
          slm = stokes_coef(i, 4) 
        end if
      end do
      v_grav = v_grav + scala*((a_e/r)**l)*plm*(clm*cos(m*lambda) + slm*sin(m*lambda)) 
    end do
  end do

  !v_grav = v_grav + temp

end function




!Program to process and prepare for plotting GRACE level 2
!data.
program legendre
  use, intrinsic :: iso_fortran_env !This is for using real128 
  implicit none
  integer, parameter :: qp = REAL128 !quadruple precision
  integer :: l=4, m=0, i, j
  real, parameter :: d2r = asin(1e0)/9e1 !declaration of a constant
  real(kind=qp), allocatable, dimension(:,:) :: ans
  real(kind=qp), allocatable, dimension(:,:) :: normalized
  integer :: nrow
  !Declaring gravitational constant and Earth's mean radius
  real(kind=qp), parameter :: G=6.6742d-11*5.9722d24, a_e = 6.37101d6  
  !Declaring Earth's Potential
  real(kind=qp) :: v_grav
  !Declaring return type of compute_v
  real(kind=qp) :: compute_v
  real(kind=qp), allocatable, dimension(:,:) :: grace_data
  integer :: in_unit = 1, nr=0, stat !Input device
  real(kind=qp) :: phi, lambda
  real(kind=qp) :: x, r
  !Computing # of row of result
  nrow = ceiling((real(l**2 + 4*l + 1))/(2))
  phi=45d0*d2r
  lambda=45d0*d2r
  x = sin(phi)
  !x=0.5

  !Reading GRACE file!
  open(unit=in_unit, file='grace_clean', action='read')
  !Here I'm counting # of record in data file
  do
    read(in_unit, *, iostat=stat)
    if(stat < 0) then
      print *, 'Reached eof'
      exit
    else if( stat > 0) then
      print *, 'Error in file format' 
      exit
    end if
    nr = nr + 1
  end do 
  !From the reference: the REWIND statement positions the file associated
  !with the specified unit to its initial point. If I don't rewind the file
  !(i.e. put the pointer on the initial position), I'll get a runtime error
  !the second time I read the file
  rewind(in_unit)
  allocate(grace_data(nr,4))
  !Here I'm reading and saving data in a 2-dim array
  do i=1, nr
    read(in_unit, *) grace_data(i, 1), grace_data(i, 2), grace_data(i, 3), grace_data(i, 4)
  end do 
  
  close(in_unit)

  !print *, 'sin(pi/2)', x
  allocate(ans(l+1, l+1))
  allocate(normalized(nrow, 3))
  call alegendre_increment_l(x, l, m, ans)
  call alegendre_increment_m(x, l, ans)
  call normalize_alegendre(ans, l, normalized)

  !G=G*5.9722d24
  r = a_e 
  v_grav = compute_v(G, a_e, r, phi, lambda, normalized, l, grace_data, nr)


  !do i=1, l
  !    print *, (ans(i,j), j=1,l+1)
  !end do

  do i=1, nrow
    print *, normalized(i, 1), ' ', normalized(i,2)&
            &, ' ', normalized(i,3)  
  end do

  print *, 'Potential: ', v_grav

  !do i=1, nr
  !  print *, grace_data(i, 1), grace_data(i, 2), grace_data(i, 3), grace_data(i, 4)
  !end do 

end program
