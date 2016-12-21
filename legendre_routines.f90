!This file contains subroutines for computing the values of associated
!Legendre polynomials through recursive relationships. 
!Author: Eugenio Thieme, eugen@thieme.it

include 'math_routines.f90'

!This subroutine is for computing the value of the polynomial of degree
!l+1 from the polynomial of degree l. Order m remains the same 
subroutine alegendre_increment_l(x, lmax, m, ans)
  use, intrinsic :: iso_fortran_env  !This is for using real128
  implicit none
  integer, parameter :: qp = REAL128  !quadruple precision
  integer :: lmax, m, l
  real(kind=qp) :: x
  real(kind=qp), dimension(lmax+1) :: ans
 
  !For avoiding confusion between array index and polynomial degree,
  !I voluntarily chose to start computing from degree 1 (P_1_m).
  !Remember that P_0_0 = 1!
  do l=1, lmax
    if(l == 1) then
      ans(l) = x
    else if(l == 2) then
      ans(l) = ((2*(l-1)+1)*x*ans(l-1)-(l-1+m))/(l-m)
    else  
      ans(l) = (((2*(l-1)+1)*x*ans(l-1))-((l-1+m)*ans(l-2)))/(l-m)
    end if
  end do
end subroutine

!This subroutine is for computing the value of the polynomial of degree
!l and order m+1, starting from P_0_0. The result is a matrix of dimension
!l*(l+1), with the element (1,1) equivalent to P_0_1. This subroutine normally
!should be called after alegendre_increment_l(), because it needs to know
!all Legendre associated polynomials with m=0, from l=1 to l=lmax
subroutine alegendre_increment_m(x, lmax, ans)
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, parameter :: qp = REAL128
  integer :: m, lmax, l
  real(kind=qp) :: x
  real(kind=qp), dimension(lmax+1, lmax+1) :: ans

  do l=1, lmax
    do m=0, l-1
      if(l == 1) then 
        ans(l, m+2) = (x**2-1)/(sqrt(1-x**2))
      else
          ans(l, m+2) = ((l-m)*x*ans(l,m+1)-(l+m)*ans(l-1, m+1))/(sqrt(1-x**2))
      end if
    end do
  end do
end subroutine


!Subroutine for computing fully normalized associated Legendre polynomials. 
!Result will be organized like GRACE level-2 data:
!first column represents degree l, second order m, and last the value of fully
!normalized Legendre. Order of rows: first all m=0 (from l=1 to l=lmax), 
!then all m=1,2,... and so on.
!For doubts, please refer to GRACE L2-User-Handbook  
subroutine normalize_alegendre(leg_com, lmax, normalized)
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, parameter :: qp = REAL128
  integer :: lmax, m, l, i=1
  integer(kind=8) :: fattoriale !These are the result of factorial
  real(kind=qp), dimension(lmax+1, lmax+1) :: leg_com 
  integer :: nrow
!Can't use allocatable objects inside subroutines unless I provide explicit
!interface... but I don't know yet what an explicit interface is :(
!  real(kind=qp), allocatable, dimension(:,:) :: normalized 
!  real(kind=qp), dimension(ceiling((real(lmax**2 + 4*lmax + 1))/(2))&
!                 &,ceiling((real(lmax**2 + 4*lmax + 1))/(2))) :: normalized 

  real(kind=qp), dimension(ceiling((real(lmax**2 + 4*lmax + 1))/(2))&
                 &,3) :: normalized 
  !Computing # of row of result
  nrow = ceiling((real(lmax**2 + 4*lmax + 1))/(2)) 
  !allocate(normalized(nrow,3))

  do m=1,(lmax+1),1
    do l=1,lmax
      if((m-1)==0) then
        normalized(i,1) = l
        normalized(i,2) = m-1
        normalized(i,3) = sqrt(real((2*l+1)))*leg_com(l,m)
        i = i+1
      else if(m <= l+1) then
        normalized(i,1) = l
        normalized(i,2) = m-1
        normalized(i,3) = sqrt(real(2*(2*l+1)*fattoriale(l-m+1))/fattoriale(l+m-1))*leg_com(l,m)
        i = i+1
      end if
    end do
  end do

end subroutine



