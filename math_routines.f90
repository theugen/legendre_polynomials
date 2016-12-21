! This is a program which computes factorial of a number

! here I try to declare a function

function fattoriale(n) result(res)
  implicit none
  integer, intent(in) :: n
  integer :: i
  ! intent(in) attribute of argument means that the argument  cannot be replaced inside function
  integer(kind=8) :: res, temp = 1
  
  if(n /= 0) then
    do i = n, 1, -1
      if(i ==  n) then
        res = n
      else 
        res = res * i
      end if
    end do
  else
    res = 1
  end if

end function fattoriale
  
