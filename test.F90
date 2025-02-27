subroutine compute_array(a, b, n)
  implicit none
  integer, intent(in) :: n
  real, intent(inout) :: a(n)
  real, intent(in) :: b(n)
  integer :: i

  do i = 1, n
    a(i) = a(i) + b(i)
  end do

end subroutine compute_array

