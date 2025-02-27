subroutine compute_array(a, b, n)
  integer, intent(in) :: n
  real, dimension(n), intent(inout) :: a
  real, dimension(n), intent(in) :: b
  integer :: i

  !$acc routine vector
  !$acc loop vector
  do i = 1, n, 1
    a(i) = a(i) + b(i)
  enddo

end subroutine compute_array
