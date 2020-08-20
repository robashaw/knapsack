module test_cases
	use iso_c_binding
	use fruit
	use constants
	implicit none
	
contains
	
	subroutine factorial_positive
		call assert_equals(1d0, factorial(1), 'factorial(1) failed')
		call assert_equals(120d0, factorial(5), 'factorial(5) failed')
		call assert_equals(3.55687428096D+14, factorial(17), 'factorial(17) failed')
	end subroutine factorial_positive
	
	subroutine factorial_zero
		call assert_equals(1d0, factorial(0), 'factorial_zero failed')
	end subroutine factorial_zero
	
	subroutine factorial_negative
		call assert_equals(1d0, factorial(-5), 'factorial_negative failed')
	end subroutine factorial_negative
	
	subroutine time_check
		character(len=100) 	:: result
		logical(sint)  		:: check
		result = trim(adjustl(format_time(14985.781d0)))
		! should give 4h: 9m:45s 
		check =  '4h: 9m:45s' .eq. result
		write(*, *) result
		call assert_true(check, 'time_check failed')
	end subroutine time_check

end module test_cases


program constants_test
  use fruit
  use test_cases
  
  implicit none

  logical(sint) :: ok
  
  call init_fruit

  ! call all tests
  call factorial_positive
  call factorial_zero
  call factorial_negative
  call time_check

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) call exit(1)
  
end program constants_test
