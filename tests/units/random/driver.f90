module test_cases
	use iso_c_binding
	use fruit
	use constants
	use randomness
	implicit none
	
contains
	
	subroutine randomint_ranges
		integer(bigint) :: result
		logical(sint)	:: test
		
		! typical range
		result = random_int(2, 21)
		test = (result .ge. 2) .and. (result .le. 21)
		call assert_true(test, 'random_int(2, 21) failed')
		
		! negative range
		result = random_int(-43, -8)
		test = (result .ge. -43) .and. (result .le. -8)
		call assert_true(test, 'random_int(-43, -8) failed')
		
		! large range
		result = random_int(-509, 12044)
		test = (result .ge. -509) .and. (result .le. 12044)
		call assert_true(test, 'random_int(-509, 12044) failed')
	end subroutine randomint_ranges
	
	subroutine randomint_special
		integer(bigint) :: result
		logical(sint)	:: test
		
		! adjacent values
		result = random_int(0, 1)
		test = result .eq. 0
		call assert_true(test, 'random_int(0, 1) failed')
		
		! equal values
		result = random_int(4, 4)
		test = result .eq. 4
		call assert_true(test, 'random_int(4, 4) failed')
		
		! reverse minv and maxv
		result = random_int(80, 21)
		test = (result .ge. 21) .and. (result .le. 80)
		call assert_true(test, 'random_int(80, 21) failed')
	end subroutine randomint_special

	subroutine nrandomint_test
		integer(bigint), dimension(32) :: results
		logical(sint)				   :: test
		integer						   :: ix
		
		results = 0
		call n_random_ints(32, 1, 100, results)
		test = .true.
		main: do ix=1,32
			test = (results(ix) .ge. 1) .and. (results(ix) .le. 100)
			if (.not. test) exit main
		end do main
		call assert_true(test, 'nrandomint_test failed')
	end subroutine nrandomint_test

end module test_cases


program random_test
  use fruit
  use test_cases
  
  implicit none

  logical(sint) :: ok
  
  call init_fruit

  ! call all tests
  call randomint_ranges
  call randomint_special
  call nrandomint_test

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) call exit(1)
  
end program random_test
