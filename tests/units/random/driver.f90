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
			test = test .and. ((results(ix) .ge. 1) .and. (results(ix) .le. 100))
		end do main
		call assert_true(test, 'nrandomint_test failed')
	end subroutine nrandomint_test
	
	subroutine nrandomint_weighted_uniform
		! n = 20, nvals = 5
		integer(bigint), dimension(2)	:: plusminus
		real(dbl), dimension(5)			:: weights
		integer(bigint), dimension(20)	:: result
		integer, dimension(5)			:: counter
		logical(sint)					:: test
		integer							:: ix, jx, bigdiff
		
		plusminus(1) = 10
		plusminus(2) = 10
		weights = 1d0 ! uniformly weighted
		
		call n_random_ints_weighted(20, 5, plusminus, weights, result)
		test = .true.
		counter = 0
		do ix=1,20
			test = test .and. ((result(ix) .ge. 1) .and. (result(ix) .le. 5))
			if (test) counter(result(ix)) = counter(result(ix)) + 1
		end do
		call assert_true(test, 'n_random_ints_weighted out of bounds error')
		
		if (test) then
			bigdiff = 0
			do ix=1,4
				do jx=ix+1,5
					bigdiff = max(abs(counter(ix)-counter(jx)), bigdiff)
				end do
			end do
			test = bigdiff .le. 5
			call assert_true(test, 'n_random_ints_weighted nonuniform weightings')
		end if
	end subroutine nrandomint_weighted_uniform
		
	subroutine nrandomint_weighted_binomial
		! n = 14, nvals = 6
		integer(bigint), dimension(2)	:: plusminus
		real(dbl), dimension(6)			:: weights
		integer(bigint), dimension(14)	:: result
		logical(sint)					:: test
		integer							:: ix
		
		plusminus(1) = 9
		plusminus(2) = 5
		weights(1) = 1d0 ! binomially weighted
		weights(2) = 5d0
		weights(3) = 10d0
		weights(4) = 10d0
		weights(5) = 5d0
		weights(6) = 1d0
		
		call n_random_ints_weighted(14, 6, plusminus, weights, result)
		test = .true.
		do ix=1,14
			test = test .and. ((result(ix) .ge. 1) .and. (result(ix) .le. 6))
		end do
		call assert_true(test, 'n_random_ints_weighted binomial weights failed')
	end subroutine nrandomint_weighted_binomial
	
	subroutine nrandomint_weighted_zero
		! n = 5, nvals = 3
		integer(bigint), dimension(2)	:: plusminus
		real(dbl), dimension(3)			:: weights
		integer(bigint), dimension(5)	:: result
		logical(sint)					:: test
		integer							:: ix
		
		plusminus(1) = 2
		plusminus(2) = 3
		weights = 0d0 ! all weights zero
		
		call n_random_ints_weighted(5, 3, plusminus, weights, result)
		test = .true.
		do ix=1,5
			test = test .and. ((result(ix) .ge. 1) .and. (result(ix) .le. 3))
		end do
		call assert_true(test, 'n_random_ints_weighted zero weights failed')
	end subroutine nrandomint_weighted_zero

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
  call nrandomint_weighted_uniform
  call nrandomint_weighted_binomial
  call nrandomint_weighted_zero

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) call exit(1)
  
end program random_test
