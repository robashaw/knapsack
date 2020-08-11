module randomness
  use constants, only : dbl, bigint
#ifdef __INTEL_COMPILER
  use ifport
#endif
  implicit none
contains
	function random_int(min, max)
		integer(bigint), intent(in) :: min
		integer(bigint), intent(in) :: max
		integer(bigint)		:: random_int
		real(dbl)			:: rn, scale 
		
		rn     = rand()
		scale  = max-min
		random_int = min + floor(rn*scale)
	end function random_int
	
	subroutine n_random_ints(n, min, max, res)
		integer(bigint), intent(in)							:: n, min, max
		integer(bigint), dimension(n), intent(out) 	:: res
		
		real(dbl)	:: rn, scale
		integer(bigint)		:: ix
		scale = max - min
		do ix=1, n
			rn = rand()
			res(ix) = min + floor(rn*scale)
		end do
	end subroutine n_random_ints
	
	subroutine n_random_ints_weighted(n, nvals, plusminus, weights, res)
		integer(bigint), intent(in)							:: n, nvals
		integer(bigint), dimension(2), intent(in)			:: plusminus
		real(dbl), dimension(nvals), intent(inout)	:: weights
		integer(bigint), dimension(n), intent(out) 	:: res
		
		integer(bigint), dimension(:), allocatable	:: expanded_ints
		integer(bigint)						:: mults(nvals), posres(n), negres(n)
		integer(bigint)								:: ix, nt, ctr, pmmin, unity
		
		pmmin = min(plusminus(1), plusminus(2))
		unity = 1
		
		! Prob distribution for adding occupations
		weights = log10(weights / minval(weights))
		mults = ceiling(weights)+1
		nt = sum(mults)
		allocate(expanded_ints(nt))
		ctr = 1
		do ix=1,nvals
			expanded_ints(ctr:ctr+mults(ix)-1) = ix
			ctr = ctr + mults(ix)
		end do 
		call n_random_ints(plusminus(1), unity, nt, posres)
		do ix=1,pmmin
			res(2*ix-1) = expanded_ints(posres(ix))
		end do
		if (pmmin .eq. plusminus(2)) res(2*pmmin+1:) = [(expanded_ints(posres(ix)), ix=pmmin+1,plusminus(1))]
		
		! Prob distribution for lowering occupations
		mults = maxval(mults) - mults + 1
		nt = sum(mults)
		deallocate(expanded_ints)
		allocate(expanded_ints(nt))
		ctr = 1
		do ix=1,nvals
			expanded_ints(ctr:ctr+mults(ix)-1) = ix
			ctr = ctr + mults(ix)
		end do 
		call n_random_ints(plusminus(2), unity, nt, negres)
		do ix=1,pmmin
			res(2*ix) = expanded_ints(negres(ix))
		end do
		if (pmmin .eq. plusminus(1)) res(2*pmmin+1:) = [(expanded_ints(negres(ix)), ix=pmmin+1,plusminus(2))]

		deallocate(expanded_ints)
	end subroutine n_random_ints_weighted
end module randomness

