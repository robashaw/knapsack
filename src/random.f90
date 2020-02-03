module random
	use constants, only : dbl
	implicit none
contains
	function random_int(min, max)
		integer, intent(in) :: min
		integer, intent(in) :: max
		integer				:: random_int
		real(dbl)			:: rn, scale 
		
		rn     = rand()
		scale  = max-min
		random_int = min + floor(rn*scale)
	end function random_int
	
	subroutine n_random_ints(n, min, max, res)
		integer, intent(in)					:: n, min, max
		integer, dimension(n), intent(out) 	:: res
		
		real(dbl)	:: rn, scale
		integer		:: ix
		scale = max - min
		do ix=1, n
			rn = rand()
			res(ix) = min + floor(rn*scale)
		end do
	end subroutine n_random_ints
	
	subroutine n_random_ints_weighted(n, nvals, plusminus, weights, res)
		integer, intent(in)							:: n, nvals
		integer, dimension(2), intent(in)			:: plusminus
		real(dbl), dimension(nvals), intent(inout)	:: weights
		integer, dimension(n), intent(out) 			:: res
		
		integer, dimension(:), allocatable	:: expanded_ints
		integer								:: mults(nvals), posres(n), negres(n)
		integer								:: ix, nt, ctr, pmmin
		
		pmmin = min(plusminus(1), plusminus(2))
		
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
		call n_random_ints(plusminus(1), 1, nt, posres)
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
		call n_random_ints(plusminus(2), 1, nt, negres)
		do ix=1,pmmin
			res(2*ix) = expanded_ints(negres(ix))
		end do
		if (pmmin .eq. plusminus(1)) res(2*pmmin+1:) = [(expanded_ints(negres(ix)), ix=pmmin+1,plusminus(2))]

		deallocate(expanded_ints)
	end subroutine n_random_ints_weighted
end module random

