module randomness
	use constants, only : dbl, bigint
#ifdef __INTEL_COMPILER
        use ifport
#endif
	implicit none
contains
	function random_int(minv, maxv)
		integer(bigint), intent(in) :: minv
		integer(bigint), intent(in) :: maxv
		integer(bigint)		:: random_int
		real(dbl)			:: rn, scale 
		
		rn     = rand()
		scale  = maxv-minv
		random_int = minv + floor(rn*scale)
	end function random_int
	
	subroutine n_random_ints(n, minv, maxv, res)
		integer(bigint), intent(in)							:: n, minv, maxv
		integer(bigint), dimension(n), intent(out) 	:: res
		
		real(dbl)	:: rn, scale
		integer(bigint)		:: ix
		scale = maxv - minv
		do ix=1, n
			rn = rand()
			res(ix) = minv + floor(rn*scale)
		end do
	end subroutine n_random_ints
	
	subroutine n_random_ints_weighted(n, nvals, plusminus, weights, res)
		integer(bigint), intent(in)					:: n, nvals
		integer(bigint), dimension(2), intent(in)	:: plusminus
		real(dbl), dimension(nvals), intent(inout)	:: weights
		integer(bigint), dimension(n), intent(out) 	:: res
		
		integer(bigint), dimension(:), allocatable	:: expanded_ints
		integer(bigint)								:: mults(nvals), posres(n), negres(n)
		integer(bigint)								:: ix, nt, ctr, pmmin, unity
		real(dbl)									:: mv
		
		pmmin = min(plusminus(1), plusminus(2))
		unity = 1
		
		! Prob distribution for adding occupations
		mv = minval(weights)
		if (mv .eq. 0) mv = 0.01d0
		weights = log10(weights / mv)
		do ix=1,nvals
			weights(ix) = max(weights(ix), 0d0)
		end do
		mults = ceiling(weights)+1
		nt = sum(mults)
		allocate(expanded_ints(nt))
		ctr = 1
		do ix=1,nvals
			expanded_ints(ctr:ctr+mults(ix)-1) = ix
			ctr = ctr + mults(ix)
		end do 
		call n_random_ints(plusminus(1), unity, nt+1, posres)
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
		call n_random_ints(plusminus(2), unity, nt+1, negres)
		do ix=1,pmmin
			res(2*ix) = expanded_ints(negres(ix))
		end do
		if (pmmin .eq. plusminus(1)) res(2*pmmin+1:) = [(expanded_ints(negres(ix)), ix=pmmin+1,plusminus(2))]

		deallocate(expanded_ints)
	end subroutine n_random_ints_weighted
end module randomness

