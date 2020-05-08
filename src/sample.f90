module sample
	use constants, only : dbl, bigint, damping, maxnchange, n_guesses, print_frequency
	use ioutil, only : energy
	use random
	use progress, only : progress_bar_time
	use hashtbl
	use nonradiative
	implicit none
contains	
	subroutine differences(n, values, res)
		integer, intent(in)						:: n
		real(dbl), dimension(n), intent(in) 	:: values
		real(dbl), dimension(n, n), intent(out)	:: res
		
		integer :: ix, jx
		res = 0d0
		do ix=1,n-1
			do jx=ix+1,n
				res(ix, jx) = values(ix) - values(jx)
				res(jx, ix) = -res(ix, jx)
			end do
		end do
	end subroutine differences
	
	subroutine add_occ(n, values, emin, emax, nsamples, noccs, occlist, enlist, tbl)
		integer, intent(in)								:: n, nsamples
		real(dbl), dimension(n), intent(in)				:: values
		integer(bigint), intent(inout)					:: noccs
		real(dbl), intent(in)							:: emin, emax
		integer, dimension(n, nsamples), intent(inout)	:: occlist
		real(dbl), dimension(nsamples), intent(inout)	:: enlist
		type(hash_tbl_sll), intent(inout)				:: tbl
			
		integer							:: outix
		real(dbl)						:: en
		character(len=:), allocatable	:: hash_str
		call get_hash(n, occlist(1:n, noccs+1), hash_str)
		call tbl%get(hash_str, outix)
		if (outix == 0) then
			en = energy(n, values, values, occlist(:, noccs+1))
			if ((en > emin) .and. (en < emax)) then
				noccs = noccs + 1
				call tbl%put(hash_str, int(noccs))
				enlist(noccs) = en
			end if
		end if
	end subroutine add_occ
	
	subroutine uniform_indices(sys, occ, plusminus, indices)
		type(sysdata), intent(in)						:: sys
		integer, dimension(sys%nlevels), intent(in)		:: occ
		integer, dimension(2), intent(out)				:: plusminus
		integer, dimension(2*maxnchange), intent(out)	:: indices
		call n_random_ints(2, 1, maxnchange, plusminus)
		call n_random_ints(sum(plusminus), 1, sys%nlevels, indices)
	end subroutine uniform_indices
	
	subroutine weighted_indices(sys, occ, plusminus, indices)
		type(sysdata), intent(in)						:: sys
		integer, dimension(sys%nlevels), intent(in)		:: occ
		integer, dimension(2), intent(out)				:: plusminus
		integer, dimension(2*maxnchange), intent(out)	:: indices

		real(dbl), dimension(sys%nlevels)	:: weights
		integer								:: ix
		call n_random_ints(2, 1, maxnchange, plusminus)
		do ix=1,sys%nlevels
			weights(ix) = sys%fcfactors(ix, occ(ix)+1)
		end do
		call n_random_ints_weighted(sum(plusminus), sys%nlevels, plusminus, weights, indices)
	end subroutine weighted_indices
	
	subroutine do_sample(sys, guesses, occs, guess_ens, emin, emax, ixfunc, tbl, occlist, enlist, noccs)
		type(sysdata), intent(in)									:: sys
		integer, dimension(sys%nlevels), intent(in)					:: occs
		integer, dimension(n_guesses, sys%nlevels), intent(in)		:: guesses
		real(dbl), intent(in)										:: emin, emax
		real(dbl), dimension(n_guesses), intent(in)					:: guess_ens
		integer, dimension(sys%nlevels, sys%nsamples), intent(out)	:: occlist
		real(dbl), dimension(sys%nsamples), intent(out)				:: enlist
		integer(bigint), intent(out)								:: noccs
		type(hash_tbl_sll), intent(out)								:: tbl
			
		interface
			subroutine ixfunc(sys, occ, plusminus, indices)
				import :: sysdata, maxnchange
				type(sysdata), intent(in)						:: sys
				integer, dimension(sys%nlevels), intent(in)		:: occ
				integer, dimension(2), intent(out)				:: plusminus
				integer, dimension(2*maxnchange), intent(out)	:: indices
			end subroutine
		end interface
		
		integer								:: jx, kx, lx, enix, occix, sn, npm
		integer, dimension(2)				:: pmix
		integer, dimension(2*maxnchange)	:: ixes
		real(dbl)							:: delta, delta_minus, delta_plus
		character(len=:), allocatable		:: hash_str
		
		if (.not. tbl%is_init) call tbl%init(sys%nsamples)
		noccs = 0
		do jx=1,n_guesses
			occlist(:, noccs+1) = guesses(jx, :)
			call add_occ(sys%nlevels, sys%energies, emin, emax, sys%nsamples, noccs, occlist, enlist, tbl)
		end do
		
		do jx=1,sys%nsamples
			if (mod(jx, print_frequency) .eq. 1) call progress_bar_time(int(jx, kind=bigint), int(sys%nsamples, kind=bigint))
			
			occix = random_int(1, int(noccs))
			call ixfunc(sys, occlist(:, occix), pmix, ixes)
			delta_minus = (emin - enlist(occix)) * damping
			delta_plus  = (emax - enlist(occix)) * damping
			delta = 0d0
			occlist(:, noccs+1) = occlist(:, occix)
			npm = sum(pmix)
			do kx=1,npm
				enix = ixes(kx)
				sn 	 = (2*mod(kx, 2) - 1) * sign(1, pmix(1)-1) * sign(1, pmix(2)-1)
				lx 	 = mod(3-2*sn, 3)
				pmix(lx) = pmix(lx) - 1
				lx = occlist(enix, noccs+1) + sn
				if ((lx .ge. 0) .and. (lx .le. occs(enix))) then
					delta = delta + sn*sys%energies(enix)
					occlist(enix, noccs+1) = lx
					if ((delta .gt. delta_minus) .and. (delta .lt. delta_plus)) then
						call add_occ(sys%nlevels, sys%energies, emin, emax, sys%nsamples, noccs, occlist, enlist, tbl)
					end if
				end if
			end do
		end do
		call progress_bar_time(int(sys%nsamples, kind=bigint), int(sys%nsamples, kind=bigint))
		write(*, *) '\n'
	end subroutine do_sample
	
	subroutine get_hash(n, occ, hash_str)
		integer, intent(in) 						:: n
		integer, dimension(n), intent(in)			:: occ
		character(len=:), allocatable, intent(out)	:: hash_str
		
		character(len=20)	:: fmt
		integer				:: ix
		write(fmt,'("(", i0, "i3)")') n
		
		if (allocated(hash_str)) deallocate(hash_str)
		allocate(character(len=3*n) :: hash_str)
		write(hash_str, fmt) (occ(ix), ix=1,n)
	end subroutine get_hash
	
	subroutine get_unique(n, nocc, occlist, tbl, nunique)
		integer, intent(in)						:: n, nocc
		integer, dimension(n, nocc), intent(in)	:: occlist
		integer, intent(out)					:: nunique
		type(hash_tbl_sll), intent(out)			:: tbl
		
		integer 						:: ix, jx, outix
		character(len=:), allocatable	:: hash_str
		
		nunique = 0
		if (.not. tbl%is_init) call tbl%init(nocc)
		do ix=1,nocc
			call get_hash(n, occlist(1:n, ix), hash_str)
			call tbl%get(hash_str, outix)
			if (outix == 0) then
				call tbl%put(hash_str, ix)
				nunique = nunique + 1
			end if
		end do
	end subroutine get_unique
end module sample