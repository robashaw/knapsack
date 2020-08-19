module sample
	use constants
	use ioutil, only : energy
	use randomness
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
	
	subroutine add_occ(sys, emin, emax, noccs, occ, enlist)
		type(sysdata), intent(inout)							:: sys
		integer(bigint), intent(inout)							:: noccs
		real(dbl), intent(in)									:: emin, emax
		integer(smallint), dimension(sys%nlevels), intent(in)	:: occ
		real(dbl), dimension(sys%mm%chunk_size), intent(inout)	:: enlist
			
		integer		:: i, outnocc
		real(dbl)	:: en
		en = energy(sys%nlevels, sys%energies, sys%energies, occ(:))
		if ((en > emin) .and. (en < emax)) then
			noccs = noccs + 1
			
			if (noccs .gt. sys%mm%chunk_size) then
				call sys%mm%block_swap(sys%energies, sys%hrfactors, noccs-1, outnocc)
				do i=1,outnocc
					enlist(i) = energy(sys%nlevels, sys%energies, sys%energies, sys%mm%current_block(:, i))
				end do
				write(*, '(1x,a,1x,i10,1x,a,1x,i10)') 'Kept top', outnocc, 'samples out of', noccs
				noccs = outnocc + 1
			end if
			
			sys%mm%current_block(:, noccs) = occ(:)
			enlist(noccs) = en
		end if
	end subroutine add_occ
	
	subroutine uniform_indices(sys, occ, plusminus, indices)
		type(sysdata), intent(in)						:: sys
		integer(smallint), dimension(sys%nlevels), intent(in)		:: occ
		integer, dimension(2), intent(out)				:: plusminus
		integer, dimension(2*maxnchange), intent(out)	:: indices
		call n_random_ints(2, 1, maxnchange, plusminus)
		call n_random_ints(sum(plusminus), 1, sys%nlevels, indices)
	end subroutine uniform_indices
	
	subroutine weighted_indices(sys, occ, plusminus, indices)
		type(sysdata), intent(in)						:: sys
		integer(smallint), dimension(sys%nlevels), intent(in)		:: occ
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
	
	subroutine do_sample(sys, guesses, occs, guess_ens, emin, emax, ixfunc, enlist, noccs)
		type(sysdata), intent(inout)										:: sys
		integer(smallint), dimension(sys%nlevels), intent(in)				:: occs
		integer(smallint), dimension(n_guesses, sys%nlevels), intent(in)	:: guesses
		real(dbl), intent(in)												:: emin, emax
		real(dbl), dimension(n_guesses), intent(in)							:: guess_ens
		real(dbl), dimension(sys%mm%chunk_size), intent(out)				:: enlist
		integer(bigint), intent(out)										:: noccs
			
		interface
			subroutine ixfunc(sys, occ, plusminus, indices)
				import :: sysdata, maxnchange, smallint
				type(sysdata), intent(in)						:: sys
				integer(smallint), dimension(sys%nlevels), intent(in)		:: occ
				integer, dimension(2), intent(out)				:: plusminus
				integer, dimension(2*maxnchange), intent(out)	:: indices
			end subroutine
		end interface
		
		integer										:: jx, kx, lx, enix, occix, sn, npm, tmp, mergerecord
		integer, dimension(2)						:: pmix
		integer, dimension(2*maxnchange)			:: ixes
		integer(smallint), dimension(sys%nlevels)	:: cocc
		real(dbl)									:: delta, delta_minus, delta_plus, pct
		character(len=100)							:: mergefile, tmpfile
		
		noccs = 0
		do jx=1,n_guesses
			call add_occ(sys, emin, emax, noccs, guesses(jx, :), enlist)
		end do
		
		do jx=1,sys%nsamples
			if (mod(jx, print_frequency/10) .eq. 0) then
				pct = real(jx)/real(sys%nsamples) * 100.0
				write(*, '(1x,a,1x,i12,1x,a,1x,f5.1,1x,a)') 'Done', jx, 'samples (', pct, '%)'
			end if
			
			occix = random_int(1, int(noccs))
			call ixfunc(sys, sys%mm%current_block(:, occix), pmix, ixes)
			delta_minus = (emin - enlist(occix)) * damping
			delta_plus  = (emax - enlist(occix)) * damping
			delta = 0d0
			cocc(:) = sys%mm%current_block(:, occix)
			npm = sum(pmix)
			do kx=1,npm
				enix = ixes(kx)
				sn 	 = (2*mod(kx, 2) - 1) * sign(1, pmix(1)-1) * sign(1, pmix(2)-1)
				lx 	 = mod(3-2*sn, 3)
				pmix(lx) = pmix(lx) - 1
				lx = cocc(enix) + sn
				if ((lx .ge. 0) .and. (lx .le. occs(enix))) then
					delta = delta + sn*sys%energies(enix)
					cocc(enix) = lx
					if ((delta .gt. delta_minus) .and. (delta .lt. delta_plus)) then
						call add_occ(sys, emin, emax, noccs, cocc, enlist)
					end if
				end if
			end do
		end do
		write(*, *) '\n'
		
		! Now we need to write the last file and merge all the files
		call sys%mm%block_swap(sys%energies, sys%hrfactors, noccs-1, tmp)
		call sys%mm%merge_all(sys%energies, sys%hrfactors, mergefile, mergerecord)
		
		! Calculate rate from merged occ file
		write(*, *) 'Calculating rate'
		call sys%calculate_from_file(mergerecord, mergefile, noccs)
		
		if (sys%do_write) then
			tmpfile = sys%mm%occprefix // '.final'
			write(*, '(1x,a,1x,a)') 'Saving final, unique occs to file as', tmpfile
			call sys%mm%move_file(mergefile, mergerecord, tmpfile)
		end if
		
		write(*, *) 'Cleaning up temporary files'
		call sys%mm%clean_up(sys%mm%occprefix, minix=1, maxix=sys%mm%current_record-1)
		if (mergefile .eq. 'merged') call sys%mm%clean_up('merged', minix=1, maxix=mergerecord)
	end subroutine do_sample
	
	subroutine get_hash(n, occ, hash_str)
		integer, intent(in) 						:: n
		integer(smallint), dimension(n), intent(in)	:: occ
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
		integer(smallint), dimension(n, nocc), intent(in)	:: occlist
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
