program main
	use constants
	use nonradiative
	use radiative
	use sample
	use knapsack 
	use ioutil
#ifdef __INTEL_COMPILER
	use ifport
#endif
        
	character(len=100)						:: arg, mergefile, tmpfile
	integer									:: ios
	type(sysdata)							:: sys
	integer(smallint), dimension(:, :), allocatable 	:: occlist
	real(dbl), dimension(:), allocatable	:: enlist
	integer(smallint), dimension(:, :), allocatable	:: guesses
	integer									:: ix, noccs, ntot, tmp, mergerecord, ncutm
	integer(bigint)							:: worst_case, bignoccs
	real(dbl)								:: emin, emax, start, finish, target_en, klag
	real(dbl), dimension(2*n_guesses)		:: guess_ens
	type(hash_tbl_sll) 						:: tbl
	
	character(8)  :: datechar
	character(10) :: timechar
	call date_and_time(datechar, timechar)
	write (*, *) 'KNAPSACK, version 1.02, Aug 2020'
	write (*, *) 'by Robert A. Shaw'
	write (*, '(1x, a, 2x, 5a, 2x, 5a)') 'Called at:', datechar(1:4), '/', datechar(5:6), '/', datechar(7:8), &
		 											   timechar(1:2), ':', timechar(3:4), ':', timechar(5:6)
	
	call get_command_argument(1, arg, status=ios)
	if (ios .eq. 0) then
		call cpu_time(start)
		call sys%from_file(arg)
		
  	 	select case (sys%calctype)
	 	case ('radiative')
	 		sys%do_rad = .true.
	 		sys%do_nonrad = .false.
	 	case ('nonradiative')
 	 		sys%do_rad = .false.
 	 		sys%do_nonrad = .true.
	 	case ('all')
 	 		sys%do_rad = .true.
 	 		sys%do_nonrad = .true.
	 	case default
 	 		sys%do_rad = .false.
 	 		sys%do_nonrad = .false.
			write (*, *) 'Nothing to do here? No calculation type set.'
	 	end select
		
		call sys%init
		
		if (sys%do_rad) then
			write (*, *)
			write (*, *) 'RADIATIVE RATE CALCULATION'
			call calculate_rate(sys%k_r, sys%radfile, sys%radunits)
			sys%k_r = sys%k_r * RAD_CONVERT * (sys%tdm**2)
			write (*, '(1x,a,1x,e14.8,1x,a)') 'k_r =', sys%k_r, 's-1' 
		end if
		
		if (sys%do_nonrad) then
			write (*, *)
			write (*, *) 'NON-RADIATIVE RATE CALCULATION'
			write (*, '(1x,a,i2,a,a)') 'Read in ', sys%nlevels, ' energy levels from file ', arg
			write (*, *)
			!if ((sys%algorithm .ne. 'brute') .and. (sys%algorithm .ne. 'stochastic')) sys%algorithm = 'hybrid'
			noccs = 1
			
			if (sys%debug_level .gt. 0) then
				write(*, *)
				write(*, '(a10,1x,a9,1x,a10)') 'ENERGY', 'MAX. OCC.', 'HR FACTOR'
				do ix=1,sys%nlevels
					write(*, '(f10.6,1x,i9,1x,f10.6)') sys%energies(ix), sys%maxoccs(ix), sys%hrfactors(ix)
				end do
			end if
			
			emin = sys%e_target * TO_EV - sys%delta_e
			emax = emin + 2*sys%delta_e
			noccs = 0
			write (*, '(1x,a,f10.6,1x,a,1x,f10.6,1x,a)') 'Searching for occs with E(n) between', emin, 'and', emax, 'eV'
			write (*, '(1x,a,a,1x,a,i2)') 'Using algorithm: ', trim(sys%algorithm), 'with ncut = ', sys%ncut
			allocate(occlist(sys%nlevels, 1))
			occlist = 0
			call sys%calculate_kic(1, occlist, init=.true., stopix=0)
			
			call lambda(sys, emax, emin, klag)
			select case (sys%algorithm)
			case ('brute')
				call brute_force(sys, emax, emin, enlist, noccs)
			case ('stochastic')
				sys%mm%keep_top = .true.
				ncutm = n_cut_guess
				do while (ncutm * sys%nlevels < sys%maxnfix)
					ncutm = ncutm + 1
				end do
				sys%maxoccs = [(min(ncutm, sys%cutoffs(ix, sys%nthresh+1)), ix=1,sys%nlevels)]
				allocate(guesses(2*n_guesses, sys%nlevels))
				write(*, '(1x,a,i2)') "Generating guesses with ncut = ", n_cut_guess
				do ix=1,n_guesses
					target_en = emin + ((real(ix)-0.5d0)/real(n_guesses)) * (emax - emin)
					call knap_n(target_en, sys%nlevels, sys%energies, sys%maxoccs, 0.05d0, sys%lambda_n,&
					 			guesses(2*ix-1, :), guess_ens(2*ix-1), .true.)
					write (*, *) guesses(2*ix-1, :), guess_ens(2*ix-1)
					call knap_n(target_en, sys%nlevels, sys%energies, sys%maxoccs, 0.05d0, sys%lambda_n,& 
								guesses(2*ix, :), guess_ens(2*ix), .false.)
					write (*, *) guesses(2*ix, :), guess_ens(2*ix)
				end do
				sys%maxoccs = [(sys%cutoffs(ix, sys%nthresh+1), ix=1,sys%nlevels)]
				write(*, '(1x,a,1x,i20,1x,a)') "\nSampling with", sys%nsamples, "samples\n"
				allocate(enlist(sys%nsamples))
				select case (sys%weighting)
				case ('uniform')
					call do_sample(sys, guesses, sys%maxoccs, guess_ens, emin, emax, uniform_indices, enlist, bignoccs)
				case default
					call do_sample(sys, guesses, sys%maxoccs, guess_ens, emin, emax, weighted_indices, enlist, bignoccs)
				end select
				bignoccs = bignoccs + 1
				noccs = sys%mm%chunk_size
			case ('fixedn')
				worst_case = sys%maxnoccs
				call screened_fixed_n(sys, emax, emin, enlist, noccs, worst_case)
			case ('fromfile')
				tmpfile = trim(adjustl(sys%mm%occprefix)) // '.*'
				write(*, '(1x,a,1x,a)') 'Reading occs from', tmpfile
				sys%mm%nrecords = sys%nrecords
				write(*, '(1x,a,1x,i4,1x,a)') 'Looking for', sys%mm%nrecords, 'records'
				do ix=1,sys%mm%nrecords
					call sys%calculate_from_file(ix, sys%mm%occprefix, ntot)
					noccs = noccs + ntot
				end do
			case ('mergefiles')	
				sys%mm%current_record = sys%nrecords + 1
				write(*, '(1x,a,1x,a,1x,a,1x,i3)') 'Merging files with prefix', sys%mm%occprefix, 'up to record', sys%mm%nrecords
				call sys%mm%merge_all(sys%energies, sys%hrfactors, mergefile, mergerecord)
				tmpfile = trim(adjustl(sys%mm%occprefix)) // '.final'
				write(*, '(1x,a,1x,a)') 'Saving final, unique occs to file as', tmpfile
				call sys%mm%move_file(mergefile, mergerecord, tmpfile)
				if (mergefile .eq. 'merged') call sys%mm%clean_up('merged', minix=1, maxix=mergerecord)
			case ('lagrange')
				sys%k_ic = 0.25 * klag * sys%gamma
			case default
				worst_case = sys%maxnoccs 
				call screened_brute_force(sys, emax, emin, enlist, noccs, worst_case)
			end select
			if ((sys%algorithm .eq. 'fixedn') .or. (sys%algorithm .eq. 'brute') .or. (sys%algorithm .eq. 'hybrid')) then
				write(*, '(1x,a,1x,i10)') 'Finished block', sys%mm%current_record
				call sys%calculate_kic(sys%mm%chunk_size, sys%mm%current_block, init=.false., stopix=noccs-1)
				if (sys%do_write) then
					call sys%mm%block_swap(sys%energies, sys%hrfactors, noccs-1, tmp)
				else
					sys%mm%current_record = sys%mm%current_record + 1
					sys%mm%current_block = 0
				end if
				
				noccs = (sys%mm%current_record-2)*sys%mm%chunk_size+noccs-1
			end if
			
			if (sys%algorithm .eq. 'stochastic') then
				write (*, '(1x,a,1x,i10,1x,a)') 'Found approximately', bignoccs, 'occs'
			else
				write (*, '(1x,a,1x,i10,1x,a)') 'Found approximately', noccs, 'occs'
			end if
			
			write (*, *)
			write (*, *) 'Calculating non-radiative rate'
			write (*, '(1x,a,1x,e14.8)') 'With calculated gamma =', sys%gamma
			sys%k_ic = 4d0 * sys%k_ic / sys%gamma 
			sys%k_ic_ht = 4d0 * sys%k_ic_ht / sys%gamma
			!call sys%calculate_kic(noccs, occlist, init=.true.)
			write (*, '(1x,a,1x,e14.8,1x,a)') 'k_ic_fc =', sys%k_ic, 's-1'
			if (sys%hessfile .ne. 'none') write (*, '(1x,a,1x,e14.8,1x,a)') 'k_ic_ht =', sys%k_ic_ht, 's-1'
			if (sys%user_gamma .gt. 0) then
				write (*, '(1x,a,1x,e14.8)') 'With user-input gamma =', sys%user_gamma
				write (*, '(1x,a,1x,e14.8,1x,a)') 'k_ic_fc =', sys%k_ic * (sys%gamma / sys%user_gamma), 's-1'
				if (sys%hessfile .ne. 'none') write (*, '(1x,a,1x,e14.8,1x,a)') 'k_ic_ht =', sys%k_ic_ht * (sys%gamma / sys%user_gamma), 's-1'
			end if
			write (*, *)
		end if
		
		if (sys%do_rad .and. sys%do_nonrad) then
			write (*, '(1x,a,1x,f12.8,1x,a)') 'PLQY =', 100d0 * sys%k_r/(sys%k_r + sys%k_ic + sys%k_ic_ht), '%'
		end if
		
		call cpu_time(finish)
		write (*, *) 'Time elapsed:', format_time(finish-start)
	else
		write (*, *) 'Could not open file', arg
	end if
end program main
