program main
	use constants, only : dbl, bigint, TO_EV, RAD_CONVERT, n_guesses, n_cut_guess, format_time
	use nonradiative
	use radiative
	use sample
	use knapsack 
	
	character(len=100)						:: arg
	integer									:: ios
	type(sysdata)							:: sys
	integer, dimension(:, :), allocatable 	:: occlist
	real(dbl), dimension(:), allocatable	:: enlist
	integer, dimension(:), allocatable		:: minoccs, maxoccs
	integer, dimension(:, :), allocatable	:: guesses
	integer									:: ix
	integer(bigint)							:: noccs, worst_case
	real(dbl)								:: emin, emax, start, finish, target_en
	real(dbl), dimension(n_guesses)			:: guess_ens
	type(hash_tbl_sll) 						:: tbl
	
	character(8)  :: date
	character(10) :: time
	call date_and_time(date, time)
	write (*, *) 'KNAPSACK, version 1.0, Dec 2019'
	write (*, *) 'by Robert A. Shaw'
	write (*, '(1x, a, 2x, 5a, 2x, 5a)') 'Called at:', date(1:4), '/', date(5:6), '/', date(7:8), &
		 											   time(1:2), ':', time(3:4), ':', time(5:6)
	
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
		
		if (sys%do_rad) then
			write (*, *)
			write (*, *) 'RADIATIVE RATE CALCULATION'
			call calculate_rate(sys%k_r, sys%radfile)
			sys%k_r = sys%k_r * RAD_CONVERT * (sys%tdm**2)
			write (*, '(1x,a,1x,e14.8,1x,a)') 'k_r =', sys%k_r, 's-1' 
		end if
		
		if (sys%do_nonrad) then
			write (*, *)
			write (*, *) 'NON-RADIATIVE RATE CALCULATION'
			write (*, '(1x,a,i2,a,a)') 'Read in ', sys%nlevels, ' energy levels from file ', arg
			write (*, *)
			if ((sys%algorithm .ne. 'brute') .and. (sys%algorithm .ne. 'stochastic')) sys%algorithm = 'hybrid'
			call sys%fc_compute
			call sys%cutoff_compute
			noccs = 1
			
			allocate(minoccs(sys%nlevels))
			allocate(maxoccs(sys%nlevels))
			minoccs = 0
			maxoccs = [(sys%cutoffs(ix, sys%nthresh+1), ix=1,sys%nlevels)]
			
			emin = sys%e_target * TO_EV - sys%delta_e
			emax = emin + 2*sys%delta_e
			noccs = 0
			write (*, '(1x,a,f10.6,1x,a,1x,f10.6,1x,a)') 'Searching for occs with E(n) between', emin, 'and', emax, 'eV'
			write (*, '(1x,a,a,1x,a,i2)') 'Using algorithm: ', trim(sys%algorithm), 'with ncut = ', sys%ncut
			select case (sys%algorithm)
			case ('brute')
				call brute_force(sys%nlevels, sys%energies, maxoccs, minoccs, emax, emin, occlist, enlist, noccs)
			case ('stochastic')
				maxoccs = [(min(n_cut_guess, sys%cutoffs(ix, sys%nthresh+1)), ix=1,sys%nlevels)]
				allocate(guesses(n_guesses, sys%nlevels))
				write(*, '(1x,a,i2)') "Generating guesses with ncut = ", n_cut_guess
				do ix=1,n_guesses
					target_en = emin + ((real(ix)-0.5d0)/real(n_guesses)) * (emax - emin)
					call knap_n(target_en, sys%nlevels, sys%energies, maxoccs, 0.05d0, guesses(ix, :), guess_ens(ix))
					write (*, *) guesses(ix, :), guess_ens(ix)
				end do
				maxoccs = [(sys%cutoffs(ix, sys%nthresh+1), ix=1,sys%nlevels)]
				write(*, *) "Sampling"
				allocate(occlist(sys%nlevels, sys%nsamples)) 
				allocate(enlist(sys%nsamples))
				select case (sys%weighting)
				case ('fc')
					call do_sample(sys, guesses, maxoccs, guess_ens, emin, emax, weighted_indices, tbl, occlist, enlist, noccs)
				case default
					call do_sample(sys, guesses, maxoccs, guess_ens, emin, emax, uniform_indices, tbl, occlist, enlist, noccs)
				end select
				noccs = noccs + 1
			case default
				worst_case = ncombinations(sys%nlevels, maxoccs, minoccs) 
				call screened_brute_force(sys, emax, emin, occlist, enlist, noccs, worst_case)
			end select
			noccs = noccs - 1
			write (*, '(1x,a,1x,i10,1x,a)') 'Found', noccs, 'occs'
			write (*, *)
			write (*, *) 'Calculating non-radiative rate'
			call sys%calculate_kic(noccs, occlist, init=.true.)
			write (*, '(1x,a,1x,e14.8,1x,a)') 'k_ic =', sys%k_ic, 's-1'
			write (*, *)
		end if
		
		if (sys%do_rad .and. sys%do_nonrad) then
			write (*, '(1x,a,1x,f12.8,1x,a)') 'PLQY =', 100d0 * sys%k_r/(sys%k_r + sys%k_ic), '%'
		end if
		
		call cpu_time(finish)
		write (*, *) 'Time elapsed:', format_time(finish-start)
	else
		write (*, *) 'Could not open file', arg
	end if
end program main
