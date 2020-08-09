module knapsack
	use constants, only : dbl, bigint, print_frequency
	use ioutil, only : ncombinations
	use progress, only : progress_bar_time
	use nonradiative
	implicit none
contains
	subroutine knap_01(capacity, n, values, res)
		integer, intent(in)					:: n
		real(dbl), intent(in) 				:: capacity
		real(dbl), dimension(n), intent(in)	:: values
		integer(smallint), dimension(n), intent(out)  :: res
		
		integer									:: maxcap, w, ix, jx, err
		real(dbl)								:: wt
		logical									:: was_added
		real(dbl), allocatable, dimension(:, :)	:: table
		
		maxcap = ceiling(capacity)+1
		allocate(table(maxcap, n+1), stat=err)
		res = 0
		table = 0d0
		
		do jx=2,n+1
			wt = values(jx-1)
			do ix=2,maxcap
				if (wt > ix-1) then
					table(ix, jx) = table(ix, jx-1)
				else
					table(ix, jx) = max(table(ix, jx-1), table(ix-int(wt), jx-1) + wt)
				end if
			end do
		end do
		
		w = maxcap
		do jx=n+1,2,-1
			was_added = (table(w, jx) .ne. table(w, jx-1))
			if (was_added) then
				res(jx-1) = 1
				wt = values(jx-1)
				w = w - int(wt)
			end if
			if (w .le. 0d0) then
				exit
			end if
		end do
		
		deallocate(table)
	end subroutine knap_01
	
	subroutine fptas(capacity, n, values, epsilon, res, sumres)
		real(dbl), intent(in)				:: capacity, epsilon
		integer, intent(in)					:: n
		real(dbl), dimension(n), intent(in) :: values
		integer(smallint), dimension(n), intent(out)  :: res
		real(dbl), intent(out)				:: sumres
		
		real(dbl)				:: K
		real(dbl), dimension(n)	:: new_values
		integer					:: ix
		
		K = epsilon * maxval(values) / real(n)
		new_values = values / K
		
		call knap_01(capacity/K, n, new_values, res)
		sumres = dot_product(values, res)
	end subroutine fptas
	
	subroutine expand(n, energies, occs, res)
		integer, intent(in)								  :: n
		real(dbl), dimension(n), intent(in) 			  :: energies
		integer(smallint), dimension(n), intent(in)  				  :: occs
		real(dbl), allocatable, dimension(:), intent(out) :: res
		
		integer	:: ix, jx, ctr, ntot
		ntot = sum(occs)
		allocate(res(ntot))
		
		ctr = 1
		do ix=1,n
			do jx=1,occs(ix)
				res(ctr) = energies(ix)
				ctr = ctr + 1
			end do
		end do
	end subroutine expand
	
	subroutine contract(nbase, ntot, base, total, res)
		integer, intent(in) 			  	   :: nbase, ntot
		integer(smallint), dimension(nbase), intent(in)  :: base
		integer(smallint), dimension(ntot), intent(in)   :: total
		integer(smallint), dimension(nbase), intent(out) :: res
		
		integer	:: ctr, ix, jx
		ctr = 1
		do ix=1,nbase
			res(ix) = sum(total(ctr:ctr+base(ix)-1))
			ctr = ctr + base(ix)
		end do
	end subroutine contract
	
	subroutine knap_n(capacity, n, energies, occs, epsilon, res, sumres)
		integer, intent(in)					:: n
		real(dbl), intent(in) 				:: capacity, epsilon
		real(dbl), dimension(n), intent(in)	:: energies
		integer(smallint), dimension(n), intent(in)	:: occs
		integer(smallint), dimension(n), intent(out)  :: res
		real(dbl), intent(out)				:: sumres
		
		integer								 :: ntot
		integer(smallint), allocatable, dimension(:) 	 :: tempres
		real(dbl), allocatable, dimension(:) :: totvalues
		
		call expand(n, energies, occs, totvalues)
		ntot = size(totvalues)
		allocate(tempres(ntot))
		call fptas(capacity, ntot, totvalues, epsilon, tempres, sumres)
		call contract(n, ntot, occs, tempres, res)
		deallocate(tempres)
		deallocate(totvalues)
	end subroutine knap_n
	
	recursive subroutine iterate(sys, noccs, emax, emin, occ, ix, enlist, occix, checked, maxnocc)
		type(sysdata), intent(inout)							:: sys
		integer, intent(in)										:: ix
		integer(bigint), intent(in)								:: noccs, maxnocc
		integer(bigint), intent(inout)							:: occix, checked
		integer(smallint), dimension(sys%nlevels), intent(inout)	:: occ
		real(dbl), dimension(noccs), intent(inout)				:: enlist
		real(dbl), intent(in)									:: emax, emin
		
		integer		:: i, n
		real(dbl) 	:: en
		
		n = sys%nlevels
		do i=sys%minoccs(ix),sys%maxoccs(ix)
			occ(ix) = i
			
			if (ix .lt. n) then
				call iterate(sys, noccs, emax, emin, occ, ix+1, enlist, occix, checked, maxnocc)
			else
				en = dot_product(sys%energies, occ)
				if ((en .lt. emax) .and. (en .gt. emin)) then
					sys%mm%current_block(1:n, occix) = occ(1:n)
					enlist(occix) = en
					occix = occix+1
					if (modulo(occix, 100000) .eq. 0) write(*, *) occix
					if (occix .gt. sys%mm%chunk_size) then
						write(*, '(1x,a,1x,i10)') 'Finished block', sys%mm%current_record
						call sys%calculate_kic(sys%mm%chunk_size, sys%mm%current_block, init=.false., stopix=sys%mm%chunk_size)
						if (sys%do_write) then
							call sys%mm%block_swap(sys%energies, sys%hrfactors, sys%mm%chunk_size)
						else
							sys%mm%current_record = sys%mm%current_record + 1
							sys%mm%current_block = 0
						end if
						occix = 1
					end if
				end if
				checked = checked + 1
				if (mod(checked, print_frequency) .eq. 0) write(*, *) 'Checked ', checked
			end if
		end do		
	end subroutine iterate
	
	subroutine brute_force(sys, emax, emin, enlist, noccs)
		type(sysdata), intent(inout)						:: sys
		real(dbl), intent(in)								:: emax, emin
		real(dbl), allocatable, dimension(:), intent(out)	:: enlist
		integer(bigint), intent(out)						:: noccs
		
		integer(bigint) :: max_nocc, occsize, checked
		integer(smallint)	:: cocc(sys%nlevels)
		integer :: ix
		cocc = 0
		cocc(1) = max(sys%minoccs(1), 1)
		max_nocc = sys%maxnoccs
		occsize = sys%mm%chunk_size
		
		allocate(enlist(occsize))
		noccs = 1
		checked = 0
		write(*, '(1x,a,1x,e12.4,1x,a)') 'Estimated', real(max_nocc), 'occupations to check'
		call iterate(sys, occsize, emax, emin, cocc, 1, enlist, noccs, checked, max_nocc)
		!call progress_bar_time(max_nocc, max_nocc)
		write(*, *)
	end subroutine brute_force
	
	subroutine screened_brute_force(sys, emax, emin, enlist, noccs, worst_case)
		type(sysdata), intent(inout)						:: sys
		real(dbl), intent(in)								:: emax, emin
		real(dbl), allocatable, dimension(:), intent(out)	:: enlist
		integer(bigint), intent(in)							:: worst_case
		integer(bigint), intent(out)						:: noccs
		
		integer(bigint)	:: max_nocc, occsize, checked
		integer(smallint)	:: cocc(sys%nlevels)
		integer :: ix
		cocc = 0
		cocc(1) = sys%cutoffs(1, sys%nthresh+1)
		max_nocc = sys%maxnoccs
		occsize = sys%mm%chunk_size
		
		allocate(enlist(occsize))
		noccs = 1
		checked = 0
		write(*, '(1x,a,1x,e12.4)') 'Worst case no. of occs to check:', real(worst_case)
		call screened_iterate(sys, occsize, emax, emin, cocc, 1, enlist, noccs, checked, max_nocc, sys%maxnfix, fixedn=.false.)
		!call progress_bar_time(max_nocc, max_nocc)
		write(*, *)
	end subroutine screened_brute_force
	
	subroutine screened_fixed_n(sys, emax, emin, enlist, noccs, worst_case)
		type(sysdata), intent(inout)						:: sys
		real(dbl), intent(in)								:: emax, emin
		real(dbl), allocatable, dimension(:), intent(out)	:: enlist
		integer(bigint), intent(in)							:: worst_case
		integer(bigint), intent(out)						:: noccs
		
		integer(bigint)	:: max_nocc, occsize, checked
		integer(smallint)	:: cocc(sys%nlevels)
		integer :: ix, nx
		max_nocc = sys%maxnoccs
		occsize = sys%mm%chunk_size
		
		allocate(enlist(occsize))
		noccs = 1
		checked = 0
		write(*, '(1x,a,1x,i5,a,1x,i5)') 'Nmin =', sys%minnfix, ', Nmax =', sys%maxnfix
		write(*, '(1x,a,1x,e12.4)') 'Worst case no. of occs to check:', real(worst_case)
		do nx=sys%minnfix, sys%maxnfix
			cocc = 0
			cocc(1) = sys%cutoffs(1, sys%nthresh+1)
			write(*, '(1x,a,1x,i3,1x,a,1x,e12.4,1x,a)') 'Doing N =', nx, 'found', real(noccs-1), 'valid occs so far'
			call screened_iterate(sys, occsize, emax, emin, cocc, 1, enlist, noccs, checked, max_nocc, nfix=nx, fixedn=.true.)
		end do
		!call progress_bar_time(max_nocc, max_nocc)
		write(*, *)
	end subroutine
	
	recursive subroutine screened_iterate(sys, noccs, emax, emin, occ, ix, enlist, occix, checked, maxnocc, nfix, fixedn)
		type(sysdata), intent(inout)								:: sys
		integer, intent(in)											:: ix, nfix
		integer(bigint), intent(in)									:: noccs, maxnocc
		integer(bigint), intent(inout)								:: occix, checked
		integer(smallint), dimension(sys%nlevels), intent(inout)	:: occ
		real(dbl), dimension(noccs), intent(inout)					:: enlist
		real(dbl), intent(in)										:: emax, emin
		logical, intent(in)											:: fixedn
		
		integer		:: i, maxix, minix, bnd, newnfix
		real(dbl) 	:: en
		maxix = min(occ(ix), nfix)
		if (ix .eq. sys%nlevels) then
			if (fixedn) then
				minix = max(nfix, 0)
				if (maxix .lt. minix) then
					continue
				end if
				maxix = minix
			else
				minix = 0
			end if
			do i=minix,maxix
				occ(ix) = i
				en = dot_product(sys%energies, occ)
				if ((en .lt. emax) .and. (en .gt. emin)) then
					sys%mm%current_block(:, occix) = occ(:)
					enlist(occix) = en
					occix = occix+1
					if (occix .gt. sys%mm%chunk_size) then
						write(*, '(1x,a,1x,i10)') 'Finished block', sys%mm%current_record
						call sys%calculate_kic(sys%mm%chunk_size, sys%mm%current_block, init=.false., stopix=sys%mm%chunk_size)
						if (sys%do_write) then
							call sys%mm%block_swap(sys%energies, sys%hrfactors, sys%mm%chunk_size)
						else
							sys%mm%current_record = sys%mm%current_record + 1
							sys%mm%current_block = 0
						end if
						occix = 1
					end if
				end if
				checked = checked + 1
				if (mod(checked, print_frequency) .eq. 0) write(*, *) 'Checked ', checked
			end do
		else
			do i=0,maxix
				occ(ix) = i
				newnfix = nfix - i
				call sys%get_next_bound(occ, ix, bnd)
				occ(ix+1) = sys%cutoffs(ix+1, bnd)
				call screened_iterate(sys, noccs, emax, emin, occ, ix+1, enlist, occix, checked, maxnocc, newnfix, fixedn)
			end do
		end if
	end subroutine screened_iterate
	
end module knapsack		
		