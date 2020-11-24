module knapsack
	use constants, only : dbl, bigint, print_frequency
	use ioutil, only : ncombinations, sort_list
	use sample, only : differences
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
	
	subroutine adjust(capacity, n, energies, epsilon, nmax, res, sumres)
		integer, intent(in)								:: n, nmax
		real(dbl), intent(in) 							:: capacity, epsilon
		real(dbl), dimension(n), intent(in)				:: energies
		integer(smallint), dimension(n), intent(inout)  :: res
		real(dbl), intent(out)							:: sumres
		
		integer	:: ntot, ndelta, ix, jx, best_ix, best_jx, iter
		integer(sint) 	  :: indices(n)
		integer(smallint) :: res_sorted(n)
		logical :: empty
		real(dbl) :: ens_sorted(n), edelta, en_diffs(n, n), best_edelta, de
		ntot = sum(res)
		ndelta = nmax - ntot
		
		if (ndelta .ne. 0) then
			call sort_list(n, energies, ens_sorted, indices)
			call differences(n, ens_sorted, en_diffs)
			do ix=1,n
				res_sorted(ix) = res(indices(ix))
			end do
			
			res_sorted(1) = res_sorted(1) + ndelta
			neg: do ix=1,n-1
				if (res_sorted(ix) .lt. 0) then
					res_sorted(ix+1) = res_sorted(ix+1) + res_sorted(ix)
					res_sorted(ix) = 0
				else 
					exit neg
				end if
			end do neg
			res_sorted(n) = max(0, res_sorted(n))
			
			edelta = dot_product(ens_sorted, res_sorted)
			edelta = capacity - edelta
			iter = 0
			do while ((abs(edelta) .gt. epsilon) .and. (iter < 100))
				best_edelta = 1d0
				do ix=1,n-1
					do jx=ix+1,n
						if (de*edelta .gt. 0) then 
							empty = res_sorted(ix) < 1
						else
							empty = res_sorted(jx) < 1
						end if 
						if (.not. empty) then
							de = abs(en_diffs(ix, jx) - edelta)
							if (de < best_edelta) then
								best_edelta = de
								best_ix = ix
								best_jx = jx
							end if
						end if
					end do
				end do
				if (de * edelta .gt. 0) then
					res_sorted(best_ix) = res_sorted(best_ix) - 1
					res_sorted(best_jx) = res_sorted(best_jx) + 1
				else 
					res_sorted(best_ix) = res_sorted(best_ix) + 1
					res_sorted(best_jx) = res_sorted(best_jx) - 1
				end if
				edelta = capacity - dot_product(ens_sorted, res_sorted)
				iter = iter + 1
			end do
			
			do ix=1,n
				res(indices(ix)) = res_sorted(ix)
			end do
			sumres = dot_product(res, energies)
		end if
	end subroutine adjust
	
	subroutine knap_n(capacity, n, energies, occs, epsilon, nmax, res, sumres, with_adjust)
		integer, intent(in)								:: n, nmax
		real(dbl), intent(in) 							:: capacity, epsilon
		real(dbl), dimension(n), intent(in)				:: energies
		integer(smallint), dimension(n), intent(in)		:: occs
		integer(smallint), dimension(n), intent(out)  	:: res
		real(dbl), intent(out)							:: sumres
		logical, intent(in)								:: with_adjust
		
		integer								 			:: ntot
		integer(smallint), allocatable, dimension(:) 	:: tempres
		real(dbl), allocatable, dimension(:) 			:: totvalues
		
		call expand(n, energies, occs, totvalues)
		ntot = size(totvalues)
		allocate(tempres(ntot))
		call fptas(capacity, ntot, totvalues, epsilon, tempres, sumres)
		call contract(n, ntot, occs, tempres, res)
		
		if (with_adjust) call adjust(capacity, n, energies, epsilon, nmax, res, sumres)
		
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
		
		integer		:: i, n, tmp
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
							call sys%mm%block_swap(sys%energies, sys%hrfactors, sys%mm%chunk_size, tmp)
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
		
		integer		:: i, maxix, minix, bnd, newnfix, tmp
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
							call sys%mm%block_swap(sys%energies, sys%hrfactors, sys%mm%chunk_size, tmp)
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
	
	subroutine lambda(sys, emax, emin, klag) 
		type(sysdata), intent(inout)	:: sys
		real(dbl), intent(in)			:: emax, emin
		real(dbl), intent(out)			:: klag
		
		real(dbl)	:: emid, lam, deltae, en, d1, d2, step, de, predde
		real(dbl)	:: rho, tau, trust, p, maxtrust, bestlam, beste
		real(dbl), dimension(3) :: derivs
		real(dbl), dimension(sys%nlevels) :: occs
		integer		:: ix

		emid = 0.5 * (emax + emin)
		
		bestlam = 0.0
		lam = 16.0
		beste = emid
		do while (beste .gt. 0.5d0)
			lam = lam - 0.5
			call lagrange(sys%nlevels, sys%hrfactors, sys%energies, lam, derivs)
			en = abs(emid - derivs(1))
			if (en .lt. beste) then
				beste = en
				bestlam = lam
			end if
		end do
		write(*, *)
		write(*, *) 'LAGRANGE MULTIPLIER RESULTS'
		write(*, '(1x,a15,1x,f10.4)') 'Starting lambda:', bestlam
		
		deltae = 1d0
		p = bestlam
		ix = 0
		maxtrust = 1.0d0
		trust = 0.2d0
		do while ((deltae .gt. 0.01d0) .and. (ix .lt. 50))
			lam = p
			call lagrange(sys%nlevels, sys%hrfactors, sys%energies, lam, derivs)
			en = emid - derivs(1)
			
			! take Newton step
			d1 = -2d0 * en * derivs(2)
			d2 = 2d0 * derivs(3)**2 + d1 * derivs(2)
			
			tau = 1d0
			if (d2 .gt. 0) then
				tau = min(abs(d1)/ (trust * d2), 1d0)
			end if 
			step = -tau * trust * d1 / abs(d1)
			
			p = lam + step
			call lagrange(sys%nlevels, sys%hrfactors, sys%energies, p, derivs)
			de = (en - emid + derivs(1))
			predde = -(step*d1 + 0.5*step*d2*step)
			rho = abs(de/predde)
			if (rho .lt. 0.25) then
				trust = 0.25*trust
			else if ((rho .gt. 0.75) .and. ((abs(step)-trust) .lt. 1E-3)) then
				trust = min(2.0*trust, maxtrust)
			else
				trust = trust 
			end if
			
			if (rho .lt. 0.1) p = lam
			
			
			deltae = abs(en)
			ix = ix + 1
		end do
		
		d1 = d1 / abs(d1)
		bestlam = lam
		beste = abs(en)
		de = beste
		deltae = -1d0
		do while (deltae .lt. 0)
			lam = lam - d1*0.005
			call lagrange(sys%nlevels, sys%hrfactors, sys%energies, lam, derivs)
			en = abs(emid - derivs(1))
			if (en .lt. beste) then
				beste = en
				bestlam = lam
			end if
			deltae = en - de
			de = en
		end do
		write(*, '(1x,a15,1x,e12.4,1x,a)') 'Best DE:', beste, 'eV'
		write(*, '(1x,a15,1x,f10.4)') 'Best Lambda:', bestlam
		
		do ix=1,sys%nlevels
			occs(ix) = sys%hrfactors(ix) * exp(-sys%energies(ix) * bestlam)
		end do
		write(*, '(1x,a15,1x,f8.2)') 'Total occ:', sum(occs)
		sys%lambda_n = nint(sum(occs))
		write(*, '(1x,a15,1x,f10.6,1x,a)') 'Total en:', dot_product(occs, sys%energies), 'eV'
		call sys%calculate_from_lagrange(occs, klag)
	end subroutine lambda
	
	subroutine lagrange(n, hrfactors, levels, lam, results)
		integer, intent(in)						:: n
		real(dbl), dimension(n), intent(in) 	:: levels, hrfactors
		real(dbl), intent(in)					:: lam
		real(dbl), dimension(3), intent(out)	:: results
	
		integer		:: ix
		real(dbl) 	:: tmp
		results = 0d0
		do ix=1,n
			tmp = hrfactors(ix) * exp(-levels(ix) * lam) * levels(ix)
			results(1) = results(1) + tmp
			results(2) = results(2) - tmp*levels(ix)
			results(3) = results(3) + tmp*(levels(ix)**2)
		end do
	end subroutine lagrange
	
end module knapsack		
		
