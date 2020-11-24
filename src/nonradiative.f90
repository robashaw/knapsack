module nonradiative
	use constants
	use ioutil
	implicit none
	
	type sysdata
		logical											:: do_rad, do_nonrad, do_write, skip_rate
		integer											:: ncut, nthresh, nlevels, nsamples, natoms, debug_level
		integer											:: maxnfix, minnfix, nrecords, lambda_n
		integer(bigint)									:: maxnoccs
		character(len=100)								:: bfile, gradfile, hessfile, algorithm, weighting
		character(len=100)								:: radfile, calctype, radunits, sortby, occprefix
		real(dbl)										:: k_ic, k_ic_ht, k_r, e_target, delta_e, gamma, tdm, memory, user_gamma
		real(dbl), dimension(:), allocatable			:: energies, hrfactors, masses, V_vq_j
		real(dbl), dimension(:, :), allocatable			:: fcfactors, V_jj, Bvqj
		integer, dimension(:, :), allocatable			:: cutoffs, bounds
		integer(sint), dimension(:), allocatable		:: energy_order
		integer(smallint), dimension(:), allocatable	:: minoccs, maxoccs
		type(memorymanager)								:: mm
	contains
		procedure	:: from_file => sysdata_from_file
		procedure	:: init => sysdata_init
		procedure	:: fc_compute => sysdata_fc_compute
		procedure	:: cutoff_compute => sysdata_cutoff_compute
		procedure	:: compute_zn => sysdata_compute_zn
		procedure	:: build_V => sysdata_build_V
		procedure	:: calculate_kic => sysdata_calc_kic
		procedure	:: calculate_from_lagrange => sysdata_calc_lagrange
		procedure	:: calculate_from_file => sysdata_calc_from_file
		procedure	:: calculate_gamma => sysdata_calc_gamma
		procedure	:: get_next_bound => sysdata_get_next_bound
		procedure	:: screened_ncombinations => sysdata_ncombinations
		procedure	:: free => sysdata_free 
	end type sysdata
	
contains
	
	subroutine sysdata_from_file(sys, filename)
		class(sysdata), intent(inout)	:: sys 
		character(len=*), intent(in)	:: filename
		
	    ! Input related variables
	    character(len=100) :: buffer, label, randostring
	    integer 	:: pos
	    integer 	:: ios = 0
	    integer 	:: line = 0
		integer 	:: ix
		real(dbl)	:: tmp_nsamples = 0d0

	    open(main_input_unit, file=filename)

	    ! ios is negative if an end of record condition is encountered or if
	    ! an endfile condition was detected.  It is positive if an error was
	    ! detected.  ios is zero otherwise.

		sys%debug_level = 0
		sys%do_write = .false.
		sys%skip_rate = .false.
		sys%memory = 0.5d0
		sys%nrecords = 0
		sys%maxnfix = -1
		sys%minnfix = -1
		sys%occprefix = 'occs'
		sys%hessfile  = 'none'
		sys%user_gamma=-1.0
	    do while (ios == 0)
	       read(main_input_unit, '(A)', iostat=ios) buffer
	       if (ios == 0) then
	          line = line + 1

	          ! Find the first instance of whitespace.  Split label and data.
	          pos = scan(buffer, ',')
	          label = trim(buffer(1:pos-1))
	          buffer = buffer(pos+1:)

	          select case (label)
			  case ('spectrum')
		  		 read(buffer, *, iostat=ios) sys%radfile
				 sys%radfile = trim(adjustl(sys%radfile))
			  case ('bfile')
			  	 read(buffer, *, iostat=ios) sys%bfile
				 sys%bfile = trim(adjustl(sys%bfile))
   			  case ('gradfile')
   			  	 read(buffer, *, iostat=ios) sys%gradfile
   				 sys%gradfile = trim(adjustl(sys%gradfile))
      		  case ('hessfile')
      		  	 read(buffer, *, iostat=ios) sys%hessfile
      			 sys%hessfile = trim(adjustl(sys%hessfile))
      		  case ('algorithm')
      			 read(buffer, *, iostat=ios) sys%algorithm
      			 sys%algorithm = trim(adjustl(sys%algorithm))
         	  case ('weighting')
         		 read(buffer, *, iostat=ios) sys%weighting
         		 sys%weighting = trim(adjustl(sys%weighting))
              case ('occprefix')
            	 read(buffer, *, iostat=ios) sys%occprefix
            	 sys%occprefix = trim(adjustl(sys%occprefix))
			  case ('calculation')
			  	 read(buffer, *, iostat=ios) sys%calctype
				 sys%calctype = trim(adjustl(sys%calctype))
   			  case ('radunits')
   			  	 read(buffer, *, iostat=ios) sys%radunits
   				 sys%radunits = trim(adjustl(sys%radunits))
      		  case ('sortby')
      			 read(buffer, *, iostat=ios) sys%sortby
      		     sys%sortby = trim(adjustl(sys%sortby))
   	          case ('debug')
   	             read(buffer, *, iostat=ios) sys%debug_level
			  case ('nrecords')
			  	 read(buffer, *, iostat=ios) sys%nrecords
			  case ('dump')
			  	 read(buffer, *, iostat=ios) sys%mm%dumplevel
	          case ('ncut')
	             read(buffer, *, iostat=ios) sys%ncut
			  case ('tofile')
			  	 read(buffer, *, iostat=ios) randostring
				 if (trim(adjustl(randostring)) .eq. 'true') then
					 sys%do_write = .true.
				 end if
   			  case ('norate')
   			  	 read(buffer, *, iostat=ios) randostring
   				 if (trim(adjustl(randostring)) .eq. 'true') then
   					 sys%skip_rate = .true.
   				 end if
   	          case ('natoms')
   	             read(buffer, *, iostat=ios) sys%natoms
	          case ('nthresh')
	             read(buffer, *, iostat=ios) sys%nthresh
			  case ('nsamples')
			  	 read(buffer, *, iostat=ios) tmp_nsamples
				 sys%nsamples = int(tmp_nsamples)
			  case ('nmax')
			  	 read(buffer, *, iostat=ios) sys%maxnfix
   			  case ('nmin')
   			  	 read(buffer, *, iostat=ios) sys%minnfix
      		  case ('memory')
      			 read(buffer, *, iostat=ios) sys%memory
   			  case ('energy')
   			  	 read(buffer, *, iostat=ios) sys%e_target
   			  case ('deltae')
   			  	 read(buffer, *, iostat=ios) sys%delta_e
			  case ('tdm')
			  	 read(buffer, *, iostat=ios) sys%tdm
      		  case ('gamma')
      			 read(buffer, *, iostat=ios) sys%user_gamma
			  case ('levels')
			  	 read(buffer, *, iostat=ios) sys%nlevels
				 allocate(sys%energies(sys%nlevels))
				 allocate(sys%hrfactors(sys%nlevels))
				 do ix=1,sys%nlevels
					 read(main_input_unit, '(A)', iostat=ios) buffer
					 pos = scan(buffer, ',')
					 label = buffer(1:pos-1)
					 buffer = buffer(pos+1:)
					 read(label, *, iostat=ios) sys%energies(ix)
					 read(buffer, *, iostat=ios) sys%hrfactors(ix)
				 end do
	          case default
	             print *, 'Skipping invalid label at line', line
	          end select
	       end if
	    end do
		close(main_input_unit)
	end subroutine sysdata_from_file
	
	subroutine sysdata_init(sys)
		class(sysdata), intent(inout)	:: sys
		
		integer	:: ix, tmp
		real(dbl)	:: nestimate
		real(dbl), dimension(sys%nlevels) :: unsorted_levels
		
		allocate(sys%energy_order(sys%nlevels))
		do ix=1,sys%nlevels
			sys%energy_order(ix) = ix
		end do
		
		if (sys%sortby .eq. 'energy') then
			unsorted_levels(:) = sys%energies(:)
			call sort_list(sys%nlevels, unsorted_levels, sys%energies, sys%energy_order)
			call reorder_list(sys%nlevels, sys%hrfactors, sys%energy_order)
		else
			! reverse sort by HR factor by default
			unsorted_levels(:) = -sys%hrfactors(:)
			call sort_list(sys%nlevels, unsorted_levels, sys%hrfactors, sys%energy_order)
			sys%hrfactors = -sys%hrfactors
			call reorder_list(sys%nlevels, sys%energies, sys%energy_order)
		end if
		
		! calculate maximum no. of occs and initialise memory manager
		if (sys%do_nonrad) then
			call sys%fc_compute
			call sys%cutoff_compute
			
			allocate(sys%minoccs(sys%nlevels))
			allocate(sys%maxoccs(sys%nlevels))
			sys%minoccs = 0
			sys%maxoccs = [(sys%cutoffs(ix, sys%nthresh+1), ix=1,sys%nlevels)]
			sys%maxnoccs = ncombinations(sys%nlevels, sys%maxoccs, sys%minoccs)
			
			if (sys%maxnfix .lt. 0) then
				sys%maxnfix = sum(sys%maxoccs)
				if (sys%algorithm .eq. 'fixedn') then
					! Estimate maxn 
					nestimate = 10**(2*sys%nthresh/sys%nlevels)
					nestimate = nestimate * exp(-sum(sys%hrfactors))
					nestimate = nestimate ** (1.0/sys%nlevels)
					sys%maxnfix = ceiling(sys%nlevels * nestimate)
				end if
			end if
			
			if (sys%minnfix .lt. 0) then
				sys%minnfix = ceiling((sys%e_target*TO_EV - sys%delta_e)/(maxval(sys%energies)))
			end if
			
			if (sys%algorithm .eq. 'stochastic') then
				sys%mm%notunique = .true.
				sys%memory = 0.5 * sys%memory
			end if
			call sys%mm%initialise(sys%nlevels, sys%maxnoccs, sys%memory, sys%nthresh, sys%occprefix)
			call sys%mm%block_swap(sys%energies, sys%hrfactors, sys%mm%chunk_size, tmp)
		end if
	end subroutine sysdata_init
	
	subroutine sysdata_fc_compute(sys)
		class(sysdata), intent(inout)	:: sys
		
		integer		:: ix, nx
		real(dbl)	:: tmp, yk
		
		if (allocated(sys%fcfactors)) deallocate(sys%fcfactors)
		allocate(sys%fcfactors(sys%nlevels, sys%ncut+1))
		
		do ix=1, sys%nlevels
			yk = sys%hrfactors(ix)
			tmp = exp(-yk)
			do nx=0, sys%ncut
				sys%fcfactors(ix, nx+1) = sqrt(tmp)
				tmp = tmp*yk/real(nx+1)
			end do
		end do
	end subroutine sysdata_fc_compute
	
	subroutine sysdata_cutoff_compute(sys)
		class(sysdata), intent(inout)	:: sys
		
		integer		:: ix, nx, jx, tmpix
		real(dbl)	:: tmp
		logical		:: found
		
		if (allocated(sys%cutoffs)) deallocate(sys%cutoffs)
		allocate(sys%cutoffs(sys%nlevels, sys%nthresh+1))
		sys%cutoffs = 1
		
		do ix=1,sys%nlevels
			tmpix = 1
			tmp = log10(sys%fcfactors(ix, tmpix))
			do nx=0,sys%nthresh
				found = .false.
				do while (.not. found)
					if ((tmp .lt. -nx) .or. (tmpix .eq. sys%ncut)) then
						sys%cutoffs(ix, nx+1) = tmpix
						found = .true.
					else
						tmpix = tmpix+1
						tmp = log10(sys%fcfactors(ix, tmpix))
					end if
				end do
			end do
		end do
		
		if (allocated(sys%bounds)) deallocate(sys%bounds)
		allocate(sys%bounds(sys%nlevels, sys%ncut+1))
		sys%bounds = sys%nthresh+1
		do ix=1,sys%nlevels
			tmpix = 1
			inner: do nx=1,sys%nthresh+1
				do jx=tmpix,sys%ncut+1
					sys%bounds(ix, jx) = nx
				end do
				tmpix = sys%cutoffs(ix, nx)+1
				if (tmpix .eq. sys%ncut+1) exit inner
			end do inner
		end do	
	
	end subroutine sysdata_cutoff_compute
	
	subroutine sysdata_free(sys)
		class(sysdata), intent(inout)	:: sys
		if (allocated(sys%energies)) deallocate(sys%energies)
		if (allocated(sys%hrfactors)) deallocate(sys%hrfactors)
		if (allocated(sys%masses)) deallocate(sys%masses)
		if (allocated(sys%fcfactors)) deallocate(sys%fcfactors)
		if (allocated(sys%V_vq_j)) deallocate(sys%V_vq_j)
		if (allocated(sys%Bvqj)) deallocate(sys%Bvqj)
		if (allocated(sys%V_jj)) deallocate(sys%V_jj)
		if (allocated(sys%cutoffs)) deallocate(sys%cutoffs)
		if (allocated(sys%bounds)) deallocate(sys%bounds)
		if (allocated(sys%energy_order)) deallocate(sys%energy_order)
		if (allocated(sys%maxoccs)) deallocate(sys%maxoccs)
		if (allocated(sys%minoccs)) deallocate(sys%minoccs)
	end subroutine sysdata_free
	
	real(dbl) function	compute_kfcn(sys, occs) result(kfcn)
		class(sysdata), intent(in) 					:: sys
		integer(smallint), dimension(sys%nlevels), intent(in)	:: occs
		
		integer :: ix
		kfcn = 1d0
		do ix = 1, sys%nlevels
			kfcn = kfcn * sys%fcfactors(ix, occs(ix)+1)
		end do
	end function compute_kfcn
	
	subroutine sysdata_compute_zn(sys, occs, res)
		class(sysdata), intent(inout)					:: sys
		integer(smallint), dimension(sys%nlevels), intent(in)		:: occs
		real(dbl), dimension(sys%nlevels), intent(out)	:: res
		
		integer		:: jx
		real(dbl)	:: kfcn, tmp
		
		kfcn = compute_kfcn(sys, occs)
		do jx=1, sys%nlevels
			res(jx) = 0.5d0 * sys%energies(jx) / sys%hrfactors(jx)
			tmp = real(occs(jx)) - sys%hrfactors(jx)
			res(jx) = sqrt(res(jx) * tmp * tmp)
			res(jx) = res(jx) * kfcn
		end do
	end subroutine sysdata_compute_zn
	
	subroutine sysdata_build_V(sys, grads, hessian)
		class(sysdata), intent(inout)					 	 	:: sys
		real(dbl), dimension(sys%natoms, 3), intent(inout)  	:: grads
		real(dbl), dimension(:, :), allocatable, intent(inout)  :: hessian
		
		integer								:: b_ios, g_ios, vx, qx, dummy1, dummy2
		character(len=1)					:: dummy_char1, dummy_char2
		real(dbl), dimension(sys%nlevels)	:: tmp
		real(dbl), dimension(3*sys%natoms, sys%nlevels) :: tmphess
		
		if (allocated(sys%masses)) deallocate(sys%masses)
		allocate(sys%masses(sys%natoms))
		if (allocated(sys%V_vq_j)) deallocate(sys%V_vq_j)
		allocate(sys%V_vq_j(sys%nlevels))
		if (allocated(sys%Bvqj)) deallocate(sys%Bvqj)
		allocate(sys%Bvqj(3*sys%natoms, sys%nlevels))
		
		if (sys%skip_rate) then
			sys%V_vq_j = 1.0/sqrt(sys%gamma)
		else	
			open(gradfile_unit, file=sys%gradfile)
			readgrads: do vx=1,sys%natoms
				do qx=1,3
					read(gradfile_unit, *, iostat=g_ios) sys%masses(vx), grads(vx, qx)
					if (g_ios .ne. 0) then
						write(*, '(1x,a,1x,i4,1x,a,1x,i4)') 'Error reading grads, ierr=', g_ios, 'line=', vx
						exit readgrads
					end if
				end do
			end do readgrads
			close(gradfile_unit)
			write(*, '(1x,a,1x,a,1x,a,1x,i4,1x,a)') 'Read grads from', trim(adjustl(sys%gradfile)), 'with', vx-1, 'atoms'
			
			open(bfile_unit, file=sys%bfile)
			read(bfile_unit, *, iostat=b_ios)
			sys%V_vq_j = 0d0
			readbfile: do vx=1,sys%natoms
				do qx=1,3
					read(bfile_unit, *, iostat=b_ios) dummy_char1, dummy_char2, dummy1, dummy2, tmp(:)
					call reorder_list(sys%nlevels, tmp, sys%energy_order)
					sys%Bvqj(3*(vx-1)+qx, :) = tmp(:) / sqrt(sys%masses(vx))
					if (b_ios .ne. 0) then
						write(*, '(1x,a,1x,i4,1x,a,1x,i4)') 'Error reading B-vectors, ierr=', g_ios, 'line=', vx
						exit readbfile
					else
						sys%V_vq_j(:) = sys%V_vq_j(:) + grads(vx, qx) * tmp(:) / sys%masses(vx)
					end if
				end do
			end do readbfile
			close(bfile_unit)
			write(*, '(1x,a,1x,a,1x,a,1x,i4,1x,a)') 'Read B-vectors from', trim(adjustl(sys%bfile)), 'with', vx-1, 'atoms'
			
			if (sys%hessfile .ne. 'none') then
				if (allocated(hessian)) deallocate(hessian)
				allocate(hessian(3*sys%natoms, 3*sys%natoms))
				if (allocated(sys%V_jj)) deallocate(sys%V_jj)
				allocate(sys%V_jj(sys%nlevels, sys%nlevels))
				
				open(hessfile_unit, file=sys%hessfile)
				readhess: do vx=1,3*sys%natoms
					do qx=1,vx
						read(hessfile_unit, *, iostat=b_ios) hessian(vx, qx)
						if (b_ios .ne. 0) then
							write(*, '(1x,a,1x,i4,1x,a,1x,i4)') 'Error reading hessian, ierr=', b_ios, 'line=', vx
							exit readhess
						end if
						hessian(qx, vx) = hessian(vx, qx)
					end do 
				end do readhess
				close(hessfile_unit)
				write(*, '(1x,a,1x,a,1x,a,1x,i4,1x,a)') 'Read Hessian from', trim(adjustl(sys%hessfile)), 'with', (vx-1)/3, 'atoms'
				
				tmphess = matmul(hessian, sys%Bvqj)
				sys%V_jj = matmul(transpose(sys%Bvqj), tmphess)		
				
				sys%V_jj = sys%V_jj * V_CONVERT / TO_S
				
			end if
		end if
		
		sys%V_vq_j = sys%V_vq_j * V_CONVERT / TO_S
		if (sys%debug_level .gt. 1) then
			write(*, *) 'Elements of V_vq_j (FC)'
			write(*, '(1x,a4,1x,a10,1x,a20)') 'J', 'ENERGY', 'V_VQ_J'
			do qx=1,sys%nlevels
				write(*, '(1x,i4,1x,f10.6,1x,e20.8)') qx, sys%energies(qx), sys%V_vq_j(qx)
			end do
			
			if ((sys%debug_level .gt. 2) .and. (sys%hessfile .ne. 'none')) then
				write(*, *) 'Elements of V_ij (HT)'
				write(*, '(1x,a4,1x,a4,1x,a20)') 'I', 'J', 'V_VQ_J'
				do qx=1,sys%nlevels
					do vx=1,qx
						write(*, '(1x,i4,1x,i4,1x,e20.8)') qx, vx, sys%V_jj(qx, vx)
					end do
				end do
			end if
		end if
	end subroutine sysdata_build_V
	
	subroutine sysdata_calc_gamma(sys)
		class(sysdata), intent(inout)	:: sys
		
		integer		:: ix
		real(dbl)	:: be, w, tmp
		sys%gamma = 0d0
		do ix=1, sys%nlevels
			w  = sys%energies(ix)
			be = exp(w / (KB*300d0)) - 1d0
			be = 1d0/be
			w  = w / PLANCK
			tmp = w * w * sys%hrfactors(ix) * (2d0 * be + 1)
			sys%gamma = sys%gamma + tmp
		end do
	
		sys%gamma = SQRT_8LN2 * sqrt(sys%gamma)
	end subroutine sysdata_calc_gamma
	
	subroutine sysdata_calc_from_file(sys, record, prefix, ntot)
		class(sysdata), intent(inout)	:: sys
		integer, intent(in)				:: record
		character(len=*), intent(in)	:: prefix
		integer, intent(out)			:: ntot
		
		integer	:: nrows, nchunk, ix, startrow, endrow
		
		call sys%mm%get_nrows(sys%nlevels, record, prefix, nrows)
		write(*, *) 'Reading from file, expecting ', nrows, ' rows'
		startrow = 1
		endrow = min(sys%mm%chunk_size, nrows)-1
		ntot = 0
		main: do while (ntot .lt. nrows)
			nchunk = endrow - startrow + 1
			call sys%mm%read_from_bin(sys%nlevels, record, sys%mm%current_block, sys%energies, sys%hrfactors, prefix, startrow, endrow)
			call sys%calculate_kic(nchunk, sys%mm%current_block, init=.false., stopix=nchunk)
			ntot = ntot + nchunk
			startrow = endrow + 1
			endrow = startrow + min(nrows - ntot, sys%mm%chunk_size) - 1
		end do main
	end subroutine sysdata_calc_from_file
	
	subroutine sysdata_calc_kic(sys, noccs, occs, init, stopix)
		class(sysdata), intent(inout)						:: sys
		integer(bigint), intent(in)							:: noccs
		integer(smallint), dimension(sys%nlevels, noccs), intent(in)	:: occs
		logical, intent(in)									:: init
		integer, intent(in)									:: stopix
		
		integer			:: nx, ix
		integer(bigint)	:: counter(20)
		real(dbl)		:: z(sys%nlevels), res(sys%nlevels), tmp
		real(dbl)		:: grads(sys%natoms, 3), sums(20), nsums(100)
		
		real(dbl), dimension(:, :), allocatable  :: hessian
		
		if (init) then
			counter = 0
			sums = 0d0
			nsums = 0d0
			sys%k_ic = 0d0
			sys%k_ic_ht = 0d0
			call sys%calculate_gamma
			if (.not. allocated(sys%V_vq_j)) call sys%build_V(grads, hessian)
		end if
		
		do nx=1,stopix
			call sys%compute_zn(occs(:, nx), res)
			tmp = dot_product(res, sys%V_vq_j)
			tmp = tmp**2
			sys%k_ic = sys%k_ic + tmp
			tmp = 4d0 * tmp / sys%gamma
			ix = int(min(max(1d0, log10(tmp)+15d0), 20d0))
			counter(ix) = counter(ix) + 1
			sums(ix) = sums(ix) + tmp
			ix = int(sum(occs(:, nx)))
			nsums(ix) = nsums(ix) + tmp
			if (sys%debug_level .gt. 3) then
				if (tmp > 1.0d0) then
					write(*, *) occs(:, nx), tmp
				end if
			end if
		end do
		
		if (sys%debug_level .gt. 2) then
			write(*, *) '\nLog10', 'Counter', 'Sum'
			do nx=1,20
				write(*, *) nx-15, counter(nx), sums(nx)
			end do
		end if
		
		if (sys%debug_level .gt. 1) then
			write(*, *) '\nNSUMS'
			do nx=1,100
				if (nsums(nx) .gt. 1e-12) then
					write(*, *) nx, nsums(nx)
				end if
			end do
		end if
		
		if (sys%hessfile .ne. 'none') then
			! Do Hertzberg-Teller calculation
			do nx=1,stopix
				call sys%compute_zn(occs(:, nx), z)
				res = matmul(sys%V_jj, z)
				do ix=1,sys%nlevels
					z(ix) = 0.5 * (occs(ix, nx) + sys%hrfactors(ix))**2
					z(ix) = z(ix) / (sys%energies(ix)*sys%hrfactors(ix))
					z(ix) = sqrt(z(ix))
				end do
				tmp = dot_product(res, z)
				tmp = tmp**2
				sys%k_ic_ht = sys%k_ic_ht + tmp
			end do
		end if
		
	end subroutine sysdata_calc_kic
	
	subroutine sysdata_calc_lagrange(sys, occs, klag)
		class(sysdata), intent(inout)					:: sys
		real(dbl), dimension(sys%nlevels), intent(in)	:: occs
		real(dbl), intent(out)							:: klag
		
		
		real(dbl) 	:: tmp, kfcn, z(sys%nlevels), res(sys%nlevels), yk, nfac
		integer		:: jx
		
		kfcn = 1d0
		do jx=1, sys%nlevels
			yk = sys%hrfactors(jx)
			nfac = factorial(nint(occs(jx)))
			tmp = (yk**occs(jx)) * exp(-yk) / nfac
			kfcn = kfcn * sqrt(tmp)
		end do

		
		do jx=1, sys%nlevels
			res(jx) = 0.5d0 * sys%energies(jx) / sys%hrfactors(jx)
			tmp = occs(jx) - sys%hrfactors(jx)
			res(jx) = sqrt(res(jx) * tmp * tmp)
			res(jx) = res(jx) * kfcn
		end do
		
		tmp = dot_product(res, sys%V_vq_j)
		tmp = tmp**2
		tmp = 4d0 * tmp / sys%gamma
		klag = tmp * 1.9652703588D+3	
		
		write(*, '(1x,a,1x,e14.8,1x,a)') 'Lagrange k_ic_fc=', klag, 's-1'
		
		if (sys%hessfile .ne. 'none') then
			z = matmul(sys%V_jj, res)
			do jx=1, sys%nlevels
				res(jx) = 0.5 * (occs(jx) + sys%hrfactors(jx))**2
				res(jx) = res(jx) / (sys%energies(jx)*sys%hrfactors(jx))
				res(jx) = sqrt(res(jx))
			end do
			tmp = dot_product(res, z)
			tmp = tmp**2
			tmp = 4d0 * tmp / sys%gamma
			tmp = tmp * 1.9652703588D+3
			
			write(*, '(1x,a,1x,e14.8,1x,a)') 'Lagrange k_ic_ht=', tmp, 's-1'
		end if	
		write(*, *)
					
	end subroutine sysdata_calc_lagrange
	
	subroutine sysdata_get_next_bound(sys, occ, ix, bnd)
		class(sysdata), intent(in)					:: sys
		integer(smallint), dimension(sys%nlevels), intent(in)	:: occ
		integer, intent(in)							:: ix
		integer, intent(out)						:: bnd
		
		integer	:: jx
		bnd = sys%nthresh + 2
		do jx=1,ix
			bnd = bnd - sys%bounds(jx, occ(jx)+1)
		end do
		bnd = max(1, bnd)

	end subroutine sysdata_get_next_bound
		
	subroutine sysdata_ncombinations(sys, maxnocc)
		class(sysdata), intent(in)		:: sys
		integer(bigint), intent(out)	:: maxnocc
		
		integer(smallint) :: occ(sys%nlevels)
		maxnocc = 0
		occ = 0
		occ(1) = sys%cutoffs(1, sys%nthresh+1)
		call next_n(sys, occ, 1, maxnocc)
	end subroutine
	
	recursive subroutine next_n(sys, occ, ix, maxnocc)
		type(sysdata), intent(in)						:: sys
		integer(bigint), intent(out) 					:: maxnocc 
		integer, intent(in)								:: ix
		integer(smallint), dimension(sys%nlevels), intent(inout)	:: occ
		
		integer :: maxix, jx, bnd
		integer(bigint) :: tmpmax
		maxix = occ(ix)
		maxnocc = 1
		if (ix .eq. sys%nlevels) then
			maxnocc = maxix + 1
		else
			do jx=0,maxix
				occ(ix) = jx
				call sys%get_next_bound(occ, ix, bnd)
				occ(ix+1) = sys%cutoffs(ix+1, bnd)
				call next_n(sys, occ, ix+1, tmpmax)
				maxnocc = maxnocc + tmpmax
			end do
		end if
	end subroutine next_n	
end module