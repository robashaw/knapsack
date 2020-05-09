module ioutil
	use constants
	use random
	implicit none
	
	type memorymanager
		real(dbl) 							  :: max_memory, threshold
		integer								  :: nrecords, current_record, nlevels, chunk_size
		integer, dimension(:, :), allocatable :: current_block
		character(len=10)					  :: screening_type	
		logical								  :: keep_top
	contains
		procedure	:: initialise => mm_init
		procedure	:: record_name => mm_record_name
		procedure	:: block_swap => mm_block_swap
		procedure	:: sort_and_screen => mm_sort_and_screen
		procedure	:: write_to_bin => mm_write_to_bin
		procedure	:: read_from_bin => mm_read_from_bin
	end type memorymanager
	
contains
	
	subroutine sort_list(n, inlist, outlist, indices)
		! implements a simple bubble sort on inlist (length n)
		! yielding the sorted outlist, and the order
		! of the indices going from out to in
		! sorts in ASCENDING order
		integer, intent(in)						:: n
		real(dbl), dimension(n), intent(in)		:: inlist 
		real(dbl), dimension(n), intent(out)	:: outlist
		integer, dimension(n), intent(out)		:: indices
		
		integer 	:: i, j, tmpix
		real(dbl) 	:: tmpval
		indices(1:n) = (/ (i, i=1,n) /)
		outlist(1:n) = inlist(1:n)
		
		do i = 1, n-1
			do j = i+1, n
				if (outlist(i) .gt. outlist(j)) then
					tmpval    = outlist(i)
					outlist(i) = outlist(j)
					outlist(j) = tmpval
					
					tmpix      = indices(i)
					indices(i) = indices(j)
					indices(j) = tmpix 
				end if
			end do
		end do
	end subroutine sort_list
	
	recursive subroutine quicksort(n, inlist, indices) 
		! NOTE: sorts in ASCENDING order of inlist 
		integer, intent(in)						:: n
		real(dbl), dimension(n), intent(inout)	:: inlist
		integer, dimension(n), intent(out)		:: indices
		
		integer 				:: i, k, p, m, j  ! indices
		real(dbl)				:: temp_real	  ! temp value
		integer					:: temp_ix		  ! temp index
		real(dbl), dimension(n) :: b, c 		  ! temp arrays
		real(dbl) 				:: random_real    ! random number
		integer					:: random_ix 	  ! random index
		
		if (n .eq. 1) then
			return
		end if
		
		! Choose pivot and swap
		random_ix = random_int(1, n)
		call swap(random_ix, n)
		
		! 3-way partition
		i = 1
		k = 1
		p = n
		main: do
			if (inlist(i) .lt. inlist(n)) then
				call swap(i, k)
				i = i + 1
				k = k + 1
			else if (inlist(i) - inlist(n) .lt. 1E-12) then
				call swap(i, p-1)
				p = p - 1
			else
				i = i + 1
			end if 
			
			! exit condition
			if (i .ge. p) then
				exit main
			end if
		end do main
		
		! move pivots to centre
		m = min(p-k, n-p+1)
		if (m .gt. 0) then
			do j=1,m
				call swap(j+k-1, n-m+j)
			end do
		end if
		
		! recursive sorts
		if (k .gt. 1) then
			call quicksort(k-1, inlist(1:k-1), indices(1:k-1))
		end if 
		
		if (m .gt. 0) then
			j = n-p+k+1
			call quicksort(p-k, inlist(j:n), indices(j:n))
		end if
		
	contains
		subroutine swap(ix, jx)
			integer, intent(in) :: ix, jx
			temp_real     = inlist(ix)
			temp_ix       = indices(ix)
			inlist(ix)    = inlist(jx)
			inlist(jx)    = temp_real
			indices(ix)   = indices(jx)
			indices(jx)   = temp_ix
		end subroutine swap
	end subroutine quicksort
	
	subroutine reorder_list(n, list, indices)
		integer, intent(in)						:: n
		real(dbl), dimension(n), intent(inout)	:: list 
		integer, dimension(n), intent(in)		:: indices
	
		integer :: ix
		real(dbl), dimension(n)	:: listcopy
		listcopy(:) = list(:)
		
		do ix=1,n
			list(ix) = listcopy(indices(ix))
		end do
	end subroutine reorder_list
	
	subroutine mm_init(mm, n, max_noccs, max_mem)
		class(memorymanager), intent(inout)		:: mm
		integer, intent(in) 					:: n
		integer(bigint), intent(in)				:: max_noccs
		real(dbl), intent(in)					:: max_mem
		
		real(dbl) :: mem_estimate
		
		mm%nlevels = n
		mm%screening_type = 'fc' ! default to Franck-Condon screening
		mm%threshold = 1D-12 ! threshold for screening 
		mm%keep_top = .false.
		
		! maximum memory availale in GB
		if (max_mem .gt. 0d0) then
			mm%max_memory = max_mem
		else
			mm%max_memory= 1.0d0 ! default to 1 GB
		end if
		
		! size in bites of a single occupation
		mem_estimate = sizeof(n) *  n 
		! maximum size if all occs enumerated, convert to GB
		mem_estimate = (mem_estimate * real(max_noccs)) / (1024d0**3)
		mm%nrecords = ceiling(mem_estimate / mm%max_memory)	
		mm%chunk_size = ceiling(real(max_noccs) / real(mm%nrecords))
		mm%current_record = 1	
		
		write(*, '(A33,F15.3,A3)') 'Initalised memory manager with:', mm%max_memory, 'GB'
		write(*, '(A33,F15.3,A3)') 'Maximum memory needed (estimate):', mem_estimate, 'GB'
		write(*, '(A12,I10,A13,I10,A12)') 'Anticipating', mm%nrecords, 'records with', mm%chunk_size, 'occs each'
	end subroutine mm_init
	
	subroutine mm_record_name(mm, filename, record)
		class(memorymanager), intent(inout)	:: mm
		character(len=*), intent(out)		:: filename
		integer, intent(in)					:: record
		
		write(filename, *) record
		filename = 'occs.' // trim(adjustl(filename))
	end subroutine mm_record_name
	
	subroutine mm_sort_and_screen(mm, levels, hrfactors, thresh, max_ix, indices, func)
		class(memorymanager), intent(inout)				:: mm
		real(dbl), intent(in)							:: thresh
		integer, intent(out)							:: max_ix
		integer, dimension(mm%chunk_size), intent(out)	:: indices
		real(dbl), dimension(mm%nlevels), intent(in)	:: levels, hrfactors
		
		interface 
			function func(n, levels, hrfactors, occs)
				import :: dbl
				integer, intent(in)					:: n
				real(dbl), dimension(n), intent(in) :: levels, hrfactors
				integer, dimension(n), intent(in) 	:: occs
				real(dbl)							:: func
			end function
		end interface
		
		! Calculate values using func
		real(dbl), dimension(mm%chunk_size) :: values
		integer :: ix
		do ix=1,mm%chunk_size
			values(ix) = -func(mm%nlevels, levels, hrfactors, mm%current_block(:, ix))
		end do
		
		! sort in descending order
		indices(1:mm%chunk_size) = (/ (ix, ix=1,mm%chunk_size) /)
		call quicksort(mm%chunk_size, values, indices)
		
		! find the cutoff index
		max_ix = 1
		main: do
			if (-values(max_ix) .gt. thresh) then
				max_ix = max_ix + 1
			else 
				exit main
			end if
		end do main
	end subroutine mm_sort_and_screen
	
	subroutine mm_block_swap(mm, levels, hrfactors)
		class(memorymanager), intent(inout)				:: mm
		real(dbl), dimension(mm%nlevels), intent(in)	:: levels, hrfactors
		
		integer	:: max_ix, ix
		integer, dimension(mm%chunk_size) :: indices
		integer, dimension(mm%nlevels)	  :: tmp_occ
		
		! check if there is currently a block
		if (allocated(mm%current_block)) then
			! sort and screen the block
			write(*, *)
			write(*, '(A15,I5,A16,1X,A10)') 'Sorting record', mm%current_record, 'screening by', mm%screening_type
			select case(mm%screening_type)
			case ('energy')
				call mm%sort_and_screen(levels, hrfactors, mm%threshold, max_ix, indices, energy)
			case default
				call mm%sort_and_screen(levels, hrfactors, mm%threshold, max_ix, indices, franck_condon)
			end select
			
			! we now only need to write the current block to file up to max_ix
			call mm%write_to_bin(mm%current_block, indices(:max_ix))
			! reset the block and increment counter
			mm%current_record = mm%current_record + 1
			if (mm%keep_top .and. (max_ix .lt. mm%chunk_size)) then
				! stochastic case need to keep the 'top' records
				mm%current_block(:, max_ix+1:) = 0
			else 
				mm%current_block = 0
			end if 
		else
			! allocate the first block
			allocate(mm%current_block(mm%nlevels, mm%chunk_size))
			mm%current_block = 0 ! initialise
		end if
	end subroutine mm_block_swap
	
	subroutine mm_write_to_bin(mm, occs, indices)
		class(memorymanager), intent(inout)		:: mm
		integer, dimension(:, :), intent(in) 	:: occs
		integer, dimension(:), intent(in)		:: indices
		
		character(len=100) :: filename
		integer	:: ios
		
		call mm%record_name(filename, mm%current_record)
		open(unit=occs_unit, file=filename, form='unformatted', access='stream')
		
		if (size(indices) .eq. 0) then
			write(occs_unit, iostat=ios) occs
		else
			write(occs_unit, iostat=ios) occs(:, indices(:))
		end if
		close(occs_unit)
		
		if (ios .ne. 0) then
			write(*, *) 'Error writing record ', trim(filename), ' with ios=', ios
		else
			write(*, *) 'Successfully wrote record ', filename
		end if
	end subroutine mm_write_to_bin
	
	subroutine mm_read_from_bin(mm, n, record, occs)
		class(memorymanager), intent(inout)					:: mm
		integer, intent(in)									:: n, record
		integer, dimension(:, :), allocatable, intent(out) 	:: occs
		
		integer	:: ios, nrows, i, filesize, firstrow(n)
		character(len=100) :: filename
		
		
		call mm%record_name(filename, record)
		open(unit=occs_unit, file=filename, form='unformatted', access='stream', iostat=ios)
		if (ios .eq. 0) then
			inquire(occs_unit, size=filesize)
			read(occs_unit) firstrow
			nrows = filesize / sizeof(firstrow)
			allocate(occs(n, nrows))
			
			occs(:, 1) = firstrow(:)
			do i=2,nrows
				read(occs_unit, iostat=ios) occs(:, i)
				if (ios .ne. 0) then
					if (i .lt. nrows) then
						write(*, *) 'Only found ', i-1, ' entries, expected ', nrows
					end if
					exit
				end if
			end do
			close(occs_unit)
			
			write(*, *) 'Successfully read record ', trim(filename), ' with dimensions ', shape(occs)
		else
			write(*, *) 'Error reading record ', filename
		end if
	end subroutine mm_read_from_bin
	
	integer(bigint) function ncombinations(n, max_occs, min_occs) result(max_nocc)
		integer, intent(in)					:: n
		integer, dimension(n), intent(in)	:: max_occs, min_occs
		integer								:: ix
		
		max_nocc = 1
		do ix=1,n
			max_nocc = max_nocc * (max_occs(ix) - min_occs(ix) + 1)
		end do
	end function ncombinations
	
	function energy(n, levels, hrfactors, occs) 
		integer, intent(in)					:: n
		real(dbl), dimension(n), intent(in) :: levels, hrfactors
		integer, dimension(n), intent(in) 	:: occs
		real(dbl)							:: energy
		
		energy = dot_product(occs, levels)
	end function energy
	
	function franck_condon(n, levels, hrfactors, occs)
		integer, intent(in)					:: n
		real(dbl), dimension(n), intent(in) :: levels, hrfactors
		integer, dimension(n), intent(in) 	:: occs
		real(dbl)							:: franck_condon
		
		integer 	:: ix
		real(dbl) 	:: tmp
		franck_condon = 1d0
		do ix=1,n
			tmp = exp(hrfactors(ix)) * (hrfactors(ix) ** occs(ix)) / factorial(occs(ix))
			franck_condon =  franck_condon * tmp
		end do
		franck_condon = sqrt(franck_condon)
	end function franck_condon
end module ioutil
