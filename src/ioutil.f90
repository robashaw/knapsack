module ioutil
	use constants
	use random
	implicit none
	
	type memorymanager
		real(dbl) 							  :: max_memory
		integer								  :: nrecords, current_record, nlevels, chunk_size
		integer, dimension(:, :), allocatable :: current_block
	contains
		procedure	:: initialise => mm_init
		procedure	:: record_name => mm_record_name
		!procedure	:: block_swap => mm_block_swap
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
			integer				:: random_ix 	  ! random index
		
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
			call quicksort(p-k, inlist(n-p+k+1:n), indices(n-p+k+1:n))
		end if
		
	contains
		subroutine swap(ix, jx)
			integer, intent(in) :: ix, jx
			temp_real = inlist(ix)
			temp_ix   = indices(ix)
			inlist(ix) = inlist(jx)
			inlist(jx) = temp_real
			indices(ix) = indices(jx)
			indices(jx) = temp_ix
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
		write(*, *) mm%chunk_size
		mm%current_record = 1	
		
		write(*, '(A33,F15.3,A3)') 'Initalised memory manager with:', mm%max_memory, 'GB'
		write(*, '(A33,F15.3,A3)') 'Maximum memory needed (estimate):', mem_estimate, 'GB'
		write(*, '(A12,I10,A8)') 'Anticipating', mm%nrecords, 'records'
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
	
	subroutine mm_write_to_bin(mm, n, occs)
		class(memorymanager), intent(inout)		:: mm
		integer, intent(in)						:: n
		integer, dimension(:, :), intent(in) 	:: occs
		
		character(len=100) :: filename
		integer	:: ios
		
		call mm%record_name(filename, mm%current_record)
		mm%current_record = mm%current_record + 1
		open(unit=occs_unit, file=filename, form='unformatted', access='stream')
		write(occs_unit, iostat=ios) occs
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
	
end module ioutil
