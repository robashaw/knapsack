module ioutil
	use constants
	use random
	implicit none
	
	type memorymanager
		real(dbl) 							  :: max_memory, threshold
		integer								  :: nrecords, current_record, nlevels, dumplevel=0
		integer(bigint)					      :: chunk_size
		integer(smallint), dimension(:, :), allocatable :: current_block
		character(len=10)					  :: screening_type	
		logical								  :: keep_top
	contains
		procedure	:: initialise => mm_init
		procedure	:: record_name => mm_record_name
		procedure	:: block_swap => mm_block_swap
		procedure	:: sort_and_screen => mm_sort_and_screen
		procedure	:: write_to_bin => mm_write_to_bin
		procedure	:: read_from_bin => mm_read_from_bin
		procedure	:: merge_files => mm_merge_files
		procedure	:: merge_all => mm_merge_all
		procedure	:: clean_up => mm_clean_up
	end type memorymanager
	
contains
	
	subroutine sort_list(n, inlist, outlist, indices)
		! implements a simple bubble sort on inlist (length n)
		! yielding the sorted outlist, and the order
		! of the indices going from out to in
		! sorts in ASCENDING order
		integer(bigint), intent(in)				:: n
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
		integer(bigint), intent(in)				:: n
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
			integer(bigint), intent(in) :: ix, jx
			temp_real     = inlist(ix)
			temp_ix       = indices(ix)
			inlist(ix)    = inlist(jx)
			inlist(jx)    = temp_real
			indices(ix)   = indices(jx)
			indices(jx)   = temp_ix
		end subroutine swap
	end subroutine quicksort
	
	subroutine reorder_list(n, list, indices)
		integer(bigint), intent(in)				:: n
		real(dbl), dimension(n), intent(inout)	:: list 
		integer, dimension(n), intent(in)		:: indices
	
		integer :: ix
		real(dbl), dimension(n)	:: listcopy
		listcopy(:) = list(:)
		
		do ix=1,n
			list(ix) = listcopy(indices(ix))
		end do
	end subroutine reorder_list
	
	subroutine mm_init(mm, n, max_noccs, max_mem, nthresh)
		class(memorymanager), intent(inout)		:: mm
		integer, intent(in) 					:: n, nthresh
		integer(bigint), intent(in)				:: max_noccs
		real(dbl), intent(in)					:: max_mem
		
		real(dbl) :: mem_estimate
		integer(smallint)	:: sizer = 0
		
		mm%nlevels = n
		mm%screening_type = 'fc' ! default to Franck-Condon screening
		mm%threshold = 10.0**(-nthresh) ! threshold for screening 
		mm%keep_top = .false.
		
		! maximum memory availale in GB
		if (max_mem .gt. 0d0) then
			mm%max_memory = max_mem
		else
			mm%max_memory= 1.0d0 ! default to 1 GB
		end if
		
		! size in bites of a single occupation
		mem_estimate = sizeof(sizer) *  n 
		! maximum size if all occs enumerated, convert to GB
		mem_estimate = (mem_estimate * real(max_noccs)) / (1024d0**3)
		mm%nrecords = ceiling(mem_estimate / mm%max_memory)	
		mm%chunk_size = ceiling(real(max_noccs) / real(mm%nrecords))
		mm%current_record = 1	
		
		write(*, '(A33,F15.3,A3)') 'Initalised memory manager with:', mm%max_memory, 'GB'
		write(*, '(A33,F15.3,A3)') 'Maximum memory needed (estimate):', mem_estimate, 'GB'
		write(*, '(A12,I10,A13,I10,A12)') 'Anticipating', mm%nrecords, 'records with', mm%chunk_size, 'occs each'
	end subroutine mm_init
	
	subroutine mm_record_name(mm, filename, record, prefix)
		class(memorymanager), intent(inout)	:: mm
		character(len=*), intent(in)		:: prefix
		character(len=*), intent(out)		:: filename
		integer, intent(in)					:: record
		
		write(filename, *) record
		filename = trim(adjustl(prefix)) // '.' // trim(adjustl(filename))
	end subroutine mm_record_name
	
	subroutine mm_clean_up(mm, prefix, minix, maxix)
		class(memorymanager), intent(inout)	:: mm
		character(len=*), intent(in)		:: prefix
		integer, intent(in)					:: minix, maxix
		
		integer				:: jx, ios
		character(len=100)	:: tmpfile
		
		do jx=minix,maxix
			call mm%record_name(tmpfile, jx, prefix)
			open(unit=occs_unit, iostat=ios, file=tmpfile, status='old')
			if (ios == 0) then
				close(occs_unit, status='delete')
				write(*, '(1x,a,1x,a)') 'Deleted', tmpfile
			end if
		end do
	end subroutine mm_clean_up
	
	subroutine mm_sort_and_screen(mm, levels, hrfactors, thresh, max_ix_in, max_ix_out, indices, func)
		class(memorymanager), intent(inout)				:: mm
		real(dbl), intent(in)							:: thresh
		integer, intent(in)								:: max_ix_in
		integer, intent(out)							:: max_ix_out
		integer, dimension(max_ix_in), intent(out)		:: indices
		real(dbl), dimension(mm%nlevels), intent(in)	:: levels, hrfactors
		
		interface 
			function func(n, levels, hrfactors, occs)
				import :: dbl, smallint
				integer, intent(in)					:: n
				real(dbl), dimension(n), intent(in) :: levels, hrfactors
				integer(smallint), dimension(n), intent(in) 	:: occs
				real(dbl)							:: func
			end function
		end interface
		
		! Calculate values using func
		real(dbl), dimension(max_ix_in) :: values
		integer(smallint), dimension(mm%nlevels) :: cocc, lastocc
		integer :: ix
		integer, dimension(max_ix_in)	:: unique_indices
		do ix=1,max_ix_in
			values(ix) = -func(mm%nlevels, levels, hrfactors, mm%current_block(:, ix))
		end do
		
		! sort in descending order
		indices(1:max_ix_in) = (/ (ix, ix=1,max_ix_in) /)
		call quicksort(max_ix_in, values, indices)
		
		! find the cutoff index
		max_ix_out = 1
		ix = 1
		lastocc = 0
		unique_indices = 0
		main: do while (ix .le. max_ix_in)
			cocc(:) = mm%current_block(:,indices(ix))
			if (-values(ix) .gt. thresh) then
				if ((sum(abs(cocc - lastocc)) .ne. 0) .and. (sum(cocc) .ne. 0)) then
					unique_indices(max_ix_out) = indices(ix)
					max_ix_out = max_ix_out + 1
					lastocc(:) = cocc(:)
				end if
				ix = ix + 1
			else 
				exit main
			end if
		end do main
		max_ix_out = max_ix_out - 1
		indices(1:max_ix_out) = unique_indices(1:max_ix_out)
	end subroutine mm_sort_and_screen
	
	subroutine mm_block_swap(mm, levels, hrfactors, max_ix_in, max_ix_out)
		class(memorymanager), intent(inout)				:: mm
		real(dbl), dimension(mm%nlevels), intent(in)	:: levels, hrfactors
		integer, intent(in)								:: max_ix_in
		integer, intent(out)							:: max_ix_out
		
		integer	:: ix
		integer, dimension(max_ix_in) 	:: indices
		integer, dimension(mm%nlevels)	:: tmp_occ
		integer(smallint), dimension(mm%nlevels, topkeep) :: thetop
		
		! check if there is currently a block
		if (allocated(mm%current_block)) then
			! sort and screen the block
			write(*, *)
			write(*, '(A15,I5,A16,1X,A10)') 'Sorting record', mm%current_record, 'screening by', mm%screening_type
			select case(mm%screening_type)
			case ('energy')
				call mm%sort_and_screen(levels, hrfactors, mm%threshold, max_ix_in, max_ix_out, indices, energy)
			case default
				call mm%sort_and_screen(levels, hrfactors, mm%threshold, max_ix_in, max_ix_out, indices, franck_condon)
			end select
			
			! we now only need to write the current block to file up to max_ix
			if (max_ix_out .eq. 0) then
				write(*, *) 'Nothing to write'
			else 
				call mm%write_to_bin(mm%current_block, indices(:max_ix_out))
				! reset the block and increment counter
				mm%current_record = mm%current_record + 1
				if (mm%keep_top .and. (max_ix_out .lt. mm%chunk_size)) then
					! stochastic case need to keep the 'top' records
					ix = min(max_ix_out, topkeep)
					thetop(:, :ix) = mm%current_block(:, indices(:ix))
					mm%current_block = 0
					mm%current_block(:, :ix) = thetop(:, :ix)
					max_ix_out = ix
				else 
					mm%current_block = 0
				end if 
			end if
		else
			! allocate the first block
			allocate(mm%current_block(mm%nlevels, mm%chunk_size))
			mm%current_block = 0 ! initialise
		end if
	end subroutine mm_block_swap
	
	subroutine mm_write_to_bin(mm, occs, indices)
		class(memorymanager), intent(inout)		:: mm
		integer(smallint), dimension(:, :), intent(in) 	:: occs
		integer, dimension(:), intent(in)		:: indices
		
		character(len=100) :: filename
		integer	:: ios
		
		call mm%record_name(filename, mm%current_record, 'occs')
		open(unit=output_unit, file=filename, form='unformatted', access='stream')
		
		if (size(indices) .eq. 0) then
			write(output_unit, iostat=ios) occs
		else
			write(output_unit, iostat=ios) occs(:, indices(:))
		end if
		close(output_unit)
		
		if (ios .ne. 0) then
			write(*, *) 'Error writing record ', trim(filename), ' with ios=', ios
		else
			write(*, *) 'Successfully wrote record ', filename
		end if
		write(*, *) "\n"
	end subroutine mm_write_to_bin
	
	subroutine mm_read_from_bin(mm, n, record, occs, levels, hrfactors, prefix)
		class(memorymanager), intent(inout)								:: mm
		integer, intent(in)												:: n, record
		integer(smallint), dimension(:, :), allocatable, intent(out) 	:: occs
		real(dbl), dimension(mm%nlevels), intent(in)					:: levels, hrfactors
		character(len=*), intent(in)									:: prefix
		
		integer	:: ios, nrows, i, filesize
		integer(smallint) :: firstrow(n)
		real(dbl) :: en, fc
		character(len=100) :: filename
		
		
		call mm%record_name(filename, record, prefix)
		open(unit=occs_unit, file=filename, form='unformatted', access='stream', iostat=ios)
		if (ios .eq. 0) then
			inquire(occs_unit, size=filesize)
			read(occs_unit) firstrow
			nrows = filesize / sizeof(firstrow)
			if (nrows .gt. mm%chunk_size) write(*, *) 'Warning: file needs more memory than is available'
			if (allocated(occs)) then
				deallocate(occs)
			end if
			allocate(occs(n, nrows))
			mm%chunk_size = nrows
			write(*, *) 'Expecting ', nrows, ' rows'
			
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
			write(*, *) "\n"
			
			if (mm%dumplevel .eq. 1) then
				write(*, *) 'OCCUPATIONS'
				do i=1,nrows
					write(*, *) occs(:, i)
				end do
			else if (mm%dumplevel .eq. 2) then
				write(*, '(1x, a12, 1x, a12)') 'ENERGY (EV)', 'FC FACTOR'
				do i=1,nrows
					en = energy(mm%nlevels, levels, hrfactors, occs(:, i))
					fc = franck_condon(mm%nlevels, levels, hrfactors, occs(:, i))
					write(*, '(1x,d12.4,2x,d12.4)') en, fc
				end do
			end if
		else
			write(*, *) 'Error reading record ', filename
		end if
	end subroutine mm_read_from_bin
	
	subroutine mm_merge_all(mm, levels, hrfactors, outfile, outrecord)
		class(memorymanager), intent(inout)				:: mm
		real(dbl), dimension(mm%nlevels), intent(in)	:: levels, hrfactors
		character(len=*), intent(out)					:: outfile
		integer, intent(out)							:: outrecord
		
		integer ::	occ_record, merge_record
		character(len=100) :: file1, file2, file3
		if (mm%current_record .eq. 2) then
			outfile = 'occs'
			outrecord = 1
		else
			if (mm%screening_type .eq. 'energy') then
				call mm%merge_files('occs.1', 'occs.2', 'merged.1', levels, hrfactors, energy)
			else
				call mm%merge_files('occs.1', 'occs.2', 'merged.1', levels, hrfactors, franck_condon)
			end if
			
			merge_record = 1
			do occ_record = 3,(mm%current_record-1)
				call mm%record_name(file1, record=merge_record, prefix='merged')
				call mm%record_name(file2, record=occ_record, prefix='occs')
				call mm%record_name(file3, record=merge_record+1, prefix='merged')
				if (mm%screening_type .eq. 'energy') then
					call mm%merge_files(file1, file2, file3, levels, hrfactors, energy)
				else
					call mm%merge_files(file1, file2, file3, levels, hrfactors, franck_condon)
				end if
				merge_record = merge_record + 1
			end do 
			
			outfile = 'merged'
			outrecord = merge_record
		end if 
	end subroutine mm_merge_all
	
	subroutine mm_merge_files(mm, file1, file2, outfile, levels, hrfactors, func)
		class(memorymanager), intent(inout)	:: mm
		character(len=*), intent(in)		:: file1, file2, outfile
		real(dbl), dimension(mm%nlevels), intent(in)	:: levels, hrfactors
		
		interface 
			function func(n, levels, hrfactors, occs)
				import :: dbl, smallint
				integer, intent(in)					:: n
				real(dbl), dimension(n), intent(in) :: levels, hrfactors
				integer(smallint), dimension(n), intent(in) 	:: occs
				real(dbl)							:: func
			end function
		end interface
		
		integer 			:: maxrecords, ios, nrows1, nrows2, i, filesize, delta
		integer				:: size1, size2, ptr1, ptr2, start1, start2, ctr1, ctr2
		integer(smallint) 	:: firstrow(mm%nlevels), sizer=0
		real(dbl)			:: val1, val2
		integer(smallint), dimension(:, :), allocatable	:: f1occs, f2occs
		
		write(*, '(1x,a,1x,a,1x,a,1x,a)') 'Merging', trim(adjustl(file1)), 'and', trim(adjustl(file2))
		
		! calculate how many records can open from each file
		maxrecords = floor(0.5 * mm%max_memory * (1024**3) / (sizeof(sizer) * mm%nlevels))
		! open the two files
		open(unit=occs_unit, file=file1, form='unformatted', access='stream', iostat=ios)
		open(unit=occs_unit_2, file=file2, form='unformatted', access='stream', iostat=ios)
		open(unit=output_unit, file=outfile, form='unformatted', access='stream')
		if (ios .eq. 0) then
			ptr1 = 1
			ptr2 = 1
			start1 = 2
			start2 = 2
			ctr1 = 0
			ctr2 = 0
			! file 1 nrecords
			inquire(occs_unit, size=filesize)
			read(occs_unit) firstrow
			nrows1 = filesize / sizeof(firstrow)
			size1 = min(nrows1, maxrecords)
			allocate(f1occs(mm%nlevels, size1))
			f1occs(:, 1) = firstrow(:)
			write(*, '(1x,a,1x,i10,1x,a)') 'File1 has', nrows1, 'rows'
			
			! file 2 nrecords
			inquire(occs_unit_2, size=filesize)
			read(occs_unit_2) firstrow
			nrows2 = filesize / sizeof(firstrow)
			size2 = min(nrows2, maxrecords)
			allocate(f2occs(mm%nlevels, size2))
			f2occs(:, 1) = firstrow(:)
			write(*, '(1x,a,1x,i10,1x,a)') 'File2 has', nrows2, 'rows'
			ptr1 = size1
			ptr2 = size2
			main: do 
				if (ptr1 .eq. size1) size1 = min(size1, nrows1-ctr1)
				if (ptr2 .eq. size2) size2 = min(size2, nrows2-ctr2)
				if (size1 .eq. 0) then
					 if (size2 .eq. 0) exit main
					 ! write the rest of file2 to output
					 do while (ctr2 .lt. nrows2)
					 	do i=ptr2,size2
							 write(unit=output_unit, iostat=ios) f2occs(:, i)
						end do
						
						size2 = min(size2, nrows2-ctr2)
						start2 = 1
						ptr2 = 1
						
						do i=start2,size2
							read(occs_unit_2, iostat=ios) f2occs(:, i)
						end do
						ctr2 = ctr2 + size2
					end do
				 	do i=ptr2,size2
						 write(unit=output_unit, iostat=ios) f2occs(:, i)
					end do			
					size2 = 0
				else if (size2 .eq. 0) then
				 	do while (ctr1 .lt. nrows1)
				 	   	do i=ptr1,size1
							 write(unit=output_unit, iostat=ios) f1occs(:, i)
						end do
					
						size1 = min(size1, nrows1-ctr1)
						start1 = 1
						ptr1 = 1
						
						do i=start1,size1
							read(occs_unit, iostat=ios) f1occs(:, i)
						end do
						ctr1 = ctr1 + size1
					end do
			 	   	do i=ptr1,size1
						 write(unit=output_unit, iostat=ios) f1occs(:, i)
					end do
					size1 = 0
				else
					if (ptr1 .ge. size1) then
						do i=start1,size1
							read(occs_unit, iostat=ios) f1occs(:, i)
						end do
						ctr1 = ctr1 + size1
						start1 = 1
						ptr1 = 1
					end if 
						
					if (ptr2 .ge. size2) then
						do i=start2,size2
							read(occs_unit_2, iostat=ios) f2occs(:, i)
						end do
						ctr2 = ctr2 + size2
						start2 = 1
						ptr2 = 1
					end if 
					
					val1 = func(mm%nlevels, levels, hrfactors, f1occs(:, ptr1))
					val2 = func(mm%nlevels, levels, hrfactors, f2occs(:, ptr2))
					do while ((ptr1 .lt. size1) .and. (ptr2 .lt. size2))
						delta = sum(abs(f1occs(:, ptr1) - f2occs(:, ptr2)))
						if (delta .eq. 0) then
							ptr2 = ptr2 + 1
							val2 = func(mm%nlevels, levels, hrfactors, f2occs(:, ptr2))
						else if (val1 .gt. val2) then
							write(unit=output_unit, iostat=ios) f1occs(:, ptr1)
							ptr1 = ptr1 + 1
							val1 = func(mm%nlevels, levels, hrfactors, f1occs(:, ptr1))
						else 
							write(unit=output_unit, iostat=ios) f2occs(:, ptr2)
							ptr2 = ptr2 + 1
							val2 = func(mm%nlevels, levels, hrfactors, f2occs(:, ptr2))
						end if
					end do
				end if
			end do main
			
			close(occs_unit)
			close(occs_unit_2)
			close(output_unit)
		else
			write(*, *) 'Error reading records'
		end if
		
	end subroutine mm_merge_files
	
	integer(bigint) function ncombinations(n, max_occs, min_occs) result(max_nocc)
		integer, intent(in)							:: n
		integer(smallint), dimension(n), intent(in)	:: max_occs, min_occs
		integer										:: ix
		
		max_nocc = 1
		do ix=1,n
			max_nocc = max_nocc * (max_occs(ix) - min_occs(ix) + 1)
		end do
	end function ncombinations
	
	function energy(n, levels, hrfactors, occs) 
		integer, intent(in)								:: n
		real(dbl), dimension(n), intent(in) 			:: levels, hrfactors
		integer(smallint), dimension(n), intent(in) 	:: occs
		real(dbl)										:: energy
		
		energy = dot_product(occs, levels)
	end function energy
	
	function franck_condon(n, levels, hrfactors, occs)
		integer, intent(in)								:: n
		real(dbl), dimension(n), intent(in) 			:: levels, hrfactors
		integer(smallint), dimension(n), intent(in) 	:: occs
		real(dbl)										:: franck_condon
		
		integer 	:: ix
		real(dbl) 	:: tmp
		franck_condon = 1d0
		do ix=1,n
			tmp = exp(hrfactors(ix)) * (hrfactors(ix) ** occs(ix)) / factorial(int(occs(ix)))
			franck_condon =  franck_condon * tmp
		end do
		franck_condon = sqrt(franck_condon)
	end function franck_condon
end module ioutil
