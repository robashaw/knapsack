module ioutil
	use constants
	implicit none
	
	type memorymanager
		real(dbl) 	:: max_memory
		integer		:: nrecords, current_record
	contains
		procedure	:: initialise => mm_init
		procedure	:: record_name => mm_record_name
		procedure	:: write_to_bin => mm_write_to_bin
		procedure	:: read_from_bin => mm_read_from_bin
	end type memorymanager
	
contains
	
	subroutine sort_list(n, inlist, outlist, indices)
		! implements a simple bubble sort on inlist (length n)
		! yielding the sorted outlist, and the order
		! of the indices going from out to in
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
	
end module ioutil
