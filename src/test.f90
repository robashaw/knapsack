program main
	use constants
	use random
	use ioutil
	
	real(dbl), dimension(:), allocatable :: testlist, outlist
	integer, dimension(:), allocatable	:: ixlist, indices
	integer   :: i, j, n, listsize, max_ix
	type(memorymanager)	:: mm
	integer(bigint)	:: maxnoccs
	real(dbl)	:: start, finish
	
	maxnoccs = 1D+9
	
	n = 8
	allocate(testlist(n), outlist(n), ixlist(n))
	
	call mm%initialise(n, maxnoccs, 1.0d0)
	do i=1, n
		testlist(i) = rand()
	end do
	
	listsize = ceiling(real(maxnoccs) / real(mm%nrecords))
	allocate(mm%current_block(n, listsize))
	mm%current_block = 0
	do j=1,listsize
		do i=1,n
			mm%current_block(i, j) = random_int(0, 5)
		end do
	end do 
	write(*, *) mm%current_block(:, 1)
	call cpu_time(start)
	call mm%write_to_bin(n, mm%current_block)
	call cpu_time(finish)
	write(*, *) 'Write time: ', finish-start
	deallocate(mm%current_block)
	
	call cpu_time(start)
	call mm%read_from_bin(n, 1, mm%current_block)
	call cpu_time(finish)
	write(*, *) 'Read time: ', finish-start
	write(*, *) mm%current_block(:, 1)
	
	allocate(indices(listsize))
	indices(1:listsize) = (/ (j, j=1,listsize) /)
	call cpu_time(start)
	call mm%sort_and_screen(testlist, testlist, 10d0, max_ix, indices, energy)
	call cpu_time(finish)
	write(*, *) 'Sort time: ', finish-start
	write(*, *) max_ix
	deallocate(indices)
	
	call sort_list(n, testlist, outlist, ixlist)
	write(*, '(A10, 20F8.4)') 'INPUT:', testlist
	write(*, '(A10, 20F8.4)') 'SORTED:', outlist
	write(*, '(A10, 20I8)') 'INDICES:', ixlist
	
	allocate(indices(n))
	indices(1:n) = (/ (j, j=1,n) /)
	call quicksort(n, testlist, indices)
	write(*, '(A10, 20F8.4)') 'INPUT:', testlist
	write(*, '(A10, 20F8.4)') 'SORTED:', outlist
	write(*, '(A10, 20I8)') 'INDICES:', indices
	
end program main