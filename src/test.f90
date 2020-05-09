program main
	use constants
	use random
	use ioutil
	
	real(dbl), dimension(:), allocatable  :: testlist, outlist
	integer, dimension(:), allocatable	  :: ixlist, indices
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
	mm%current_record = 2
	call mm%write_to_bin(mm%current_block, indices)
	call cpu_time(finish)
	write(*, *) 'Write time: ', finish-start
	deallocate(mm%current_block)
	
	call cpu_time(start)
	call mm%read_from_bin(n, 2, mm%current_block)
	call cpu_time(finish)
	write(*, *) 'Read time: ', finish-start
	write(*, *) mm%current_block(:, 1)
	
	allocate(indices(listsize))
	indices(1:listsize) = (/ (j, j=1,listsize) /)
	call cpu_time(start)
	call mm%sort_and_screen(testlist, testlist, 6d0, max_ix, indices, energy)
	call cpu_time(finish)
	write(*, *) 'Sort time: ', finish-start
	write(*, *) 'No. of non-negligible records: ', max_ix
	
	mm%current_record = 1
	mm%screening_type = 'energy'
	mm%threshold = 6d0
	write(*, *) mm%current_block(:, 1004)
	write(*, *) mm%current_block(:, indices(1004))
	call cpu_time(start)
	call mm%block_swap(testlist, testlist)
	call cpu_time(finish)
	write(*, *) mm%current_block(:, 1004)
	write(*, *) 'Current record: ', mm%current_record
	write(*, *) 'Block swap time: ', finish-start
	
	call cpu_time(start)
	call mm%read_from_bin(n, 1, mm%current_block)
	call cpu_time(finish)
	write(*, *) 'Read time: ', finish-start
	write(*, *) mm%current_block(:, 1004)
	deallocate(indices)
	
	call sort_list(n, testlist, outlist, ixlist)
	write(*, '(A10, 20F8.4)') 'INPUT:', testlist
	write(*, '(A10, 20F8.4)') 'BUBBLESORT:', outlist
	write(*, '(A10, 20I8)') 'INDICES:', ixlist
	
	allocate(indices(n))
	indices(1:n) = (/ (j, j=1,n) /)
	call quicksort(n, testlist, indices)
	write(*, '(A10, 20F8.4)') 'QUICKSORT:', testlist
	write(*, '(A10, 20I8)') 'INDICES:', indices
	
end program main