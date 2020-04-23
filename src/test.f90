program main
	use constants
	use ioutil
	
	real(dbl), dimension(:), allocatable :: testlist, outlist
	integer, dimension(:), allocatable	:: ixlist
	integer   :: i, n, listsize
	integer, dimension(:, :), allocatable :: occs
	type(memorymanager)	:: mm
	integer(bigint)	:: maxnoccs
	
	maxnoccs = 1D+9
	
	n = 20
	allocate(testlist(n), outlist(n), ixlist(n))
	
	call mm%initialise(n, maxnoccs, 0.5d0)
	do i=1, n
		testlist(i) = rand()
	end do
	
	listsize = ceiling(real(maxnoccs) / real(mm%nrecords))
	allocate(occs(n, listsize))
	occs = 0
	call mm%write_to_bin(n, occs)
	deallocate(occs)
	
	call mm%read_from_bin(n, 1, occs)
	
	call sort_list(n, testlist, outlist, ixlist)
	write(*, '(A10, 20F8.4)') 'INPUT:', testlist
	write(*, '(A10, 20F8.4)') 'SORTED:', outlist
	write(*, '(A10, 20I8)') 'INDICES:', ixlist
	
end program main