module hashtbl
	use constants, only : dbl
	implicit none
	integer, parameter :: tbl_size = 50
	
	type sllist
		type(sllist), pointer 		  :: child => null()
		character(len=:), allocatable :: key
		integer						  :: val
	contains
		procedure :: put  => put_sll
		procedure :: get  => get_sll
		procedure :: free => free_sll
	end type sllist
	
	type hash_tbl_sll
		type(sllist), dimension(:), allocatable :: vec
		integer									:: vec_len = 0
		logical									:: is_init = .false.
	contains
		procedure :: init => init_hash_tbl_sll
		procedure :: put  => put_hash_tbl_sll
		procedure :: get  => get_hash_tbl_sll
		procedure :: free => free_hash_tbl_sll
	end type hash_tbl_sll
	
	public :: hash_tbl_sll
contains
	recursive subroutine put_sll(list, key, val)
		class(sllist), intent(inout)	:: list
		character(len=*), intent(in)	:: key
		integer, intent(in)				:: val
		integer							:: keylen
		
		keylen = len(key)
		if (allocated(list%key)) then
			if (list%key /= key) then
				if (.not. associated(list%child)) allocate(list%child)
				call put_sll(list%child, key, val)
			end if
		else
			if (.not. allocated(list%key)) allocate(character(len=keylen) :: list%key)
			list%key = key
			list%val = val
		end if
	end subroutine put_sll
	
	recursive subroutine get_sll(list, key, val)
		class(sllist), intent(in)		:: list
		character(len=*), intent(in)	:: key
		integer, intent(out)			:: val
	
		if (allocated(list%key) .and. (list%key == key)) then
			val = list%val
		else if (associated(list%child)) then
			call get_sll(list%child, key, val)
		else
			val = 0
		end if 
	end subroutine get_sll
	
	recursive subroutine free_sll(list)
		class(sllist), intent(inout)	:: list
		if (associated(list%child)) then
			call free_sll(list%child)
			deallocate(list%child)
		end if
		list%child => null()
		if (allocated(list%key)) deallocate(list%key)
	end subroutine free_sll
	
	subroutine init_hash_tbl_sll(tbl, tbl_len)
		class(hash_tbl_sll), intent(inout)	:: tbl
		integer, optional, intent(in)		:: tbl_len
		
		if (allocated(tbl%vec)) deallocate(tbl%vec)
		if (present(tbl_len)) then
			allocate(tbl%vec(0:tbl_len-1))
			tbl%vec_len = tbl_len
		else
			allocate(tbl%vec(0:tbl_size-1))
			tbl%vec_len = tbl_size
		end if
		tbl%is_init = .true.
	end subroutine init_hash_tbl_sll
	
	integer function gen_hash(str) result(hash)
		character(len=*), intent(in) :: str
		integer :: i
		
		hash = 5381
		do i=1,len(str)
			hash = (ishft(hash,5) + hash) + ichar(str(i:i))
		end do
	end function gen_hash
	
	subroutine put_hash_tbl_sll(tbl, key, val)
		class(hash_tbl_sll), intent(inout)	:: tbl
		character(len=*), intent(in)		:: key
		integer, intent(in)					:: val
		integer								:: hash
		
		hash = abs(mod(gen_hash(key), tbl%vec_len))
		call tbl%vec(hash)%put(key=key, val=val)
	end subroutine put_hash_tbl_sll
	
	subroutine get_hash_tbl_sll(tbl, key, val)
		class(hash_tbl_sll), intent(in)		:: tbl
		character(len=*), intent(in)		:: key
		integer, intent(out)				:: val
		integer								:: hash
		
		hash = abs(mod(gen_hash(key), tbl%vec_len))
		call tbl%vec(hash)%get(key=key, val=val)
	end subroutine get_hash_tbl_sll
	
	subroutine free_hash_tbl_sll(tbl)
		class(hash_tbl_sll), intent(inout) :: tbl
		integer							   :: i, low, high
		
		low  = lbound(tbl%vec, dim=1)
		high = ubound(tbl%vec, dim=1)
		if (allocated(tbl%vec)) then
			do i=low,high
				call tbl%vec(i)%free()
			end do
			deallocate(tbl%vec)
		end if
		tbl%is_init = .false.
	end subroutine free_hash_tbl_sll
	
end module hashtbl