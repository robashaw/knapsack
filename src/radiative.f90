module radiative
	use constants
	implicit none
	
contains
	
	subroutine trapezium(value, x, f, N)
		integer, intent(in) 			    :: N
		real(dbl), dimension(N), intent(in)	:: x, f
		real(dbl), intent(out)				:: value
		
		integer		:: i
		real(dbl)	:: width, tmp
		value = 0d0
		do i=1, N-1
			width = x(i+1) - x(i)
			tmp = 0.5 * width * (f(i) + f(i+1))
			value = value + tmp
		end do
	end subroutine trapezium
	
	subroutine calculate_rate(value, filename, eunits)
		real(dbl), intent(out)			:: value
		character(len=*), intent(in)	:: filename
		
		integer 	:: nlines, nvalues, firstline, ios, line, i
		real(dbl)	:: x, f, norm, unit_conversion
		real(dbl), allocatable, dimension(:)  :: xvals, fvals
		character(len=100)	:: iom
		
		write(*, *) 'Calculating radiative rate'
		write(*, *) 'Reading spectral data from ', filename
		write(*, *) 'Assuming energy units of ', eunits
		
		select case(eunits)
		case ('kj')
			unit_conversion = KJ_TO_EV
		case ('ha')
			unit_conversion = HA_TO_EV
		case ('kcal')
			unit_conversion = KCAL_TO_EV
		case ('cm')
			unit_conversion = CM_TO_EV
		case default
			unit_conversion = 1d0 ! assume in eV already
		end select	
		
		open(radiative_unit, file=filename, iostat=ios)
		firstline = 0
		f = 0d0
		nlines = -1
		do while (ios == 0)
			nlines = nlines + 1
			if (firstline == 0) then
				if (abs(f) > tolint) firstline = nlines
			end if
			read(radiative_unit, *, iostat=ios, iomsg=iom) x, f
		end do
		nvalues = nlines - firstline + 1
		if (nvalues < 1) then
			write(*, '(3A)') 'File ', filename, ' is empty or all zero'
		else
			rewind(radiative_unit)
			allocate(xvals(nvalues))
			allocate(fvals(nvalues))
			
			i = 1
			do line = 1, nlines
				read(radiative_unit, *, iostat=ios, iomsg=iom) x, f
				if (line .ge. firstline) then 
					xvals(i) = x * unit_conversion ! make sure in eV
					fvals(i) = f
					i = i + 1
				end if
			end do
			
			call trapezium(norm, xvals, fvals, nvalues)
			write(*, *) 'Spectral density norm:', norm
			if (norm .gt. tolint) then
				fvals = fvals / norm
				do i = 1, nvalues
					fvals(i) = fvals(i) * (xvals(i)**3)
				end do
				call trapezium(value, xvals, fvals, nvalues)
			else
				write(*, *) 'Aborting radiative calculation, spectral density has negligible norm'
			end if 
			
			deallocate(xvals)
			deallocate(fvals)
		end if
		close(radiative_unit)
		
	end subroutine 
end module radiative