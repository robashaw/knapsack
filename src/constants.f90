module constants
	use iso_fortran_env, only: real64, int64
	
	implicit none
	integer, parameter 		:: dbl = real64
	integer, parameter		:: bigint = int64
	integer, parameter		:: maxnchange = 10
	integer, parameter		:: main_input_unit = 15
	integer, parameter		:: bfile_unit = 16
	integer, parameter		:: gradfile_unit = 17
	integer, parameter		:: radiative_unit = 18
	integer, parameter		:: occs_unit = 19
	integer, parameter		:: n_guesses = 10
	integer, parameter		:: n_cut_guess = 4	
	integer, parameter		:: print_frequency = 100000
	real(dbl), parameter	:: tolint = 1D-12
	real(dbl), parameter 	:: damping = 1.0
	real(dbl), parameter	:: PLANCK = 4.135667696D-15 ! in eV.s
	real(dbl), parameter	:: SQRT_8LN2 = 2.354820045
	real(dbl), parameter	:: KB = 8.617333262D-5 ! in eV/K
	real(dbl), parameter	:: TO_ME = 1.822888486D+4
	real(dbl), parameter	:: TO_EV = 27.21138602 ! Ha to eV
	real(dbl), parameter	:: TO_S  = 2.418884326509D-17 ! a.u. to s
	real(dbl), parameter	:: TO_CM = 8065.5 ! eV to cm-1
	real(dbl), parameter	:: KJ_TO_EV = 0.0103641 ! kJ/mol to eV
	real(dbl), parameter	:: CM_TO_EV = 1.23981D-4 ! cm-1 to eV
	real(dbl), parameter	:: KCAL_TO_EV = 0.0433634 ! kcal/ol to eV
	real(dbl), parameter	:: HA_TO_EV = 27.2107 ! Hartree to eV
	real(dbl), parameter	:: V_CONVERT = 1.49198625108D+5 ! TO_ME * sqrt(TO_ME/TO_EV)
	real(dbl), parameter	:: RAD_CONVERT = 1.0781987879D+6
	
contains
	real(dbl) function factorial(n) result(fac)
	! Only really appropriate up to about n=20, but unlikely to need higher
		integer, intent(in) :: n
		
		integer :: ix
		fac = 1d0
		if (n .gt. 0) then
			do ix = 1, n
				fac = fac * real(n)
			end do
		end if
	end function factorial
	
	function format_time(seconds) result(t)
		real(dbl), intent(in)	:: seconds 
		character(len=20)		:: t
		
		integer	:: secs, mins, hrs
		hrs = floor(seconds / 3600d0)
		mins = floor((seconds - hrs*3600d0) / 60d0)
		secs = floor(seconds - hrs*3600d0 - mins*60d0)
		write(t, '(i4,a2,i2,a2,i2,a1)') hrs, 'h:', mins, 'm:', secs, 's'
	end function format_time
		
end module