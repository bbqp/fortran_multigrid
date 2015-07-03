! precision_specs.f08
!
! Integer parameters for precision specification.
! This file contains single, double, and quad
! precision when available.

module precision_specs
	use, intrinsic :: iso_fortran_env
	implicit none
	
	integer, parameter :: R32 = REAL32;
	integer, parameter :: R64 = REAL64;
	integer, parameter :: R128 = REAL128;	! WARNING: May only be a 10-byte extended real
											! on some systems, normally 16 bytes when available.
end module precision_specs