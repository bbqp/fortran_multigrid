! precision_specs.f08
!
! Integer parameters for precision specification.
! This file contains single, double, and quad
! precision when available.

module precision_specs

	integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = selected_real_kind(2 * precision(1.0_sp))
    integer, parameter :: qp = selected_real_kind(2 * precision(1.0_dp))
	
end module precision_specs