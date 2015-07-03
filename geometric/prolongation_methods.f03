module prolongation_methods
	use structured_grids
	use unstructured_grids
	use precision_specs
	implicit none
	
	interface linear_interp
		module procedure linear_interp_2d
		module procedure linear_interp_3d
	end interface linear_interp
	
	interface cubic_interp
		module procedure cubic_interp_2d
		module procedure cubic_interp_3d
	end interface cubic_interp
	
	contains
	
	subroutine linear_interp_2d(grid, UC, UF)
		type(Grid2D), intent(in) :: grid
		real(R64), dimension(1:grid%nx, 1:grid%ny), intent(in) :: UC
		real(R64), dimension(1:(2 * grid%nx - 1), 1:(2 * grid%ny - 1)), intent(inout) :: UF
		
		! Map the course points to the fine points.
		UF(1:size(UF, 1):2, 1:size(UF, 2):2) = UC
		
		! Interpolate the points for each x-line.
		UF(2:(size(UF, 1) - 1):2, 3:(size(UF, 2) - 2):2) = 0.5 *	( &
																		UF(1:(size(UF, 1) - 2):2, 3:(size(UF, 2) - 2):2) + &
																		UF(3:size(UF, 1):2, 3:(size(UF, 2) - 2):2) &
																	)
		
		! Interpolate the points for each y-line.
		UF(3:(size(UF, 1) - 2):2, 2:(size(UF, 2) - 1):2) = 0.5 *	( &
																		UF(3:(size(UF, 1) - 2):2, 1:(size(UF, 2) - 2):2) + &
																		UF(3:(size(UF, 1) - 2):2, 3:size(UF, 2):2) &
																	)
		
		! Interpolate the interior points.
		UF(2:(size(UF, 1) - 1):2, 2:(size(UF, 2) - 1):2) = 0.25 *	( &
																		UF(1:(size(UF, 1) - 2):2, 1:(size(UF, 2) - 2):2) + &
																		UF(1:(size(UF, 1) - 2):2, 3:size(UF, 2):2) + &
																		UF(3:size(UF, 1):2, 1:(size(UF, 2) - 2):2) + &
																		UF(3:size(UF, 1):2, 3:size(UF, 2):2) &
																	)
		
		! TODO: Interpolation of boundaries.
	end subroutine linear_interp_2d
	
	subroutine linear_interp_3d(grid, UC, UF)
		type(Grid3D), intent(in) :: grid
		real(R64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz), intent(in) :: UC
		real(R64), dimension(1:(2 * grid%nx - 1), 1:(2 * grid%ny - 1), 1:(2 * grid%nz - 1)), intent(inout) :: UF
		
		! Map the course points to the fine points.
		UF(1:size(UF, 1):2, 1:size(UF, 2):2, 1:size(UF, 3):2) = UC
		
		! Interpolate points for each coarse x plane.
		xplane_linear_interp_3d(grid, UC, UF)
		
		! Interpolate points for each coarse y plane.
		UF(3:(size(UF, 1) - 2):2, 2:(size(UF, 2) - 1):2) = 0.5 *	( &
																		UF(3:(size(UF, 1) - 2):2, 1:(size(UF, 2) - 2):2) + &
																		UF(3:(size(UF, 1) - 2):2, 3:size(UF, 2):2) &
																	)
		
		! Interpolate points for each coarse z plane.
		
		! Interpolate the interior points.
		UF(2:(size(UF, 1) - 1):2, 2:(size(UF, 2) - 1):2, 2:(size(UF, 3) - 1):2) = &
		
		0.125 * ( &
			UF(1:(size(UF, 1) - 2):2, 1:(size(UF, 2) - 2):2, 1:(size(UF, 3) - 2):2) + &
			UF(1:(size(UF, 1) - 2):2, 3:size(UF, 2):2, 1:(size(UF, 3) - 2):2) + &
			UF(3:size(UF, 1):2, 1:(size(UF, 2) - 2):2, 1:(size(UF, 3) - 2):2) + &
			UF(3:size(UF, 1):2, 3:size(UF, 2):2, 1:(size(UF, 3) - 2):2) + &
			
			UF(1:(size(UF, 1) - 2):2, 1:(size(UF, 2) - 2):2, 3:size(UF, 3):2) + &
			UF(1:(size(UF, 1) - 2):2, 3:size(UF, 2):2, 3:size(UF, 3):2) + &
			UF(3:size(UF, 1):2, 1:(size(UF, 2) - 2):2, 3:size(UF, 3):2) + &
			UF(3:size(UF, 1):2, 3:size(UF, 2):2, 3:size(UF, 3):2) &
		)
		
		! TODO: Interpolation of boundaries.
	end subroutine linear_interp_3d
	
	!
	! Subroutine which interpolates the coarse planes of the solution. 
	!
	subroutine xplane_linear_interp_3d(grid, UC, UF)
		type(Grid3D), intent(in) :: grid
		real(R64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz), intent(in) :: UC
		real(R64), dimension(1:(2 * grid%nx - 1), 1:(2 * grid%ny - 1), 1:(2 * grid%nz - 1)), intent(inout) :: UF
		
		!
		! For each z-level, interpolate the x-lines, y-lines, and interior points.
		!
		
		! First compute the fine x-line values.
		UF(2:(SIZE(UF, 1) - 1):2, 3:(SIZE(UF, 2) - 2):2, 3:(SIZE(UF, 3) - 2):2) = &
		
		0.5 * ( &
			UF(1:(SIZE(UF, 1) - 2):2, 3:(SIZE(UF, 2) - 2):2, 3:(SIZE(UF, 3) - 2):2) + &
			UF(3:SIZE(UF, 1):2, 3:(SIZE(UF, 2) - 2):2, 3:(SIZE(UF, 3) - 2):2) &
		)
		
		! Now compute the fine y-line values.
		UF(3:(SIZE(UF, 1) - 2):2, 2:(SIZE(UF, 2) - 1):2, 3:(SIZE(UF, 3) - 2):2) = &
		
		0.5 * ( &
			UF(3:(SIZE(UF, 1) - 2):2, 2:(SIZE(UF, 2) - 1):2, 3:(SIZE(UF, 3) - 2):2) + &
			UF(3:(SIZE(UF, 1) - 2):2, 2:(SIZE(UF, 2) - 1):2, 3:(SIZE(UF, 3) - 2):2) &
		)
		
		! Now compute the interior points
		UF(2:(SIZE(UF, 1) - 1):2, 2:(SIZE(UF, 1) - 1):2, 3:(SIZE(UF, 3) - 2):2) = &
		
		0.25 * ( &
			UF(i - 1, j - 1, k) + &
			UF(i - 1, j + 1, k) + &
			UF(i + 1, j - 1, k) + &
			UF(i + 1, j + 1, k) + &
		)
	end subroutine xplane_linear_interp_3d
	
end module













