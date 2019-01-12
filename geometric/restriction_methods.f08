module restriction_methods

	use precision_specs
	use structured_grids
	implicit none
	
	interface inject
		module procedure inject_2d
		module procedure inject_3d
	end interface inject

	interface half_weight
		module procedure half_weight2d
		module procedure half_weight3d
	end interface half_weight
	
	interface full_weight
		module procedure full_weight2d
		module procedure full_weight3d
	end interface full_weight

	contains
	
	subroutine inject_2d(grid, UC, UF)
		type(Grid2D), intent(in) :: grid
		real(R64), dimension(1:((grid%nx + 1) / 2), 1:((grid%ny + 1) / 2)), intent(in) :: UC
		real(R64), dimension(1:grid%nx, 1:grid%ny), intent(in) :: UF
		
		! Set the coinciding fine and coarse points.
		UC = UF(1:size(UF, 1):2, 1:size(UF, 2):2)
	end subroutine
	
	subroutine inject_3d(grid, UC, UF)
		type(Grid3D), intent(in) :: grid
		real(R64), dimension(1:((grid%nx + 1) / 2), 1:((grid%ny + 1) / 2), 1:((grid%nz + 1) / 2)), intent(in) :: UC
		real(R64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz), intent(in) :: UF
		
		! Set the coinciding fine and coarse points.
		UC = UF(1:size(UF, 1):2, 1:size(UF, 2):2)
	end subroutine
	
	!
	! Half-weight restruction.
	!
	subroutine half_weight_2d(grid, UC, UF)
		type(Grid3D), intent(in) :: grid
		real(R64), dimension(1:((grid%nx + 1) / 2), 1:((grid%ny + 1) / 2), 1:((grid%nz + 1) / 2)), intent(in) :: UC
		real(R64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz), intent(in) :: UF
		
		integer :: i, j
		
		! First restrict the interior points.
		do j = 2, SIZE(UC, 2) - 1
			do i = 2, SIZE(UC, 1) - 1
				UC(i, j) = 0.125 * ( &
					4.0 * UF(2*i - 1, 2*j - 1) + &
					
					UF(2*i - 2, 2*j - 1) + &
					UF(2*i, 2*j - 1) + &
					
					UF(2*i - 1, 2*j - 2) + &
					UF(2*i - 1, 2*j) &
				)
			end do
		end do
		
		! First restrict the x-valued boundaries for fixed y-values.
		do i = 2, SIZE(UC, 1)
			UC(i, 1) = 0.25 * ( UF(2*i - 2, 1) +  2.0 * UF(2*i - 1, 1) + UF(2*i, 1) )
			UC(i, SIZE(UC, 2)) = 0.25 * ( UF(2*i - 2, SIZE(UF, 2)) +  2.0 * UF(2*i - 1, SIZE(UF, 2)) + UF(2*i, SIZE(UF, 2)) )
		end do
		
		! Now restrict the y-valued boundaries for fixed x-values.
		do j = 2, SIZE(UC, 2)
			UC(1, j) = 0.25 * ( UF(1, 2*j - 2) +  2.0 * UF(1, 2*j - 1) + UF(1, 2*j) )
			UC(SIZE(UC, 1), j) = 0.25 * ( UF(1, 2*j - 2) +  2.0 * UF(1, 2*j - 1) + UF(1, 2*j) )
		end do
		
	end subroutine half_weight_2d
	
	!
	! Full weighting restriction in 2 variables.
	!
	subroutine full_weight2d(grid, UC, UF)
		type(Grid2D), intent(in) :: grid
		real(R64), dimension(1:((grid%nx + 1) / 2), 1:((grid%ny + 1) / 2)), intent(in) :: UC
		real(R64), dimension(1:grid%nx, 1:grid%ny), intent(in) :: UF
		
		UC = 0.0
		
		! First weight the interior points.
		do j = 2, SIZE(UC, 2) - 1
			do i = 2, SIZE(UC, 1) - 1
				UC(i, j) =	0.0625 *	( &
											4.0 * UF(2*i - 1, 2*j - 1) + &
											
											2.0 * ( &
												UF(2*i - 2, 2*j - 1) + &
												UF(2*i, 2*j - 1) + &
												UF(2*i - 1, 2*j - 2) + &
												UF(2*i - 1, 2*j) + &
											) + &
											
											UF(2*i - 2, 2*j - 2) + &
											UF(2*i - 2, 2*j) + &
											UF(2*i, 2*j - 2) + &
											UF(2*i, 2*j) + &
										)
			end do
		end do
		
		! First restrict the x-valued boundaries for fixed y-values.
		do i = 2, SIZE(UC, 1)
			UC(i, 1) = 0.25 * ( UF(2*i - 2, 1) +  2.0 * UF(2*i - 1, 1) + UF(2*i, 1) )
			UC(i, SIZE(UC, 2)) = 0.25 * ( UF(2*i - 2, SIZE(UF, 2)) +  2.0 * UF(2*i - 1, SIZE(UF, 2)) + UF(2*i, SIZE(UF, 2)) )
		end do
		
		! Now restrict the y-valued boundaries for fixed x-values.
		do j = 2, SIZE(UC, 2)
			UC(1, j) = 0.25 * ( UF(1, 2*j - 2) +  2.0 * UF(1, 2*j - 1) + UF(1, 2*j) )
			UC(SIZE(UC, 1), j) = 0.25 * ( UF(1, 2*j - 2) +  2.0 * UF(1, 2*j - 1) + UF(1, 2*j) )
		end do
		
	end subroutine full_weight2d
		
	!
	! A vectorized version of the above subroutine (i.e. without do loops).
	!
	subroutine vectorized_full_weight_2d(grid, UC, UF)
		type(Grid2D), intent(in) :: grid
		real(R64), dimension(1:((grid%nx + 1) / 2), 1:((grid%ny + 1) / 2)), intent(in) :: UC
		real(R64), dimension(1:grid%nx, 1:grid%ny), intent(in) :: UF
		
		! Initialize the course restriction.
		UC = 0.0;
		
		! Now perform the restriction of the interior points.
		UC(2:(SIZE(UC, 1) - 1), 2:(SIZE(UC, 2) - 1)) = &
			0.0625 * ( &
				4.0 * UF(3:2:(SIZE(UF, 1) - 2), 3:2:(SIZE(UF, 2) - 2)) + &
				
				2.0 * ( &
					UF(2:(SIZE(UF, 1) - 3):2, 3:(SIZE(UF, 2) - 2):2) + &
					UF(4:(SIZE(UF, 1) - 1):2, 3:(SIZE(UF, 2) - 2):2) + &
					UF(3:(SIZE(UF, 1) - 2):2, 2:(SIZE(UF, 2) - 3):2) + &
					UF(3:(SIZE(UF, 1) - 2):2, 4:(SIZE(UF, 2) - 1):2) &
				) + &
				
				UF(2:(SIZE(UF, 1) - 3):2, 2:(SIZE(UF, 2) - 3):2) + &
				UF(2:(SIZE(UF, 1) - 3):2, 4:(SIZE(UF, 2) - 1):2) + &
				UF(4:(SIZE(UF, 1) - 1):2, 2::2(SIZE(UF, 2) - 3):2) + &
				UF(4:(SIZE(UF, 1) - 1):2, 4:(SIZE(UF, 2) - 1):2) &
			)
		
		! Now restrict the boundaries. In theory, dirichlet boundaries imply
		! that the residual is equal to zero, but we'll do the restriction anyway.
		
		! Restrict the x-valued portion first.
		UC(:, 1) = 0.5 * ( UF(2:2:(SIZE(UF, 1) - 3), 1) + UF(4:2:(SIZE(UF, 1) - 1), 1) )
		UC(:, SIZE(UC, 2)) = ( UF(2:2:(SIZE(UF, 1) - 3), SIZE(UC, 2)) + UF(4:2:(SIZE(UF, 1) - 1), SIZE(UC, 2)) )
		
		! Now restrict the y-valued portion.
		UC(1, :) = 0.5 * ( UF(1, 2:(SIZE(UF, 2) - 3):2) + UF(1, 4:(SIZE(UF, 2) - 1):2) )
		UC(SIZE(UC, 1), :) = 0.5 * ( UF(SIZE(UC, 1), 2:(SIZE(UF, 2) - 3):2) + UF(SIZE(UC, 1), 4:(SIZE(UF, 2) - 1):2) )
		
	end subroutine vectorized_full_weight_2d
	
	!
	! Full weighting restriction for a function of three variables.
	!
	subroutine full_weight3d(grid, UC, UF)
		type(Grid3D), intent(in) :: grid
		real(R64), dimension(1:((grid%nx + 1) / 2), 1:((grid%ny + 1) / 2), 1:((grid%nz + 1) / 2)), intent(in) :: UC
		real(R64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz), intent(in) :: UF
		
		integer :: i, j, k
		
		! Perform full-weighting on the interior points.
		do k = 2, SIZE(UC, 3) - 1
			do j = 2, SIZE(UC, 2) - 1
				do i = 2, SIZE(UC, 1) - 1
					
					UC(i, j, k) = ( &
						4.0 * UF(2*i - 1, 2*j - 1, 2*k - 1) + &
						
						2.0 * ( &
								UF(2*i - 2, 2*j - 1) + &
								UF(2*i, 2*j - 1) + &
								UF(2*i - 1, 2*j - 2) + &
								UF(2*i - 1, 2*j) + &
						) + &
						
						UF(2*i - 2, 2*j - 2) + &
						UF(2*i - 2, 2*j) + &
						UF(2*i, 2*j - 2) + &
						UF(2*i, 2*j) + &
					) / 64.0
					
				end do
			end do
		end do
	end subroutine full_weight3d
	
	!
	! Full weighting restriction for a function of three variables,
	! but sans the do loops of the above subroutine
	!
	subroutine vetorized_full_weight3d(grid, UC, UF)
		type(Grid3D), intent(in) :: grid
		real(R64), dimension(1:((grid%nx + 1) / 2), 1:((grid%ny + 1) / 2), 1:((grid%nz + 1) / 2)), intent(in) :: UC
		real(R64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz), intent(in) :: UF
		
		
	end subroutine vetorized_full_weight3d
end module











