! structured_grids.f03
!
! Module that contains user-defined types and subroutines
! for 2D and 3D grids.

module structured_grids
    use, intrinsic iso_fortran_env
    implicit none
    
    ! Parameters for identifying boundary types.
    integer, parameter :: DIRICHLET = 0
    integer, parameter :: NEUMANN = 1
    integer, parameter :: ROBIN = 2
    
    contains
        
        !-----------------------------------------------------------!
        ! Type definition and subroutines for 2D rectangular grids. !
        !-----------------------------------------------------------!
        
        type :: Grid2D
            ! Coordinates of the bounding domain.
            real(REAL64) :: x0, xn, y0, yn
            
            ! The initial number of points in each direction.
            integer :: numx0, numy0
            
            ! The number of points at the current grid level.
            integer :: numx, numy
            
            ! The current step size on the given grid.
            integer :: hx, hy
            
            ! The current grid level and the total number of grids.
            integer :: current_grid_level
            integer :: num_grids
            
            ! Array for the boundary types.
            integer, dimension(1:4) :: boundary_flags
        end type Grid2D
        
        subroutine grid2d_init(grid, x0, xn, y0, yn, nx0, ny0, ngrids, bflags)
            type(Grid2D), intent(inout) :: grid
            real(REAL64), intent(in) :: x0, xn, y0, yn
            integer, intent(in) :: nx0, ny0, ngrids
            integer, dimension(1:4), intent(in) :: bflags
            
            ! Integers for looping over and initializing grid indices.
            integer :: i, j
            
            ! Set the endpoints.
            grid%x0 = x0
            grid%xn = xn
            grid%y0 = y0
            grid%yn = yn
            
            ! Set the initial number of points in each direction.
            grid%numx0 = nx0
            grid%numy0 = ny0
            
            ! Set the current grid number to the maximum
            ! number of grids, and the number of x- and y-
            ! points accordingly.
            grid%current_grid_level = ngrids - 1
            grid%num_grids = ngrids
            grid%numx = 2**(ngrids - 1) * (nx0 - 1) + 1
            grid%numy = 2**(ngrids - 1) * (ny0 - 1) + 1
            
            ! Now set the step sizes.
            grid%hx = (grid%xn - grid%x0) / (grid%numx - 1)
            grid%hy = (grid%yn - grid%y0) / (grid%numy - 1)
            
            ! Set the boundary nodes to dirichlet by default.
            grid%boundary_flags = bflags
        end subroutine init_grid2d
        
        subroutine set_boundary_types_2d(grid, node_types, houndary_num)
            type(Grid2D), intent(inout) :: grid
            integer, dimension(:), intent(in) :: node_types
            integer, dimension(:), intent(in) :: houndary_nums
            
            integer :: k
            
            if(allocated(grid%node_flags)) then
                do k = 1, size(boundary_nums)
                    select case (boundary_nums(k))
                        case (1)
                            grid%node_flags(grid%nx, :) = node_types(k)
                        case (2)
                            grid%node_flags(:, grid%ny) = node_types(k)
                        case (3)
                            grid%node_flags(1, :) = node_types(k)
                        case (4)
                            grid%node_flags(:, 1) = node_types(k)
                        case default
                            print *, 'Error: Array index out of bounds for the number of boundaries.'
                    end select
                end do
            end if
            
        end subroutine set_boundary_types_2d
        
        subroutine grid2d_refine(grid)
            type(Grid2D), intent(inout) :: grid
            
            grid2d_jump_to_level(grid, grid%current_level + 1)
        end subroutine grid2d_refine
        
        subroutine grid2d_coarsen(grid)
            type(Grid2D), intent(inout) :: grid
            
            grid2d_jump_to_level(grid, grid%current_level - 1)
        end subroutine grid2d_coarsen
        
        subroutine grid2d_jump_to_level(grid, k)
            type(Grid2D), intent(inout) :: grid
            integer, intent(in) :: k
            
            ! Integers for looping over and initializing grid indices.
            integer :: i, j
            
            ! We need to have a valid grid number;
            ! otherwise, we do nothing and exit.
            if(0 <= k .and. k < grid%num_grids - 1) then
                ! Set the current grid level.
                grid%current_level = k
                
                ! Set the number of x- and y-points.
                grid%numx = 2**k * (grid%numx0 - 1) + 1
                grid%numy = 2**k * (grid%numy0 - 1) + 1
                
                ! Now set the step sizes.
                grid%hx = (grid%xn - grid%x0) / (grid%numx - 1)
                grid%hy = (grid%yn - grid%y0) / (grid%numy - 1)
            endif
        end subroutine update_grid2d_info
        
        !-----------------------------------------------------------!
        ! Type definition and subroutines for 3D rectangular grids. !
        !-----------------------------------------------------------!
        
        type :: Grid3D
            ! Coordinates of the bounding domain.
            real(REAL64) :: x0, xn, y0, yn, z0, zn
            
            ! The initial number of points in each direction.
            integer :: nx0, ny0, nz0
            
            ! The number of points at the current grid level.
            integer :: nx, ny, nz
            
            ! The current step size on the given grid.
            integer :: hx, hy, hz
            
            ! The current (kth) grid and the total number of grids m.
            integer :: k, m
            
            ! Array for the boundary types.
            integer, dimension(1:8) :: boundary_flags = DIRICHLET_NODE
        end type Grid3D
        
        subroutine init_grid3d(grid, x0, xn, y0, yn, z0, zn, nx0, ny0, nz0, num_grids)
            type(Grid3D), intent(inout) :: grid
            real(REAL64), intent(in) :: x0, xn, y0, yn, z0, zn
            integer, intent(in) :: nx0, ny0, nz0, num_grids
            
            ! Set the endpoints.
            grid%x0 = x0
            grid%xn = xn
            grid%y0 = y0
            grid%yn = yn
            grid%z0 = z0
            grid%zn = zn
            
            ! Set the initial number of points in each direction.
            grid%nx0 = nx0
            grid%ny0 = ny0
            grid%nz0 = nz0
            
            ! Set the current grid number to the maximum
            ! number of grids, and the number of x- and y-
            ! points accordingly.
            grid%k = num_grids
            grid%m = num_grids
            grid%nx = 2**(grid%m - 1) * (grid%nx0 - 1) + 1
            grid%ny = 2**(grid%m - 1) * (grid%ny0 - 1) + 1
            grid%nz = 2**(grid%m - 1) * (grid%nz0 - 1) + 1
            
            ! Now set the step sizes.
            grid%hx = (grid%xn - grid%x0) / (grid%nx - 1)
            grid%hy = (grid%yn - grid%y0) / (grid%ny - 1)
            grid%hz = (grid%zn - grid%z0) / (grid%nz - 1)
        end subroutine init_grid3d
        
        subroutine update_grid3d_info(grid, k)
            type(Grid3D), intent(inout) :: grid
            integer, intent(in) :: k
            
            ! We need to have a valid grid number;
            ! otherwise, we do nothing and exit.
            if(k .lt. 1 .or. k .gt. grid%m) then
                print *, 'Invalid grid number.'
                
                return;
            endif
            
            ! Set the current grid level.
            grid%k = k
            
            ! Set the number of x- and y-points.
            grid%nx = 2**(grid%k - 1) * (grid%nx0 - 1) + 1
            grid%ny = 2**(grid%k - 1) * (grid%ny0 - 1) + 1
            grid%nz = 2**(grid%k - 1) * (grid%nz0 - 1) + 1
            
            ! Now set the step sizes.
            grid%hx = (grid%xn - grid%x0) / (grid%nx - 1)
            grid%hy = (grid%yn - grid%y0) / (grid%ny - 1)
            grid%hz = (grid%zn - grid%z0) / (grid%nz - 1)
        end subroutine update_grid3d_info
        
end module structured_grids






