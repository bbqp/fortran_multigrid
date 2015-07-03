! unstructured_grids.f03

! Module that contains user-defined types and
! subroutines for 2D and 3D unstructured grids.

module grid_types
	use precision_specs
	
	implicit none
	
	! Parameters for identifying node types.
	integer, parameter :: FREE_NODE = 0			! Also refers to an interior node for which we are solving.
	integer, parameter :: DIRICHLET_NODE = 1	! A node with a fixed value on the boundary.
	integer, parameter :: NEUMANN_NODE = 2		! A node corresponding to flux on the boundary.
	integer, parameter :: ROBIN_NODE = 3		! A node whose value is a combination of the previous two.
	
	! Integer parameters for edge types.
	integer, parameter :: FREE_EDGE	= 4			! An edge in the interior of the domain.
	integer, parameter :: DIRICHLET_EDGE = 5	! An edge with a fixed (function) value on the boundary.
	integer, parameter :: NEUMANN_EDGE = 6		! An edge for flux on the boundary
	integer, parameter :: ROBIN_EDGE = 7		! An edge with a combination of the previous two conditions.
	
	! Integer parameters for edge types.
	integer, parameter :: FREE_FACE	= 8			! An edge in the interior of the domain.
	integer, parameter :: DIRICHLET_FACE = 9	! An edge with a fixed (function) value on the boundary.
	integer, parameter :: NEUMANN_FACE = 10		! An edge for flux on the boundary
	integer, parameter :: ROBIN_FACE = 11		! An edge with a combination of the previous two conditions.
	
	contains
		
		!-------------------------------------------------!
		! Type definition and subroutines for arbitrary   !
		! (unstructured) grids. These grids will be       !
		! particularly useful for finite element methods. !
		! in combination with algebraic multigrid.		  !
		!-------------------------------------------------!
		
		type UnstructuredGrid(space_dim)
			real(R64), dimension(:,:), allocatable :: nodes
			integer, dimension(:,:), allocatable :: edges
			integer, dimension(:,:), allocatable :: elems
			integer, dimension(:), allocatable :: node_ordering
			integer, dimension(:), allocatable :: node_types
			
			integer :: num_nodes
			integer :: num_edges
			integer :: num_elems
		end type UnstructuredGrid
		
		interface init_unstructured_grid
		
			subroutine init_unstructured_default_grid(grid)
				implicit none
				
				type(UnstructuredGrid), intent(in) :: grid
			end subroutine init_unstructured_default_grid
			
			subroutine init_unstructured_grid_from_file(grid, gridfile)
				implicit none
				
				type(UnstructuredGrid), intent(inout) :: grid
				character(:), intent(in) :: gridfile
			end subroutine init_unstructured_grid_from_file
			
		end interface init_unstructured_grid
		
		! Now define the methods contained in the interface for
		! initializing the unstructured grids.
		
		subroutine init_unstructured_default_grid(grid)
			implicit none
			
			type(UnstructuredGrid), intent(in) :: grid
			
			! Looping variables for adding points and setting elements, edges.
			integer :: i, j
			
			! Condition to indicate a point is on a boundary.
			logical :: on_boundary
			
			! Here, our "unstructured" domain will be the unit square
			! with a default step size of 1/8 (default of 9 interior
			! points in each direction).
			real(R64) :: x0, xn
			real(R64) :: y0, yn
			real(R64) :: x, y
			integer :: nx, ny
			real(R64) :: hx, hy
			
			x0 = 0
			xn = 1
			y0 = 0
			yn = 1
			
			nx = 33
			ny = 33
			
			hx = (xn - x0) / (nx - 1)
			hy = (yn - y0) / (ny - 1)
			
			! Now set the fields in the grid.
			grid%num_nodes = nx * ny
			grid%num_edges = 3 * (nx - 1) * (ny - 1) + (nx - 1) + (ny - 1)
			grid%num_elems = 2 * (nx - 1) * (ny - 1)
			
			ALLOCATE(grid%nodes(grid%num_nodes, 2))
			ALLOCATE(grid%edges(grid%num_edges, 2))
			ALLOCATE(grid%elems(grid%num_elems, 2))
			ALLOCATE(grid%node_ordering(grid%num_nodes))
			ALLOCATE(grid%node_types(grid%num_nodes))
			
			! Initialize the nodes and ordering.
			do j = 1, ny
				y = y0 + (j - 1) * h
				
				do i = 1, nx
					x = x0 + (i - 1) * h
					
					grid%nodes(nx * (j - 1) + i, :) = (/ x, y /)
					
					! Check if we are on a boundary.
					on_boundary = 	(
										(abs(x - x0) .lt. 1.0D-15) .or.
										(abs(x - xn) .lt. 1.0D-15)
									)
									.and.
									(
										(abs(y - y0) .lt. 1.0D-15) .and.
										(abs(y - yn) .lt. 1.0D-15)
									)
									
					if(on_boundary) then
						grid%node_types(nx * j + i) = DIRICHLET_NODE
					else
						grid%node_types(nx * j + i) = FREE_NODE
					end if
				end do
			end do
			
			! Now initialize the edges.
			do j = 1, ny - 1
				do i = 1, nx - 1
					grid%edges
				end do
			end do
			
			
		end subroutine init_unstructured_default_grid
		
		subroutine init_unstructured_grid_from_file(grid, gridfile)
			implicit none
			
			! Declare the intent of the variables.
			type(UnstructuredGrid), intent(inout) :: grid
			character(:), intent(in) :: gridfile
			
			! Looping variables for going through nodes, edges, and elements.
			integer :: i, j
			
			! Open the file and read in the number of nodes.
			open(unit = 1, gridfile, form = 'formatted')
			
			! Read in the number of nodes, 
			read(unit = 1, 'I5') grid%num_nodes
			
			! Allocate space for all of the nodes,
			! then read them in from file.
			allocate(grid%nodes(1:grid%num_nodes))
			read(unit = 1, '(F15.7, F15.7)') grid%nodes
			
			close(unit = 1)
		end subroutine init_unstructured_grid_from_file
		
		
end module grids






