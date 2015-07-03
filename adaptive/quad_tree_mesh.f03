module quad_tree_mesh
	use precision_specs
	implicit none
	
	type :: qtree_node
		real(dp), dimension(1:2, 1:2) :: bounds
		real(dp), dimension(1:2) :: centroid
		real(dp) :: hx, hy
		
		! Neighbors for the given node.
		type(qtree_node), dimension(:), allocatable :: neighbors
		
		! Children for the given node.
		type(qtree_node), dimension(:), allocatable :: children
	end type qtree_node

	type :: qtree_mesh
		type(qtree_node) :: root
		integer :: tree_size = 0
		
		integer, dimension(1:4) :: boundary_flags
		procedure, pointer :: eps
		procedure, pointer :: f
	end type qtree_mesh
	
	contains
	
	subroutine qtree_refine_global(qtree, xcoord, ycoord)
		type(qtree_mesh), intent(inout) :: qtree
		real(dp), intent(in) :: xcoord
		real(dp), intent(in) :: ycoord
		
		type(qtree_node), pointer :: next, prev
		
		if(associated(qtree%root)) then
			prev => null()
			next => qtree%root
			
			while(associated(next))
				
				! Once we've hit a leaf, refine it, then move up.
				if(qtree_is_leaf(next)) then
					
					call qtnode_refine(next)
					
					! Set the pointers to move back up the tree.
					next%
				else
				
				end if
				
			end while
			
			
		end if
	end subroutine qtree_refine_global
	
	subroutine qtree_refine_local(qtree, xcoord, ycoord)
		type(qtree_mesh), intent(inout) :: qtree
		real(dp), intent(in) :: xcoord
		real(dp), intent(in) :: ycoord
		
		type(qtree_node), pointer :: next, prev
		
		if(associated(qtree%root)) then
			prev => null()
			next => qtree%root
			
			while(associated(next))
				
				if(qtree_is_leaf(next)) then
					
					
				else
				
				end if
				
			end while
			
			
		end if
	end subroutine qtree_refine_local
	
end module

















