module quad_tree_mesh
	use precision_specs
	implicit none
	
	type :: qtree_node
		real(dp), dimension(1:2, 1:2) :: bounds
		real(dp) :: 
		integer :: num_children = 0
		
		type(qtree_node), pointer :: north
		type(qtree_node), pointer :: south
		type(qtree_node), pointer :: east
		type(qtree_node), pointer :: west
	end type qtree_node

	type :: qtree_mesh
		type(qtree_node) :: root
		integer :: tree_size = 0
	end type qtree_mesh
	
	contains
	
	subroutine qtree_refine_global(qtree, xcoord, ycoord)
		type(qtree_mesh), intent(inout) :: qtree
		real(dp), intent(in) :: xcoord
		real(dp), intent(in) :: ycoord
		
		
	end subroutine qtree_refine_global
	
	subroutine qtree_refine_local(qtree, xcoord, ycoord)
		type(qtree_mesh), intent(inout) :: qtree
		real(dp), intent(in) :: xcoord
		real(dp), intent(in) :: ycoord
		
		
	end subroutine qtree_refine_local
	
end module

















