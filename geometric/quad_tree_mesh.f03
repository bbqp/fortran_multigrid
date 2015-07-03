module quad_tree_mesh
	use precision_specs
	implicit none
	
	type, private :: qtree_node
		real(R64) :: x = 0.0
		real(R64) :: y = 0.0
		integer :: num_children = 0
		
		type(qtree_node), pointer :: child1
		type(qtree_node), pointer :: child2
		type(qtree_node), pointer :: child3
		type(qtree_node), pointer :: child4
	end type qtree_node

	type :: qtree_mesh
		type(qtree_node), private :: root
		integer :: tree_size = 0
		
		contains
		
		procedure :: add => qtree_add_point
		procedure :: get => qtree_get_point
		procedure :: clear => qtree_clear
	end type qtree_mesh

	interface
		module procedure qtree_add_point
		module procedure qtree_get_point
		module procedure qtree_clear
	end interface
	
	contains
	
	subroutine qtree_add_point(qtree, xcoord, ycoord)
		type(qtree_mesh), intent(inout) :: qtree
		real(R64), intent(in) :: xcoord
		real(R64), intent(in) :: ycoord
		
		
	end subroutine qtree_add_point

end module

















