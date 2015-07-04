module quad_tree_mesh
	use precision_specs
	implicit none
	
	type :: qtree_node
		real(dp), dimension(1:2, 1:2) :: bounds
		real(dp), dimension(1:2) :: centroid
		real(dp) :: h
		integer :: idx
		logical :: need_to_refine
				
		real(dp) :: previous_val
		real(dp) :: current_val
		
		! Pointer to the parent node.
		type(qtree_node), pointer :: parent
		
		! Neighbors for the given node.
		type(qtree_node), pointer :: north
		type(qtree_node), pointer :: south
		type(qtree_node), pointer :: east
		type(qtree_node), pointer :: west
		
		! Children for the given node.
		type(qtree_node), pointer :: nechild
		type(qtree_node), pointer :: nwchild
		type(qtree_node), pointer :: swchild
		type(qtree_node), pointer :: sechild
	end type qtree_node

	type :: qtree_mesh
		type(qtree_node), pointer :: root
		integer :: tree_size = 0
		
		integer, dimension(1:4) :: boundary_flags
		
		! Pointers for the coefficient and right-hand-side functions.
		procedure, pointer :: coefficient
		procedure, pointer :: right_hand_side
	end type qtree_mesh
	
	contains
	
	!--------------------------------------------------------------!
	! Subroutines for working with the individual quad tree nodes. !
	!--------------------------------------------------------------!
	
	subroutine qtree_node_init(node, parent, bnds, step_size, indx)
		type(qtree_node), intent(inout) :: node
		type(qtree_node), intent(in) :: parent
		real(dp), intent(in), dimension(1:2, 1:2) :: bnds
		real(dp), intent(in) :: step_size
		integer, intent(in) :: indx
		
		node%bounds = bnds
		node%centroid = 0.5 * sum(bnds, 2)
		node%h = step_size
		node%idx = indx
		
		node%previous_val = 0.0_dp
		node%current_val = 0.0_dp
		
		! Set the parent for the given node.
		node%parent => parent
		
		! Neighbors for the given node.
		node%north = NULL()
		node%south = NULL()
		node%east = NULL()
		node%west = NULL()
		
		! Children for the given node.
		node%nechild = NULL()
		node%nwchild = NULL()
		node%swchild = NULL()
		node%sechild = NULL()
	end subroutine qtree_node_init
	
	subroutine qtree_node_refine(node)
		type(qtree_node), intent(inout) :: node
		real(dp) :: bnds
		
		! Allocate and initialize all of the children.
		allocate(qtree_node%children)
		
		! Initialize the child in the northeast quadrant.
		bnds(1, :) = (/ node%bounds(1, 1) + node%h / 2.0, node%bounds(1, 2) /)
		bnds(2, :) = (/ node%bounds(2, 1) + node%h / 2.0, node%bounds(2, 2) /)
		qtree_node_init(node%nechild, node, bnds, node%h / 2.0, 4 * node%idx - 2)
		
		! Initialize the child in the northwest quadrant.
		bnds(1, :) = (/ node%bounds(1, 1), node%bounds(1, 1) + node%h / 2.0 /) 
		bnds(2, :) = (/ node%bounds(2, 1) + node%h / 2.0, node%bounds(2, 2) /) 
		qtree_node_init(node%nwchild, node, bnds, node%h / 2.0, 4 * node%idx - 1)
		
		! Initialize the child in the southwest quadrant.
		bnds(1, :) = (/ node%bounds(1, 1), node%bounds(1, 1) + node%h / 2.0 /) 
		bnds(2, :) = (/ node%bounds(2, 1), node%bounds(2, 1) + node%h / 2.0 /) 
		qtree_node_init(node%swchild, node, bnds, node%h / 2.0, 4 * node%idx)
		
		! Initialize the child in the southeast quadrant.
		bnds(1, :) = (/ node%bounds(1, 1) + node%h / 2.0, node%bounds(1, 2) /)
		bnds(2, :) = (/ node%bounds(2, 1) + node%h / 2.0, node%bounds(2, 2) /)  
		qtree_node_init(node%sechild, node, bnds, node%h / 2.0, 4 * node%idx + 1)
	end subroutine qtree_node_refine
	
	!
	! Subroutine that checks whether or not a node has children. Quad tree
	! nodes either have no children or all four, so we can aribitrarily check
	! whether or not once of the children is associated to satisfy this.
	!
	function qtree_node_is_leaf(node) result(is_leaf)
		type(qtree_node), intent(in) :: node
		logical :: is_leaf
		
		is_leaf = .not. associated(node%nechild)
	end function qtnode_is_leaf
	
	!---------------------------------------------------------!
	! Subroutines for working with the quad tree mesh proper. !
	!---------------------------------------------------------!
	
	subroutine qtree_mesh_refine_global(qtree)
		type(qtree_mesh), intent(inout) :: qtree
		type(qtree_node), pointer :: next, prev
		
		if(associated(qtree%root)) then
			prev => null()
			next => qtree%root
			
			while(associated(next))
				
				if(qtree_is_leaf(next)) then
					! Once we've hit a leaf, refine it, then move up.
					call qtree_node_refine(next)
					
					! Set the pointers to move back up the tree.
					prev => next
					next => next%parent
				else
					! In this case, we will either move to the next available child or up the tree.
					if(associated(prev, node%nechild)) then
						prev => next
						next => node%nwchild
					else if(associated(prev, node%nwchild)) then
						prev => next
						next => node%swchild
					else if(associated(prev, node%swchild)) then
						prev => next
						next => node%sechild
					else
						prev = next
						next => next%parent
					end if
				end if
				
			end while
		end if
	end subroutine qtree_refine_global
	
	subroutine qtree_refine_local(qtree)
		type(qtree_mesh), intent(inout) :: qtree
		type(qtree_node), pointer :: next, prev
		
		if(associated(qtree%root)) then
			prev => null()
			next => qtree%root
			
			while(associated(next))
				
				if(qtree_node_is_leaf(next) .and. next%need_to_refine) then
					! Once we've hit a leaf that is flagged for refinement, refine it, then move up.
					call qtree_node_refine(next)
					
					! Set the pointers to move back up the tree.
					prev => next
					next => next%parent
				else
					! In this case, we will either move to the next available child or up the tree.
					if(associated(prev, node%nechild)) then
						prev => next
						next => node%nwchild
					else if(associated(prev, node%nwchild)) then
						prev => next
						next => node%swchild
					else if(associated(prev, node%swchild)) then
						prev => next
						next => node%sechild
					else
						prev = next
						next => next%parent
					end if
				end if
				
			end while
		end if
	end subroutine qtree_refine_local
	
end module

















