module line_smoothers
	use banded_matrices
	use precision_specs
	use structured_grids
	implicit none
	
	interface assemble_laplacian
		module procedure assemble_laplacian1d
		module procedure assemble_laplacian2d
	end interface assemble_laplacian
	
	contains
	
	function assemble_laplacian1d(grid) result(A)
		type(Grid2D), intent(in) :: grid
		
		! The output of our method for line relaxation.
		type(banded_matrix), intent(out) :: A
		
		bmatrix(
		
		return
	end function assemble_laplacian2d
	
end module