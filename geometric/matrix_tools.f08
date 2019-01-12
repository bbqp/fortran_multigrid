module matrix_tools
    use, intrinsic iso_fortran_env
    use banded_matrices
    use structured_grids
    
    implicit none
    
    contains
        
        ! This function assembles the corresponding laplacian for a given 2D grid.
        ! When assembling the laplacian, we assume the corresponding vector of
        ! unknowns is laid out so that the rows corresponding to varying x-values
        ! (fixed y-values) are laid out end-to-end, starting from the smallest x-value.
        subroutine assemble_laplacian_2d(grid, A)
            type(Grid2D), intent(in) :: grid
            type(banded_matrix), intent(out) :: A
            
            integer, parameter :: M = grid%nx
            integer, parameter :: N = grid%ny
            real(R64), parameter :: dx = grid%hx
            real(R64), parameter :: dy = grid%hy
            
            integer, dimension(1:5) :: band_nums = (/ -M, -1, 0, 1, M /)
            integer :: i
            
            ! Initialize the band matrix.
            call bmat_init(A, .true., .true., M*N, M*N, 2, 2, band_nums)
            
            ! Set the individual bands of the matrix.
            call bmat_set_band(A, -2, -1 / dy**2)
            call bmat_set_band(A, -1, -1 / dx**2)
            call bmat_set_band(A, 0, 2.0 * ((1 / dx**2) + (1 / dy**2)))
            call bmat_set_band(A, 1, -1 / dx**2)
            call bmat_set_band(A, 2, -1 / dy**2)
        end subroutine assemble_laplacian_2d
        
        subroutine reshape_unknowns2d(grid, U, V)
            real(R64), dimension(:,:), intent(in) :: U
            real(R64), dimension(:), intent(inout) :: V
            
            V = reshape(U, SIZE(U(2:SIZE(U, 1), 2:SIZE(U, 2)), 1) * SIZE(U(2:SIZE(U, 1), 2:SIZE(U, 2)), 2))
        end subroutine reshape_unknowns2d+
        
end module matrix_tools
