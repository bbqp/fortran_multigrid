! point_smoothers.f08
!
! Point smoothers for geometric multigrid.

module point_smoothers
    use, intrinsic iso_fortran_env
    use structured_grids
    implicit none
    
    ! Interface for the weighted jacobi method.
    interface wjacobi
        module procedure wjacobi_2d
        module procedure wjacobi_3d
    end interface wjacobi
    
    ! Interface for the guass seidel method.
    interface gauss_seidel
        module procedure gauss_seidel_2d
        module procedure gauss_seidel_3d
    end interface gauss_seidel
    
    ! Interface for the (weighted) SOR method.
    interface wsor
        module procedure wsor_2d
        module procedure wsor_3d
    end interface wsor
    
    contains
    
    !
    ! Subroutine definitions for weighted jacobi in 2D and 3D
    !
    
    subroutine wjacobi_2d(grid, U, F, w)
        type(Grid2D), intent(in) :: grid                                ! The grid that contains step size information, etc.
        real(REAL64), dimension(1:grid%nx, 1:grid%ny), intent(inout) :: U  ! The (discrete) function we are relaxing.
        real(REAL64), dimension(1:grid%nx, 1:grid%ny), intent(in) :: F     ! The right-hand-side of the PDE, in discrete form.
        real(REAL64), intent(in) :: w                                      ! The weight for the method.
        
        ! A temporary array that holds the next guess.
        real(REAL64), dimension(1:grid%nx, 1:grid%ny) :: temp
        
        ! Counters for the loops.
        integer :: i, j, k
        
        ! Temporary x- and y-coordinates.
        real(REAL64) :: x, y
        
        ! A temporary variable that stores the ratio of the square of the step sizes.
        real(REAL64) :: alpha
        
        ! Store the ratio of the step sizes.
        alpha = (grid%hx / grid%hy)**2
        
        ! First compute the interior points.
        !$omp parallel do private(i, j) collapse(2)
        do j = 2, grid%ny - 1
            do i = 2, grid%nx - 1
                temp(i, j) =    w * &
                                ( &
                                    F(i, j) + &
                                    U(i + 1, j) + U(i - 1, j) + &
                                    alpha * ( U(i, j + 1) + U(i, j - 1) )
                                ) / (2 * (1 + alpha))
                                + &
                                (1 - w) * U(i, j)
            end do
        end do
        !$omp end parallel do

        ! Set the new version of U to the current value of the temporary array.
        U = temp
    end subroutine wjacobi_2d
    
    subroutine wjacobi_3d(grid, U, F, w)
        type(Grid3D), intent(in) :: grid
        real(REAL64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz), intent(inout) :: U
        real(REAL64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz), intent(in) :: F
        real(REAL64), intent(in) :: w
        
        ! A temporary array that holds the next guess.
        real(REAL64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz) :: temp
        
        ! Counters for the loops.
        integer :: i, j, k, m
        
        ! Temporary x- and y-coordinates.
        real(REAL64) :: x, y, z
        
        ! A temporary variable that stores the ratio of the square of the step sizes.
        real(REAL64) :: alpha, beta
        
        ! Store the ratio of the step sizes.
        alpha = (grid%hx / grid%hy)**2
        beta = (grid%hx / grid%hz)**2
        
        ! First compute the interior points.
        !$omp parallel do private(i, j, k) collapse(3)
        do k = 2, grid%nz - 1
            do j = 2, grid%ny - 1
                do i = 2, grid%nx - 1
                    temp(i, j, k) = w * &
                                    ( &
                                        F(i, j, k) + &
                                        U(i + 1, j, k) + U(i - 1, j, k) + &
                                        alpha * ( U(i, j + 1, k) + U(i, j - 1, k) ) + &
                                        beta * (U(i, j, k + 1) + U(i, j, k - 1) ) &
                                    ) / (2 * (1 + alpha)) &
                                    + &
                                    (1 - w) * U(i, j, k)
                end do
            end do
        end do
        !$omp end parallel do
       
        ! Set the new version of U to the current value of the temporary array.
        U = temp
    end subroutine wjacobi_3d
    
    !
    ! Subroutine definitions for gauss seidel in 2D and 3D
    !
    
    subroutine gauss_seidel_2d(grid, U, F)
        type(Grid2D), intent(in) :: grid
        real(REAL64), dimension(1:grid%nx, 1:grid%ny), intent(inout) :: U
        real(REAL64), dimension(1:grid%nx, 1:grid%ny), intent(in) :: F
        
        ! A parameter to hold the ratio of our step sizes.
        real(REAL64), parameter :: alpha = (grid%hx / grid%hy)**2
        
        ! Counter variables for looping over the interior points and the boundaries.
        integer :: i, j, k
        
        ! Temporary x- and y-coordinates.
        real(REAL64) :: x, y
        
        ! Compute the new interior points.
        !$omp parallel do private(i, j) collapse(2)
        do j = 2, grid%ny - 1
            do i = 2, grid%nx - 1
                U(i, j) =   ( &
                                F(i, j) + &
                                U(i - 1, j) + U(i + 1, j) + &
                                alpha * ( &
                                    U(i, j - 1) + U(i, j + 1) &
                                ) &
                            ) / (2.0 * (1 + alpha))
            end do
        end do
        !$omp end parallel do
    end subroutine gauss_seidel_2d
    
    subroutine gauss_seidel_3d(grid, U)
        type(Grid3D), intent(in) :: grid
        real(REAL64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz), intent(inout) :: U
        
        real, parameter :: alpha = (grid%hx / grid%hy)**2
        real, parameter :: beta = (grid%hx / grid%hz)**2
        
        integer :: i, j, k
        
        ! First relax the interior points.
        !$omp parallel do private(i, j, k) collapse(3)
        do k = 2, grid%nz - 1
            do j = 2, grid%ny - 1
                do i = 2, grid%nx - 1
                    U(i, j, k) =    ( &
                                        F(i, j, k) + &
                                        
                                        U(i - 1, j, k) + U(i + 1, j, k) + &
                                        
                                        alpha * ( &
                                            U(i, j - 1, k) + U(i, j + 1, k)
                                        ) + &
                                        
                                        beta * ( &
                                            U(i, j, k - 1) + U(i, j, k + 1)
                                        )
                                    ) / (2.0 * (1 + alpha + beta))
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine gauss_seidel_3d
    
    !
    ! Subroutine definitions for (weighted) successive overrelaxation in 2D and 3D.
    !
    
    subroutine wsor_2d(grid, U, F)
        type(Grid2D), intent(in) :: grid
        real(REAL64), dimension(1:grid%nx, 1:grid%ny), intent(inout) :: U
        real(REAL64), dimension(1:grid%nx, 1:grid%ny), intent(in) :: F
        
        real(REAL64), parameter :: alpha = (grid%hx / grid%hy)**2
        
        integer :: i, j, k
        
        ! First relax the interior points.
        !$omp parallel do private(i, j) collapse(2)
        do j = 2, grid%ny - 1
            do i = 2, grid%nx - 1
                U(i, j) =   w * ( &
                                    F(i, j) + &
                                    
                                    U(i - 1, j) + U(i + 1, j) + &
                                    
                                    alpha * ( &
                                        U(i, j - 1) + U(i, j + 1) &
                                    ) &
                                ) / (2.0 * (1 + alpha)) + &
                                
                            (1 - w) * U(i, j)
            end do
        end do
        !$omp end parallel do
    end subroutine wsor_2d
    
    subroutine wsor_3d(grid, U, F)
        type(Grid3D), intent(in) :: grid
        real(REAL64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz), intent(inout) :: U
        real(REAL64), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz), intent(in) :: F
        
        real(REAL64), parameter :: alpha = (grid%hx / grid%hy)**2
        real(REAL64), parameter :: beta = (grid%hx / grid%hz)**2
        
        integer :: i, j, k, m
        
        ! First relax the interior points.
        !$omp parallel do private(i, j, k) collapse(3)
        do k = 2, grid%nz - 1
            do j = 2, grid%ny - 1
                do i = 2, grid%nx - 1
                    U(i, j, k) =    w * ( &
                                        F(i, j, k) + &
                                        
                                        U(i - 1, j, k) + U(i + 1, j, k) + &
                                        
                                        alpha * ( &
                                            U(i, j - 1, k) + U(i, j + 1, k) &
                                        ) + &
                                        
                                        beta * ( &
                                            U(i, j, k - 1) + U(i, j, k + 1) &
                                        )
                                    ) / (2.0 * (1 + alpha + beta)) + &
                                    
                                    (1 - w) * U(i, j, k)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine wsor_3d

end module point_smoothers























