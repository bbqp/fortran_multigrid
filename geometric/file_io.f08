module file_io
    implicit none
    use, intrinsic iso_fortran_env
    use structured_grids

    interface write_to_file
        module procedure write_to_file2d
        module procedure write_to_file3d
    end interface write_to_file

    contains
    
    subroutine write_to_file2d(grid, U, filename)
        type(Grid2D), intent(in) :: grid
        type(REAL64), intent(in), dimension(1:grid%nx, 1:grid%ny) :: U
        character(len = *), intent(in) :: filename
        
        integer, parameter :: file_unit = 1
        integer :: error_code
        
        ! Open the file for writing in stream mode.
        open(unit = file_unit, file = filename, status = 'replace', access = 'stream', error = error_code)
        
        ! Write the size of the array in the first twelve bytes.
        !
        ! TODO: May need to worry about architecture/compiler implementation of integer types.
        write(file_unit) grid%nx
        write(file_unit) grid%ny
        
        ! Write the array itself.
        write(file_unit) U
        
        ! Close the file.
        close(file_unit)
    end subroutine
    
    subroutine write_to_file3d(grid, U, filename)
        type(Grid2D), intent(in) :: grid
        type(REAL64), intent(in), dimension(1:grid%nx, 1:grid%ny, 1:grid%nz) :: U
        character(len = *), intent(in) :: filename
        
        integer, parameter :: file_unit = 1
        integer :: error_code
        
        ! Open the file for writing in stream mode.
        open(unit = file_unit, file = filename, status = 'replace', error = error_code)
        
        ! Write the size of the array in the first twelve bytes.
        !
        ! TODO: May need to worry about architecture/compiler implementation of integer types.
        write(file_unit) grid%nx
        write(file_unit) grid%ny
        write(file_unit) grid%nz
        
        ! Write the array itself.
        write(file_unit) U
        
        ! Close the file.
        close(file_unit)
    end subroutine
    
end module
