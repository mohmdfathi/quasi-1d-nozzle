program nozzle_1D
    use parameters
    implicit none

    integer :: i
    
    real, allocatable :: xi(:)    ! x-location of grid points 
    real, allocatable :: Ai(:)    ! nozzle cross-sectional area at each grid point
    real, allocatable :: dAdxi(:) ! derivative of area dA/dx at each grid point


    ! allocate arrays and generate 1D grid
    allocate(xi(0:Imax), Ai(0:Imax), dAdxi(0:Imax))
    do i = 0, Imax
        xi(i) = x0 + i * dx
    end do
    write(*,*) "Grid generated."
    
 
    ! compute geometry using elemental functions
    Ai    = A(xi)
    dAdxi = dAdx(xi)
    write(*,*) "Nozzle area and gradients computed."

    
    ! cleanup
    deallocate(xi, Ai, dAdxi)
    write(*,*) "Done. Memory deallocated."

end program nozzle_1D