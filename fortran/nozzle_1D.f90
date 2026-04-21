program nozzle_1D
    use parameters
    implicit none

    integer :: i
    
    real, allocatable :: xi(:)    ! x-location of grid points 
    real, allocatable :: Ai(:)    ! nozzle cross-sectional area at each grid point
    real, allocatable :: dAdxi(:) ! derivative of area dA/dx at each grid point

    real, allocatable :: T(:)     ! static temperature
    real, allocatable :: p(:)     ! static pressure
    real, allocatable :: rho(:)   ! density
    real, allocatable :: u(:)     ! axial velocity


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

    
    ! allocate primitive variables
    allocate(T(0:Imax), p(0:Imax), rho(0:Imax), u(0:Imax))

    ! initialize primitive variables
    ! based on isentropic relations (M = 0.5)
    associate( fac => 1.0 + 0.5*(k-1.0)*0.5**2 )
    T   = T0 / fac             ! static temperature from isentropic relations
    P   = P0 / fac**(k/(k-1.)) ! static pressure from isentropic relations
    rho = P / ( R * T )        ! density from ideal gas law
    u   = 0.5 * sqrt( kR * T ) ! velocity magnitude (from Mach definition)
    end associate
    write(*,*) "Primitive variables initialized."
 

    
    ! cleanup
    deallocate(xi, Ai, dAdxi, T, p, rho, u)
    write(*,*) "Done. Memory deallocated."

end program nozzle_1D