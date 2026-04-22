module params
implicit none

! ==============================================================
!  CONTENTS
! ==============================================================
! 1. Grid & solver parameters    ::  Imax, iter_max, CFL, Cx, dx, domain limits
! 2. Physical constants          ::  Gas properties (R, k, cv), inlet/outlet conditions (P0, T0, Pb)
! 3. Geometry functions          ::  A(x) → nozzle cross-sectional area, dAdx(x) → spatial derivative of area
! 4. Grid arrays                 ::  xi → grid coordinates, Ai → area distribution on grid, dAdxi → area gradient on grid
! 5. gridGeneration              ::  Allocates and builds computational grid and geometry fields
! ==============================================================

! ---- solver settings
integer, parameter :: Imax = 12800         ! Grid size (number of cells in x-direction)
integer, parameter :: iter_max = 500000    ! Maximum number of iterations for time marching

! ---- physical constants
real, parameter    :: R   = 286.9865       ! Gas constant [J/(kg·K)]
real, parameter    :: k   = 1.4            ! Specific heat ratio of air
real, parameter    :: cv  = 718.0          ! Specific heat at constant volume [J/(kg·K)]
real, parameter    :: kR  = k * R          ! Derived constant (k * R)
real, parameter    :: P0  = 25000.0        ! Inlet total pressure [Pa]
real, parameter    :: T0  = 300.0          ! Inlet total temperature [K]
real, parameter    :: Pb  = 15000.0        ! Outlet back pressure [Pa]

! ---- coefficients 
real, parameter    :: CFL = 0.4            ! CFL number for time-step stability
real, parameter    :: Cx  = 5.0            ! Numerical dissipation coefficient 

! ---- geometry info
real, parameter    :: x0  = -0.25           ! Domain start location [m]
real, parameter    :: x1  =  1.03           ! Domain end location [m]
real, parameter    :: dx  = ( x1 - x0 ) / Imax ! Grid spacing

real, allocatable :: xi(:)                  ! grid coordinates
real, allocatable :: Ai(:)                  ! cross-sectional area A(x)
real, allocatable :: dAdxi(:)               ! area gradient dA/dx

contains

! ==============================================================
! nozzle area A(x) 
! ==============================================================
elemental real function A(x)
    implicit none
    real, intent(in) :: x

    if( x <= 0.0 )then
        A = 0.03150 - 0.0318 * x
    else
        A = 0.03150 + 0.0510 * x
    endif

end function

! ==============================================================
! derivative dA/dx(x)
! ==============================================================
elemental real function dAdx(x)
    implicit none
    real, intent(in) :: x

    if( x <= 0.0 )then
        dAdx = -0.0318
    else
        dAdx = 0.0510
    endif

end function

! ==============================================================
! grid generation subroutine
! ==============================================================
subroutine gridGeneration()
    implicit none 
    integer :: i 
    
    allocate( xi(0:Imax), Ai(0:Imax), dAdxi(0:Imax) )

    do i = 0, Imax
        xi(i) = x0 + i * dx
    end do
    Ai    = A(xi)      ! area distribution
    dAdxi = dAdx(xi)   ! area gradient
end subroutine

end module