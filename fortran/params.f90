module parameters
implicit none

    ! grid & solver settings
    integer, parameter :: Imax = 7500          ! Grid size (number of cells in x-direction)
    integer, parameter :: iter_max = 100000    ! Maximum number of iterations for time marching

    ! physical constants
    real, parameter    :: R   = 286.9865       ! Gas constant [J/(kg·K)]
    real, parameter    :: k   = 1.4            ! Specific heat ratio of air
    real, parameter    :: cv  = 718.0          ! Specific heat at constant volume [J/(kg·K)]
    real, parameter    :: kR  = k * R          ! Derived constant (k * R)
    real, parameter    :: P0  = 25000.0        ! Inlet total pressure [Pa]
    real, parameter    :: T0  = 300.0          ! Inlet total temperature [K]
    real, parameter    :: Pb  = 15000.0        ! Outlet back pressure [Pa]

    ! coefficients 
    real, parameter    :: CFL = 0.4            ! CFL number for time-step stability
    real, parameter    :: Cx  = 0.3            ! Numerical dissipation coefficient 
 
    ! geometry info
    real, parameter    :: x0  = -0.25           ! Domain start location [m]
    real, parameter    :: x1  =  0.50           ! Domain end location [m]
    real, parameter    :: dx  = ( x1 - x0 ) / Imax ! Grid spacing

contains

    ! ----------------------------
    ! nozzle area A(x) 
    elemental real function A(x)
        real, intent(in) :: x

        if( x <= -0.20 )then
            A = 0.03786

        elseif( x <= 0.0 )then
            A = 0.03150 + (0.03786 - 0.03150) / (-0.20) * x

        else
            A = 0.03150 + (0.05700 - 0.03150) / (0.50) * x
        endif

    end function


    ! ----------------------------
    ! derivative dA/dx 
    elemental real function dAdx(x)
        real, intent(in) :: x

        if (x <= -0.20) then
            dAdx = 0.0

        elseif (x <= 0.0) then
            dAdx = (0.03786 - 0.03150) / (-0.20)

        else
            dAdx = (0.05700 - 0.03150) / (0.50)
        endif

    end function

end module parameters