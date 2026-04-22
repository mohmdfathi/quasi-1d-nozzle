program nozzle_1D
    use params 
    use routines
    implicit none

    ! ---- control variables
    integer           :: iter       ! iteration counter
    integer           :: i          ! grid index
    real              :: dt         ! time step (CFL)

    ! ---- primitive variables
    real, allocatable :: p(:)       ! pressure
    real, allocatable :: rho(:)     ! density
    real, allocatable :: u(:)       ! velocity

    ! ---- solution vectors and fluxes
    real, allocatable :: Q(:,:)     ! conservative variables
    real, allocatable :: Qn(:,:)    ! updated conservative variables
    real, allocatable :: F(:,:)     ! flux vector
    real, allocatable :: H(:,:)     ! source term
    real, allocatable :: dQdt(:,:)  ! time derivative

    
    ! ---- grid generation
    call gridGeneration()

    ! ---- primitive variables initialization 
    allocate(p(0:Imax), rho(0:Imax), u(0:Imax))
    associate( fac => 1.0 + 0.5*(k-1.0)*0.5**2 )
        p   = P0 / fac**(k/(k-1.0))       ! pressure
        rho = p / (R * T0 / fac)          ! density
        u   = 0.5 * sqrt( k * p / rho )   ! velocity (M = 0.5)
    end associate

    ! ---- allocate solver arrays
    allocate(Q(3,0:Imax), Qn(3,0:Imax),F(3,0:Imax), H(3,0:Imax), dQdt(3,0:Imax))

    ! ---- time-marching loop (MacCormack scheme)
    do iter = 1, iter_max

        ! ---- CFL time step
        dt = CFL * dx / maxval( abs(u) + sqrt(k * p / rho) )

        ! ---- build conservative forms
        call primitivesToConservatives(rho, u, p, Q)
        call primitivesToFluxesSources(rho, u, p, F, H)

        ! ---- predictor step (backward difference)
        do i = 1, Imax-1
            dQdt(:,i) = (F(:,i) - F(:,i-1)) / dx + H(:,i)
            Qn(:,i) = Q(:,i) - dt * dQdt(:,i) 
        end do
        call addArtificialViscosity(p, Q, Qn)

        call conservativesToPrimitives(Qn, rho, u, p)
        call primitivesBoundaries(rho, u, p)
        call primitivesToFluxesSources(rho, u, p, F, H)
        
        ! ---- corrector step (forward difference)
        do i = 1, Imax-1
            dQdt(:,i) = 0.5 * (dQdt(:,i) + (F(:,i+1) - F(:,i)) / dx + H(:,i))
            Qn(:,i) = Q(:,i) - dt * dQdt(:,i)
        end do
        call addArtificialViscosity(p, Q, Qn)

        ! ---- Convergence monitor: 
        if( mod(iter, 1000) == 0 )then
            block; real :: mdot_ref, imbalance

            mdot_ref = P0 * A(0.0) * sqrt( k / (R * T0) * (2.0 / (k + 1.0))**((k + 1.0) / (k - 1.0)) )
            imbalance = ( Qn(1,Imax-1) - Qn(1,1) ) / mdot_ref

            write(*,'(A,I10,3X,A,ES12.4)') "iter =", iter, "mass_imbalance =", imbalance
            end block
        end if

        call conservativesToPrimitives(Qn, rho, u, p)
        call primitivesBoundaries(rho, u, p)

    end do

    ! -------- output to CSV --------
    ! TODO..

    !  ---- cleanup
    deallocate(xi, Ai, dAdxi, p, rho, u, Q, Qn, F, H, dQdt)
    write(*,*) "Simulation finished."

end program nozzle_1D