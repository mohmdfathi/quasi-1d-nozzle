program nozzle_1D
    use params 
    use omp_lib
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

    ! ---- wall time 
    integer           :: count_start, count_end, count_rate
    real              :: elapsed_time

    ! ---- threads info
    integer           :: is, ie
    integer           :: tid, nthreads

    
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

    call system_clock(count_start, count_rate)

    ! ---- time-marching loop (MacCormack scheme)
    do iter = 1, iter_max

        ! ---- CFL time step
        dt = CFL * dx / maxval( abs(u) + sqrt(k * p / rho) )

        !$omp parallel default(shared) private(is,ie,tid)
        tid = omp_get_thread_num()
        nthreads = omp_get_num_threads()

        call getThreadBounds(tid, nthreads, is, ie)
       
        ! ---- build conservative forms
        call primitivesToConservatives(is, ie, rho, u, p, Q)
        call primitivesToFluxesSources(is, ie, rho, u, p, F, H)
        
        !$omp barrier
        ! ---- predictor step (backward difference)
        do i = max(1,is), min(Imax-1,ie)
            dQdt(:,i) = (F(:,i) - F(:,i-1)) / dx + H(:,i)
            Qn(:,i) = Q(:,i) - dt * dQdt(:,i) 
        end do
        call addArtificialViscosity(is, ie, p, Q, Qn)
        call conservativesToPrimitives(is, ie, Qn, rho, u, p)
        call primitivesBoundaries(is, ie, rho, u, p)
        call primitivesToFluxesSources(is, ie, rho, u, p, F, H)
        
        !$omp barrier
        ! ---- corrector step (forward difference)
        do i = max(1,is), min(Imax-1,ie)
            dQdt(:,i) = 0.5 * (dQdt(:,i) + (F(:,i+1) - F(:,i)) / dx + H(:,i))
            Qn(:,i) = Q(:,i) - dt * dQdt(:,i)
        end do
        call addArtificialViscosity(is, ie, p, Q, Qn)
        call conservativesToPrimitives(is, ie, Qn, rho, u, p)
        call primitivesBoundaries(is, ie, rho, u, p)

        !$omp end parallel

        ! ---- Convergence monitor
        if( mod(iter, 1000) == 0 )then
            block; real :: mdot_ref, imbalance

            mdot_ref = P0 * A_throat * sqrt( k / (R * T0) * (2.0 / (k + 1.0))**((k + 1.0) / (k - 1.0)) )
            imbalance = ( Qn(1,Imax-1) - Qn(1,1) ) / mdot_ref

            write(*,'(A,I10,3X,A,ES12.4)') "iter =", iter, "mass_imbalance =", imbalance
            end block
        end if


    end do

    call system_clock(count_end)
    elapsed_time = real( count_end - count_start )/real( count_rate )
    write(*,'(A,I0,A)') "Total wall time = ", nint(elapsed_time), " s"

    ! -------- output to CSV --------
    if( nthreads == 1 )then 
        call writeResults("nozzle_maccormack_serial.csv")
    else
        call writeResults("nozzle_maccormack_parallel.csv")
    end if

    !  ---- cleanup
    deallocate(xi, Ai, dAdxi, p, rho, u, Q, Qn, F, H, dQdt)
    write(*,*) "Simulation finished."

contains

! ==============================================================
! write results to a CSV file
! ==============================================================
subroutine writeResults(filename)
    implicit none
    character(len=*), intent(in) :: filename

    integer    :: i
    real       :: rho_ref, a0, mdot_ref

    open(unit=10, file=filename, status="replace")

    write(10,'(A)') "x/L,p/p0,Mach,mdot/mdot_ref,A"

    ! ---- reference values ----
    rho_ref = P0 / (R * T0)
    a0      = sqrt(k * R * T0)
    mdot_ref = rho_ref * a0 * A_throat

    do i = 0, Imax
        write(10,'(5E16.8)') &
            xi(i)/(x1 - x0), &                         ! x/L
            p(i)/P0, &                                 ! p/p0
            u(i)/sqrt(k*p(i)/rho(i)), &                ! Mach
            (rho(i)*u(i)*Ai(i)) / mdot_ref, &          ! mdot/mdot_ref
            Ai(i)                                      ! A(x)
    end do

    close(10)

end subroutine

! ==============================================================
! compute contiguous index range [is, ie] for each thread
! over a 1D domain (0 ... Imax), with load-balanced partitioning
! ==============================================================
subroutine getThreadBounds(tid, nthreads, is, ie)
    implicit none
    integer, intent(in )  :: tid, nthreads
    integer, intent(out)  :: is, ie 

    integer               :: n, chunk, rem
   
    ! total number of points (0 ... Imax)
    n = Imax + 1 

    ! uniform chunk size + leftover points
    chunk = n / nthreads
    rem   = mod(n, nthreads)

    ! distribute remainder among first threads
    if( tid < rem )then
        is = tid * ( chunk + 1 )
        ie = is + chunk
    else
        is = tid * chunk + rem
        ie = is + chunk - 1
    end if

end subroutine
    
end program nozzle_1D