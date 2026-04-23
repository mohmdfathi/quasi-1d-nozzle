module routines
implicit none

! ==============================================================
!  CONTENTS
! ==============================================================
! 1. primitivesToConservatives   ::  Convert primitive variables (ρ, u, p) → conservative form Q
! 2. conservativesToPrimitives   ::  Convert conservative variables Q → primitive variables
! 3. primitivesBoundaries        ::  Apply inlet/outlet boundary conditions
! 4. primitivesToFluxesSources   ::  Compute flux vector F and geometric source term H
! 5. addArtificialViscosity      ::  Jameson-type artificial dissipation for shock capturing
! ==============================================================

contains

! ==============================================================
! primitives → conservatives
! ==============================================================
subroutine primitivesToConservatives(is, ie, rho, u, p, Q)
    use params,     only : Imax, Ai, k
    implicit none
    integer, intent(in) :: is, ie
    real, intent(in)    :: rho(0:Imax), u(0:Imax), p(0:Imax)
    real, intent(out)   :: Q(3,0:Imax)

    integer            :: i
    real               :: E 

    do concurrent( i = is : ie )

    E = p(i) / ((k - 1.0) * rho(i)) + 0.5 * u(i)**2

    Q(1,i) = rho(i) * Ai(i)
    Q(2,i) = rho(i) * u(i) * Ai(i)
    Q(3,i) = rho(i) * E * Ai(i)
    end do 

end subroutine


! ==============================================================
! conservatives → primitives
! ==============================================================
subroutine conservativesToPrimitives(is, ie, Q, rho, u, p)
    use params,     only : Imax, Ai, k
    implicit none
    integer, intent(in) :: is, ie
    real, intent(in)    :: Q(3,0:Imax)
    real, intent(out)   :: rho(0:Imax), u(0:Imax), p(0:Imax)

    integer           :: i
    real              :: E 

    do concurrent( i = is : ie )

        rho(i) = Q(1,i) / Ai(i)
        u(i)   = Q(2,i) / Q(1,i)
        E      = Q(3,i) / (rho(i) * Ai(i))

        p(i) = (k - 1.0) * rho(i) * (E - 0.5 * u(i)**2)
    end do 


end subroutine


! ==============================================================
! boundary conditions
! ==============================================================
subroutine primitivesBoundaries(is, ie, rho, u, p)
    use params,     only : Imax, T0, P0, Pb, R, kR, k
    implicit none
    integer, intent(in) :: is, ie
    real, intent(inout) :: rho(0:Imax), u(0:Imax), p(0:Imax)

    real                :: M2, T

    ! inlet (subsonic)
    if( is == 0 )then
        p(0) = p(1)

        M2 = ((P0/p(0))**((k-1)/k) - 1.0) * (2.0/(k-1.0))
        T  = T0 /( 1.0 + 0.5 * ( k - 1.0 ) * M2 )

        rho(0) = p(0) / (R*T)
        u(0)   = sqrt(M2 * kR * T)
    end if

    ! outlet (subsonic)
    if( ie == Imax )then
        p(Imax)   = Pb
        rho(Imax) = 2. * rho(Imax-1) - rho(Imax-2)
        u(Imax)   = 2. * u(Imax-1)   - u(Imax-2)
    end if
end subroutine


! ==============================================================
! fluxes + source term
! ==============================================================
subroutine primitivesToFluxesSources(is, ie, rho, u, p, F, H)
    use params,     only : Imax, Ai, dAdxi, k
    implicit none
    integer, intent(in) :: is, ie
    real, intent(in)    :: rho(0:Imax), u(0:Imax), p(0:Imax) 
    real, intent(out)   :: F(3,0:Imax), H(3,0:Imax)

    integer             :: i
    real                :: E 

    do concurrent( i = is : ie )
        E = p(i) / ((k - 1.0) * rho(i)) + 0.5 * u(i)**2

        F(1,i) = rho(i) * u(i) * Ai(i)
        F(2,i) = (rho(i) * u(i)**2 + p(i)) * Ai(i)
        F(3,i) = u(i) * (rho(i) * E + p(i)) * Ai(i)

        H(1,i) = 0.0
        H(2,i) = -p(i) * dAdxi(i)
        H(3,i) = 0.0
    end do 


end subroutine


! ==============================================================
! artificial viscosity (Jameson style)
! ==============================================================
subroutine addArtificialViscosity(is, ie, p, Q, Qn)
    use params,     only : Imax, Cx
    implicit none
    integer, intent(in) :: is, ie
    real, intent(in)    :: p(0:Imax)
    real, intent(in)    :: Q(3,0:Imax)
    real, intent(inout) :: Qn(3,0:Imax) 

    integer             :: i
    real                :: nu
    real                :: D2(3)

    do concurrent( i = max(is,1) : min(ie, Imax-1) )

        nu = abs(p(i+1) - 2.0*p(i) + p(i-1)) / &
                (p(i+1) + 2.0*p(i) + p(i-1))

        D2 = Q(:,i+1) - 2.0*Q(:,i) + Q(:,i-1)

        Qn(:,i) = Qn(:,i) + Cx * nu * D2
    end do 

end subroutine

end module