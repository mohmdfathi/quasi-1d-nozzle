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
subroutine primitivesToConservatives(rho, u, p, Q)
    use params,   only : Imax, Ai, k
    implicit none

    real, intent(in)  :: rho(0:Imax), u(0:Imax), p(0:Imax)
    real, intent(out) :: Q(3,0:Imax)

    real :: E(0:Imax)

    E = p / ((k - 1.0) * rho) + 0.5 * u**2

    Q(1,:) = rho * Ai
    Q(2,:) = rho * u * Ai
    Q(3,:) = rho * E * Ai
end subroutine


! ==============================================================
! conservatives → primitives
! ==============================================================
subroutine conservativesToPrimitives(Q, rho, u, p)
    use params,   only : Imax, Ai, k
    implicit none

    real, intent(in)  :: Q(3,0:Imax)
    real, intent(out) :: rho(0:Imax), u(0:Imax), p(0:Imax)

    real :: E(0:Imax)

    rho = Q(1,:) / Ai
    u   = Q(2,:) / Q(1,:)
    E   = Q(3,:) / (rho * Ai)

    p = (k - 1.0) * rho * (E - 0.5 * u**2)
end subroutine


! ==============================================================
! boundary conditions
! ==============================================================
subroutine primitivesBoundaries(rho, u, p)
    use params,     only : Imax, T0, P0, Pb, R, kR, k
    implicit none

    real, intent(inout) :: rho(0:Imax), u(0:Imax), p(0:Imax)

    real :: M2, T

    ! inlet (subsonic)
    p(0) = p(1)

    M2 = ((P0/p(0))**((k-1)/k) - 1.0) * (2.0/(k-1.0))
    T  = T0 / (1.0 + 0.5*(k-1.0)*M2)

    rho(0) = p(0) / (R*T)
    u(0)   = sqrt(M2 * kR * T)

    ! outlet (extrapolation)
    p(Imax)   = Pb
    rho(Imax) = 2*rho(Imax-1) - rho(Imax-2)
    u(Imax)   = 2*u(Imax-1)   - u(Imax-2)
end subroutine


! ==============================================================
! fluxes + source term
! ==============================================================
subroutine primitivesToFluxesSources(rho, u, p, F, H)
    use params,   only : Imax, Ai, dAdxi, k
    implicit none

    real, intent(in)  :: rho(0:Imax), u(0:Imax), p(0:Imax) 
    real, intent(out) :: F(3,0:Imax), H(3,0:Imax)

    real :: E(0:Imax)

    E = p / ((k - 1.0) * rho) + 0.5 * u**2

    F(1,:) = rho * u * Ai
    F(2,:) = (rho * u**2 + p) * Ai
    F(3,:) = u * (rho * E + p) * Ai

    H(1,:) = 0.0
    H(2,:) = -p * dAdxi
    H(3,:) = 0.0
end subroutine


! ==============================================================
! artificial viscosity (Jameson style)
! ==============================================================
subroutine addArtificialViscosity(p, Q, Qn)
    use params,     only : Imax, Cx
    implicit none

    real, intent(in)    :: p(0:Imax)
    real, intent(in)    :: Q(3,0:Imax)
    real, intent(inout) :: Qn(3,0:Imax) 

    integer :: i
    real :: nu
    real :: D2(3)

    do i = 2, Imax-2

        nu = abs(p(i+1) - 2.0*p(i) + p(i-1)) / &
             (p(i+1) + 2.0*p(i) + p(i-1))

        D2 = Q(:,i+1) - 2.0*Q(:,i) + Q(:,i-1)

        Qn(:,i) = Qn(:,i) + Cx * nu * D2
    end do

end subroutine

end module