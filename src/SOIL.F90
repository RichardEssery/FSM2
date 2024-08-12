!-----------------------------------------------------------------------
! Update soil temperatures
!-----------------------------------------------------------------------
subroutine SOIL(csoil,dt,Gsoil,ksoil,Tsoil)

use LAYERS, only: &
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsoil               ! Number of soil layers

implicit none

real, intent(in) :: &
  dt,                &! Timestep (s)
  Gsoil,             &! Heat flux into soil (W/m^2)
  csoil(Nsoil),      &! Areal heat capacity of soil layers (J/K/m^2)
  ksoil(Nsoil)        ! Thermal conductivity of soil layers (W/m/K)

real, intent(inout) :: &
  Tsoil(Nsoil)        ! Soil layer temperatures (K)

integer :: n          ! Soil layer counter

real :: &
  a(Nsoil),          &! Below-diagonal matrix elements
  b(Nsoil),          &! Diagonal matrix elements
  c(Nsoil),          &! Above-diagonal matrix elements
  dTs(Nsoil),        &! Temperature increments (K)
  U(Nsoil),          &! Thermal transmittance between layers (W/m^2/K)
  rhs(Nsoil)          ! Matrix equation rhs

U(1) = 2 / (Dzsoil(1)/ksoil(1) + Dzsoil(2)/ksoil(2))
a(1) = 0
b(1) = csoil(1) + U(1)*dt
c(1) = - U(1)*dt
rhs(1) = (Gsoil - U(1)*(Tsoil(1) - Tsoil(2)))*dt
do n = 2, Nsoil - 1
  U(n) = 2 / (Dzsoil(n)/ksoil(n) + Dzsoil(n+1)/ksoil(n+1))
  a(n) = c(n-1)
  b(n) = csoil(n) + (U(n-1) + U(n))*dt
  c(n) = - U(n)*dt
  rhs(n) = U(n-1)*(Tsoil(n-1) - Tsoil(n))*dt  &
           + U(n)*(Tsoil(n+1) - Tsoil(n))*dt 
end do
n = Nsoil
U(n) = ksoil(n)/Dzsoil(n)
a(n) = c(n-1)
b(n) = csoil(n) + (U(n-1) + U(n))*dt
c(n) = 0
rhs(n) = U(n-1)*(Tsoil(n-1) - Tsoil(n))*dt
call TRIDIAG(Nsoil,Nsoil,a,b,c,rhs,dTs)
do n = 1, Nsoil
  Tsoil(n) = Tsoil(n) + dTs(n)
end do

end subroutine SOIL
