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

integer :: k          ! Soil layer counter

real :: &
  a(Nsoil),          &! Below-diagonal matrix elements
  b(Nsoil),          &! Diagonal matrix elements
  c(Nsoil),          &! Above-diagonal matrix elements
  dTs(Nsoil),        &! Temperature increments (k)
  gs(Nsoil),         &! Thermal conductivity between layers (W/m^2/k)
  rhs(Nsoil)          ! Matrix equation rhs

do k = 1, Nsoil - 1
  gs(k) = 2 / (Dzsoil(k)/ksoil(k) + Dzsoil(k+1)/ksoil(k+1))
end do
a(1) = 0
b(1) = csoil(1) + gs(1)*dt
c(1) = - gs(1)*dt
rhs(1) = (Gsoil - gs(1)*(Tsoil(1) - Tsoil(2)))*dt
do k = 2, Nsoil - 1
  a(k) = c(k-1)
  b(k) = csoil(k) + (gs(k-1) + gs(k))*dt
  c(k) = - gs(k)*dt
  rhs(k) = gs(k-1)*(Tsoil(k-1) - Tsoil(k))*dt  &
           + gs(k)*(Tsoil(k+1) - Tsoil(k))*dt 
end do
k = Nsoil
gs(k) = ksoil(k)/Dzsoil(k)
a(k) = c(k-1)
b(k) = csoil(k) + (gs(k-1) + gs(k))*dt
c(k) = 0
rhs(k) = gs(k-1)*(Tsoil(k-1) - Tsoil(k))*dt
call TRIDIAG(Nsoil,Nsoil,a,b,c,rhs,dTs)
do k = 1, Nsoil
  Tsoil(k) = Tsoil(k) + dTs(k)
end do

end subroutine SOIL
