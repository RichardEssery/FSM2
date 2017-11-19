!-----------------------------------------------------------------------
! Surface exchange coefficients
!-----------------------------------------------------------------------
subroutine SFEXCH(fsnow,KH,KHsurf,KHveg)

#include "OPTS.h"


use CONSTANTS, only: &
  g,                 &! Acceleration due to gravity (m/s^2)
  vkman               ! Von Karman constant

use DRIVING, only: &
  Ta,                &! Air temperature (K)
  Ua,                &! Wind speed (m/s)
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind measurement height (m)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only : &
  bstb,              &! Stability slope parameter
  z0sn                ! Snow roughness length (m)

use PARAMMAPS, only: &
  hcan,              &! Canopy height (m)
  fveg,              &! Canopy cover fraction
  VAI,               &! Vegetation area index
  z0sf                ! Snow-free surface roughness length (m)

use STATE_VARIABLES, only : &
  Tcan,              &! Canopy air space temperature (K)
  Tsurf               ! Surface skin temperature (K)

implicit none

real, intent(in) :: &
  fsnow(Nx,Ny)        ! Snow cover fraction

real, intent(out) :: &
  KH(Nx,Ny),         &! Eddy diffusivity for heat fluxes (m/s)
  KHsurf(Nx,Ny),     &! Surface eddy diffusivity (m/s)
  KHveg(Nx,Ny)        ! Vegetation eddy diffusivity (m/s)

integer :: &
  i,j                 ! Point counters

real :: &
  CD,                &! Drag coefficient
  dh,                &! Displacement height (m)
  fh,                &! Stability factor
  RiB,               &! Bulk Richardson number
  ustar,             &! Friction velocity (m/s)
  z0,                &! Roughness length for momentum (m)
  z0surf,            &! Ground surface roughness length (m)
  z0veg               ! Vegetation roughness length (m)

do j = 1, Ny
do i = 1, Nx

! Roughness lengths and friction velocity
  z0surf = (z0sn**fsnow(i,j)) * (z0sf(i,j)**(1 - fsnow(i,j)))
  z0veg = 0.1*hcan(i,j)
  z0 = (z0veg**fveg(i,j)) * (z0surf**(1 - fveg(i,j)))
  dh = 0.67*hcan(i,j)
  CD = (vkman / log((zU - dh)/z0))**2
  ustar = sqrt(CD)*Ua(i,j)

#if EXCHNG == 0
! No stability adjustment
  fh = 1
#endif
#if EXCHNG == 1
! Stability adjustment (Louis et al. 1982, quoted by Beljaars 1992)
  RiB = g*(Ta(i,j) - Tsurf(i,j))*zU**2 / (zT*Ta(i,j)*Ua(i,j)**2)
#if CANMOD == 1
  if (VAI(i,j) > 0) RiB = g*(Ta(i,j) - Tcan(i,j))*(zU - dh)**2 / ((zT - dh)*Ta(i,j)*Ua(i,j)**2)
#endif
  if (RiB > 0) then 
    fh = 1/(1 + 3*bstb*RiB*sqrt(1 + bstb*RiB))
  else
    fh = 1 - 3*bstb*RiB / (1 + 3*bstb**2*CD*sqrt(-RiB*zU/z0))
  end if
#endif

! Eddy diffusivities
  KH(i,j) = fh*vkman*ustar / log((zT - dh)/z0)
  KHsurf(i,j) = vkman*ustar*((1 - fveg(i,j))/log(10.) + 0.01*fveg(i,j))
  KHveg(i,j) = 0.05*VAI(i,j)*sqrt(ustar)

end do
end do

end subroutine SFEXCH
