!-----------------------------------------------------------------------
! Mass baance of canopy snow
!-----------------------------------------------------------------------
subroutine CANOPY(Eveg,unload)

#include "OPTS.h"

use CONSTANTS, only: &
  Tm                  ! Melting point (K)

use DRIVING, only: &
  dt,                &! Timestep (s)
  Sf                  ! Snowfall rate (kg/m2/s)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  tcnc,              &! Canopy unloading time scale for cold snow (s)
  tcnm                ! Canopy unloading time scale for melting snow (s)

use PARAMMAPS, only: &
  fveg,              &! Canopy cover fraction
  scap                ! Canopy snow capacity (kg/m^2)

use STATE_VARIABLES, only: &
  Sveg,              &! Canopy snow mass (kg/m^2)
  Tveg                ! Vegetation temperature (K)

implicit none

real, intent(in) :: &
  Eveg(Nx,Ny)         ! Moisture flux from vegetation (kg/m^2/s)

real, intent(out) :: &
  unload(Nx,Ny)       ! Snow mass unloaded from canopy (kg/m^2)

real :: &
  intcpt,            &! Canopy interception (kg/m^2)
  Evegs,             &! Canopy snow sublimation rate (kg/m^2/s)
  tunl                ! Canopy snow unloading timescale (s)

integer :: & 
  i,j                 ! Grid coordinates

do j = 1, Ny
do i = 1, Nx
  unload(i,j) = 0
  if (fveg(i,j) > 0) then

  ! interception
    intcpt = (scap(i,j) - Sveg(i,j))*(1 - exp(-fveg(i,j)*Sf(i,j)*dt/scap(i,j)))
    Sveg(i,j) = Sveg(i,j) + intcpt
    Sf(i,j) = Sf(i,j) - intcpt/dt

  ! sublimation
    Evegs = 0
    if (Sveg(i,j) > 0 .or. Tveg(i,j) < Tm) Evegs = Eveg(i,j)
    Sveg(i,j) = Sveg(i,j) - Evegs*dt
    Sveg(i,j) = max(Sveg(i,j), 0.)

  ! unloading
    tunl = tcnc
    if (Tveg(i,j) >= Tm) tunl = tcnm
    tunl = max(tunl, dt)
    unload(i,j) = Sveg(i,j)*dt/tunl
    Sveg(i,j) = Sveg(i,j) - unload(i,j)

  end if
end do
end do

end subroutine CANOPY

