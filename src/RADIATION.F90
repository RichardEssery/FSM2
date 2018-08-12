!-----------------------------------------------------------------------
! Surface and canopy net shortwave radiation
!-----------------------------------------------------------------------
subroutine RADIATION(alb,fsnow,SWsrf,SWveg)

#include "OPTS.h"

use CONSTANTS, only: &
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm                  ! Melting point (K)

use DRIVING, only: &
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Sf,                &! Snowfall rate (kg/m2/s)
  SW,                &! Incoming shortwave radiation (W/m2)
  Ta                  ! Air temperature (K)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  asmx,              &! Maximum albedo for fresh snow
  asmn,              &! Minimum albedo for melting snow
  avg0,              &! Snow-free vegetation albedo
  avgs,              &! Snow-covered vegetation albedo
  hfsn,              &! Snowcover fraction depth scale (m)
  kext,              &! Canopy radiation extinction coefficient
  Salb,              &! Snowfall to refresh albedo (kg/m^2)
  Talb,              &! Albedo decay temperature threshold (C)
  tcld,              &! Cold snow albedo decay timescale (s)
  tmlt                ! Melting snow albedo decay timescale (s)

use PARAMMAPS, only: &
  alb0,              &! Snow-free ground albedo
  fsky,              &! Sky view fraction
  fveg,              &! Canopy cover fraction
  scap,              &! Canopy snow capacity (kg/m^2)
  trcn,              &! Canopy transmissivity
  VAI                 ! Vegetation area index

use STATE_VARIABLES, only: &
  albs,              &! Snow albedo
  Ds,                &! Snow layer thicknesses (m)
  Sveg,              &! Snow mass on vegetation (kg/m^2)
  Tsrf                ! Surface skin temperature (K)
 
implicit none

real, intent(out) :: &
  alb(Nx,Ny),        &! Albedo
  fsnow(Nx,Ny),      &! Ground snowcover fraction
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny)        ! Net SW radiation absorbed by vegetation (W/m^2)

integer :: &
  i,j                 ! Point counters

real :: &
  alim,              &! Limiting snow albedo
  acan,              &! Canopy albedo
  asrf,              &! Surface albedo
  aveg,              &! Vegetation albedo
  fcans,             &! Canopy snowcover fraction
  rt,                &! Reciprocal timescale for albedo adjustment (1/s)
  snowdepth,         &! Snow depth (m)
  tau                 ! Snow albedo decay timescale (s)

! Nonlocal shading
do j = 1, Ny
do i = 1, Nx
  LW(i,j) = fsky(i,j)*LW(i,j) + (1 - fsky(i,j))*sb*Ta(i,j)**4
  SW(i,j) = fsky(i,j)*SW(i,j)
end do
end do

! Snow albedo
do j = 1, Ny
do i = 1, Nx
#if ALBEDO == 0
! Diagnostic
  albs(i,j) = asmn + (asmx - asmn)*(Tsrf(i,j) - Tm) / Talb
#endif
#if ALBEDO == 1
! Prognostic
  tau = tcld
  if (Tsrf(i,j) >= Tm) tau = tmlt
  rt = 1/tau + Sf(i,j)/Salb
  alim = (asmn/tau + Sf(i,j)*asmx/Salb)/rt
  albs(i,j) = alim + (albs(i,j) - alim)*exp(-rt*dt)
#endif
  if (albs(i,j) < min(asmx, asmn)) albs(i,j) = min(asmx, asmn)
  if (albs(i,j) > max(asmx, asmn)) albs(i,j) = max(asmx, asmn)
end do
end do

! Surface and canopy net shortwave radiation
do j = 1, Ny
do i = 1, Nx
! Partial snowcover on ground
  snowdepth = sum(Ds(:,i,j))
  fsnow(i,j) = tanh(snowdepth/hfsn)
  asrf = fsnow(i,j)*albs(i,j) + (1 - fsnow(i,j))*alb0(i,j)
! Partial snowcover on canopy
  acan = 0
  if (scap(i,j) > 0) then 
    fcans = Sveg(i,j) / scap(i,j)
    aveg = (1 - fcans)*avg0 + fcans*avgs
    acan = fveg(i,j)*aveg
  end if
! Effective albedo and net radiation
  alb(i,j) = acan + (1 - acan)*asrf*trcn(i,j)**2
  SWsrf(i,j) = (1 - acan)*(1 - asrf)*trcn(i,j)*SW(i,j)
  SWveg(i,j) = (1 - acan)*(1 + asrf*trcn(i,j))*(1 - trcn(i,j))*SW(i,j)
end do
end do

end subroutine RADIATION
