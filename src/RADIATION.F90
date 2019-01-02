!-----------------------------------------------------------------------
! Surface and canopy net shortwave radiation
!-----------------------------------------------------------------------
subroutine RADIATION(alb,fsnow,SWsrf,SWveg)

#include "OPTS.h"

use CMOR, only : &
  albedo,            &! Surface albedo
  albsn,             &! Snow albedo
  rsus,              &! Surface upwelling shortwave radiation (W/m^2)
  snc                 ! Snow area fraction	

use CONSTANTS, only: &
  I0,                &! Solar constant (W/m^2)
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm                  ! Melting point (K)

use DRIVING, only: &
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Sdif,              &! Diffuse shortwave radiation (W/m^2)
  Sdir,              &! Direct-beam shortwave radiation (W/m^2)
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
  Gcn1,              &! Leaf angle distribution parameter
  Gcn2,              &! Leaf angle distribution parameter
  hfsn,              &! Snowcover fraction depth scale (m)
  kdif,              &! Diffuse radiation extinction coefficient
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
  azim,              &! Solar azimuth (radians)
  dfrac,             &! Diffuse fraction of shortwave radiation
  elev,              &! Solar elevation (radians)
  fcans,             &! Canopy snowcover fraction
  Kt,                &! Clearness parameter
  rt,                &! Reciprocal timescale for albedo adjustment (1/s)
  snowdepth,         &! Snow depth (m)
  tau,               &! Snow albedo decay timescale (s)
  tdif,              &! Canopy transmissivity for diffuse radiation
  tdir                ! Canopy transmissivity for direct-beam radiation

#if SWPART == 1
call SOLARPOS(azim,elev)
#endif

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
#if SNFRAC == 0
  fsnow(i,j) = snowdepth / (snowdepth + hfsn)
#endif
#if SNFRAC ==1
  fsnow(i,j) = tanh(snowdepth/hfsn)
#endif
  asrf = fsnow(i,j)*albs(i,j) + (1 - fsnow(i,j))*alb0(i,j)
! Partial snowcover on canopy
  fcans = 0
  if (scap(i,j) > epsilon(scap)) fcans = Sveg(i,j) / scap(i,j)
  aveg = (1 - fcans)*avg0 + fcans*avgs
  acan = fveg(i,j)*aveg
#if SWPART == 0
! SW radiation assumed to be diffuse
  Sdif = SW(i,j)
  Sdir = 0
#endif
#if SWPART == 1
! Global SW radiation partitioned into diffuse and direct components
  Kt = 0
  if (elev > 0) Kt = SW(i,j) / (I0*sin(elev))
  dfrac = 1 - 0.09*Kt
  if (Kt > 0.22) dfrac = 0.95 - 0.16*Kt + 4.39*Kt**2 - 16.64*Kt**3 + 12.34*Kt**4 
  if (Kt > 0.8)  dfrac = 0.165
  Sdif = dfrac*SW(i,j)
  Sdir = (1 - dfrac)*SW(i,j)
#endif
#if SWPART == 2
  SW(i,j) = Sdif + Sdir
#if DRIV1D != 1
  stop 'DRIV1D=1 required with SWPART=2'
#endif
#endif
  Sdif = fsky(i,j)*Sdif
  tdif = trcn(i,j)
  tdir = 0
  if (elev > 0) tdir = exp(-(Gcn1 + Gcn2*sin(elev))*VAI(i,j)/sin(elev))
! Effective albedo and net radiation
  alb(i,j) = acan + (1 - acan)*asrf*tdif**2
  if (Sdif > epsilon(Sdif))  &
    alb(i,j) = acan + (1 - acan)*asrf*tdif*(tdif*Sdif + tdir*Sdir) / (Sdif + Sdir)
  SWsrf(i,j) = (1 - acan)*(1 - asrf)*(tdif*Sdif + tdir*Sdir)
  SWveg(i,j) = (1 - acan)*(1 - (1 - asrf)*tdif - asrf*tdif*tdif)*Sdif +  &
               (1 - acan)*(1 - (1 - asrf)*tdir - asrf*tdif*tdir)*Sdir
end do
end do

! Add thermal emissions from surroundings
do j = 1, Ny
do i = 1, Nx
  LW(i,j) = fsky(i,j)*LW(i,j) + (1 - fsky(i,j))*sb*Ta(i,j)**4
end do
end do

#if TXTOUT == 1
albedo = alb(1,1)
albsn = albs(1,1)
rsus = (1 - albedo)*SW(1,1)
snc = fsnow(1,1)
#endif

end subroutine RADIATION
