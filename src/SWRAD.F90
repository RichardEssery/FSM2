!-----------------------------------------------------------------------
! Surface and canopy net shortwave radiation
!-----------------------------------------------------------------------
subroutine SWRAD(alb,fcans,fsnow,SWsrf,SWveg)

#include "OPTS.h"

use CONSTANTS, only: &
  Tm                  ! Melting point (K)

use DRIVING, only: &
  month,             &! Month of year
  day,               &! Day of month
  hour,              &! Hour of day
  dt,                &! Timestep (s)
  lat,               &! Latitude (radians)
  noon,              &! Local offset from solar noon (hours)
  Sf,                &! Snowfall rate (kg/m2/s)
  SW                  ! Incoming shortwave radiation (W/m2)

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
  VAI                 ! Vegetation area index

use STATE_VARIABLES, only: &
  albs,              &! Snow albedo
  Ds,                &! Snow layer thicknesses (m)
  Sveg,              &! Snow mass on vegetation (kg/m^2)
  Tsrf                ! Surface skin temperature (K)
 
implicit none

real, intent(out) :: &
  alb(Nx,Ny),        &! Albedo
  fcans(Nx,Ny),      &! Canopy snowcover fraction
  fsnow(Nx,Ny),      &! Ground snowcover fraction
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny)        ! Net SW radiation absorbed by vegetation (W/m^2)

integer :: &
  i,j                 ! Point counters

real, parameter :: &
  pi = 3.14159        ! pi

real :: &
  SWdif(Nx,Ny),      &! Incoming diffuse shortwave radiation (W/m2)
  SWdir(Nx,Ny)        ! Incoming direct shortwave radiation (W/m2)

real :: &
  alim,              &! Limiting snow albedo
  acan,              &! Canopy albedo
  asrf,              &! Surface albedo
  aveg,              &! Vegetation albedo
  dang,              &! Day angle (radians)
  doy,               &! Day of year
  diff,              &! Diffuse fraction
  decl,              &! Solar declination (radians)
  sinelev,           &! Sine of solar elevation
  eqtm,              &! Equation of time (hours)
  hang,              &! Hour angle (radians)
  kt,                &! Atmospheric transmissivity
  rt,                &! Reciprocal timescale for albedo adjustment (1/s)
  snowdepth,         &! Snow depth (m)
  tau,               &! Snow albedo decay timescale (s)
  tdif,              &! Canopy transmission of diffuse SW radiation
  tdir                ! Canopy transmission of direct SW radiation

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

#if SWPART == 0
! SW radiation assumed to be diffuse
sinelev = 0
SWdif = SW(:,:)
SWdir = 0
#endif
#if SWPART == 1
! Calculate direct and diffuse SW radiation
doy = 275*month/9 - 3*((1 + (month - 9)/7)/100 + 1)/4 &
                  - 7*(1 + (month + 9)/12)/4 + day - 29
dang = 2*pi*(doy - 1)/365
decl = 0.006918 - 0.399912*cos(dang) + 0.070257*sin(dang)       &
		- 0.006758*cos(2*dang) + 0.000907*sin(2*dang)   &
		- 0.002697*cos(3*dang) + 0.001480*sin(3*dang)
eqtm = (0.000075 + 0.001868*cos(dang) - 0.032077*sin(dang)      &
		 - 0.014615*cos(2*dang) - 0.040849*sin(2*dang)) &
       *(229.18/60)
hang = (pi/12)*(12 + noon - hour - eqtm)
sinelev = sin(lat)*sin(decl) + cos(lat)*cos(decl)*cos(hang)
!azim = asin(cos(decl)*sin(hang)/cos(elev))
!if (sin(elev) < sin(decl)/sin(lat)) then
!  if (azim < 0) azim = azim + 2*pi
!  azim = pi - azim
!end if
do j = 1, Ny
do i = 1, Nx
  diff = 1
  if (sinelev > epsilon(sinelev)) then
    kt = SW(i,j) / (1367*sinelev)
    diff = 0.165
    if (kt < 0.22) diff = 1 - 0.09*kt
    if (kt < 0.8)  diff = 0.9511 - 0.1604*kt + 4.388*kt**2   & 
                                 - 16.638*kt**3 + 12.336*kt**4
  end if
  SWdif(i,j) = diff*SW(i,j)
  SWdir(i,j) = (1 - diff)*SW(i,j)
end do
end do
#endif

! Surface and canopy net shortwave radiation
do j = 1, Ny
do i = 1, Nx
! Partial snowcover on ground
  snowdepth = sum(Ds(:,i,j))
  fsnow(i,j) = tanh(snowdepth/hfsn)
  asrf = fsnow(i,j)*albs(i,j) + (1 - fsnow(i,j))*alb0(i,j)
! Partial snowcover on canopy
  acan = 0
  fcans(i,j) = 0
  if (scap(i,j) > 0) then 
    fcans(i,j) = Sveg(i,j) / scap(i,j)
    aveg = (1 - fcans(i,j))*avg0 + fcans(i,j)*avgs
    acan = fveg(i,j)*aveg
  end if
! Effective albedo and net radiation
  tdif = fsky(i,j)
  tdir = 0
  if (sinelev > epsilon(sinelev)) tdir = exp(-0.5*VAI(i,j)/sinelev)
  alb(i,j) = acan + (1 - acan)*asrf*tdif**2
  if (SW(i,j) > epsilon(SW))  &
    alb(i,j) = acan + (1 - acan)*asrf*tdif*(tdif*SWdif(i,j) + tdir*SWdir(i,j)) / SW(i,j)
  SWsrf(i,j) = (1 - acan)*(1 - asrf)*(tdif*SWdif(i,j) + tdir*SWdir(i,j)) 
  SWveg(i,j) = (1 - acan)*(1 + asrf*tdif)*(1 - tdif)*SWdif(i,j) +    &
               (1 - acan)*(1 - tdir + (1 - tdif)*tdir*asrf)*SWdir(i,j)
end do
end do

end subroutine SWRAD
