!-----------------------------------------------------------------------
! Read and downscale meteorological driving data
!-----------------------------------------------------------------------
subroutine FSM2_DRIVE(Ncols,Nrows,fsky,lat,noon,                       &
                      year,month,day,hour,elev,EoF,                    &
                      LW,Ps,Qa,Rf,Sdif,Sdir,Sf,Ta,Ua)  
  
#include "OPTS.h"

use CONSTANTS, only: &
  eps,               &! Ratio of molecular weights of water and dry air
  e0,                &! Saturation vapour pressure at Tm (Pa)
  I0,                &! Solar constant (W/m^2)
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm                  ! Melting point (K)

use IOUNITS, only: &
  umet                ! Meteorological driving file unit number

implicit none

integer,intent(in) :: &
  Ncols,             &! Number of columns in grid
  Nrows               ! Number of rows in grid 

real, intent(in) :: &
  lat,               &! Latitude (radians)
  noon,              &! Time of solar noon (hour)
  fsky(Ncols,Nrows)   ! Skyview not obstructed by remote vegetation

integer,intent(out) :: &
  year,              &! Year
  month,             &! Month of year
  day                 ! Day of month

logical, intent(out) :: &
  EoF                 ! End-of-run flag

real, intent(out) :: &
  elev,              &! Solar elevation (radians)
  hour,              &! Hour of day
  LW(Ncols,Nrows),   &! Incoming longwave radiation (W/m^2)
  Ps(Ncols,Nrows),   &! Surface pressure (Pa)
  Qa(Ncols,Nrows),   &! Specific humidity (kg/kg)
  Rf(Ncols,Nrows),   &! Rainfall rate (kg/m^2/s)
  Sdif(Ncols,Nrows), &! Diffuse shortwave radiation (W/m^2)
  Sdir(Ncols,Nrows), &! Direct-beam shortwave radiation (W/m^2)
  Sf(Ncols,Nrows),   &! Snowfall rate (kg/m^2/s)
  Ta(Ncols,Nrows),   &! Air temperature (K)
  Ua(Ncols,Nrows)     ! Wind speed (m/s)

real :: &
  azim,              &! Solar azimuth (radians)
  dfrac,             &! Diffuse fraction of shortwave radiation
  es,                &! Saturation vapour pressure (Pa)
  Kt,                &! Sky clearness parameter
  RH,                &! Relative humidity (%)
  SW,                &! Incoming shortwave radiation (W/m^2)
  Tc                  ! Temperature (C)

#if DRIV1D == 1
! FSM driving data
read(umet,*,end=1) year,month,day,hour,SW,LW(1,1),Sf(1,1),Rf(1,1),     &
                   Ta(1,1),RH,Ua(1,1),Ps(1,1)
! Convert relative to specific humidity
Tc = Ta(1,1) - Tm
es = e0*exp(17.5043*Tc/(241.3 + Tc))
Qa(1,1) = (RH/100)*eps*es/Ps(1,1)
#endif
#if DRIV1D == 2
! ESM-SnowMIP driving data
read(umet,*,end=1) year,month,day,hour,SW,LW(1,1),Rf(1,1),Sf(1,1),     &
                   Ta(1,1),Qa(1,1),RH,Ua(1,1),Ps(1,1)
#endif

! Lower limit on wind speed
Ua(1,1) = max(Ua(1,1), 0.1)

! Copy point meteorological data to grid without downscaling
LW = LW(1,1)
Ps = Ps(1,1)
Qa = Qa(1,1)
Rf = Rf(1,1)
Sf = Sf(1,1)
Ta = Ta(1,1)
Ua = Ua(1,1)

#if SWPART == 0
! All SW radiation assumed to be diffuse
elev = 0
Sdif = SW
Sdir = 0
#endif
#if SWPART == 1
call SOLARPOS(year,month,day,hour,lat,noon,azim,elev)
! Global SW radiation partitioned into diffuse and direct components
Kt = 0
if (elev > 0) Kt = SW / (I0*sin(elev))
dfrac = 1 - 0.09*Kt
if (Kt > 0.22) dfrac = 0.95 - 0.16*Kt + 4.39*Kt**2 - 16.64*Kt**3 + 12.34*Kt**4 
if (Kt > 0.8)  dfrac = 0.165
Sdif = dfrac*SW
Sdir = (1 - dfrac)*SW
#endif

! Shading by remote vegetation
LW = fsky*LW + (1 - fsky)*sb*Ta**4
Sdif = fsky*Sdif
Sdir = fsky*Sdir  ! Replace with direct-beam transmittance

return

! End of driving data file
1 EoF = .true.

end subroutine FSM2_DRIVE
