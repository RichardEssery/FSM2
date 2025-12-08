!-----------------------------------------------------------------------
! Read meteorological driving data and partition SW radiation
!-----------------------------------------------------------------------
subroutine FSM2_DRIVE(lat,noon,                                        &
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
  
use PARAMETERS, only: &
  Pmlt,              &! Precipitation multiplier for ensemble generation
  Tadd                ! Temperature offset for ensemble generation (K)

implicit none

real, intent(in) :: &
  lat,               &! Latitude (radians)
  noon                ! Time of solar noon (hour)

integer,intent(out) :: &
  year,              &! Year
  month,             &! Month of year
  day                 ! Day of month

logical, intent(out) :: &
  EoF                 ! End-of-run flag

real, intent(out) :: &
  elev,              &! Solar elevation (radians)
  hour,              &! Hour of day
  LW,                &! Incoming longwave radiation (W/m^2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Rf,                &! Rainfall rate (kg/m^2/s)
  Sdif,              &! Diffuse shortwave radiation (W/m^2)
  Sdir,              &! Direct-beam shortwave radiation (W/m^2)
  Sf,                &! Snowfall rate (kg/m^2/s)
  Ta,                &! Air temperature (K)
  Ua                  ! Wind speed (m/s)

real :: &
  azim,              &! Solar azimuth (radians)
  dfrac,             &! Diffuse fraction of shortwave radiation
  es,                &! Saturation vapour pressure (Pa)
  fs,                &! Snowfall fraction
  Kt,                &! Sky clearness parameter
  Pr,                &! Total precipitation (kg/m^2/s)
  Qs,                &! Saturation specific humidity
  RH,                &! Relative humidity (%)
  SW,                &! Incoming shortwave radiation (W/m^2)
  Tc                  ! Temperature (C)

#if DRIV1D == 1
! FSM driving data
read(umet,*,end=1) year,month,day,hour,SW,LW,Sf,Rf,Ta,RH,Ua,Ps
! Convert relative to specific humidity
Tc = Ta - Tm
es = e0*exp(17.5043*Tc/(241.3 + Tc))
Qa = (RH/100)*eps*es/Ps
#endif
#if DRIV1D == 2
! ESM-SnowMIP driving data
read(umet,*,end=1) year,month,day,hour,SW,LW,Rf,Sf,Ta,Qa,RH,Ua,Ps
#endif

! Lower limit on wind speed
Ua = max(Ua, 0.1)

#if ENSMBL == 1
! Perturbed driving data
Ta = Ta + Tadd
Pr = Pmlt*(Rf + Sf)
Tc = Ta - Tm
es = e0*exp(17.5043*Tc/(241.3 + Tc))
Qs = eps*es/Ps
Qa = min(Qa,Qs)
fs = 1 - 0.5*Tc
fs = min(1.,max(fs,0.))
Rf = (1 - fs)*Pr
Sf = fs*Pr
#endif
 
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

return

! End of driving data file
1 EoF = .true.

end subroutine FSM2_DRIVE
