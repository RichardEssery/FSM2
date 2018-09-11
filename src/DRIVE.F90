!-----------------------------------------------------------------------
! Read point driving data
!-----------------------------------------------------------------------
subroutine DRIVE(EoR)

#include "OPTS.h"

use CONSTANTS, only: &
  eps,               &! Ratio of molecular weights of water and dry air
  e0,                &! Saturation vapour pressure at Tm (Pa)
  Tm                  ! Melting point (K)

use DRIVING, only: &
  year,              &! Year
  month,             &! Month of year
  day,               &! Day of month
  hour,              &! Hour of day
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Rf,                &! Rainfall rate (kg/m2/s)
  Sf,                &! Snowfall rate (kg/m2/s)
  SW,                &! Incoming shortwave radiation (W/m2)
  Ta,                &! Air temperature (K)
  Ua,                &! Wind speed (m/s)
  Pscl,              &! Precipitation adjustment scale (1/km)
  Tlps,              &! Temperature lapse rate (K/km)
  Tsnw,              &! Snow threshold temperature (K)
  zaws                ! Weather station elevation for downscaling (m)

use GRID, only: &
  Nx,Ny,             &! Grid dimensions
  ztop                ! Land surface elevations (m)

use IOUNITS, only: &
  umet                ! Driving file unit number

implicit none

logical, intent(out) :: &
  EoR                 ! End-of-run flag

integer :: & 
  i,j                 ! Point counters

! Point driving data
real :: &
  LWp,               &! Incoming longwave radiation (W/m2)
  Psp,               &! Surface pressure (Pa)
  Qap,               &! Specific humidity (kg/kg)
  Rfp,               &! Rainfall rate (kg/m2/s)
  RHp,               &! Relative humidity (%)
  Sfp,               &! Snowfall rate (kg/m2/s)
  SWp,               &! Incoming shortwave radiation (W/m2)
  Tap,               &! Air temperature (K)
  Uap                 ! Wind speed (m/s)

real :: &
  es,                &! Saturation vapour pressure (Pa)
  Tc                  ! Temperature (C)

#if DRIV1D == 0
! FSM driving data
read(umet,*,end=1) year,month,day,hour,SWp,LWp,Sfp,Rfp,Tap,RHp,Uap,Psp
Tc = Tap - Tm
es = e0*exp(17.5043*Tc/(241.3 + Tc))
Qap = (RHp/100)*eps*es/Psp
#endif

#if DRIV1D == 1
! ESM-SnowMIP driving data
read(umet,*,end=1) year,month,day,hour,SWp,LWp,Rfp,Sfp,Tap,Qap,RHp,Uap,Psp
#endif

Uap = max(Uap, 0.1)

LW(:,:) = LWp
PS(:,:) = Psp
Qa(:,:) = Qap
Rf(:,:) = Rfp
Sf(:,:) = Sfp
SW(:,:) = SWp 
Ta(:,:) = Tap
Ua(:,:) = Uap

#if DOWNSC == 1
do j = 1, Ny
do i = 1, Nx
  Ta(i,j) = Tap - Tlps*(ztop(i,j) - zaws)
  Rf(i,j) = (Rfp + Sfp)*(1 + Pscl*(ztop(i,j) - zaws))/(1 - Pscl*(ztop(i,j) - zaws))
  Sf(i,j) = 0
  if (Ta(i,j) < Tsnw) then
    Sf(i,j) = Rf(i,j)
    Rf(i,j) = 0
  end if
end do
end do
#endif

return

! End of driving data file
1 EoR = .true.

end subroutine DRIVE
