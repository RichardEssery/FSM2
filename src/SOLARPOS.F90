!-----------------------------------------------------------------------
! Azimuth and elevation angles of the sun
!-----------------------------------------------------------------------
subroutine SOLARPOS(year,month,day,hour,lat,noon,azim,elev)

use CONSTANTS, only: &
  pi                  ! pi
 
implicit none

integer, intent(in) :: &
  year,              &! Year
  month,             &! Month of year
  day                 ! Day of month

real, intent(in) :: &
  hour,              &! Hour of day
  lat,               &! Latitude (radians)
  noon                ! Time of solar noon (hours)

real, intent(out) :: &
  azim,              &! Solar azimuth (radians)
  elev                ! Solar elevation (radians)

integer :: &
  DoY                 ! Day of year

real :: &
  dangle,            &! Day angle (radians)
  declin,            &! Solar declination (radians)
  eqtime,            &! Equation of time (hours)
  hangle              ! Hour angle (radians)

DoY = (7*year)/4 - 7*(year+(month+9)/12)/4 + (275*month)/9 + day - 30
dangle = 2*pi*(DoY - 1)/365
declin = 0.006918 - 0.399912*cos(dangle)   + 0.070257*sin(dangle)      &
                  - 0.006758*cos(2*dangle) + 0.000907*sin(2*dangle)    &
                  - 0.002697*cos(3*dangle) + 0.001480*sin(3*dangle)
eqtime = (0.000075 + 0.001868*cos(dangle)   - 0.032077*sin(dangle)     &
                   - 0.014615*cos(2*dangle) - 0.04089*sin(2*dangle))   &
         *(12/pi)
hangle = (pi/12)*(noon - hour - eqtime)
elev = asin(sin(declin)*sin(lat) + cos(declin)*cos(lat)*cos(hangle))
azim = acos((sin(elev)*sin(lat) - sin(declin))/(cos(elev)*cos(lat)))
if (hangle < 0) azim = - azim

end subroutine SOLARPOS
