!-----------------------------------------------------------------------
! Monin-Obukhov stability functions
!-----------------------------------------------------------------------

real function psim(z,rL)  ! Stability function for momentum
use CONSTANTS, only: &
  pi                  ! pi
implicit none
real, intent(in) :: &
  rL,                &! Reciprocal of Obukhov length (1/m)
  z                   ! Height (m)
real :: &
  x,                 &! (1 - 16*z/L)^(1/4)
  zeta                ! z/L
zeta = z*rL
zeta = max(min(zeta,1.),-2.)
if (zeta > 0) then
  psim = -5*zeta
else
  x = (1 - 16*zeta)**0.25
  psim = 2*log((1 + x)/2) + log((1 + x**2)/2) - 2*atan(x) + pi/2
end if
end function psim

real function psih(z,rL)  ! Stability function for heat
implicit none
real, intent(in) :: &
  rL,                &! Reciprocal of Obukhov length (1/m)
  z                   ! Height (m)
real :: &
  x,                 &! (1 - 16*z/L)^(1/4)
  zeta                ! z/L
zeta = z*rL
zeta = max(min(zeta,1.),-2.)
if (zeta > 0) then
  psih = -5*zeta
else
  x = (1 - 16*zeta)**0.25
  psih = 2*log((1 + x**2)/2)
end if
end function psih

