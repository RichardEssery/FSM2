!-----------------------------------------------------------------------
! Two-stream reflectance and transmittance of a canopy layer
!-----------------------------------------------------------------------
subroutine TWOSTREAM(elev,fcans,lveg,fdir,rdif,rdir,tdif,tdir)

use PARAMETERS, only: &
  avg0,              &! Canopy element reflectivity
  avgs,              &! Canopy snow reflectivity
  kext                ! Vegetation light extinction coefficient

implicit none

real, intent(in) :: &
  elev,              &! Solar elevation (radians)
  fcans,             &! Canopy layer snowcover fraction
  lveg                ! Canopy layer vegetation area index

real, intent(out) :: &
  fdir,              &! Forward-scattered fraction of direct beam
  rdif,              &! Canopy layer diffuse reflectance
  rdir,              &! Canopy layer direct-beam reflectance
  tdif,              &! Canopy layer diffuse transmittance
  tdir                ! Canopy layer direct-beam transmittance

real :: &
  a1,a2,             &! Meador-Weaver alpha coefficients
  b1,b2,b3,          &! Direct-beam numerator terms
  beta,              &! Diffuse upscatter parameter
  beta0,             &! Direct-beam upscatter parameter
  D,                 &! Denominator
  g1,g2,g3,g4,       &! Meador-Weaver gamma coefficients
  k,                 &! Extinction coefficient
  mu,                &! Sine of solar elevation
  omega,             &! Scattering coefficient
  tau                 ! Optical thickness

omega = (1 - fcans)*avg0 + fcans*avgs
beta = 0.67
g1 = 2*(1 - (1 - beta)*omega)
g2 = 2*beta*omega
k = sqrt(g1**2 - g2**2)
tau = kext*lveg

! diffuse
D = k + g1 + (k - g1)*exp(-2*k*tau)
rdif = (g2/D)*(1 - exp(-2*k*tau))
tdif = 2*(k/D)*exp(-k*tau)

! direct beam
fdir = 0
rdir = 0
tdir = 0
if (elev > 0) then
  mu = sin(elev)
  beta0 = (0.5 + mu)*(1 - mu*log((1 + mu)/mu))
  g3 = beta0
  g4 = 1 - beta0
  a1 = g1*g4 + g2*g3
  a2 = g1*g3 + g2*g4
  D = (1 - k**2*mu**2)*((k + g1)*exp(k*tau) + (k - g1)*exp(-k*tau))
  b1 = (1 - k*mu)*(a2 + k*g3)*exp(k*tau)
  b2 = (1 + k*mu)*(a2 - k*g3)*exp(-k*tau)
  b3 = 2*k*(g3 - a2*mu)*exp(-tau/mu)
  rdir = (omega/d)*(b1 - b2 - b3)
  b1 = (1 + k*mu)*(a1 + k*g4)*exp(k*tau)
  b2 = (1 - k*mu)*(a1 - k*g4)*exp(-k*tau)
  b3 = 2*k*(g4 + a1*mu)*exp(tau/mu)
  tdir = exp(-tau/mu)
  if (tau > 30*mu) then
    fdir = 2*k*(omega/d)*(g4 + a1*mu)
  else
    fdir = exp(-tau/mu)*(omega/d)*(b2 + b3 - b1)
  end if
end if

end subroutine TWOSTREAM

