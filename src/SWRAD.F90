!-----------------------------------------------------------------------
! Surface and vegetation net shortwave radiation
!-----------------------------------------------------------------------
subroutine SWRAD(alb0,Dsnw,dt,elev,fcans,lveg,Sdif,Sdir,Sf,Tsrf,       &
                 albs,fsnow,SWout,SWsrf,SWsub,SWveg,tdif)

#include "OPTS.h"

use CONSTANTS, only: &
  Tm                  ! Melting point (K)

use LAYERS, only: &
  Ncnpy,             &! Number of canopy layers
  Nsmax               ! Maximum number of snow layers

use PARAMETERS, only: &
  acn0,              &! Snow-free dense canopy albedo
  acns,              &! Snow-covered dense canopy albedo
  asmx,              &! Maximum albedo for fresh snow
  asmn,              &! Minimum albedo for melting snow
  hfsn,              &! Snowcover fraction depth scale (m)
  kext,              &! Vegetation light extinction coefficient
  Salb,              &! Snowfall to refresh albedo (kg/m^2)
  Talb,              &! Snow albedo decay temperature threshold (C)
  tcld,              &! Cold snow albedo decay time scale (s)
  tmlt                ! Melting snow albedo decay time scale (s)
 
implicit none

real, intent(in) :: &
  alb0,              &! Snow-free ground albedo
  dt,                &! Timestep (s)
  elev,              &! Solar elevation (radians)
  Sdif,              &! Diffuse shortwave radiation (W/m^2)
  Sdir,              &! Direct-beam shortwave radiation (W/m^2)
  Sf,                &! Snowfall rate (kg/m2/s)
  Tsrf,              &! Snow/ground surface temperature (K)
  Dsnw(Nsmax),       &! Snow layer thicknesses (m)
  fcans(Ncnpy),      &! Canopy layer snowcover fractions
  lveg(Ncnpy)         ! Canopy layer vegetation area indices

real, intent(inout) :: &
  albs                ! Snow albedo

real, intent(out) :: &
  fsnow,             &! Ground snowcover fraction
  SWout,             &! Outgoing SW radiation (W/m^2)
  SWsrf,             &! SW absorbed by snow/ground surface (W/m^2)
  SWsub,             &! Subcanopy downward SW radiation (W/m^2)
  SWveg(Ncnpy),      &! SW absorbed by vegetation layers (W/m^2)
  tdif(Ncnpy)         ! Canopy layer diffuse transmittances

integer :: n          ! Canopy layer counter

real :: &
  alim,              &! Limiting snow albedo
  asrf,              &! Snow/ground surface albedo
  snd,               &! Snow depth (m)
  tdec                ! Snow albedo decay time scale (s)

real :: &
  A(2*Ncnpy+1,2*Ncnpy+1),   &! Canopy radiative transfer matrix
  b(2*Ncnpy+1),      &! Canopy layer boundary SW fluxes (W/m^2)
  x(2*Ncnpy+1),      &! Canopy SW sources (W/m^2)
  acan(Ncnpy),       &! Dense canopy albedo
  fdir(Ncnpy),       &! Forward-scattered fraction of direct beam
  rdif(Ncnpy),       &! Canopy layer diffuse reflectance
  rdir(Ncnpy),       &! Canopy layer direct-beam reflectance
  tdir(Ncnpy)         ! Canopy layer direct-beam transmittance

#if ALBEDO == 1
! Diagnostic snow albedo
albs = asmn + (asmx - asmn)*(Tsrf - Tm) / Talb
#elif ALBEDO == 2
! Prognostic snow albedo
tdec = tcld
if (Tsrf >= Tm) tdec = tmlt
alim = (asmn/tdec + asmx*Sf/Salb)/(1/tdec + Sf/Salb)
albs = alim + (albs - alim)*exp(-(1/tdec + Sf/Salb)*dt)
#else
stop 'Unknown option ALBEDO'
#endif
albs = max(min(albs,asmx),asmn)

! Partial snowcover on ground
snd = sum(Dsnw(:))
#if SNFRAC == 1
fsnow = min(snd/hfsn, 1.)
#elif SNFRAC == 2
fsnow = tanh(snd/hfsn)
#elif SNFRAC == 3
fsnow = snd / (snd + hfsn)
#else
stop 'Unknown option SNFRAC'
#endif

! Surface and vegetation net shortwave radiation
asrf = (1 - fsnow)*alb0 + fsnow*albs
SWveg(:) = 0
SWout = asrf*(Sdif + Sdir)
SWsub = Sdif + Sdir
tdif(:) = 0
tdir(:) = 0
if (lveg(1) > 0) then
#if CANRAD == 1
  acan(:) = (1 - fcans(:))*acn0 + fcans(:)*acns
  fdir(:) = 0
  tdif(:) = exp(-1.6*kext*lveg(:))
  tdir(:) = tdif(:)
  if (elev > 0) tdir(:) = exp(-kext*lveg(:)/sin(elev))
  rdif(:) = (1 - tdif(:))*acan(:)
  rdir(:) = (1 - tdir(:))*acan(:)
#elif CANRAD == 2
  do n =1, Ncnpy
    call TWOSTREAM(elev,fcans(n),lveg(n),fdir(n),rdif(n),rdir(n),tdif(n),tdir(n))
  end do
#else
  stop 'Unknown option CANRAD'
#endif
  A(:,:) = 0
  do n = 1, 2*Ncnpy + 1
    A(n,n) = 1
  end do
#if CANMOD == 1
  A(1,2) = -rdif(1)
  A(2,1) = -asrf
  A(3,2) = -tdif(1)
  b(1) = tdif(1)*Sdif + fdir(1)*Sdir
  b(2) = asrf*tdir(1)*Sdir
  b(3) = rdif(1)*Sdif + rdir(1)*Sdir
  call LUDCMP(3,A,b,x)
  SWout = x(3)
  SWveg(1) = Sdif - x(1) + x(2) - x(3) + (1 - tdir(1))*Sdir
  SWsub = x(1) + tdir(1)*Sdir
#elif CANMOD == 2
  A(1,4) = -rdif(1)
  A(2,1) = -tdif(2)
  A(2,3) = -rdif(2)
  A(3,2) = -asrf
  A(4,1) = -rdif(2)
  A(4,3) = -tdif(2)
  A(5,4) = -tdif(1)
  b(1) = tdif(1)*Sdif + fdir(1)*Sdir
  b(2) = fdir(2)*tdir(1)*Sdir
  b(3) = asrf*tdir(1)*tdir(2)*Sdir
  b(4) = rdir(2)*tdir(1)*Sdir
  b(5) = rdif(1)*Sdif + rdir(1)*Sdir
  call LUDCMP(5,A,b,x)
  SWout = x(5)
  SWveg(1) = Sdif - x(1) + x(4) - x(5) + (1 - tdir(1))*Sdir
  SWveg(2) = x(1) - x(2) + x(3) - x(4) + tdir(1)*(1 - tdir(2))*Sdir
  SWsub = x(2) + tdir(1)*tdir(2)*Sdir
#else
  stop 'Unknown option CANMOD'
#endif
end if
SWsrf = (1 - asrf)*SWsub

! Diffuse transmittance for LW radiation
tdif(:) = exp(-1.6*kext*lveg(:))

end subroutine SWRAD
