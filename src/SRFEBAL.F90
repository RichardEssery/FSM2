!-----------------------------------------------------------------------
! Surface and canopy energy balance
!-----------------------------------------------------------------------
subroutine SRFEBAL(cveg,Ds1,dt,fcans,fsnow,gs1,ks1,lveg,LW,Ps,Qa,      &
                   SWsrf,Sveg,SWveg,Ta,tdif,Ts1,Tveg0,Ua,VAI,vegh,     &
                   zT,zU,Tsrf,Qcan,Sice,Tcan,Tveg,                     &
                   Esrf,Eveg,Gsrf,H,LE,LWout,LWsub,Melt,subl,Tsub,Usub)

#include "OPTS.h"

use CONSTANTS, only : &
  cp,                &! Specific heat capacity of air (J/K/kg)
  g,                 &! Acceleration due to gravity (m/s^2)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  Lv,                &! Latent heat of vapourisation (J/kg)
  Rair,              &! Gas constant for air (J/K/kg)
  Rwat,              &! Gas constant for water vapour (J/K/kg)
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm,                &! Melting point (K)
  vkman               ! Von Karman constant

use LAYERS, only: &
  Ncnpy,             &! Number of canopy layers
  Nsmax,             &! Maximum number of snow layers
  fvg1,              &! Fraction of vegetation in upper canopy layer
  zsub                ! Subcanopy wind speed diagnostic height (m)

use PARAMETERS, only: &
  gsnf,              &! Snow-free vegetation moisture conductance (m/s)
  hbas,              &! Canopy base height (m)
  kext,              &! Vegetation light extinction coefficient
  leaf,              &! Leaf boundary resistance (s/m)^(1/2)
  wcan,              &! Canopy wind decay coefficient
  z0sf,              &! Snow-free surface roughness length (m)
  z0sn                ! Snow roughness length (m)

implicit none

real, intent(in) :: &
  Ds1,               &! Surface layer thickness (m)
  dt,                &! Timestep (s)
  fsnow,             &! Ground snowcover fraction
  gs1,               &! Surface moisture conductance (m/s)
  ks1,               &! Surface layer thermal conductivity (W/m/K)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  SWsrf,             &! SW absorbed by snow/ground surface (W/m^2)
  Ta,                &! Air temperature (K)
  Ts1,               &! Surface layer temperature (K)
  Ua,                &! Wind speed (m/s)
  VAI,               &! Vegetation area index
  vegh,              &! Canopy height (m)
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind speed measurement height (m)

real, intent(in) :: &
  cveg(Ncnpy),       &! Vegetation heat capacities (J/K/m^2)
  fcans(Ncnpy),      &! Canopy layer snowcover fractions
  lveg(Ncnpy),       &! Canopy layer vegetation area indices
  Sveg(Ncnpy),       &! Snow mass on vegetation layers (kg/m^2)
  SWveg(Ncnpy),      &! SW absorbed by vegetation layers (W/m^2)
  tdif(Ncnpy),       &! Canopy layer diffuse transmittances
  Tveg0(Ncnpy)        ! Vegetation temperatures at start of timestep (K)

real, intent(inout) :: &
  Tsrf,              &! Snow/ground surface temperature (K)
  Qcan(Ncnpy),       &! Canopy air space humidities
  Sice(Nsmax),       &! Ice content of snow layers (kg/m^2)
  Tcan(Ncnpy),       &! Canopy air space temperatures (K)
  Tveg(Ncnpy)         ! Vegetation layer temperatures (K)

real, intent(out) :: &
  Esrf,              &! Moisture flux from the surface (kg/m^2/s)
  Gsrf,              &! Heat flux into snow/ground surface (W/m^2)
  H,                 &! Sensible heat flux to the atmosphere (W/m^2)
  LE,                &! Latent heat flux to the atmosphere (W/m^2)
  LWout,             &! Outgoing LW radiation (W/m^2)
  LWsub,             &! Subcanopy downward LW radiation (W/m^2)
  Melt,              &! Surface melt rate (kg/m^2/s)
  subl,              &! Sublimation rate (kg/m^2/s)
  Tsub,              &! Subcanopy air temperature (K)
  Usub,              &! Subcanopy wind speed (m/s)
  Eveg(Ncnpy)         ! Moisture flux from vegetation layers (kg/m^2/s)

integer :: &
  n,                 &! Canopy layer counter
  ne                  ! Energy balance iteration counter

real :: &
  d,                 &! Displacement height (m)
  Dsrf,              &! dQsat/dT at ground surface temperature (1/K)
  dEs,               &! Change in surface moisture flux (kg/m^2/s)
  dGs,               &! Change in surface heat flux (kg/m^2/s)
  dHs,               &! Change in surface sensible heat flux (kg/m^2/s)
  dTs,               &! Change in surface temperature (K)
  E,                 &! Moisture flux to the atmosphere (kg/m^2/s)
  ebal,              &! Surface energy balance closure (W/m^2)
  Ecan,              &! Within-canopy moisture flux (kg/m^2/s)
  fveg,              &! Vegetation fraction
  ga,                &! Aerodynamic conductance to the atmosphere (m/s)
  gc,                &! Conductance within canopy air space (m/s)
  gs,                &! Surface to canopy air space conductance (m/s)
  Hcan,              &! Within-canopy sensible heat flux (W/m^2)
  Hsrf,              &! Sensible heat flux from the surface (W/m^2)
  Kh,                &! Eddy diffusivity at canopy top (m^2/s)
  Lsrf,              &! Latent heat for phase change on ground (J/kg)
  psih,              &! Stability function for heat
  psim,              &! Stability function for momentum
  Qsrf,              &! Saturation humidity at surface temperature
  rd,                &! Dense vegetation aerodynamic resistance (s/m)
  rho,               &! Air density (kg/m^3)
  rL,                &! Reciprocal of Obukhov length (1/m)
  ro,                &! Open aerodynamic resistance (s/m)
  Rsrf,              &! Net radiation absorbed by the surface (W/m^2)
  Ssub,              &! Mass of snow available for sublimation (kg/m^2)
  Uc,                &! Within-canopy wind speed (m/s)
  Uh,                &! Wind speed at canopy top (m/s)
  ustar,             &! Friction velocity (m/s)
  wsrf,              &! Surface water availability factor
  zT1,               &! Temperature measurement height with offset (m)
  zU1,               &! Wind measurement height with offset (m)
  z0g,               &! Snow/ground surface roughness length (m)
  z0h,               &! Roughness length for heat (m)
  z0v                 ! Vegetation roughness length (m)

real :: &
  dEv(Ncnpy),        &! Change in vegetation moisture flux (kg/m^2/s)
  dHv(Ncnpy),        &! Change in veg sensible heat flux (kg/m^2/s)
  dQc(Ncnpy),        &! Change in canopy air humidity (kg/kg)
  dTv(Ncnpy),        &! Change in vegetation temperature (K)
  dTc(Ncnpy),        &! Change in canopy air temperature (K)
  Dveg(Ncnpy),       &! dQsat/dT at vegetation layer temperature (1/K)
  gv(Ncnpy),         &! Vegetation to canopy air space conductance (m/s)
  Hveg(Ncnpy),       &! Sensible heat flux from vegetation (W/m^2)
  Lcan(Ncnpy),       &! Latent heat for canopy water phase change (J/kg)
  Qveg(Ncnpy),       &! Saturation humidity at vegetation temperature
  Rveg(Ncnpy),       &! Net radiation absorbed by vegetation (W/m^2)
  wveg(Ncnpy),       &! Vegetation water availability factor
  zh(Ncnpy)           ! Vegetation layer heights (m)

real :: &
  J(3*Ncnpy+1,3*Ncnpy+1),&! Jacobian of energy and mass balance equations
  f(3*Ncnpy+1),          &! Residuals of energy and mass balance equations
  x(3*Ncnpy+1)            ! Temperature and humidity increments

#if ZOFFST == 0
! Heights specified above ground
zU1 = zU
zT1 = zT
#endif
#if ZOFFST == 1
! Heights specified above canopy top
zU1 = zU + vegh
zT1 = zT + vegh
#endif
#if CANMOD == 1
zh(1) = hbas + 0.5*(vegh - hbas)
#endif
#if CANMOD == 2
zh(1) = (1 - 0.5*fvg1)*vegh
zh(2) = 0.5*(1 - fvg1)*vegh
#endif

! Roughness lengths
fveg = 1 - exp(-kext*VAI)
d = 0.67*vegh
z0g = (z0sn**fsnow) * (z0sf**(1 - fsnow))
z0h = 0.1*z0g
z0v = 0.1*vegh

! Saturation humidity and air density
call QSAT(Ps,Tsrf,Qsrf)
Lsrf = Ls
if (Tsrf > Tm) Lsrf = Lv
Dsrf = Lsrf*Qsrf/(Rwat*Tsrf**2)
rho = Ps/(Rair*Ta)

if (VAI == 0) then  ! open
Eveg(:) = 0
Hveg(:) = 0
Lcan(:) = 0
LWsub = LW
rL = 0
ustar = vkman*Ua/log(zU1/z0g)
ga = vkman*ustar/log(zT1/z0h)

do ne = 1, 10
#if EXCHNG == 0
  ! No stability adjustment
#elif EXCHNG == 1
  ! Stability adjustment
  if (ne<8) rL = -vkman*g*ga*(Tsrf - Ta)/(Ta*ustar**3)
  ustar = vkman*Ua/(log(zU1/z0g) - psim(zU1,rL) + psim(z0g,rL))
  ga = vkman*ustar/(log(zT1/z0h) - psih(zT1,rL) + psih(z0h,rL))
#else
  stop 'Unknown option EXCHNG'
#endif

  ! Surface water availability
  if (Qa > Qsrf) then
    wsrf = 1
  else
    wsrf = fsnow + (1 - fsnow)*gs1/(gs1 + ga)
  end if

  ! Explicit fluxes
  Esrf = rho*wsrf*ga*(Qsrf - Qa)
  Eveg = 0
  Gsrf = 2*ks1*(Tsrf - Ts1)/Ds1
  Hsrf = cp*rho*ga*(Tsrf - Ta)
  Hveg = 0
  Melt = 0
  Rsrf = SWsrf + LW - sb*Tsrf**4

  ! Surface energy balance increments without melt
  dTs = (Rsrf - Gsrf - Hsrf - Lsrf*Esrf) /  &
        (4*sb*Tsrf**3 + 2*ks1/Ds1 + rho*(cp + Lsrf*Dsrf*wsrf)*ga)
  dEs = rho*wsrf*ga*Dsrf*dTs
  dGs = 2*ks1*dTs/Ds1 
  dHs = cp*rho*ga*dTs
  ! Surface melting
  if (Tsrf + dTs > Tm .and. Sice(1) > 0) then
    Melt = sum(Sice) / dt
    dTs = (Rsrf - Gsrf - Hsrf - Lsrf*Esrf - Lf*Melt) /  &
          (4*sb*Tsrf**3 + 2*ks1/Ds1 + rho*(cp + Ls*Dsrf*wsrf)*ga)
    dEs = rho*wsrf*ga*Dsrf*dTs
    dGs = 2*ks1*dTs/Ds1
    dHs = cp*rho*ga*dTs
    if (Tsrf + dTs < Tm) then
      call QSAT(Ps,Tm,Qsrf)
      Esrf = rho*wsrf*ga*(Qsrf - Qa)  
      Gsrf = 2*ks1*(Tm - Ts1)/Ds1
      Hsrf = cp*rho*ga*(Tm - Ta)
      Rsrf = SWsrf + LW - sb*Tm**4 
      Melt = (Rsrf - Gsrf - Hsrf - Lsrf*Esrf)/Lf
      Melt = max(Melt, 0.)
      dEs = 0
      dGs = 0
      dHs = 0
      dTs = Tm - Tsrf
    end if
  end if
  
  ! Update surface temperature and fluxes
  Esrf = Esrf + dEs
  Gsrf = Gsrf + dGs
  Hsrf = Hsrf + dHs
  Tsrf = Tsrf + dTs
  ebal = SWsrf + LW - sb*Tsrf**4 - Gsrf - Hsrf - Lsrf*Esrf - Lf*Melt
  if (ne>4 .and. abs(ebal)<0.01) exit
  
end do

else ! forest
rL = 0
ustar = fveg*vkman*Ua/log((zU1-d)/z0v) + (1 - fveg)*vkman*Ua/log(zU1/z0g)
Kh = vkman*ustar*(vegh - d)
rd = log((zT1-d)/(vegh-d))/(vkman*ustar) + vegh*(exp(wcan*(1 - zh(1)/vegh) - 1))/(wcan*Kh)
ro = log(zT1/zh(1))/(vkman*ustar)
ga = fveg/rd + (1 - fveg)/ro
do ne = 1, 10
  ! Aerodynamic resistance
#if EXCHNG == 1
  ustar = fveg*vkman*Ua/(log((zu1-d)/z0v) - psim(zU1-d,rL) + psim(z0v,rL)) + &
          (1 - fveg)*vkman*Ua/(log(zU1/z0g) - psim(zU1,rL) + psim(z0g,rL))
  if (ne<8) rL = -vkman*g*ga*(Tcan(1) - Ta)/(Ta*ustar**3)
  if (rL > 0) then
    Kh = vkman*ustar*(vegh - d)/(1 + 5*(vegh - d)*rL)
  else
    Kh = vkman*ustar*(vegh - d)*sqrt(1 - 16*(vegh - d)*rL)
  end if
  rd = (log((zT1-d)/(vegh-d)) - psih(zT1-d,rL) + psih(vegh-d,rL))/(vkman*ustar) +  &
       vegh*(exp(wcan*(1 - zh(1)/vegh) - 1))/(wcan*Kh)
  ro = (log(zT1/zh(1)) - psih(zT1,rL) + psih(zh(1),rL))/(vkman*ustar)
  ga = fveg/rd + (1 - fveg)/ro
#endif
  Uh = (ustar/vkman)*(log((vegh-d)/z0v) - psim(vegh-d,rl) + psim(z0v,rl))
  do n =1, Ncnpy
    Uc = fveg*exp(wcan*(zh(n)/vegh - 1))*Uh  +   &
         (1 - fveg)*(ustar/vkman)*(log(zh(n)/z0g) - psim(zh(n),rL) + psim(z0g,rL))
    gv(n) = sqrt(Uc)*lveg(n)/leaf
  end do
  rd = vegh*exp(wcan)*(exp(-wcan*zh(2)/vegh) - exp(-wcan*zh(1)/vegh))/(wcan*Kh)
  ro = (log(zh(1)/zh(2)) - psih(zh(1),rL) + psih(zh(2),rL))/(vkman*ustar)
  gc = fveg/rd + (1 - fveg)/ro
  n = Ncnpy
  Uc = exp(wcan*(hbas/vegh - 1))*Uh
  rd = log(hbas/z0g)*log(hbas/z0h)/(vkman**2*Uc) +  &
       vegh*exp(wcan)*(exp(-wcan*hbas/vegh) - exp(-wcan*zh(n)))/(wcan*Kh)
  ro = (log(zh(n)/z0h) - psih(zh(n),rL) + psih(z0h,rL))/(vkman*ustar)
  gs = fveg/rd + (1 - fveg)/ro

  ! Saturation humidity
  do n =1, Ncnpy
    call QSAT(Ps,Tveg(n),Qveg(n))
    Lcan(n) = Ls
    if (Tveg(n) > Tm) Lcan(n) = Lv
  end do
  Dveg(:) = Lcan(:)*Qveg(:)/(Rwat*Tveg(:)**2)

  ! Water availability
  if (Qcan(Ncnpy) > Qsrf) then
    wsrf = 1
  else
    wsrf = fsnow + (1 - fsnow)*gs1/(gs1 + gs)
  end if
  do n =1, Ncnpy
    if (Qcan(n) > Qveg(n)) then
      wveg(n) = 1
    else
      wveg(n) = fcans(n) + (1 - fcans(n))*gsnf/(gsnf + gv(n))
    end if
  end do

#if CANMOD == 1
! 1-layer canopy model

  ! Explicit fluxes
  E = rho*ga*(Qcan(1) - Qa)
  Esrf = rho*wsrf*gs*(Qsrf - Qcan(1))
  Eveg(1) = rho*wveg(1)*gv(1)*(Qveg(1) - Qcan(1))
  Gsrf = 2*ks1*(Tsrf - Ts1)/Ds1
  H = rho*cp*ga*(Tcan(1) - Ta)
  Hsrf = rho*cp*gs*(Tsrf - Tcan(1))
  Hveg(1) = rho*cp*gv(1)*(Tveg(1) - Tcan(1))
  Melt = 0
  Rsrf = SWsrf + tdif(1)*LW - sb*Tsrf**4 + (1 - tdif(1))*sb*Tveg(1)**4
  Rveg(1) = SWveg(1) + (1 - tdif(1))*(LW + sb*Tsrf**4 - 2*sb*Tveg(1)**4)

  ! Surface energy balance increments without melt
  J(1,1) = -rho*gs*(cp + Lsrf*Dsrf*wsrf) - 4*sb*Tsrf**3 - 2*ks1/Ds1
  J(1,2) = Lsrf*rho*wsrf*gs
  J(1,3) = rho*cp*gs
  J(1,4) = 4*(1 - tdif(1))*sb*Tveg(1)**3
  J(2,1) = 4*(1 - tdif(1))*sb*Tsrf**3
  J(2,2) = Lcan(1)*rho*wveg(1)*gv(1)
  J(2,3) = rho*cp*gv(1)
  J(2,4) = - rho*gv(1)*(cp + Lcan(1)*Dveg(1)*wveg(1))        &
           - 8*(1 - tdif(1))*sb*Tveg(1)**3 - cveg(1)/dt
  J(3,1) = -gs
  J(3,2) = 0
  J(3,3) = ga + gs + gv(1)
  J(3,4) = -gv(1)
  J(4,1) = -Dsrf*wsrf*gs
  J(4,2) = ga + wsrf*gs + wveg(1)*gv(1)
  J(4,3) = 0
  J(4,4) = -Dveg(1)*wveg(1)*gv(1)
  f(1)   = -(Rsrf - Gsrf - Hsrf - Lsrf*Esrf)
  f(2)   = -(Rveg(1) - Hveg(1) - Lcan(1)*Eveg(1) -           &
             cveg(1)*(Tveg(1) - Tveg0(1))/dt)
  f(3)   = -(H - Hveg(1) - Hsrf) / (rho*cp)
  f(4)   = -(E - Eveg(1) - Esrf) / rho
  call LUDCMP(4,J,f,x)
  dTs = x(1)
  dQc(1) = x(2)
  dTc(1) = x(3)
  dTv(1) = x(4)
  dEs = rho*wsrf*gs*(Dsrf*dTs - dQc(1))
  dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
  dGs = 2*ks1*dTs/Ds1
  dHs = rho*cp*gs*(dTs - dTc(1))
  dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))

  ! Surface melting
  if (Tsrf + dTs > Tm .and. Sice(1) > 0) then
    Melt = sum(Sice) / dt
    f(1) = f(1) + Lf*Melt
    call LUDCMP(4,J,f,x)
    dTs = x(1)
    dQc(1) = x(2)
    dTc(1) = x(3)
    dTv(1) = x(4)
    dEs = rho*wsrf*gs*(Dsrf*dTs - dQc(1))
    dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
    dGs = 2*ks1*dTs/Ds1
    dHs = rho*cp*gs*(dTs - dTc(1))
    dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))
    if (Tsrf + dTs < Tm) then
      call QSAT(Ps,Tm,Qsrf)
      Esrf = rho*wsrf*gs*(Qsrf - Qcan(1))
      Gsrf = 2*ks1*(Tm - Ts1)/Ds1
      Hsrf = rho*cp*gs*(Tm - Tcan(1))
      Rsrf = SWsrf + tdif(1)*LW - sb*Tm**4 + (1 - tdif(1))*sb*Tveg(1)**4
      Rveg(1) = SWveg(1) + (1 - tdif(1))*(LW + sb*Tm**4 - 2*sb*Tveg(1)**4) 
      J(1,1) = -1
      J(2,1) = 0
      J(3,1) = 0
      J(4,1) = 0
      f(1)   = -(Rsrf - Gsrf - Hsrf - Lsrf*Esrf)
      f(2)   = -(Rveg(1) - Hveg(1) - Lcan(1)*Eveg(1) -  &
                 cveg(1)*(Tveg(1) - Tveg0(1))/dt)
      f(3)   = -(H - Hveg(1) - Hsrf)/(rho*cp)
      f(4)   = -(E - Eveg(1) - Esrf)/rho
      call LUDCMP(4,J,f,x)
      Melt = x(1)/Lf
      dQc(1) = x(2)
      dTc(1) = x(3)
      dTv(1) = x(4)
      dTs = Tm - Tsrf
      dEs = 0
      dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
      dGs = 0
      dHs = 0
      dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))
    end if
  end if
  LWout = (1 - tdif(1))*sb*Tveg(1)**4 + tdif(1)*sb*Tsrf**4
  LWsub = tdif(1)*LW + (1 - tdif(1))*sb*Tveg(1)**4 
#endif

#if CANMOD == 2
! 2-layer canopy model

  ! Explicit fluxes
  E = rho*ga*(Qcan(1) - Qa)
  Ecan = rho*gc*(Qcan(2) - Qcan(1))
  Esrf = rho*wsrf*gs*(Qsrf - Qcan(2))
  Eveg(1) = rho*wveg(1)*gv(1)*(Qveg(1) - Qcan(1))
  Eveg(2) = rho*wveg(2)*gv(2)*(Qveg(2) - Qcan(2))
  Gsrf = 2*ks1*(Tsrf - Ts1)/Ds1
  H = rho*cp*ga*(Tcan(1) - Ta)
  Hcan = rho*cp*gc*(Tcan(2) - Tcan(1))
  Hsrf = rho*cp*gs*(Tsrf - Tcan(2))
  Hveg(1) = rho*cp*gv(1)*(Tveg(1) - Tcan(1))
  Hveg(2) = rho*cp*gv(2)*(Tveg(2) - Tcan(2))
  Melt = 0
  Rsrf = SWsrf + tdif(1)*tdif(2)*LW +              &
         (1 - tdif(1))*tdif(2)*sb*Tveg(1)**4 +     &
         (1 - tdif(2))*sb*Tveg(2)**4 - sb*Tsrf**4 
  Rveg(1) = SWveg(1) + (1 - tdif(1))*(LW - 2*sb*Tveg(1)**4    & 
            + (1 - tdif(2))*sb*Tveg(2)**4 + tdif(2)*sb*Tsrf**4) 
  Rveg(2) = SWveg(2) +                                               & 
            (1 - tdif(2))*(tdif(1)*LW + (1 - tdif(1))*sb*Tveg(1)**4  & 
               - 2*sb*Tveg(2)**4 + sb*Tsrf**4)

! Surface energy balance increments without melt
  J(1,1) = - rho*gs*(cp + Lsrf*Dsrf*wsrf) - 4*sb*Tsrf**3 - 2*ks1/Ds1
  J(1,2) = 0
  J(1,3) = 0
  J(1,4) = 4*(1 - tdif(1))*tdif(2)*sb*Tveg(1)**3
  J(1,5) = Lsrf*rho*wsrf*gs
  J(1,6) = rho*cp*gs
  J(1,7) = 4*(1 - tdif(2))*sb*Tveg(2)**3
  J(2,1) = 4*(1 - tdif(1))*tdif(2)*sb*Tveg(2)**3 
  J(2,2) = Lcan(1)*rho*wveg(1)*gv(1)
  J(2,3) = rho*cp*gv(1)
  J(2,4) = - rho*gv(1)*(cp + Lcan(1)*Dveg(1)*wveg(1))    & 
           - 8*(1 - tdif(1))*sb*Tveg(1)**3 - cveg(1)/dt
  J(2,5) = 0
  J(2,6) = 0
  J(2,7) = 4*(1 - tdif(1))*(1 - tdif(2))*sb*Tveg(2)**3
  J(3,1) = 4*(1 - tdif(2))*sb*Tsrf**3 
  J(3,2) = 0
  J(3,3) = 0
  J(3,4) = 4*(1 - tdif(1))*(1 - tdif(2))*sb*Tveg(1)**3
  J(3,5) = Lcan(2)*rho*wveg(2)*gv(2)
  J(3,6) = rho*cp*gv(2)
  J(3,7) = - rho*gv(2)*(cp + Lcan(2)*Dveg(2)*wveg(2))    &
           - 8*(1 - tdif(2))*sb*Tveg(2)**3 - cveg(2)/dt
  J(4,1) = 0
  J(4,2) = 0
  J(4,3) = ga + gc + gv(1)
  J(4,4) = -gv(1)
  J(4,5) = 0
  J(4,6) = -gc
  J(4,7) = 0
  J(5,1) = -gs
  J(5,2) = 0
  J(5,3) = -gc
  J(5,4) = 0
  J(5,5) = 0
  J(5,6) = gc + gs + gv(2)
  J(5,7) = -gv(2)
  J(6,1) = 0
  J(6,2) = ga + gc + wveg(1)*gv(1)
  J(6,3) = 0
  J(6,4) = -Dveg(1)*wveg(1)*gv(1)
  J(6,5) = -gc
  J(6,6) = 0
  J(6,7) = 0
  J(7,1) = -Dsrf*wsrf*gs
  J(7,2) = -gc
  J(7,3) = 0
  J(7,4) = 0
  J(7,5) = gc + wsrf*gs + wveg(2)*gv(2)
  J(7,6) = 0
  J(7,7) = -Dveg(2)*wveg(2)*gv(2)
  f(1)   = -(Rsrf - Gsrf - Hsrf - Lsrf*Esrf)
  f(2)   = -(Rveg(1) - Hveg(1) - Lcan(1)*Eveg(1) -  & 
             cveg(1)*(Tveg(1) - Tveg0(1))/dt)
  f(3)   = -(Rveg(2) - Hveg(2) - Lcan(2)*Eveg(2) -  &
             cveg(2)*(Tveg(2) - Tveg0(2))/dt)
  f(4)   = -(H - Hcan - Hveg(1))/(rho*cp)
  f(5)   = -(Hcan - Hsrf - Hveg(2))/(rho*cp)
  f(6)   = -(E - Ecan - Eveg(1))/rho
  f(7)   = -(Ecan - Esrf - Eveg(2))/rho
  call LUDCMP(7,J,f,x)
  dTs    = x(1)
  dQc(1) = x(2)
  dTc(1) = x(3)
  dTv(1) = x(4)
  dQc(2) = x(5)
  dTc(2) = x(6)
  dTv(2) = x(7)
  dEs = rho*wsrf*gs*(Dsrf*dTs - dQc(2))
  dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
  dEv(2) = rho*wveg(2)*gv(2)*(Dveg(2)*dTv(2) - dQc(2))
  dGs = 2*ks1*dTs/Ds1
  dHs = rho*cp*gs*(dTs - dTc(2))
  dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))
  dHv(2) = rho*cp*gv(2)*(dTv(2) - dTc(2))

  ! Surface melting
  if (Tsrf + dTs > Tm .and. Sice(1) > 0) then
    Melt = sum(Sice)/dt
    f(1) = f(1) + Lf*Melt
    call LUDCMP(7,J,f,x)
    dTs    = x(1)
    dQc(1) = x(2)
    dTc(1) = x(3)
    dTv(1) = x(4)
    dQc(2) = x(5)
    dTc(2) = x(6)
    dTv(2) = x(7)
    dEs = rho*wsrf*gs*(Dsrf*dTs - dQc(2))
    dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
    dEv(2) = rho*wveg(2)*gv(2)*(Dveg(2)*dTv(2) - dQc(2))
    dGs = 2*ks1*dTs/Ds1
    dHs = rho*cp*gs*(dTs - dTc(2))
    dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))
    dHv(2) = rho*cp*gv(2)*(dTv(2) - dTc(2))
    if (Tsrf + dTs < Tm) then
      call QSAT(Ps,Tm,Qsrf)
      Esrf = rho*wsrf*gs*(Qsrf - Qcan(2))
      Hsrf = rho*cp*gs*(Tm - Tcan(2))
      Rsrf = Rsrf + sb*Tsrf**4 - sb*Tm**4
      Rveg(1) = Rveg(1) + (1 - tdif(1))*tdif(2)*sb*(Tm**4 - Tsrf**4) 
      Rveg(2) = Rveg(2) + (1 - tdif(2))*sb*(Tm**4 - Tsrf**4)
      J(1,1) = -1
      J(2,1) = 0
      J(3,1) = 0
      J(4,1) = 0
      J(5,1) = 0
      J(6,1) = 0
      J(7,1) = 0
      f(1) = -(Rsrf - Gsrf - Hsrf - Lsrf*Esrf)
      f(2) = -(Rveg(1) - Hveg(1) - Lcan(1)*Eveg(1) -  & 
               cveg(1)*(Tveg(1) - Tveg0(1))/dt)
      f(3) = -(Rveg(2) - Hveg(2) - Lcan(2)*Eveg(2) -  &
               cveg(2)*(Tveg(2) - Tveg0(2))/dt)
      f(4) = -(H - Hcan - Hveg(1))/(rho*cp)
      f(5) = -(Hcan - Hsrf - Hveg(2))/(rho*cp)
      f(6) = -(E - Ecan - Eveg(1))/rho
      f(7) = -(Ecan - Esrf - Eveg(2))/rho
      call LUDCMP(7,J,f,x)
      Melt = x(1)/Lf
      dQc(1) = x(2)
      dTc(1) = x(3)
      dTv(1) = x(4)
      dQc(2) = x(5)
      dTc(2) = x(6)
      dTv(2) = x(7)
      dTs = Tm - Tsrf
      dEs = 0
      dEv(1) = rho*wveg(1)*gv(1)*(Dveg(1)*dTv(1) - dQc(1))
      dEv(2) = rho*wveg(2)*gv(2)*(Dveg(2)*dTv(2) - dQc(2))
      dGs = 0
      dHs = 0
      dHv(1) = rho*cp*gv(1)*(dTv(1) - dTc(1))
      dHv(2) = rho*cp*gv(2)*(dTv(2) - dTc(2))
    end if
  end if
  LWout = (1 - tdif(1))*sb*Tveg(1)**4 +          &
          (1 - tdif(2))*tdif(1)*sb*Tveg(1)**4 +  &
          tdif(1)*tdif(2)*sb*Tsrf**4
  LWsub = tdif(1)*tdif(2)*LW +                   &
          (1 - tdif(1))*tdif(2)*sb*Tveg(1)**4 +  &
          (1 - tdif(2))*sb*Tveg(2)**4
#endif

  ! Update vegetation temperatures and fluxes
  Eveg(:) = Eveg(:) + dEv(:)
  Hveg(:) = Hveg(:) + dHv(:)
  Qcan(:) = Qcan(:) + dQc(:)
  Tcan(:) = Tcan(:) + dTc(:)
  Tveg(:) = Tveg(:) + dTv(:)

  ! Update surface temperature and fluxes
  Esrf = Esrf + dEs
  Gsrf = Gsrf + dGs
  Hsrf = Hsrf + dHs
  Tsrf = Tsrf + dTs
  ebal = SWsrf + LWsub - sb*Tsrf**4 - Gsrf - Hsrf - Lsrf*Esrf - Lf*Melt
  
if (ne>4 .and. abs(ebal)<0.01) exit
end do
end if  ! forest

#if EBDIAG == 1
if (Melt > 0) then
  if (sum(Sice(:)) - Esrf*dt - Melt*dt > 0) then
    write(71,100) SWsrf,LWsub-sb*Tsrf**4,Gsrf,Hsrf,Lsrf*Esrf,Lf*Melt
  end if
end if
100 format(6(f16.8))
#endif

! Sublimation limited by available snow
subl = 0
Ssub = sum(Sice(:)) - Melt*dt
if (Ssub > 0 .or. Tsrf<Tm) then
  Esrf = min(Esrf, Ssub/dt)
  subl = Esrf
end if
if (VAI>0) then
  do n =1, Ncnpy
    if (Sveg(n)>0 .or. Tveg(n)<Tm) then
      Eveg(n) = min(Eveg(n), Sveg(n)/dt)
      subl = subl + Eveg(n)
    end if
  end do
end if

! Fluxes to the atmosphere
E = Esrf + sum(Eveg(:))
H = Hsrf + sum(Hveg(:))
LE = Lsrf*Esrf + sum(Lcan(:)*Eveg(:))

! Diagnostics
if (VAI==0) then
  LWout = sb*Tsrf**4
  ustar = vkman*Ua/(log(zU1/z0g) - psim(zU1,rL) + psim(z0g,rL))
  Usub = (ustar/vkman)*(log(zsub/z0g) - psim(zsub,rL) + psim(z0g,rL))
  gs = vkman*ustar/(log(zsub/z0h) - psih(zsub,rL) + psih(z0h,rL))
else
  Uc = exp(wcan*(hbas/vegh - 1))*Uh
  Usub = fveg*Uc*log(zsub/z0g)/log(hbas/z0g) +  &
         (1 - fveg)*Ua*(log(zsub/z0g) - psim(zsub,rL) + psim(z0g,rL)) / &
                       (log(zU/z0g) - psim(zU,rL) + psim(z0g,rL))
  rd = log(zsub/z0g)*log(zsub/z0h)/(vkman**2*Uc)
  ro = (log(zsub/z0h) - psih(zsub,rL) + psih(z0h,rL))/(vkman*ustar)
  gs = fveg/rd + (1 - fveg)/ro
end if
Tsub = Tsrf - Hsrf/(cp*rho*gs)

end subroutine SRFEBAL

