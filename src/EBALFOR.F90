!-----------------------------------------------------------------------
! Surface and forest canopy energy balance
!-----------------------------------------------------------------------
subroutine EBALFOR(Ds1,KHa,KHg,KHv,KWg,KWv,ks1,SWsrf,SWveg,Ts1,Tveg0,  &
                   Esrf,Eveg,G,H,Hsrf,LE,LEsrf,Melt,Rnet,Rsrf)

#include "OPTS.h"

use CMOR, only : &
  rlus,              &! Surface upwelling longwave radiation (W/m^2)
  tcs,               &! Vegetation canopy temperature (K)
  ts                  ! Surface temperature (K)

use CONSTANTS, only : &
  cp,                &! Specific heat capacity of air (J/K/kg)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  Lv,                &! Latent heat of vapourisation (J/kg)
  Rair,              &! Gas constant for air (J/K/kg)
  Rwat,              &! Gas constant for water vapour (J/K/kg)
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm                  ! Melting point (K)

use DRIVING, only: &
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Ta                  ! Air temperature (K)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMMAPS, only: &
  canh,              &! Canopy heat capacity (J/K/m^2)
  fveg,              &! Canopy cover fraction
  trcn                ! Canopy transmissivity

use STATE_VARIABLES, only : &
  Qcan,              &! Canopy air space humidity
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sveg,              &! Snow mass on vegetation (kg/m^2)
  Tcan,              &! Canopy air space temperature (K)
  Tsrf,              &! Surface temperature (K)
  Tveg                ! Vegetation temperature (K)

implicit none

real, intent(in) :: &
  Ds1(Nx,Ny),        &! Surface layer thickness (m)
  KHa(Nx,Ny),        &! Eddy diffusivity from the canopy air space (m/s)
  KHg(Nx,Ny),        &! Eddy diffusivity for heat from the ground (m/s)
  KHv(Nx,Ny),        &! Eddy diffusivity for heat from vegetation (m/s)
  KWg(Nx,Ny),        &! Eddy diffusivity for water from the ground (m/s)
  KWv(Nx,Ny),        &! Eddy diffusivity for water from vegetation (m/s)
  ks1(Nx,Ny),        &! Surface layer thermal conductivity (W/m/K)
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny),      &! Net SW radiation absorbed by vegetation (W/m^2)
  Ts1(Nx,Ny),        &! Surface layer temperature (K)
  Tveg0(Nx,Ny)        ! Vegetation temperature at start of timestep (K)

real, intent(out) :: &
  Esrf(Nx,Ny),       &! Moisture flux from the surface (kg/m^2/s)
  Eveg(Nx,Ny),       &! Moisture flux from vegetation (kg/m^2/s)
  G(Nx,Ny),          &! Heat flux into surface (W/m^2)
  H(Nx,Ny),          &! Sensible heat flux to the atmosphere (W/m^2)
  Hsrf(Nx,Ny),       &! Sensible heat flux from the surface (W/m^2)
  LE(Nx,Ny),         &! Latent heat flux to the atmosphere (W/m^2)
  LEsrf(Nx,Ny),      &! Latent heat flux from the surface (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny),       &! Net radiation (W/m^2)
  Rsrf(Nx,Ny)         ! Net radiation absorbed by the surface (W/m^2)

integer :: & 
  i,j                 ! Point counters

real :: &
  A(4,4),            &! Jacobian of energy and mass balance equations
  b(4),              &! Residuals of energy and mass balance equations
  x(4)                ! Temperature and humidity increments

real :: &
  Dsrf,              &! dQsat/dT at surface temperature (1/K)
  Dveg,              &! dQsat/dT at vegetation temperature (1/K)
  dEs,               &! Change in surface moisture flux (kg/m^2/s)
  dEv,               &! Change in vegetation moisture flux (kg/m^2/s)
  dGs,               &! Change in surface heat flux (kg/m^2/s)
  dHs,               &! Change in surface sensible heat flux (kg/m^2/s)
  dHv,               &! Change in vegetation sensible heat flux (kg/m^2/s)
  dQc,               &! Change in canopy air humidity (kg/kg)
  dTc,               &! Change in canopy air temperature (K)
  dTs,               &! Change in surface temperature (K)
  dTv,               &! Change in vegetation temperature (K)
  E,                 &! Moisture flux to the atmosphere (kg/m^2/s)
  Hveg,              &! Sensible heat flux from vegetation (W/m^2)
  Lsrf,              &! Latent heat for phase change on the ground (J/kg)
  Lveg,              &! Latent heat for phase change on vegetation (J/kg)
  Qsrf,              &! Saturation humidity at surface temperature
  Qveg,              &! Saturation humidity at vegetation temperature
  rho,               &! Air density (kg/m^3)
  Rveg,              &! Net radiation absorbed by vegetation (W/m^2)model
  Ssub                ! Mass of snow available for sublimation (kg/m^2)

! 1-layer canopy model
do j = 1, Ny
do i = 1, Nx
  if (fveg(i,j) > 0) then

    ! Saturation humidity and density of air
    call QSAT(Ps(i,j),Tsrf(i,j),Qsrf)
    Lsrf = Ls
    if (Tsrf(i,j) > Tm) Lsrf = Lv
    Dsrf = Lsrf*Qsrf / (Rwat*Tsrf(i,j)**2)
    call QSAT(Ps(i,j),Tveg(i,j),Qveg)
    Lveg = Ls
    if (Tveg(i,j) > Tm) Lveg = Lv
    Dveg = Lveg*Qveg / (Rwat*Tveg(i,j)**2)
    rho = Ps(i,j) / (Rair*Ta(i,j))

    ! Explicit fluxes
    E = rho*KHa(i,j)*(Qcan(i,j) - Qa(i,j))
    Esrf(i,j) = rho*KWg(i,j)*(Qsrf - Qcan(i,j))
    Eveg(i,j) = rho*KWv(i,j)*(Qveg - Qcan(i,j))
    G(i,j) = 2*ks1(i,j)*(Tsrf(i,j) - Ts1(i,j))/Ds1(i,j)
    H(i,j) = rho*cp*KHa(i,j)*(Tcan(i,j) - Ta(i,j))
    Hsrf(i,j) = rho*cp*KHg(i,j)*(Tsrf(i,j) - Tcan(i,j))
    Hveg = rho*cp*KHv(i,j)*(Tveg(i,j) - Tcan(i,j))
    LE(i,j) = Lsrf*Esrf(i,j) + Lveg*Eveg(i,j)
    Melt(i,j) = 0
    Rsrf(i,j) = SWsrf(i,j) + trcn(i,j)*LW(i,j) - sb*Tsrf(i,j)**4 + (1 - trcn(i,j))*sb*Tveg(i,j)**4
    Rveg = SWveg(i,j) + (1 - trcn(i,j))*(LW(i,j) + sb*Tsrf(i,j)**4 - 2*sb*Tveg(i,j)**4) 

    ! Surface energy balance increments without melt
    A(1,1) = 0
    A(1,2) = - (KHa(i,j) + KHv(i,j) + KHg(i,j))
    A(1,3) = KHg(i,j)
    A(1,4) = KHv(i,j)
    b(1)   = (H(i,j) - Hveg - Hsrf(i,j)) / (rho*cp)
    A(2,1) = - (KHa(i,j) + KWv(i,j) + KWg(i,j))
    A(2,2) = 0
    A(2,3) = Dsrf*KWg(i,j)
    A(2,4) = Dveg*KWv(i,j)
    b(2)   = (E - Eveg(i,j) - Esrf(i,j)) / rho
    A(3,1) = - Lsrf*rho*KWg(i,j)
    A(3,2) = - rho*cp*KHg(i,j)
    A(3,3) = rho*(cp*KHg(i,j) + Lsrf*Dsrf*KWg(i,j)) + 4*sb*Tsrf(i,j)**3 + 2*ks1(i,j)/Ds1(i,j)
    A(3,4) = - 4*(1 - trcn(i,j))*sb*Tveg(i,j)**3
    b(3)   = Rsrf(i,j) - Hsrf(i,j) - Lsrf*Esrf(i,j) - G(i,j)
    A(4,1) = - Lveg*rho*KWv(i,j)
    A(4,2) = - rho*cp*KHv(i,j)
    A(4,3) = -4*(1 - trcn(i,j))*sb*Tsrf(i,j)**3
    A(4,4) = canh(i,j)/dt + rho*(cp*KHv(i,j) + Lveg*Dveg*KWv(i,j)) + 8*(1 - trcn(i,j))*sb*Tveg(i,j)**3
    b(4)   = Rveg - Hveg - Lveg*Eveg(i,j) - canh(i,j)*(Tveg(i,j) - Tveg0(i,j))/dt
    call LUDCMP(4,A,b,x)
    dQc = x(1)
    dTc = x(2)
    dTs = x(3)
    dTv = x(4)
    dEs = rho*KWg(i,j)*(Dsrf*dTs - dQc)
    dEv = rho*KWv(i,j)*(Dveg*dTs - dQc)
    dGs = 2*ks1(i,j)*dTs/Ds1(i,j)
    dHs = rho*cp*KHg(i,j)*(dTs - dTc)
    dHv = rho*cp*KHv(i,j)*(dTv - dTc)

    ! Surface melting
    if (Tsrf(i,j) + dTs > Tm .and. Sice(1,i,j) > 0) then
      Melt(i,j) = sum(Sice(:,i,j)) / dt
      b(3) = Rsrf(i,j) - Hsrf(i,j) - Lsrf*Esrf(i,j) - G(i,j) - Lf*Melt(i,j)
      call LUDCMP(4,A,b,x)
      dQc = x(1)
      dTc = x(2)
      dTs = x(3)
      dTv = x(4)
      dEs = rho*KWg(i,j)*(Dsrf*dTs - dQc)
      dEv = rho*KWv(i,j)*(Dveg*dTs - dQc)
      dGs = 2*ks1(i,j)*dTs/Ds1(i,j)
      dHs = rho*cp*KHg(i,j)*(dTs - dTc)
      dHv = rho*cp*KHv(i,j)*(dTv - dTc)
      if (Tsrf(i,j) + dTs < Tm) then
        call QSAT(Ps(i,j),Tm,Qsrf)
        Esrf(i,j) = rho*KWg(i,j)*(Qsrf - Qcan(i,j))
        G(i,j) = 2*ks1(i,j)*(Tm - Ts1(i,j))/Ds1(i,j)
        Hsrf(i,j) = rho*cp*KHg(i,j)*(Tm - Tcan(i,j))
        Rsrf(i,j) = SWsrf(i,j) + trcn(i,j)*LW(i,j) - sb*Tm**4 + (1 - trcn(i,j))*sb*Tveg(i,j)**4
        Rveg = SWveg(i,j) + (1 - trcn(i,j))*(LW(i,j) + sb*Tm**4 - 2*sb*Tveg(i,j)**4) 
        A(1,3) = 0
        b(1)   = (H(i,j) - Hveg - Hsrf(i,j)) / (rho*cp)
        A(2,3) = 0
        b(2)   = (E - Eveg(i,j) - Esrf(i,j)) / rho
        A(3,3) = 1
        b(3)   = Rsrf(i,j) - Hsrf(i,j) - Lsrf*Esrf(i,j) - G(i,j)
        A(4,3) = 0
        b(4)   = Rveg - Hveg - Lveg*Eveg(i,j) - canh(i,j)*(Tveg(i,j) - Tveg0(i,j))/dt
        call LUDCMP(4,A,b,x)
        dQc = x(1)
        dTc = x(2)
        Melt(i,j) = x(3) / Lf
        dTv = x(4)
        dTs = Tm - Tsrf(i,j)
        dEs = 0
        dEv = rho*KWv(i,j)*(Dveg*dTs - dQc)
        dGs = 2*ks1(i,j)*dTs/Ds1(i,j)
        dHs = 0
        dHv = rho*cp*KHv(i,j)*(dTv - dTc)
      end if
    end if

    ! Update temperatures and fluxes
    Qcan(i,j) = Qcan(i,j) + dQc
    Tcan(i,j) = Tcan(i,j) + dTc
    Tsrf(i,j) = Tsrf(i,j) + dTs
    Tveg(i,j) = Tveg(i,j) + dTv
    Esrf(i,j) = Esrf(i,j) + dEs
    Eveg(i,j) = Eveg(i,j) + dEv
    Hsrf(i,j) = Hsrf(i,j) + dHs
    Hveg = Hveg + dHv
    G(i,j) = G(i,j) + dGs
    H(i,j) = Hsrf(i,j) + Hveg
    LE(i,j) = Lsrf*Esrf(i,j) + Lveg*Eveg(i,j)
    LEsrf(i,j) = Lsrf*Esrf(i,j)
    Rnet(i,j) = SWsrf(i,j) + SWveg(i,j) +  &
                LW(i,j) - trcn(i,j)*sb*Tsrf(i,j)**4 - (1 - trcn(i,j))*sb*Tveg(i,j)**4

    ! Sublimation limited by amount of snow after melt
    Ssub = sum(Sice(:,i,j)) - Melt(i,j)*dt
    if (Ssub > 0 .and. Esrf(i,j)*dt > Ssub) then
      Esrf(i,j) = Ssub / dt
      LEsrf(i,j) = Ls*Esrf(i,j)
      Hsrf(i,j) = Rnet(i,j) - G(i,j) - LEsrf(i,j) - Lf*Melt(i,j)
    end if

#if TXTOUT == 1
   rlus = trcn(i,j)*sb*Tsrf(i,j)**4 + (1 - trcn(i,j))*sb*Tveg(i,j)**4
   tcs = Tveg(1,1)
   ts = (rlus/sb)**0.25
#endif

  end if
end do
end do

end subroutine EBALFOR

