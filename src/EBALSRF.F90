!-----------------------------------------------------------------------
! Surface energy balance in open areas or zero-layer forest canopy model
!-----------------------------------------------------------------------
subroutine EBALSRF(Ds1,KH,KHa,KHv,KWg,KWv,ks1,SWsrf,SWveg,Ts1, &
                   Esrf,Eveg,G,H,Hsrf,LE,LEsrf,Melt,Rnet,Rsrf)

#include "OPTS.h"

use CONSTANTS, only: &
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
  fveg,              &! Canopy cover fraction
  trcn                ! Canopy transmissivity

use STATE_VARIABLES, only: &
  Sice,              &! Ice content of snow layers (kg/m^2)
  Tsrf,              &! Surface temperature (K)
  Tveg                ! Vegetation temperature (K)

implicit none

real, intent(in) :: &
  Ds1(Nx,Ny),        &! Surface layer thickness (m)
  KH(Nx,Ny),         &! Eddy diffusivity for heat to the atmosphere (m/s)
  KHa(Nx,Ny),        &! Eddy diffusivity from the canopy air space (m/s)
  KHv(Nx,Ny),        &! Eddy diffusivity for heat from vegetation (m/s)
  KWg(Nx,Ny),        &! Eddy diffusivity for water from the ground (m/s)
  KWv(Nx,Ny),        &! Eddy diffusivity for water from vegetation (m/s)
  ks1(Nx,Ny),        &! Surface layer thermal conductivity (W/m/K)
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny),      &! Net SW radiation absorbed by vegetation (W/m^2)
  Ts1(Nx,Ny)          ! Surface layer temperature (K)

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
  D,                 &! dQsat/dT (1/K)
  dE,                &! Change in surface moisture flux (kg/m^2/s)
  dG,                &! Change in surface heat flux (W/m^2)
  dH,                &! Change in sensible heat flux (W/m^2)
  dR,                &! Change in net radiation (W/m^2)
  dTs,               &! Change in surface skin temperatures (K)
  E,                 &! Moisture flux to the atmosphere (kg/m^2/s)
  Lh,                &! Latent heat (J/kg)
  Qs,                &! Saturation humidity
  rho,               &! Air density (kg/m^3)
  Ssub                ! Mass of snow available for sublimation (kg/m^2)

do j = 1, Ny
do i = 1, Nx

#if CANMOD == 1
! Surface energy balance in forests handled by subroutine EBALFOR
  if (fveg(i,j) == 0) then
#endif

    Tveg(i,j) = Ta(i,j) 

    ! Saturation humidity and density of air
    call QSAT(Ps(i,j),Tsrf(i,j),Qs)
    Lh = Lv
    if (Tsrf(i,j) < Tm .or. Sice(1,i,j) > 0) Lh = Ls
    D = Lh*Qs/(Rwat*Tsrf(i,j)**2)
    rho = Ps(i,j) / (Rair*Ta(i,j))

    ! Explicit fluxes
    Esrf(i,j) = rho*KWg(i,j)*(Qs - Qa(i,j))
    G(i,j) = 2*ks1(i,j)*(Tsrf(i,j) - Ts1(i,j))/Ds1(i,j)
    H(i,j) = cp*rho*KH(i,j)*(Tsrf(i,j) - Ta(i,j))
    LE(i,j) = Lh*Esrf(i,j)
    Melt(i,j) = 0
    Rnet(i,j) = SWsrf(i,j) + trcn(i,j)*LW(i,j) - sb*Tsrf(i,j)**4  &
                           + (1 - trcn(i,j))*sb*Tveg(i,j)**4

    ! Surface energy balance increments without melt
    dTs = (Rnet(i,j) - G(i,j) - H(i,j) - LE(i,j)) /  &
          (4*sb*Tsrf(i,j)**3 + 2*ks1(i,j)/Ds1(i,j) + rho*(cp*KH(i,j) + Lh*D*KWg(i,j)))
    dE = rho*KWg(i,j)*D*dTs
    dG = 2*ks1(i,j)*dTs/Ds1(i,j) 
    dH = cp*rho*KH(i,j)*dTs
    dR = -4*sb*Tsrf(i,j)**3*dTs

    ! Surface melting
    if (Tsrf(i,j) + dTs > Tm .and. Sice(1,i,j) > 0) then
      Melt(i,j) = sum(Sice(:,i,j))/dt
      dTs = (Rnet(i,j) - G(i,j) - H(i,j) - LE(i,j) - Lf*Melt(i,j)) /  &
            (4*sb*Tsrf(i,j)**3 + 2*ks1(i,j)/Ds1(i,j) + rho*(cp*KH(i,j) + Ls*D*KWg(i,j)))
      dE = rho*KWg(i,j)*D*dTs
      dG = 2*ks1(i,j)*dTs/Ds1(i,j)
      dH = cp*rho*KH(i,j)*dTs
      dR = -4*sb*Tsrf(i,j)**3*dTs
      if (Tsrf(i,j) + dTs < Tm) then
          call QSAT(Ps(i,j),Tm,Qs)
          Esrf(i,j) = rho*KWg(i,j)*(Qs - Qa(i,j))  
          G(i,j) = 2*ks1(i,j)*(Tm - Ts1(i,j))/Ds1(i,j)
          H(i,j) = cp*rho*KH(i,j)*(Tm - Ta(i,j))
          LE(i,j) = Ls*Esrf(i,j)
          Rnet(i,j) = SWsrf(i,j) + trcn(i,j)*LW(i,j) - sb*Tm**4  &
                                 + (1 - trcn(i,j))*sb*Tveg(i,j)**4
          Melt(i,j) = (Rnet(i,j) - H(i,j) - LE(i,j) - G(i,j)) / Lf
          Melt(i,j) = max(Melt(i,j), 0.)
          dE = 0
          dG = 0
          dH = 0
          dR = 0
          dTs = Tm - Tsrf(i,j)
      end if

    end if

    ! Update surface temperature and fluxes
    Esrf(i,j) = Esrf(i,j) + dE
    G(i,j) = G(i,j) + dG
    H(i,j) = H(i,j) + dH
    LE(i,j) = Lh*Esrf(i,j)
    Rnet(i,j) = Rnet(i,j) + dR
    Tsrf(i,j) = Tsrf(i,j) + dTs

    ! Sublimation limited by amount of snow after melt
    Ssub = sum(Sice(:,i,j)) - Melt(i,j)*dt
    if (Ssub > 0 .and. Esrf(i,j)*dt > Ssub) then
      Esrf(i,j) = Ssub / dt
      LE(i,j) = Ls*Esrf(i,j)
      H(i,j) = Rnet(i,j) - G(i,j) - LE(i,j) - Lf*Melt(i,j)
    end if
    Hsrf(i,j) = H(i,j)
    LEsrf(i,j) = LE(i,j)
    Rsrf(i,j) = Rnet(i,j)

#if CANMOD == 0
    ! Add fluxes from canopy in zero-layer model
    Eveg(i,j) = 0
    if (fveg(i,j) > 0) then
      Eveg(i,j) = - KWv(i,j)*Esrf(i,j) / (KHa(i,j) + KWv(i,j))
      H(i,j) = KHa(i,j)*H(i,j) / (KHa(i,j) + KHv(i,j))
      Lh = Ls
      if (Tveg(i,j) > Tm) Lh = Lv
      LE(i,j) = LE(i,j) + Lh*Eveg(i,j)
      Rnet(i,j) = Rnet(i,j) + SWveg(i,j) +  &
                  (1 - trcn(i,j))*(LW(i,j) + sb*Tsrf(i,j)**4 - 2*sb*Tveg(i,j)**4)
    end if
#endif
#if CANMOD == 1
  end if
#endif

end do
end do

end subroutine EBALSRF
