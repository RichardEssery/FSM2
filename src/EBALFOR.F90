!-----------------------------------------------------------------------
! Solve surface energy balance for forested points
!-----------------------------------------------------------------------
subroutine EBALFOR(Dz1,gevap,KH,KHsurf,KHveg,ksurf,SWsurf,SWveg,Ts1, &
                   Esurf,Eveg,Gsurf,Hatmo,Latmo,Melt,Rnet)

#include "OPTS.h"

use CONSTANTS, only : &
  cp,                &! Specific heat capacity of air (J/K/kg)
  Lc,                &! Latent heat of condensation (J/kg)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  Rgas,              &! Gas constant for dry air (J/K/kg)
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
  fsky,              &! Sky view fraction
  fveg,              &! Canopy cover fraction
  scap                ! Canopy snow capacity (kg/m^2)

use STATE_VARIABLES, only : &
  Qcan,              &! Canopy air space humidity
  Sveg,              &! Snow mass on vegetation (kg/m^2)
  Sice,              &! Ice content of snow layers (kg/m^2)
  Tcan,              &! Canopy air space temperature (K)
  Tsurf,             &! Surface temperature (K)
  Tveg                ! Vegetation temperature (K)

implicit none

real, intent(in) :: &
  Dz1(Nx,Ny),        &! Surface layer thickness (m)
  gevap(Nx,Ny),      &! Surface moisture conductance (m/s)
  KH(Nx,Ny),         &! Eddy diffusivity for heat fluxes (m/s)
  KHsurf(Nx,Ny),     &! Surface eddy diffusivity (m/s)
  KHveg(Nx,Ny),      &! Vegetation eddy diffusivity (m/s)
  ksurf(Nx,Ny),      &! Surface layer thermal conductivity (W/m/K)
  SWsurf(Nx,Ny),     &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny),      &! Net SW radiation absorbed by vegetation (W/m^2)
  Ts1(Nx,Ny)          ! Surface layer temperature (K)

real, intent(out) :: &
  Esurf(Nx,Ny),      &! Moisture flux from the surface (kg/m^2/s)
  Eveg(Nx,Ny),       &! Moisture flux from vegetation (kg/m^2/s)
  Gsurf(Nx,Ny),      &! Heat flux into surface (W/m^2)
  Hatmo(Nx,Ny),      &! Sensible heat flux to the atmosphere (W/m^2)
  Latmo(Nx,Ny),      &! Latent heat flux to the atmosphere (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny)         ! Net radiation (W/m^2)

integer :: & 
  i,j                 ! Point counters

real :: &
  A(4,4),            &! Matrix of coefficients
  b(4),              &! Vector of residuals
  x(4)                ! Solution of matrix equation

real :: &
  Dsurf,             &! dQsat/dT at surface temperature (1/K)
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
  Eatmo,             &! Moisture flux to the atmosphere (kg/m^2/s)
  Hsurf,             &! Sensible heat flux from the surface (W/m^2)
  Hveg,              &! Sensible heat flux from vegetation (W/m^2)
  Lsurf,             &! Latent heat for phase change on the ground (J/kg)
  Lveg,              &! Latent heat for phase change on vegetation (J/kg)
  psis,              &! Surface moisture availability factor
  psiv,              &! Canopy moisture availability factor
  Qsurf,             &! Saturation humidity at surface temperature
  Qveg,              &! Saturation humidity at vegetation temperature
  rho,               &! Air density (kg/m^3)
  Rsurf,             &! Net radiation absorbed by the surface (W/m^2)
  Rveg                ! Net radiation absorbed by vegetation (W/m^2)model

#if CANMOD == 0
! 0-layer canopy 
do j = 1, Ny
do i = 1, Nx
  if (fveg(i,j) > 0) then

    Tveg(i,j) = Ta(i,j) 

  ! Saturation humidity, moisture availability and air density
    call QSAT(Ps(i,j),Tsurf(i,j),Qsurf)
    Lsurf = Ls
    if (Tsurf(i,j) > Tm) Lsurf = Lc
    Dsurf = Lsurf*Qsurf / (Rwat*Tsurf(i,j)**2)
    psis = gevap(i,j) / (gevap(i,j) + KH(i,j))
    if (Qsurf < Qa(i,j) .or. Sice(1,i,j) > 0) psis = 1
    call QSAT(Ps(i,j),Tveg(i,j),Qveg)
    Lveg = Ls
    if (Tveg(i,j) > Tm) Lveg = Lc
    psiv = Sveg(i,j) / Scap(i,j)
    if (Qveg < Qa(i,j) .or. psiv > 1) psiv = 1
    rho = Ps(i,j) / (Rgas*Ta(i,j))

  ! Explicit fluxes
    Esurf(i,j) = psis*rho*KH(i,j)*(Qsurf - Qa(i,j))
    Eveg(i,j)  = psiv*rho*KH(i,j)*(Qveg - Qa(i,j))
    Gsurf(i,j) = 2*ksurf(i,j)*(Tsurf(i,j) - Ts1(i,j))/Dz1(i,j)
    Hatmo(i,j) = cp*rho*KH(i,j)*(Tsurf(i,j) - Ta(i,j))
    Latmo(i,j) = Lsurf*Esurf(i,j)
    Melt(i,j) = 0
    Rsurf = SWsurf(i,j) + fsky(i,j)*LW(i,j) - sb*Tsurf(i,j)**4  &
                        + (1 - fsky(i,j))*sb*Tveg(i,j)**4

  ! Surface energy balance increments without melt
    dTs = (Rsurf - Hatmo(i,j) - Latmo(i,j) - Gsurf(i,j)) /  &
          ((cp + Lsurf*psis*Dsurf)*rho*KH(i,j) + 2*ksurf(i,j)/Dz1(i,j) + 4*sb*Tsurf(i,j)**3)
    dEs = psis*rho*KH(i,j)*Dsurf*dTs
    dGs = 2*ksurf(i,j)*dTs/Dz1(i,j)
    dHs = cp*rho*KH(i,j)*dTs

  ! Surface melting
    if (Tsurf(i,j) + dTs > Tm .and. Sice(1,i,j) > 0) then
      Melt(i,j) = sum(Sice(:,i,j))/dt
      dTs = (Rsurf - Hatmo(i,j) - Latmo(i,j) - Gsurf(i,j) - Lf*Melt(i,j)) /  &
            ((cp + Ls*Dsurf)*rho*KH(i,j) + 2*ksurf(i,j)/Dz1(i,j) + 4*sb*Tsurf(i,j)**3)
      dEs = rho*KH(i,j)*Dsurf*dTs
      dGs = 2*ksurf(i,j)*dTs/Dz1(i,j)
      dHs = cp*rho*KH(i,j)*dTs
      if (Tsurf(i,j) + dTs < Tm) then
        call QSAT(Ps(i,j),Tm,Qsurf)
        Esurf(i,j) = rho*KH(i,j)*(Qsurf - Qa(i,j))  
        Gsurf(i,j) = 2*ksurf(i,j)*(Tm - Ts1(i,j))/Dz1(i,j)
        Hatmo(i,j) = cp*rho*KH(i,j)*(Tm - Ta(i,j))
        Latmo(i,j) = Ls*Esurf(i,j)
        Rsurf = SWsurf(i,j) + fsky(i,j)*LW(i,j) - sb*Tm**4  &
                            + (1 - fsky(i,j))*sb*Tveg(i,j)**4
        Melt(i,j) = (Rsurf - Hatmo(i,j) - Latmo(i,j) - Gsurf(i,j)) / Lf
        Melt(i,j) = max(Melt(i,j), 0.)
        dEs = 0
        dGs = 0
        dHs = 0
        dTs = Tm - Tsurf(i,j)
      end if
    end if

  ! Update surface temperature and fluxes
    Tsurf(i,j) = Tsurf(i,j) + dTs
    Esurf(i,j) = Esurf(i,j) + dEs
    Gsurf(i,j) = Gsurf(i,j) + dGs
    Hatmo(i,j) = Hatmo(i,j) + dHs
    Latmo(i,j) = Lsurf*Esurf(i,j) + Lveg*Eveg(i,j)
    Rnet(i,j) = SWsurf(i,j) + SWveg(i,j) + LW(i,j) - fsky(i,j)*sb*Tsurf(i,j)**4  &
                            - (1 - fsky(i,j))*sb*Tveg(i,j)**4

  end if
end do
end do
#endif

#if CANMOD == 1
! 1-layer canopy model
do j = 1, Ny
do i = 1, Nx
  if (fveg(i,j) > 0) then

  ! Saturation humidity, moisture availability and air density
  call QSAT(Ps(i,j),Tsurf(i,j),Qsurf)
  Lsurf = Ls
  if (Tsurf(i,j) > Tm) Lsurf = Lc
  Dsurf = Lsurf*Qsurf / (Rwat*Tsurf(i,j)**2)
  psis = gevap(i,j) / (KHsurf(i,j) + gevap(i,j))
  if (Qsurf < Qcan(i,j) .or. Sice(1,i,j) > 0) psis = 1
  call QSAT(Ps(i,j),Tveg(i,j),Qveg)
  Lveg = Ls
  if (Tveg(i,j) > Tm) Lveg = Lc
  Dveg = Lveg*Qveg / (Rwat*Tveg(i,j)**2)
  psiv = Sveg(i,j) / Scap(i,j)
  if (Qveg < Qa(i,j) .or. psiv > 1) psiv = 1
  rho = Ps(i,j) / (Rgas*Ta(i,j))

  ! Explicit fluxes
  Eatmo = rho*KH(i,j)*(Qcan(i,j) - Qa(i,j))
  Esurf(i,j) = psis*rho*KHsurf(i,j)*(Qsurf - Qcan(i,j))
  Eveg(i,j) = psiv*rho*KHveg(i,j)*(Qveg - Qcan(i,j))
  Gsurf(i,j) = 2*ksurf(i,j)*(Tsurf(i,j) - Ts1(i,j))/Dz1(i,j)
  Hatmo(i,j) = rho*cp*KH(i,j)*(Tcan(i,j) - Ta(i,j))
  Hsurf = rho*cp*KHsurf(i,j)*(Tsurf(i,j) - Tcan(i,j))
  Hveg = rho*cp*KHveg(i,j)*(Tveg(i,j) - Tcan(i,j))
  Latmo(i,j) = Lsurf*Esurf(i,j) + Lveg*Eveg(i,j)
  Melt(i,j) = 0
  Rsurf = SWsurf(i,j) + fsky(i,j)*LW(i,j) - sb*Tsurf(i,j)**4 + (1 - fsky(i,j))*sb*Tveg(i,j)**4
  Rveg = SWveg(i,j) + (1 - fsky(i,j))*(LW(i,j) + sb*Tsurf(i,j)**4 - 2*sb*Tveg(i,j)**4) 

  ! Surface energy balance increments without melt
  A(1,1) = 0
  A(1,2) = - (KH(i,j) + KHveg(i,j) + KHsurf(i,j))
  A(1,3) = KHsurf(i,j)
  A(1,4) = KHveg(i,j)
  b(1)   = (Hatmo(i,j) - Hveg - Hsurf) / (rho*cp)
  A(2,1) = - (KH(i,j) + psiv*KHveg(i,j) + psis*KHsurf(i,j))
  A(2,2) = 0
  A(2,3) = psis*Dsurf*KHsurf(i,j)
  A(2,4) = psiv*Dveg*KHveg(i,j)
  b(2)   = (Eatmo - Eveg(i,j) - Esurf(i,j)) / rho
  A(3,1) = - Lsurf*psis*rho*KHsurf(i,j)
  A(3,2) = - rho*cp*KHsurf(i,j)
  A(3,3) = (cp + Lsurf*psis*Dsurf)*rho*KHsurf(i,j) + 4*sb*Tsurf(i,j)**3 + 2*ksurf(i,j)/Dz1(i,j)
  A(3,4) = - 4*(1 - fsky(i,j))*sb*Tveg(i,j)**3
  b(3)   = Rsurf - Hsurf - Lsurf*Esurf(i,j) - Gsurf(i,j)
  A(4,1) = - Lveg*psiv*rho*KHveg(i,j)
  A(4,2) = - rho*cp*KHveg(i,j)
  A(4,3) = -4*(1 - fsky(i,j))*sb*Tsurf(i,j)**3
  A(4,4) = canh(i,j)/dt + (cp + Lveg*psiv*Dveg)*rho*KHveg(i,j) + 8*(1 - fsky(i,j))*sb*Tveg(i,j)**3
  b(4)   = Rveg - Hveg - Lveg*Eveg(i,j)
  call LUDCMP(4,A,b,x)
  dQc = x(1)
  dTc = x(2)
  dTs = x(3)
  dTv = x(4)
  dEs = psis*rho*KHsurf(i,j)*(Dsurf*dTs - dQc)
  dEv = psiv*rho*KHveg(i,j)*(Dveg*dTs - dQc)
  dGs = 2*ksurf(i,j)*dTs/Dz1(i,j)
  dHs = rho*cp*KHsurf(i,j)*(dTs - dTc)
  dHv = rho*cp*KHveg(i,j)*(dTv - dTc)

  ! Surface melting
  if (Tsurf(i,j) + dTs > Tm .and. Sice(1,i,j) > 0) then
    Melt(i,j) = sum(Sice(:,i,j)) / dt
    b(3) = Rsurf - Hsurf - Lsurf*Esurf(i,j) - Gsurf(i,j) - Lf*Melt(i,j)
    call LUDCMP(4,A,b,x)
    dQc = x(1)
    dTc = x(2)
    dTs = x(3)
    dTv = x(4)
    dEs = rho*KHsurf(i,j)*(Dsurf*dTs - dQc)
    dEv = psiv*rho*KHveg(i,j)*(Dveg*dTs - dQc)
    dGs = 2*ksurf(i,j)*dTs/Dz1(i,j)
    dHs = rho*cp*KHsurf(i,j)*(dTs - dTc)
    dHv = rho*cp*KHveg(i,j)*(dTv - dTc)
    if (Tsurf(i,j) + dTs < Tm) then
      call QSAT(Ps(i,j),Tm,Qsurf)
      Esurf(i,j) = rho*KHsurf(i,j)*(Qsurf - Qcan(i,j))
      Gsurf(i,j) = 2*ksurf(i,j)*(Tm - Ts1(i,j))/Dz1(i,j)
      Hsurf = rho*cp*KHsurf(i,j)*(Tm - Tcan(i,j))
      Rsurf = SWsurf(i,j) + fsky(i,j)*LW(i,j) - sb*Tm**4 + (1 - fsky(i,j))*sb*Tveg(i,j)**4
      Rveg = SWveg(i,j) + (1 - fsky(i,j))*(LW(i,j) + sb*Tm**4 - 2*sb*Tveg(i,j)**4) 
      A(1,3) = 0
      b(1)   = (Hatmo(i,j) - Hveg - Hsurf) / (rho*cp)
      A(2,3) = 0
      b(2)   = (Eatmo - Eveg(i,j) - Esurf(i,j)) / rho
      A(3,3) = 1
      b(3)   = Rsurf - Hsurf - Lsurf*Esurf(i,j) - Gsurf(i,j)
      A(4,3) = 0
      b(4)   = Rveg - Hveg - Lveg*Eveg(i,j)
      call LUDCMP(4,A,b,x)
      dQc = x(1)
      dTc = x(2)
      Melt = x(3) / Lf
      dTv = x(4)
      dTs = Tm - Tsurf(i,j)
      dEs = 0
      dEv = psiv*rho*KHveg(i,j)*(Dveg*dTs - dQc)
      dGs = 0
      dHs = 0
      dHv = rho*cp*KHveg(i,j)*(dTv - dTc)
    end if
  end if

  ! Update temperatures and fluxes
  Qcan(i,j) = Qcan(i,j) + dQc
  Tcan(i,j) = Tcan(i,j) + dTc
  Tsurf(i,j) = Tsurf(i,j) + dTs
  Tveg(i,j) = Tveg(i,j) + dTv
  Esurf(i,j) = Esurf(i,j) + dEs
  Eveg(i,j) = Eveg(i,j) + dEv
  Gsurf(i,j) = Gsurf(i,j) + dGs
  Hsurf = Hsurf + dHs
  Hveg = Hveg + dHv
  Hatmo(i,j) = Hsurf + Hveg
  Latmo(i,j) = Lsurf*Esurf(i,j) + Lveg*Eveg(i,j)
  Rnet = SWsurf(i,j) + SWveg(i,j) + LW(i,j) - fsky(i,j)*sb*Tsurf(i,j)**4 - (1 - fsky(i,j))*sb*Tveg(i,j)**4

  end if
end do
end do
#endif

end subroutine EBALFOR

