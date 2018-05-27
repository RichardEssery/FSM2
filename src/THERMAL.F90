!-----------------------------------------------------------------------
! Thermal properties of snow and soil
!-----------------------------------------------------------------------
subroutine THERMAL(csoil,Ds1,gs1,ks1,ksnow,ksoil,Ts1,Tveg0)

#include "OPTS.h"

use CONSTANTS, only: &
  grav,              &! Acceleration due to gravity (m/s^2)
  hcap_ice,          &! Specific heat capacity of ice (J/K/kg)
  hcap_wat,          &! Specific heat capacity of water (J/K/kg)
  hcon_air,          &! Thermal conductivity of air (W/m/K)
  hcon_ice,          &! Thermal conducivity of ice (W/m/K)
  hcon_wat,          &! Thermal conductivity of water (W/m/K)
  Lf,                &! Latent heat of fusion (J/kg)
  rho_ice,           &! Density of ice (kg/m^3)
  rho_wat,           &! Density of water (kg/m^3)
  Tm                  ! Melting point (K)

use GRID, only: &
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  bthr,              &! Snow thermal conductivity exponent
  gsat,              &! Surface conductance for saturated soil (m/s)
  hfsn,              &! Snow cover fraction depth scale (m)
  kfix,              &! Thermal conductivity at fixed snow density (W/m/K)
  rhof                ! Fresh snow density (kg/m^3)

use SOILPARAMS, only: &
  b,                 &! Clapp-Hornberger exponent
  hcap_soil,         &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil,         &! Thermal conductivity of dry soil (W/m/K)
  sathh,             &! Saturated soil water pressure (m)
  Vcrit,             &! Volumetric soil moisture content at critical point
  Vsat                ! Volumetric soil moisture content at saturation

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  theta,             &! Volumetric soil moisture content
  Tsnow,             &! Snow layer temperatures (K)
  Tsoil,             &! Soil layer temperatures (K)
  Tveg                ! Vegetation temperature (K)

implicit none

real, intent(out) :: &
  Ds1(Nx,Ny),        &! Surface layer thickness (m)
  gs1(Nx,Ny),        &! Surface moisture conductance (m/s)
  ks1(Nx,Ny),        &! Surface thermal conductivity (W/m/K)
  Ts1(Nx,Ny),        &! Surface layer temperature (K)
  Tveg0(Nx,Ny),      &! Vegetation temperature at start of timestep (K)
  csoil(Nsoil,Nx,Ny),&! Areal heat capacity of soil (J/K/m^2)
  ksnow(Nsmax,Nx,Ny),&! Thermal conductivity of snow (W/m/K)
  ksoil(Nsoil,Nx,Ny)  ! Thermal conductivity of soil (W/m/K)

integer :: &
  i,j,               &! Point counters
  k                   ! Level counter

real :: &
  dPsidT,            &! Rate of change of ice potential with temperature (m/K)
  dthudT,            &! Rate of change of unfrozen soil moisture with temperature (1/K)
  hcon_sat,          &! Thermal conductivity of saturated soil (W/m/K)
  Mf,                &! Frozen moisture content of soil layers (kg/m^2)
  Mu,                &! Unfrozen moisture content of soil layers (kg/m^2) 
  rhos,              &! Snow density (kg/m^3)
  Smf,               &! Fractional frozen soil moisture content
  Smu,               &! Fractional unfrozen soil moisture conctent
  snowdepth,         &! Snow depth (m)
  sthf,              &! Frozen soil moisture content
  sthu,              &! Unfrozen soil moisure content
  Tc,                &! Soil temperature (C)
  thice,             &! Soil ice saturation at current liquid / ice ratio
  thwat,             &! Soil water saturation at current liquid / ice ratio
  Tmax                ! Maximum temperature for frozen soil moisture (K)

! Thermal conductivity of snow
do j = 1, Ny
do i = 1, Nx
! Fixed
  ksnow(:,i,j) = kfix
#if CONDCT == 1
! Density function
  do k = 1, Nsnow(i,j)
    rhos = rhof
#if DENSTY == 1
    if (Ds(k,i,j) > epsilon(Ds)) rhos = (Sice(k,i,j) + Sliq(k,i,j)) / Ds(k,i,j)
#endif
    ksnow(k,i,j) = hcon_ice*(rhos/rho_ice)**bthr
  end do
#endif
end do
end do

! Heat capacity and thermal conductivity of soil
dPsidT = - rho_ice*Lf/(rho_wat*grav*Tm)
do j = 1, Ny
do i = 1, Nx
  do k = 1, Nsoil
    csoil(k,i,j) = hcap_soil(i,j)*Dzsoil(k)
    ksoil(k,i,j) = hcon_soil(i,j)
    if (theta(k,i,j) > epsilon(theta)) then
      dthudT = 0
      sthu = theta(k,i,j)
      sthf = 0
      Tc = Tsoil(k,i,j) - Tm
      Tmax = Tm + (sathh(i,j)/dPsidT)*(Vsat(i,j)/theta(k,i,j))**b(i,j)
      if (Tsoil(k,i,j) < Tmax) then
        dthudT = (-dPsidT*Vsat(i,j)/(b(i,j)*sathh(i,j))) *  &
                 (dPsidT*Tc/sathh(i,j))**(-1/b(i,j) - 1)
        sthu = Vsat(i,j)*(dPsidT*Tc/sathh(i,j))**(-1/b(i,j))
        sthu = min(sthu, theta(k,i,j))
        sthf = (theta(k,i,j) - sthu)*rho_wat/rho_ice
      end if
      Mf = rho_ice*Dzsoil(k)*sthf
      Mu = rho_wat*Dzsoil(k)*sthu
      csoil(k,i,j) = hcap_soil(i,j)*Dzsoil(k) + hcap_ice*Mf + hcap_wat*Mu  &
                     + rho_wat*Dzsoil(k)*((hcap_wat - hcap_ice)*Tc + Lf)*dthudT
      Smf = rho_ice*sthf/(rho_wat*Vsat(i,j))
      Smu = sthu/Vsat(i,j)
      thice = 0
      if (Smf > 0) thice = Vsat(i,j)*Smf/(Smu + Smf) 
      thwat = 0
      if (Smu > 0) thwat = Vsat(i,j)*Smu/(Smu + Smf)
      hcon_sat = hcon_soil(i,j)*(hcon_wat**thwat)*(hcon_ice**thice)/(hcon_air**Vsat(i,j))
      ksoil(k,i,j) = (hcon_sat - hcon_soil(i,j))*(Smf + Smu) + hcon_soil(i,j)
      if (k == 1) gs1(i,j) = gsat*max((Smu*Vsat(i,j)/Vcrit(i,j))**2, 1.)
    end if
  end do
end do
end do

! Surface layer
do j = 1, Ny
do i = 1, Nx
  Ds1(i,j) = max(Dzsoil(1), Ds(1,i,j))
  Ts1(i,j) = Tsoil(1,i,j) + (Tsnow(1,i,j) - Tsoil(1,i,j))*Ds(1,i,j)/Dzsoil(1)
  ks1(i,j) = Dzsoil(1) / (2*Ds(1,i,j)/ksnow(1,i,j) + (Dzsoil(1) - 2*Ds(1,i,j))/ksoil(1,i,j))
  if (Ds(1,i,j) > 0.5*Dzsoil(1)) ks1(i,j) = ksnow(1,i,j)
  if (Ds(1,i,j) > Dzsoil(1)) Ts1(i,j) = Tsnow(1,i,j)
  Tveg0(i,j) = Tveg(i,j)
end do
end do

end subroutine THERMAL
