!-----------------------------------------------------------------------
! Thermal properties of snow and soil
!-----------------------------------------------------------------------
subroutine THERMAL(Dsnw,Nsnow,Sice,Sliq,Tsnow,Tsoil,Vsmc,              &
                   csoil,Ds1,gs1,ksnow,ksoil,ks1,Ts1)

#include "OPTS.h"

use CONSTANTS, only: &
  g,                 &! Acceleration due to gravity (m/s^2)
  hcap_ice,          &! Specific heat capacity of ice (J/K/kg)
  hcap_wat,          &! Specific heat capacity of water (J/K/kg)
  hcon_air,          &! Thermal conductivity of air (W/m/K)
  hcon_ice,          &! Thermal conducivity of ice (W/m/K)
  hcon_wat,          &! Thermal conductivity of water (W/m/K)
  Lf,                &! Latent heat of fusion (J/kg)
  rho_ice,           &! Density of ice (kg/m^3)
  rho_wat,           &! Density of water (kg/m^3)
  Tm                  ! Melting point (K)

use LAYERS, only: &
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsmax,             &! Maximum number of snow layers
  Nsoil               ! Number of soil layers

use PARAMETERS, only: &
  gsat,              &! Surface conductance for saturated soil (m/s)
  kfix,              &! Fixed thermal conductivity of snow (W/m/K)
  rhof                ! Fresh snow density (kg/m^3)

use SOILPROPS, only: &
  b,                 &! Clapp-Hornberger exponent
  hcap_soil,         &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil,         &! Thermal conductivity of dry soil (W/m/K)
  sathh,             &! Saturated soil water pressure (m)
  Vcrit,             &! Volumetric soil moisture at critical point
  Vsat                ! Volumetric soil moisture at saturation

implicit none

integer, intent(in) :: &
  Nsnow               ! Number of snow layers

real, intent(in) :: &
  Dsnw(Nsmax),       &! Snow layer thicknesses (m)
  Sice(Nsmax),       &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmax),       &! Liquid content of snow layers (kg/m^2)
  Tsnow(Nsmax),      &! Snow layer temperatures (K)
  Tsoil(Nsoil),      &! Soil layer temperatures (K)
  Vsmc(Nsoil)         ! Volumetric soil moisture content in layers

real, intent(out) :: &
  Ds1,               &! Surface layer thickness (m)
  gs1,               &! Surface moisture conductance (m/s)
  ks1,               &! Surface layer thermal conductivity (W/m/K)
  Ts1,               &! Surface layer temperature (K)
  csoil(Nsoil),      &! Areal heat capacity of soil layers (J/K/m^2)
  ksnow(Nsmax),      &! Thermal conductivity of snow layers (W/m/K)
  ksoil(Nsoil)        ! Thermal conductivity of soil layers (W/m/K)

integer :: &
  k                   ! Level counter

real :: &
  dPsidT,            &! d(ice potential)/dT (m/K)
  dthudT,            &! d(unfrozen soil moisture)/dT (1/K)
  hcon_sat,          &! Thermal conductivity of saturated soil (W/m/K)
  Mf,                &! Frozen moisture content of soil layer (kg/m^2)
  Mu,                &! Unfrozen moisture content of soil layer (kg/m^2)
  rhos,              &! Snow density (kg/m^3)
  Smf,               &! Fractional frozen soil moisture content
  Smu,               &! Fractional unfrozen soil moisture content
  snd,               &! Snow depth (m)
  sthf,              &! Frozen soil moisture content
  sthu,              &! Unfrozen soil moisure content
  Tc,                &! Soil temperature (C)
  thice,             &! Soil ice saturation at current liquid/ice ratio
  thwat,             &! Soil water saturation at current liquid/ice ratio
  Tmax                ! Maximum temperature for frozen soil moisture (K)

! Thermal conductivity of snow
ksnow = kfix
#if CONDCT == 1
do k = 1, Nsnow
  rhos = rhof
#if DENSTY == 1
  if (Dsnw(k) > epsilon(Dsnw)) rhos = (Sice(k) + Sliq(k)) / Dsnw(k)
#endif
  ksnow(k) = 2.224*(rhos/rho_wat)**1.885
end do
#endif

! Heat capacity and thermal conductivity of soil
dPsidT = - rho_ice*Lf/(rho_wat*g*Tm)
do k = 1, Nsoil
  csoil(k) = hcap_soil*Dzsoil(k)
  ksoil(k) = hcon_soil
  if (Vsmc(k) > epsilon(Vsmc)) then
    dthudT = 0
    sthu = Vsmc(k)
    sthf = 0
    Tc = Tsoil(k) - Tm
    Tmax = Tm + (sathh/dPsidT)*(Vsat/Vsmc(k))**b
    if (Tsoil(k) < Tmax) then
      dthudT = (-dPsidT*Vsat/(b*sathh)) * (dPsidT*Tc/sathh)**(-1/b - 1)
      sthu = Vsat*(dPsidT*Tc/sathh)**(-1/b)
      sthu = min(sthu, Vsmc(k))
      sthf = (Vsmc(k) - sthu)*rho_wat/rho_ice
    end if
    Mf = rho_ice*Dzsoil(k)*sthf
    Mu = rho_wat*Dzsoil(k)*sthu
    csoil(k) = hcap_soil*Dzsoil(k) + hcap_ice*Mf + hcap_wat*Mu +       &
               rho_wat*Dzsoil(k)*((hcap_wat - hcap_ice)*Tc + Lf)*dthudT
    Smf = rho_ice*sthf/(rho_wat*Vsat)
    Smu = sthu/Vsat
    thice = 0
    if (Smf > 0) thice = Vsat*Smf/(Smu + Smf) 
    thwat = 0
    if (Smu > 0) thwat = Vsat*Smu/(Smu + Smf)
    hcon_sat = hcon_soil*(hcon_wat**thwat)*(hcon_ice**thice) /         &
              (hcon_air**Vsat)
    ksoil(k) = (hcon_sat - hcon_soil)*(Smf + Smu) + hcon_soil
    if (k == 1) gs1 = gsat*max((Smu*Vsat/Vcrit)**2, 1.)
  end if
end do

! Surface layer
Ds1 = max(Dzsoil(1), Dsnw(1))
Ts1 = Tsoil(1) + (Tsnow(1) - Tsoil(1))*Dsnw(1)/Dzsoil(1)
ks1 = Dzsoil(1)/(2*Dsnw(1)/ksnow(1) + (Dzsoil(1) - 2*Dsnw(1))/ksoil(1))
snd = sum(Dsnw)
if (snd > 0.5*Dzsoil(1)) ks1 = ksnow(1)
if (snd > Dzsoil(1)) Ts1 = Tsnow(1)

end subroutine THERMAL
