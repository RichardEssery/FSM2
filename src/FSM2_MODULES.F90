#include "OPTS.h"

!-----------------------------------------------------------------------
! Physical constants
!-----------------------------------------------------------------------
module CONSTANTS
real, parameter :: &
  cp = 1005,         &! Specific heat capacity of air (J/K/kg)
  eps = 0.622,       &! Ratio of molecular weights of water and dry air
  e0 = 611.213,      &! Saturation vapour pressure at Tm (Pa)
  g = 9.81,          &! Acceleration due to gravity (m/s^2)
  hcap_ice = 2100,   &! Specific heat capacity of ice (J/K/kg)
  hcap_wat = 4180,   &! Specific heat capacity of water (J/K/kg)
  hcon_air = 0.025,  &! Thermal conductivity of air (W/m/K)
  hcon_clay = 1.16,  &! Thermal conductivity of clay (W/m/K)
  hcon_ice = 2.24,   &! Thermal conducivity of ice (W/m/K)
  hcon_sand = 1.57,  &! Thermal conductivity of sand (W/m/K)
  hcon_wat = 0.56,   &! Thermal conductivity of water (W/m/K)
  I0 = 1367,         &! Solar constant (W/m^2)
  Lf = 0.334e6,      &! Latent heat of fusion (J/kg)
  Lv = 2.501e6,      &! Latent heat of vapourisation (J/kg)
  Ls = Lf + Lv,      &! Latent heat of sublimation (J/kg)
  mu_wat = 1.78e-3,  &! Dynamic viscosity of water (kg/m/s)
  pi = 3.14159,      &! pi
  Rair = 287,        &! Gas constant for air (J/K/kg)
  Rwat = 462,        &! Gas constant for water vapour (J/K/kg)
  rho_ice = 917,     &! Density of ice (kg/m^3)
  rho_wat = 1000,    &! Density of water (kg/m^3)
  sb = 5.67e-8,      &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm = 273.15,       &! Melting point (K)
  vkman = 0.4         ! von Karman constant
end module CONSTANTS

!-----------------------------------------------------------------------
! Input / output file unit numbers
!-----------------------------------------------------------------------
module IOUNITS
integer, parameter :: &
  ucan = 11,         &! Subcanopy diagnostics file unit number
  udmp = 12,         &! Start / dump file unit number
  uflx = 13,         &! Flux output file unit number
  umap = 14,         &! Map input file unit number
  umet = 15,         &! Meteorological driving file unit number
  usta = 16           ! State output file unit number
end module IOUNITS

!-----------------------------------------------------------------------
! Canopy, snow and soil layers
!-----------------------------------------------------------------------
module LAYERS
integer :: &
  Ncnpy,             &! Number of canopy layers
  Nsmax,             &! Maximum number of snow layers
  Nsoil               ! Number of soil layers
real, allocatable :: &
  Dzsnow(:),         &! Minimum snow layer thicknesses (m)
  Dzsoil(:)           ! Soil layer thicknesses (m)
real :: &
  fvg1,              &! Fraction of vegetation in upper canopy layer
  zsub                ! Subcanopy wind speed diagnostic height (m)
end module LAYERS

!-----------------------------------------------------------------------
! Soil properties
!-----------------------------------------------------------------------
module SOILPROPS
real :: &
  b,                 &! Clapp-Hornberger exponent
  hcap_soil,         &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil,         &! Thermal conductivity of dry soil (W/m/K)
  sathh,             &! Saturated soil water pressure (m)
  Vcrit,             &! Volumetric soil moisture concentration at critical point
  Vsat                ! Volumetric soil moisture concentration at saturation
end module SOILPROPS
