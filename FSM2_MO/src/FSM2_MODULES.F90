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
! Parameters
!-----------------------------------------------------------------------
module PARAMETERS

! Vegetation parameters
real :: &
  acn0,              &! Snow-free dense canopy albedo
  acns,              &! Snow-covered dense canopy albedo
  avg0,              &! Canopy element reflectivity
  avgs,              &! Canopy snow reflectivity
  cvai,              &! Vegetation heat capacity per unit VAI (J/K/m^2)
  gsnf,              &! Snow-free vegetation moisture conductance (m/s)
  hbas,              &! Canopy base height (m)
  kext,              &! Vegetation light extinction coefficient
  leaf,              &! Leaf boundary resistance (s/m)^(1/2)
  svai,              &! Intercepted snow capacity per unit VAI (kg/m^2)
  tunl,              &! Canopy snow unloading time scale (s)
  wcan                ! Canopy wind decay coefficient

! Snow parameters
real :: &
  asmn,              &! Minimum albedo for melting snow
  asmx,              &! Maximum albedo for fresh snow
  eta0,              &! Reference snow viscosity (Pa s)
  hfsn,              &! Snowcover fraction depth scale (m)
  kfix,              &! Fixed thermal conductivity of snow (W/m/K)
  rcld,              &! Maximum density for cold snow (kg/m^3)
  rfix,              &! Fixed snow density (kg/m^3)
  rgr0,              &! Fresh snow grain radius (m)
  rhof,              &! Fresh snow density (kg/m^3)
  rhow,              &! Wind-packed snow density (kg/m^3)
  rmlt,              &! Maximum density for melting snow (kg/m^3)
  Salb,              &! Snowfall to refresh albedo (kg/m^2)
  snda,              &! Thermal metamorphism parameter (1/s)
  Talb,              &! Snow albedo decay temperature threshold (C)
  tcld,              &! Cold snow albedo decay time scale (s)
  tmlt,              &! Melting snow albedo decay time scale (s)
  trho,              &! Snow compaction timescale (s)
  Wirr,              &! Irreducible liquid water content of snow
  z0sn                ! Snow roughness length (m)

! Ground surface and soil parameters
real :: &
  fcly,              &! Soil clay fraction
  fsnd,              &! Soil sand fraction
  gsat,              &! Surface conductance for saturated soil (m/s)
  z0sf                ! Snow-free surface roughness length (m)

end module PARAMETERS

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
