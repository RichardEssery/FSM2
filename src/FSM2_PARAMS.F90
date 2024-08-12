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
  eunl,              &! Exponential unloading time scale (s)
  gsnf,              &! Snow-free vegetation moisture conductance (m/s)
  hbas,              &! Canopy base height (m)
  kext,              &! Vegetation light extinction coefficient
  leaf,              &! Leaf boundary resistance (s/m)^(1/2)
  munl,              &! Melt unloading fraction
  svai,              &! Intercepted snow capacity per unit VAI (kg/m^2)
  Tunl,              &! Temperature unloading parameter (K s)
  Uunl,              &! Wind unloading parameter (m)
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
  
! Ensemble parameters
real :: &
  Pmlt,              &! Precipitation multiplier for ensemble generation
  Tadd                ! Temperature offset for ensemble generation (K)

end module PARAMETERS

!-----------------------------------------------------------------------
! Set default parameter values or read from namelist
!-----------------------------------------------------------------------
subroutine FSM2_PARAMS

#include "OPTS.h"

use PARAMETERS

implicit none

namelist /params/ acn0,acns,avg0,avgs,cvai,gsnf,hbas,kext,leaf,svai,   &
                  Tunl,Uunl,wcan,Pmlt,Tadd,                            &
                  asmn,asmx,eta0,hfsn,kfix,rcld,rfix,rgr0,rhof,rmlt,   &
                  Salb,snda,Talb,tcld,tmlt,trho,Wirr,z0sn,fcly,fsnd,   &
                  gsat,z0sf 

! Vegetation parameters
acn0 = 0.1            ! Snow-free dense canopy albedo
acns = 0.4            ! Snow-covered dense canopy albedo
avg0 = 0.27           ! Canopy element reflectivity
avgs = 0.77           ! Canopy snow reflectivity
cvai = 3.6e4          ! Vegetation heat capacity per unit VAI (J/K/m^2)
eunl = 240*3600       ! Exponential unloading time scale (s)
gsnf = 0.01           ! Snow-free vegetation moisture conductance (m/s)
hbas = 2              ! Canopy base height (m)
kext = 0.5            ! Vegetation light extinction coefficient
leaf = 20             ! Leaf boundary resistance (s/m)^(1/2)
munl = 0.4            ! Melt unloading fraction
svai = 4.4            ! Intercepted snow capacity per unit VAI (kg/m^2)
Tunl = 1.87e5         ! Temperature unloading parameter (K s)
Uunl = 1.56e5         ! Wind unloading parameter (m)
wcan = 2.5            ! Canopy wind decay coefficient

! Snow parameters
asmn = 0.5            ! Minimum albedo for melting snow
asmx = 0.85           ! Maximum albedo for fresh snow
eta0 = 3.7e7          ! Reference snow viscosity (Pa s)
hfsn = 0.1            ! Snowcover fraction depth scale (m)
kfix = 0.24           ! Fixed thermal conductivity of snow (W/m/K)
rcld = 300            ! Maximum density for cold snow (kg/m^3)
rfix = 300            ! Fixed snow density (kg/m^3)
rgr0 = 5e-5           ! Fresh snow grain radius (m)
rhof = 100            ! Fresh snow density (kg/m^3)
rhow = 300            ! Wind-packed snow density (kg/m^3)
rmlt = 500            ! Maximum density for melting snow (kg/m^3)
Salb = 10             ! Snowfall to refresh albedo (kg/m^2)
snda = 2.8e-6         ! Thermal metamorphism parameter (1/s)
Talb = -2             ! Snow albedo decay temperature threshold (C)
tcld = 1000*3600      ! Cold snow albedo decay time scale (s)
tmlt = 100*3600       ! Melting snow albedo decay time scale (s)
trho = 200*3600       ! Snow compaction timescale (s)
Wirr = 0.03           ! Irreducible liquid water content of snow
z0sn = 0.001          ! Snow roughness length (m)

! Ground surface and soil parameters
fcly = 0.3            ! Soil clay fraction
fsnd = 0.6            ! Soil sand fraction
gsat = 0.01           ! Surface conductance for saturated soil (m/s)
z0sf = 0.1            ! Snow-free surface roughness length (m)

! Ensemble parameters
Pmlt = 1              ! Precipitation multiplier for ensemble generation
Tadd = 0              ! Temperature offset for ensemble generation (K)

read(5,params)
#if DENSTY == 0
rhof = rfix
#endif

end subroutine FSM2_PARAMS
