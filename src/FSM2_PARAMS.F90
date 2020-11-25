!-----------------------------------------------------------------------
! Set default parameter values or read from namelist
!-----------------------------------------------------------------------
subroutine FSM2_PARAMS

#include "OPTS.h"

use PARAMETERS

implicit none

#if SETPAR == 1
namelist /params/ acn0,acns,avg0,avgs,cvai,gsnf,hbas,kext,leaf,svai,   &
                  tunl,wcan,                                           &
                  asmn,asmx,eta0,hfsn,kfix,rcld,rfix,rgr0,rhof,rmlt,   &
                  Salb,snda,Talb,tcld,tmlt,trho,Wirr,z0sn,             &
                  fcly,fsnd,gsat,z0sf 

! Vegetation parameters
acn0 = 0.1            ! Snow-free dense canopy albedo
acns = 0.4            ! Snow-covered dense canopy albedo
avg0 = 0.21           ! Canopy element reflectivity
avgs = 0.6            ! Canopy snow reflectivity
cvai = 3.6e4          ! Vegetation heat capacity per unit VAI (J/K/m^2)
gsnf = 0.01           ! Snow-free vegetation moisture conductance (m/s)
hbas = 2              ! Canopy base height (m)
kext = 0.5            ! Vegetation light extinction coefficient
leaf = 20             ! Leaf boundary resistance (s/m)^(1/2)
svai = 4.4            ! Intercepted snow capacity per unit VAI (kg/m^2)
tunl = 240*3600       ! Canopy snow unloading time scale (s)
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
tcld = 3.6e6          ! Cold snow albedo decay time scale (s)
tmlt = 3.6e5          ! Melting snow albedo decay time scale (s)
trho = 200*3600       ! Snow compaction timescale (s)
Wirr = 0.03           ! Irreducible liquid water content of snow
z0sn = 0.001          ! Snow roughness length (m)

! Ground surface and soil parameters
fcly = 0.3            ! Soil clay fraction
fsnd = 0.6            ! Soil sand fraction
gsat = 0.01           ! Surface conductance for saturated soil (m/s)
z0sf = 0.1            ! Snow-free surface roughness length (m)

read(5,params)
#if DENSTY == 0
rhof = rfix
#endif

#endif
end subroutine FSM2_PARAMS
