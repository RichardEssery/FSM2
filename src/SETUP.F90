!-----------------------------------------------------------------------
! Set parameter values and initialize prognostic variables
!-----------------------------------------------------------------------
subroutine SETUP

#include "OPTS.h"

use CONSTANTS, only: &
  hcon_air,          &! Thermal conductivity of air (W/m/K)
  hcon_clay,         &! Thermal conductivity of clay (W/m/K)
  hcon_sand,         &! Thermal conductivity of sand (W/m/K)
  rho_wat,           &! Density of water (kg/m^3)
  Tm                  ! Melting point (K)

use DIAGNOSTICS, only: &
  diags,             &! Averaged diagnostics
  Nave,              &! Number of timesteps in average outputs
  Ndiags,            &! Number of averaged diagnostics
  Nsmp,              &! Timestep for sample outputs
  SWin,              &! Cumulated incoming solar radiation (J/m^2)
  SWout               ! Cumulated reflected solar radiation (J/m^2)

use DRIVING, only: &
  dt,                &! Timestep (s)
  lat,               &! Latitude (radians)
  noon,              &! Local offset from solar noon (hours)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Rf,                &! Rainfall rate (kg/m2/s)
  Sf,                &! Snowfall rate (kg/m2/s)
  SW,                &! Incoming shortwave radiation (W/m2)
  Ta,                &! Air temperature (K)
  Ua,                &! Wind speed (m/s)
  zT,                &! Temperature measurement height (m)
  zU                  ! Wind measurement height (m)

use GRID, only: &
  Dzsnow,            &! Minimum snow layer thicknesses (m)
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions

use IOUNITS, only: &
  udmp,              &! Dump file unit number
  uave,              &! Average output file unit number
  umet,              &! Driving file unit number
  umta,              &! Metadata file unit number
  usmp,              &! Sample output file unit number
  ustr                ! Start file unit number

use PARAMETERS, only: &
  asmx,              &! Maximum albedo for fresh snow
  asmn,              &! Minimum albedo for melting snow
  avg0,              &! Snow-free vegetation albedo
  avgs,              &! Snow-covered vegetation albedo
  bstb,              &! Stability slope parameter
  bthr,              &! Snow thermal conductivity exponent
  canc,              &! Canopy snow capacity per unit LAI (kg/m^2)
  cunc,              &! Canopy unloading time scale for cold snow (s)
  cunm,              &! Canopy unloading time scale for melting snow (s)
  gsat,              &! Surface conductance for saturated soil (m/s)
  hfsn,              &! Snow cover fraction depth scale (m)
  kext,              &! Canopy radiation extinction coefficient
  kfix,              &! Thermal conductivity at fixed snow density (W/m/K)
  rho0,              &! Fixed snow density (kg/m^3)
  rhof,              &! Fresh snow density (kg/m^3)
  rcld,              &! Maximum density for cold snow (kg/m^3)
  rmlt,              &! Maximum density for melting snow (kg/m^3)
  Salb,              &! Snowfall to refresh albedo (kg/m^2)
  Talb,              &! Albedo decay temperature threshold (C)
  tcld,              &! Cold snow albedo decay time scale (s)
  tmlt,              &! Melting snow albedo decay time scale (s)
  trho,              &! Snow compaction time scale (s)
  Wirr,              &! Irreducible liquid water content of snow
  z0sn                ! Snow roughness length (m)

use PARAMMAPS, only: &
  alb0,              &! Snow-free ground albedo
  canh,              &! Canopy heat capacity (J/K/m^2)
  fcly,              &! Soil clay fraction
  fsnd,              &! Soil sand fraction
  fsky,              &! Sky view fraction
  fveg,              &! Canopy cover fraction
  hcan,              &! Canopy height (m)
  scap,              &! Canopy snow capacity (kg/m^2)
  VAI,               &! Vegetation area index
  z0sf                ! Snow-free roughness length (m)

use SOILPARAMS, only: &
  b,                 &! Clapp-Hornberger exponent
  hcap_soil,         &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil,         &! Thermal conductivity of dry soil (W/m/K)
  sathh,             &! Saturated soil water pressure (m)
  Vcrit,             &! Volumetric soil moisture concentration at critical point
  Vsat                ! Volumetric soil moisture concentration at saturation

use STATE_VARIABLES, only: &
  albs,              &! Snow albedo
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers 
  Qcan,              &! Canopy air space humidity
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Sveg,              &! Canopy snow mass (kg/m^2)
  Tcan,              &! Canopy air space temperature (K)
  theta,             &! Volumetric moisture content of soil layers
  Tsnow,             &! Snow layer temperatures (K)
  Tsoil,             &! Soil layer temperatures (K)
  Tsurf,             &! Surface skin temperature (K)
  Tveg                ! Vegetation temperature (K)

implicit none

character(len=70) :: &
  ave_file,          &! Average output file name
  dmp_file,          &! Dump file name
  met_file,          &! Driving file name
  runid,             &! Run identifier
  smp_file,          &! Sample output file name
  start_file          ! Start file name

character(len=70) :: &
  alb0_file,         &! Snow-free ground albedo map file name
  canh_file,         &! Canopy heat capacity map file name
  fcly_file,         &! Soil clay fraction map file name
  fsnd_file,         &! Soil sand fraction map file name
  fsky_file,         &! Sky view fraction map file name
  fveg_file,         &! Canopy cover fraction map file name
  hcan_file,         &! Canopy height map file name
  scap_file,         &! Canopy snow capacity map file name
  VAI_file,          &! Vegetation area index map file name
  z0sf_file           ! Snow-free roughness length map file name

integer :: & 
  i,j,               &! Point counters
  k,                 &! Level counter
  time(8)             ! Date and time

real :: &
  hcon_min            ! Thermal conductivity of soil minerals (W/m/K)

real, allocatable :: &
  fsat(:),           &! Initial moisture content of soil layers as fractions of saturation
  Tprof(:)            ! Initial soil layer temperatures (K)

namelist /drive/ met_file,dt,lat,noon,zT,zU
namelist /gridpnts/ Nsmax,Nsoil,Nx,Ny
namelist /gridlevs/ Dzsnow,Dzsoil
namelist /initial/ fsat,Tprof,start_file
namelist /maps/ alb0,canh,fcly,fsnd,fsky,fveg,hcan,scap,VAI,z0sf,  &
                alb0_file,canh_file,fcly_file,fsnd_file,fsky_file, &
                fveg_file,hcan_file,scap_file,VAI_file,z0sf_file 
namelist /outputs/ Nave,Nsmp,ave_file,dmp_file,smp_file,runid
namelist /params/ asmx,asmn,avg0,avgs,bstb,bthr,canc,cunc,cunm,gsat,hfsn,kext,kfix,  &
                  rho0,rhof,rcld,rmlt,Salb,Talb,tcld,tmlt,trho,Wirr,z0sn

! Grid parameters
Nx = 1
Ny = 1
Nsmax = 3
Nsoil = 4
read(5, gridpnts)
allocate(Dzsnow(Nsmax))
allocate(Dzsoil(Nsoil))
if (Nsmax == 3) Dzsnow = (/0.1, 0.2, 0.4/)
if (Nsoil == 4) Dzsoil = (/0.1, 0.2, 0.4, 0.8/)
read(5, gridlevs)

! Driving data
met_file = 'met'
dt = 1
lat = 0
noon = 0
zT = 2
zU = 10
read(5,drive)
lat = (3.14159/180)*lat  ! convert latitude to radians
open(umet, file = met_file)
allocate(LW(Nx,Ny))
allocate(Ps(Nx,Ny))
allocate(Qa(Nx,Ny))
allocate(Rf(Nx,Ny))
allocate(Sf(Nx,Ny))
allocate(SW(Nx,Ny))
allocate(Ta(Nx,Ny))
allocate(Ua(Nx,Ny))         

! Defaults for canopy parameters
avg0 = 0.1
avgs = 0.4
canc = 4.4
cunc = 240
cunm = 2.4
kext = 0.5

! Defaults for snow parameters
asmx = 0.8
asmn = 0.5
bstb = 5
bthr = 2
hfsn = 0.1
kfix = 0.24
rho0 = 300
rhof = 100
rcld = 300
rmlt = 500
Salb = 10
Talb = -2
tcld = 1000
tmlt = 100
trho = 200
Wirr = 0.03
z0sn = 0.01

! Defaults for surface parameters
bstb = 5
gsat = 0.01

! Read parameter namelist and overwrite defaults
read(5,params)
#if DENSTY == 0
rhof = rho0
#endif

! Surface data from defaults, namelist or named map files
allocate(alb0(Nx,Ny))
allocate(canh(Nx,Ny))
allocate(fcly(Nx,Ny))
allocate(fsnd(Nx,Ny))
allocate(fsky(Nx,Ny))
allocate(fveg(Nx,Ny))
allocate(hcan(Nx,Ny))
allocate(scap(Nx,Ny))
allocate(VAI(Nx,Ny))
allocate(z0sf(Nx,Ny))
alb0_file = 'none'
canh_file = 'none'
fcly_file = 'none'
fsnd_file = 'none'
fsky_file = 'none'
fveg_file = 'none'
hcan_file = 'none'
scap_file = 'none'
VAI_file  = 'none'
z0sf_file = 'none'
alb0(:,:) = 0.2
canh(:,:) = -1
fcly(:,:) = 0.3
fsnd(:,:) = 0.6
fsky(:,:) = -1
fveg(:,:) = -1
hcan(:,:) = 0
scap(:,:) = -1
VAI(:,:)  = 0
z0sf(:,:) = 0.1
read(5,maps)
call READMAPS(alb0_file,alb0)
call READMAPS(canh_file,canh)
call READMAPS(fcly_file,fcly)
call READMAPS(fsnd_file,fsnd)
call READMAPS(fsky_file,fsky)
call READMAPS(fveg_file,fveg)
call READMAPS(hcan_file,hcan)
call READMAPS(scap_file,scap)
call READMAPS(VAI_file,VAI)
call READMAPS(z0sf_file,z0sf)
if (canh(1,1) < 0) canh(:,:) = 2500*VAI(:,:)
if (fsky(1,1) < 0) fsky(:,:) = exp(-kext*VAI(:,:))
if (fveg(1,1) < 0) fveg(:,:) = 1 - exp(-kext*VAI(:,:))
if (scap(1,1) < 0) scap(:,:) = canc*VAI(:,:)

! Derived soil parameters
allocate(b(Nx,Ny))
allocate(hcap_soil(Nx,Ny))
allocate(hcon_soil(Nx,Ny))
allocate(sathh(Nx,Ny))
allocate(Vsat(Nx,Ny))
allocate(Vcrit(Nx,Ny))
do j = 1, Ny
do i = 1, Nx
  if (fcly(i,j) + fsnd(i,j) > 1) fcly(i,j) = 1 - fsnd(i,j)
  b(i,j) = 3.1 + 15.7*fcly(i,j) - 0.3*fsnd(i,j)
  hcap_soil(i,j) = (2.128*fcly(i,j) + 2.385*fsnd(i,j))*1E6 / (fcly(i,j) + fsnd(i,j))
  sathh(i,j) = 10**(0.17 - 0.63*fcly(i,j) - 1.58*fsnd(i,j))
  Vsat(i,j) = 0.505 - 0.037*fcly(i,j) - 0.142*fsnd(i,j)
  Vcrit(i,j) = Vsat(i,j)*(sathh(i,j)/3.364)**(1/b(i,j))
  hcon_min = (hcon_clay**fcly(i,j)) * (hcon_sand**(1 - fcly(i,j)))
  hcon_soil(i,j) = (hcon_air**Vsat(i,j)) * (hcon_min**(1 - Vsat(i,j)))
end do
end do

! Convert time scales from hours to seconds
dt = 3600*dt
cunc = max(3600*cunc, dt)
cunm = max(3600*cunm, dt)
tcld = 3600*tcld
tmlt = 3600*tmlt
trho = 3600*trho

! Allocate state variables
allocate(albs(Nx,Ny))
allocate(Ds(Nsmax,Nx,Ny))
allocate(Nsnow(Nx,Ny))
allocate(Qcan(Nx,Ny))
allocate(Sice(Nsmax,Nx,Ny))
allocate(Sliq(Nsmax,Nx,Ny))
allocate(Sveg(Nx,Ny))
allocate(Tcan(Nx,Ny))
allocate(theta(Nsoil,Nx,Ny))
allocate(Tsnow(Nsmax,Nx,Ny))
allocate(Tsoil(Nsoil,Nx,Ny))
allocate(Tsurf(Nx,Ny))
allocate(Tveg(Nx,Ny))

! Default initialization of state variables
albs(:,:)    = 0.8
Ds(:,:,:)    = 0
Nsnow(:,:)   = 0
Qcan(:,:)    = 0
Sice(:,:,:)  = 0
Sliq(:,:,:)  = 0
Sveg(:,:)    = 0
Tcan(:,:)    = 285
Tsnow(:,:,:) = Tm
Tsoil(:,:,:) = 285
Tveg(:,:)    = 285

! Initial soil profiles from namelist
allocate(fsat(Nsoil))
allocate(Tprof(Nsoil))
fsat(:)  = 0.5
Tprof(:) = 285
start_file = 'none'
read(5, initial)
do k = 1, Nsoil
  theta(k,:,:) = fsat(k)*Vsat(:,:)
  Tsoil(k,:,:) = Tprof(k)
end do
Tsurf(:,:) = Tsoil(1,:,:)

! Initialize state variables from a named start file
if (start_file /= 'none') then
  open(ustr, file = start_file)
  read(ustr,*) albs(:,:)
  read(ustr,*) Ds(:,:,:)
  read(ustr,*) Nsnow(:,:)
  read(ustr,*) Qcan(:,:)
  read(ustr,*) Sice(:,:,:)
  read(ustr,*) Sliq(:,:,:)
  read(ustr,*) Sveg(:,:)
  read(ustr,*) Tcan(:,:)
  read(ustr,*) theta(:,:,:)
  read(ustr,*) Tsnow(:,:,:)
  read(ustr,*) Tsoil(:,:,:)
  read(ustr,*) Tsurf(:,:)
  read(ustr,*) Tveg(:,:)
  close(ustr)
end if

! Outputs
allocate(diags(Nx,Ny,Ndiags))
allocate(SWin(Nx,Ny))
allocate(SWout(Nx,Ny))
diags(:,:,:) = 0
SWin(:,:) = 0
SWout(:,:) = 0
Nave = 24
Nsmp = 12
ave_file = 'ave'
dmp_file = 'dump'
smp_file = 'smp'
runid = 'none'
read(5, outputs)
if (runid == 'none') then
  runid = ''
endif
open(uave, file = trim(runid) // '_' // trim(ave_file))
open(udmp, file = trim(runid) // '_' // trim(dmp_file))
open(usmp, file = trim(runid) // '_' // trim(smp_file))

! Write options and namelists to metadata file
open(umta, file = trim(runid) // '_' // 'runinfo')
write(umta,*) '##################################'
write(umta,*) '#                                #'
write(umta,*) '# Flexible Snow Model FSM 2.0    #'
write(umta,*) '#                                #'
write(umta,*) '##################################'
write(umta,*) 
call date_and_time(VALUES=time)
write(umta,100) time(5),time(6),time(3),time(2),time(1)
100 format (' Run at ',i2,':',i2,1x,i2,'/',i2,'/',i4)
write(umta,*) 
write(umta,*) 'PHYSICS OPTIONS'
write(umta,*) 'snow albedo        ', ALBEDO_OPT
write(umta,*) 'canopy model       ', CANMOD_OPT
write(umta,*) 'snow conductivity  ', CONDCT_OPT
write(umta,*) 'snow density       ', DENSTY_OPT
write(umta,*) 'turbulent exchange ', EXCHNG_OPT
write(umta,*) 'snow hydraulics    ', HYDROL_OPT
write(umta,*)
write(umta,*) 'NAMELISTS'
write(umta,drive)
write(umta,gridpnts)
write(umta,gridlevs)
write(umta,maps)
write(umta,params)
write(umta,initial)
write(umta,outputs)

end subroutine SETUP
