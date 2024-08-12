!----------------------------------------------------------------------!
! Flexible Snow Model (FSM version 2.1.0)                              !
!                                                                      !
! Richard Essery                                                       !
! School of GeoSciences                                                !
! University of Edinburgh                                              !
!----------------------------------------------------------------------!
program FSM2

#include "OPTS.h"

#if PROFNC == 1
use netcdf
#endif

use CONSTANTS, only: &
  hcon_air,          &! Thermal conductivity of air (W/m/K)
  hcon_clay,         &! Thermal conductivity of clay (W/m/K)
  hcon_sand           ! Thermal conductivity of sand (W/m/K)

use IOUNITS, only: &
  ucan,              &! Subcanopy diagnostics file unit number
  udmp,              &! Start / dump file unit number
  uflx,              &! Flux output file unit number
  umet,              &! Meteorological driving file unit number
  usta                ! State output file unit number

use LAYERS, only: &
  Dzsnow,            &! Minimum snow layer thicknesses (m)
  Dzsoil,            &! Soil layer thicknesses (m)
  fvg1,              &! Fraction of vegetation in upper canopy layer
  Ncnpy,             &! Number of canopy layers
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  zsub                ! Subcanopy wind speed diagnostic height (m)

use PARAMETERS, only: &
  fcly,              &! Soil clay fraction
  fsnd,              &! Soil sand fraction
  rgr0                ! Fresh snow grain radius (m)

use SOILPROPS, only: &
  b,                 &! Clapp-Hornberger exponent
  hcap_soil,         &! Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil,         &! Thermal conductivity of dry soil (W/m/K)
  sathh,             &! Saturated soil water pressure (m)
  Vcrit,             &! Volumetric soil moisture at critical point
  Vsat                ! Volumetric soil moisture at saturation

implicit none

! Grid dimensions
integer :: &
  n,                 &! Layer/point counter
  Npnts               ! Number of points

! Site characteristics
real :: &
  lat,               &! Latitude (radians)
  noon                ! Time of solar noon (hour)

! Meteorological driving data
character(len=70) :: &
  met_file            ! Meteorological driving file name
integer :: &
  year,              &! Year
  month,             &! Month of year
  day                 ! Day of month
logical :: EoF        ! End-of-file flag
real :: &
  dt,                &! Timestep (s)
  elev,              &! Solar elevation (radians)
  hour,              &! Hour of day
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind speed measurement height (m)
real :: &
  LW,                &! Incoming longwave radiation (W/m^2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Rf,                &! Rainfall rate (kg/m^2/s)
  Sdif,              &! Diffuse shortwave radiation (W/m^2)
  Sdir,              &! Direct-beam shortwave radiation (W/m^2)
  Sf,                &! Snowfall rate (kg/m^2/s)
  Ta,                &! Air temperature (K)
  trans,             &! Wind-blown snow transport rate (kg/m^2/s)
  Ua                  ! Wind speed (m/s)

! Model state variables  
integer, allocatable :: &
  Nsnow(:)            ! Number of snow layers
real, allocatable :: &
  albs(:),           &! Snow albedo
  Tsrf(:),           &! Snow/ground surface temperature (K)
  Dsnw(:,:),         &! Snow layer thicknesses (m)
  Qcan(:,:),         &! Canopy air space humidities
  Rgrn(:,:),         &! Snow layer grain radii (m)
  Sice(:,:),         &! Ice content of snow layers (kg/m^2)
  Sliq(:,:),         &! Liquid content of snow layers (kg/m^2)
  Sveg(:,:),         &! Snow mass on vegetation layers (kg/m^2)
  Tcan(:,:),         &! Canopy air space temperatures (K)
  Tsnow(:,:),        &! Snow layer temperatures (K)
  Tsoil(:,:),        &! Soil layer temperatures (K)
  Tveg(:,:),         &! Vegetation layer temperatures (K)
  Vsmc(:,:),         &! Volumetric moisture content of soil layers
  fsat(:),           &! Initial soil layer moisture/saturation
  Tprf(:)             ! Initial soil layer temperatures (K)

! Diagnostics
real, allocatable :: &
  fsnow(:),          &! Ground snowcover fraction
  H(:),              &! Sensible heat flux to the atmosphere (W/m^2)
  LE(:),             &! Latent heat flux to the atmosphere (W/m^2)
  LWout(:),          &! Outgoing LW radiation (W/m^2)
  LWsub(:),          &! Subcanopy downward LW radiation (W/m^2)
  Melt(:),           &! Surface melt rate (kg/m^2/s)
  Roff(:),           &! Runoff from snow (kg/m^2/s)
  snd(:),            &! Snow depth (m)
  snw(:),            &! Total snow mass on ground (kg/m^2) 
  subl(:),           &! Sublimation rate (kg/m^2/s)
  svg(:),            &! Total snow mass on vegetation (kg/m^2)
  SWout(:),          &! Outgoing SW radiation (W/m^2)
  SWsub(:),          &! Subcanopy downward SW radiation (W/m^2)
  Tsub(:),           &! Subcanopy air temperature (K)
  Usub(:),           &! Subcanopy wind speed (m/s)
  Wflx(:,:)           ! Water flux into snow layer (kg/m^2/s)

! Vegetation characteristics
character(len=70) :: &
  alb0_file,         &! Snow-free ground albedo map file name
  vegh_file,         &! Canopy height map file name
  VAI_file            ! Vegetation area index map file name
real, allocatable :: &
  alb0(:),           &! Snow-free ground albedo
  vegh(:),           &! Canopy height (m)
  VAI(:)              ! Vegetation area index

! Start and dump file names
character(len=70) :: &
  dump_file,         &! End dump file name
  runid,             &! Run identifier
  start_file          ! Start file name

! NetCDF variables
integer :: &
  ncid,              &! Dataset ID
  rec,               &! Record number
  status,            &! Error status
  varid(17)           ! Variable IDs

namelist    /drive/ met_file,dt,lat,noon,zT,zU
namelist /gridpnts/ Npnts,Nsmax,Nsoil
namelist /gridlevs/ Dzsnow,Dzsoil,fvg1,zsub
namelist  /initial/ fsat,Tprf,start_file
namelist  /outputs/ dump_file,runid
namelist      /veg/ alb0,vegh,VAI,alb0_file,vegh_file,VAI_file
                    
call FSM2_PARAMS

! Grid dimensions
#if CANMOD == 1
Ncnpy = 1
#endif
#if CANMOD == 2
Ncnpy = 2
#endif
Npnts = 1
Nsmax = 3
Nsoil = 4
read(5,gridpnts)

! Canopy, snow and soil layers
fvg1 = 0.5
zsub = 1.5
allocate(Dzsnow(Nsmax))
allocate(Dzsoil(Nsoil))
if (Nsmax == 3) Dzsnow = (/0.1, 0.2, 0.4/)
if (Nsoil == 4) Dzsoil = (/0.1, 0.2, 0.4, 0.8/)
read(5, gridlevs)

! Site and driving data characteristics
met_file = 'met'
dt = 3600
lat = 0
noon = 12
zT = 2
zU = 10
read(5,drive)
open(umet, file = met_file)
read(umet,*) year,month,day,hour
rewind(umet)
lat = (3.14159/180)*lat
trans = 0  ! no snow transport

! Vegetation characteristics from defaults, namelist or named files
allocate(alb0(Npnts))
allocate(vegh(Npnts))
allocate(VAI(Npnts))
alb0_file = 'none'
vegh_file = 'none'
VAI_file  = 'none'
alb0 = 0.2
vegh = 0
VAI  = 0
read(5,veg)
if (alb0_file /= 'none') call FSM2_VEG(alb0_file,Npnts,alb0)
if (vegh_file /= 'none') call FSM2_VEG(vegh_file,Npnts,vegh)
if (VAI_file  /= 'none') call FSM2_VEG(VAI_file,Npnts,VAI)

! Soil properties
b = 3.1 + 15.7*fcly - 0.3*fsnd
hcap_soil = (2.128*fcly + 2.385*fsnd)*1e6 / (fcly + fsnd)
sathh = 10**(0.17 - 0.63*fcly - 1.58*fsnd)
Vsat = 0.505 - 0.037*fcly - 0.142*fsnd
Vcrit = Vsat*(sathh/3.364)**(1/b)
hcon_soil = (hcon_air**Vsat) * ((hcon_clay**fcly)*(hcon_sand**(1 - fcly))**(1 - Vsat))

! Allocate state variable arrays
allocate(albs(Npnts))
allocate(Nsnow(Npnts))
allocate(Tsrf(Npnts))
allocate(Dsnw(Nsmax,Npnts))
allocate(Qcan(Ncnpy,Npnts))
allocate(Rgrn(Nsmax,Npnts))
allocate(Sice(Nsmax,Npnts))
allocate(Sliq(Nsmax,Npnts))
allocate(Sveg(Ncnpy,Npnts))
allocate(Tcan(Ncnpy,Npnts))
allocate(Tsnow(Nsmax,Npnts))
allocate(Tsoil(Nsoil,Npnts))
allocate(Tveg(Ncnpy,Npnts))
allocate(Vsmc(Nsoil,Npnts))

! Default initialization of state variables
albs  = 0.8
Dsnw  = 0
Nsnow = 0
Qcan  = 0
Rgrn  = rgr0
Sice  = 0
Sliq  = 0
Sveg  = 0
Tcan  = 285
Tsnow = 273
Tveg  = 285
! Missing values for vegetation at non-forest points
do n = 1, Ncnpy
  where(VAI==0) Sveg(n,:) = -999./Ncnpy
  where(VAI==0) Tveg(n,:) = -999
end do

! Initial soil profiles from namelist
allocate(fsat(Nsoil))
allocate(Tprf(Nsoil))
fsat = 0.5
Tprf = 285
start_file = 'none'
read(5,initial)
do n = 1, Nsoil
  Tsoil(n,:) = Tprf(n)
  Vsmc(n,:) = fsat(n)*Vsat
end do
Tsrf = Tsoil(1,:)

! Initialize state variables from a named start file
if (start_file /= 'none') then
  open(udmp,file = start_file)
  read(udmp,*) albs
  read(udmp,*) Dsnw
  read(udmp,*) Nsnow
  read(udmp,*) Qcan
  read(udmp,*) Rgrn
  read(udmp,*) Sice
  read(udmp,*) Sliq
  read(udmp,*) Sveg
  read(udmp,*) Tcan
  read(udmp,*) Tsnow
  read(udmp,*) Tsoil
  read(udmp,*) Tsrf
  read(udmp,*) Tveg
  read(udmp,*) Vsmc
  close(udmp)
end if

! Allocate diagnostic output arrays
allocate(fsnow(Npnts))
allocate(H(Npnts))
allocate(LE(Npnts))
allocate(LWout(Npnts))
allocate(LWsub(Npnts))
allocate(Melt(Npnts))
allocate(Roff(Npnts))
allocate(snd(Npnts))
allocate(snw(Npnts))
allocate(subl(Npnts))
allocate(svg(Npnts))
allocate(SWout(Npnts))
allocate(SWsub(Npnts))
allocate(Tsub(Npnts))
allocate(Usub(Npnts))
allocate(Wflx(Nsmax,Npnts))

! Output files
dump_file = 'dump'
runid = 'none'
read(5,outputs)
if (runid == 'none') runid = ''
#if PROFNC == 1
if (Npnts>1) stop 'NetCDF output only available for Npnts = 1'
call FSM2_PREPNC(runid,year,month,day,hour,ncid,rec,varid)
#else
open(ucan, file = trim(runid)//'subc.txt')
open(uflx, file = trim(runid)//'flux.txt')
open(usta, file = trim(runid)//'stat.txt')
#endif

! Run the model
EoF = .false.
do
  call FSM2_DRIVE(lat,noon,year,month,day,hour,elev,EoF,               &
                  LW,Ps,Qa,Rf,Sdif,Sdir,Sf,Ta,Ua)  
  if (EoF) goto 1
  do n = 1, Npnts
    call FSM2_TIMESTEP(                                                &
                       ! Driving variables                             &
                       dt,elev,zT,zU,LW,Ps,Qa,Rf,Sdif,Sdir,Sf,Ta,      &
                       trans,Ua,                                       &
                       ! Vegetation characteristics                    &
                       alb0(n),vegh(n),VAI(n),                         &
                       ! State variables                               &
                       albs(n),Tsrf(n),Dsnw(:,n),Nsnow(n), Qcan(:,n),  &
                       Rgrn(:,n),Sice(:,n),Sliq(:,n),Sveg(:,n),        &
                       Tcan(:,n),Tsnow(:,n),Tsoil(:,n),Tveg(:,n),      &
                       Vsmc(:,n),                                      &
                       ! Diagnostics                                   &
                       fsnow(n),H(n),LE(n),LWout(n),LWsub(n),Melt(n),  &
                       Roff(n),snd(n),snw(n),subl(n),svg(n),SWout(n),  &
                       SWsub(n),Tsub(n),Usub(n),Wflx(:,n)              )
  end do
#if PROFNC == 1
  call FSM2_WRITENC(Dsnw(:,1),dt,H,LE,LWout,Melt,ncid,Nsnow,           &
                    Rgrn(:,1),Roff,Sice(:,1),Sliq(:,1),snd,snw,SWout,  &
                    Tsnow(:,1),Tsoil(:,1),Tsrf,varid,Wflx(:,1),rec) 
#else
  call FSM2_OUTPUT(Npnts,year,month,day,hour,                          &
                   H,LE,LWout,LWsub,Melt,Roff,snd,snw,subl,svg,SWout,  &
                   SWsub,Tsoil,Tsrf,Tsub,Tveg,Usub)
#endif
end do
1 continue

! Write out state variables at end of run
open(udmp,file = trim(runid) // trim(dump_file))
write(udmp,*) albs
write(udmp,*) Dsnw
write(udmp,*) Nsnow
write(udmp,*) Qcan
write(udmp,*) Rgrn
write(udmp,*) Sice
write(udmp,*) Sliq
write(udmp,*) Sveg
write(udmp,*) Tcan
write(udmp,*) Tsnow
write(udmp,*) Tsoil
write(udmp,*) Tsrf
write(udmp,*) Tveg
write(udmp,*) Vsmc
close(udmp)

#if PROFNC == 1
status = nf90_close(ncid) 
#endif

end program FSM2
