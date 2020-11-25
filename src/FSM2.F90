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
  Ncols,             &! Number of columns in grid
  Nrows               ! Number of rows in grid 

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
real, allocatable :: &
  LW(:,:),           &! Incoming longwave radiation (W/m^2)
  Ps(:,:),           &! Surface pressure (Pa)
  Qa(:,:),           &! Specific humidity (kg/kg)
  Rf(:,:),           &! Rainfall rate (kg/m^2/s)
  Sdif(:,:),         &! Diffuse shortwave radiation (W/m^2)
  Sdir(:,:),         &! Direct-beam shortwave radiation (W/m^2)
  Sf(:,:),           &! Snowfall rate (kg/m^2/s)
  Ta(:,:),           &! Air temperature (K)
  trans(:,:),        &! Wind-blown snow transport rate (kg/m^2/s)
  Ua(:,:)             ! Wind speed (m/s)

! Model state variables  
integer, allocatable :: &
  Nsnow(:,:)          ! Number of snow layers
real, allocatable :: &
  albs(:,:),         &! Snow albedo
  Tsrf(:,:),         &! Snow/ground surface temperature (K)
  Dsnw(:,:,:),       &! Snow layer thicknesses (m)
  Qcan(:,:,:),       &! Canopy air space humidities
  Rgrn(:,:,:),       &! Snow layer grain radii (m)
  Sice(:,:,:),       &! Ice content of snow layers (kg/m^2)
  Sliq(:,:,:),       &! Liquid content of snow layers (kg/m^2)
  Sveg(:,:,:),       &! Snow mass on vegetation layers (kg/m^2)
  Tcan(:,:,:),       &! Canopy air space temperatures (K)
  Tsnow(:,:,:),      &! Snow layer temperatures (K)
  Tsoil(:,:,:),      &! Soil layer temperatures (K)
  Tveg(:,:,:),       &! Vegetation layer temperatures (K)
  Vsmc(:,:,:),       &! Volumetric moisture content of soil layers
  fsat(:),           &! Initial soil layer moisture/saturation
  Tprf(:)             ! Initial soil layer temperatures (K)

! Diagnostics
real, allocatable :: &
  H(:,:),            &! Sensible heat flux to the atmosphere (W/m^2)
  LE(:,:),           &! Latent heat flux to the atmosphere (W/m^2)
  LWout(:,:),        &! Outgoing LW radiation (W/m^2)
  LWsub(:,:),        &! Subcanopy downward LW radiation (W/m^2)
  Melt(:,:),         &! Surface melt rate (kg/m^2/s)
  Roff(:,:),         &! Runoff from snow (kg/m^2/s)
  snd(:,:),          &! Snow depth (m)
  snw(:,:),          &! Total snow mass on ground (kg/m^2) 
  subl(:,:),         &! Sublimation rate (kg/m^2/s)
  svg(:,:),          &! Total snow mass on vegetation (kg/m^2)
  SWout(:,:),        &! Outgoing SW radiation (W/m^2)
  SWsub(:,:),        &! Subcanopy downward SW radiation (W/m^2)
  Usub(:,:),         &! Subcanopy wind speed (m/s)
  Wflx(:,:,:)         ! Water flux into snow layer (kg/m^2/s)

! Vegetation characteristics
character(len=70) :: &
  alb0_file,         &! Snow-free ground albedo map file name
  fsky_file,         &! Skyview fraction map file name
  vegh_file,         &! Canopy height map file name
  VAI_file            ! Vegetation area index map file name
real, allocatable :: &
  alb0(:,:),         &! Snow-free ground albedo
  fsky(:,:),         &! Skyview not obstructed by remote vegetation
  vegh(:,:),         &! Canopy height (m)
  VAI(:,:)            ! Vegetation area index

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

! Counters
integer :: &
  i,                 &! Grid row counter
  j,                 &! Grid column counter
  k                   ! Soil layer counter

namelist    /drive/ met_file,dt,lat,noon,zT,zU
namelist /gridpnts/ Ncols,Nrows,Nsmax,Nsoil
namelist /gridlevs/ Dzsnow,Dzsoil,fvg1,zsub
namelist  /initial/ fsat,Tprf,start_file
namelist  /outputs/ dump_file,runid
namelist      /veg/ alb0,fsky,vegh,VAI,  &
                    alb0_file,fsky_file,vegh_file,VAI_file

#if SETPAR == 1
call FSM2_PARAMS
#endif

! Grid dimensions
Ncols = 1
Nrows = 1
Nsmax = 3
Nsoil = 4
read(5,gridpnts)

! Canopy, snow and soil layers
#if CANMOD == 1
Ncnpy = 1
#endif
#if CANMOD == 2
Ncnpy = 2
#endif
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

! Allocate driving data arrays
allocate(LW(Ncols,Nrows))
allocate(Ps(Ncols,Nrows))
allocate(Qa(Ncols,Nrows))
allocate(Rf(Ncols,Nrows))
allocate(Sdif(Ncols,Nrows))
allocate(Sdir(Ncols,Nrows))
allocate(Sf(Ncols,Nrows))
allocate(Ta(Ncols,Nrows))
allocate(trans(Ncols,Nrows))
allocate(Ua(Ncols,Nrows))
trans(:,:) = 0

! Vegetation characteristics from defaults, namelist or named map files
allocate(alb0(Ncols,Nrows))
allocate(fsky(Ncols,Nrows))
allocate(vegh(Ncols,Nrows))
allocate(VAI(Ncols,Nrows))
alb0_file = 'none'
fsky_file = 'none'
vegh_file = 'none'
VAI_file  = 'none'
alb0 = 0.2
fsky = 1
vegh = 0
VAI  = 0
read(5,veg)
if (alb0_file /= 'none') call FSM2_MAP(alb0_file,Ncols,Nrows,alb0)
if (fsky_file /= 'none') call FSM2_MAP(fsky_file,Ncols,Nrows,fsky)
if (vegh_file /= 'none') call FSM2_MAP(vegh_file,Ncols,Nrows,vegh)
if (VAI_file  /= 'none') call FSM2_MAP(VAI_file,Ncols,Nrows,VAI)

! Soil properties
b = 3.1 + 15.7*fcly - 0.3*fsnd
hcap_soil = (2.128*fcly + 2.385*fsnd)*1e6 / (fcly + fsnd)
sathh = 10**(0.17 - 0.63*fcly - 1.58*fsnd)
Vsat = 0.505 - 0.037*fcly - 0.142*fsnd
Vcrit = Vsat*(sathh/3.364)**(1/b)
hcon_soil = (hcon_air**Vsat) * ((hcon_clay**fcly)*(hcon_sand**(1 - fcly))**(1 - Vsat))

! Allocate state variable arrays
allocate(albs(Ncols,Nrows))
allocate(Nsnow(Ncols,Nrows))
allocate(Tsrf(Ncols,Nrows))
allocate(Dsnw(Nsmax,Ncols,Nrows))
allocate(Qcan(Ncnpy,Ncols,Nrows))
allocate(Rgrn(Nsmax,Ncols,Nrows))
allocate(Sice(Nsmax,Ncols,Nrows))
allocate(Sliq(Nsmax,Ncols,Nrows))
allocate(Sveg(Ncnpy,Ncols,Nrows))
allocate(Tcan(Ncnpy,Ncols,Nrows))
allocate(Tsnow(Nsmax,Ncols,Nrows))
allocate(Tsoil(Nsoil,Ncols,Nrows))
allocate(Tveg(Ncnpy,Ncols,Nrows))
allocate(Vsmc(Nsoil,Ncols,Nrows))

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
do k = 1, Ncnpy
  where(VAI==0) Sveg(k,:,:) = -999./Ncnpy
  where(VAI==0) Tveg(k,:,:) = -999
end do

! Initial soil profiles from namelist
allocate(fsat(Nsoil))
allocate(Tprf(Nsoil))
fsat = 0.5
Tprf = 285
start_file = 'none'
read(5,initial)
do k = 1, Nsoil
  Tsoil(k,:,:) = Tprf(k)
  Vsmc(k,:,:) = fsat(k)*Vsat
end do
Tsrf = Tsoil(1,:,:)

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
allocate(H(Ncols,Nrows))
allocate(LE(Ncols,Nrows))
allocate(LWout(Ncols,Nrows))
allocate(LWsub(Ncols,Nrows))
allocate(Melt(Ncols,Nrows))
allocate(Roff(Ncols,Nrows))
allocate(snd(Ncols,Nrows))
allocate(snw(Ncols,Nrows))
allocate(subl(Ncols,Nrows))
allocate(svg(Ncols,Nrows))
allocate(SWout(Ncols,Nrows))
allocate(SWsub(Ncols,Nrows))
allocate(Usub(Ncols,Nrows))
allocate(Wflx(Nsmax,Ncols,Nrows))

! Output files
dump_file = 'dump'
runid = 'none'
read(5,outputs)
if (runid == 'none') runid = ''
#if PROFNC == 1
if (Ncols*Nrows>1) stop 'NetCDF output only available for Nrows = Ncols = 1'
call FSM2_PREPNC(runid,year,month,day,hour,ncid,rec,varid)
#else
if (maxval(VAI) > 0) open(ucan, file = trim(runid)//'subc.txt')
open(uflx, file = trim(runid)//'flux.txt')
open(usta, file = trim(runid)//'stat.txt')
#endif

! Run the model
EoF = .false.
do
  call FSM2_DRIVE(Ncols,Nrows,fsky,lat,noon,                           &
                  year,month,day,hour,elev,EoF,                        &
                  LW,Ps,Qa,Rf,Sdif,Sdir,Sf,Ta,Ua)  
  if (EoF) goto 1
  do i = 1, Nrows
  do j = 1, Ncols
    call FSM2_TIMESTEP(                                                &
                       ! Driving variables                             &
                       dt,elev,zT,zU,LW(j,i),Ps(j,i),Qa(j,i),          &
                       Rf(j,i),Sdif(j,i),Sdir(j,i),Sf(j,i),            &
                       Ta(j,i),trans(j,i),Ua(j,i),                     &
                       ! Vegetation characteristics                    &
                       alb0(j,i),vegh(j,i),VAI(j,i),                   &
                       ! State variables                               &
                       albs(j,i),Tsrf(j,i),Dsnw(:,j,i),Nsnow(j,i),     &
                       Qcan(:,j,i),Rgrn(:,j,i),Sice(:,j,i),            &
                       Sliq(:,j,i),Sveg(:,j,i),Tcan(:,j,i),            &
                       Tsnow(:,j,i),Tsoil(:,j,i),Tveg(:,j,i),          &
                       Vsmc(:,j,i),                                    &
                       ! Diagnostics                                   &
                       H(j,i),LE(j,i),LWout(j,i),LWsub(j,i),           &
                       Melt(j,i),Roff(j,i),snd(j,i),snw(j,i),          &
                       subl(j,i),svg(j,i),SWout(j,i),SWsub(j,i),       &
                       Usub(j,i),Wflx(:,j,i)                           )
  end do
  end do
#if PROFNC == 1
  call FSM2_WRITENC(Dsnw(:,1,1),dt,H(1,1),LE(1,1),LWout(1,1),Melt(1,1),&
                    ncid,Nsnow(1,1),Rgrn(:,1,1),Roff(1,1),Sice(:,1,1), &
                    Sliq(:,1,1),snd(1,1),snw(1,1),SWout(1,1),          &
                    Tsnow(:,1,1),Tsoil(:,1,1),Tsrf(1,1),varid,         &
                    Wflx(:,1,1),rec) 
#else
  call FSM2_OUTPUT(Ncols,Nrows,year,month,day,hour,                    &
                   H,LE,LWout,LWsub,Melt,Roff,snd,snw,subl,svg,SWout,  &
                   SWsub,Tsoil,Tsrf,Tveg,Usub,VAI)
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
