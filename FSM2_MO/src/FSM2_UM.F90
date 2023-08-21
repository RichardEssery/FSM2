!----------------------------------------------------------------------!
! Drive FSM2 with UM driving data                                      !
!                                                                      !
! Richard Essery                                                       !
! School of GeoSciences                                                !
! University of Edinburgh                                              !
!----------------------------------------------------------------------!
program FSM2_GLOBAL

#include "OPTS.h"

use netcdf

use CONSTANTS, only: &
  hcon_air,          &! Thermal conductivity of air (W/m/K)
  hcon_clay,         &! Thermal conductivity of clay (W/m/K)
  hcon_sand           ! Thermal conductivity of sand (W/m/K)

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

! Meteorological driving data
real :: &
  dt,                &! Timestep (s)
  elev,              &! Solar elevation (radians)
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind speed measurement height (m)
real, allocatable :: &
  LW(:,:),           &! Incoming longwave radiation (W/m^2)
  Ps(:,:),           &! Surface pressure (Pa)
  Qa(:,:),           &! Specific humidity (kg/kg)
  Rf(:,:),           &! Rainfall rate (kg/m^2/s)
  RH(:,:),           &! Relative humidity (%)
  Sdif(:,:),         &! Diffuse shortwave radiation (W/m^2)
  Sdir(:,:),         &! Direct-beam shortwave radiation (W/m^2)
  Sf(:,:),           &! Snowfall rate (kg/m^2/s)
  Ta(:,:),           &! Air temperature (K)
  trans(:),          &! Wind-blown snow transport rate (kg/m^2/s)
  Ua(:,:)             ! Wind speed (m/s)

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
  Vsmc(:,:)           ! Volumetric moisture content of soil layers
logical :: &
  start_file          ! True if start file exists

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
  SWup(:),           &! Cumulated outgoing SW radiation (J/m^2)
  Tsub(:),           &! Subcanopy air temperature (K)
  Usub(:),           &! Subcanopy wind speed (m/s)
  Wflx(:,:)           ! Water flux into snow layer (kg/m^2/s)

! Surface characteristics
real, allocatable :: &
  alb0(:),           &! Snow-free ground albedo
  vegh(:),           &! Canopy height (m)
  VAI(:)              ! Vegetation area index

! NetCDF variables
character(80) :: &
  timeunits           ! Time units
integer :: &
  sdimid,            &! Snow layer dimension ID
  tdimid,            &! Time dimension ID
  xdimid,            &! x dimension ID
  ydimid,            &! y dimension ID
  hour,              &! Time (hours)
  ncid,              &! Dataset ID
  Nout,              &! Number of timesteps between outputs
  Ntime,             &! Number of timesteps
  rec,               &! Record number
  status,            &! Error status
  varid,varids(14)    ! Variable IDs

! Simulation domain
integer :: &
  Nland               ! Number of land points
integer, allocatable :: &
  x(:),              &! x indices of land points
  y(:)                ! y indices of land points
!integer, parameter :: &
!  Nx = 800,          &! Number of x points in simulation domain
!  Ny = 800,          &! Number of y points in simulation domain
!  x0 = 200,          &! First x point in simulation domain
!  y0 = 200            ! First y point in simulation domain
 
integer, parameter :: &
  Nx = 100,          &! Number of x points in simulation domain
  Ny = 200,          &! Number of y points in simulation domain
  x0 = 820,          &! First x point in simulation domain
  y0 = 780            ! First y point in simulation domain 
  
real :: &
  lats(Nx,Ny),       &! Latitudes (degrees) 
  lons(Nx,Ny),       &! Longitudes (degrees) 
  var2d(Nx,Ny)        ! Variable on 2D grid
real, allocatable :: &
  var3d(:,:,:)        ! Variable on 3D grid

! Counters
integer :: &
  i,                 &! x point counter
  j,                 &! y point counter
  l,                 &! Land point counter
  n                   ! Timestep counter
  
call FSM2_PARAMS

! Grid dimensions
#if CANMOD == 1
Ncnpy = 1
#endif
#if CANMOD == 2
Ncnpy = 2
#endif
Nsmax = 3
Nsoil = 4

! Canopy, snow and soil layers
fvg1 = 0.5
zsub = 1.5
allocate(Dzsnow(Nsmax))
allocate(Dzsoil(Nsoil))
Dzsnow = (/0.1, 0.2, 0.4/)
Dzsoil = (/0.1, 0.2, 0.4, 0.8/)
allocate(var3d(Nx,Ny,Nsmax))

! Land mask and number of points
status = nf90_open('ancils.nc',nf90_nowrite,ncid)
status = nf90_inq_varid(ncid,'lat',varid)
status = nf90_get_var(ncid,varid,lats,start=(/x0,y0/),count=(/Nx,Ny/))
status = nf90_inq_varid(ncid,'lon',varid)
status = nf90_get_var(ncid,varid,lons,start=(/x0,y0/),count=(/Nx,Ny/))
status = nf90_inq_varid(ncid,'landmask',varid)
status = nf90_get_var(ncid,varid,var2d,start=(/x0,y0/),count=(/Nx,Ny/))
Nland = count(var2d == 1)
allocate(y(Nland))
allocate(x(Nland))
l = 0
do i = 1, Nx
do j = 1, Ny
  if (var2d(i,j) == 1) then
    l = l + 1
    x(l) = i
    y(l) = j
  end if
end do
end do

! Surface characteristics
allocate(alb0(Nland))
allocate(vegh(Nland))
allocate(VAI(Nland))
status = nf90_inq_varid(ncid,'alb0',varid)
status = nf90_get_var(ncid,varid,var2d,start=(/x0,y0/),count=(/Nx,Ny/))
do l = 1, Nland
  alb0(l) = var2d(x(l),y(l))
end do
status = nf90_inq_varid(ncid,'vegh',varid)
status = nf90_get_var(ncid,varid,var2d,start=(/x0,y0/),count=(/Nx,Ny/))
do l = 1, Nland
  vegh(l) = var2d(x(l),y(l))
end do
status = nf90_inq_varid(ncid,'VAI',varid)
status = nf90_get_var(ncid,varid,var2d,start=(/x0,y0/),count=(/Nx,Ny/))
do l = 1, Nland
  VAI(l) = var2d(x(l),y(l))
end do
status = nf90_close(ncid)

! Soil properties
b = 3.1 + 15.7*fcly - 0.3*fsnd
hcap_soil = (2.128*fcly + 2.385*fsnd)*1e6 / (fcly + fsnd)
sathh = 10**(0.17 - 0.63*fcly - 1.58*fsnd)
Vsat = 0.505 - 0.037*fcly - 0.142*fsnd
Vcrit = Vsat*(sathh/3.364)**(1/b)
hcon_soil = (hcon_air**Vsat) * ((hcon_clay**fcly)*(hcon_sand**(1 - fcly))**(1 - Vsat))

! Driving data characteristics from an input file
dt = 3600
Nout = 24
zT = 2
zU = 10
status = nf90_open('LWdown.nc', nf90_nowrite, ncid)
status = nf90_inq_dimid(ncid,'time',tdimid)
status = nf90_inquire_dimension(ncid,tdimid,len=Ntime)
status = nf90_inq_varid(ncid,'time',varid)
status = nf90_get_att(ncid,varid,'units',timeunits)
status = nf90_get_var(ncid,varid,hour)
status = nf90_close(ncid)

! Allocate driving data arrays
allocate(LW(Nland,Ntime))
allocate(Ps(Nland,Ntime))
allocate(Qa(Nland,Ntime))
allocate(Rf(Nland,Ntime))
allocate(RH(Nland,Ntime))
allocate(Sf(Nland,Ntime))
allocate(Sdif(Nland,Ntime))
allocate(Sdir(Nland,Ntime))
allocate(Ta(Nland,Ntime))
allocate(trans(Nland))
allocate(Ua(Nland,Ntime))
where(Ua < 0.1) Ua = 0.1
trans(:) = 0

! Read driving data
call READNC_LAND(Nland,Ntime,x0,y0,Nx,Ny,x,y,'LWdown',LW)
call READNC_LAND(Nland,Ntime,x0,y0,Nx,Ny,x,y,'Psurf',Ps)
call READNC_LAND(Nland,Ntime,x0,y0,Nx,Ny,x,y,'RelHum',RH)
call READNC_LAND(Nland,Ntime,x0,y0,Nx,Ny,x,y,'Rainf',Rf)
call READNC_LAND(Nland,Ntime,x0,y0,Nx,Ny,x,y,'Snowf',Sf)
call READNC_LAND(Nland,Ntime,x0,y0,Nx,Ny,x,y,'SWdown',Sdif)
call READNC_LAND(Nland,Ntime,x0,y0,Nx,Ny,x,y,'Tair',Ta)
call READNC_LAND(Nland,Ntime,x0,y0,Nx,Ny,x,y,'Wind',Ua)
where(RH>100) RH = 100 
Qa = (RH/100)*0.622*(611.2/Ps)*exp(17.5043*(Ta - 273.15)/(Ta - 31.85))
Sdir(:,:) = 0
elev = 0

! Allocate state variable arrays
allocate(albs(Nland))
allocate(Nsnow(Nland))
allocate(Tsrf(Nland))
allocate(Dsnw(Nsmax,Nland))
allocate(Qcan(Ncnpy,Nland))
allocate(Rgrn(Nsmax,Nland))
allocate(Sice(Nsmax,Nland))
allocate(Sliq(Nsmax,Nland))
allocate(Sveg(Ncnpy,Nland))
allocate(Tcan(Ncnpy,Nland))
allocate(Tsnow(Nsmax,Nland))
allocate(Tsoil(Nsoil,Nland))
allocate(Tveg(Ncnpy,Nland))
allocate(Vsmc(Nsoil,Nland))

! Initialize state variables
inquire(file='start.txt',exist=start_file)
if (start_file) then  ! Initialize from start file if it exists
  open(8,file='start.txt')
  read(8,*) albs
  read(8,*) Dsnw
  read(8,*) Nsnow
  read(8,*) Qcan
  read(8,*) Rgrn
  read(8,*) Sice
  read(8,*) Sliq
  read(8,*) Sveg
  read(8,*) Tcan
  read(8,*) Tsnow
  read(8,*) Tsoil
  read(8,*) Tsrf
  read(8,*) Tveg
  read(8,*) Vsmc
  close(8)
else                  ! Cold start
  albs = 0.8
  Dsnw = 0
  Nsnow = 0
  Qcan = 0
  Rgrn = rgr0
  Sice = 0
  Sliq = 0
  Sveg = 0
  Tcan = 285
  Tsnow = 273
  Tsoil = 285
  Tsrf = 285
  Tveg = 285
  Vsmc = 0.5*Vsat
end if

! Allocate diagnostic output arrays
allocate(fsnow(Nland))
allocate(H(Nland))
allocate(LE(Nland))
allocate(LWout(Nland))
allocate(LWsub(Nland))
allocate(Melt(Nland))
allocate(Roff(Nland))
allocate(snd(Nland))
allocate(snw(Nland))
allocate(subl(Nland))
allocate(svg(Nland))
allocate(SWout(Nland))
allocate(SWsub(Nland))
allocate(SWup(Nland))
allocate(Tsub(Nland))
allocate(Usub(Nland))
allocate(Wflx(Nsmax,Nland))

! Prepare NetCDF output file
status = nf90_create('FSM2out.nc',nf90_clobber,ncid)
status = nf90_def_dim(ncid,'x',Nx,xdimid)
status = nf90_def_dim(ncid,'y',Ny,ydimid)
status = nf90_def_dim(ncid,'snow_layer',Nsmax,sdimid)
status = nf90_def_dim(ncid,'time',NF90_UNLIMITED,tdimid)
status = nf90_def_var(ncid,'lon',NF90_REAL,(/xdimid,ydimid/),varids(1))
status = nf90_put_att(ncid,varids(1),'units','degrees_east')
status = nf90_def_var(ncid,'lat',NF90_REAL,(/xdimid,ydimid/),varids(2))
status = nf90_put_att(ncid,varids(2),'units','degrees_north')
status = nf90_def_var(ncid,'time',NF90_REAL,tdimid,varids(3))
status = nf90_put_att(ncid,varids(3),'units',timeunits)
status = nf90_def_var(ncid,'snc',NF90_REAL,(/xdimid,ydimid,tdimid/),varids(4))
status = nf90_put_att(ncid,varids(4),'long_name','snowcover fraction')
status = nf90_put_att(ncid,varids(4),'units','-')
status = nf90_def_var(ncid,'snd',NF90_REAL,(/xdimid,ydimid,tdimid/),varids(5))
status = nf90_put_att(ncid,varids(5),'long_name','snowdepth')
status = nf90_put_att(ncid,varids(5),'units','m')
status = nf90_def_var(ncid,'snw',NF90_REAL,(/xdimid,ydimid,tdimid/),varids(6))
status = nf90_put_att(ncid,varids(6),'long_name','mass of snowpack')
status = nf90_put_att(ncid,varids(6),'units','kg m-2')
status = nf90_def_var(ncid,'Nsnow',NF90_REAL,(/xdimid,ydimid,tdimid/),varids(7))
status = nf90_put_att(ncid,varids(7),'long_name','number of snow layers')
status = nf90_put_att(ncid,varids(7),'units','-')
status = nf90_def_var(ncid,'tsl',NF90_REAL,(/xdimid,ydimid,tdimid/),varids(8))
status = nf90_put_att(ncid,varids(8),'long_name','temperature of first soil layer')
status = nf90_put_att(ncid,varids(8),'units','K')
status = nf90_def_var(ncid,'SWup',NF90_REAL,(/xdimid,ydimid,tdimid/),varids(9))
status = nf90_put_att(ncid,varids(9),'long_name','upwelling SW radiation')
status = nf90_put_att(ncid,varids(9),'units','J m-2')
#if SMRTIO == 1
status = nf90_def_var(ncid,'Dsnw',NF90_REAL,(/xdimid,ydimid,sdimid,tdimid/),varids(10))
status = nf90_put_att(ncid,varids(10),'long_name','thickness of snow layers')
status = nf90_put_att(ncid,varids(10),'units','m')
status = nf90_def_var(ncid,'lqsn',NF90_REAL,(/xdimid,ydimid,sdimid,tdimid/),varids(11))
status = nf90_put_att(ncid,varids(11),'long_name','mass fraction of liquid water in snow layers')
status = nf90_put_att(ncid,varids(11),'units','-')
status = nf90_def_var(ncid,'snowrho',NF90_REAL,(/xdimid,ydimid,sdimid,tdimid/),varids(12))
status = nf90_put_att(ncid,varids(12),'long_name','density of snow layers')
status = nf90_put_att(ncid,varids(12),'units','kg m-3')
status = nf90_def_var(ncid,'rgrn',NF90_REAL,(/xdimid,ydimid,sdimid,tdimid/),varids(13))
status = nf90_put_att(ncid,varids(13),'long_name','grain radius in snow layers')
status = nf90_put_att(ncid,varids(13),'units','m')
status = nf90_def_var(ncid,'tsnl',NF90_REAL,(/xdimid,ydimid,sdimid,tdimid/),varids(14))
status = nf90_put_att(ncid,varids(14),'units','K')
status = nf90_put_att(ncid,varids(14),'long_name','temperature of snow layers')
#endif
do varid = 1, 14
  status = nf90_put_att(ncid,varids(varid),'missing_value',1.e+20)
end do
status = nf90_enddef(ncid)
status = nf90_put_var(ncid,varids(1),lons)
status = nf90_put_var(ncid,varids(2),lats) 

! Run the model and write daily output
rec = 1
snd(:) = 0
snw(:) = 0
SWup(:) = 0
do n = 1, Ntime
  write(6,fmt='(a,i3,a)',advance='no') achar(13),int(100*n/float(Ntime)),'% complete'
  do l = 1, Nland
    call FSM2_TIMESTEP(                                                &
                       ! Driving variables                             &
                       dt,elev,zT,zU,LW(l,n),Ps(l,n),Qa(l,n),Rf(l,n),  &
                       Sdif(l,n),Sdir(l,n),Sf(l,n),Ta(l,n),trans(l),   &
                       Ua(l,n),                                        &
                       ! Vegetation characteristics                    &
                       alb0(l),vegh(l),VAI(l),                         &
                       ! State variables                               &
                       albs(l),Tsrf(l),Dsnw(:,l),Nsnow(l),Qcan(:,l),   &
                       Rgrn(:,l),Sice(:,l),Sliq(:,l),Sveg(:,l),        &
                       Tcan(:,l),Tsnow(:,l),Tsoil(:,l),Tveg(:,l),      &
                       Vsmc(:,l),                                      &
                       ! Diagnostics                                   &
                       fsnow(l),H(l),LE(l),LWout(l),LWsub(l),Melt(l),  &
                       Roff(l),snd(l),snw(l),subl(l),svg(l),SWout(l),  &
                       SWsub(l),Tsub(l),Usub(l),Wflx(:,l)              )   
  end do
  SWup(:) = SWup(:) + SWout(:)*dt
  if (mod(n,Nout)==0) then
    status = nf90_put_var(ncid,varids(3),hour,start=(/rec/))
    var2d(:,:) = 1.e+20
    do l = 1, Nland
      var2d(x(l),y(l)) = fsnow(l)
    end do
    status = nf90_put_var(ncid,varids(4),var2d,start=(/1,1,rec/),      &
                          count=(/Nx,Ny,1/))
    do l = 1, Nland
      var2d(x(l),y(l)) = snd(l)
    end do
    status = nf90_put_var(ncid,varids(5),var2d,start=(/1,1,rec/),      &
                          count=(/Nx,Ny,1/))
    do l = 1, Nland
      var2d(x(l),y(l)) = snw(l)
    end do
    status = nf90_put_var(ncid,varids(6),var2d,start=(/1,1,rec/),      &
                          count=(/Nx,Ny,1/))
    do l = 1, Nland
      var2d(x(l),y(l)) = Nsnow(l)
    end do
    status = nf90_put_var(ncid,varids(7),var2d,start=(/1,1,rec/),          &
                          count=(/Nx,Ny,1/))
    do l = 1, Nland
      var2d(x(l),y(l)) = Tsoil(1,l)
    end do
    status = nf90_put_var(ncid,varids(8),var2d,start=(/1,1,rec/),          &
                          count=(/Nx,Ny,1/))
    do l = 1, Nland
      var2d(x(l),y(l)) = SWup(l)
    end do
    status = nf90_put_var(ncid,varids(9),var2d,start=(/1,1,rec/),      &
                          count=(/Nx,Ny,1/))
    SWup(:) = 0
#if SMRTIO == 1
    var3d(:,:,:) = 1.e+20
    do l = 1, Nland
      var3d(x(l),y(l),:) = Dsnw(:Nsnow(l),l)
    end do
    status = nf90_put_var(ncid,varids(10),var3d,start=(/1,1,1,rec/),       &
                          count=(/Nx,Ny,Nsmax,1/))
    do l = 1, Nland  ! liquid fraction
      var3d(x(l),y(l),:Nsnow(l)) = Sliq(:Nsnow(l),l) /                     &
                                   (Sice(:Nsnow(l),l) + Sliq(:Nsnow(l),l))
    end do
    status = nf90_put_var(ncid,varids(11),var3d,start=(/1,1,1,rec/),       &
                          count=(/Nx,Ny,Nsmax,1/))
    do l = 1, Nland  ! density
      var3d(x(l),y(l),:Nsnow(l)) = (Sice(:Nsnow(l),l) + Sliq(:Nsnow(l),l)) &
                                   / Dsnw(:Nsnow(l),l)
    end do
    status = nf90_put_var(ncid,varids(12),var3d,start=(/1,1,1,rec/),       &
                          count=(/Nx,Ny,Nsmax,1/))
    do l = 1, Nland
      var3d(x(l),y(l),:) = Rgrn(:Nsnow(l),l)
    end do
    status = nf90_put_var(ncid,varids(13),var3d,start=(/1,1,1,rec/),       &
                          count=(/Nx,Ny,Nsmax,1/))
    do l = 1, Nland         
      var3d(x(l),y(l),:) = Tsnow(:Nsnow(l),l)
    end do
    status = nf90_put_var(ncid,varids(14),var3d,start=(/1,1,1,rec/),       &
                          count=(/Nx,Ny,Nsmax,1/))
#endif
    rec = rec + 1
  end if
  hour = hour + dt/3600
end do
status=nf90_close(ncid) 
write(6,*) 'minimum and maximum SWE',minval(snw),maxval(snw)

! Write out state variables at end of run
open(8,file='dump.txt')
write(8,*) albs
write(8,*) Dsnw
write(8,*) Nsnow
write(8,*) Qcan
write(8,*) Rgrn
write(8,*) Sice
write(8,*) Sliq
write(8,*) Sveg
write(8,*) Tcan
write(8,*) Tsnow
write(8,*) Tsoil
write(8,*) Tsrf
write(8,*) Tveg
write(8,*) Vsmc
close(8)

end program FSM2_GLOBAL
