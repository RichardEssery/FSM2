!-----------------------------------------------------------------------
! Set parameters, initialize prognostic variables and write metadata
!-----------------------------------------------------------------------
subroutine SETUP

#include "OPTS.h"
use CONSTANTS
use DIAGNOSTICS
use DRIVING
use GRID
use IOUNITS
use PARAMETERS
use PARAMMAPS
use SOILPARAMS 
use STATE_VARIABLES 

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
namelist /params/ asmx,asmn,avg0,avgs,bstb,bthr,cden,cvai,cveg,eta0,etaa,etab,gsat,hfsn,  &
                  kext,kfix,rchd,rchz,rho0,rhoc,rhof,rcld,rmlt,Salb,snda,sndb,sndc,Talb,  &
                  tcnc,tcnm,tcld,tmlt,trho,Wirr,z0sn,z0zh

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
dt = 3600
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
cden = 0.004
cvai = 4.4
cveg = 20
kext = 0.5
rchd = 0.67
rchz = 0.1
tcnc = 240
tcnm = 2.4

! Defaults for snow parameters
asmx = 0.8
asmn = 0.5
bstb = 5
bthr = 2
eta0 = 3.7e7
etaa = 0.081
etab = 0.018
hfsn = 0.1
kfix = 0.24
rho0 = 300
rhoc = 150
rhof = 100
rcld = 300
rmlt = 500
Salb = 10
snda = 2.8e-6
sndb = 0.042
sndc = 0.046
Talb = -2
tcld = 1000
tmlt = 100
trho = 200
Wirr = 0.03
z0sn = 0.01

! Defaults for ground surface parameters
bstb = 5
gsat = 0.01
z0zh = 10

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
if (fveg(1,1) < 0) fveg(:,:) = 1 - exp(-VAI(:,:))
if (scap(1,1) < 0) scap(:,:) = cvai*VAI(:,:)

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
  hcap_soil(i,j) = (2.128*fcly(i,j) + 2.385*fsnd(i,j))*1e6 / (fcly(i,j) + fsnd(i,j))
  sathh(i,j) = 10**(0.17 - 0.63*fcly(i,j) - 1.58*fsnd(i,j))
  Vsat(i,j) = 0.505 - 0.037*fcly(i,j) - 0.142*fsnd(i,j)
  Vcrit(i,j) = Vsat(i,j)*(sathh(i,j)/3.364)**(1/b(i,j))
  hcon_min = (hcon_clay**fcly(i,j)) * (hcon_sand**(1 - fcly(i,j)))
  hcon_soil(i,j) = (hcon_air**Vsat(i,j)) * (hcon_min**(1 - Vsat(i,j)))
end do
end do

! Convert time scales from hours to seconds
tcnc = 3600*tcnc
tcnm = 3600*tcnm
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
allocate(Tsrf(Nx,Ny))
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
Tsrf(:,:) = Tsoil(1,:,:)

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
  read(ustr,*) Tsrf(:,:)
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
if (runid == 'none') runid = ''
open(uave, file = trim(runid) // trim(ave_file))
open(udmp, file = trim(runid) // trim(dmp_file))
open(usmp, file = trim(runid) // trim(smp_file))

! Write options and namelists to metadata file
#if ALBEDO == 0
#define ALBEDO_OPT 'diagnostic'
#elif ALBEDO == 1
#define ALBEDO_OPT 'prognostic'
#endif
#if CANMOD == 0
#define CANMOD_OPT 'zero layer'
#elif CANMOD == 1
#define CANMOD_OPT 'one layer'
#endif
#if CONDCT == 0
#define CONDCT_OPT 'constant'
#elif CONDCT == 1
#define CONDCT_OPT 'Yen (1981)'
#endif
#if DENSTY == 0
#define DENSTY_OPT 'constant'
#elif DENSTY == 1
#define DENSTY_OPT 'Verseghy (1991)'
#elif DENSTY == 2
#define DENSTY_OPT 'Anderson (1976)'
#endif
#if EXCHNG == 0
#define EXCHNG_OPT 'constant'
#elif EXCHNG == 1
#define EXCHNG_OPT 'Louis (1979)'
#endif
#if HYDROL == 0
#define HYDROL_OPT 'free draining'
#elif HYDROL == 1
#define HYDROL_OPT 'bucket'
#endif
open(umta, file =  trim(runid) // 'runinfo')
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
