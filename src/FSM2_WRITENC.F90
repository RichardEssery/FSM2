!-----------------------------------------------------------------------
! Write variabls for one timestep to the NetCDF output file
!-----------------------------------------------------------------------
subroutine FSM2_WRITENC(Dsnw,dt,H,LE,LWout,Melt,ncid,Nsnow,Rgrn,Roff,  &
                        Sice,Sliq,snd,snw,SWout,Tsnow,Tsoil, Tsrf,     &
                        varid,Wflx,rec)

use netcdf

use LAYERS, only: &
  Nsmax,             &! Maximum number of snow layers
  Nsoil               ! Number of soil layers

implicit none

integer, intent(in) :: &
  ncid,              &! Dataset ID
  Nsnow,             &! Number of snow layers
  varid(17)           ! Variable IDs

real, intent(in) :: &
  dt,                &! Timestep (s)
  H,                 &! Sensible heat flux to the atmosphere (W/m^2)
  LE,                &! Latent heat flux to the atmosphere (W/m^2)
  LWout,             &! Outgoing LW radiation (W/m^2)
  Melt,              &! Surface melt rate (kg/m^2/s)
  Roff,              &! Runoff from snow (kg/m^2/s)
  snd,               &! Snow depth (m)
  snw,               &! Total snow mass on ground (kg/m^2) 
  SWout,             &! Outgoing SW radiation (W/m^2)
  Tsrf,              &! Snow/ground surface temperature (K)
  Dsnw(Nsmax),       &! Snow layer thicknesses (m)
  Rgrn(Nsmax),       &! Snow layer grain radii (m)
  Sice(Nsmax),       &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmax),       &! Liquid content of snow layers (kg/m^2)
  Tsnow(Nsmax),      &! Snow layer temperatures (K)
  Tsoil(Nsoil),      &! Soil layer temperatures (K)
  Wflx(Nsmax)         ! Water flux into snow layer (kg/m^2/s)

integer,intent(inout) :: &
  rec                 ! Record number

integer :: &
  count(2),          &! Number of points to write
  start(2),          &! Starting point
  status              ! Error status

real :: &
  time,              &! Time (hours)
  var(Nsmax)          ! Variable on snow layers

! Time
time = (rec - 1)*dt/3600.
status = nf90_put_var(ncid,varid(1),time,start=(/rec/))

! Surface variables
status = nf90_put_var(ncid,varid(2),LE,start=(/rec/))
status = nf90_put_var(ncid,varid(3),H,start=(/rec/))
status = nf90_put_var(ncid,varid(4),LWout,start=(/rec/))
status = nf90_put_var(ncid,varid(5),SWout,start=(/rec/))
status = nf90_put_var(ncid,varid(6),Melt,start=(/rec/))
status = nf90_put_var(ncid,varid(7),Tsrf,start=(/rec/))

! Bulk snow variables
status = nf90_put_var(ncid,varid(8),snd,start=(/rec/))
status = nf90_put_var(ncid,varid(9),Roff,start=(/rec/))
status = nf90_put_var(ncid,varid(10),snw,start=(/rec/))

! Snow layer variables
start = (/1,rec/)
count = (/Nsmax,1/)
var(:) = 0
var(:Nsnow) = Dsnw(:Nsnow)
status = nf90_put_var(ncid,varid(11),var,start=start,count=count)
var(:Nsnow) = Sliq(:Nsnow) / (Sice(:Nsnow) + Sliq(:Nsnow))
var(Nsnow:) = var(Nsnow)
status = nf90_put_var(ncid,varid(12),var,start=start,count=count)
var(:Nsnow) = (Sice(:Nsnow) + Sliq(:Nsnow)) / Dsnw(:Nsnow)
var(Nsnow:) = var(Nsnow)
status = nf90_put_var(ncid,varid(13),var,start=start,count=count)
var(:Nsnow) = Rgrn(:Nsnow)
var(Nsnow:) = var(Nsnow)
status = nf90_put_var(ncid,varid(14),var,start=start,count=count)
var(:Nsnow) = Tsnow(:Nsnow)
var(Nsnow:) = var(Nsnow)
status = nf90_put_var(ncid,varid(15),var,start=start,count=count)
var(:) = 0
var(:Nsnow) = Wflx(:Nsnow)
if (Nsnow>0) var(Nsnow:) = Roff
status = nf90_put_var(ncid,varid(16),var,start=start,count=count)

! Soil layer variables
start = (/1,rec/)
count = (/Nsoil,1/)
status = nf90_put_var(ncid,varid(17),Tsoil,start=start,count=count)

rec = rec + 1

end subroutine FSM2_WRITENC
