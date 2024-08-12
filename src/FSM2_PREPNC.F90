!-----------------------------------------------------------------------
! Prepare NetCDF output file
!-----------------------------------------------------------------------
subroutine FSM2_PREPNC(runid,year,month,day,hour,ncid,rec,varid)

use netcdf

use LAYERS, only: &
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsmax,             &! Maximum number of snow layers
  Nsoil               ! Number of soil layers

implicit none

character(len=*), intent(in) :: &
  runid               ! Run identifier

integer, intent(in) :: &
  year,              &! Year
  month,             &! Month of year
  day                 ! Day of month

real, intent(in) :: &
  hour                ! Hour of day

integer,intent(out) :: &
  ncid,              &! Dataset ID
  rec,               &! Record number
  varid(17)           ! Variable IDs

character(len=31) :: &
  units               ! NetCDF time units

integer :: &
  nvar,              &! Variable number
  snow_dimid,        &! Snow layer dimension ID
  soil_dimid,        &! Soil layer dimension ID
  time_dimid,        &! Time dimension ID
  status              ! Error status

! Start time written to time units
write(units,'(A12,i4,A1,i2.2,A1,i2.2,A1,i2.2,A6)')  &
            'hours since ',year,'-',month,'-',day,' ',int(hour),':00:00'

! Create NetCDF file
status = nf90_create(trim(runid)//'FSM2out.nc',nf90_clobber,ncid)

! Dimensions
status = nf90_def_dim(ncid,'snow_layer',Nsmax,snow_dimid)
status = nf90_def_dim(ncid,'soil_layer',Nsoil,soil_dimid)
status = nf90_def_dim(ncid,'time',NF90_UNLIMITED,time_dimid)

! Define time variable
nvar = 1
status = nf90_def_var(ncid,'time',NF90_REAL,time_dimid,varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','time')
status = nf90_put_att(ncid,varid(nvar),'units',units)
status = nf90_put_att(ncid,varid(nvar),'calendar','standard')

! Define surface variables

nvar = nvar + 1
status = nf90_def_var(ncid,'hfls',NF90_REAL,time_dimid,varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','surface upward latent heat flux')
status = nf90_put_att(ncid,varid(nvar),'units','W m-2')

nvar = nvar + 1
status = nf90_def_var(ncid,'hfss',NF90_REAL,time_dimid,varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','surface upward sensible heat flux')
status = nf90_put_att(ncid,varid(nvar),'units','W m-2')

nvar = nvar + 1
status = nf90_def_var(ncid,'rlus',NF90_REAL,time_dimid,varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','surface upwelling longwave radiation')
status = nf90_put_att(ncid,varid(nvar),'units','W m-2')

nvar = nvar + 1
status = nf90_def_var(ncid,'rsus',NF90_REAL,time_dimid,varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','surface upwelling shortwave radiation')
status = nf90_put_att(ncid,varid(nvar),'units','W m-2')

nvar = nvar + 1
status = nf90_def_var(ncid,'snm',NF90_REAL,time_dimid,varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','surface snow melt')
status = nf90_put_att(ncid,varid(nvar),'units','kg m-2 s-1')

nvar = nvar + 1
status = nf90_def_var(ncid,'ts',NF90_REAL,time_dimid,varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','surface temperature')
status = nf90_put_att(ncid,varid(nvar),'units','K')

! Define bulk snow variables

nvar = nvar + 1
status = nf90_def_var(ncid,'snd',NF90_REAL,time_dimid,varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','snow depth')
status = nf90_put_att(ncid,varid(nvar),'units','m')

nvar = nvar + 1
status = nf90_def_var(ncid,'snmsl',NF90_REAL,time_dimid,varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','water flowing out of snowpack')
status = nf90_put_att(ncid,varid(nvar),'units','kg m-2 s-1')

nvar = nvar + 1
status = nf90_def_var(ncid,'snw',NF90_REAL,time_dimid,varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','surface snow mass')
status = nf90_put_att(ncid,varid(nvar),'units','kg m-2')

! Define snow layer variables

nvar = nvar + 1
status = nf90_def_var(ncid,'Dsnw',NF90_REAL,(/snow_dimid,time_dimid/),varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','thickness of snow layers')
status = nf90_put_att(ncid,varid(nvar),'units','m')

nvar = nvar + 1
status = nf90_def_var(ncid,'lqsn',NF90_REAL,(/snow_dimid,time_dimid/),varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','mass fraction of liquid water in snow layers')
status = nf90_put_att(ncid,varid(nvar),'units','-')

nvar = nvar + 1
status = nf90_def_var(ncid,'snowrho',NF90_REAL,(/snow_dimid,time_dimid/),varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','density of snow layers')
status = nf90_put_att(ncid,varid(nvar),'units','kg m-3')

nvar = nvar + 1
status = nf90_def_var(ncid,'rgrn',NF90_REAL,(/snow_dimid,time_dimid/),varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'long_name','grain radius in snow layers')
status = nf90_put_att(ncid,varid(nvar),'units','m')

nvar = nvar + 1
status = nf90_def_var(ncid,'tsnl',NF90_REAL,(/snow_dimid,time_dimid/),varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'units','K')
status = nf90_put_att(ncid,varid(nvar),'long_name','temperature of snow layers')

nvar = nvar + 1
status = nf90_def_var(ncid,'wflx',NF90_REAL,(/snow_dimid,time_dimid/),varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'units','kg m-2 s-1')
status = nf90_put_att(ncid,varid(nvar),'long_name','water flux into snow layers')

! Define soil layer variables

nvar = nvar + 1
status = nf90_def_var(ncid,'tsl',NF90_REAL,(/soil_dimid,time_dimid/),varid(nvar))
status = nf90_put_att(ncid,varid(nvar),'units','K')
status = nf90_put_att(ncid,varid(nvar),'long_name','temperature of soil layers')

status = nf90_def_var(ncid,'Dzsoil',NF90_REAL,soil_dimid,nvar)
status = nf90_put_att(ncid,nvar,'long_name','soil layer thickness')
status = nf90_put_att(ncid,nvar,'units','m')

status = nf90_enddef(ncid)

status = nf90_put_var(ncid,nvar,Dzsoil)

rec = 1

end subroutine FSM2_PREPNC
