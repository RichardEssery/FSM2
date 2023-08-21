
!-----------------------------------------------------------------------
! Read and flatten a (x,y,time) netCDF data file
!-----------------------------------------------------------------------
subroutine READNC_LAND(Nland,Ntime,x0,y0,Nx,Ny,x,y,varname,var)

use netcdf

implicit none

character(len=*), intent(in) :: &
  varname             ! Variable name

integer, intent(in) :: &
  Nland,             &! Number of land points
  Ntime,             &! Number of timesteps
  Nx,                &! Number of x points in simulation domain
  Ny,                &! Number of y points in simulation domain
  x0,                &! First x point in simulation domain
  y0,                &! First y point in simulation domain
  x(Nland),          &! x indices of land points
  y(Nland)            ! y indices of land points

real, intent(out) :: &
  var(Nland,Ntime)    ! Variable on land points

integer :: &
  l,                 &! Land point counter
  ncid,              &! NetCDF dataset ID
  varid,             &! NetCDF variable ID
  status              ! NetCDF error status

real :: &
  var2d(Nx,Ny,Ntime) ! Variable on grid
  
print*,'Reading ',varname
status = nf90_open(varname//'.nc',nf90_nowrite,ncid)
status = nf90_inq_varid(ncid,varname,varid)
status = nf90_get_var(ncid,varid,var2d,start=(/x0,y0,1/),  &
                                 count=(/Nx,Ny,Ntime/))
status = nf90_close(ncid)

do l = 1, Nland
  var(l,:) = var2d(x(l),y(l),:)
end do

end subroutine READNC_LAND

