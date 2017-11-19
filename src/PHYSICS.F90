!-----------------------------------------------------------------------
! Call physics subroutines
!-----------------------------------------------------------------------
subroutine PHYSICS

use GRID, only: &
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions

implicit none

! Eddy diffusivities
real :: &
  KH(Nx,Ny),         &! Eddy diffusivity for heat fluxes (m/s)
  KHsurf(Nx,Ny),     &! Surface eddy diffusivity (m/s)
  KHveg(Nx,Ny)        ! Vegetation eddy diffusivity (m/s)

! Surface properties
real :: &
  alb(Nx,Ny),        &! Albedo
  Dz1(Nx,Ny),        &! Surface layer thickness (m)
  fsnow(Nx,Ny),      &! Snowcover fraction
  gevap(Nx,Ny),      &! Surface moisture conductance (m/s)
  ksurf(Nx,Ny),      &! Surface thermal conductivity (W/m/K)
  Ts1(Nx,Ny)          ! Surface layer temperature (K)

! Snow properties
real :: &
  ksnow(Nsmax,Nx,Ny)  ! Thermal conductivity of snow (W/m/K)

! Soil properties
real :: &
  csoil(Nsoil,Nx,Ny),&! Areal heat capacity of soil (J/K/m^2)
  ksoil(Nsoil,Nx,Ny)  ! Thermal conductivity of soil (W/m/K)

! Fluxes
real :: &
  Esurf(Nx,Ny),      &! Moisture flux from the surface (kg/m^2/s)
  Eveg(Nx,Ny),       &! Moisture flux from vegetation (kg/m^2/s)
  Gsurf(Nx,Ny),      &! Heat flux into surface (W/m^2)
  Gsoil(Nx,Ny),      &! Heat flux into soil (W/m^2)
  Hatmo(Nx,Ny),      &! Sensible heat flux to the atmosphere (W/m^2)
  Latmo(Nx,Ny),      &! Latent heat flux to the atmosphere (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny),       &! Net radiation (W/m^2)
  Roff(Nx,Ny),       &! Runoff from snow (kg/m^2)
  SWsurf(Nx,Ny),     &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny)        ! Net SW radiation absorbed by vegetation (W/m^2)

call SWRAD(alb,fsnow,SWsurf,SWveg)

call THERMAL(csoil,Dz1,gevap,ksnow,ksoil,ksurf,Ts1)

call SFEXCH(fsnow,KH,KHsurf,KHveg)

call EBALFOR(Dz1,gevap,KH,KHsurf,KHveg,ksurf,SWsurf,SWveg,Ts1, &
             Esurf,Eveg,Gsurf,Hatmo,Latmo,Melt,Rnet)

call EBALOPN(Dz1,gevap,KH,ksurf,SWsurf,Ts1,Esurf,Gsurf,Hatmo,Latmo,Melt,Rnet)

call CANOPY(Eveg)

call SNOW(Esurf,Gsurf,ksnow,ksoil,Melt,Gsoil,Roff)

call SOIL(csoil,Gsoil,ksoil)

call CUMULATE(alb,Gsurf,Hatmo,Latmo,Melt,Rnet,Roff)

end subroutine PHYSICS
