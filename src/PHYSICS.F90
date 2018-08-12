!-----------------------------------------------------------------------
! Call physics subroutines
!-----------------------------------------------------------------------
subroutine PHYSICS

use GRID, only: &
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  Nitr                ! Number of iterations in energy balance calulation

implicit none

! Eddy diffusivities
real :: &
  KH(Nx,Ny),         &! Eddy diffusivity for heat to the atmosphere (m/s)
  KHa(Nx,Ny),        &! Eddy diffusivity from the canopy air space (m/s)
  KHg(Nx,Ny),        &! Eddy diffusivity for heat from the ground (m/s)
  KHv(Nx,Ny),        &! Eddy diffusivity for heat from vegetation (m/s)
  KWg(Nx,Ny),        &! Eddy diffusivity for water from the ground (m/s)
  KWv(Nx,Ny)          ! Eddy diffusivity for water from vegetation (m/s)

! Surface properties
real :: &
  alb(Nx,Ny),        &! Albedo
  Ds1(Nx,Ny),        &! Surface layer thickness (m)
  fsnow(Nx,Ny),      &! Ground snowcover fraction
  gs1(Nx,Ny),        &! Surface moisture conductance (m/s)
  ks1(Nx,Ny),        &! Surface thermal conductivity (W/m/K)
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
  Esrf(Nx,Ny),       &! Moisture flux from the surface (kg/m^2/s)
  Eveg(Nx,Ny),       &! Moisture flux from vegetation (kg/m^2/s)
  G(Nx,Ny),          &! Heat flux into surface (W/m^2)
  Gsoil(Nx,Ny),      &! Heat flux into soil (W/m^2)
  H(Nx,Ny),          &! Sensible heat flux to the atmosphere (W/m^2)
  Hsrf(Nx,Ny),       &! Sensible heat flux from the surface (W/m^2)
  LE(Nx,Ny),         &! Latent heat flux to the atmosphere (W/m^2)
  LEsrf(Nx,Ny),      &! Latent heat flux from the surface (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny),       &! Net radiation (W/m^2)
  Roff(Nx,Ny),       &! Runoff from snow (kg/m^2)
  Rsrf(Nx,Ny),       &! Net radiation absorbed by the surface (W/m^2)
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny),      &! Net SW radiation absorbed by vegetation (W/m^2)
  Tveg0(Nx,Ny),      &! Vegetation temperature at start of timestep (K)
  unload(Nx,Ny)       ! Snow mass unloaded from canopy (kg/m^2)

integer :: & 
  n                   ! Iteration counter

call RADIATION(alb,fsnow,SWsrf,SWveg)

call THERMAL(csoil,Ds1,gs1,ks1,ksnow,ksoil,Ts1,Tveg0)

do n = 1, Nitr

  call SFEXCH(fsnow,gs1,KH,KHa,KHg,KHv,KWg,KWv)

  call EBALFOR(Ds1,KHa,KHg,KHv,KWg,KWv,ks1,SWsrf,SWveg,Ts1,Tveg0, &
               Esrf,Eveg,G,H,Hsrf,LE,LEsrf,Melt,Rnet,Rsrf)

  call EBALSRF(Ds1,KH,KHa,KHv,KWg,KWv,ks1,SWsrf,SWveg,Ts1, &
               Esrf,Eveg,G,H,Hsrf,LE,LEsrf,Melt,Rnet,Rsrf)

end do

call CANOPY(Eveg,unload)

call SNOW(Esrf,G,ksnow,ksoil,Melt,unload,Gsoil,Roff)

call SOIL(csoil,Gsoil,ksoil)

call CUMULATE(alb,G,H,Hsrf,LE,LEsrf,Melt,Rnet,Roff,Rsrf)

end subroutine PHYSICS
