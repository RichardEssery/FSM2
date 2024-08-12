!-----------------------------------------------------------------------
! Call FSM2 physics subroutines for one timestep at one point
!-----------------------------------------------------------------------
subroutine FSM2_TIMESTEP(dt,elev,zT,zU,                                &
                         LW,Ps,Qa,Rf,Sdif,Sdir,Sf,Ta,trans,Ua,         &
                         alb0,vegh,VAI,                                &
                         albs,Tsrf,Dsnw,Nsnow,Qcan,Rgrn,Sice,Sliq,     &
                         Sveg,Tcan,Tsnow,Tsoil,Tveg,Vsmc,              &
                         fsnow,H,LE,LWout,LWsub,Melt,Roff,snd,snw,     &
                         subl,svg,SWout,SWsub,Tsub,Usub,Wflx)

use LAYERS, only: &
  Ncnpy,             &! Number of canopy layers
  Nsmax,             &! Maximum number of snow layers
  Nsoil               ! Number of soil layers

implicit none

! Site characteristics
real, intent(in) :: &!
  dt,                &! Timestep (s)
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind speed measurement height (m)

! Meteorological variables
real, intent(in) :: &
  elev,              &! Solar elevation (radians)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Rf,                &! Rainfall rate (kg/m2/s)
  Sdif,              &! Diffuse shortwave radiation (W/m^2)
  Sdir,              &! Direct-beam shortwave radiation (W/m^2)
  Ta,                &! Air temperature (K)
  trans,             &! Wind-blown snow transport rate (kg/m^2/s)
  Ua                  ! Wind speed (m/s)
real, intent(inout) :: &
  Sf                  ! Snowfall rate (kg/m2/s)

! Vegetation characteristics
real, intent(in) :: &
  alb0,              &! Snow-free ground albedo
  vegh,              &! Canopy height (m)
  VAI                 ! Vegetation area index

! State variables
integer, intent(inout) :: &
  Nsnow               ! Number of snow layers
real, intent(inout) :: &
  albs,              &! Snow albedo
  Tsrf,              &! Snow/ground surface temperature (K)
  Dsnw(Nsmax),       &! Snow layer thicknesses (m)
  Qcan(Ncnpy),       &! Canopy air space humidities
  Rgrn(Nsmax),       &! Snow layer grain radii (m)
  Sice(Nsmax),       &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmax),       &! Liquid content of snow layers (kg/m^2)
  Sveg(Ncnpy),       &! Snow mass on vegetation layers (kg/m^2)
  Tcan(Ncnpy),       &! Canopy air space temperatures (K)
  Tsnow(Nsmax),      &! Snow layer temperatures (K)
  Tsoil(Nsoil),      &! Soil layer temperatures (K)
  Tveg(Ncnpy),       &! Vegetation layer temperatures (K)
  Vsmc(Nsoil)         ! Volumetric moisture content of soil layers

! Diagnostics
real, intent(out) :: &
  fsnow,             &! Ground snowcover fraction
  H,                 &! Sensible heat flux to the atmosphere (W/m^2)
  LE,                &! Latent heat flux to the atmosphere (W/m^2)
  LWout,             &! Outgoing LW radiation (W/m^2)
  LWsub,             &! Subcanopy downward LW radiation (W/m^2)
  Melt,              &! Surface melt rate (kg/m^2/s)
  Roff,              &! Runoff from snow (kg/m^2/s)
  snd,               &! Snow depth (m)
  snw,               &! Total snow mass on ground (kg/m^2) 
  subl,              &! Sublimation rate (kg/m^2/s)
  svg,               &! Total snow mass on vegetation (kg/m^2)
  SWout,             &! Outgoing SW radiation (W/m^2)
  SWsub,             &! Subcanopy downward SW radiation (W/m^2)
  Tsub,              &! Subcanopy air temperature (K)
  Usub,              &! Subcanopy wind speed (m/s)
  Wflx(Nsmax)         ! Water flux into snow layer (kg/m^2/s)

! Vegetation properties
real :: &
  cveg(Ncnpy),       &! Vegetation heat capacities (J/K/m^2)
  fcans(Ncnpy),      &! Canopy layer snowcover fractions
  lveg(Ncnpy),       &! Canopy layer vegetation area indices
  Scap(Ncnpy),       &! Canopy layer snow capacities (kg/m^2)
  tdif(Ncnpy),       &! Canopy layer diffuse transmittances
  Tveg0(Ncnpy)        ! Vegetation temperatures at start of timestep (K)

! Snow properties
real :: &
  ksnow(Nsmax)        ! Thermal conductivities of snow layers (W/m/K)

! Surface properties
real :: &
  Ds1,               &! Surface layer thickness (m)
  gs1,               &! Surface moisture conductance (m/s)
  ks1,               &! Surface layer thermal conductivity (W/m/K)
  Ts1                 ! Surface layer temperature (K)

! Soil properties
real :: &
  csoil(Nsoil),      &! Areal heat capacity of soil layers (J/K/m^2)
  ksoil(Nsoil)        ! Thermal conductivity of soil layers (W/m/K)

! Fluxes
real :: &
  drip,              &! Melt water drip from vegetation (kg/m^2)
  Esrf,              &! Moisture flux from the surface (kg/m^2/s)
  Gsrf,              &! Heat flux into snow/ground surface (W/m^2)
  Gsoil,             &! Heat flux into soil (W/m^2)
  SWsrf,             &! SW absorbed by snow/ground surface (W/m^2)
  unload,            &! Snow mass unloaded from vegetation (kg/m^2)
  Eveg(Ncnpy),       &! Moisture flux from vegetation layers (kg/m^2/s)
  SWveg(Ncnpy)        ! SW absorbed by vegetation layers (W/m^2)

call CANOPY(Sveg,Tveg,VAI,cveg,fcans,lveg,Scap,Tveg0)

call SWRAD(alb0,Dsnw,dt,elev,fcans,lveg,Sdif,Sdir,Sf,Tsrf,             &
           albs,fsnow,SWout,SWsrf,SWsub,SWveg,tdif)

call THERMAL(Dsnw,Nsnow,Sice,Sliq,Tsnow,Tsoil,Vsmc,csoil,Ds1,          &
             gs1,ksnow,ksoil,ks1,Ts1)

call SRFEBAL(cveg,Ds1,dt,fcans,fsnow,gs1,ks1,lveg,LW,Ps,Qa,            &
             SWsrf,Sveg,SWveg,Ta,tdif,Ts1,Tveg0,Ua,VAI,vegh,zT,zU,     &
             Tsrf,Qcan,Sice,Tcan,Tveg,                                 &
             Esrf,Eveg,Gsrf,H,LE,LWout,LWsub,Melt,subl,Tsub,Usub)

call INTERCEPT(dt,cveg,Eveg,lveg,Scap,Ta,Ua,Sf,Sveg,Tveg,drip,svg,unload)

call SNOW(dt,drip,Esrf,Gsrf,ksnow,ksoil,Melt,Rf,Sf,Ta,trans,Tsrf,unload, &
          Nsnow,Dsnw,Rgrn,Sice,Sliq,Tsnow,Tsoil,Gsoil,Roff,snd,snw,Wflx)

call SOIL(csoil,dt,Gsoil,ksoil,Tsoil)

end subroutine FSM2_TIMESTEP
