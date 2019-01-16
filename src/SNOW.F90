!-----------------------------------------------------------------------
! Snow thermodynamics and hydrology
!-----------------------------------------------------------------------
subroutine SNOW(Esrf,G,ksnow,ksoil,Melt,unload,Gsoil,Roff)

#include "OPTS.h"

use CONSTANTS, only: &
  grav,              &! Acceleration due to gravity (m/s^2)
  hcap_ice,          &! Specific heat capacity of ice (J/K/kg)
  hcap_wat,          &! Specific heat capacity of water (J/K/kg)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  rho_ice,           &! Density of ice (kg/m^3)
  rho_wat,           &! Density of water (kg/m^3)
  Tm                  ! Melting point (K)

use DRIVING, only: &
  dt,                &! Timestep (s)
  Rf,                &! Rainfall rate (kg/m^2/s)
  Sf,                &! Snowfall rate (kg/m^2/s)
  Ta,                &! Air temperature (K)
  Ua                  ! Wind speed (m/s)

use GRID, only: &
  Dzsnow,            &! Minimum snow layer thicknesses (m)
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  eta0,              &! Reference snow viscosity (Pa s)
  rgr0,              &! Fresh snow grain radius (m)
  rho0,              &! Fixed snow density (kg/m^3)
  rhob,              &! Temperature factor in fresh snow density (kg/m^3/K)
  rhoc,              &! Wind factor in fresh snow density (kg s^0.5/m^3.5)
  rhof,              &! Fresh snow density (kg/m^3)
  rcld,              &! Maximum density for cold snow (kg/m^3)
  rmlt,              &! Maximum density for melting snow (kg/m^3)
  snda,              &! Thermal metamorphism parameter (1/s)
  trho,              &! Snow compaction timescale (s)
  Wirr                ! Irreducible liquid water content of snow

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers
  rgrn,              &! Snow layer grain radius (m)
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Tsnow,             &! Snow layer temperatures (K)
  Tsoil,             &! Soil layer temperatures (K)
  Tsrf                ! Surface skin temperature (K)

implicit none

real, intent(in) :: &
  Esrf(Nx,Ny),       &! Moisture flux from the surface (kg/m^2/s)
  G(Nx,Ny),          &! Heat flux into surface (W/m^2)
  ksnow(Nsmax,Nx,Ny),&! Thermal conductivity of snow (W/m/K)
  ksoil(Nsoil,Nx,Ny),&! Thermal conductivity of soil (W/m/K)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  unload(Nx,Ny)       ! Snow mass unloaded from canopy (kg/m^2)

real, intent(out) :: &
  Gsoil(Nx,Ny),      &! Heat flux into soil (W/m^2)
  Roff(Nx,Ny)         ! Runoff from snow (kg/m^2)

integer :: &
  i,j,               &! Point counters
  k,                 &! Level counter
  knew,              &! New snow layer pointer
  kold,              &! Old snow layer pointer
  Nold                ! Previous number of snow layers

real :: &
  coldcont,          &! Layer cold content (J/m^2)
  dnew,              &! New snow layer thickness (m)
  dSice,             &! Change in layer ice content (kg/m^2)
  Esnow,             &! Snow sublimation rate (kg/m^2/s)
  ggr,               &! Grain area growth rate (m^2/s)
  mass,              &! Mass of overlying snow (kg/m^2)
  phi,               &! Porosity
  rhonew,            &! Density of new snow (kg/m^3)
  rhos,              &! Density of snow layer (kg/m^3)
  SliqMax,           &! Maximum liquid content for layer (kg/m^2)
  snowdepth,         &! Snow depth (m)
  wt                  ! Layer weighting

real :: &
  a(Nsmax),          &! Below-diagonal matrix elements
  b(Nsmax),          &! Diagonal matrix elements
  c(Nsmax),          &! Above-diagonal matrix elements
  csnow(Nsmax),      &! Areal heat capacity of snow (J/K/m^2)
  dTs(Nsmax),        &! Temperature increments (k)
  D(Nsmax),          &! Layer thickness before adjustment (m)
  E(Nsmax),          &! Energy contents before adjustment (J/m^2)
  Gs(Nsmax),         &! Thermal conductivity between layers (W/m^2/k)
  rhs(Nsmax),        &! Matrix equation rhs
  R(Nsmax),          &! Snow grain radii before adjustment (kg/m^2)
  S(Nsmax),          &! Ice contents before adjustment (kg/m^2)
  U(Nsmax),          &! Layer internal energy contents (J/m^2)
  W(Nsmax)            ! Liquid contents before adjustment (kg/m^2)

Gsoil(:,:) = G(:,:)
Roff(:,:) = Rf(:,:)*dt

! Points with existing snowpack
do j = 1, Ny
do i = 1, Nx
  if (Nsnow(i,j) > 0) then  

  ! Heat conduction
    do k = 1, Nsnow(i,j)
      csnow(k) = Sice(k,i,j)*hcap_ice + Sliq(k,i,j)*hcap_wat
    end do
    if (Nsnow(i,j) == 1) then
      Gs(1) = 2 / (Ds(1,i,j)/ksnow(1,i,j) + Dzsoil(1)/ksoil(1,i,j))
      dTs(1) = (G(i,j) + Gs(1)*(Tsoil(1,i,j) - Tsnow(1,i,j)))*dt /  &
               (csnow(1) + Gs(1)*dt)
    else
      do k = 1, Nsnow(i,j) - 1
        Gs(k) = 2 / (Ds(k,i,j)/ksnow(k,i,j) + Ds(k+1,i,j)/ksnow(k+1,i,j))
      end do
      a(1) = 0
      b(1) = csnow(1) + Gs(1)*dt
      c(1) = - Gs(1)*dt
      rhs(1) = (G(i,j) - Gs(1)*(Tsnow(1,i,j) - Tsnow(2,i,j)))*dt
      do k = 2, Nsnow(i,j) - 1
        a(k) = c(k-1)
        b(k) = csnow(k) + (Gs(k-1) + Gs(k))*dt
        c(k) = - Gs(k)*dt
        rhs(k) = Gs(k-1)*(Tsnow(k-1,i,j) - Tsnow(k,i,j))*dt  &
                 + Gs(k)*(Tsnow(k+1,i,j) - Tsnow(k,i,j))*dt 
      end do
      k = Nsnow(i,j)
      Gs(k) = 2 / (Ds(k,i,j)/ksnow(k,i,j) + Dzsoil(1)/ksoil(1,i,j))
      a(k) = c(k-1)
      b(k) = csnow(k) + (Gs(k-1) + Gs(k))*dt
      c(k) = 0
      rhs(k) = Gs(k-1)*(Tsnow(k-1,i,j) - Tsnow(k,i,j))*dt  &
               + Gs(k)*(Tsoil(1,i,j) - Tsnow(k,i,j))*dt
      call TRIDIAG(Nsnow(i,j),Nsmax,a,b,c,rhs,dTs)
    end if 
    do k = 1, Nsnow(i,j)
      Tsnow(k,i,j) = Tsnow(k,i,j) + dTs(k)
    end do
    k = Nsnow(i,j)
    Gsoil(i,j) = Gs(k)*(Tsnow(k,i,j) - Tsoil(1,i,j))

  ! Convert melting ice to liquid water
    dSice = Melt(i,j)*dt
    do k = 1, Nsnow(i,j)
      coldcont = csnow(k)*(Tm - Tsnow(k,i,j))
      if (coldcont < 0) then
        dSice = dSice - coldcont/Lf
        Tsnow(k,i,j) = Tm
      end if
      if (dSice > 0) then
        if (dSice > Sice(k,i,j)) then  ! Layer melts completely
          dSice = dSice - Sice(k,i,j)
          Ds(k,i,j) = 0
          Sliq(k,i,j) = Sliq(k,i,j) + Sice(k,i,j)
          Sice(k,i,j) = 0
        else                       ! Layer melts partially
          Ds(k,i,j) = (1 - dSice/Sice(k,i,j))*Ds(k,i,j)
          Sice(k,i,j) = Sice(k,i,j) - dSice
          Sliq(k,i,j) = Sliq(k,i,j) + dSice
          dSice = 0                ! Melt exhausted
        end if
      end if
    end do

  ! Remove snow by sublimation 
    dSice = max(Esrf(i,j), 0.)*dt
    if (dSice > 0) then
      do k = 1, Nsnow(i,j)
        if (dSice > Sice(k,i,j)) then  ! Layer sublimates completely
          dSice = dSice - Sice(k,i,j)
          Ds(k,i,j) = 0
          Sice(k,i,j) = 0
        else                       ! Layer sublimates partially
          Ds(k,i,j) = (1 - dSice/Sice(k,i,j))*Ds(k,i,j)
          Sice(k,i,j) = Sice(k,i,j) - dSice
          dSice = 0                ! Sublimation exhausted
        end if
      end do
    end if

  ! Snow hydraulics
#if HYDROL == 0
  ! Free-draining snow 
    do k = 1, Nsnow(i,j)
      Roff(i,j) = Roff(i,j) + Sliq(k,i,j)
      Sliq(k,i,j) = 0
    end do
#endif
#if HYDROL == 1
  ! Bucket storage 
    do k = 1, Nsnow(i,j)
      phi = 0
      if (Ds(k,i,j) > epsilon(Ds)) phi = 1 - Sice(k,i,j)/(rho_ice*Ds(k,i,j))
      SliqMax = rho_wat*Ds(k,i,j)*phi*Wirr
      Sliq(k,i,j) = Sliq(k,i,j) + Roff(i,j)
      Roff(i,j) = 0
      if (Sliq(k,i,j) > SliqMax) then       ! Liquid capacity exceeded
        Roff(i,j) = Sliq(k,i,j) - SliqMax   ! so drainage to next layer
        Sliq(k,i,j) = SliqMax
      end if
      coldcont = csnow(k)*(Tm - Tsnow(k,i,j))
      if (coldcont > 0) then       ! Liquid can freeze
        dSice = min(Sliq(k,i,j), coldcont/Lf)
        Sliq(k,i,j) = Sliq(k,i,j) - dSice
        Sice(k,i,j) = Sice(k,i,j) + dSice
        Tsnow(k,i,j) = Tsnow(k,i,j) + Lf*dSice/csnow(k)
      end if
    end do
#endif

  ! Snow compaction
#if DENSTY == 0
  ! Fixed snow density
    do k = 1, Nsnow(i,j)
      Ds(k,i,j) = (Sice(k,i,j) + Sliq(k,i,j)) / rho0
    end do
#endif
#if DENSTY == 1
  ! Snow compaction with age
    do k = 1, Nsnow(i,j)
      if (Ds(k,i,j) > epsilon(Ds)) then
        rhos = (Sice(k,i,j) + Sliq(k,i,j)) / Ds(k,i,j)
        if (Tsnow(k,i,j) >= Tm) then
            if (rhos < rmlt) rhos = rmlt + (rhos - rmlt)*exp(-dt/trho)
        else
            if (rhos < rcld) rhos = rcld + (rhos - rcld)*exp(-dt/trho)
        end if
        Ds(k,i,j) = (Sice(k,i,j) + Sliq(k,i,j)) / rhos
      end if
    end do
#endif
#if DENSTY == 2
  ! Snow compaction by overburden
    mass = 0
    do k = 1, Nsnow(i,j)
      mass = mass + 0.5*(Sice(k,i,j) + Sliq(k,i,j)) 
      if (Ds(k,i,j) > epsilon(Ds)) then
        rhos = (Sice(k,i,j) + Sliq(k,i,j)) / Ds(k,i,j)
        rhos = rhos + (rhos*grav*mass*dt/(eta0*exp(-(Tsnow(k,i,j) - Tm)/12.4 + rhos/55.6))   &
                    + dt*rhos*snda*exp((Tsnow(k,i,j) - Tm)/23.8 - max(rhos - 150, 0.)/21.7))
        Ds(k,i,j) = (Sice(k,i,j) + Sliq(k,i,j)) / rhos
      end if
      mass = mass + 0.5*(Sice(k,i,j) + Sliq(k,i,j))
    end do
#endif

  ! Snow grain growth
!#if SGRAIN == 0
    do k = 1, Nsnow(i,j)
      ggr = 2e-13
      if (Tsnow(k,i,j) < Tm) then
        if (rgrn(k,i,j) < 1.50e-4) then
          ggr = 2e-14
        else
          ggr = 7.3e-8*exp(-4600/Tsnow(k,i,j))
        end if
      end if
      rgrn(k,i,j) = rgrn(k,i,j) + dt*ggr/rgrn(k,i,j)
    end do
!#endif

  end if  ! Existing snowpack
end do
end do

! Add snow, recalculate layers
do j = 1, Ny
do i = 1, Nx

! Add snowfall and frost to layer 1 with fresh snow density and grain size
  Esnow = 0
  if (Esrf(i,j) < 0 .and. Tsrf(i,j) < Tm) Esnow = Esrf(i,j)
  dSice = (Sf(i,j) - Esnow)*dt
#if DENSTY == 0
  rhonew = rho0
#else
  rhonew = max(rhof + rhob*(Ta(i,j) - Tm) + rhoc*Ua(i,j)**0.5, 50.)
#endif
  Ds(1,i,j) = Ds(1,i,j) + dSice / rhonew
  if (Sice(1,i,j) + dSice > epsilon(Sice))  &
    rgrn(1,i,j) = (Sice(1,i,j)*rgrn(1,i,j) + dSice*rgr0) / (Sice(1,i,j) + dSice)
  Sice(1,i,j) = Sice(1,i,j) + dSice

! Add canopy unloading to layer 1 with bulk snow density and grain size
  rhos = rhof
  mass = sum(Sice(:,i,j)) + sum(Sliq(:,i,j))
  snowdepth = sum(Ds(:,i,j))
  if (snowdepth > snowdepth) rhos = mass / snowdepth
  Ds(1,i,j) = Ds(1,i,j) + unload(i,j) / rhos
  Sice(1,i,j) = Sice(1,i,j) + unload(i,j)

! New snowpack
  if (Nsnow(i,j) == 0 .and. Sice(1,i,j) > 0) then
    Nsnow(i,j) = 1
    Tsnow(1,i,j) = min(Ta(i,j), Tm)
  end if
  snowdepth = sum(Ds(:,i,j))

! Store state of old layers
  D(:) = Ds(:,i,j)
  R(:) = rgrn(:,i,j)
  S(:) = Sice(:,i,j)
  W(:) = Sliq(:,i,j)
  do k = 1, Nsnow(i,j)
    csnow(k) = Sice(k,i,j)*hcap_ice + Sliq(k,i,j)*hcap_wat
    E(k) = csnow(k)*(Tsnow(k,i,j) - Tm)
  end do
  Nold = Nsnow(i,j)

! Initialise new layers
  Ds(:,i,j) = 0
  rgrn(:,i,j) = 0
  Sice(:,i,j) = 0
  Sliq(:,i,j) = 0
  Tsnow(:,i,j) = Tm
  U(:) = 0
  Nsnow(i,j) = 0

  if (snowdepth > 0) then  ! Existing or new snowpack

  ! Re-assign and count snow layers
    dnew = snowdepth
    Ds(1,i,j) = dnew
    k = 1
    if (Ds(1,i,j) > Dzsnow(1)) then 
      do k = 1, Nsmax
        Ds(k,i,j) = Dzsnow(k)
        dnew = dnew - Dzsnow(k)
        if (dnew <= Dzsnow(k) .or. k == Nsmax) then
          Ds(k,i,j) = Ds(k,i,j) + dnew
          exit
        end if
      end do
    end if
    Nsnow(i,j) = k

  ! Fill new layers from the top downwards
    knew = 1
    dnew = Ds(1,i,j)
    do kold = 1, Nold
      do
        if (D(kold) < dnew) then
        ! All snow from old layer partially fills new layer
          rgrn(knew,i,j) = rgrn(knew,i,j) + S(kold)*R(kold)
          Sice(knew,i,j) = Sice(knew,i,j) + S(kold)
          Sliq(knew,i,j) = Sliq(knew,i,j) + W(kold)
          U(knew) = U(knew) + E(kold)
          dnew = dnew - D(kold)
          exit
        else
        ! Some snow from old layer fills new layer
          wt = dnew / D(kold)
          rgrn(knew,i,j) = rgrn(knew,i,j) + wt*S(kold)*R(kold)
          Sice(knew,i,j) = Sice(knew,i,j) + wt*S(kold) 
          Sliq(knew,i,j) = Sliq(knew,i,j) + wt*W(kold)
          U(knew) = U(knew) + wt*E(kold)
          D(kold) = (1 - wt)*D(kold)
          E(kold) = (1 - wt)*E(kold)
          S(kold) = (1 - wt)*S(kold)
          W(kold) = (1 - wt)*W(kold)
          knew = knew + 1
          if (knew > Nsnow(i,j)) exit
          dnew = Ds(knew,i,j)
        end if
      end do
    end do

  ! Diagnose snow layer temperatures
    do k = 1, Nsnow(i,j)
     csnow(k) = Sice(k,i,j)*hcap_ice + Sliq(k,i,j)*hcap_wat
     Tsnow(k,i,j) = Tm + U(k) / csnow(k)
     rgrn(k,i,j) = rgrn(k,i,j) / Sice(k,i,j)
    end do

  end if ! Existing or new snowpack

end do
end do

end subroutine SNOW
