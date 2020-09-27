!-----------------------------------------------------------------------
! Snow thermodynamics and hydrology
!-----------------------------------------------------------------------
subroutine SNOW(dt,drip,Esrf,Gsrf,ksnow,ksoil,Melt,Rf,Sf,Ta,trans,     &
                Tsrf,unload,Nsnow,Dsnw,Rgrn,Sice,Sliq,Tsnow,Tsoil,     &
                Gsoil,Roff,snd,snw,Wflx)

#include "OPTS.h"

use CONSTANTS, only: &
  g,                 &! Acceleration due to gravity (m/s^2)
  hcap_ice,          &! Specific heat capacity of ice (J/K/kg)
  hcap_wat,          &! Specific heat capacity of water (J/K/kg)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  mu_wat,            &! Dynamic viscosity of water (kg/m/s)
  rho_ice,           &! Density of ice (kg/m^3)
  rho_wat,           &! Density of water (kg/m^3)
  Tm                  ! Melting point (K)

use LAYERS, only: &
  Dzsnow,            &! Minimum snow layer thicknesses (m)
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsmax,             &! Maximum number of snow layers
  Nsoil               ! Number of soil layers

use PARAMETERS, only: &
  eta0,              &! Reference snow viscosity (Pa s)
  rcld,              &! Maximum density for cold snow (kg/m^3)
  rfix,              &! Fixed snow density (kg/m^3)
  rgr0,              &! Fresh snow grain radius (m)
  rhof,              &! Fresh snow density (kg/m^3)
  rhow,              &! Wind-packed snow density (kg/m^3)
  rmlt,              &! Maximum density for melting snow (kg/m^3)
  snda,              &! Thermal metamorphism parameter (1/s)
  trho,              &! Snow compaction timescale (s)
  Wirr                ! Irreducible liquid water content of snow

implicit none

real, intent(in) :: &
  dt,                &! Timestep (s)
  drip,              &! Melt water drip from vegetation (kg/m^2)
  Esrf,              &! Moisture flux from the surface (kg/m^2/s)
  Gsrf,              &! Heat flux into snow/ground surface (W/m^2)
  Melt,              &! Surface melt rate (kg/m^2/s)
  Rf,                &! Rainfall rate (kg/m^2/s)
  Sf,                &! Snowfall rate (kg/m^2/s)
  Ta,                &! Air temperature (K)
  trans,             &! Wind-blown snow transport rate (kg/m^2/s)
  Tsrf,              &! Snow/ground surface temperature (K)
  unload,            &! Snow mass unloaded from vegetation (kg/m^2)
  ksnow(Nsmax),      &! Thermal conductivity of snow layers (W/m/K)
  ksoil(Nsoil)        ! Thermal conductivity of soil layers (W/m/K)

integer, intent(inout) :: &
  Nsnow               ! Number of snow layers

real, intent(inout) :: &
  Dsnw(Nsmax),       &! Snow layer thicknesses (m)
  Rgrn(Nsmax),       &! Snow layer grain radius (m)
  Sice(Nsmax),       &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmax),       &! Liquid content of snow layers (kg/m^2)
  Tsnow(Nsmax),      &! Snow layer temperatures (K)
  Tsoil(Nsoil)        ! Soil layer temperatures (K)

real, intent(out) :: &
  Gsoil,             &! Heat flux into soil (W/m^2)
  Roff,              &! Runoff from snow (kg/m^2/s)
  snd,               &! Snow depth (m)
  snw,               &! Total snow mass on ground (kg/m^2)
  Wflx(Nsmax)         ! Water flux into snow layer (kg/m^2/s)

integer :: &
  i,j,               &! Hydrology iteration counters
  k,                 &! Snow layer counter
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
  rhos,              &! Density of snow layer (kg/m^3)
  SliqMax,           &! Maximum liquid content for layer (kg/m^2)
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
  phi(Nsmax),        &! Porosity of snow layers
  rhs(Nsmax),        &! Matrix equation rhs
  R(Nsmax),          &! Snow grain radii before adjustment (kg/m^2)
  S(Nsmax),          &! Ice contents before adjustment (kg/m^2)
  U(Nsmax),          &! Layer internal energy contents (J/m^2)
  W(Nsmax)            ! Liquid contents before adjustment (kg/m^2)

real :: &
  dth,               &! Hydrology timestep (s)
  dtheta(Nsmax),     &! Change in liquid water content
  ksat(Nsmax),       &! Saturated hydraulic conductivity (m/s)
  thetar(Nsmax),     &! Irreducible water content
  thetaw(Nsmax),     &! Volumetric liquid water content
  theta0(Nsmax),     &! Liquid water content at start of timestep
  Qw(Nsmax+1)         ! Water flux at snow layer boundaruess (m/s)

! No snow
Gsoil = Gsrf
Roff = Rf + drip/dt
Wflx(:) = 0

! Existing snowpack
if (Nsnow > 0) then 

  ! Heat conduction
  do k = 1, Nsnow
    csnow(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
  end do
  if (Nsnow == 1) then
    Gs(1) = 2 / (Dsnw(1)/ksnow(1) + Dzsoil(1)/ksoil(1))
    dTs(1) = (Gsrf + Gs(1)*(Tsoil(1) - Tsnow(1)))*dt /  &
             (csnow(1) + Gs(1)*dt)
  else
    do k = 1, Nsnow - 1
      Gs(k) = 2 / (Dsnw(k)/ksnow(k) + Dsnw(k+1)/ksnow(k+1))
    end do
    a(1) = 0
    b(1) = csnow(1) + Gs(1)*dt
    c(1) = - Gs(1)*dt
    rhs(1) = (Gsrf - Gs(1)*(Tsnow(1) - Tsnow(2)))*dt
    do k = 2, Nsnow - 1
      a(k) = c(k-1)
      b(k) = csnow(k) + (Gs(k-1) + Gs(k))*dt
      c(k) = - Gs(k)*dt
      rhs(k) = Gs(k-1)*(Tsnow(k-1) - Tsnow(k))*dt  &
               + Gs(k)*(Tsnow(k+1) - Tsnow(k))*dt 
    end do
    k = Nsnow
    Gs(k) = 2 / (Dsnw(k)/ksnow(k) + Dzsoil(1)/ksoil(1))
    a(k) = c(k-1)
    b(k) = csnow(k) + (Gs(k-1) + Gs(k))*dt
    c(k) = 0
    rhs(k) = Gs(k-1)*(Tsnow(k-1) - Tsnow(k))*dt  &
             + Gs(k)*(Tsoil(1) - Tsnow(k))*dt
    call TRIDIAG(Nsnow,Nsmax,a,b,c,rhs,dTs)
  end if 
  do k = 1, Nsnow
    Tsnow(k) = Tsnow(k) + dTs(k)
  end do
  k = Nsnow
  Gsoil = Gs(k)*(Tsnow(k) - Tsoil(1))

  ! Convert melting ice to liquid water
  dSice = Melt*dt
  do k = 1, Nsnow
    coldcont = csnow(k)*(Tm - Tsnow(k))
    if (coldcont < 0) then
      dSice = dSice - coldcont/Lf
      Tsnow(k) = Tm
    end if
    if (dSice > 0) then
      if (dSice > Sice(k)) then  ! Layer melts completely
        dSice = dSice - Sice(k)
        Dsnw(k) = 0
        Sliq(k) = Sliq(k) + Sice(k)
        Sice(k) = 0
      else                       ! Layer melts partially
        Dsnw(k) = (1 - dSice/Sice(k))*Dsnw(k)
        Sice(k) = Sice(k) - dSice
        Sliq(k) = Sliq(k) + dSice
        dSice = 0
      end if
    end if
  end do

  ! Remove snow by sublimation 
  dSice = Esrf*dt
  if (dSice > 0) then
    do k = 1, Nsnow
      if (dSice > Sice(k)) then  ! Layer sublimates completely
        dSice = dSice - Sice(k)
        Dsnw(k) = 0
        Sice(k) = 0
      else                       ! Layer sublimates partially
        Dsnw(k) = (1 - dSice/Sice(k))*Dsnw(k)
        Sice(k) = Sice(k) - dSice
        dSice = 0
      end if
    end do
  end if

  ! Remove wind-trasported snow 
  dSice = trans*dt
  if (dSice > 0) then
    do k = 1, Nsnow
      if (dSice > Sice(k)) then  ! Layer completely removed
        dSice = dSice - Sice(k)
        Dsnw(k) = 0
        Sice(k) = 0
      else                       ! Layer partially removed
        Dsnw(k) = (1 - dSice/Sice(k))*Dsnw(k)
        Sice(k) = Sice(k) - dSice
        dSice = 0
      end if
    end do
  end if

#if DENSTY == 0
  ! Fixed snow density
  do k = 1, Nsnow
    if (Dsnw(k) > epsilon(Dsnw)) Dsnw(k) = (Sice(k) + Sliq(k)) / rfix
  end do
#endif
#if DENSTY == 1
  ! Snow compaction with age
  do k = 1, Nsnow
    if (Dsnw(k) > epsilon(Dsnw)) then
      rhos = (Sice(k) + Sliq(k)) / Dsnw(k)
      if (Tsnow(k) >= Tm) then
          if (rhos < rmlt) rhos = rmlt + (rhos - rmlt)*exp(-dt/trho)
      else
          if (rhos < rcld) rhos = rcld + (rhos - rcld)*exp(-dt/trho)
      end if
      Dsnw(k) = (Sice(k) + Sliq(k)) / rhos
    end if
  end do
#endif
#if DENSTY == 2
  ! Snow compaction by overburden
    mass = 0
    do k = 1, Nsnow
      mass = mass + 0.5*(Sice(k) + Sliq(k)) 
      if (Dsnw(k) > epsilon(Dsnw)) then
        rhos = (Sice(k) + Sliq(k)) / Dsnw(k)
        rhos = rhos + (rhos*g*mass*dt/(eta0*exp(-(Tsnow(k) - Tm)/12.4 + rhos/55.6))  &
               + dt*rhos*snda*exp((Tsnow(k) - Tm)/23.8 - max(rhos - 150, 0.)/21.7))
        Dsnw(k) = (Sice(k) + Sliq(k)) / rhos
      end if
      mass = mass + 0.5*(Sice(k) + Sliq(k))
    end do
#endif

  ! Snow grain growth
  do k = 1, Nsnow
    ggr = 2e-13
    if (Tsnow(k) < Tm) then
      if (Rgrn(k) < 1.50e-4) then
        ggr = 2e-14
      else
        ggr = 7.3e-8*exp(-4600/Tsnow(k))
      end if
    end if
    Rgrn(k) = Rgrn(k) + dt*ggr/Rgrn(k)
  end do

end if  ! Existing snowpack

! Add snowfall and frost to layer 1 with fresh snow density and grain size
Esnow = 0
if (Esrf < 0 .and. Tsrf < Tm) Esnow = Esrf
dSice = (Sf - Esnow)*dt
Dsnw(1) = Dsnw(1) + dSice / rhof
if (Sice(1) + dSice > epsilon(Sice)) then
  Rgrn(1) = (Sice(1)*Rgrn(1) + dSice*rgr0) / (Sice(1) + dSice)
end if
Sice(1) = Sice(1) + dSice

! Add canopy unloading to layer 1 with bulk snow density and grain size
rhos = rhof
snw = sum(Sice(:)) + sum(Sliq(:))
snd = sum(Dsnw(:))
if (snd > epsilon(snd)) rhos = snw / snd
Dsnw(1) = Dsnw(1) + unload / rhos
if (Sice(1) + unload > epsilon(Sice)) then
  Rgrn(1) = (Sice(1)*Rgrn(1) + unload*rgr0) / (Sice(1) + unload)
end if
Sice(1) = Sice(1) + unload

! Add wind-blown snow to layer 1 with wind-packed density and fresh grain size
dSice = - trans*dt
if (dSice > 0) then
  Dsnw(1) = Dsnw(1) + dSice / rhow
  Rgrn(1) = (Sice(1)*Rgrn(1) + dSice*rgr0) / (Sice(1) + dSice)
  Sice(1) = Sice(1) + dSice
end if

! New snowpack
if (Nsnow == 0 .and. Sice(1) > 0) then
  Nsnow = 1
  Rgrn(1) = rgr0
  Tsnow(1) = min(Ta, Tm)
end if

! Store state of old layers
D(:) = Dsnw(:)
R(:) = Rgrn(:)
S(:) = Sice(:)
W(:) = Sliq(:)
do k = 1, Nsnow
  csnow(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
  E(k) = csnow(k)*(Tsnow(k) - Tm)
end do
Nold = Nsnow
snd = sum(Dsnw(:))

! Initialise new layers
Dsnw(:) = 0
Rgrn(:) = 0
Sice(:) = 0
Sliq(:) = 0
Tsnow(:) = Tm
U(:) = 0
Nsnow = 0

if (snd > 0) then  ! Existing or new snowpack
 
  ! Re-assign and count snow layers
  dnew = snd
  Dsnw(1) = dnew
  k = 1
  if (Dsnw(1) > Dzsnow(1)) then 
    do k = 1, Nsmax
      Dsnw(k) = Dzsnow(k)
      dnew = dnew - Dzsnow(k)
      if (dnew <= Dzsnow(k) .or. k == Nsmax) then
        Dsnw(k) = Dsnw(k) + dnew
        exit
      end if
    end do
  end if
  Nsnow = k

  ! Fill new layers from the top downwards
  knew = 1
  dnew = Dsnw(1)
  do kold = 1, Nold
    do
      if (D(kold) < dnew) then
        ! All snow from old layer partially fills new layer
        Rgrn(knew) = Rgrn(knew) + S(kold)*R(kold)
        Sice(knew) = Sice(knew) + S(kold)
        Sliq(knew) = Sliq(knew) + W(kold)
        U(knew) = U(knew) + E(kold)
        dnew = dnew - D(kold)
        exit
      else
        ! Some snow from old layer fills new layer
        wt = dnew / D(kold)
        Rgrn(knew) = Rgrn(knew) + wt*S(kold)*R(kold)
        Sice(knew) = Sice(knew) + wt*S(kold) 
        Sliq(knew) = Sliq(knew) + wt*W(kold)
        U(knew) = U(knew) + wt*E(kold)
        D(kold) = (1 - wt)*D(kold)
        E(kold) = (1 - wt)*E(kold)
        S(kold) = (1 - wt)*S(kold)
        W(kold) = (1 - wt)*W(kold)
        knew = knew + 1
        if (knew > Nsnow) exit
        dnew = Dsnw(knew)
      end if
    end do
  end do

  ! Diagnose snow layer temperatures
  do k = 1, Nsnow
    csnow(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
    Tsnow(k) = Tm + U(k) / csnow(k)
    Rgrn(k) = Rgrn(k) / Sice(k)
  end do

  ! Drain, retain or freeze snow in layers
#if HYDROL == 0
  ! Free-draining snow, no retention or freezing 
  Wflx(1) = Roff
  do k = 1, Nsnow
    Roff = Roff + Sliq(k) / dt
    Sliq(k) = 0
    if (k < Nsnow) Wflx(k+1) = Roff
  end do
#endif
#if HYDROL == 1
  ! Bucket storage 
  if (maxval(Sliq)>0 .or. Rf>0) then
  do k = 1, Nsnow
    phi(k) = 1 - Sice(k)/(rho_ice*Dsnw(k))
    SliqMax = rho_wat*Dsnw(k)*phi(k)*Wirr
    Sliq(k) = Sliq(k) + Roff*dt
    Wflx(k) = Roff
    Roff = 0
    if (Sliq(k) > SliqMax) then       ! Liquid capacity exceeded
      Roff = (Sliq(k) - SliqMax)/dt   ! so drainage to next layer
      Sliq(k) = SliqMax
    end if
    csnow(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
    coldcont = csnow(k)*(Tm - Tsnow(k))
    if (coldcont > 0) then            ! Liquid can freeze
      dSice = min(Sliq(k), coldcont/Lf)
      Sliq(k) = Sliq(k) - dSice
      Sice(k) = Sice(k) + dSice
      Tsnow(k) = Tsnow(k) + Lf*dSice/csnow(k)
    end if
  end do
  end if
#endif
#if HYDROL == 2
  ! Gravitational drainage 
  if (maxval(Sliq)>0 .or. Rf>0) then
  Qw(:) = 0
  Qw(1) = Rf/rho_wat
  Roff = 0
  do k = 1, Nsnow
    ksat(k) = 0.31*(rho_wat*g/mu_wat)*Rgrn(k)**2*exp(-7.8*Sice(k)/(rho_wat*Dsnw(k)))
    phi(k) = 1 - Sice(k)/(rho_ice*Dsnw(k))
    thetar(k) = Wirr*phi(k)
    thetaw(k) = Sliq(k)/(rho_wat*Dsnw(k))
    if (thetaw(k)>phi(k)) then
      Roff = Roff + rho_wat*Dsnw(k)*(thetaw(k) - phi(k))/dt
      thetaw(k) = phi(k)
    end if
  end do
  dth = 0.1*dt
  do i = 1, 10  ! subdivide timestep
    theta0(:) = thetaw(:)
    do j = 1, 10  ! Newton-Raphson iteration
      a(:) = 0
      b(:) = 1/dth 
      if (thetaw(1) > thetar(1)) then
        b(1) = 1/dth + 3*ksat(1)*(thetaw(1) - thetar(1))**2/(phi(1) - thetar(1))**3/Dsnw(1)
        Qw(2) = ksat(1)*((thetaw(1) - thetar(1))/(phi(1) - thetar(1)))**3
      end if
      rhs(1) = (thetaw(1) - theta0(1))/dth + (Qw(2) - Qw(1))/Dsnw(1)
      do k = 2, Nsnow
        if (thetaw(k-1) > thetar(k-1))  &
          a(k) = - 3*ksat(k-1)*(thetaw(k-1) - thetar(k-1))**2/(phi(k-1) - thetar(k-1))**3/Dsnw(k-1)
        if (thetaw(k) > thetar(k)) then
          b(k) = 1/dth + 3*ksat(k)*(thetaw(k) - thetar(k))**2/(phi(k) - thetar(k))**3/Dsnw(k)
          Qw(k+1) = ksat(k)*((thetaw(k) - thetar(k))/(phi(k) - thetar(k)))**3
        end if
        rhs(k) = (thetaw(k) - theta0(k))/dth + (Qw(k+1) - Qw(k))/Dsnw(k)
      end do
      dtheta(1) = - rhs(1)/b(1)
      do k = 2, Nsnow
        dtheta(k) = - (a(k)*dtheta(k-1) + rhs(k))/b(k)
      end do 
      do k = 1, Nsnow
        thetaw(k) = thetaw(k) + dtheta(k)
        if (thetaw(k) > phi(k)) then
          Qw(k+1) = Qw(k+1) + (thetaw(k) - phi(k))*Dsnw(k)/dth
          thetaw(k) = phi(k)
        endif
      end do
    end do
    Wflx(:) = Wflx(:) + rho_wat*Qw(1:Nsmax)/10
    Roff = Roff + rho_wat*Qw(Nsnow+1)/10
  end do
  Sliq(:) = rho_wat*Dsnw(:)*thetaw(:)
  do k = 1, Nsnow
    csnow(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
    coldcont = csnow(k)*(Tm - Tsnow(k))
    if (coldcont > 0) then            ! Liquid can freeze
      dSice = min(Sliq(k), coldcont/Lf)
      Sliq(k) = Sliq(k) - dSice
      Sice(k) = Sice(k) + dSice
      Tsnow(k) = Tsnow(k) + Lf*dSice/csnow(k)
    end if
  end do
  end if
#endif
snw = sum(Sice(:)) + sum(Sliq(:))
end if ! Existing or new snowpack

end subroutine SNOW
