!-----------------------------------------------------------------------
! Snow thermodynamics and hydrology
!-----------------------------------------------------------------------
subroutine SNOW(dt,drip,Esrf,Gsrf,ksnow,ksoil,Melt,Rf,Sf,Ta,trans,     &
                Tsrf,unload,Nsnow,Dsnw,Rgrn,Sice,Sliq,Tsnow,Tsoil,     &
                Gsoil,Roff,snd,snw,Wflx)

#include "OPTS.h"

use CONSTANTS, only: &
  e0,                &! Saturation vapour pressure at Tm (Pa)
  g,                 &! Acceleration due to gravity (m/s^2)
  hcap_ice,          &! Specific heat capacity of ice (J/K/kg)
  hcap_wat,          &! Specific heat capacity of water (J/K/kg)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  mu_wat,            &! Dynamic viscosity of water (kg/m/s)
  rho_ice,           &! Density of ice (kg/m^3)
  rho_wat,           &! Density of water (kg/m^3
  Rwat,              &! Gas constant for water vapour (J/K/kg)
  Tm                  ! Melting point (n)

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
  Ta,                &! Air temperature (n)
  trans,             &! Wind-blown snow transport rate (kg/m^2/s)
  Tsrf,              &! Snow/ground surface temperature (n)
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
  Tsnow(Nsmax),      &! Snow layer temperatures (n)
  Tsoil(Nsoil)        ! Soil layer temperatures (n)

real, intent(out) :: &
  Gsoil,             &! Heat flux into soil (W/m^2)
  Roff,              &! Runoff from snow (kg/m^2/s)
  snd,               &! Snow depth (m)
  snw,               &! Total snow mass on ground (kg/m^2)
  Wflx(Nsmax)         ! Water flux into snow layer (kg/m^2/s)

integer :: &
  i,j,               &! Hydrology iteration counters
  n,                 &! Snow layer counter
  new,               &! New snow layer pointer
  Nold                ! Previous number of snow layers

real :: &
  coldcont,          &! Layer cold content (J/m^2)
  dnew,              &! New snow layer thickness (m)
  dpdT,              &! Derivative of saturation vapour density (kg/m^3/K)
  dSice,             &! Change in layer ice content (kg/m^2)
  dTdz,              &! Snow temperature gradient (K/m)
  Esnow,             &! Snow sublimation rate (kg/m^2/s)
  ggr,               &! Grain area growth rate (m^2/s)
  mass,              &! Mass of overlying snow (kg/m^2)
  qv,                &! Vertical vapour flux in snow (kg/m^2/s)
  rhos,              &! Density of snow layer (kg/m^3)
  SliqMax,           &! Maximum liquid content for layer (kg/m^2)
  Tl,                &! Snow layer lower boundary temperature (n)
  Tu,                &! Snow layer upper boundary temperature (n)
  wt                  ! Layer weighting

real :: &
  a(Nsmax),          &! Below-diagonal matrix elements
  b(Nsmax),          &! Diagonal matrix elements
  c(Nsmax),          &! Above-diagonal matrix elements
  csnow(Nsmax),      &! Areal heat capacity of snow (J/K/m^2)
  dTs(Nsmax),        &! Temperature increments (n)
  D(Nsmax),          &! Layer thickness before adjustment (m)
  E(Nsmax),          &! Energy contents before adjustment (J/m^2)
  phi(Nsmax),        &! Porosity of snow layers
  rhs(Nsmax),        &! Matrix equation rhs
  R(Nsmax),          &! Snow grain radii before adjustment (kg/m^2)
  S(Nsmax),          &! Ice contents before adjustment (kg/m^2)
  U(Nsmax),          &! Thermal transmittance between layers (W/m^2/K)
                      ! / Layer internal energy contents (J/m^2)
  W(Nsmax)            ! Liquid contents before adjustment (kg/m^2)

real :: &
  dth,               &! Hydrology timestep (s)
  dtheta(Nsmax),     &! Change in liquid water content
  ksat(Nsmax),       &! Saturated hydraulic conductivity (m/s)
  thetar(Nsmax),     &! Irreducible water content
  thetaw(Nsmax),     &! Volumetric liquid water content
  theta0(Nsmax),     &! Liquid water content at start of timestep
  Qw(Nsmax+1)         ! Water flux at snow layer boundaries (m/s)

! No snow
Gsoil = Gsrf
Roff = Rf + drip/dt
Wflx(:) = 0

! Existing snowpack
if (Nsnow > 0) then 

  ! Heat conduction
  do n = 1, Nsnow
    csnow(n) = Sice(n)*hcap_ice + Sliq(n)*hcap_wat
  end do
  if (Nsnow == 1) then
    U(1) = 2 / (Dsnw(1)/ksnow(1) + Dzsoil(1)/ksoil(1))
    dTs(1) = (Gsrf + U(1)*(Tsoil(1) - Tsnow(1)))*dt /  &
             (csnow(1) + U(1)*dt)
  else
    U(1) = 2 / (Dsnw(1)/ksnow(1) + Dsnw(2)/ksnow(2))
    a(1) = 0
    b(1) = csnow(1) + U(1)*dt
    c(1) = - U(1)*dt
    rhs(1) = (Gsrf - U(1)*(Tsnow(1) - Tsnow(2)))*dt
    do n = 2, Nsnow - 1
      U(n) = 2 / (Dsnw(n)/ksnow(n) + Dsnw(n+1)/ksnow(n+1))
      a(n) = c(n-1)
      b(n) = csnow(n) + (U(n-1) + U(n))*dt
      c(n) = - U(n)*dt
      rhs(n) = U(n-1)*(Tsnow(n-1) - Tsnow(n))*dt  &
               + U(n)*(Tsnow(n+1) - Tsnow(n))*dt 
    end do
    n = Nsnow
    U(n) = 2 / (Dsnw(n)/ksnow(n) + Dzsoil(1)/ksoil(1))
    a(n) = c(n-1)
    b(n) = csnow(n) + (U(n-1) + U(n))*dt
    c(n) = 0
    rhs(n) = U(n-1)*(Tsnow(n-1) - Tsnow(n))*dt  &
             + U(n)*(Tsoil(1) - Tsnow(n))*dt
    call TRIDIAG(Nsnow,Nsmax,a,b,c,rhs,dTs)
  end if 
  do n = 1, Nsnow
    Tsnow(n) = Tsnow(n) + dTs(n)
  end do
  n = Nsnow
  Gsoil = U(n)*(Tsnow(n) - Tsoil(1))

  ! Convert melting ice to liquid water
  dSice = Melt*dt
  do n = 1, Nsnow
    coldcont = csnow(n)*(Tm - Tsnow(n))
    if (coldcont < 0) then
      dSice = dSice - coldcont/Lf
      Tsnow(n) = Tm
    end if
    if (dSice > 0) then
      if (dSice > Sice(n)) then  ! Layer melts completely
        dSice = dSice - Sice(n)
        Dsnw(n) = 0
        Sliq(n) = Sliq(n) + Sice(n)
        Sice(n) = 0
      else                       ! Layer melts partially
        Dsnw(n) = (1 - dSice/Sice(n))*Dsnw(n)
        Sice(n) = Sice(n) - dSice
        Sliq(n) = Sliq(n) + dSice
        dSice = 0
      end if
    end if
  end do

  ! Remove snow by sublimation 
  dSice = Esrf*dt
  if (dSice > 0) then
    do n = 1, Nsnow
      if (dSice > Sice(n)) then  ! Layer sublimates completely
        dSice = dSice - Sice(n)
        Dsnw(n) = 0
        Sice(n) = 0
      else                       ! Layer sublimates partially
        Dsnw(n) = (1 - dSice/Sice(n))*Dsnw(n)
        Sice(n) = Sice(n) - dSice
        dSice = 0
      end if
    end do
  end if

  ! Remove wind-transported snow 
  dSice = trans*dt
  if (dSice > 0) then
    do n = 1, Nsnow
      if (dSice > Sice(n)) then  ! Layer completely removed
        dSice = dSice - Sice(n)
        Dsnw(n) = 0
        Sice(n) = 0
      else                       ! Layer partially removed
        Dsnw(n) = (1 - dSice/Sice(n))*Dsnw(n)
        Sice(n) = Sice(n) - dSice
        dSice = 0
      end if
    end do
  end if

#if DENSTY == 0
  ! Fixed snow density
  do n = 1, Nsnow
    if (Dsnw(n) > epsilon(Dsnw)) Dsnw(n) = (Sice(n) + Sliq(n)) / rfix
  end do
#elif DENSTY == 1
  ! Snow compaction with age
  do n = 1, Nsnow
    if (Dsnw(n) > epsilon(Dsnw)) then
      rhos = (Sice(n) + Sliq(n)) / Dsnw(n)
      if (Tsnow(n) >= Tm) then
          if (rhos < rmlt) rhos = rmlt + (rhos - rmlt)*exp(-dt/trho)
      else
          if (rhos < rcld) rhos = rcld + (rhos - rcld)*exp(-dt/trho)
      end if
      Dsnw(n) = (Sice(n) + Sliq(n)) / rhos
    end if
  end do
#elif DENSTY == 2
  ! Snow compaction by overburden
    mass = 0
    do n = 1, Nsnow
      mass = mass + 0.5*(Sice(n) + Sliq(n)) 
      if (Dsnw(n) > epsilon(Dsnw)) then
        rhos = (Sice(n) + Sliq(n)) / Dsnw(n)
        rhos = rhos + (rhos*g*mass*dt/(eta0*exp(-(Tsnow(n) - Tm)/12.4 + rhos/55.6))  &
               + dt*rhos*snda*exp((Tsnow(n) - Tm)/23.8 - max(rhos - 150, 0.)/21.7))
        Dsnw(n) = (Sice(n) + Sliq(n)) / rhos
      end if
      mass = mass + 0.5*(Sice(n) + Sliq(n))
    end do
#else
    stop 'Unknown option DENSTY'
#endif

#if SGRAIN == 0
  ! Snow grain growth not represented
  Rgrn(:) = rgr0
#elif SGRAIN == 1
  ! Temperature dependent snow grain growth
  do n = 1, Nsnow
    ggr = 2e-13
    if (Tsnow(n) < Tm) then
      if (Rgrn(n) < 1.50e-4) then
        ggr = 2e-14
      else
        ggr = 7.3e-8*exp(-4600/Tsnow(n))
      end if
    end if
    Rgrn(n) = Rgrn(n) + dt*ggr/Rgrn(n)
  end do
#elif SGRAIN == 2
  ! Temperature gradient dependent snow grain growth
  do n = 1, Nsnow
    Tu = Tsrf
    if (n > 1) Tu = (Dsnw(n-1)*Tsnow(n) + Dsnw(n)*Tsnow(n-1)) / (Dsnw(n) + Dsnw(n-1))
    Tl = (Dzsoil(1)*Tsnow(n) + Dsnw(n)*Tsoil(1)) / (Dsnw(n) + Dzsoil(1))
    if (n < Nsnow) Tl = (Dsnw(n+1)*Tsnow(n) + Dsnw(n)*Tsnow(n+1)) / (Dsnw(n) + Dsnw(n+1))
    dTdz = abs(Tu - Tl) / Dsnw(n)
    thetaw(n) = Sliq(n)/(rho_wat*Dsnw(n))
    if (thetaw(n) < 1e-4) then
      dpdT = (e0/(Rwat*Tsnow(n)**2))*(Ls/(Rwat*Tsnow(n)) - 1)*exp((Ls/Rwat)*(1/Tm - 1/Tsnow(n)))
      qv = 9.2e-5*(Tsnow(n)/Tm)**6*dpdT*dTdz
      ggr = 1.25e-7*min(qv, 1e-6)
    else
      ggr = 1e-12*min(thetaw(n) + 0.05, 0.14)
    end if
    Rgrn(n) = Rgrn(n) + dt*ggr/Rgrn(n)
  end do
#else
  stop 'Unknown option SGRAIN'
#endif

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
do n = 1, Nsnow
  csnow(n) = Sice(n)*hcap_ice + Sliq(n)*hcap_wat
  E(n) = csnow(n)*(Tsnow(n) - Tm)
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
  n = 1
  if (Dsnw(1) > Dzsnow(1)) then 
    do n = 1, Nsmax
      Dsnw(n) = Dzsnow(n)
      dnew = dnew - Dzsnow(n)
      if (dnew <= Dzsnow(n) .or. n == Nsmax) then
        Dsnw(n) = Dsnw(n) + dnew
        exit
      end if
    end do
  end if
  Nsnow = n

  ! Fill new layers from the top downwards
  new = 1
  dnew = Dsnw(1)
  do n = 1, Nold
    do
      if (D(n) < dnew) then
        ! All snow from old layer partially fills new layer
        Rgrn(new) = Rgrn(new) + S(n)*R(n)
        Sice(new) = Sice(new) + S(n)
        Sliq(new) = Sliq(new) + W(n)
        U(new) = U(new) + E(n)
        dnew = dnew - D(n)
        exit
      else
        ! Some snow from old layer fills new layer
        wt = dnew / D(n)
        Rgrn(new) = Rgrn(new) + wt*S(n)*R(n)
        Sice(new) = Sice(new) + wt*S(n) 
        Sliq(new) = Sliq(new) + wt*W(n)
        U(new) = U(new) + wt*E(n)
        D(n) = (1 - wt)*D(n)
        E(n) = (1 - wt)*E(n)
        S(n) = (1 - wt)*S(n)
        W(n) = (1 - wt)*W(n)
        new = new + 1
        if (new > Nsnow) exit
        dnew = Dsnw(new)
      end if
    end do
  end do

  ! Diagnose snow layer temperatures
  do n = 1, Nsnow
    csnow(n) = Sice(n)*hcap_ice + Sliq(n)*hcap_wat
    Tsnow(n) = Tm + U(n) / csnow(n)
    Rgrn(n) = Rgrn(n) / Sice(n)
  end do

  ! Drain, retain or freeze snow in layers
#if HYDROL == 0
  ! Free-draining snow, no retention or freezing 
  Wflx(1) = Roff
  do n = 1, Nsnow
    Roff = Roff + Sliq(n) / dt
    Sliq(n) = 0
    if (n < Nsnow) Wflx(n+1) = Roff
  end do
#elif HYDROL == 1
  ! Bucket storage 
  if (maxval(Sliq)>0 .or. Rf>0) then
  do n = 1, Nsnow
    phi(n) = 1 - Sice(n)/(rho_ice*Dsnw(n))
    if (phi(n) < 0) phi(n) = 0
    SliqMax = rho_wat*Dsnw(n)*phi(n)*Wirr
    Sliq(n) = Sliq(n) + Roff*dt
    Wflx(n) = Roff
    Roff = 0
    if (Sliq(n) > SliqMax) then       ! Liquid capacity exceeded
      Roff = (Sliq(n) - SliqMax)/dt   ! so drainage to next layer
      Sliq(n) = SliqMax
    end if
    csnow(n) = Sice(n)*hcap_ice + Sliq(n)*hcap_wat
    coldcont = csnow(n)*(Tm - Tsnow(n))
    if (coldcont > 0) then            ! Liquid can freeze
      dSice = min(Sliq(n), coldcont/Lf)
      Sliq(n) = Sliq(n) - dSice
      Sice(n) = Sice(n) + dSice
      Tsnow(n) = Tsnow(n) + Lf*dSice/csnow(n)
    end if
  end do
  end if
#elif HYDROL == 2
  ! Gravitational drainage 
  if (maxval(Sliq)>0 .or. Rf>0) then
  Qw(:) = 0
  Qw(1) = Rf/rho_wat
  Roff = 0
  thetaw(:) = 0
  do n = 1, Nsnow
    ksat(n) = 0.31*(rho_wat*g/mu_wat)*Rgrn(n)**2*exp(-7.8*Sice(n)/(rho_wat*Dsnw(n)))
    phi(n) = 1 - Sice(n)/(rho_ice*Dsnw(n))
    thetar(n) = Wirr*phi(n)
    thetaw(n) = Sliq(n)/(rho_wat*Dsnw(n))
    if (thetaw(n)>phi(n)) then
      Roff = Roff + rho_wat*Dsnw(n)*(thetaw(n) - phi(n))/dt
      thetaw(n) = phi(n)
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
      do n = 2, Nsnow
        if (thetaw(n-1) > thetar(n-1))  &
          a(n) = - 3*ksat(n-1)*(thetaw(n-1) - thetar(n-1))**2/(phi(n-1) - thetar(n-1))**3/Dsnw(n-1)
        if (thetaw(n) > thetar(n)) then
          b(n) = 1/dth + 3*ksat(n)*(thetaw(n) - thetar(n))**2/(phi(n) - thetar(n))**3/Dsnw(n)
          Qw(n+1) = ksat(n)*((thetaw(n) - thetar(n))/(phi(n) - thetar(n)))**3
        end if
        rhs(n) = (thetaw(n) - theta0(n))/dth + (Qw(n+1) - Qw(n))/Dsnw(n)
      end do
      dtheta(1) = - rhs(1)/b(1)
      do n = 2, Nsnow
        dtheta(n) = - (a(n)*dtheta(n-1) + rhs(n))/b(n)
      end do 
      do n = 1, Nsnow
        thetaw(n) = max(thetaw(n) + dtheta(n), 0.)
        if (thetaw(n) > phi(n)) then
          Qw(n+1) = Qw(n+1) + (thetaw(n) - phi(n))*Dsnw(n)/dth
          thetaw(n) = phi(n)
        endif
      end do
    end do
    Wflx(:) = Wflx(:) + rho_wat*Qw(1:Nsmax)/10
    Roff = Roff + rho_wat*Qw(Nsnow+1)/10
  end do
  Sliq(:) = rho_wat*Dsnw(:)*thetaw(:)
  do n = 1, Nsnow
    csnow(n) = Sice(n)*hcap_ice + Sliq(n)*hcap_wat
    coldcont = csnow(n)*(Tm - Tsnow(n))
    if (coldcont > 0) then            ! Liquid can freeze
      dSice = min(Sliq(n), coldcont/Lf)
      Sliq(n) = Sliq(n) - dSice
      Sice(n) = Sice(n) + dSice
      Tsnow(n) = Tsnow(n) + Lf*dSice/csnow(n)
    end if
  end do
  end if
#else
  stop 'Unknown option HYDROL'
#endif
end if ! Existing or new snowpack
snw = sum(Sice(:)) + sum(Sliq(:))

end subroutine SNOW
