!-----------------------------------------------------------------------
! Mass balance of snow intercepted by vegetation
!-----------------------------------------------------------------------
subroutine INTERCEPT(dt,cveg,Eveg,lveg,Scap,Ta,Ua,Sf,Sveg,Tveg,        &
                     drip,svg,unload)
                     
#include "OPTS.h"

use CONSTANTS, only: &
  Lf,                &! Latent heat of fusion (J/kg)
  Tm                  ! Melting point (K)

use LAYERS, only: &
  Ncnpy,             &! Number of canopy layers
  fvg1                ! Fraction of vegetation in upper canopy layer

use PARAMETERS, only: &
  eunl,              &! Exponential unloading time scale (s)
  kext,              &! Vegetation light extinction coefficient
  munl,              &! Melt unloading fraction
  Tunl,              &! Temperature unloading parameter (K s)
  Uunl                ! Wind unloading parameter (m)

implicit none

real, intent(in) :: &
  dt,                &! Timestep (s)
  cveg(Ncnpy),       &! Vegetation heat capacities (J/K/m^2)
  Eveg(Ncnpy),       &! Moisture flux from vegetation layers (kg/m^2/s)
  lveg(Ncnpy),       &! Canopy layer vegetation area indices
  Scap(Ncnpy),       &! Vegetation layer snow capacities (kg/m^2)
  Ta,                &! Air temperature (K)
  Ua                  ! Wind speed (m/s)

real, intent(inout) :: &
  Sf,                &! Snowfall rate (kg/m2/s)
  Sveg(Ncnpy),       &! Snow mass on vegetation layers (kg/m^2)
  Tveg(Ncnpy)         ! Vegetation layer temperatures (K)

real, intent(out) :: &
  drip,              &! Melt water drip from vegetation (kg/m^2)
  svg,               &! Total snow mass on vegetation (kg/m^2)
  unload              ! Snow mass unloaded from vegetation (kg/m^2)

integer :: & 
  n                   ! Canopy layer counter
  
real :: &
  dsvg,              &! Change in canopy snow mass (kg/m^2)
  fveg,              &! Vegetation fraction
  melt                ! Canopy snow melt (kg/m^2)

drip = 0
svg = 0
unload = 0
if (lveg(1) > epsilon(lveg)) then
  do n = 1, Ncnpy
  
    ! Interception of falling snow 
    fveg = 1 - exp(-kext*lveg(n)) 
#if CANINT == 1
    dsvg = fveg*Sf*dt
#elif CANINT == 2
    dsvg = (Scap(n) - Sveg(n))*(1 - exp(-fveg*Sf*dt/Scap(n)))
#else
    stop 'Unknown option CANINT'
#endif
    if (Sveg(n) + dsvg > Scap(n)) dsvg = Scap(n) - Sveg(n)
    Sveg(n) = Sveg(n) + dsvg
    Sf = Sf - dsvg/dt 
    
    ! Sublimation of canopy snow
    if (Eveg(n) > 0) then
      if (Sveg(n) > 0) Sveg(n) = Sveg(n) - Eveg(n)*dt
      if (Sveg(n) < 0) Sveg(n) = 0
    else
      if (Tveg(n) < Tm) then  ! frost on canopy
        Sveg(n) = Sveg(n) - Eveg(n)*dt
        if (Sveg(n) > Scap(n)) then
          unload = unload + Sveg(n) - Scap(n)
          Sveg(n) = Scap(n)
        end if
      end if
    end if

    ! Melting of canopy snow
    melt = 0
    if (Tveg(n) > Tm) then
      melt = cveg(n)*(Tveg(n) - Tm)/Lf
      if (melt > Sveg(n)) melt = Sveg(n)
      drip = drip + melt
      Sveg(n) = Sveg(n) - melt
      Tveg(n) = Tveg(n) - Lf*melt/cveg(n)
    end if
 
    ! Unloading of canopy snow   
#if CANUNL == 1
    dsvg = Sveg(n)*dt/eunl + munl*melt
#elif CANUNL == 2 
    dsvg = (max(Tveg(n) - 270.15, 0.)/Tunl + Ua/Uunl)*dt*Sveg(n)
#else
    stop 'Unknown option CANUNL'
#endif
    if (dsvg > Sveg(n)) dsvg = Sveg(n)
    Sveg(n) = Sveg(n) - dsvg
    unload = unload + dsvg
    
    Sveg(n) = min(max(Sveg(n),0.),Scap(n))
    svg = svg + Sveg(n)
  end do
end if

end subroutine INTERCEPT
