!-----------------------------------------------------------------------
! Mass balance of snow intercepted by vegetation
!-----------------------------------------------------------------------
subroutine INTERCEPT(dt,cveg,Eveg,Scap,Sf,Sveg,Tveg,drip,svg,unload)

use CONSTANTS, only: &
  Lf,                &! Latent heat of fusion (J/kg)
  Tm                  ! Melting point (K)

use LAYERS, only: &
  Ncnpy               ! Number of canopy layers

use PARAMETERS, only: &
  tunl                ! Canopy snow unloading time scale (s)

implicit none

real, intent(in) :: &
  dt,                &! Timestep (s)
  cveg(Ncnpy),       &! Vegetation heat capacities (J/K/m^2)
  Eveg(Ncnpy),       &! Moisture flux from vegetation layers (kg/m^2/s)
  Scap(Ncnpy)         ! Vegetation layer snow capacities (kg/m^2)

real, intent(inout) :: &
  Sf,                &! Snowfall rate (kg/m2/s)
  Sveg(Ncnpy),       &! Snow mass on vegetation layers (kg/m^2)
  Tveg(Ncnpy)         ! Vegetation layer temperatures (K)

real, intent(out) :: &
  drip,              &! Melt water drip from vegetation (kg/m^2)
  svg,               &! Total snow mass on vegetation (kg/m^2)
  unload              ! Snow mass unloaded from vegetation (kg/m^2)

integer :: k          ! Vegetation layer counter

real :: &
  intcpt,            &! Vegetation layer snow interception (kg/m^2)
  melt                ! Vegetation layer snow melt (kg/m^2)

drip = 0
svg = 0
unload = 0
if (Scap(1) > 0) then
  do k = 1, Ncnpy

    ! Interception of falling snow
    intcpt = (Scap(k) - Sveg(k))*(1 - exp(-Sf*dt/Scap(k)))
    Sveg(k) = Sveg(k) + intcpt
    Sf = Sf - intcpt/dt

    ! Sublimation of canopy snow
    if (Eveg(k)>0) then
      if (Sveg(k)>0) Sveg(k) = Sveg(k) - Eveg(k)*dt
    else
      if (Tveg(k)<Tm) Sveg(k) = Sveg(k) - Eveg(k)*dt
    end if
    Sveg(k) = max(Sveg(k), 0.)

    ! Unloading of canopy snow
    unload = unload + Sveg(k)*dt/tunl
    Sveg(k) = (1 - dt/tunl)*Sveg(k)

    ! Melting of canopy snow
    if (Tveg(k) > Tm) then
      melt = cveg(k)*(Tveg(k) - Tm)/Lf
      if (melt > Sveg(k)) melt = Sveg(k)
      drip = drip + melt
      Sveg(k) = Sveg(k) - melt
      Tveg(k) = Tveg(k) - Lf*melt/cveg(k)
    end if

  end do
end if
svg = sum(Sveg)

end subroutine INTERCEPT
