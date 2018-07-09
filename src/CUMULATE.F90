!-----------------------------------------------------------------------
! Cumulate diagnostics
!-----------------------------------------------------------------------
subroutine CUMULATE(alb,G,H,Hsrf,LE,LEsrf,Melt,Rnet,Roff,Rsrf)

use CONSTANTS, only: &
  Tm                  ! Melting point (K)

use DIAGNOSTICS, only: &
  Nave,              &! Number of timesteps in average outputs
  diags,             &! Cumulated diagnostics
  SWin,              &! Cumulated incoming solar radiation (J/m^2)
  SWout               ! Cumulated reflected solar radiation (J/m^2)

use DRIVING, only: &
  dt,                &! Timestep (s)
  SW                  ! Incoming shortwave radiation (W/m^2)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use STATE_VARIABLES, only: &
  Tsoil,             &! Soil layer temperatures (K)
  Tsrf                ! Surface skin temperature (K)

implicit none

real, intent(in) :: &
  alb(Nx,Ny),        &! Albedo
  G(Nx,Ny),          &! Heat flux into surface (W/m^2)
  H(Nx,Ny),          &! Sensible heat flux to the atmosphere (W/m^2)
  Hsrf(Nx,Ny),       &! Sensible heat flux from the surface (W/m^2)
  LE(Nx,Ny),         &! Latent heat flux to the atmosphere (W/m^2)
  LEsrf(Nx,Ny),      &! Latent heat flux from the surface (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny),       &! Net radiation (W/m^2)
  Roff(Nx,Ny),       &! Runoff from snow (kg/m^2)
  Rsrf(Nx,Ny)         ! Net radiation absorbed by the surface (W/m^2)

SWin(:,:) = SWin (:,:)+ SW(:,:)*dt
SWout(:,:) = SWout(:,:) + alb(:,:)*SW(:,:)*dt
diags(:,:,1) = diags(:,:,1) + G(:,:)
diags(:,:,2) = diags(:,:,2) + H(:,:)
diags(:,:,3) = diags(:,:,3) + Hsrf(:,:)
diags(:,:,4) = diags(:,:,4) + LE(:,:)
diags(:,:,5) = diags(:,:,5) + LEsrf(:,:)
diags(:,:,6) = diags(:,:,6) + Melt(:,:) * dt * Nave
diags(:,:,7) = diags(:,:,7) + Rnet(:,:)
diags(:,:,8) = diags(:,:,8) + Roff(:,:) * Nave
diags(:,:,9) = diags(:,:,9) + Rsrf(:,:)
diags(:,:,10) = diags(:,:,10) + Tsrf(:,:) - Tm
diags(:,:,11) = diags(:,:,11) + Tsoil(2,:,:) - Tm

end subroutine CUMULATE
