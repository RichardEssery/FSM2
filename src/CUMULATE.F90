!-----------------------------------------------------------------------
! Cumulate fluxes
!-----------------------------------------------------------------------
subroutine CUMULATE(alb,G,Gsoil,H,Hsrf,LE,LEsrf,Melt,Rnet,Roff,Rsrf)

use CONSTANTS, only: &
  Lf,                &! Latent heat of fusion (J/kg)
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

implicit none

real, intent(in) :: &
  alb(Nx,Ny),        &! Albedo
  G(Nx,Ny),          &! Heat flux into surface (W/m^2)
  Gsoil(Nx,Ny),      &! Heat flux into soil (W/m^2)
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
diags(:,:,2) = diags(:,:,2) + Gsoil(:,:)
diags(:,:,3) = diags(:,:,3) + H(:,:)
diags(:,:,4) = diags(:,:,4) + Hsrf(:,:)
diags(:,:,5) = diags(:,:,5) + LE(:,:)
diags(:,:,6) = diags(:,:,6) + LEsrf(:,:)
diags(:,:,7) = diags(:,:,7) + Lf*Melt(:,:)
diags(:,:,8) = diags(:,:,8) + Rnet(:,:)
diags(:,:,9) = diags(:,:,9) + Roff(:,:) * Nave
diags(:,:,10) = diags(:,:,10) + Rsrf(:,:)

end subroutine CUMULATE
