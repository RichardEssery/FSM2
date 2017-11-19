!-----------------------------------------------------------------------
! Solve surface energy balance
!-----------------------------------------------------------------------
subroutine EBALOPN(Dz1,gevap,KH,ksurf,SWsurf,Ts1,   &
                   Esurf,Gsurf,Hatmo,Latmo,Melt,Rnet)


use CONSTANTS, only: &
  cp,                &! Specific heat capacity of air (J/K/kg)
  Lc,                &! Latent heat of condensation (J/kg)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  Rgas,              &! Gas constant for dry air (J/K/kg)
  Rwat,              &! Gas constant for water vapour (J/K/kg)
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm                  ! Melting point (K)

use DRIVING, only: &
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Ta                  ! Air temperature (K)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMMAPS, only: &
  VAI                 ! Vegetation area index

use STATE_VARIABLES, only: &
  Sice,              &! Ice content of snow layers (kg/m^2)
  Tsurf               ! Surface skin temperature (K)

implicit none

real, intent(in) :: &
  Dz1(Nx,Ny),        &! Surface layer thickness (m)
  gevap(Nx,Ny),      &! Surface moisture conductance (m/s)
  KH(Nx,Ny),         &! Eddy diffusivity (m/s)
  ksurf(Nx,Ny),      &! Surface layer thermal conductivity (W/m/K)
  SWsurf(Nx,Ny),     &! Net SW radiation absorbed by the surface (W/m^2)
  Ts1(Nx,Ny)          ! Surface layer temperature (K)

real, intent(out) :: &
  Esurf(Nx,Ny),      &! Moisture flux from the surface (kg/m^2/s)
  Gsurf(Nx,Ny),      &! Heat flux into surface (W/m^2)
  Hatmo(Nx,Ny),      &! Sensible heat flux to the atmosphere (W/m^2)
  Latmo(Nx,Ny),      &! Latent heat flux to the atmosphere (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny)         ! Net radiation (W/m^2)

integer :: & 
  i,j                 ! Point counters

real :: &
  D,                 &! dQsat/dT (1/K)
  dE,                &! Change in surface moisture flux (kg/m^2/s)
  dG,                &! Change in surface heat flux (W/m^2)
  dH,                &! Change in sensible heat flux (W/m^2)
  dTs,               &! Change in surface skin temperatures (K)
  Lh,                &! Latent heat (J/kg)
  psi,               &! Moisture availability factor
  Qs,                &! Saturation humidity at surface layer temperature
  rho                 ! Air density (kg/m^3)

do j = 1, Ny
do i = 1, Nx

  if (VAI(i,j) == 0) then

    ! Saturation humidity, moisture availability and air density
    call QSAT(Ps(i,j),Tsurf(i,j),Qs)
    Lh = Ls
    if (Tsurf(i,j) > Tm) Lh = Lc
    D = Lh*Qs/(Rwat*Tsurf(i,j)**2)
    psi = gevap(i,j) / (gevap(i,j) + KH(i,j))
    if (Qs < Qa(i,j) .or. Sice(1,i,j) > 0) psi = 1
    rho = Ps(i,j) / (Rgas*Ta(i,j))

    ! Explicit fluxes
    Esurf(i,j) = psi*rho*KH(i,j)*(Qs - Qa(i,j))
    Gsurf(i,j) = 2*ksurf(i,j)*(Tsurf(i,j) - Ts1(i,j))/Dz1(i,j)
    Hatmo(i,j) = cp*rho*KH(i,j)*(Tsurf(i,j) - Ta(i,j))
    Latmo(i,j) = Lh*Esurf(i,j)
    Melt(i,j) = 0
    Rnet(i,j) = SWsurf(i,j) + LW(i,j) - sb*Tsurf(i,j)**4

    ! Surface energy balance increments without melt
    dTs = (Rnet(i,j) - Hatmo(i,j) - Latmo(i,j) - Gsurf(i,j)) /  &
          ((cp + Lh*psi*D)*rho*KH(i,j) + 2*ksurf(i,j)/Dz1(i,j) + 4*sb*Tsurf(i,j)**3)
    dE = psi*rho*KH(i,j)*D*dTs
    dG = 2*ksurf(i,j)*dTs/Dz1(i,j) 
    dH = cp*rho*KH(i,j)*dTs

    ! Surface melting
    if (Tsurf(i,j)  + dTs > Tm .and. Sice(1,i,j) > 0) then
      Melt(i,j) = sum(Sice(:,i,j))/dt
      dTs = (Rnet(i,j) - Hatmo(i,j) - Latmo(i,j) - Gsurf(i,j) - Lf*Melt(i,j)) /  &
          ((cp + Lh*psi*D)*rho*KH(i,j) + 2*ksurf(i,j)/Dz1(i,j) + 4*sb*Tsurf(i,j)**3)
      dE = rho*KH(i,j)*D*dTs
      dG = 2*ksurf(i,j)*dTs/Dz1(i,j)
      dH = cp*rho*KH(i,j)*dTs
      if (Tsurf(i,j) + dTs < Tm) then
          call QSAT(Ps(i,j),Tm,Qs)
          Esurf(i,j) = rho*KH(i,j)*(Qs - Qa(i,j))  
          Gsurf(i,j) = 2*ksurf(i,j)*(Tm - Ts1(i,j))/Dz1(i,j)
          Hatmo(i,j) = cp*rho*KH(i,j)*(Tm - Ta(i,j))
          Latmo(i,j) = Ls*Esurf(i,j)
          Rnet(i,j) = SWsurf(i,j) + LW(i,j) - sb*Tm**4
          Melt(i,j) = (Rnet(i,j) - Hatmo(i,j) - Latmo(i,j) - Gsurf(i,j)) / Lf
          Melt(i,j) = max(Melt(i,j), 0.)
          dE = 0
          dG = 0
          dH = 0
          dTs = Tm - Tsurf(i,j)
      end if
    end if

    ! Update surface temperature and fluxes
    Tsurf(i,j) = Tsurf(i,j) + dTs
    Esurf(i,j) = Esurf(i,j) + dE
    Gsurf(i,j) = Gsurf(i,j) + dG
    Hatmo(i,j) = Hatmo(i,j) + dH
    Latmo(i,j) = Lh*Esurf(i,j)
    Rnet(i,j) = Hatmo(i,j) + Latmo(i,j) + Gsurf(i,j) + Lf*Melt(i,j)

  end if

end do
end do

end subroutine EBALOPN
