!-----------------------------------------------------------------------
! Write output
!-----------------------------------------------------------------------
subroutine OUTPUT(type)

#include "OPTS.h"

use DRIVING, only: &
  year,              &! Year
  month,             &! Month of year
  day,               &! Day of month
  hour                ! Hour of day

use IOUNITS, only: &
  uave,              &! Average output file unit number
  uprf,              &! Profile output file unit number
  usmp                ! Sample output file unit number

use DIAGNOSTICS, only: &
  diags,             &! Cumulated diagnostics
  Nave,              &! Number of timesteps in average outputs
  Ndiags,            &! Number of averaged diagnostics
  SWin,              &! Cumulated incoming solar radiation (J/m^2)
  SWout               ! Cumulated reflected solar radiation (J/m^2)

use GRID, only: &
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers
  rgrn,              &! Snow layer grain radius (m)
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Sveg,              &! Snow mass on vegetation (kg/m^2)
  theta,             &! Volumetric moisture content of soil layers
  Tsnow,             &! Snow layer temperatures (K)
  Tsoil,             &! Soil layer temperatures (K)
  Tsrf,              &! Surface skin temperature (K)
  Tveg                ! Vegetation temperature (K)

implicit none

character(len=3), intent(in) :: &
  type                ! Output type

integer :: &
  i,j,               &! Point counters
  k                   ! Level counter

real :: &
  alb(Nx,Ny),        &! Effective albedo
  snowdepth(Nx,Ny),  &! Snow depth (m)
  SWE(Nx,Ny)          ! Snow water equivalent (kg/m^2) 

real :: &
  zl                  ! Height of layer above ground (m)

! Output state variable samples and profiles
if (type == 'smp') then

  do j = 1, Ny
  do i = 1, Nx
    snowdepth(i,j) = sum(Ds(:,i,j))
    SWE(i,j) = sum(Sice(:,i,j)) + sum(Sliq(:,i,j))
  end do
  end do
  write(usmp,100) year,month,day,hour,snowdepth(:,:),SWE(:,:),Sveg(:,:),  &
                                      Tsoil(2,:,:),Tsrf(:,:),Tveg(:,:)

  if (Nx*Ny == 1) then
    write(uprf,'(5i4)') year,month,day,Nsnow(1,1),Nsoil
    zl = snowdepth(1,1)
    do k = 1, Nsnow(1,1)
      zl = zl - 0.5*Ds(k,1,1)
      write(uprf,'(6e14.5)') zl,Ds(k,1,1),Tsnow(k,1,1),rgrn(k,1,1),Sice(k,1,1),Sliq(k,1,1)
      zl = zl - 0.5*Ds(k,1,1)
    end do
    zl = 0
    do k = 1, Nsoil
      zl = zl - 0.5*Dzsoil(k)
      write(uprf,'(4e14.5)') zl,Dzsoil(k),Tsoil(k,1,1),theta(k,1,1)
      zl = zl - 0.5*Dzsoil(k)
    end do
  else
    stop 'Profile output only available for 1D simulations'
  end if

end if

! Output and reset averages
if (type == 'ave') then
  do j = 1, Ny
  do i = 1, Nx
    if (SWin(i,j) > 0) then
      alb(i,j) = SWout(i,j) / SWin(i,j)
    else
      alb(i,j) = -9
    end if
  end do
  end do
  diags(:,:,:) = diags(:,:,:) / Nave
  write(uave,100) year,month,day,hour,alb,diags(:,:,:)
  diags(:,:,:) = 0
  SWin(:,:) = 0
  SWout(:,:) = 0
end if

100 format(3(i4),*(f12.3))

end subroutine OUTPUT
