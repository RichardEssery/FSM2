!-----------------------------------------------------------------------
! Write output
!-----------------------------------------------------------------------
subroutine OUTPUT(type)

use DRIVING, only: &
  year,              &! Year
  month,             &! Month of year
  day,               &! Day of month
  hour                ! Hour of day

use IOUNITS, only: &
  uave,              &! Average output file unit number
  usmp                ! Sample output file unit number

use DIAGNOSTICS, only: &
  diags,             &! Cumulated diagnostics
  Nave,              &! Number of timesteps in average outputs
  Ndiags,            &! Number of averaged diagnostics
  SWin,              &! Cumulated incoming solar radiation (J/m^2)
  SWout               ! Cumulated reflected solar radiation (J/m^2)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Sveg                ! Snow mass on vegetation (kg/m^2)

implicit none

character(len=3), intent(in) :: &
  type                ! Output type

integer :: &
  i,j,               &! Point counters
  N                   ! Diagnostic counter

real :: &
  alb(Nx,Ny),        &! Effective albedo
  snowdepth(Nx,Ny),  &! Snow depth (m)
  SWE(Nx,Ny)          ! Snow water equivalent (kg/m^2) 

! Output samples
if (type == 'smp') then
  do j = 1, Ny
  do i = 1, Nx
    snowdepth(i,j) = sum(Ds(:,i,j))
    SWE(i,j) = sum(Sice(:,i,j)) + sum(Sliq(:,i,j))
  end do
  end do
  if (Nx == 1 .or. Ny == 1) then
    write(usmp,100) year,month,day,hour,snowdepth(:,:),SWE(:,:),Sveg(:,:)
  else
    write(usmp,'(4i4)') year,month,day,hour
    write(usmp,*) snowdepth(:,:)
    write(usmp,*) SWE(:,:)
    write(usmp,*) Sveg(:,:)
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
  if (Nx == 1 .or. Ny == 1) then
    write(uave,100) year,month,day,hour,alb,diags(:,:,:)
  else
    write(uave,'(3i4)') year,month,day
    do n = 1, Ndiags
      write(uave,*) diags(:,:,n)
    end do
  end if
  diags(:,:,:) = 0
  SWin(:,:) = 0
  SWout(:,:) = 0
end if

100 format(3(i4),*(f12.3))

end subroutine OUTPUT
