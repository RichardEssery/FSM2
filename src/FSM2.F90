!-----------------------------------------------------------------------
! Flexible Snow Model (FSM version 2.0)
!
! Richard Essery
! School of GeoSciences
! University of Edinburgh
!-----------------------------------------------------------------------
program FSM2

use DIAGNOSTICS, only: &
  Nave,              &! Number of timesteps in average outputs
  Nsmp                ! Timestep for sample outputs

implicit none

logical :: EoR        ! End-of-run flag

integer :: i          ! Timestep counter

call SETUP

! Loop over timesteps
EoR = .false.
do
  do i = 1, Nave
    call DRIVE(EoR)
    if (EoR) goto 1
    call PHYSICS
    if (i == Nsmp) call OUTPUT('smp')
  end do
  call OUTPUT('ave')
end do

1 continue

call DUMP

end program FSM2
