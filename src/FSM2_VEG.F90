!-----------------------------------------------------------------------
! Read vegetation data from a file
!-----------------------------------------------------------------------
subroutine FSM2_VEG(map_file,Npnts,var)

use IOUNITS, only: &
  umap                ! Map input file unit number

implicit none

character(len=*) :: &
  map_file            ! Map file name

integer, intent(in) :: &
  Npnts               ! Number of points

real, intent(out) :: &
  var(Npnts)          ! Mapped variable

open(umap,file=map_file)
read(umap,*) var(:)
close(umap)

end subroutine FSM2_VEG
