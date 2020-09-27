!-----------------------------------------------------------------------
! Read data from a map file
!-----------------------------------------------------------------------
subroutine FSM2_MAP(map_file,Ncols,Nrows,var)

use IOUNITS, only: &
  umap                ! Map input file unit number

implicit none

character(len=*) :: &
  map_file            ! Map file name

integer, intent(in) :: &
  Ncols,             &! Number of columns in grid
  Nrows               ! Number of rows in grid

real, intent(out) :: &
  var(Nrows,Ncols)    ! Mapped variable

integer :: &
  i                   ! Grid row counter

open(umap,file=map_file)
do i = 1, Nrows
  read(umap,*) var(i,:)
end do
close(umap)

end subroutine FSM2_MAP
