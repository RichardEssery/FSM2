!-----------------------------------------------------------------------
! Read data from a map file if named
!-----------------------------------------------------------------------
subroutine READMAPS(map_file,map_var)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use IOUNITS, only: &
  umap                ! Map file unit number

implicit none

character(len=*) :: &
  map_file            ! Map file name

real, intent(out) :: &
  map_var(Nx,Ny)      ! Mapped variable

if (map_file /= 'none') then
  open(umap,file=map_file)
  read(umap,*) map_var
  close(umap)
end if

end subroutine READMAPS

