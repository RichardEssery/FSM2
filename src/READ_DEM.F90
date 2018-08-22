!-----------------------------------------------------------------------
! Read elevations from a DEM file with an ESRI header
!-----------------------------------------------------------------------
subroutine READ_DEM(ztop_file)

use GRID, only: &
  Nx,Ny,             &! Grid dimensions
  ztop                ! Land surface elevations (m)

use IOUNITS, only: &
  umap                ! Map file unit number

implicit none

character(len=*), intent(in) :: &
  ztop_file           ! DEM file name

character(len=12) :: &
  hdr                 ! Header text

open(umap,file=ztop_file)
read(umap,*) hdr, Nx
read(umap,*) hdr, Ny
read(umap,*)
read(umap,*)
read(umap,*)
read(umap,*)
allocate(ztop(Nx,Ny))
read(umap,*) ztop
close(umap)

end subroutine READ_DEM

