!-----------------------------------------------------------------------
! Write out state variables at end of run
!-----------------------------------------------------------------------
subroutine DUMP

use IOUNITS, only : &
  udmp                ! Dump file unit number

use STATE_VARIABLES, only : &
  albs,              &! Snow albedo
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers 
  Qcan,              &! Canopy air space humidity
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Sveg,              &! Canopy snow mass (kg/m^2)
  Tcan,              &! Canopy air space temperature (K)
  theta,             &! Volumetric moisture content of soil layers
  Tsnow,             &! Snow layer temperatures (K)
  Tsoil,             &! Soil layer temperatures (K)
  Tsrf,              &! Surface skin temperature (K)
  Tveg                ! Vegetation temperature (K)

implicit none

write(udmp,*) albs(:,:)
write(udmp,*) Ds(:,:,:)
write(udmp,*) Nsnow(:,:)
write(udmp,*) Qcan(:,:)
write(udmp,*) Sice(:,:,:)
write(udmp,*) Sliq(:,:,:)
write(udmp,*) Sveg(:,:)
write(udmp,*) Tcan(:,:)
write(udmp,*) theta(:,:,:)
write(udmp,*) Tsnow(:,:,:)
write(udmp,*) Tsoil(:,:,:)
write(udmp,*) Tsrf(:,:)
write(udmp,*) Tveg(:,:)
close(udmp)

end subroutine DUMP
