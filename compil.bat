::-------------------------------------------------------------------------------------------
:: Flexible Snow Model DOS compilation script
::
:: Richard Essery
:: School of GeoSciences
:: University of Edinburgh
::-------------------------------------------------------------------------------------------

cd src
(
echo /* Process options                                : Possible values */
echo #define ALBEDO 1   /* snow albedo                 : 0, 1            */
echo #define CANMOD 1   /* forest canopy               : 0, 1            */
echo #define CONDCT 1   /* snow thermal conductivity   : 0, 1            */
echo #define DENSTY 1   /* snow density                : 0, 1, 2         */
echo #define EXCHNG 1   /* turbulent exchange          : 0, 1            */
echo #define HYDROL 1   /* snow hydraulics             : 0, 1            */
echo #define SNFRAC 1   /* snow cover fraction         : 0, 1            */
echo /* Driving data options                           : Possible values */
echo #define DRIV1D 0   /* 1D driving data format      : 0, 1, 2         */
echo #define DOWNSC 0   /* 1D driving data downscaling : 0, 1            */
echo #define DEMHDR 0   /* DEM header                  : 0, 1            */
echo #define SWPART 0   /* SW radiation partition      : 0, 1, 2         */
echo #define ZOFFST 0   /* Measurement height offset   : 0, 1            */
) > OPTS.h

set mods= MODULES.f90
set routines= CANOPY.F90 CUMULATE.F90 DRIVE.F90 DUMP.F90 EBALFOR.F90 ^
EBALSRF.F90 FSM2.F90 LUDCMP.F90 OUTPUT.F90 PHYSICS.F90 QSAT.F90 ^
RADIATION.F90 READMAPS.F90 READ_DEM.F90 SETUP.F90 SNOW.F90 SOIL.F90 ^
SFEXCH.F90 SOLARPOS.F90 THERMAL.F90 TRIDIAG.F90
gfortran %mods% %routines% -cpp -o FSM2
del *.mod
move FSM2.exe ../FSM2.exe
cd ..

