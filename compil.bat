::-------------------------------------------------------------------------------------------
:: Flexible Snow Model DOS compilation script
::
:: Richard Essery
:: School of GeoSciences
:: University of Edinburgh
::-------------------------------------------------------------------------------------------

cd src
(
echo /* Process options */
echo #define ALBEDO 1   /* snow albedo: 0 - diagnostic, 1 - prognostic             */
echo #define CANMOD 1   /* forest canopy: 0 - zero layer, 1 - one layer            */
echo #define CONDCT 1   /* snow thermal conductivity: 0 - constant, 1 - Yen 1981   */
echo #define DENSTY 1   /* snow density: 0 - constant, 1 - Verseghy 1991           */
echo #define EXCHNG 1   /* turbulent exchange: 0 - constant, 1 - Louis 1979        */
echo #define HYDROL 1   /* snow hydraulics: 0 - free draining, 1 - bucket          */
echo /* Driving data options */
echo #define DRIV1D 0   /* 1D driving data format: 0 - FSM, 1 - ESM-SnowMIP           */
echo #define SWPART 0   /* SW radiation: 0 - total, 1 - direct and diffuse calculated */
) > OPTS.h

set mods= MODULES.f90
set routines= CANOPY.F90 CUMULATE.F90 DRIVE.F90 DUMP.F90 EBALFOR.F90 ^
EBALSRF.F90 FSM2.F90 LUDCMP.F90 OUTPUT.F90 PHYSICS.F90 QSAT.F90 ^
READMAPS.F90 SETUP.F90 SNOW.F90 SOIL.F90 SFEXCH.F90 SWRAD.F90 ^
THERMAL.F90 TRIDIAG.F90
gfortran %mods% %routines% -cpp -o FSM2
del *.mod
move FSM2.exe ../FSM2.exe
cd ..

