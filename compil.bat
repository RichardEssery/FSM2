::-------------------------------------------------------------------------
:: Flexible Snow Model DOS compilation script
::
:: Richard Essery
:: School of GeoSciences
:: University of Edinburgh
::-------------------------------------------------------------------------

cd src
(
echo /* Process options                                : Possible values */
echo #define ALBEDO 2   /* snow albedo                 : 1, 2            */
echo #define CANMOD 1   /* forest canopy layers        : 1, 2            */
echo #define CANRAD 1   /* canopy radiative properties : 1, 2            */
echo #define CONDCT 1   /* snow thermal conductivity   : 0, 1            */
echo #define DENSTY 1   /* snow density                : 0, 1, 2         */
echo #define EXCHNG 1   /* turbulent exchange          : 0, 1            */
echo #define HYDROL 1   /* snow hydraulics             : 0, 1, 2         */
echo #define SNFRAC 1   /* snow cover fraction         : 1, 2            */
echo /* Driving data options                           : Possible values */
echo #define DRIV1D 0   /* 1D driving data format      : 1, 2            */
echo #define SETPAR 1   /* parameter inputs            : 0, 1            */
echo #define SWPART 0   /* SW radiation partition      : 0, 1            */
echo #define ZOFFST 0   /* Measurement height offset   : 0, 1            */
echo /* Output options                                 : Possible values */
echo #define PROFNC 0   /* netCDF output               : 0, 1            */
) > OPTS.h

set mods= FSM2_MODULES.f90

set routines= FSM2.F90 FSM2_DRIVE.F90 FSM2_MAP.F90 FSM2_OUTPUT.F90 ^
FSM2_PARAMS.F90 FSM2_TIMESTEP.F90 CANOPY.F90 INTERCEPT.F90 LUDCMP.F90 ^
PSIMH.F90 QSAT.F90 SNOW.F90 SOIL.F90 SOLARPOS.F90 SRFEBAL.F90 SWRAD.F90 ^
THERMAL.F90 TRIDIAG.F90

gfortran %mods% %routines% -cpp -o FSM2
del *.mod
move FSM2.exe ../FSM2.exe
cd ..

