::-------------------------------------------------------------------------------------------
:: Factorial Snow Model DOS compilation script
::
:: Richard Essery
:: School of GeoSciences
:: University of Edinburgh
::-------------------------------------------------------------------------------------------

cd src
set mods= MODULES.f90
set routines= CANOPY.F90 CUMULATE.f90 DRIVE.f90 DUMP.f90 EBALFOR.F90 EBALOPN.F90 FSM2.F90 LUDCMP.F90 OUTPUT.F90 PHYSICS.F90 QSAT.F90 READMAPS.F90 SETUP.F90 SFEXCH.F90 SNOW.F90 SOIL.F90 SWRAD.F90 THERMAL.F90 TRIDIAG.F90
gfortran %mods% %routines% -cpp -o FSM.exe
del *.mod
move FSM.exe ../FSM.exe
cd ..

