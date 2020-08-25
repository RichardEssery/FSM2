#############################################################################
# Flexible Snow Model compilation script
#
# Richard Essery
# School of GeoSciences
# University of Edinburgh
############################################################################
# Make sure line seperator is Unix-style (LF)
# see https://www.jetbrains.com/help/pycharm/configuring-line-endings-and-line-separators.html
#############################################################################
FC=gfortran
cd src

cat > OPTS.h << EOF
/* Process options                                : Possible values */
#define ALBEDO 1   /* snow albedo                 : 0, 1            */
#define CANMOD 1   /* forest canopy               : 0, 1            */
#define CONDCT 1   /* snow thermal conductivity   : 0, 1            */
#define DENSTY 1   /* snow density                : 0, 1, 2         */
#define EXCHNG 1   /* turbulent exchange          : 0, 1            */
#define HYDROL 1   /* snow hydraulics             : 0, 1            */
#define SNFRAC 1   /* snow cover fraction         : 0, 1            */

/* Driving data options                           : Possible values */
#define DRIV1D 0   /* 1D driving data format      : 0, 1, 2         */
#define DOWNSC 0   /* 1D driving data downscaling : 0, 1            */
#define DEMHDR 0   /* DEM header                  : 0, 1            */
#define SWPART 0   /* SW radiation partition      : 0, 1, 2         */
#define ZOFFST 0   /* Measurement height offset   : 0, 1            */
EOF

$FC -cpp -o FSM2 -O3  \
MODULES.F90 CANOPY.F90 CUMULATE.F90 DRIVE.F90 DUMP.F90 EBALFOR.F90  \
EBALSRF.F90 FSM2.F90 LUDCMP.F90 OUTPUT.F90 PHYSICS.F90 QSAT.F90     \
RADIATION.F90 READMAPS.F90 READ_DEM.F90 SETUP.F90 SNOW.F90 SOIL.F90 \
SFEXCH.F90 SOLARPOS.F90 THERMAL.F90 TRIDIAG.F90
mv FSM2 ../FSM2
rm *.mod
cd ..

