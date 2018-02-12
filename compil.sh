#############################################################################
# Flexible Snow Model compilation script
#
# Richard Essery
# School of GeoSciences
# University of Edinburgh
#############################################################################
FC=gfortran
cd src

cat > OPTS.h << EOF
/* Process options */
#define ALBEDO 1   /* snow albedo: 0 - diagnostic, 1 - prognostic                          */
#define CANMOD 0   /* forest canopy: 0 - zero layer, 1 - one layer                         */
#define CONDCT 1   /* snow thermal conductivity: 0 - constant, 1 - Yen (1981)              */
#define DENSTY 1   /* snow density: 0 - constant, 1 - Verseghy (1991), 2 - Anderson (1976) */
#define EXCHNG 1   /* turbulent exchange: 0 - constant, 1 - Louis (1979)                   */
#define HYDROL 1   /* snow hydraulics: 0 - free draining, 1 - bucket                       */

/* Driving data options */
#define DRIV1D 0   /* 1D driving data format: 0 - FSM, 1 - ESM-SnowMIP           */
#define SWPART 0   /* SW radiation: 0 - total, 1 - direct and diffuse calculated */
EOF

$FC -cpp -o FSM2 -O3 \
MODULES.F90 CANOPY.F90 CUMULATE.F90 DRIVE.F90 DUMP.F90 EBALFOR.F90  \
EBALOPN.F90 FSM2.F90 LUDCMP.F90 OUTPUT.F90 PHYSICS.F90 QSAT.F90     \
READMAPS.F90 SETUP.F90 SNOW.F90 SOIL.F90 SFEXCH.F90 SWRAD.F90       \
THERMAL.F90 TRIDIAG.F90
mv FSM2 ../FSM2
rm *.mod
cd ..

