########################################################################
# Compile FSM2 for UM driving data
#
# Richard Essery
# School of GeoSciences
# University of Edinburgh
########################################################################

NDIR=/usr/include

cd src

cat > OPTS.h << EOF
/* Process options                                  : Possible values */
#define ALBEDO 2   /* snow albedo                   : 1, 2            */
#define CANMOD 1   /* forest canopy layers          : 1, 2            */
#define CANRAD 1   /* canopy radiative properties   : 1, 2            */
#define CONDCT 1   /* snow thermal conductivity     : 0, 1            */
#define DENSTY 1   /* snow density                  : 0, 1            */
#define EXCHNG 1   /* turbulent exchange            : 0, 1            */
#define HYDROL 1   /* snow hydraulics               : 0, 1            */
#define SGRAIN 2   /* snow grain growth             : 1, 2            */
#define SNFRAC 2   /* snow cover fraction           : 1, 2            */
/* Driving data options                             : Possible values */
#define SETPAR 0   /* parameter inputs              : 0, 1            */
#define SWPART 0   /* SW radiation partition        : 0, 1            */
#define ZOFFST 1   /* measurement height offset     : 0, 1            */
/* Output data options                              : Possible values */
#define SMRTIO 1   /* snow layer outputs for SMRT   : 0, 1            */
EOF

gfortran -cpp -O3 -o ../FSM2 -I$NDIR                                   \
FSM2_MODULES.F90 FSM2_UM.F90 FSM2_PARAMS.F90 FSM2_TIMESTEP.F90         \
READNC.F90 CANOPY.F90 INTERCEPT.F90 LUDCMP.F90 PSIMH.F90 QSAT.F90      \
SNOW.F90 SOIL.F90 SRFEBAL.F90 SWRAD.F90 THERMAL.F90 TRIDIAG.F90        \
TWOSTREAM.F90 -lnetcdff -lnetcdf
rm *.mod
cd ..


