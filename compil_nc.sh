########################################################################
# Flexible Snow Model compilation script for netCDF output             #
#                                                                      #
# Richard Essery                                                       #
# School of GeoSciences                                                #
# University of Edinburgh                                              #
########################################################################
FC=gfortran             # compiler
NDIR=/usr/include/      # NetCDF include path

cd src

cat > OPTS.h << EOF
/* Process options                                  : Possible values */
#define ALBEDO 2   /* snow albedo                   : 1, 2            */
#define CANINT 1   /* canopy interception of snow   : 1, 2            */
#define CANMOD 1   /* forest canopy layers          : 1, 2            */
#define CANRAD 1   /* canopy radiative properties   : 1, 2            */
#define CANUNL 1   /* unloading of canopy           : 1, 2            */
#define CONDCT 1   /* snow thermal conductivity     : 0, 1            */
#define DENSTY 1   /* snow density                  : 0, 1, 2         */
#define EXCHNG 1   /* turbulent exchange            : 0, 1            */
#define HYDROL 1   /* snow hydraulics               : 0, 1, 2         */
#define SGRAIN 1   /* snow grain growth             : 1, 2            */
#define SNFRAC 1   /* snow cover fraction           : 1, 2, 3         */
/* Driving data options                             : Possible values */
#define DRIV1D 1   /* 1D driving data format        : 1, 2            */
#define SWPART 0   /* SW radiation partition        : 0, 1            */
#define ZOFFST 0   /* measurement height offset     : 0, 1            */
/* Output options                                   : Possible values */
#define PROFNC 1   /* NetCDF output                 : 0, 1            */
EOF

$FC -cpp -O3 -o FSM2 -I$NDIR                                           \
FSM2_MODULES.F90 FSM2_PARAMS.F90 FSM2.F90 FSM2_DRIVE.F90               \
FSM2_PREPNC.F90 FSM2_VEG.F90 FSM2_WRITENC.F90 FSM2_TIMESTEP.F90        \
CANOPY.F90 INTERCEPT.F90 LUDCMP.F90 PSIMH.F90 QSAT.F90 SNOW.F90        \
SOIL.F90 SOLARPOS.F90 SRFEBAL.F90 SWRAD.F90 THERMAL.F90 TRIDIAG.F90    \
TWOSTREAM.F90 -lnetcdff -lnetcdf 
mv FSM2 ../FSM2
rm *.mod
cd ..

