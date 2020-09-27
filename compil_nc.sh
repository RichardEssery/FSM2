########################################################################
# Flexible Snow Model compilation script for netCDF output             #
#                                                                      #
# Richard Essery                                                       #
# School of GeoSciences                                                #
# University of Edinburgh                                              #
########################################################################
FC=gfortran
cd src

cat > OPTS.h << EOF
/* Process options                                  : Possible values */
#define ALBEDO 2   /* snow albedo                   : 1, 2            */
#define CANMOD 1   /* forest canopy layers          : 1, 2            */
#define CANRAD 1   /* canopy radiative properties   : 1, 2            */
#define CONDCT 1   /* snow thermal conductivity     : 0, 1            */
#define DENSTY 1   /* snow density                  : 0, 1, 2         */
#define EXCHNG 1   /* turbulent exchange            : 0, 1            */
#define HYDROL 1   /* snow hydraulics               : 0, 1, 2         */
#define SNFRAC 1   /* snow cover fraction           : 1, 2            */
/* Driving data options                             : Possible values */
#define DRIV1D 1   /* 1D driving data format        : 1, 2            */
#define SETPAR 1   /* parameter inputs              : 0, 1            */
#define SWPART 0   /* SW radiation partition        : 0, 1            */
#define ZOFFST 0   /* measurement height offset     : 0, 1            */
/* Output options                                   : Possible values */
#define PROFNC 1   /* NetCDF output                 : 0, 1            */
EOF

files=(FSM2_MODULES.F90 FSM2.F90 FSM2_DRIVE.F90 FSM2_MAP.F90          \
       FSM2_PARAMS.F90 FSM2_PREPNC.F90 FSM2_TIMESTEP.F90              \
       FSM2_WRITENC.F90                                               \
       CANOPY.F90 INTERCEPT.F90 LUDCMP.F90 PSIMH.F90 QSAT.F90         \
       SNOW.F90 SOIL.F90 SOLARPOS.F90 SRFEBAL.F90 SWRAD.F90           \
       THERMAL.F90 TRIDIAG.F90)
> FSM2_temp.f90
for file in ${files[*]}
do
  $FC -cpp -E -P -o out $file
  cat --squeeze-blank out >> FSM2_temp.f90
done

$FC -O3 -c -I/usr/lib64/gfortran/modules FSM2_temp.f90
$FC -o FSM2 FSM2_temp.o -L/usr/lib64 -lnetcdff
mv FSM2 ../FSM2

rm out
rm *.mod
rm *.o
rm FSM2_temp.f90
cd ..

