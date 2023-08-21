########################################################################
# Compile and run FSM2 with UM driving data
#
# Richard Essery
# School of GeoSciences
# University of Edinburgh
########################################################################

module load gcc/8.1.0 mpi/mpich/3.2.1/gnu hdf5/1.8.20/gnu/8.1.0 netcdf/4.6.1/gnu/8.1.0
NDIR=/project/ukmo/rhel7/fortran/opt/gfortran/packages/gnu/8.1.0/netcdf/4.6.1/include
LNDIR=/project/ukmo/rhel7/fortran/opt/gfortran/packages/gnu/8.1.0/netcdf/4.6.1/lib
MDIR=

mkdir -p output
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
#define SGRAIN 1   /* snow grain growth             : 1, 2            */
#define SNFRAC 2   /* snow cover fraction           : 1, 2            */
/* Driving data options                             : Possible values */
#define SETPAR 0   /* parameter inputs              : 0, 1            */
#define SWPART 0   /* SW radiation partition        : 0, 1            */
#define ZOFFST 1   /* measurement height offset     : 0, 1            */
/* Output data options                              : Possible values */
#define SMRTIO 1   /* snow layer outputs for SMRT   : 0, 1            */
EOF

gfortran -cpp -O3 -o ../FSM2 -I$NDIR -L$LNDIR                          \
FSM2_MODULES.F90 FSM2_UM.F90 FSM2_PARAMS.F90 FSM2_TIMESTEP.F90         \
READNC.F90 CANOPY.F90 INTERCEPT.F90 LUDCMP.F90 PSIMH.F90 QSAT.F90      \
SNOW.F90 SOIL.F90 SRFEBAL.F90 SWRAD.F90 THERMAL.F90 TRIDIAG.F90        \
TWOSTREAM.F90 -lnetcdff -lnetcdf
rm *.mod
cd ..

months=(20171002-20171130 20171130-20180131 20180131-20180401)
vars=(LWdown Psurf Rainf RelHum Snowf SWdown Tair Wind)
rm -f start.txt
for m in ${months[*]}
do
    echo $m
    for v in ${vars[*]}
    do
      # ln -s -f $MDIR'/'$v'_'$m'.nc' $v.nc
      ln -s -f $MDIR'/NetCDF/'$v'_'$m'.nc' $v.nc
    done
    ./FSM2
    mv dump.txt start.txt
    mv FSM2out.nc 'output/FSM2out_'$y$m'.nc'
done
rm -f start.txt

