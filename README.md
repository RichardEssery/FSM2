# FSM2 quickstart guide

The Flexible Snow Model (FSM2) is a multi-physics energy balance model of snow accumulation and melt, extending the Factorial Snow Model [(Essery, 2015)](#Essery2015) with additional physics, driving and output options. FSM2 adds forest canopy model options and the possibility of running simulations for more than one point at the same time. For greater efficiency than FSM, which selects physics options when it is run, FSM2 options are selected when the model is compiled. Otherwise, FSM2 is built and run in the same way as FSM; for details, see the user guide in docs.

## Building the model

FSM2 is coded in Fortran and consists of subroutines and modules contained in the src directory. A linux executable FSM2 is produced by running script compil.sh, which uses the [gfortran](https://gcc.gnu.org/wiki/GFortran) compiler. Physics and driving data configurations are selected in the compilation script by defining options that are copied to a preprocessor file before compilation.

## Running the model

FSM2 requires meteorological driving data and namelists to set options and parameters. An example can be run with the command

    ./FSM2 < nlst_Sod_1314

which run simulations for the winter of 2013-2014 at Sodankylä, Finland [(Essery et al, 2016)](#Essery2016). Two points are simulated: one with forest cover and one without.

## References

<a name="Essery2015"></a> Essery (2015). A Factorial Snowpack Model (FSM 1.0). *Geoscientific Model Development*, **8**, 3867-3876, [doi:10.5194/gmd-8-3867-2015](http://www.geosci-model-dev.net/8/3867/2015/)

<a name="Essery2016"></a> Essery et al. (2016). A 7-year dataset for driving and evaluating snow models at an Arctic site (Sodankylä, Finland). *Geosci. Instrum. Method. Data Syst.*, **5**, 219-227, [doi:10.5194/gi-5-219-2016](https://www.geosci-instrum-method-data-syst.net/5/219/2016/)


