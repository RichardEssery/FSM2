# FSM2 quickstart guide

The Flexible Snow Model (FSM2) is a multi-physics energy balance model of snow accumulation and melt, extending the Factorial Snow Model (Essery, 2015) with additional physics, driving and output options. FSM2 adds forest canopy model options and the possibility of running simulations for more than one point at the same time. For greater efficiency than FSM, which selects physics options when it is run, FSM2 options are selected when the model is compiled. Otherwise, FSM2 is built and run in the same way as FSM; for details, see the user guide in docs.

## Building the model

FSM2 is coded in Fortran and consists of subroutines and modules contained in the src directory. A linux executable FSM2 is produced by running script compil.sh, which uses the [gfortran](https://gcc.gnu.org/wiki/GFortran) compiler by default. Physics and driving data configurations are selected in the compilation script by defining options that are copied to a preprocessor file before compilation.

## Running the model

FSM2 requires meteorological driving data and namelists to set options and parameters. An example can be run with the command

    ./FSM2 < nlst_Alptal.txt

which runs simulations for the winter of 2004-2005 at Alptal, Switzerland (Stähli and Gustafsson, 2006). Two points are simulated: one with forest cover and one without.

## References

Essery (2015). A Factorial Snowpack Model (FSM 1.0). *Geoscientific Model Development*, **8**, 3867-3876, [doi:10.5194/gmd-8-3867-2015](http://www.geosci-model-dev.net/8/3867/2015/)

Stähli and Gustafsson (2006). The role of snow interception in winter-time radiation processes of a coniferous sub-alpine forest.
*Hydrological Processes*, **23**, 2498–2512, [doi:10.1002/hyp.7180](https://onlinelibrary.wiley.com/doi/abs/10.1002/hyp.7180)


