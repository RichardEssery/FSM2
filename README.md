# FSM2

The Factorial Snow Model (FSM) is a multi-physics energy balance model of accumulation and melt of snow on the ground. Version 2 adds forest canopy model options and the possibility of running simulations for more than one point at the same time. In contrast with FSM1, which selects options when it is run, FSM2 options are selected when it is compiled for greater efficiency. A full description is given by [Essery (2015)](#Essery). FSM2 is built and run in the same way as FSM1.

## Building the model

FSM2 is coded in Fortran. A linux executable `FSM2` or a Windows executable `FSM2.exe` is produced by running the script `compil.sh` or the batch file `compil.bat`. Both use the [gfortran](https://gcc.gnu.org/wiki/GFortran) compiler but could be edited to use other compilers. 

Model options are selected by defining variables in `src/OPTS.h` before compilation

| Variable | Description                  | Options                        |
|----------|------------------------------|--------------------------------|
| ALBEDO   | Snow albedo options          | <ul><li>diagnosed</li><li>prognostic</li></ul> |
| CANMOD   | Canopy model options         | 0 zero-layer <br> 1 - one-layer)                |
| CONDCT   | Thermal conductivity options |(0 - fixed, 1 - density function) |
| DENSTY   | Snow density options         |(0 - fixed, 1 - prognostic)               |
| EXCHNG   | Surface exchange options     |(0 - fixed, 1 - stability adjusted)   |
| HYDROL   | Snow hydrology options       | (0 - free draining, 1 - bucket)      | 



## Running the model

FSM2 requires meteorological driving data and namelists to set options and parameters. The model is run with the command

    ./FSM2 < nlst.txt

or

    FSM2.exe < nlst.txt

where `nlst.txt` is a text file containing seven namelists described below; `nlst_CdP_0506.txt` gives an example to run FSM for the winter of 2005-2006 at Col de Porte ([Morin et al. 2011](#Morin)). All of the namelists have to be present in the same order as in the example, but any or all of the namelist variables listed in the tables below can be omitted; defaults are then used.

### Driving data

Meteorological driving data are read from a text file named in namelist `&drive`. A driving data file has 12 columns containing the variables listed in the table below. Each row of the file corresponds with a timestep. Driving data for the Col de Porte example are given in file `data/met_CdP_0506.txt`.

| Variable | Units  | Description       |
|----------|--------|-------------------|
| year     | years  | Year              |
| month    | months | Month of the year |
| day      | days   | Day of the month  |
| hour     | hours  | Hour of the day   |
| SW       | W m<sup>-2</sup> | Incoming shortwave radiation  |
| LW       | W m<sup>-2</sup> | Incoming longwave radiation   |
| Sf       | kg m<sup>-2</sup> s<sup>-1</sup> | Snowfall rate |
| Rf       | kg m<sup>-2</sup> s<sup>-1</sup> | Rainfall rate |
| Ta       | K      | Air temperature      |
| RH       | RH     | Relative humidity    |
| Ua       | m s<sup>-1</sup> | Wind speed |
| Ps       | Pa     | Surface air pressure |

### Grid dimensions namelist 

`&gridpnts`

| Variable | Default | Description |
|----------|---------|-------------|
| Nsmax    | 3       | Maximum number of snow layers                    |
| Nsoil    | 4       | Number of soil layers                            |
| Nx       | 1       | Number of grid points in x direction or sequence |
| Ny       | 1       | Number of grid points in y direction             |

FSM2 can be run at a point, at a sequence of points or on a rectangular grid.


### Model levels namelist 

`&gridlevs`

| Variable | Default          | Units | Description |
|----------|------------------|-------|-------------|
| Dzsnow   | 0.1 0.2 0.4      | m     | Minimum snow layer thicknesses  |
| Dzsoil   | 0.1 0.2 0.4 0.8  | m     | Soil layer thicknesses          |

If layer thicknesses are specified in `&gridlevs`, they must match the numbers of layers specfied in `&gridpnts`; this is not checked automatically.

### Driving data namelist 

`&drive`

| Variable | Default | Units | Description |
|----------|---------|-------|-------------|
| met_file | 'met.txt' | string  | Driving file name              |
| dt       | 3600      | s       | Time step                      |
| zT       | 2         | m       | Temperature measurement height |
| zU       | 10        | m       | Wind speed measurement height  |

Measurement heights have to be above the canopy height.

### Parameter namelist 

`&params`

| Variable | Default | Units | Description |
|----------|---------|-------|-------------|
| asmx | 0.8  | -    | Maximum albedo for fresh snow                                   |
| asmn | 0.5  | -    | Minimum albedo for melting snow                                 |
| avgs | 0.2  | -    | Snow-covered vegetation albedo                                  |
| bstb | 5    | -    | Atmospheric stability adjustment parameter (if n<sub>e</sub>=1) |
| bthr | 2    | -    | Thermal conductivity exponent (if n<sub>c</sub>=1)              |
| gsat | 0.01 | m s<sup>-1</sup>  | Surface conductance for saturated soil             |
| canc | 4.4  | kg m<sup>-2</sup> | Canopy snow capacity per unit vegetation area      |
| cmlt | 240  | days | Melting canopy snow unloading time scale                        | 
| cunl | 2.4  | days | Cold canopy snow unloading time scale                           | 
| hfsn | 0.1  | m    | Snow cover fraction depth scale                                 |
| kext | 0.5  | -    | Canopy radiation extinction coefficient                         |
| kfix | 0.24 | W m<sup>-1</sup> K<sup>-1</sup> | Fixed thermal conductivity (if n<sub>c</sub>=0) |
| rho0 | 300  | kg m<sup>-3</sup> | Fixed snow density (if n<sub>d</sub>=0)               |
| rhof | 100  | kg m<sup>-3</sup> | Fresh snow density (if n<sub>d</sub>=1)               |
| rcld | 300  | kg m<sup>-3</sup> | Maximum density for cold snow (if n<sub>d</sub>=1)    |
| rmlt | 500  | kg m<sup>-3</sup> | Maximum density for melting snow (if n<sub>d</sub>=1) |
| Salb | 10   | kg m<sup>-2</sup> | Snowfall to refresh albedo (if n<sub>a</sub>=1)       |
| Talb | -2   | &deg;C| Albedo decay temperature threshold (if n<sub>a</sub>=0)           |
| tcld | 1000 | h    | Cold snow albedo decay time scale (if n<sub>a</sub>=1)             |
| tmlt | 100  | h    | Melting snow albedo decay time scale (if n<sub>a</sub>=1)          |
| trho | 200  | h    | Compaction time scale (if n<sub>d</sub>=1)                         |
| Wirr | 0.03 | -    | Irreducible liquid water content (if n<sub>w</sub>=1)              |
| z0sn | 0.01 | m    | Snow roughness length                                              |


### Site characteristics 

`&maps`

| Variable | Default | Units | Description |
|----------|---------|-------|-------------|
| alb0 | 0.2  | -    | Snow-free ground albedo                                         |
| fcly | 0.3  | -    | Soil clay fraction                                              |
| fsnd | 0.6  | -    | Soil sand fraction                                              |
| fsky | exp(-kext*VAI) | -    | Sky view fraction                                     |
| hcan | 0    | m    | Canopy height                                                   |
| scap | canc*VAI    | kg m<sup>-2</sup> | Canopy snow capacity                        |
| VAI  | 0    | -    | Vegetation area index                                           |
| z0sf | 0.1  | m    | Snow-free roughness length                                      |

Site characteristics can either be left as default values, set to a sequence of values in the namelist or read from a named map file. e.g. for a simulation with 10 points, the snow-free ground albedo can be reset to a constant value of 0.1 by including

    alb0 = 10*0.1

in the namelist, set to sequence by including

    alb0 = 0.2 0.2 0.1 0.1 0.1 0.1 0.1 0.1 0.2 0.2

or read from a file albedo.txt containing 10 values by including

    alb0_file = 'albedo.txt'

Note that sky view can be set independently of vegetation cover to allow for grid cells shaded by topography or vegetation in neighbouring cells.



### Initial values namelist 

`&initial`

| Variable   | Default | Units  | Description |
|------------|---------|--------|-------------|
| start_file | 'none'  | string | Start file  |
| fsat       | 4 * 0.5 | -      | Initial moisture content of soil layers as fractions of saturation |
| Tsoil      | 4 * 285 | K      | Initial temperatures of soil layers                                |

Soil temperature and moisture content are taken from the namelist and FSM is initialized in a snow-free state by default. If a start file is named, it should be a text file containing initial values for each of the state variables in order:

| Variable    |  Units             | Description |
|-------------|--------------------|-------------|
|  albs       |  -                 | Snow albedo |
|  Ds(1:3)    |  m                 | Snow layer thicknesses                     |
|  Nsnow      |  -                 | Number of snow layers                      | 
|  Sice(1:3)  |  kg m<sup>-2</sup> | Ice content of snow layers                 |
|  Sliq(1:3)  |  kg m<sup>-2</sup> | Liquid content of snow layers              |
|  theta(1:4) |  -                 | Volumetric moisture content of soil layers |
|  Tsnow(1:3) |  K                 | Snow layer temperatures                    | 
|  Tsoil(1:4) |  K                 | Soil layer temperatures                    |
|  Tsurf      |  K                 | Surface skin temperature                   |

Snow and soil layers are numbered from the top downwards. File `data/dump.txt` is an example of a dump file produced at the end of a run which could be used to restart the model.

### Output namelist 

`&outputs`

| Variable  | Default    | Description |
|-----------|------------|-------------|
| Nave      | 24         | Number of timesteps in averaged outputs |
| out_file  | 'out.txt'  | Output file name |
| dump_file | 'dump.txt' | Dump file name   |

At present a simple fixed output format is used. The output text file has 10 columns:

| Variable | Units  | Description       |
|----------|--------|-------------------|
| year     | years  | Year              |
| month    | months | Month of the year |
| day      | days   | Day of the month  |
| hour     | hours  | Hour of the day   |
| alb      | -      | Effective albedo  |
| Rof      | kg m<sup>-2</sup> | Cumulated runoff from snow    |
| snd      | m      | Average snow depth                       |
| SWE      | kg m<sup>-2</sup> | Average snow water equivalent |
| Tsf      | &deg;C | Average surface temperature              |
| Tsl      | &deg;C | Average soil temperature at 20 cm depth  |

Example file `data/out_CdP_0506.txt` contains output from a run of FSM configuration 31; `data/obs_CdP_0506.txt` contains daily observations of the same variables, with -99 indicating missing data.

At the end of a run, the state variables are written to a dump file with the same format as the start file.
 

## References

<a name="Essery"></a> Essery (2015). A Factorial Snowpack Model (FSM 1.0). *Geoscientific Model Development*, **8**, 3867-3876, [doi:10.5194/gmd-8-3867-2015](http://www.geosci-model-dev.net/8/3867/2015/)

<a name="Morin"></a> Morin et al. (2012). A 18-yr long (1993-2011) snow and meteorological dataset from a mid-altitude mountain site (Col de Porte, France, 1325 m alt.) for driving and evaluating snowpack models. *Earth System Science Data*, **4**(1), 13-21, [doi:10.5194/essd-4-13-2012](http://www.earth-syst-sci-data.net/4/13/2012/essd-4-13-2012.html)

