# FSM2

The Flexible Snow Model (FSM2) is a multi-physics energy balance model of snow accumulation and melt, extending the Factorial Snow Model [(Essery, 2015)](#Essery2015) with additional physics, driving and output options. FSM2 adds forest canopy model options and the possibility of running simulations for more than one point at the same time. In contrast with FSM, which selects physics options when it is run, FSM2 options are selected when compiled for greater efficiency. Otherwise, FSM2 is built and run in the same way as FSM.

## Building the model

FSM2 is coded in Fortran. A linux executable FSM2 or a Windows executable FSM2.exeis produced  by running the script compil.sh or the batch file compil.bat. Both use the [gfortran](https://gcc.gnu.org/wiki/GFortran) compiler but could be edited to use other compilers. Physics and driving data configurations are selected in the compilation script by defining options that are copied to a preprocessor file src/OPTS.h before compilation.

### Physics options
| Option  | Description                  | Possible choices                      |
|---------|------------------------------|---------------------------------------|
| ALBEDO  | Snow albedo options          | 0 - diagnostic <br> 1 - prognostic    |
| CANMOD  | Canopy model options         | 0 - zero layer <br> 1 - one layer     |
| CONDCT  | Thermal conductivity options | 0 - fixed <br> 1 - density function   |
| DENSTY  | Snow density options         | 0 - fixed <br> 1 - Verseghy (1991) <br> 2 - Anderson (1976) |
| EXCHNG  | Surface exchange options     | 0 - fixed <br> 1 - stability adjusted |
| HYDROL  | Snow hydrology options       | 0 - free draining <br> 1 - bucket     | 

### Driving data options
| Option  | Description                  | Possible choices                      |
|---------|------------------------------|---------------------------------------|
| DRIV1D  | 1D driving data format       | 0 - FSM format <br> 1 - ESM-SnowMIP format                            | 
| SWPART  | SW radiation partition       | 0 - Total SW radiation used <br> 1 - Diffuse and direct SW calculated | 


## Running the model

FSM2 requires meteorological driving data and namelists to set options and parameters. The model is run with the command

    ./FSM2 < nlst

or

    FSM2.exe < nlst

where nlst is a text file containing eight namelists described below. All of the namelists have to be present in the order given below, but any or all of the namelist variables listed in the tables can be omitted; defaults are then used.

### Grid dimensions namelist &gridpnts

| Variable | Default | Description |
|----------|---------|-------------|
| Nsmax    | 3       | Maximum number of snow layers                    |
| Nsoil    | 4       | Number of soil layers                            |
| Nx       | 1       | Number of grid points in x direction or sequence |
| Ny       | 1       | Number of grid points in y direction             |

FSM2 can be run at a point, at a sequence of points or on a rectangular grid by selecting values for Nx and Ny.

### Model levels namelist &gridlevs

| Variable | Default          | Units | Description |
|----------|------------------|-------|-------------|
| Dzsnow   | 0.1 0.2 0.4      | m     | Minimum snow layer thicknesses  |
| Dzsoil   | 0.1 0.2 0.4 0.8  | m     | Soil layer thicknesses          |

Snow and soil layers are numbered from the top downwards. If layer thicknesses are specified in &gridlevs, they must match the numbers of layers specified in &gridpnts; this is not checked automatically.

### Driving data namelist &drive and data files

| Variable | Default | Units   | Description                    |
|----------|---------|---------|--------------------------------|
| met_file | 'met'   | string  | Driving file name              |
| dt       | 3600    | s       | Time step                      |
| lat      | 0       | degrees | Latitude                       |
| noon     | 0       | hours   | Local offset from solar noon   |
| zT       | 2       | m       | Temperature measurement height |
| zU       | 10      | m       | Wind speed measurement height  |

Measurement heights have to be above the canopy height.

For simulations at a point or for a set of nearby points with common meteorology, 1D driving data are read from the text file named in namelist &drive. The file should have 12 columns containing the variables listed in the table below. Each row of the file corresponds with a timestep.

| Variable | Units  | Description                                    |
|----------|--------|------------------------------------------------|
| year     | years  | Year                                           |
| month    | months | Month of the year                              |
| day      | days   | Day of the month                               |
| hour     | hours  | Hour of the day                                |
| SW       | W m<sup>-2</sup> | Incoming shortwave radiation         |
| LW       | W m<sup>-2</sup> | Incoming longwave radiation          |
| Sf       | kg m<sup>-2</sup> s<sup>-1</sup> | Snowfall rate        |
| Rf       | kg m<sup>-2</sup> s<sup>-1</sup> | Rainfall rate        |
| Ta       | K      | Air temperature                                |
| RH       | RH     | Relative humidity                              |
| Ua       | m s<sup>-1</sup> | Wind speed                           |
| Ps       | Pa     | Surface air pressure                           |

### Parameter namelist &params

| Variable | Default | Units | Description |
|----------|---------|-------|-------------|
| asmx | 0.8  | -    | Maximum albedo for fresh snow                               |
| asmn | 0.5  | -    | Minimum albedo for melting snow                             |
| avgs | 0.2  | -    | Snow-covered vegetation albedo                              |
| bstb | 5    | -    | Atmospheric stability adjustment parameter (if EXCHNG=1)    |
| bthr | 2    | -    | Thermal conductivity exponent (if CONDCT=1)                 |
| gsat | 0.01 | m s<sup>-1</sup>  | Surface conductance for saturated soil         |
| canc | 4.4  | kg m<sup>-2</sup> | Canopy snow capacity per unit vegetation area  |
| cunc | 240  | hours| Canopy unloading time scale for cold snow                   | 
| cunm | 2.4  | hours| Canopy unloading time scale for melting snow <br> (immediate unloading if < dt) | 
| eta0 | 3.7e7| Pa s                      | Reference snow viscosity (if DENSTY=2) |
| etaa | 0.081| K<sup>-1</sup>            | Snow viscosity parameter (if DENSTY=2) |
| etab | 0.018| m<sup>3</sup>kg<sup>-1</sup> | Snow viscosity parameter (if DENSTY=2) |
| hfsn | 0.1  | m    | Snow cover fraction depth scale                             |
| kext | 0.5  | -    | Canopy radiation extinction coefficient                     |
| kfix | 0.24 | W m<sup>-1</sup> K<sup>-1</sup> | Fixed thermal conductivity (if CONDCT=0)|
| rho0 | 300  | kg m<sup>-3</sup> | Fixed snow density (if DENSTY=0)               |
| rhoc | 150  | kg m<sup>-3</sup> | Critical snow density (if DENSTY=2)            |
| rhof | 100  | kg m<sup>-3</sup> | Fresh snow density (if DENSTY=1)               |
| rcld | 300  | kg m<sup>-3</sup> | Maximum density for cold snow (if DENSTY=1)    |
| rmlt | 500  | kg m<sup>-3</sup> | Maximum density for melting snow (if DENSTY=1) |
| Salb | 10   | kg m<sup>-2</sup> | Snowfall to refresh albedo (if ALBEDO=1)       |
| snda | 2.8e-6| s<sup>-1</sup>   | Snow densification parameter (if DENSTY=2)     |
| sndb | 0.042| K<sup>-1</sup>    | Snow densification parameter (if DENSTY=2)     |
| sndc | 0.046| m<sup>3</sup> kg<sup>-1</sup> | Snow densification parameter (if DENSTY=2) |
| Talb | -2   | &deg;C| Albedo decay temperature threshold (if ALBEDO=0)           |
| tcld | 1000 | h    | Cold snow albedo decay time scale (if ALBEDO=1)             |
| tmlt | 100  | h    | Melting snow albedo decay time scale (if ALBEDO=1)          |
| trho | 200  | h    | Compaction time scale (if DENSTY=1)                         |
| Wirr | 0.03 | -    | Irreducible liquid water content (if HYDROL=1)              |
| z0sn | 0.01 | m    | Snow roughness length                                       |

### Site characteristics namelist &maps and map files

| Variable | Default | Units | Description |
|----------|---------|-------|-------------|
| alb0     | 0.2     | -     | Snow-free ground albedo                             |
| canh     | 2500*VAI| J/K/m<sup>-2</sup> | Canopy heat capacity (if CANMOD=1)     |
| fcly     | 0.3     | -     | Soil clay fraction                                  |
| fsnd     | 0.6     | -     | Soil sand fraction                                  |
| fsky     |     exp(-kext*VAI) | -    | Sky view fraction                         |
| fveg     | 1 - exp(-kext*VAI) | -    | Canopy cover fraction                     |
| hcan     | 0       | m     | Canopy height                                       |
| scap     | canc*VAI| kg m<sup>-2</sup> | Canopy snow capacity                    |
| VAI      | 0       | -     | Vegetation area index                               |
| z0sf     | 0.1     | m     | Snow-free roughness length                          |

Site characteristics can either be left as default values, set to a sequence of values in the namelist or read from a named map file. e.g. for a simulation with 10 points, the snow-free ground albedo can be reset to a constant value of 0.1 by including

    alb0 = 10*0.1

in namelist &maps, set to sequence by including

    alb0 = 0.2 0.2 0.1 0.1 0.1 0.1 0.1 0.1 0.2 0.2

or read from a file albedo.txt containing 10 values by including

    alb0_file = 'albedo.txt'

Note that sky view can be set independently of vegetation cover to allow for grid cells shaded by topography or vegetation in neighbouring cells.

### Initial values namelist &initial

| Variable   | Default | Units  | Description |
|------------|---------|--------|-------------|
| start_file | 'none'  | string | Start file  |
| fsat       | 4 * 0.5 | -      | Initial moisture content of soil layers as fractions of saturation |
| Tsoil      | 4 * 285 | K      | Initial temperatures of soil layers                                |

Soil temperature and moisture content are taken from the namelist and FSM2 is initialized in a snow-free state by default. If a start file is named, it should be a text file containing initial values for each of the state variables in order:

| Variable    |  Units             | Description |
|-------------|--------------------|-------------|
| albs(1:Nx,1:Ny)          |  -                  | Snow albedo                                |
| Ds(1:Nsmax,1:Nx,1:Ny)    |  m                  | Snow layer thicknesses                     |
| Nsnow(1:Nx,1:Ny)         |  -                  | Number of snow layers                      | 
| Qcan(1:Nx,1:Ny)          |  -                  | Canopy air space humidity                  |
| Sice(1:Nsmax,1:Nx,1:Ny)  |  kg m<sup>-2</sup>  | Ice content of snow layers                 |
| Sliq(1:Nsmax,1:Nx,1:Ny)  |  kg m<sup>-2</sup>  | Liquid content of snow layers              |
| Sveg(1:Nsmax,1:Nx,1:Ny)  |  kg m<sup>-2</sup>  | Canopy snow mass                           |
| Tcan(1:Nx,1:Ny)          |  K                  | Canopy air space temperature               |
| theta(1:Nsoil,1:Nx,1:Ny) |  -                  | Volumetric moisture content of soil layers |
| Tsnow(1:Nsmax,1:Nx,1:Ny) |  K                  | Snow layer temperatures                    | 
| Tsoil(1:Nsoil,1:Nx,1:Ny) |  K                  | Soil layer temperatures                    |
| Tsurf(1:Nx,1:Ny)         |  K                  | Surface skin temperature                   |
| Tveg(1:Nx,1:Ny)          |  K                  | Vegetation temperature                     |

The easiest way to generate a start file is to spin up the model by running it for a whole number of years without a start file and then rename the dump file produced at the end of the run as a start file for a new run.

### Output namelist &outputs and output files

Although still simple, FSM2 has more flexible output options than FSM. 

| Variable  | Default    | Description |
|-----------|------------|-------------|
| Nave      | 24         | Number of timesteps in averaged outputs |
| Nsmp      | 12         | Timestep of sample outputs, <= Nave     |
| runid     | none       | Run identifier string                   |
| dump_file | 'dump'     | Dump file name                          |

Flux variable are averaged over Nave timesteps and written to file ave, and state variables are written to file smp at timestep number Nsmp during every averaging period. For the default output frequencies, daily averages and samples at noon will be produced if the driving data has a one-hour timestep and starts at 01:00. Full timeseries are written if Nave = 1 and Nsmp = 1. At the end of a run, the state variables are written to a dump file with the same format as the start file. A run identifier, if specified, is prefixed on all output file names. If the run identifier includes a directory name (e.g. runid = 'output/'), the directory has to exist before the model is run.

The sample file has 4 + 3Nx*Ny columns:

| Variable       | Units              | Description           |
|----------------|--------------------|-----------------------|
| year           | years              | Year                  |
| month          | months             | Month of the year     |
| day            | days               | Day of the month      |
| hour           | hours              | Hour of the day       |
| snd(1:Nx*Ny)   | m                  | Snow depth            |
| SWE(1:Nx*Ny)   | kg m<sup>-2</sup>  | Snow water equivalent |
| Sveg(1:Nx*Ny)  | kg m<sup>-2</sup>  | Canopy snow mass      |


The average file has 3 + 7Nx*Ny columns:

| Variable       | Units              | Description           |
|----------------|--------------------|-----------------------|
| year           | years              | Year                  |
| month          | months             | Month of the year     |
| day            | days               | Day of the month      |
| alb(1:Nx*Ny)   | -                  | Flux-weighted albedo  |
| Gsurf(1:Nx*Ny) | W m <sup>-2</sup>  | Ground heat flux      |
| Hatmo(1:Nx*Ny) | W m <sup>-2</sup>  | Sensible heat flux    |
| Latmo(1:Nx*Ny) | W m <sup>-2</sup>  | Latent heat flux      |
| Melt(1:Nx*Ny)  | kg m <sup>-2</sup> | Cumulated melt        |
| Rnet(1:Nx*Ny)  | W m <sup>-2</sup>  | Net radiation         |
| Roff(1:Nx*Ny)  | kg m <sup>-2</sup> | Cumulated runoff      |
| Tsurf(1:Nx*Ny) | C                  | Surface temperature   |
| Tsoil(1:Nx*Ny) | C                  | 20 cm soil temperature|

A metadata file runifo is produce containing copies of all the namelists and the physics options for the run.
 

## References

<a name="Essery2015"></a> Essery (2015). A Factorial Snowpack Model (FSM 1.0). *Geoscientific Model Development*, **8**, 3867-3876, [doi:10.5194/gmd-8-3867-2015](http://www.geosci-model-dev.net/8/3867/2015/)

<a name="Essery2016"></a> Essery et al. (2016). A 7-year dataset for driving and evaluating snow models at an Arctic site (Sodankylä, Finland). *Geosci. Instrum. Method. Data Syst.*, **5**, 219-227, [doi:10.5194/gi-5-219-2016](https://www.geosci-instrum-method-data-syst.net/5/219/2016/)

