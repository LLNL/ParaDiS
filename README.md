
# ParaDiS

**IMPORTANT: ParaDiS development is no longer supported as the code is being replaced by our new open-source project [OpenDiS](https://github.com/OpenDiS/OpenDiS), which includes the high-performance [ExaDiS](https://github.com/LLNL/exadis) core library that runs on GPU. More information about OpenDiS is available at the [OpenDiS documentation](https://opendis.github.io/OpenDiS).**

ParaDiS Public Release Version 4.0

ParaDiS (Parallel Dislocation Simulator) is a simulation tool that performs direct numerical simulation of dislocation ensembles, the carriers of plasticity, to predict the strength in crystalline materials from the fundamental physics of defect motion, evolution, and interaction.

The code has been successfully deployed on high performance computing architectures and used to study the origins of strength and strain hardening for cubic crystals, the strength of micro-pillars, and irradiated materials at LLNL. The ParaDiS code has been successfully deployed on more than one hundred thousand CPU's with over ten million active degrees of freedom.

## Installation

### Quick start

The public distribution of ParaDiS is built using conventional make files. 

1. Select a target system (e.g. `SYS=gcc` at line 52) and compilation mode (`MODE=SERIAL` or `MODE=PARALLEL` at lines 66-67) in file `makefile.setup`
2. Compile the code from the root of the distribution using the `make` command
3. Successful compilation will produce binary executables in the `./bin` directory
4. Test your installation by running an example
```
./bin/paradis tests/frank_read_src.ctrl
```

### Detailed instructions

The build system is comprised of a number of makefiles:

* `makefile`           : overall makefile to build ParaDiS and supporting utilities
* `makefile.setup`     : use to enable/disable various application features
* `makefile.sys`       : system specific build settings for various supported systems
* `makefile.srcs`      : complete list of all source files used to build system
* `src/makefile`       : makefile that builds ParaDiS 
* `src/makefile.dep`   : makefile dependencies (auto-generated via make depend)
* `tests/makefile`     : makefile for cleaning up and managing the test directory
* `utils/makefile`     : makefile to building the various utility applications
* `ext/makefile`       : makefile for building any external dependencies

If you want or need to customize the build for your machine, you may need to adjust the two 
files - `makefile.sys` and `makefile.setup` before executing `make`.

Executing make from the main directory will build the main ParaDiS executable
and all the supporting tools. If you only want to build the main ParaDiS application, 
execute `make` from the main source directory (`src/`).

#### makefile.sys

Prior to building ParaDiS, you will need to identify the target system 
and environment for the application. Several target systems are preset and
are listed below:

```
linux.intel     Linux systems using native Intel compilers.
linux           Generic linux system
gcc             Generic system build using gnu compilers
aix             IBM aix systems using native compilers (LLNL's ice, berg, purple, um, up, uv...)
mac             MacBook Pro 
bgp             LC BlueGene/P systems (dawn, dawndev)
bgq             LC BlueGene/Q systems (sequoia, rzuseq)
mc-cc           Stanford ME Linux system using intel compilers
wcr             Stanford ME Linux system using intel compilers
cygwin          Stanford Linux emulator for Windows PC
xt4             Cray XT4 systems (NERSC's franklin)
```

If you are attempting to build the ParaDiS software on a system that is not listed above,
you may need to copy and/or adjust an existing system in the `makefile.sys` file and
set the `SYS=` parameter in `makefile.setup`.

#### makefile.setup

All the user-specific features and settings are enabled/disabled via the 
`makefile.setup` file.  Here you can specify whether the application will 
run serially or parallel (via MPI), set the optimization level, enable the
X-Windows display, etc.  Details on all the switches are documented in the
file.

The main parameters you need to select are:

```
# required parameter, must identify the target host machine
SYS=[linux.intel | linux | gcc | aix| mac | bgp| bgq | mc-cc| wcr| cygwin | xt4 ]

MODE=SERIAL     # sets the execution mode to serial
MODE=PARALLEL   # sets the execution mode to parallel (default, requires an MPI installation) 

XLIB_MODE=ON    # enables  the X-Window visualization (do not use for production runs)
XLIB_MODE=OFF   # disables the X-Window visualization (default)
```

All other parameters and settings are detailed in the `makefile.setup` file.

## Directory Structure

Brief description of the directories within this distribution:

* `./bin`      : executable applications, created during build
* `./src`      : C/C++ source files (*.cc)
* `./include`  : C/C++ include files (*.h)
* `./obj`    : object file directories (parallel and serial builds)
* `./docs`     : support documentation 
* `./inputs`   : generic input files (Rijm tables, FMM tables, gnuplot files, and X-Windows defaults)
* `./tests`    : example tests (`*.ctrl`, `*.data`, `*.sh`)
* `./utils`    : support utilities
* `./tools`    : support tools

## Applications and Tools

The following applications and support utilities are created when building 
the ParaDiS application from the root location:

* `./bin/paradis`           : main ParaDiS simulator
* `./bin/paradisconvert`    : conversion utility for older (unsupported) control and data files
* `./bin/calcdensity`       : dislocation density calculator
* `./bin/ctablegen`        : utility for creating PBC image correction tables
* `./bin/ctablegenp`        : utility for creating PBC image correction tables (parallel/MPI version)
* `./bin/paradisgen`        : utility for creating initial, random dislocation networks
* `./bin/paradisrepart`     : utility for repartitioning existing domain decompositions
* `./bin/stresstablegen`    : utility for creating far-field stress tables

## Simulation Examples

There are numerous examples of simulation input files located in the 
`tests/` directory. ParaDiS is mainly controlled by two input files, a `.data` file and a `.ctrl` file.
The `data` file specifies the domain size and initial dislocation network.
The `ctrl` file specifies the materials parameters, loading procedure, and numerical parameters of the simulation.

Most of the control files in the `./tests` directory are setup for parallel execution by 8 processors. This enables
the tests to be run on a contemporary multicore desktop. The total number of MPI processes
is control by the domain parameters:

```
numXdoms = 2  # number of X-axis domains
numYdoms = 2  # number of Y-axis domains
numZdoms = 2  # number of Z-axis domains
```

The total domain count (processes) = (`numXdoms` * `numYdoms` * `numZdoms`)
Note that the total domain count identified in the control file needs to be consistent with the 
number of MPI tasks/processes that are launched when the application is run. The command line 
for launching an MPI-enabled run of ParaDiS using the above domain configuration would be,

```
mpirun -n 8 ./bin/paradis ./tests/mg-cAxis.ctrl
```

Important note: several of the test files use precomputed FMM tables.  Those tables are located
in the `./inputs` directory. The paths to these files are controlled in the ParaDiS control files
assume that the ParaDiS application is started from the main ParaDiS directory (e.g. `./bin/paradis`).
If you do not execute ParaDiS in this way, you must specify the paths to the input control, data,
and any auxiliary files relative to where you are running the application.

Some of the examples include Fast Multipole Method (FMM) expansion tables.  You can create FMM image 
correction tables using the `ctablegen` utility. 

To generate an FMM table for an isotropic simulation: 
```
  ./bin/ctablegen -nu 3.327533e-01 -mu 6.488424e+10 -mporder 2 -torder 5 -outfile inputs/fmm-ctab.data
```
See the ParaDiS users guide for more details on building FMM correction tables.

Initial dislocation networks can be created using the `paradisgen` application.

## Citation

```
@article{arsenlis2007enabling,
  title={Enabling strain hardening simulations with dislocation dynamics},
  author={Arsenlis, Athanasios and Cai, Wei and Tang, Meijie and Rhee, Moono and Oppelstrup, Tomas and Hommes, Gregg and Pierce, Tom G and Bulatov, Vasily V},
  journal={Modelling and Simulation in Materials Science and Engineering},
  volume={15},
  number={6},
  pages={553},
  year={2007},
  publisher={IOP Publishing}
}
```

## License

ParaDiS is released under the BSD-3 license. See [LICENSE](LICENSE) for details.

LLNL-CODE-853453
