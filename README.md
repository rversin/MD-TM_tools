# MD-TM_tools

> Studies transmembrane helix properties.

## Features

calc_angles.py can :
- calculate the crossing angle between two transmembrane helices for a simple structure, or a trajectory
- calculate the tilt angle for a simple structure, or a trajectory

calc_resi-dist.py can :
- calculate all the distances for a defined cutoff between the atoms of two transmembrane helices
  - usefull to determine a cutoff in order to then build a contact map
- calculate the frequencies of each found couple of residues over a trajectory
  - usefull to build a contact map

## Requirements

Python >= 3.6 is mandatory for running the scripts.

The script are written in Python 3 and need the following modules :
  - argparse
  - numpy
  - math
  - sys
  - MDAnalysis.

## Installation (development)

1. Install conda (either with Miniconda or Anaconda, we recommend Miniconda)

2. Clone this GitHub repository:
```
$ git clone https://github.com/rversin/MD-TM_tools.git
```

3. Create conda environment:
```
$ conda env create -f MD-TM_tools.yml
$ conda activate MD-TM_tools
```

## Usage

### For calc_angles.py

```
usage: calc_angles.py [-h] [-xtc XTC] [-mb] [-bh1 BH1] [-eh1 EH1] [-bh2 BH2] [-eh2 EH2] [-z Z] gro

This script reads a structure file (pdb or gro file) and a gromacs trajectory file
(.xtc, is optional), and calculate either the crossing-angle between two helices,
or the tilt angle.
The results will be shown in the shell.

Get the help :
python3 calc_angles.py -h

Usage for crossing angle:
1) for a whole trajectory
python3 calc_angles.py structure.gro -xtc trajectory.xtc -bh1 1 -eh1 20 -bh2 21 -eh2 42 > CA.dat
2) for only one structure
python3 calc_angles.py structure.gro -bh1 1 -eh1 20 -bh2 21 -eh2 42 > CA.dat

Usage for tilt angle:
1) for a whole trajectory
python3 calc_angles.py structure.gro -xtc trajectory.xtc -bh1 1 -eh1 20 -mb > TILT.dat
python3 calc_angles.py structure.gro -xtc trajectory.xtc -bh1 21 -eh1 42 -mb -z -20 > TILT.dat
2) for only one structure
python3 calc_angles.py structure.gro -bh1 1 -eh1 20 -mb > TILT.dat
python3 calc_angles.py structure.gro -bh1 21 -eh1 42 -mb -z -20 > TILT.dat
WARNING : the tilt angle calculation assumes that the bilayer normal is aligned along Z.

positional arguments:
  gro         .gro file

optional arguments:
  -h, --help  show this help message and exit
  -xtc XTC    .xtc file, if "none" only the gro file is taken into account
  -mb         If in command line, calculates the angle of the helix with the membrane. If not,
              calculates the crossing angle between the two helix.
  -bh1 BH1    Number of the first residue for the first helix
  -eh1 EH1    Number of the last residue for the first helix
  -bh2 BH2    Number of the first residue for the second helix
  -eh2 EH2    Number of the last residue for the second helix
  -z Z        Define the Z component of the axis normal to the membrane

```

### For calc_resi-dist.py

```
usage: calc_resi-dist.py [-h] [-xtc XTC] [-d D] [-time_min TIME_MIN] [-time_max TIME_MAX] [-p] [-a]
                         [-ele {bb,all}] [-bh1 BH1] [-eh1 EH1] [-bh2 BH2] [-eh2 EH2]
                         gro

This script reads a structure file (pdb or gro file) and a gromacs trajectory file
(.xtc, is optional), and can either calculate the distances for a definded cutoff
between all the defined atoms of the helix or calculates the contact statistics
used then to build contact maps.

To build a contact maps :
python3 calc_resi-dist.py structure.gro -xtc trajectory.xtc -a -d 7 -bh1 5 -eh1 25 -bh2 36 -eh2 56 -a

positional arguments:
  gro                 .gro file

optional arguments:
  -h, --help          show this help message and exit
  -xtc XTC            .xtc file, if "none" only the gro file is taken into account
  -d D                distance cutoff
  -time_min TIME_MIN  minimal time cutoff
  -time_max TIME_MAX  maximum time cutoff, the default value is at10 microsec
  -p                  prints all the distances
  -a                  make contact analysis
  -ele {bb,all}       elements checked for distances, bb id for backbone atomsall is for all the
                      atoms
  -bh1 BH1            Number of the first residue for the first helix
  -eh1 EH1            Number of the last residue for the first helix
  -bh2 BH2            Number of the first residue for the second helix
  -eh2 EH2            Number of the last residue for the second helix

```

## Contributors

  - Raphaelle Versini
  - Patrick Fuchs
  - Antoine Taly
