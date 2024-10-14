# VASP_PDOS_PMG
Reads and plots density of states from a VASP calculation from command line. Works for both non/ and spin-polarised calulations. 
Reads data from the ```vasprun.xml``` file using pymatgen [[1]](#References) and outputs plot as Portable Network Graphic (PNG).

See below for usage or type ```pdos.py -h``` or ```--help```.




## Usage:

```
pdos.py [-h] [-x XLIM XLIM] [-y YLIM YLIM] [-r RES] [-g GAUSS] [-f] [-s SPECIES [SPECIES ...] | -o ORBITALS [ORBITALS ...] | -a ATOMS [ATOMS ...]] vasprun
```
Plots the (projected) density of states with fermi energy shifted to 0. If no elements or orbital types are specified for projection, it will plot the total
density of states. Requires POSCAR/CONTCAR and vasprun.xml files from a VASP calculation.

*positional arguments:*

 * ```vasprun``` - (**required**) filename or path to vasprun XML-file.

*optional arguments:*

  * ```-h, --help```:             show help message and exit

  * ```-s SPECIES [SPECIES ...], --species SPECIES [SPECIES ...]```:
                         list of species for projections (no orbital projection) e.g. ```-s Ca Al``` etc.

  * ```-o ORBITALS [ORBITALS ...], --orbitals ORBITALS [ORBITALS ...]```:
                         list of specific species and orbital projection e.g. ```-o Ca s Al s Al p ```etc.

  * ```-a ATOMS [ATOMS ...], --atoms ATOMS [ATOMS ...]```:
                         list of specific atomic sites and orbital projection e.g. ```-o Ca1 s Al3 s Al3 p``` etc.

  * ```-x XLIM XLIM, --xlim XLIM XLIM```: 
                         upper and lower x-limits (energy) of plot window e.g. ```-x -10 5``` etc.

  * ```-y YLIM YLIM, --ylim YLIM YLIM```
                         upper and lower y-limits (DOS) of plot window e.g. ```-y 0 50``` etc.

  * ```-r RES, --res RES  ```: 
                         resolution of png image (in dpi).

  * ```-g GAUSS, --gauss GAUSS```:
                         sigma for gaussian smearing (in eV). Defaults to no smearing.

  * ```-f, --fill ```: 
                         optional flag to fill the area under the curve. Stacks pDOS's in the order that they are provided.

## Important notes: 


Plotting any orbital-projected DOS requires a VASP calculation with ```LORBIT = 11```.

This code has not been tested for non-collinear calulations yet.

Pymatgen (up to) version 2023.12.18 is required.                

## References:

[1] Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier,
Michael Kocher, Shreyas Cholia, Dan Gunter, Vincent Chevrier, Kristin A.
Persson, Gerbrand Ceder. *Python Materials Genomics (pymatgen) : A Robust,
Open-Source Python Library for Materials Analysis.* Computational Materials
Science, 2013, 68, 314–319. https://doi.org/10.1016/j.commatsci.2012.10.028.

Pymatgen is released under the MIT license. See https://pymatgen.org/ for more info.