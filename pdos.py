#!/usr/bin/env python3

"""
Author: P.Dirk

Reads and plots density of states from a VASP calculation from command line. Works for both non/ and spin-polarised calulations. 
Reads data from the vasprun.xml file using pymatgen[1] and outputs plot as Portable Network Graphic (PNG).

See below for usage or type pdos.py -h or --help.

Important notes: 
---------------------------------

Plotting any orbital-projected DOS requires a VASP calculation with LORBIT = 11

This code has not been tested for non-collinear calulations yet.

Pymatgen (up to) version 2023.12.18 is required.

--------------------------------- September 2024


usage: pdos.py [-h] [-s SPECIES [SPECIES ...] | -o ORBITALS [ORBITALS ...] | -a ATOMS [ATOMS ...]] [-x XLIM XLIM] [-y YLIM YLIM] vasprun

plots the (projected) density of states with fermi energy shifted to 0. If no elements or orbital types are specified for projection, it will plot the total
density of states. Requires POSCAR and vasprun.xml files from a SCF calculation.

positional arguments:

  vasprun - provide name of xml-file (usually just vasprun.xml or path)

optional arguments:

  -h, --help            show this help message and exit

  -s SPECIES [SPECIES ...], --species SPECIES [SPECIES ...]
                        <optional> list of species for projections (no orbital projection) e.g. -s Ca Al etc.

  -o ORBITALS [ORBITALS ...], --orbitals ORBITALS [ORBITALS ...]
                        <optional> list of specific species and orbital projection e.g. -o Ca s Al s Al p etc.

  -a ATOMS [ATOMS ...], --atoms ATOMS [ATOMS ...]
                        <optional> list of specific atomic sites and orbital projection e.g. -o Ca1 s Al3 s Al3 p etc.

  -x XLIM XLIM, --xlim XLIM XLIM
                        <optional> upper and lower x-limits of DOS plot window e.g. -x -10 5 etc.

  -y YLIM YLIM, --ylim YLIM YLIM
                        <optional> upper and lower y-limits of DOS plot window e.g. -y 0 50 etc.

  -r RES, --res RES     
                        <optional> resolution of png image (in dpi)

  -g GAUSS, --gauss GAUSS
                        <optional> sigma for gaussian smearing (in eV)

  -f, --fill  
                        <optional> fill the area under the curve
                        

References:

[1] Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier,
Michael Kocher, Shreyas Cholia, Dan Gunter, Vincent Chevrier, Kristin A.
Persson, Gerbrand Ceder. *Python Materials Genomics (pymatgen) : A Robust,
Open-Source Python Library for Materials Analysis.* Computational Materials
Science, 2013, 68, 314–319. https://doi.org/10.1016/j.commatsci.2012.10.028.

Pymatgen is released under the MIT license. See https://pymatgen.org/ for more info.                       
"""

from pymatgen.io.vasp.outputs import Vasprun # type: ignore
from pymatgen.electronic_structure.core import OrbitalType # type: ignore
from pymatgen.electronic_structure.plotter import DosPlotter # type: ignore
import pymatgen.core.structure as mg # type: ignore
from pymatgen.core.periodic_table import Element # type: ignore
import os
import argparse
import re
from collections import defaultdict
import matplotlib as mpl # type: ignore


def check_species(species,check) -> None:
    if any(sp in species for sp in check) == False: # check if species names are valid
        print('at least one requested species not found in POSCAR')
        return SystemExit(1)
    else:
        return None

def check_orbitals(orbitals,check) -> None:
    if any(orb in orbitals for orb in check) == False: # check if orbital type names are valid
        print('at least one requested orbital type not valid')
        return SystemExit(1)
    else:
        return None
    
def check_labels(labels,check) -> None:
    if any(l in labels for l in check) == False: # check if site labels are valid
        print('at least one requested site label not valid')
        return SystemExit(1)
    else:
        return None



def enumerate_site_labels(labels) -> list:
    # count occurrences of each label
    counts = defaultdict(int)
    for l in labels:
        counts[l] += 1

    current_index = defaultdict(int) # track current index for each string
    result = []

    for l in labels:
        if counts[l] > 1:
            current_index[l] += 1
            new_label = f"{l}{current_index[l]}"

        result.append(new_label)
    
    return result


def get_usr_args():
    parser = argparse.ArgumentParser(
        description="plots the (projected) density of states with fermi energy shifted to 0. If no elements or orbital types are specified for projection, it will plot the total density of states. Requires POSCAR and vasprun.xml files from a SCF calculation."
        )
    parser.add_argument(
        "vasprun",
        help = '<required> filename or path to vasprun XML-file.',
        type = str
        )
    parser.add_argument(
        "-x",
        "--xlim",
        nargs = 2,
        help = '<optional> upper and lower x-limits of DOS plot window e.g. -x -10 5 etc.',
        default = (None, None),
        type = float
        )
    parser.add_argument(
        "-y",
        "--ylim",
        nargs = 2,
        help = '<optional> upper and lower y-limits of DOS plot window e.g. -y 0 50 etc.',
        default = (None, None),
        type = float
        )
    parser.add_argument(
        "-r",
        "--res",
        nargs = 1,
        help = '<optional> resolution of png image (in dpi)',
        default = [400],
        type = float
        )
    parser.add_argument(
        "-g",
        "--gauss",
        nargs = 1,
        help = '<optional> sigma for gaussian smearing (in eV)',
        default = None,
        type = float        
    )
    parser.add_argument(
        "-f",
        "--fill",
        action = 'store_true',
        help = '<optional> fill the area under the curve',
        default = False     
    )
    # options for which data to plot
    group = parser.add_mutually_exclusive_group(
        required = False
        )
    group.add_argument(
        "-s",
        "--species",
        nargs = '+',
        help = '<optional> list of species for projections (no orbital projection) e.g. -s Ca Al etc.'
        )
    group.add_argument(
        "-o",
        "--orbitals", 
        nargs = '+',
        help = '<optional> list of specific species and orbital projection e.g. -o Ca s Al s Al p etc.'
        )
    group.add_argument(
        "-a",
        "--atoms", 
        nargs = '+', 
        help = '<optional> list of specific atomic sites and orbital projection e.g. -o Ca1 s Al3 s Al3 p etc.'
        )
    
    args = parser.parse_args()

    vasprun_file : str = args.vasprun
    species_request : list | None = args.species
    orbitals_request : list | None = args.orbitals
    sites_request : list | None = args.atoms
    resolution : list[float] = args.res
    sigma : float | None = args.gauss
    fill : bool = args.fill
    limits_x : tuple[float,float] = (
        args.xlim[0],
        args.xlim[1]
        )
    limits_y : tuple[float,float] = (
        args.ylim[0],
        args.ylim[1]
        )

    return vasprun_file, species_request, orbitals_request, sites_request, resolution, limits_x, limits_y, sigma, fill



def parse_vasprun(vasprun_file: str, complete_dos: bool = True):

    if not os.path.isfile(str(vasprun_file)): # check if vasprun file exists

        print("The file \"{}\" doesn't exists".format(vasprun_file))
        raise SystemExit(1)

    else: # extract DOS data from vasprun file

        dos_vasprun=Vasprun(str(vasprun_file),parse_eigen=True)

        if complete_dos == False:
            total_dos=dos_vasprun.tdos
            return total_dos
        else:
            dos_data=dos_vasprun.complete_dos
            return dos_data


def get_struct_data():

    if not os.path.isfile("CONTCAR"): # check if CONTCAR file exists

        print("WARNING: The CONTCAR file doesn't exists, using POSCAR instead. This will probably crash for ionic relaxation runs!")

        if not os.path.isfile("POSCAR"): # check if POSCAR file exists 

            print("The POSCAR file doesn't exists.")
            raise SystemExit(1)

        else:  # read site positions and labels from POSCAR (this fails for relaxation runs, since positions are incompatible with the ones in vasprun.xml)

            species_list = mg.Structure.from_file("POSCAR").symbol_set
            sites=mg.IStructure.from_file("POSCAR").sites
            labels_allowed=enumerate_site_labels(mg.IStructure.from_file("POSCAR").labels)

    else: # extract spcies names from CONTCAR file

        species_list = mg.Structure.from_file("CONTCAR").symbol_set
        sites=mg.IStructure.from_file("CONTCAR").sites
        labels_allowed=enumerate_site_labels(mg.IStructure.from_file("CONTCAR").labels)
    
    return species_list, sites, labels_allowed


orbitals_allowed = { 
    's': OrbitalType.s,
    'p': OrbitalType.p,
    'd': OrbitalType.d,
    'f': OrbitalType.f
    }


if __name__ == '__main__':

    vasprun_file, species_request, orbitals_request, sites_request, resolution, limits_x, limits_y, smearing, fill = get_usr_args()

    species_list, sites, labels_allowed = get_struct_data() # TODO add option to view structural/calculation info before plotting

    # initialise plotter
    plotter=DosPlotter(
        zero_at_efermi = True,
        stack = fill,
        sigma = smearing
    )

    if species_request == None and orbitals_request == None and sites_request == None: # -s,-o and -a flags not specified, plot only total DOS

        print('No projections specified. Plotting total DOS...')
        dos_data = parse_vasprun(vasprun_file,complete_dos=False)
        plotter.add_dos('total', dos_data)

    else:
        if species_request != None: # only -s flags specified (element projected DOS)

            print('plotting element-projected DOS...')
            check_species(species_request,species_list)
            dos_data = parse_vasprun(vasprun_file)
            pdos = dos_data.get_element_dos()

            for el in species_request: # plot pDOS
                plotter.add_dos(el, pdos[Element(el)])

        elif orbitals_request != None: # only -o flags specified (element-spd projected DOS)

            print('plotting element-orbital-projected DOS...')
            elements = orbitals_request[::2]
            orbitals = orbitals_request[1::2]

            check_species(elements,species_list)
            check_orbitals(orbitals,list(orbitals_allowed.keys()))

            dos_data = parse_vasprun(vasprun_file)

            for i,el in enumerate(elements): # plot pDOS

                pdos = dos_data.get_element_spd_dos(el)
                plotter.add_dos('{}({})'.format(el,orbitals[i]), pdos[orbitals_allowed[orbitals[i]]])

        elif sites_request != None: # only -a flags specified (atom-site spd-projected DOS)

            print('plotting site-resolved orbital-projected DOS')
            site_labels = sites_request[::2]
            orbitals = sites_request[1::2]

            # separate element info drom site number info
            elements = ["".join(re.split("[^a-zA-Z]*", site)) for site in site_labels]
            site_numbers = [int("".join(re.split("[^1-9]*", site))) for site in site_labels]

            check_species(elements,species_list)
            check_orbitals(orbitals,list(orbitals_allowed.keys()))
            check_labels(site_labels,labels_allowed)

            dos_data = parse_vasprun(vasprun_file)

            for i,l in enumerate(site_labels): # plot pDOS
                site_index = labels_allowed.index(site_labels[i])
                pdos = dos_data.get_site_spd_dos(sites[site_index])
                plotter.add_dos('{}({})'.format(l,orbitals[i]), pdos[orbitals_allowed[orbitals[i]]])


    # print fermi energy, band gap, and conduction and valence band positions
    fermi_energy=dos_data.efermi
    print('fermi energy: ', fermi_energy)
    cbm,vbm=dos_data.get_cbm_vbm()
    print('band gap (eV): %.5f'%dos_data.get_gap())
    print('VBM (eV) : ',vbm, ' CBM (eV) : ',cbm)

    # produce plot and save as png file
    mpl.rcParams['figure.dpi'] = resolution[0] # TODO make plots look prettier
    plotter.get_plot()
    plotter.save_plot('DOS.png',ylim=limits_y,xlim=limits_x) # TODO alternative file formats?
