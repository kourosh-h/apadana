import mbuild as mb
from mbuild.lib.recipes import Monolayer
from mbuild.lib.recipes import Alkane
from mbuild.lib.moieties import Silane
from mbuild.lib.atoms import H
from mbuild.lib.recipes import TiledCompound
from mbuild import Port
from foyer import Forcefield
import parmed as pmd
from parmed import gromacs, amber, unit as u
#import openbabel
# import pybel
import py3Dmol
# import nglview as nv
#from ase import Atoms
#from ase.io.trajectory import Trajectory
#from ase.io import read
import os
import random
import numpy as np
import pandas as pd
from IPython.display import display, Latex
from .fileReading import *
from .fileWriting import *
from .forcefields import *
from .functionalizer import *
from .modellingMolecules import *
from pathlib import Path

# Setup package directories
apadana_package_directory = Path(__file__).parent
interface_charmm_xml = apadana_package_directory / 'EMAMI_CHARMM.xml'
ADDITIONAL_FILES = apadana_package_directory / "additional_files"



class AlphaQuartz():
    """
    Creates or loads a silica surface (Alpha-Quartz).
    Supports three modes:
    1. 'cell': Build from default unit cell (tiling enabled).
    2. 'slab': Load default relaxed slab from package (no tiling).
    3. 'user': Load custom GRO/PDB files provided by the user.
    """
    def __init__(self, nameSilicaSur = 'SilicaSUR',
                 x_blocks = 1,
                 y_blocks = 1,
                 freespace = 6,
                 res_name = 'SUR',
                 SIinsilica = "SI",
                 bOinsilica = "OS",
                 sOinsilica = "OD",
                 sHinsilica = "HD"):
        """
        :param nameSilicaSur: name of the silica surface
        :param x_blocks: number of repeats in x dimension
        :param y_blocks:  number of repeats in y dimension
        :param freespace: setting zhi and zlo; the distance in nm
        :param res_name: residue name of silica surface; str
        :param SIinsilica: Si atom name
        :param bOinsilica: O atom name
        :param sOinsilica: OH atom name
        :param sHinsilica: HO atom name
        """
        self.x_blocks = x_blocks
        self.y_blocks = y_blocks
        self.ressilica = res_name
        self.freespace = freespace * 2
        self.SIinsilica = SIinsilica
        self.bOinsilica = bOinsilica
        self.sOinsilica = sOinsilica
        self.sHinsilica = sHinsilica
        self.nameSilicaSur = nameSilicaSur
        
        # Default internal package file paths
        self.pkg_cell_gro = ADDITIONAL_FILES / "unitcell.gro"
        self.pkg_cell_pdb = ADDITIONAL_FILES / "unitcell.pdb"
        self.pkg_slab_gro = ADDITIONAL_FILES / "slab.gro"
        self.pkg_slab_pdb = ADDITIONAL_FILES / "slab.pdb"


    def get_dimensions(self, pdb_path):
        """
        Reads dimensions from a PDB file. 
        - lenx, leny, lenz are taken from the CRYST1 line.
        - lenzS is calculated from the range of the Z-coordinates.
        """
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"PDB file not found at: {pdb_path}")

        lenx = leny = lenz = 0.0
        z_coords = []

        with open(pdb_path, 'r') as f:
            for line in f:
                # Extract Box Dimensions from CRYST1 line
                if line.startswith("CRYST1"):
                    parts = line.split()
                    # PDB CRYST1 values are in Angstroms
                    lenx = float(parts[1])
                    leny = float(parts[2])
                    lenz = float(parts[3])

                # Extract Z-coordinates from ATOM or HETATM lines
                elif line.startswith("ATOM") or line.startswith("HETATM"):
                    # In PDB format, Z-coordinate is typically from index 46 to 54
                    # Or split and take index 7 (if the file is well-spaced)
                    z_val = float(line[46:54].strip())
                    z_coords.append(z_val)

        if not z_coords:
            raise ValueError(f"No coordinate data found in {pdb_path}")

        # Calculate thickness (lenzS) same as the old function
        max_z = max(z_coords)
        min_z = min(z_coords)
        lenzS = round(abs(max_z - min_z), 3)

        # Calculate empty space (emp) for centering logic
        # (lenz - thickness) / 2 / 10 to convert back to nm for translation logic
        emp = (lenz - lenzS) / 2 / 10 

        return lenx, leny, lenz, emp, lenzS


    def make_unitcell_pdb(self, input_pdb, input_gro):

        """
        Internal helper: Renames atoms in the unit cell PDB to match Forcefield names.
        """
        outputPDBname = 'unitcell.pdb'
        strs = [f"{x}" for x in range(100)]

        lenx, leny, lenz, emp, lenzS = self.get_dimensions(input_pdb)

        with open(input_pdb,
                  'r') as pdb_file:
            with open(f'{outputPDBname}', 'w+') as new_pdb:
                new_pdb.writelines(
                    'CRYST1 {:>.3f}   {:>.3f}   {:>.3f}  90.00  90.00  90.00 P 1           1 \n'.format(lenx, leny,
                                                                                                            lenz))
                for line in pdb_file:
                    segments = line.split()
                    if segments[0] == 'ATOM':
                        segments[4] = '1'
                        segments[3] = self.ressilica
                        seg_n = 'X'
                        if segments[2][0] == "S":
                            segments[2] = self.SIinsilica
                            seg_ele = 'Si'
                        if segments[2][:2] == "OS": # and (segments[2][1]) in strs:
                            segments[2] = self.bOinsilica
                            seg_ele = 'O'
                        if segments[2][:2] == "OD":
                            segments[2] = self.sOinsilica
                            seg_ele = 'O'
                        if segments[2][:2] == "HD":
                            segments[2] = self.sHinsilica
                            seg_ele = 'H'

                        new_pdb.writelines(
                            "{:>s}{:>7s}{:>5s}{:>4s}{:>2s}{:>4s}{:>12s}{:>8s}{:>8s}{:>6s}{:>6s}{:>12s}\n".format(
                                segments[0], segments[1], segments[2], segments[3], seg_n, segments[4], segments[5],
                                segments[6], segments[7], segments[8], segments[9], seg_ele, ))
                new_pdb.writelines("END\n")

    def make_silica(self, mode='cell', user_gro=None, user_pdb=None):
        """
        Main method to provide the silica surface object.
        - mode='cell': Default small unit cell tiled by x_blocks and y_blocks.
        - mode='slab': Pre-existing large relaxed slab from package files.
        - mode='user': Custom files provided via user_gro and user_pdb.
        """
        
        if mode == 'cell':
            self.make_unitcell_pdb(self.pkg_cell_pdb, self.pkg_cell_pdb)
            outputPDBname = 'unitcell.pdb'
            lenx, leny, lenz, emp, lenzS = self.get_dimensions(self.pkg_cell_pdb)
            silicasurface = mb.Compound()
            silicasurface.name = self.ressilica
            silicasurface.resname = self.ressilica

            # the unitcell that we built before is loaded here
            UnitCell = mb.load(outputPDBname, compound=silicasurface)

            silicasurface.periodicity = (True, True, True)

            silicasurface.translate(np.array([0.0, 0.0, -emp]))
            silicasurface.box = mb.Box([lenx / 10, leny / 10, (lenzS / 10) + (self.freespace)])
            silicasurface.translate(np.array([0.0, 0.0, (self.freespace / 2) + (lenzS / 10 / 20)]))

            silicasurface = TiledCompound(silicasurface, n_tiles=(self.x_blocks, self.y_blocks, 1))

            silicasurface.generate_bonds(name_a=self.SIinsilica, name_b=self.bOinsilica, dmin=0.0,
                                         dmax=0.19)  
            silicasurface.generate_bonds(name_a=self.SIinsilica, name_b=self.sOinsilica, dmin=0.0,
                                         dmax=0.19)  
            silicasurface.generate_bonds(name_a=self.sOinsilica, name_b=self.sHinsilica, dmin=0.0,
                                         dmax=0.120)  

            #return silicasurface

        if mode == 'slab':
            self.make_unitcell_pdb(self.pkg_slab_pdb, self.pkg_slab_pdb)
            outputPDBname = 'unitcell.pdb'
            lenx, leny, lenz, emp, lenzS = self.get_dimensions(self.pkg_slab_pdb)
            silicasurface = mb.Compound()
            silicasurface.name = self.ressilica
            silicasurface.resname = self.ressilica

            # the unitcell that we built before is loaded here
            UnitCell = mb.load(outputPDBname, compound=silicasurface)

            silicasurface.periodicity = (True, True, True)

            silicasurface.translate(np.array([0.0, 0.0, -emp]))
            silicasurface.box = mb.Box([lenx / 10, leny / 10, (lenzS / 10) + (self.freespace)])
            silicasurface.translate(np.array([0.0, 0.0, (self.freespace / 2) + (lenzS / 10 / 20)]))

            silicasurface = TiledCompound(silicasurface, n_tiles=(self.x_blocks, self.y_blocks, 1))

            silicasurface.generate_bonds(name_a=self.SIinsilica, name_b=self.bOinsilica, dmin=0.0,
                                         dmax=0.19)  
            silicasurface.generate_bonds(name_a=self.SIinsilica, name_b=self.sOinsilica, dmin=0.0,
                                         dmax=0.19)  
            silicasurface.generate_bonds(name_a=self.sOinsilica, name_b=self.sHinsilica, dmin=0.0,
                                         dmax=0.120) 

            #return silicasurface

        if mode == 'user':
            self.make_unitcell_pdb(user_pdb, user_gro)
            outputPDBname = 'unitcell.pdb'
            lenx, leny, lenz, emp, lenzS = self.get_dimensions(user_pdb)
            silicasurface = mb.Compound()
            silicasurface.name = self.ressilica
            silicasurface.resname = self.ressilica

            # the unitcell that we built before is loaded here
            UnitCell = mb.load(outputPDBname, compound=silicasurface)

            silicasurface.periodicity = (True, True, True)

            silicasurface.translate(np.array([0.0, 0.0, -emp]))
            silicasurface.box = mb.Box([lenx / 10, leny / 10, (lenzS / 10) + (self.freespace)])
            silicasurface.translate(np.array([0.0, 0.0, (self.freespace / 2) + (lenzS / 10 / 20)]))

            silicasurface = TiledCompound(silicasurface, n_tiles=(self.x_blocks, self.y_blocks, 1))

            silicasurface.generate_bonds(name_a=self.SIinsilica, name_b=self.bOinsilica, dmin=0.0,
                                         dmax=0.19)  
            silicasurface.generate_bonds(name_a=self.SIinsilica, name_b=self.sOinsilica, dmin=0.0,
                                         dmax=0.19)  
            silicasurface.generate_bonds(name_a=self.sOinsilica, name_b=self.sHinsilica, dmin=0.0,
                                         dmax=0.120)  

        return silicasurface

