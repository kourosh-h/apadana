import mbuild as mb
from mbuild.lib.recipes import Monolayer
from mbuild.lib.recipes import Alkane
from mbuild.lib.moieties import Silane
from mbuild.lib.atoms import H
from mbuild.lib.recipes import TiledCompound
from mbuild import Port
#from foyer import Forcefield
# from foyer.examples.utils import example_file_path
import parmed as pmd
from parmed import gromacs, amber, unit as u
#import openbabel
# import pybel
#import py3Dmol
# import nglview as nv
#from ase import Atoms
#from ase.io.trajectory import Trajectory
#from ase.io import read
# import psi4
# import resp
import os
import random
import numpy as np
import pandas as pd
from IPython.display import display, Latex
# import textract
# import pdfquery
# from pdfquery import PDFQuery
# import PyPDF2


class Hydroxyl(mb.Compound):
    def __init__(self, i=1, j=1):
        super(Hydroxyl, self).__init__()
        oxygen = mb.Particle(name=f'O{i}', element='O', pos=[0, 0, 0])
        self.add(oxygen, label='OH')
        hydrogen = mb.Particle(name=f'H{j}', element='H', pos=[0.120, 0, 0])
        self.add(hydrogen, label='HO')
        self.add_bond((self[0], self['HO'][0]))
        port = mb.Port(anchor=self[0], orientation=[0, 1, 0], )
        self.add(port, label='up')
        self['up'].translate([0, 0.06, 0])

        # port2 = mb.Port(anchor=self[0], orientation=[0,1,0],)
        # self.add(port2, label='down')
        # self['down'].rotate(np.pi/2, [1, 0, 0])
        # self['down'].translate([0, -0.06, 0])


# ----------------------------------------------------------------------
class SiOH2old(mb.Compound):
    def __init__(self):
        super(SiOH2, self).__init__()
        self.add(mb.Particle(name='SI1', element='Si', charge=0.78, pos=[0, 0, 0]), label='SI[$]')
        self.add(mb.Particle(name='O1', element='O', pos=[0.190, 0, 0]), label='OH[$]')
        self.add(mb.Particle(name='O2', element='O', pos=[-0.190, 0, 0]), label='OH[$]')
        self.add(mb.Particle(name='H1', element='H', pos=[0.190, 0, 0.120]), label='HO[$]')
        self.add(mb.Particle(name='H2', element='H', pos=[-0.190, 0, -0.120]), label='HO[$]')

        self.add_bond((self[0], self['OH'][0]))
        self.add_bond((self[0], self['OH'][1]))
        self.add_bond((self['OH'][0], self['HO'][0]))
        self.add_bond((self['OH'][1], self['HO'][1]))
        theta = 0.5 * (180 - 109.5) * np.pi / 180
        self.add(mb.Port(anchor=self[0]), label='up')

        self['up'].translate([0, 0.095, 0])
        self['up'].rotate(theta, around=[1, 0, 0])

        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, -0.095, 0])
        self['down'].rotate(np.pi, [0, 1, 0])
        # self['down'].rotate(-theta, around=[1, 0, 0])


# ----------------------------------------------------------------------

class SiOH2(mb.Compound):
    def __init__(self):
        super(SiOH2, self).__init__()
        self.add(mb.Particle(name='SI1', element='Si', charge=0.78, pos=[0, 0, 0]), label='SI[$]')
        self.add(mb.Particle(name='O1', element='O', pos=[0.180, 0, 0]), label='OH[$]')
        self.add(mb.Particle(name='O2', element='O', pos=[-0.180, 0, 0]), label='OH[$]')
        self.add(mb.Particle(name='H1', element='H', pos=[0.192, 0, -0.11]), label='HO[$]')
        self.add(mb.Particle(name='H2', element='H', pos=[-0.192, 0, -0.11]), label='HO[$]')

        self.add_bond((self[0], self['OH'][0]))
        self.add_bond((self[0], self['OH'][1]))
        self.add_bond((self['OH'][0], self['HO'][0]))
        self.add_bond((self['OH'][1], self['HO'][1]))
        theta = 0.5 * (180 - 109.5) * np.pi / 180
        self['OH'][0].rotate(theta, around=[0, 1, 0.4])
        self['OH'][1].rotate(-theta, around=[0, 1, 0.4])

        self['HO'][0].rotate(theta, around=[0, 1, 0.4])
        self['HO'][1].rotate(-theta, around=[0, 1, 0.4])

        self.add(mb.Port(anchor=self[0]), label='up')

        self['up'].translate([0, 0.098, 0])
        self['up'].rotate(theta, around=[1, 0, 0])

        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, -0.095, 0])
        self['down'].rotate(np.pi, [0, 1, 0])
        # self['down'].rotate(-theta, around=[1, 0, 0])


# ----------------------------------------------------------------------
class CH2(mb.Compound):
    def __init__(self, i=1, j=1):
        super(CH2, self).__init__()
        # Add carbon
        self.add(mb.Particle(name=f'C{i}', element='C', pos=[0, 0, 0]), label='CT[$]')
        # self['CT'][0].element= 'C'

        # Add hydrogens
        self.add(mb.Particle(name=f'H{j}', element='H', pos=[-0.109, 0, 0.0]), label='HC[$]')
        self.add(mb.Particle(name=f'H{j + 1}', element='H', pos=[0.109, 0, 0.0]), label='HC[$]')

        theta = 0.5 * (180 - 109.5) * np.pi / 180
        self['HC'][0].rotate(-theta, around=[0, 1, 0])
        self['HC'][1].rotate(theta, around=[0, 1, 0])

        # Add bonds between the atoms
        self.add_bond((self['CT'][0], self['HC'][0]))
        self.add_bond((self['CT'][0], self['HC'][1]))

        # Add ports
        self.add(mb.Port(anchor=self[0]), label='up')
        self['up'].translate([0, 0.154 / 2, 0])
        self['up'].rotate(theta, around=[1, 0, 0])

        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, -0.154 / 2, 0])
        self['down'].rotate(np.pi, [0, 1, 0])
        self['down'].rotate(-theta, around=[1, 0, 0])
    # ----------------------------------------------------------------------


class CF2(mb.Compound):
    def __init__(self):
        super(CF2, self).__init__()
        # Add carbon
        self.add(mb.Particle(name='CF', element='C', pos=[0, 0, 0]), label='CT[$]')

        # Add hydrogens
        self.add(mb.Particle(name='F', element='F', pos=[-0.109, 0, 0.0]), label='FC[$]')
        self.add(mb.Particle(name='F', element='F', pos=[0.109, 0, 0.0]), label='FC[$]')

        theta = 0.5 * (180 - 109.5) * np.pi / 180
        self['FC'][0].rotate(-theta, around=[0, 1, 0])
        self['FC'][1].rotate(theta, around=[0, 1, 0])

        # Add bonds between the atoms
        self.add_bond((self['CF'][0], self['FC'][0]))
        self.add_bond((self['CF'][0], self['FC'][1]))

        # Add ports
        self.add(mb.Port(anchor=self[0]), label='up')
        self['up'].translate([0, 0.154 / 2, 0])
        self['up'].rotate(theta, around=[1, 0, 0])

        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, -0.154 / 2, 0])
        self['down'].rotate(np.pi, [0, 1, 0])
        self['down'].rotate(-theta, around=[1, 0, 0])
    # ----------------------------------------------------------------------


class CH3(mb.Compound):
    def __init__(self, i=1, j=1):
        super(CH3, self).__init__()
        theta = 0.5 * (180 - 109.5) * np.pi / 180

        carbon = mb.Particle(name=f'C{i}', element='C', pos=[0, 0, 0])
        self.add(carbon, label='CT[$]')

        self.add(mb.Particle(name=f'H{j}', element='H', pos=[0.1, 0, -0.07]), label='HC[$]')
        self.add(mb.Particle(name=f'H{j + 1}', element='H', pos=[-0.1, 0, -0.07]), label='HC[$]')
        self.add(mb.Particle(name=f'H{j + 2}', element='H', pos=[0, 0.1, 0.07]), label='HC[$]')

        self.add_bond((self[0], self['HC'][0]))
        self.add_bond((self[0], self['HC'][1]))
        self.add_bond((self[0], self['HC'][2]))

        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, -0.154 / 2, 0])
        self['down'].rotate(np.pi, [0, 1, 0])
        self['down'].rotate(-theta, around=[1, 0, 0])
    # ----------------------------------------------------------------------


class CF3(mb.Compound):
    def __init__(self):
        super(CF3, self).__init__()
        theta = 0.5 * (180 - 109.5) * np.pi / 180

        carbon = mb.Particle(name='CF', element='C', pos=[0, 0, 0])
        self.add(carbon, label='CF[$]')

        self.add(mb.Particle(name='F', element='F', pos=[0.1, 0, -0.07]), label='FC[$]')
        self.add(mb.Particle(name='F', element='F', pos=[-0.1, 0, -0.07]), label='FC[$]')
        self.add(mb.Particle(name='F', element='F', pos=[0, 0.1, 0.07]), label='FC[$]')

        self.add_bond((self[0], self['FC'][0]))
        self.add_bond((self[0], self['FC'][1]))
        self.add_bond((self[0], self['FC'][2]))

        self.add(mb.Port(anchor=self[0]), label='down')
        self['down'].translate([0, -0.154 / 2, 0])
        self['down'].rotate(np.pi, [0, 1, 0])
        self['down'].rotate(-theta, around=[1, 0, 0])
    # ----------------------------------------------------------------------


class AlkanePolymer(mb.Compound):
    def __init__(self, chain_length, i=1, j=4):
        super(AlkanePolymer, self).__init__()
        theta = 180 * np.pi / 180
        # i=1
        # j=4
        last_monomer = CH2(i, j)
        # last_monomer[0].name  = 'C1'
        # last_monomer['HC'][0].name  = 'H1'
        # last_monomer['HC'][1].name  = 'H2'
        self.add(last_monomer)
        self.add(last_monomer['down'], 'down', containment=False)
        for monomer in range(chain_length - 1):
            i += 1
            j += 2
            current_monomer = CH2(i, j)
            # current_monomer[0].name  = f'C{i+2}'
            # current_monomer['HC'][0].name  = f'H{i+3}'
            # current_monomer['HC'][1].name  = f'H{i+4}'
            last_monomer.rotate(theta, around=[0, 1, 0])
            mb.force_overlap(move_this=current_monomer,
                             from_positions=current_monomer['down'],
                             to_positions=last_monomer['up'])
            current_monomer.rotate(theta, around=[0, 1, 0])
            self.add(current_monomer)
            last_monomer = current_monomer
        last_monomer.rotate(theta, around=[0, 1, 0])

        self.add(last_monomer['up'], 'up', containment=False)
    # ----------------------------------------------------------------------


class FCPolymer(mb.Compound):
    def __init__(self, chain_length):
        super(FCPolymer, self).__init__()
        theta = 180 * np.pi / 180
        last_monomer = CH2()
        self.add(last_monomer)
        self.add(last_monomer['down'], 'down', containment=False)
        for i in range(chain_length - 1):
            if i < 1:
                current_monomer = CH2()
            else:
                current_monomer = CF2()
            last_monomer.rotate(theta, around=[0, 1, 0])
            mb.force_overlap(move_this=current_monomer,
                             from_positions=current_monomer['down'],
                             to_positions=last_monomer['up'])
            current_monomer.rotate(theta, around=[0, 1, 0])
            self.add(current_monomer)
            last_monomer = current_monomer
        last_monomer.rotate(theta, around=[0, 1, 0])

        self.add(last_monomer['up'], 'up', containment=False)
    # ----------------------------------------------------------------------


class Alkysilane0(mb.Compound):
    def __init__(self, chainLength=6, cap=True):
        super(Alkysilane0, self).__init__()
        if cap == True:
            silane = SiOH2()
            # silane[0].name ='SI1'
            silane.name = 'SIT'
            backbone = AlkanePolymer(chain_length=chainLength, i=1, j=4)
            backbone.name = 'BAK'
            # tail = H()
            i = chainLength + 1
            j = (chainLength * 2) + 1 + 3
            tail = CH3(i, j)
            tail.name = 'TIL'
            cap = Hydroxyl(i=3, j=3)
            cap.name = 'CAP'
            self.add(silane)
            mb.force_overlap(move_this=cap, from_positions=cap['up'], to_positions=silane['down'])
            self.add(cap)
            mb.force_overlap(move_this=backbone, from_positions=backbone['down'], to_positions=silane['up'])
            self.add(backbone)
            mb.force_overlap(move_this=tail, from_positions=tail['down'], to_positions=backbone['up'])
            self.add(tail)

        if cap == False:
            silane = SiOH2()
            # silane[0].name ='SI1'
            silane.name = 'SIT'
            backbone = AlkanePolymer(chain_length=chainLength, i=1, j=3)
            backbone.name = 'BAK'
            # tail = H()
            i = chainLength + 1
            j = (chainLength * 2) + 1 + 2
            tail = CH3(i, j)
            tail.name = 'TIL'
            self.add(silane)
            mb.force_overlap(move_this=backbone, from_positions=backbone['down'], to_positions=silane['up'])
            self.add(backbone)
            mb.force_overlap(move_this=tail, from_positions=tail['down'], to_positions=backbone['up'])
            self.add(tail)
            self.add(silane['down'], 'down', containment=False)
        # ----------------------------------------------------------------------


class Alkysilane():
    def __init__(self, alk_length=6, res_name='OTS'):
        self.alk_length = alk_length
        self.res_name = res_name
        # self.cap = cap

    def complete(self):
        Alkylsilane_cap = Alkysilane0(self.alk_length, cap=True)
        Alkylsilane_cap.name = self.res_name
        Alkylsilane_cap.resname = self.res_name
        # Alkylsilane_pybel = Alkylsilane_cap.to_pybel(residues=self.res_name)
        # Alkylsilane_pybel.localopt()
        # Alkylsilane_pybel1 = mb.Compound()
        # Alkylsilane_pybel1.from_pybel(pybel_mol= Alkylsilane_pybel)
        return Alkylsilane_cap  # Alkylsilane_pybel1

    def no_H(self):
        Alkylsilane_noH = Alkysilane0(self.alk_length, cap=True)
        Alkylsilane_noH.name = self.res_name
        Alkylsilane_noH.resname = self.res_name
        HL = Alkylsilane_noH.particles_by_name('H3')
        Alkylsilane_noH.remove(HL)
        pybel1 = Alkylsilane_noH.to_pybel(residues=self.res_name)
        pybel1.localopt()
        Alkylsilane_pybel1 = mb.Compound()
        Alkylsilane_pybel1.from_pybel(pybel_mol=pybel1)
        return Alkylsilane_pybel1

    def graft_Top(self):
        Alkylsilane_grafting = Alkysilane0(self.alk_length, cap=False)
        Alkylsilane_grafting.name = self.res_name
        Alkylsilane_grafting.resname = self.res_name
        return Alkylsilane_grafting

    def graft_Botm(self):
        Alkylsilane_grafting_B = self.graft_Top()
        Alkylsilane_grafting_B.spin(np.pi, [1, 0, 0])
        return Alkylsilane_grafting_B

    def savefiles(self, molecule, file_name):
        molecule.save(f'{file_name}.pdb', overwrite=True, residues=self.res_name)
        molecule.save(f'{file_name}.mol2', overwrite=True, residues=self.res_name)
        molecule.save(f'{file_name}.xyz', overwrite=True, residues=self.res_name)

    # ----------------------------------------------------------------------


