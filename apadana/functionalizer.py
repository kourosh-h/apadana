from .fileReading import *
from .fileWriting import *
from .forcefields import *
from .modellingSurface import *
from .modellingMolecules import *
from mbuild.lib.atoms import H

import pandas as pd
from parmed import gromacs
import random
import numpy as np
from pathlib import Path

# Set up package paths
apadana_package_directory = Path(__file__).parent
interface_charmm_xml = apadana_package_directory / 'EMAMI_CHARMM.xml'

def cal_coverage(coverage_num: float, surface_object):
    """Calculate the number of chains based on surface area and fractional coverage."""
    number_of_chains = coverage_num * (surface_object.box.Lx * surface_object.box.Ly)
    return round(number_of_chains)

def report_sizes(surface_object, chainL, cover, n_grafted):
    """Write simulation box details and grafting results to a text file."""
    with open("report_sizes.txt" , 'w+') as rep:
        rep.writelines(f"boxX = {surface_object.box.Lx}\n")
        rep.writelines(f"boxY = {surface_object.box.Ly}\n")
        rep.writelines(f"boxZ = {surface_object.box.Lz}\n")
        rep.writelines(f"coverage = {cover}\n")
        rep.writelines(f"chainLength = {chainL}\n")
        rep.writelines(f"number of grafted chains = {n_grafted}\n")
    return None

def functionalizer(x_size,
                   y_size,
                   chainLength,
                   coverage,
                   free_Space,
                   width_mid,
                   si_res_name,
                   alkylsilaneresname,
                   namefuncsilica,
                   surface_mode='slab', 
                   user_gro=None, 
                   user_pdb=None):
    
    print("\n" + "-"*50)
    print("STEP 1: INITIALIZING MOLECULES AND SURFACE")
    print("-"*50)
    
    # Initialize alkylsilane builders for top and bottom grafting
    ALK = Alkysilane(alk_length=(chainLength - 1), res_name=alkylsilaneresname)
    OTS_GRAFT = ALK.graft_Top()
    OTS_GRAFT_B = ALK.graft_Botm()
    hydrogen = H()
    
    # Initialize the Surface Builder
    SSur = AlphaQuartz(x_blocks=x_size, y_blocks=y_size, freespace=free_Space, res_name=si_res_name)
    
    # Load or build the silica surface
    print(f"Loading surface in mode: {surface_mode}...")
    Silica_Surface = SSur.make_silica(mode=surface_mode, user_gro=user_gro, user_pdb=user_pdb)
    
    # Calculate coverage requirements
    n_grafted_chains = cal_coverage(coverage, Silica_Surface)
    print(f"Targeting {coverage*100}% coverage.")
    print(f"Number of grafted chains needed per side: {n_grafted_chains}")

    print("\n" + "-"*50)
    print("STEP 2: PREPARING TOP SURFACE FOR GRAFTING")
    print("-"*50)
    
    # Identify surface Hydrogen atoms (HD) on the upper half of the box
    h_top_atoms = []
    for atom1 in list(Silica_Surface.particles()):
        if atom1.name == "HD" and atom1.pos[2] > Silica_Surface.box.Lz / 2:
            h_top_atoms.append(atom1)
    print(f"Found {len(h_top_atoms)} available hydroxyl sites on Top surface.")

    # Randomly select sites and replace them with grafting ports (OP)
    print("Selecting random sites and adding Top ports...")
    for atom in random.sample(h_top_atoms, k=n_grafted_chains):
        for oxygen in Silica_Surface.bond_graph.neighbors(atom):
            oxygen.name = 'OP'
        Silica_Surface.remove(atom)

    # Attach mbuild Ports to the selected dangling oxygens
    for atom in list(Silica_Surface.particles()):
        if atom.name == "OP" and atom.pos[2] > Silica_Surface.box.Lz / 2:
            port1 = Port(anchor=atom)
            port1.spin(np.pi / 2, [1, 0, 0])
            port1.translate(np.array([0.0, 0.0, 0.1]))
            Silica_Surface.add(port1, f"port1_{len(Silica_Surface.referenced_ports())}")

    # Build the first monolayer (Top side)
    print("Applying Top monolayer...")
    m1Surface = Monolayer(chains=[OTS_GRAFT], backfill=hydrogen, surface=Silica_Surface, tile_x=1, tile_y=1)

    print("\n" + "-"*50)
    print("STEP 3: PREPARING BOTTOM SURFACE FOR GRAFTING")
    print("-"*50)

    # Identify surface Hydrogen atoms (HD) on the lower half of the box
    h_bottom_atoms = []
    for atom2 in list(m1Surface.particles()):
        if atom2.name == "HD" and atom2.pos[2] < Silica_Surface.box.Lz / 2:
            h_bottom_atoms.append(atom2)
    print(f"Found {len(h_bottom_atoms)} available hydroxyl sites on Bottom surface.")

    # Randomly select sites and replace them with grafting ports (OP)
    print("Selecting random sites and adding Bottom ports...")
    for atom in random.sample(h_bottom_atoms, k=n_grafted_chains):
        for oxygen in m1Surface.bond_graph.neighbors(atom):
            oxygen.name = 'OP'
        m1Surface.remove(atom)

    # Attach mbuild Ports to the selected dangling oxygens
    for atom in list(m1Surface.particles()):
        if atom.name == "OP" and atom.pos[2] < Silica_Surface.box.Lz / 2:
            port2 = Port(anchor=atom)
            port2.spin(-np.pi / 2, [1, 0, 0])
            port2.translate(np.array([0.0, 0.0, -0.1]))
            m1Surface.add(port2, f"port2_{len(m1Surface.referenced_ports())}")

    # Build the second monolayer (Bottom side)
    print("Applying Bottom monolayer...")
    func_surface = Monolayer(chains=[OTS_GRAFT_B], backfill=hydrogen, surface=m1Surface, tile_x=1, tile_y=1)

    print("\n" + "-"*50)
    print("STEP 4: SAVING COORDINATES AND TOPOLOGY")
    print("-"*50)

    # Save coordinate files (.pdb, .mol2, .gro)
    print(f"Exporting files: {namefuncsilica}.pdb/mol2/gro")
    func_surface.save(f'{namefuncsilica}.pdb', overwrite=True, residues=[alkylsilaneresname, SSur.ressilica])
    func_surface.save(f'{namefuncsilica}.mol2', overwrite=True, residues=[alkylsilaneresname, SSur.ressilica])
    func_surface.save(f'{namefuncsilica}.gro', overwrite=True, residues=[alkylsilaneresname, SSur.ressilica])
    
    # Save GROMACS topology using foyer and the CHARMM XML
    print("Applying forcefield and generating GROMACS .top file...")
    func_surface.save(f'{namefuncsilica}.top',
                      forcefield_files=interface_charmm_xml, overwrite=True,
                      box=func_surface.box, foyer_kwargs={"assert_bond_params": True, "assert_angle_params": True,
                                                          "assert_dihedral_params": True}, forcefield_debug=True,
                      residues=[alkylsilaneresname, SSur.ressilica])

    # Call custom writers for ITP and primary TOP files
    itpwriter(namefuncsilica)
    topologywriter(namefuncsilica)

    print("\n" + "-"*50)
    print("STEP 5: POST-PROCESSING INDEXES AND RESTRAINTS")
    print("-"*50)

    # Mapping atom indexes for GROMACS index groups
    ots_names = [i.name for i in OTS_GRAFT.particles()]
    surface_names = [i.name for i in Silica_Surface.particles()]
    itp = make_dic(f'{namefuncsilica}.itp')
    SUR = []
    OTS = []
    OTSS = [[] for k in range(n_grafted_chains * 2)]

    print("Mapping residues and atom indexes...")
    for i, l in enumerate(itp['[ atoms ]']):
        if l.split()[4] in surface_names:
            SUR.append(int(l.split()[5]))
            resSUR = int(l.split()[2])
        if l.split()[4] in ots_names:
            OTS.append(int(l.split()[5]))
            resOTS = (int(l.split()[2]) - resSUR) - 1
            OTSS[resOTS].append(int(l.split()[5]))
            
    # Locate the middle layer of the silica for position restraints
    domain = (width_mid/2)*10
    gmx_gro = gromacs.GromacsGroFile.parse(f'{namefuncsilica}.gro')
    gro = pd.DataFrame({'index':[gmx_gro.atoms[i].idx for i in range(len(gmx_gro.atoms))],
                        'name': [gmx_gro.atoms[i].name for i in range(len(gmx_gro.atoms))], 
                        'res': [gmx_gro.atoms[i].residue.name for i in range(len(gmx_gro.atoms))],
                       'x_axis': [gmx_gro.atoms[i].xx for i in range(len(gmx_gro.atoms))], 
                       'y_axis': [gmx_gro.atoms[i].xy for i in range(len(gmx_gro.atoms))], 
                       'z_axis': [gmx_gro.atoms[i].xz for i in range(len(gmx_gro.atoms))]})

    surz = [l for i,l in enumerate(gro['z_axis']) if gro['res'][i]=='SUR']
    midpoint = ((max(surz)-min(surz))/2)+min(surz)
    mid_layer = [i for i in gro['index'] if gro['res'][i]=='SUR' and midpoint-domain < gro['z_axis'][i] < midpoint+domain]

    # Write final simulation support files
    print("Writing index file, posre file, and gromacs shell script...")
    indexwriter(SUR, OTS, OTSS, mid_layer)
    posreswriter(SUR, OTS, mid_layer)
    gromacssh(namefuncsilica)
    report_sizes(func_surface, chainLength, coverage, n_grafted_chains)
    
    print("\n" + "="*50)
    print(f"FUNCTIONALIZATION SUCCESSFUL: {namefuncsilica} is ready!")
    print("="*50 + "\n")