"""
APADANA: Surface Functionalization Test Script
----------------------------------------------
This script demonstrates how to graft alkylsilane chains onto a silica surface
using the Apadana package. It allows users to choose between building a 
new surface or using a pre-relaxed one.
"""

import sys
from pathlib import Path

# --- 1. SETUP ENVIRONMENT ---
# Ensure the program can find the 'apadana' package in the current directory
current_dir = Path(__file__).parent
sys.path.append(str(current_dir))

from apadana.functionalizer import functionalizer

# --- 2. DEFINE SIMULATION PARAMETERS ---
# Modify these variables to change the output of your functionalization
params = {
    "x_size": 1,              # Number of unit cell repeats in X (only for 'cell' mode)
    "y_size": 1,              # Number of unit cell repeats in Y (only for 'cell' mode)
    "chainLength": 4,         # Length of the alkyl chain (number of carbons)
    "coverage": 0.5,          # Surface coverage (0.5 = 50% of available sites)
    "free_Space": 6,          # Vacuum space in nm above/below the silica slab
    "width_mid": 1,           # Thickness of the middle layer (Angstrom) for restraints
    "si_res_name": "SUR",     # Residue name for the silica surface in the .gro file
    "alkylsilaneresname": "OTS", # Residue name for the grafted molecules
    "namefuncsilica": "func" # Base name for output files (func.gro, func.top, etc.)
}

# --- 3. CHOOSE SURFACE MODE ---
# Option 'cell': Build a new surface from the default unit cell (tiling is enabled)
# Option 'slab': Use the pre-relaxed slab included in the package (no tiling)
# Option 'user': Load your own custom .gro and .pdb files
mode_selection = 'slab' 

# Paths for 'user' mode (only used if surface_mode='user')
user_files = {
    "user_pdb": "./my_custom_surface.pdb"
}

# --- 4. RUN FUNCTIONALIZATION ---
print("="*60)
print(f"STARTING APADANA TEST: Mode = {mode_selection.upper()}")
print("="*60)

try:
    functionalizer(
        x_size=params["x_size"],
        y_size=params["y_size"],
        chainLength=params["chainLength"],
        coverage=params["coverage"],
        free_Space=params["free_Space"],
        width_mid=params["width_mid"],
        si_res_name=params["si_res_name"],
        alkylsilaneresname=params["alkylsilaneresname"],
        namefuncsilica=params["namefuncsilica"],
        surface_mode=mode_selection
    )
    #user_pdb=user_files["user_pdb"]
    print("\nSUCCESS: All files have been generated in the current directory.")
    print(f"Check '{params['namefuncsilica']}.gro' for the final structure.")

except Exception as e:
    print(f"\nERROR: Simulation failed. See details below:")
    print(f"--> {e}")

print("="*60)