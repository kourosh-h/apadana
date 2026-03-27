# 🧬 Apadana v3.0

A specialized tool for automated construction of functionalized silica surfaces.

## 🛠 Quick Installation
1. conda create -n apadana_env python=3.10 -y
2. conda activate apadana_env
3. conda install -c mosdef -c conda-forge mbuild=0.16.0 foyer=0.12.1 parmed openbabel -y
4. pip install setuptools==69.5.1
5. pip install .

## 📖 Features
- **Surface Modes:** 'cell' (unit cell tiling), 'slab' (default relaxed), 'user' (custom PDB/GRO).
- **Automation:** Random grafting, port placement, and GROMACS topology generation.

## 📁 Usage Example
from apadana.functionalizer import functionalizer
functionalizer(chainLength=4, coverage=0.5, surface_mode='slab', namefuncsilica='test')