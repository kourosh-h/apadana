from setuptools import setup, find_packages


setup(
    name='apadana',
    version='3.0',
    packages=find_packages(),
    package_data={
        'apadana': [
            'README.md',
            'EMAMI_CHARMM.xml',
            'EMAMI.xml',
            'LICENSE',
            'charmm-gui-0720347578/*', 'charmm-gui-0720347578/gromacs/*', 'charmm-gui-0720347578/openmm/*',
            'charmm36-jul2022.ff/*',
            'additional_files/*',
        ],
    },
    include_package_data=True,
    description='Apadana makes functionalized silica surface with alkylsilanes',
    author='Kourosh Hasheminejad',
    author_email='kourosh.hasheminejad@gmail.com',
    #url='https://github.com/yourusername/mypackage',
    install_requires=[
        "mbuild==0.16.0",
        "foyer==0.12.1",
        "setuptools==69.5.1",
        "parmed",
        "pandas",
        "numpy",
        "py3Dmol"
    ],
    python_requires=">=3.9, <3.11",
)
